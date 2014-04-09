#!/usr/bin/env perl

# PODNAME: mark_duplicates.pl
# ABSTRACT: Mark duplicates for a BAM file

## Author         : is1
## Maintainer     : is1
## Created        : 2014-04-02
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Getopt::Long;
use Pod::Usage;
use File::Slurp;
use File::Spec;
use File::Path qw( make_path );
use DETCT::Misc::BAM;
use DETCT::Misc::SAMtools;
use DETCT::Misc::Picard;

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $dir    = q{.};
my $method = 'native';
my $input_bam_file;
my $output_bam_file;
my $consider_tags;
my $fix_mate_info;
my $samtools_binary     = 'samtools';
my $java_binary         = 'java';
my $sort_bam_jar        = 'picard-tools-1.110-detct/SortSam.jar';
my $mark_duplicates_jar = 'picard-tools-1.110-detct/MarkDuplicates.jar';
my $fix_mate_info_jar   = 'picard-tools-1.110-detct/FixMateInformation.jar';
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Ensure working directory exists
make_path($dir);

# Intermediate files
my $sorted_bam_file         = File::Spec->catfile( $dir, 'sorted.bam' );
my $fixmateinfo_bam_file    = File::Spec->catfile( $dir, 'fixmateinfo.bam' );
my $markduplicates_bam_file = File::Spec->catfile( $dir, 'markduplicates.bam' );

my $metrics;

# Delete intermediate files if necessary
for
  my $file ( $sorted_bam_file, $fixmateinfo_bam_file, $markduplicates_bam_file )
{
    if ( -e $file ) {
        unlink $file;
    }
}

if ( $method eq 'native' ) {

    # Sort BAM file by read name
    DETCT::Misc::Picard::sort_bam(
        {
            dir             => $dir,
            input_bam_file  => $input_bam_file,
            output_bam_file => $sorted_bam_file,
            java_binary     => $java_binary,
            sort_bam_jar    => $sort_bam_jar,
            sort_order      => 'queryname',
        }
    );

    if ($fix_mate_info) {

        # Fix mate information
        DETCT::Misc::Picard::fix_mate_info(
            {
                dir               => $dir,
                input_bam_file    => $sorted_bam_file,
                output_bam_file   => $fixmateinfo_bam_file,
                java_binary       => $java_binary,
                fix_mate_info_jar => $fix_mate_info_jar,
            }
        );
    }

    my $markduplicates_input_bam_file =
      $fix_mate_info ? $fixmateinfo_bam_file : $sorted_bam_file;

    # Mark duplicates
    $metrics = DETCT::Misc::BAM::mark_duplicates(
        {
            input_bam_file  => $markduplicates_input_bam_file,
            output_bam_file => $markduplicates_bam_file,
            consider_tags   => $consider_tags,
        }
    );

    # Sort BAM file by coordinate
    DETCT::Misc::SAMtools::sort_bam(
        {
            dir             => $dir,
            input_bam_file  => $markduplicates_bam_file,
            output_bam_file => $output_bam_file,
            samtools_binary => $samtools_binary,
            sort_order      => 'coordinate',
        }
    );
}
elsif ( $method eq 'picard' ) {
    if ($fix_mate_info) {

        # Fix mate information
        DETCT::Misc::Picard::fix_mate_info(
            {
                dir               => $dir,
                input_bam_file    => $input_bam_file,
                output_bam_file   => $fixmateinfo_bam_file,
                java_binary       => $java_binary,
                fix_mate_info_jar => $fix_mate_info_jar,
            }
        );
    }

    my $markduplicates_input_bam_file =
      $fix_mate_info ? $fixmateinfo_bam_file : $input_bam_file;

    DETCT::Misc::Picard::mark_duplicates(
        {
            dir                 => $dir,
            input_bam_file      => $markduplicates_input_bam_file,
            output_bam_file     => $output_bam_file,
            metrics_file        => $output_bam_file . '.metrics',
            java_binary         => $java_binary,
            mark_duplicates_jar => $mark_duplicates_jar,
            consider_tags       => $consider_tags,
        }
    );

    $metrics = DETCT::Misc::Picard::extract_mark_duplicates_metrics(
        {
            metrics_file => $output_bam_file . '.metrics',
        }
    );
}

# Index BAM file
DETCT::Misc::SAMtools::make_index(
    {
        dir             => $dir,
        bam_file        => $output_bam_file,
        samtools_binary => $samtools_binary,
    }
);

# Stats
DETCT::Misc::SAMtools::flagstats(
    {
        dir             => $dir,
        bam_file        => $output_bam_file,
        samtools_binary => $samtools_binary,
        output_file     => $output_bam_file . '.flagstat',
    }
);

# Delete intermediate files
unlink $sorted_bam_file, $fixmateinfo_bam_file, $markduplicates_bam_file;

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'dir=s'                 => \$dir,
        'method=s'              => \$method,
        'input_bam_file=s'      => \$input_bam_file,
        'output_bam_file=s'     => \$output_bam_file,
        'consider_tags'         => \$consider_tags,
        'fix_mate_info'         => \$fix_mate_info,
        'samtools_binary=s'     => \$samtools_binary,
        'java_binary=s'         => \$java_binary,
        'sort_bam_jar=s'        => \$sort_bam_jar,
        'mark_duplicates_jar=s' => \$mark_duplicates_jar,
        'fix_mate_info_jar=s'   => \$fix_mate_info_jar,
        'debug'                 => \$debug,
        'help'                  => \$help,
        'man'                   => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$input_bam_file ) {
        pod2usage("--input_bam_file must be specified\n");
    }
    if ( !$output_bam_file ) {
        pod2usage("--output_bam_file must be specified\n");
    }
    if ( $method ne 'native' && $method ne 'picard' ) {
        pod2usage("--method must be 'native' or 'picard'\n");
    }

    return;
}

=head1 USAGE

    mark_duplicates.pl
        [--dir directory]
        [--method string]
        [--input_bam_file file]
        [--output_bam_file file]
        [--consider_tags]
        [--fix_mate_info]
        [--samtools_binary file]
        [--java_binary file]
        [--sort_bam_jar file]
        [--mark_duplicates_jar file]
        [--fix_mate_info_jar file]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--dir DIRECTORY>

Working directory for analysis.

=item B<--method native|picard>

Mark duplicates method to use - native (default) or picard.

=item B<--input_bam_file FILE>

Input BAM file.

=item B<--output_bam_file FILE>

Output BAM file.

=item B<--consider_tags>

Consider tag (part after # in read name) when deciding if reads are duplicates.

=item B<--fix_mate_info>

Fix mate information before marking duplicates.

=item B<--samtools_binary FILE>

SAMtools binary.

=item B<--java_binary FILE>

Java binary.

=item B<--sort_bam_jar FILE>

Picard SortSam JAR file.

=item B<--mark_duplicates_jar FILE>

Picard MarkDuplicates JAR file.

=item B<--fix_mate_info_jar FILE>

Picard FixMateInformation JAR file.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
