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
my $samtools_binary     = 'samtools';
my $java_binary         = 'java';
my $mark_duplicates_jar = 'picard-tools-1.110-detct/MarkDuplicates.jar';
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Ensure working directory exists
make_path($dir);

# Intermediate files
my $sorted_bam_file         = File::Spec->catfile( $dir, 'sorted.bam' );
my $markduplicates_bam_file = File::Spec->catfile( $dir, 'markduplicates.bam' );

if ( $method eq 'native' ) {

    # Delete intermediate files if necessary
    if ( -e $sorted_bam_file ) {
        unlink $sorted_bam_file;
    }
    if ( -e $markduplicates_bam_file ) {
        unlink $markduplicates_bam_file;
    }

    # Sort BAM file by read name
    DETCT::Misc::SAMtools::sort_bam(
        {
            dir             => $dir,
            input_bam_file  => $input_bam_file,
            output_bam_file => $sorted_bam_file,
            samtools_binary => $samtools_binary,
            sort_order      => 'queryname',
        }
    );

    # Mark duplicates
    DETCT::Misc::BAM::mark_duplicates(
        {
            input_bam_file  => $sorted_bam_file,
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

    # Delete intermediate files
    unlink $sorted_bam_file, $markduplicates_bam_file;
}
elsif ( $method eq 'picard' ) {
    DETCT::Misc::Picard::mark_duplicates(
        {
            dir                 => $dir,
            input_bam_file      => $input_bam_file,
            output_bam_file     => $output_bam_file,
            metrics_file        => $output_bam_file . '.metrics',
            java_binary         => $java_binary,
            mark_duplicates_jar => $mark_duplicates_jar,
            consider_tags       => $consider_tags,
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

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'dir=s'                 => \$dir,
        'method=s'              => \$method,
        'input_bam_file=s'      => \$input_bam_file,
        'output_bam_file=s'     => \$output_bam_file,
        'consider_tags'         => \$consider_tags,
        'samtools_binary=s'     => \$samtools_binary,
        'java_binary=s'         => \$java_binary,
        'mark_duplicates_jar=s' => \$mark_duplicates_jar,
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
        [--samtools_binary file]
        [--java_binary file]
        [--mark_duplicates_jar file]
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

=item B<--samtools_binary FILE>

SAMtools binary.

=item B<--java_binary FILE>

Java binary.

=item B<--mark_duplicates_jar FILE>

Picard MarkDuplicates JAR file.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
