#!/usr/bin/env perl

# PODNAME: filter_by_tag.pl
# ABSTRACT: Filter BAM file by tag(s)

## Author         : is1
## Maintainer     : is1
## Created        : 2014-09-19
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
use DETCT::Misc::BAM;

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $input_bam_file;
my $output_bam_file;
my @read_tags;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

DETCT::Misc::BAM::filter_by_tag(
    {
        source_bam_file => $input_bam_file,
        target_bam_file => $output_bam_file,
        tags            => \@read_tags,
    }
);

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'input_bam_file=s'  => \$input_bam_file,
        'output_bam_file=s' => \$output_bam_file,
        'read_tags=s@{1,}'  => \@read_tags,
        'debug'             => \$debug,
        'help'              => \$help,
        'man'               => \$man,
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
    if ( !@read_tags ) {
        pod2usage("--read_tags must be specified\n");
    }

    return;
}

=head1 USAGE

    filter_by_tag.pl
        [--input_bam_file file]
        [--output_bam_file file]
        [--read_tags tags...]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--input_bam_file FILE>

Input BAM file.

=item B<--output_bam_file FILE>

Output BAM file.

=item B<--read_tags TAGS>

Read tags.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
