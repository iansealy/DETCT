#!/usr/bin/env perl

# PODNAME: detag_fastq.pl
# ABSTRACT: Extract tags from transcript counting FASTQ files and process files

## Author         : is1
## Maintainer     : is1
## Created        : 2012-12-15
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
use DETCT::Misc::Tag;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Default options
## no critic (ProhibitMagicNumbers)
my $fastq_read1_input;
my $fastq_read2_input;
my $fastq_output_prefix;
my $read1_required_length    = 30;
my $read2_required_length    = 54;
my $read1_5prime_trim_length = 0;
my $read2_5prime_trim_length = 0;
my $polyt_trim_length        = 14;
my $polyt_min_length         = 10;
my @read_tags;
my $no_pair_suffix        = 0;
my $treat_n_in_polyt_as_t = 0;
my ( $help, $man );
## use critic

# Get and check command line options
get_and_check_options();

DETCT::Misc::Tag::detag_trim_fastq(
    {
        fastq_read1_input        => $fastq_read1_input,
        fastq_read2_input        => $fastq_read2_input,
        fastq_output_prefix      => $fastq_output_prefix,
        read1_required_length    => $read1_required_length,
        read2_required_length    => $read2_required_length,
        read1_5prime_trim_length => $read1_5prime_trim_length,
        read2_5prime_trim_length => $read2_5prime_trim_length,
        polyt_trim_length        => $polyt_trim_length,
        polyt_min_length         => $polyt_min_length,
        read_tags                => \@read_tags,
        no_pair_suffix           => $no_pair_suffix,
        treat_n_in_polyt_as_t    => $treat_n_in_polyt_as_t,
    }
);

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'fastq_read1_input=s'        => \$fastq_read1_input,
        'fastq_read2_input=s'        => \$fastq_read2_input,
        'fastq_output_prefix=s'      => \$fastq_output_prefix,
        'read1_required_length=i'    => \$read1_required_length,
        'read2_required_length=i'    => \$read2_required_length,
        'read1_5prime_trim_length=i' => \$read1_5prime_trim_length,
        'read2_5prime_trim_length=i' => \$read2_5prime_trim_length,
        'polyt_trim_length=i'        => \$polyt_trim_length,
        'polyt_min_length=i'         => \$polyt_min_length,
        'read_tags=s@{1,}'           => \@read_tags,
        'no_pair_suffix'             => \$no_pair_suffix,
        'treat_n_in_polyt_as_t'      => \$treat_n_in_polyt_as_t,
        'help'                       => \$help,
        'man'                        => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$fastq_read1_input ) {
        pod2usage("--fastq_read1_input must be specified\n");
    }
    if ( !$fastq_read2_input ) {
        pod2usage("--fastq_read2_input must be specified\n");
    }
    if ( !$fastq_output_prefix ) {
        pod2usage("--fastq_output_prefix must be specified\n");
    }
    if ( !@read_tags ) {
        pod2usage("--read_tags must be specified\n");
    }

    return;
}

=head1 USAGE

    detag_fastq.pl
        [--fastq_read1_input file]
        [--fastq_read2_input file]
        [--fastq_output_prefix prefix]
        [--read1_required_length int]
        [--read2_required_length int]
        [--read1_5prime_trim_length int]
        [--read2_5prime_trim_length int]
        [--polyt_trim_length int]
        [--polyt_min_length int]
        [--read_tags tags...]
        [--no_pair_suffix]
        [--treat_n_in_polyt_as_t]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--fastq_read1_input FILE>

Input FASTQ file for read 1.

=item B<--fastq_read2_input FILE>

Input FASTQ file for read 2.

=item B<--fastq_output_prefix FILE>

Prefix for output FASTQ files.

=item B<--read1_required_length INT>

Length of read 1 after detagging and optional trimming.

=item B<--read2_required_length INT>

Length of read 2 after optional trimming.

=item B<--read1_5prime_trim_length INT>

Bases to trim from 5' end of read 1 after detagging and before any 3' trimming.

=item B<--read2_5prime_trim_length INT>

Bases to trim from 5' end of read 2 before any 3' trimming.

=item B<--polyt_trim_length INT>

Length of (largely) polyT to be trimmed.

=item B<--polyt_min_length INT>

Minimum number of consecutive Ts in length of polyT.

=item B<--read_tags TAGS>

Read tags.

=item B<--no_pair_suffix>

Input FASTQ file don't have pair suffixes.

=item B<--treat_n_in_polyt_as_t>

Assume any Ns in polyT are Ts.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
