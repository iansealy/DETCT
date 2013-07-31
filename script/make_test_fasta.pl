#!/usr/bin/env perl

# PODNAME: make_test_fasta.pl
# ABSTRACT: Make transcript counting test file in FASTA format

## Author         : is1
## Maintainer     : is1
## Created        : 2012-11-12
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

=head1 DESCRIPTION

This script generates test transcript counting FASTA files. The number and
maximum length of chromosomes can be varied.

=head1 EXAMPLES

    # Generate random FASTA file using default values
    perl script/make_test_fasta.pl > test.fa

    # Generate FASTA file with reproducible chromosomes using default values
    perl script/make_test_fasta.pl --seed 1 > test.fa

    # Generate FASTA file with 25 chromosomes (each up to 50 Mbp long)
    perl script/make_test_fasta.pl \
        --seq_region_count 25 \
        --seq_region_max_length 50_000_000 \
        > test.fa

=cut

# Default options
## no critic (ProhibitMagicNumbers)
my $seed;
my $seq_region_count      = 1;
my $seq_region_max_length = 1_000_000;
my ( $help, $man );
## use critic

# Get and check command line options
get_and_check_options();

# Ensure reproducible chromosome lengths if seed set
if ( defined $seed ) {
    srand $seed;
}

# Make each chromosome of random length
my %length_of;
foreach my $seq_region ( 1 .. $seq_region_count ) {
    my $length = int rand( $seq_region_max_length + 1 );
    $length_of{$seq_region} = $length;
}

# Ensure sequences are always random
srand;

# Generate sequence for each chromosome one by one
foreach my $seq_region ( 1 .. $seq_region_count ) {
    printf ">%s\n", $seq_region;
    my $length_required = $length_of{$seq_region};
    my $length_printed  = 0;
    while ($length_required) {
        ## no critic (ProhibitMagicNumbers)
        print qw( A G C T a g c t ) [ int rand 8 ];
        ## use critic
        $length_required--;
        $length_printed++;

        # Wrap every 80 bases
        ## no critic (ProhibitMagicNumbers)
        if ( !( $length_printed % 80 ) ) {
            ## use critic
            print "\n";
        }
    }

    # Final new line if haven't just printed one
    ## no critic (ProhibitMagicNumbers)
    if ( $length_printed % 80 ) {
        ## use critic
        print "\n";
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'seed=i'                  => \$seed,
        'seq_region_count=i'      => \$seq_region_count,
        'seq_region_max_length=i' => \$seq_region_max_length,
        'help'                    => \$help,
        'man'                     => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$seq_region_count ) {
        pod2usage("--seq_region_count must be a positive integer\n");
    }
    if ( !$seq_region_max_length ) {
        pod2usage("--seq_region_max_length must be a positive integer\n");
    }

    return;
}

=head1 USAGE

    make_test_fasta.pl
        [--seed seed]
        [--seq_region_count int]
        [--seq_region_max_length int]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--seed INT>

Random seed (to get reproducible chromosome lengths).

=item B<--seq_region_count INT>

Number of seq regions (default to 1).

=item B<--seq_region_max_length INT>

Maximum length of each seq region (defaults to 1,000,000 bp).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
