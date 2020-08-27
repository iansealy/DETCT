#!/usr/bin/env perl

# PODNAME: make_test_fastq.pl
# ABSTRACT: Make transcript counting test files in FASTQ format

## Author         : is1
## Maintainer     : is1
## Created        : 2013-01-08
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
use DETCT::Misc qw( write_or_die );

=head1 DESCRIPTION

This script generates test transcript counting FASTQ files.

=head1 EXAMPLES

    # Generate random FASTQ files using default values
    perl -Ilib script/make_test_fastq.pl --read_tags NNNNCTACCA

    # Generate FASTQ files with reproducible reads using default values
    perl -Ilib script/make_test_fastq.pl --seed 1 

    # Generate random FASTQ files with 1000 read pairs and 54 bp reads
    perl -Ilib script/make_test_fastq.pl \
        --read_tags NNNNCTACCA \
        --read_pair_count 1000 \
        --read_length 54

=cut

# Default options
## no critic (ProhibitMagicNumbers)
my $seed;
my $output_prefix   = 'test';
my $read_pair_count = 100;
my @read_tags;
my $read_length  = 75;
my $polyt_length = 14;
my ( $help, $man );
## use critic

# Get and check command line options
get_and_check_options();

# Assume all tags are same length
my $tag_length = length $read_tags[0];

# Add dummy tag for reads that don't match a real tag
push @read_tags, q{X} x $tag_length;

# Ensure reproducible FASTQ files if seed set
if ( defined $seed ) {
    srand $seed;
}

# Generate start of each read name
## no critic (ProhibitMagicNumbers)
my $read_name_base = 'HS';
$read_name_base .= ( int rand 50 ) + 1;        # Instrument name
$read_name_base .= q{_};
$read_name_base .= ( int rand 20_000 ) + 1;    # Run
$read_name_base .= q{:};
$read_name_base .= ( int rand 8 ) + 1;         # Flowcell lane
$read_name_base .= q{:};
## use critic

my %tag_count;

## no critic (RequireBriefOpen)
open my $fh1, '>', $output_prefix . '_1.fastq';
open my $fh2, '>', $output_prefix . '_2.fastq';
## use critic
foreach ( 1 .. $read_pair_count ) {
    my $read_name = get_read_name($read_name_base);

    my $tag = @read_tags[ int rand scalar @read_tags ];    # Random tag

    # 20% of read 1s have no polyT
    my $has_polyt = int rand 5 ? 1 : 0;    ## no critic (ProhibitMagicNumbers)

    write_or_die( $fh1, q{@}, $read_name, '/1', "\n" );
    write_or_die( $fh1, get_read1_seq( $read_length, $tag, $has_polyt ), "\n" );
    write_or_die( $fh1, "+\n" );
    write_or_die( $fh1, q{~} x $read_length, "\n" );
    write_or_die( $fh2, q{@}, $read_name, '/2', "\n" );
    write_or_die( $fh2, get_read2_seq($read_length), "\n" );
    write_or_die( $fh2, "+\n" );
    write_or_die( $fh2, q{~} x $read_length, "\n" );

    if ( !$has_polyt ) {
        $tag = $read_tags[-1];             # Dummy tag
    }
    $tag_count{$tag}++;
}
close $fh1;
close $fh2;

# Display tag counts
foreach my $read_tag (@read_tags) {
    print $output_prefix, "\t", $read_tag, ":\t",
      ( $tag_count{$read_tag} || 0 ), "\n";
}

# Construct read name
sub get_read_name {
    my ( $read_name, ) = @_;

    ## no critic (ProhibitMagicNumbers)
    $read_name .= ( int rand 3_000 ) + 1;      # Tile number
    $read_name .= q{:};
    $read_name .= ( int rand 20_000 ) + 1;     # Cluster x coordinate
    $read_name .= q{:};
    $read_name .= ( int rand 200_000 ) + 1;    # Cluster y coordinate
    ## use critic

    return $read_name;
}

# Get read 1 sequence (just random but with tag)
sub get_read1_seq {
    my ( $read_len, $tag, $has_polyt ) = @_;

    my $is_dummy_tag = $tag =~ m/X/xms ? 1 : 0;    # If tag is X then no tag

    # Replace IUPAC codes in tag with random bases
    $tag =~ s/ N / qw( A G C T )[ int rand 4 ] /xmsge;
    $tag =~ s/ B / qw( G C T   )[ int rand 3 ] /xmsge;
    $tag =~ s/ D / qw( A G T   )[ int rand 3 ] /xmsge;
    $tag =~ s/ H / qw( A C T   )[ int rand 3 ] /xmsge;
    $tag =~ s/ V / qw( A G C   )[ int rand 3 ] /xmsge;
    $tag =~ s/ R / qw( A G     )[ int rand 2 ] /xmsge;
    $tag =~ s/ Y / qw( C T     )[ int rand 2 ] /xmsge;
    $tag =~ s/ K / qw( G T     )[ int rand 2 ] /xmsge;
    $tag =~ s/ M / qw( A C     )[ int rand 2 ] /xmsge;
    $tag =~ s/ S / qw( G C     )[ int rand 2 ] /xmsge;
    $tag =~ s/ W / qw( A T     )[ int rand 2 ] /xmsge;
    $tag =~ s/ X / qw( A G C T )[ int rand 4 ] /xmsge;    # Not actually IUPAC

    # Make the last two bases be Ns so should never match a real tag
    if ( $is_dummy_tag && $has_polyt ) {                  # No need if not polyT
        substr $tag, -2, 2, 'NN';    ## no critic (ProhibitMagicNumbers)
    }

    # 20% of reads have a single mismatch somewhere in the tag
    if ( int rand 5 ) {              ## no critic (ProhibitMagicNumbers)
        my $mismatch_base = int rand length $tag;
        my $base = substr $tag, $mismatch_base, 1;
        $base =~ tr/AGCT/TCGA/;
        substr $tag, $mismatch_base, 1, $base;
    }

    # Read begins with tag then polyT (or add polyA if doesn't have polyT)
    my $seq = $tag;
    $seq .= $has_polyt ? q{T} x $polyt_length : q{A} x $polyt_length;
    $read_len -= length $seq;

    # Rest of read is random
    ## no critic (ProhibitMagicNumbers)
    $seq .= join q{}, map { qw( A G C T ) [ int rand 4 ] } 1 .. $read_len;
    ## use critic

    return $seq;
}

# Get read 2 sequence (just random)
sub get_read2_seq {
    my ($read_len) = @_;

    ## no critic (ProhibitMagicNumbers)
    return join q{}, map { qw( A G C T ) [ int rand 4 ] } 1 .. $read_len;
    ## use critic
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'seed=i'            => \$seed,
        'output_prefix=s'   => \$output_prefix,
        'read_pair_count=i' => \$read_pair_count,
        'read_tags=s@{1,}'  => \@read_tags,
        'read_length=i'     => \$read_length,
        'polyt_length=i'    => \$polyt_length,
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
    if ( !$output_prefix ) {
        pod2usage("--output_prefix must be specified\n");
    }
    if ( !$read_pair_count ) {
        pod2usage("--read_pair_count must be a positive integer\n");
    }
    if ( !$read_length ) {
        pod2usage("--read_length must be a positive integer\n");
    }
    if ( !@read_tags ) {
        pod2usage("--read_tags must be specified\n");
    }

    return;
}

=head1 USAGE

    make_test_fastq.pl
        [--seed seed]
        [--output_prefix prefix]
        [--read_pair_count int]
        [--read_tags tags...]
        [--read_length int]
        [--polyt_length int]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--seed INT>

Random seed (to get reproducible chromosome lengths).

=item B<--output_prefix FILE>

Prefix for output FASTQ files.

=item B<--read_pair_count INT>

Number of read pairs aligned to each seq region (defaults to 100).

=item B<--read_tags TAGS>

Read tags.

=item B<--read_length INT>

Length of reads (defaults to 75 bp).

=item B<--polyt_length INT>

Length of polyT in read 1.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
