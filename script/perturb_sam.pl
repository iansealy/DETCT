#!/usr/bin/env perl

# PODNAME: perturb_sam.pl
# ABSTRACT: Randomly perturb SAM file by moving reads

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-05
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
use DETCT::Misc qw( print_or_die printf_or_die );

=head1 DESCRIPTION

This script reads a SAM file from STDIN and outputs the same SAM file but with
random perturbations to the position of each read. The resulting SAM file will
not be valid because PNEXT is not updated.

=head1 EXAMPLES

    # Perturb BAM file by maximum of 200 bp for each read
    samtools view -h test1.bam \
        | perl -Ilib script/perturb_sam.pl \
        | samtools view -bS - | samtools sort - test2

    # Reproducibly perturb BAM file by maximum of 50 bp for each read
    samtools view -h test1.bam \
        | perl -Ilib script/perturb_sam.pl --seed 1 --max_move 50 \
        | samtools view -bS - | samtools sort - test2

=cut

# Default options
## no critic (ProhibitMagicNumbers)
my $seed;
my $max_move = 200;
my ( $help, $man );
## use critic

# Get and check command line options
get_and_check_options();

# Ensure reproducible perturbations if seed set
if ( defined $seed ) {
    srand $seed;
}

# Construct command line
my @cl = ('perturb_sam.pl');
if ($seed) {
    push @cl, '--seed', $seed;
}
push @cl, '--max_move', $max_move;
my $cl = join q{ }, @cl;

# Read SAM from STDIN and output perturbed SAM to STDOUT
my %length_of;
my $seen_pg    = 0;
my $printed_pg = 0;
while ( my $line = <> ) {
    chomp $line;

    # Add @PG header line after existing ones
    if ( $line =~ m/\@PG/xms ) {
        $seen_pg = 1;
    }
    if ( !$printed_pg && $seen_pg && $line !~ m/\@PG/xms ) {
        ## no critic (RequireInterpolationOfMetachars)
        printf_or_die( "%s\n", join "\t", '@PG', 'ID:2', 'PN:perturb_sam.pl',
            "CL:$cl" );
        ## use critic
        $printed_pg = 1;
    }

    if ( $line =~ m/\@SQ/xms ) {

        # @SQ header line
        my ( undef, $seq_region, $length ) = split /\t/xms, $line;
        ( undef, $seq_region ) = split /:/xms, $seq_region;
        ( undef, $length )     = split /:/xms, $length;
        $length_of{$seq_region} = $length;
        print_or_die( $line, "\n" );
    }
    elsif ( $line =~ m/\@/xms ) {

        # Other header line
        print_or_die( $line, "\n" );
    }
    else {
        my @fields     = split /\t/xms, $line;
        my $seq_region = $fields[2];
        my $pos = $fields[3];    ## no critic (ProhibitMagicNumbers)
        my $move =
          ( int rand( $max_move * 2 + 1 ) ) - $max_move;    # e.g. -200 to 200
        $pos += $move;
        if ( $pos < 1 ) {
            $pos = 1;
        }
        elsif ( $pos > $length_of{$seq_region} ) {
            $pos = $seq_region;
        }
        $fields[3] = $pos;    ## no critic (ProhibitMagicNumbers)
        printf_or_die( "%s\n", join "\t", @fields );
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'seed=i'     => \$seed,
        'max_move=i' => \$max_move,
        'help'       => \$help,
        'man'        => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

=head1 USAGE

    perturb_sam.pl
        [--seed seed]
        [--max_move int]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--seed INT>

Random seed (to get reproducible perturbations.

=item B<--max_move INT>

Maximum number of bases each read can be moved (default to 200).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
