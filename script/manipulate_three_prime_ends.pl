#!/usr/bin/env perl

# PODNAME: manipulate_three_prime_ends.pl
# ABSTRACT: Manipulate 3' ends file

## Author         : is1
## Maintainer     : is1
## Created        : 2016-08-04
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
use Path::Tiny;
use DETCT::Misc qw( printf_or_die );
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Default options
my $all_file;
my $ends_file;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

output_ends( $ends_file, get_annotation($all_file) );

sub get_annotation {
    my ($file) = @_;

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %annotation_for;

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID from chr, start or end (depending on strand) and strand
        ## no critic (ProhibitMagicNumbers)
        my $strand       = $fields[4];
        my $start_or_end = $strand > 0 ? 1 : 2;
        my $region       = join q{:}, @fields[ 0, $start_or_end, 4 ];
        my @annotation   = @fields[ 9 .. 13 ];
        ## use critic
        $annotation_for{$region} = \@annotation;
    }
    close $fh;

    return \%annotation_for;
}

sub output_ends {
    my ( $file, $annotation_for ) = @_;

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my @header = (
        '#Chr',
        'Region start',
        'Region end',
        q{3' end strand},
        q{3' end position},
        q{3' end read count},
        'polyA?',
        '14 bp upstream',
        '14 bp downstream',
        'Distance Hexamer Upstream (up to 40 bp)',
        'Hexamer',
        'Transposon distance',
        'Transposon position',
        'Continuous RNA-Seq transcripts',
        'Ensembl Gene ID',
        'Gene type',
        'Ensembl Transcript ID',
        'Transcript type',
        'Gene name',
        'Region aligment count',
        'Region alignments',
    );
    printf_or_die( "%s\n", join "\t", @header );

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID from chr, start or end (depending on strand) and strand
        ## no critic (ProhibitMagicNumbers)
        my $strand            = $fields[3];
        my $start_or_end      = $strand > 0 ? 1 : 2;
        my @region            = @fields[ 0 .. 3 ];
        my $region_alignments = $fields[4];
        my $ends              = $fields[5];
        my $annotation_region = join q{:}, @region[ 0, $start_or_end, 3 ];
        my @annotation        = @{ $annotation_for->{$annotation_region} };
        my $region            = join q{:}, @region[ 0 .. 3 ];
        ## use critic
        my ( $alignment_count, $other_alignments ) =
          remove_region_alignment( $region, $region_alignments );

        my @ends;
        if ($ends) {
            foreach my $end ( split /,/xms, $ends ) {
                my @end_data = split /[:\/]/xms, $end;
                splice @end_data, 1, 1;    # Remove redundant strand
                if ( scalar @end_data == 9 )
                {                          ## no critic (ProhibitMagicNumbers)
                    push @end_data, q{};    # Add missing last field
                }
                push @ends, \@end_data;
            }
            @ends = sort { $a->[0] <=> $b->[0] } @ends;
            if ( $strand < 0 ) {
                @ends = reverse @ends;
            }
        }
        else {
            @ends = ( [ q{} x 10 ] );       ## no critic (ProhibitMagicNumbers)
        }

        foreach my $end (@ends) {
            my @output = @region;
            push @output, @{$end};
            push @output, @annotation;
            push @output, $alignment_count, $other_alignments;
            @output = map { defined $_ && length $_ ? $_ : q{-} } @output;
            printf_or_die( "%s\n", join "\t", @output );
        }
    }
    close $fh;

    return;
}

# Exclude region aligning to itself
sub remove_region_alignment {
    my ( $region, $alignments ) = @_;

    return 0, undef if !$alignments;    # Region not even aligned to itself

    my ( $region_chr, $region_start, $region_end, $region_strand ) =
      split /:/xms, $region;

    my @other_alignments;
    my @alignments = split /,/xms, $alignments;
    foreach my $alignment (@alignments) {
        my ( $chr, $coords, $strand ) = split /:/xms, $alignment;
        my ( $start, $end ) = split /-/xms, $coords;
        $strand = $strand eq q{+} ? 1 : -1;  ## no critic (ProhibitMagicNumbers)
        next
          if $chr eq $region_chr
          && $start == $region_start
          && $end == $region_end
          && $strand == $region_strand;
        push @other_alignments, $alignment;
    }

    my $other_alignments = join q{,}, @other_alignments;

    return scalar @alignments, $other_alignments;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'all_file=s'  => \$all_file,
        'ends_file=s' => \$ends_file,
        'help'        => \$help,
        'man'         => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$all_file ) {
        pod2usage("--all_file must be specified\n");
    }
    if ( !$ends_file ) {
        pod2usage("--ends_file must be specified\n");
    }

    return;
}

=head1 USAGE

    manipulate_three_prime_ends.pl
        [--all_file file]
        [--ends_file file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--all_file FILE>

DETCT output file (containing all regions).

=item B<--ends_file FILE>

DETCT 3' ends file.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
