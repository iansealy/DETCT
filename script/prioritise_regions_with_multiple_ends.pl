#!/usr/bin/env perl

# PODNAME: prioritise_regions_with_multiple_ends.pl
# ABSTRACT: Prioritise regions with multiple 3' ends

## Author         : is1
## Maintainer     : is1
## Created        : 2018-03-08
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
use Sort::Naturally;
use List::Util qw( max sum );
use List::MoreUtils qw( uniq );
use DETCT::Misc qw( printf_or_die );
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Default options
my @files;
my $count_threshold = 0;
my $sig_level       = 0.05;    ## no critic (ProhibitMagicNumbers)
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Get regions with mean counts above threshold
my $is_region = get_regions_above_threshold(@files);

# Get genes with multiple 3' ends arbitrarily from first file
my $is_multiple_end_gene =
  get_genes_with_multiple_ends( $files[0], $is_region );

# Get regions for each gene with multiple 3' ends
my ( $regions_for, $header_fields ) =
  get_regions_for_genes( $files[0], $is_region, $is_multiple_end_gene );

# Get significance matrix from all files
my $significance_matrix =
  get_significance_matrix( \@files, $is_region, $is_multiple_end_gene );

# Score each gene according to how many comparisons are different across ends
my $score_for = score_genes( $significance_matrix, $regions_for );

# Output matrix ordered by score
my @comparisons = remove_common_prefix_suffix(@files);
push @{$header_fields}, @comparisons;
output_matrix( $significance_matrix, $regions_for, $score_for, $header_fields );

sub get_regions_above_threshold {
    my (@files) = @_;    ## no critic (ProhibitReusedNames)

    my %is_region;
    my %counts_for;
    foreach my $file (@files) {
        my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
        if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
            confess sprintf '%s is not .csv or .tsv file', $file;
        }

        open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
        my $header = <$fh>;
        my @header_fields =
          DETCT::Misc::Output::parse_line( $header, $extension );
        my %sample_for;
        my $i = 0;
        foreach my $field (@header_fields) {
            if ( $field =~ m/\A ([\w.-]+) \s normalised \s count \z/xms ) {
                my $sample = $1;
                $sample_for{$i} = $sample;
            }
            $i++;
        }
        while ( my $line = <$fh> ) {
            my @fields = DETCT::Misc::Output::parse_line( $line, $extension );
            ## no critic (ProhibitMagicNumbers)
            # Get region ID by joining chr, start, end and strand
            my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
            ## use critic
            $i = 0;
            foreach my $field (@fields) {
                if ( exists $sample_for{$i} ) {
                    push @{ $counts_for{$region}{ $sample_for{$i} } }, $field;
                }
                $i++;
            }
        }
        close $fh;
    }

    foreach my $region ( keys %counts_for ) {
        my @mean_counts;
        foreach my $sample ( keys %{ $counts_for{$region} } ) {
            my $mean = sum( @{ $counts_for{$region}{$sample} } ) /
              scalar @{ $counts_for{$region}{$sample} };
            push @mean_counts, $mean;
        }
        my $mean = sum(@mean_counts) / scalar @mean_counts;
        next if $mean < $count_threshold;
        $is_region{$region} = 1;
    }

    return \%is_region;
}

sub get_genes_with_multiple_ends {
    my ( $file, $is_region ) = @_;    ## no critic (ProhibitReusedNames)

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %region_count_for;
    my %end_count_for;

    open my $fh, '<', $file;          ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        ## no critic (ProhibitMagicNumbers)
        # Get region ID by joining chr, start, end and strand
        my $region              = join q{:}, @fields[ 0, 1, 2, 4 ];
        my $genes               = $fields[9];
        my @three_prime_end_pos = split /,/xms, $fields[3];
        ## use critic
        next if !$genes;
        next if !exists $is_region->{$region};
        foreach my $gene ( split /,/xms, $genes ) {
            $region_count_for{$gene}++;
            foreach my $three_prime_end_pos (@three_prime_end_pos) {

                # Get 3' end ID by join chr and end position
                my $three_prime_end = join q{:}, $fields[0],
                  $three_prime_end_pos;
                $end_count_for{$gene}{$three_prime_end}{$region} = 1;
            }
        }
    }
    close $fh;

    my %is_multiple_end_gene;
  GENE: foreach my $gene ( keys %region_count_for ) {
        next if $region_count_for{$gene} == 1;
        foreach my $end ( keys %{ $end_count_for{$gene} } ) {
            next GENE
              if scalar
              keys %{ $end_count_for{$gene}{$end} } == $region_count_for{$gene};
        }
        $is_multiple_end_gene{$gene} = 1;
    }

    return \%is_multiple_end_gene;
}

sub get_regions_for_genes {
    ## no critic (ProhibitReusedNames)
    my ( $file, $is_region, $is_multiple_end_gene ) = @_;
    ## use critic

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %regions_for;

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header        = <$fh>;
    my @header_fields = DETCT::Misc::Output::parse_line( $header, $extension );
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );
        @fields = map { defined $_ && length $_ > 0 ? $_ : q{-} } @fields;

        ## no critic (ProhibitMagicNumbers)
        # Get region ID by joining chr, start, end and strand
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        my $genes  = $fields[9];
        ## use critic
        next if !$genes;
        next if !exists $is_region->{$region};
        foreach my $gene ( split /,/xms, $genes ) {
            next if !exists $is_multiple_end_gene->{$gene};
            ## no critic (ProhibitMagicNumbers)
            $regions_for{$gene}{$region} = [ @fields[ 0 .. 14 ] ];
            ## use critic
        }
    }
    close $fh;

    ## no critic (ProhibitMagicNumbers)
    return \%regions_for, [ @header_fields[ 0 .. 14 ] ];
    ## use critic
}

sub get_significance_matrix {
    ## no critic (ProhibitReusedNames)
    my ( $files, $is_region, $is_multiple_end_gene ) = @_;
    ## use critic

    my %significance_matrix;

    my $i = 0;
    foreach my $file ( @{$files} ) {
        my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
        if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
            confess sprintf '%s is not .csv or .tsv file', $file;
        }
        open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
        my $header = <$fh>;
        while ( my $line = <$fh> ) {
            my @fields = DETCT::Misc::Output::parse_line( $line, $extension );
            ## no critic (ProhibitMagicNumbers)
            my $genes = $fields[9];

            # Get region ID by joining chr, start, end and strand
            my $region          = join q{:}, @fields[ 0, 1, 2, 4 ];
            my $adjusted_pvalue = $fields[7];
            ## use critic

            # Skip if no gene ID or below threshold or not multiple 3' ends
            next if !$genes;
            next if !exists $is_region->{$region};
            foreach my $gene ( split /,/xms, $genes ) {
                next if !exists $is_multiple_end_gene->{$gene};
                push @{ $significance_matrix{$region} },
                    $adjusted_pvalue eq 'NA'      ? 'NA'
                  : $adjusted_pvalue < $sig_level ? 1
                  :                                 0;
                last;
            }
        }
        close $fh;
        $i++;
    }

    return \%significance_matrix;
}

sub score_genes {
    ## no critic (ProhibitReusedNames)
    my ( $significance_matrix, $regions_for ) = @_;
    ## use critic

    my %score_for;

    foreach my $gene ( keys %{$regions_for} ) {
        my @regions = keys %{ $regions_for->{$gene} };
        confess "Only one region for $gene" if scalar @regions == 1;
        my @scores = (0);
        foreach my $region1_idx ( 0 .. scalar @regions - 2 ) {
            my $region1 = $regions[$region1_idx];
            my @sigs =
              grep { $_ ne 'NA' } @{ $significance_matrix->{$region1} };
            next if uniq(@sigs) == 1;
            foreach my $region2_idx ( $region1_idx + 1 .. scalar @regions - 1 )
            {
                my $region2 = $regions[$region2_idx];
                @sigs =
                  grep { $_ ne 'NA' } @{ $significance_matrix->{$region2} };
                next if uniq(@sigs) == 1;
                my $score = 0;
                foreach my $i (
                    0 .. scalar @{ $significance_matrix->{$region1} } - 1 )
                {
                    next if $significance_matrix->{$region1}->[$i] eq 'NA';
                    next if $significance_matrix->{$region2}->[$i] eq 'NA';
                    if ( $significance_matrix->{$region1}->[$i] !=
                        $significance_matrix->{$region2}->[$i] )
                    {
                        $score++;
                    }
                }
                push @scores, $score;
            }
        }
        $score_for{$gene} = max(@scores);
    }

    return \%score_for;
}

sub remove_common_prefix_suffix {
    my @strings = @_;

    foreach ( 0 .. 1 ) {
        @strings = map { scalar reverse } @strings;
        my $common_prefix_idx = 0;
      CHAR: while (1) {
            my $char = substr $strings[0], $common_prefix_idx, 1;
            foreach my $string (@strings) {
                last CHAR if substr( $string, $common_prefix_idx, 1 ) ne $char;
            }
            $common_prefix_idx++;
        }
        if ($common_prefix_idx) {
            @strings = map { substr $_, $common_prefix_idx } @strings;
        }
    }

    return @strings;
}

sub output_matrix {
    ## no critic (ProhibitReusedNames)
    my ( $significance_matrix, $regions_for, $score_for, $header_fields ) = @_;
    ## use critic

    printf_or_die( "%s\n", join "\t", @{$header_fields} );

    foreach my $gene (
        reverse sort { $score_for->{$a} <=> $score_for->{$b} || $a cmp $b }
        keys %{$score_for}
      )
    {
        foreach my $region ( nsort( keys %{ $regions_for->{$gene} } ) ) {
            printf_or_die(
                "%s\n", join "\t",
                @{ $regions_for->{$gene}->{$region} },
                @{ $significance_matrix->{$region} }
            );
        }
    }

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'files=s@{2,}'      => \@files,
        'count_threshold=i' => \$count_threshold,
        'sig_level=i'       => \$sig_level,
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
    if ( !@files ) {
        pod2usage("--files must be specified\n");
    }

    return;
}

=head1 USAGE

    prioritise_regions_with_multiple_ends.pl
        [--files file...]
        [--count_threshold int]
        [--sig_level float]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--files FILE...>

Differential expression pipeline output files (e.g. all.tsv).

=item B<--count_threshold INT>

Mean count threshold for removing regions (defaults to 0).

=item B<--sig_level FLOAT>

Adjusted p-value significance level (defaults to 0.05).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
