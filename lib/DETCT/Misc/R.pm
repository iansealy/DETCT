## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::R;
## use critic

# ABSTRACT: Miscellaneous functions for running R

## Author         : is1
## Maintainer     : is1
## Created        : 2012-11-21
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use English qw( -no_match_vars );
use POSIX qw( WIFEXITED);
use File::Slurp;
use File::Spec;
use File::Path qw( make_path );
use Sort::Naturally;
use List::Util qw( sum );
use List::MoreUtils qw( uniq );

use base qw( Exporter );
our @EXPORT_OK = qw(
  run_deseq
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func run_deseq

  Usage       : my $regions_ref = DETCT::Misc::R::run_deseq( {
                    dir          => '.',
                    regions      => $regions_hash_ref,
                    samples      => $samples_ary_ref,
                    r_binary     => 'R',
                    deseq_script => 'bin/run_deseq.R',
                } );
  Purpose     : Run DESeq
  Returns     : Arrayref [
                    Arrayref [
                        String (region sequence name),
                        Int (region start),
                        Int (region end),
                        Int (region maximum read count),
                        Float (region log probability sum),
                        String (3' end sequence name) or undef,
                        Int (3' end position) or undef,
                        Int (3' end strand) or undef,
                        Int (3' end read count) or undef,
                        Arrayref [
                            Int (count)
                            ...
                        ],
                        Arrayref [
                            Int (normalised count)
                            ...
                        ],
                        Int (p value) or undef,
                        Int (adjusted p value) or undef,
                        Arrayref [
                            Int (condition fold change) or undef,
                            Int (log2 condition fold change) or undef,
                        ],
                        Arrayref [
                            Arrayref [
                                Int (group fold change) or undef,
                                Int (log2 group fold change) or undef,
                            ],
                            ... (groups)
                        ]
                    ],
                    ... (regions)
                ]
  Parameters  : Hashref {
                    dir          => String (the working directory),
                    regions      => Hashref (of arrayrefs of regions),
                    samples      => Arrayref (of samples)
                    r_binary     => String (the R binary),
                    deseq_script => String (the DESeq script),
                }
  Throws      : If directory is missing
                If regions are missing
                If samples are missing
                If R binary is missing
                If DESeq script is missing
                If command line can't be run
  Comments    : None

=cut

sub run_deseq {
    my ($arg_ref) = @_;

    confess 'No directory specified'    if !defined $arg_ref->{dir};
    confess 'No regions specified'      if !defined $arg_ref->{regions};
    confess 'No samples specified'      if !defined $arg_ref->{samples};
    confess 'No R binary specified'     if !defined $arg_ref->{r_binary};
    confess 'No DESeq script specified' if !defined $arg_ref->{deseq_script};

    # Get conditions and groups
    my @samples = @{ $arg_ref->{samples} };
    my @conditions = uniq( nsort( map { $_->condition } @samples ) );
    my @groups = uniq( nsort( map { $_->group } grep { $_->group } @samples ) );
    @groups = grep { defined $_ } @groups;

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Write regions to input file
    my $input_file = File::Spec->catfile( $arg_ref->{dir}, 'input.txt' );
    my @sample_names = map { $_->name } @samples;
    ## no critic (RequireBriefOpen)
    open my $input_fh, '>', $input_file;
    print {$input_fh} ( join "\t", q{}, @sample_names ), "\n";
    foreach my $seq_name ( nsort( keys %{ $arg_ref->{regions} } ) ) {
        foreach my $region ( @{ $arg_ref->{regions}->{$seq_name} } ) {
            my $counts = $region->[-1];
            my $start  = $region->[0];
            my $end    = $region->[1];
            ## no critic (ProhibitMagicNumbers)
            my $strand = $region->[6];
            ## use critic
            my $region_text = join q{:}, $seq_name, $start, $end, $strand;
            print {$input_fh} ( join "\t", $region_text, @{$counts} ), "\n";
        }
    }
    close $input_fh;
    ## use critic

    # Write samples to input file
    my $samples_file = File::Spec->catfile( $arg_ref->{dir}, 'samples.txt' );
    my $last_col_to_print = @groups > 1 ? 2 : 1;
    my $condition_prefix  = condition_prefix( \@samples );
    my $group_prefix      = group_prefix( \@samples );
    my $samples_text;
    foreach my $sample (@samples) {
        my @row = (
            $sample->name,
            $condition_prefix . $sample->condition,
            $group_prefix . ( $sample->group || q{} ),
        )[ 0 .. $last_col_to_print ];
        $samples_text .= ( join "\t", @row ) . "\n";
    }
    open my $samples_fh, '>', $samples_file;
    my @header = ( q{}, 'condition', 'group' )[ 0 .. $last_col_to_print ];
    print {$samples_fh} ( join "\t", @header ), "\n";
    print {$samples_fh} $samples_text;
    close $samples_fh;

    my $output_file = File::Spec->catfile( $arg_ref->{dir}, 'output.txt' );
    my $size_factors_file =
      File::Spec->catfile( $arg_ref->{dir}, 'size_factors.txt' );
    my $power_file  = File::Spec->catfile( $arg_ref->{dir}, 'power.txt' );
    my $qc_pdf_file = File::Spec->catfile( $arg_ref->{dir}, 'qc.pdf' );
    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'deseq.o' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'deseq.e' );

    my $cmd = join q{ }, $arg_ref->{r_binary}, '--slave', '--args',
      $input_file, $samples_file, $output_file, $size_factors_file,
      $qc_pdf_file, $power_file, '<', $arg_ref->{deseq_script};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    # Get size factors for each sample
    my @size_factors = read_file($size_factors_file);
    chomp @size_factors;

    # Get output
    my %pval_for;
    my %padj_for;
    foreach my $line ( read_file($output_file) ) {
        chomp $line;
        my ( $region_text, $pval, $padj ) = split /\t/xms, $line;
        $pval_for{$region_text} = $pval;
        $padj_for{$region_text} = $padj;
    }

    # Reformat output into array of arrayrefs
    my @output;
    foreach my $seq_name ( nsort( keys %{ $arg_ref->{regions} } ) ) {
        foreach my $region ( @{ $arg_ref->{regions}->{$seq_name} } ) {
            my $counts = $region->[-1];
            my $start  = $region->[0];
            my $end    = $region->[1];
            ## no critic (ProhibitMagicNumbers)
            my $strand = $region->[6];
            ## use critic
            my $region_text = join q{:}, $seq_name, $start, $end, $strand;

            # Add sequence name to region
            unshift @{$region}, $seq_name;

            # Normalise counts and store for fold change calculation
            my @normalised_counts;
            my %counts_for_condition;
            my %counts_for_group_condition;
            my $sample_index = 0;
            foreach my $sample (@samples) {
                my $normalised_count =
                  $counts->[$sample_index] / $size_factors[$sample_index];
                push @normalised_counts, $normalised_count;
                if ( scalar @conditions == 2 ) {
                    push @{ $counts_for_condition{ $sample->condition } },
                      $normalised_count;
                }
                if ( scalar @conditions == 2 && scalar @groups > 1 ) {
                    push @{ $counts_for_group_condition{ $sample->group }
                          { $sample->condition } }, $normalised_count;
                }
                $sample_index++;
            }
            push @{$region}, \@normalised_counts;

            # Add p value and adjusted p value
            push @{$region}, $pval_for{$region_text}, $padj_for{$region_text};

            # Add condition fold change
            push @{$region},
              calc_condition_fold_change( \@conditions,
                \%counts_for_condition );

            #
            push @{$region},
              calc_group_fold_changes( \@conditions, \@groups,
                \%counts_for_group_condition );

            push @output, $region;
        }
    }

    return \@output;
}

=func condition_prefix

  Usage       : condition_prefix( \@samples );
  Purpose     : Return a prefix if all conditions are numeric
  Returns     : String (the prefix)
  Parameters  : Arrayref (of samples)
  Throws      : No exceptions
  Comments    : Numeric conditions won't be treated as R factors

=cut

sub condition_prefix {
    my ($samples) = @_;

    return prefix( 'condition', $samples );
}

=func group_prefix

  Usage       : group_prefix( \@samples );
  Purpose     : Return a prefix if all groups are numeric
  Returns     : String (the prefix)
  Parameters  : Arrayref (of samples)
  Throws      : No exceptions
  Comments    : Numeric groups won't be treated as R factors

=cut

sub group_prefix {
    my ($samples) = @_;

    return prefix( 'group', $samples );
}

=func prefix

  Usage       : prefix( 'condition', \@samples );
  Purpose     : Return a prefix if all of chosen attribute are numeric
  Returns     : String (the prefix)
  Parameters  : String (the attribute)
                Arrayref (of samples)
  Throws      : No exceptions
  Comments    : Numeric attributes won't be treated as R factors

=cut

sub prefix {
    my ( $attribute, $samples ) = @_;

    foreach my $sample ( @{$samples} ) {
        if ( $sample->$attribute && $sample->$attribute !~ m/\A [\d.]+ \z/xms )
        {
            return q{};
        }
    }

    return $attribute;
}

=func calc_condition_fold_change

  Usage       : calc_condition_fold_change( \@conditions, \%counts );
  Purpose     : Calculate the condition fold change if two conditions
  Returns     : Arrayref [
                    Int (condition fold change) or undef,
                    Int (log2 condition fold change) or undef,
                ],
  Parameters  : Arrayref (of conditions)
                Hashref (of counts)
  Throws      : No exceptions
  Comments    : None

=cut

sub calc_condition_fold_change {
    my ( $conditions, $counts_for_condition ) = @_;

    my $fold_change;
    my $log2_fold_change;

    if ( scalar @{$conditions} == 2 ) {
        ( $fold_change, $log2_fold_change ) = calc_fold_change(
            $counts_for_condition->{ $conditions->[0] },
            $counts_for_condition->{ $conditions->[1] }
        );
    }

    return [ $fold_change, $log2_fold_change ];
}

=func calc_group_fold_changes

  Usage       : calc_group_fold_changes( \@conditions, \@groups, \%counts );
  Purpose     : Calculate fold change for each group if two conditions
  Returns     : Arrayref [
                    Arrayref [
                        Int (group fold change) or undef,
                        Int (log2 group fold change) or undef,
                    ],
                    ... (groups)
                ]
  Parameters  : Arrayref (of conditions)
                Arrayref (of groups)
                Hashref (of counts)
  Throws      : No exceptions
  Comments    : None

=cut

sub calc_group_fold_changes {
    my ( $conditions, $groups, $counts_for_group_condition ) = @_;

    my @group_fold_changes;

    if ( scalar @{$conditions} == 2 && scalar @{$groups} > 1 ) {
        foreach my $group ( @{$groups} ) {
            my ( $group_fold_change, $group_log2_fold_change ) =
              calc_fold_change(
                $counts_for_group_condition->{$group}{ $conditions->[0] },
                $counts_for_group_condition->{$group}{ $conditions->[1] }
              );
            push @group_fold_changes,
              [ $group_fold_change, $group_log2_fold_change ];
        }
    }

    return \@group_fold_changes;
}

=func calc_fold_change

  Usage       : ($fold_change, $log2_fold_change)
                    = calc_fold_change(\@array1, \@array2);
  Purpose     : Calculate the fold change in mean value of two arrays
  Returns     : Int (fold change)
                Int (log2 fold change)
  Parameters  : Arrayref
                Arrayref
  Throws      : No exceptions
  Comments    : None

=cut

sub calc_fold_change {
    my ( $array1_ref, $array2_ref ) = @_;

    my $fold_change;
    my $log2_fold_change;
    my $mean1 = sum( @{$array1_ref} ) / scalar @{$array1_ref};
    my $mean2 = sum( @{$array2_ref} ) / scalar @{$array2_ref};
    if ( $mean1 && $mean2 ) {
        $fold_change      = $mean1 / $mean2;             # e.g. mutant / sibling
        $log2_fold_change = log($fold_change) / log 2;
    }

    return $fold_change, $log2_fold_change;
}

1;
