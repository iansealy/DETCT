## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::R;
## VERSION
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
use Path::Tiny;
use File::Spec;
use File::Path qw( make_path );
use Sort::Naturally;
use List::Util qw( sum );
use List::MoreUtils qw( uniq );
use DETCT::Misc qw( write_or_die );

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
                    analysis     => $analysis,
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
                        Int (3' end position) or arrayref or undef,
                        Int (3' end strand) or undef,
                        Int (3' end read count) or arrayref or undef,
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
                    dir                  => String (the working directory),
                    regions              => Hashref (of arrayrefs of regions),
                    analysis             => DETCT::Analysis,
                    r_binary             => String (the R binary),
                    deseq_script         => String (the DESeq script),
                    filter_percentile    => Int (the filter percentile) or undef,
                    spike_prefix         => String (the spike prefix) or undef,
                    normalisation_method => String (the normalisation method),
                    deseq_model          => String (the DESeq model),
                }
  Throws      : If directory is missing
                If regions are missing
                If analysis is missing
                If R binary is missing
                If DESeq script is missing
                If command line can't be run
  Comments    : None

=cut

sub run_deseq {    ## no critic (ProhibitExcessComplexity)
    my ($arg_ref) = @_;

    confess 'No directory specified'    if !defined $arg_ref->{dir};
    confess 'No regions specified'      if !defined $arg_ref->{regions};
    confess 'No analysis specified'     if !defined $arg_ref->{analysis};
    confess 'No R binary specified'     if !defined $arg_ref->{r_binary};
    confess 'No DESeq script specified' if !defined $arg_ref->{deseq_script};

    # Normalise based on spikes?
    my $normalisation_method = $arg_ref->{normalisation_method};
    my $spike_prefix         = $arg_ref->{spike_prefix};
    confess 'No spike prefix specified'
      if $normalisation_method eq 'spike' && !$spike_prefix;

    # Get conditions and groups
    my @samples    = @{ $arg_ref->{analysis}->get_all_samples() };
    my @conditions = $arg_ref->{analysis}->list_all_conditions();
    my @groups     = $arg_ref->{analysis}->list_all_groups();

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Write regions to input file
    my $input_file = File::Spec->catfile( $arg_ref->{dir}, 'input.txt' );
    my $spike_input_file =
      File::Spec->catfile( $arg_ref->{dir}, 'spike_input.txt' );
    my @sample_names = map { $_->name } @samples;
    if ( $normalisation_method ne 'spike' ) {
        write_deseq_input( $input_file, \@sample_names, $arg_ref->{regions} );
    }
    else {
        write_deseq_input( $input_file, \@sample_names, $arg_ref->{regions},
            $spike_prefix, 0 );
        write_deseq_input( $spike_input_file, \@sample_names,
            $arg_ref->{regions}, $spike_prefix, 1 );
    }

    # Write samples to input file
    my $samples_file = File::Spec->catfile( $arg_ref->{dir}, 'samples.txt' );
    my $condition_prefix = condition_prefix( \@samples );
    my $group_prefix     = group_prefix( \@samples );
    my $samples_text;
    foreach my $sample (@samples) {
        my @row = (
            $sample->name,
            $condition_prefix . $sample->condition,
            map { $group_prefix . $_ } @{ $sample->groups },
        );
        $samples_text .= ( join "\t", @row ) . "\n";
    }
    my @header = ( q{}, 'condition' );
    my $num_groups = scalar @{ $samples[0]->groups };
    if ( $num_groups == 1 ) {
        push @header, 'group';
    }
    elsif ( $num_groups > 1 ) {
        push @header, map { 'group' . $_ } ( 1 .. $num_groups );
    }
    open my $samples_fh, '>', $samples_file;
    write_or_die( $samples_fh, ( join "\t", @header ), "\n" );
    write_or_die( $samples_fh, $samples_text );
    close $samples_fh;

    my $control_condition = q{-};
    if ( $arg_ref->{analysis}->control_condition ) {
        $control_condition =
          $condition_prefix . $arg_ref->{analysis}->control_condition;
    }
    elsif ( scalar @conditions == 2 ) {
        $control_condition = $condition_prefix . $conditions[1];
    }
    my $experimental_condition = q{-};
    if ( $arg_ref->{analysis}->experimental_condition ) {
        $experimental_condition =
          $condition_prefix . $arg_ref->{analysis}->experimental_condition;
    }
    elsif ( scalar @conditions == 2 ) {
        $experimental_condition = $condition_prefix . $conditions[0];
    }

    if (   $arg_ref->{analysis}->control_condition
        && $arg_ref->{analysis}->experimental_condition )
    {
        @conditions = ( $control_condition, $experimental_condition );
    }

    my $output_file = File::Spec->catfile( $arg_ref->{dir}, 'output.txt' );
    my $size_factors_file =
      File::Spec->catfile( $arg_ref->{dir}, 'size_factors.txt' );
    my $qc_pdf_file       = File::Spec->catfile( $arg_ref->{dir}, 'qc.pdf' );
    my $filter_percentile = $arg_ref->{filter_percentile};
    my $stdout_file       = File::Spec->catfile( $arg_ref->{dir}, 'deseq.o' );
    my $stderr_file       = File::Spec->catfile( $arg_ref->{dir}, 'deseq.e' );

    my $cmd = join q{ }, $arg_ref->{r_binary}, '--slave', '--args',
      $input_file, $samples_file, $output_file, $size_factors_file,
      $qc_pdf_file, $filter_percentile, $normalisation_method,
      $arg_ref->{deseq_model}, $spike_input_file, $control_condition,
      $experimental_condition, '<', $arg_ref->{deseq_script};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    # Get size factors for each sample
    my @size_factors = path($size_factors_file)->lines( { chomp => 1 } );

    # Get output
    my %pval_for;
    my %padj_for;
    foreach my $line ( path($output_file)->lines( { chomp => 1 } ) ) {
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
                    foreach my $group ( @{ $sample->groups } ) {
                        push @{ $counts_for_group_condition{$group}
                              { $sample->condition } }, $normalised_count;
                    }
                }
                $sample_index++;
            }
            push @{$region}, \@normalised_counts;

            # Add p value and adjusted p value
            push @{$region}, $pval_for{$region_text} || 'NA',
              $padj_for{$region_text} || 'NA';

            # Add condition fold change
            push @{$region},
              calc_condition_fold_change( \@conditions,
                \%counts_for_condition );

            # Add group fold changes
            push @{$region},
              calc_group_fold_changes( \@conditions, \@groups,
                \%counts_for_group_condition );

            push @output, $region;
        }
    }

    return \@output;
}

=func write_deseq_input

  Usage       : write_deseq_input( 'input.txt', $sample_names, $regions );
  Purpose     : Write DESeq input to a file
  Returns     : undef
  Parameters  : String (the input filename)
                Arrayref of strings (the sample names)
                Hashref (of arrayrefs of regions)
                String (the spike prefix) or undef
                Boolean (whether to match spike prefix) or undef
  Throws      : No exceptions
  Comments    : None

=cut

sub write_deseq_input {
    my ( $input_file, $sample_names, $regions, $spike_prefix, $match_prefix ) =
      @_;

    ## no critic (RequireBriefOpen)
    open my $input_fh, '>', $input_file;
    write_or_die( $input_fh, ( join "\t", q{}, @{$sample_names} ), "\n" );
    foreach my $seq_name ( nsort( keys %{$regions} ) ) {
        next
          if $spike_prefix
          && ( ( $match_prefix && $seq_name !~ m/\A $spike_prefix /xms )
            || ( !$match_prefix && $seq_name =~ m/\A $spike_prefix /xms ) );
        foreach my $region ( @{ $regions->{$seq_name} } ) {
            my $counts = $region->[-1];
            my $start  = $region->[0];
            my $end    = $region->[1];
            ## no critic (ProhibitMagicNumbers)
            my $strand = $region->[6];
            ## use critic
            my $region_text = join q{:}, $seq_name, $start, $end, $strand;
            write_or_die( $input_fh, ( join "\t", $region_text, @{$counts} ),
                "\n" );
        }
    }
    close $input_fh;
    ## use critic

    return;
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

    foreach my $sample ( @{$samples} ) {
        if ( $sample->condition !~ m/\A [\d.]+ \z/xms ) {
            return q{};
        }
    }

    return 'condition';
}

=func group_prefix

  Usage       : group_prefix( \@samples );
  Purpose     : Return a prefix if all of a group are numeric
  Returns     : String (the prefix)
  Parameters  : Arrayref (of samples)
  Throws      : No exceptions
  Comments    : Numeric groups won't be treated as R factors

=cut

sub group_prefix {
    my ($samples) = @_;

    my @not_numeric = (0) * scalar @{ $samples->[0]->groups };

    foreach my $sample ( @{$samples} ) {
        my $i = 0;
        foreach my $group ( @{ $sample->groups } ) {
            if ( $group !~ m/\A [\d.]+ \z/xms ) {
                $not_numeric[$i] = 1;
            }
            $i++;
        }
    }

    my $prefix = 'group';
    if ( ( scalar grep { $_ } @not_numeric ) == scalar @not_numeric ) {
        $prefix = q{};
    }

    return $prefix;
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

    if ( scalar @{$conditions} == 2 && scalar @{$groups} > 0 ) {
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

    # Can't calculate fold change with empty arrays
    if ( !defined $array1_ref || !defined $array2_ref ) {
        return ( undef, undef );
    }

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
