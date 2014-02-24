## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::PeakHMM;
## use critic

# ABSTRACT: Miscellaneous functions for running peaks HMM

## Author         : is1
## Maintainer     : is1
## Created        : 2012-10-29
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
use POSIX qw( WIFEXITED ceil );
use File::Slurp;
use File::Spec;
use File::Path qw( make_path );
use Memoize qw( memoize flush_cache );
use Cwd;
use DETCT::Misc qw( write_or_die );

use base qw( Exporter );
our @EXPORT_OK = qw(
  merge_read_peaks
  summarise_read_peaks
  run_peak_hmm
  join_hmm_bins
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func merge_read_peaks

  Usage       : my $peaks_ref = DETCT::Misc::PeakHMM::merge_read_peaks( {
                    peak_buffer_width => 100,
                    seq_name          => '1',
                    peaks             => $peaks_ary_ref,
                } );
  Purpose     : Merge read peaks (overlapping reads)
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
                        Arrayref [
                            Int (peak start),
                            Int (peak end),
                            Int (peak read count)
                        ],
                        ... (peaks)
                    ]
                }
  Parameters  : Hashref {
                    peak_buffer_width => Int (the peak buffer size),
                    seq_name          => String (the sequence name),
                    peaks             => Arrayref (of peaks),
                }
  Throws      : If peak buffer width is missing
                If sequence name is missing
                If peaks are missing
  Comments    : None

=cut

sub merge_read_peaks {
    my ($arg_ref) = @_;

    confess 'No peak buffer width specified'
      if !defined $arg_ref->{peak_buffer_width};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No peaks specified'         if !defined $arg_ref->{peaks};

    my $peaks_ref = $arg_ref->{peaks};

    # Sort peaks by start then end
    @{$peaks_ref} =
      sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$peaks_ref};

    # Peak variables
    my @merged_peaks;
    my $current_merged_peak_read_count;
    my $current_merged_peak_start;
    my $current_merged_peak_end;

    # Merge peaks
    foreach my $peak ( @{$peaks_ref} ) {
        my ( $peak_start, $peak_end, $peak_read_count ) = @{$peak};

        # We're starting the first merged peak
        if ( !defined $current_merged_peak_start ) {
            $current_merged_peak_start      = $peak_start;
            $current_merged_peak_end        = $peak_end;
            $current_merged_peak_read_count = $peak_read_count;
            next;
        }

        # Extend or finish current merged peak?
        if ( $peak_start - $current_merged_peak_end <
            $arg_ref->{peak_buffer_width} )
        {
            # Extend current merged peak
            $current_merged_peak_end = $peak_end;
            $current_merged_peak_read_count += $peak_read_count;
        }
        else {
            # Finish current merged peak
            push @merged_peaks,
              [
                $current_merged_peak_start, $current_merged_peak_end,
                $current_merged_peak_read_count
              ];

            # Start new merged peak
            $current_merged_peak_start      = $peak_start;
            $current_merged_peak_end        = $peak_end;
            $current_merged_peak_read_count = $peak_read_count;
        }
    }

    # Finish last merged peak
    if ($current_merged_peak_read_count) {
        push @merged_peaks,
          [
            $current_merged_peak_start, $current_merged_peak_end,
            $current_merged_peak_read_count
          ];
    }

    return { $arg_ref->{seq_name} => \@merged_peaks };
}

=func summarise_read_peaks

  Usage       : my $summary_ref = DETCT::Misc::PeakHMM::summarise_read_peaks( {
                    peak_buffer_width => 100,
                    hmm_sig_level     => 0.001,
                    total_bp          => 10_000_000,
                    read_length       => 54,
                    peaks             => $peaks_ary_ref,
                } );
  Purpose     : Summarise read peak distribution for HMM
  Returns     : Hashref {
                    total_read_count_per_mb     => Float,
                    total_sig_read_count_per_mb => Float,
                    total_sig_peak_width_in_mb  => Float,
                    median_sig_peak_width       => Int,
                    total_sig_peaks             => Int,
                    peak_buffer_width           => Int,
                    read_threshold              => Int,
                }
  Parameters  : Hashref {
                    peak_buffer_width => Int (the peak buffer size),
                    hmm_sig_level     => Float (the HMM significance level),
                    total_bp          => Int (the total bp),
                    read_length       => Int (the read length),
                    peaks             => Arrayref (of peaks),
                }
  Throws      : If peak buffer width is missing
                If HMM significance level is missing
                If total bp is missing
                if read length is missing
                If peaks are missing
  Comments    : Source of logic is summary.pl from
                http://www.sph.umich.edu/csg/qin/HPeak/

=cut

sub summarise_read_peaks {
    my ($arg_ref) = @_;

    confess 'No peak buffer width specified'
      if !defined $arg_ref->{peak_buffer_width};
    confess 'No HMM significance level specified'
      if !defined $arg_ref->{hmm_sig_level};
    confess 'No total bp specified'    if !defined $arg_ref->{total_bp};
    confess 'No read length specified' if !defined $arg_ref->{read_length};
    confess 'No peaks specified'       if !defined $arg_ref->{peaks};

    my $total_peaks = scalar @{ $arg_ref->{peaks} };

    if ( !$total_peaks ) {

        # No peaks so won't be running HMM
        return {};
    }

    # Get total read count
    my $total_read_count = 0;
    foreach my $peak ( @{ $arg_ref->{peaks} } ) {
        my ( $start, $end, $read_count ) = @{$peak};
        $total_read_count += $read_count;
    }

    # Get avg reads/bp
    my $avg_reads_per_bp = $total_read_count / $arg_ref->{total_bp};

    # Identify significant peaks
    memoize('_calc_log_sum');
    my @sig_peak_widths;
    my $total_sig_read_count = 0;
    my $total_sig_peak_width = 0;
    foreach my $peak ( @{ $arg_ref->{peaks} } ) {
        my ( $start, $end, $read_count ) = @{$peak};
        my $width         = $end - $start + 1;
        my $avg_reads     = $avg_reads_per_bp * $width;
        my $log_avg_reads = log $avg_reads;
        my $exp_avg_reads = exp $avg_reads;

        # Gather info for significant peaks
        my $sum = 1;
        my $i   = 1;
        while ( $i < $read_count ) {
            $sum += exp _calc_log_sum( $i, $log_avg_reads );
            last if $sum >= $exp_avg_reads;
            $i++;
        }
        my $prob = 1 - exp( -$avg_reads ) * $sum;
        if ( $prob < $arg_ref->{hmm_sig_level} / $total_peaks ) {
            push @sig_peak_widths, $width;
            $total_sig_read_count += $read_count;
            $total_sig_peak_width += $width;
        }

        # Expire Memoize cache for each peak
        flush_cache('_calc_log_sum');
    }

    # Calculate hit threshold
    my $proportion_bp_in_peaks =
      $total_read_count * $arg_ref->{read_length} / $arg_ref->{total_bp};
    my $read_threshold = 0;
    my $prob           = 1;
    my $sum            = 1;
    while ( $prob > $arg_ref->{hmm_sig_level} / $total_peaks ) {
        $read_threshold++;
        my $log_sum = 0;
        foreach my $i ( 1 .. $read_threshold ) {
            $log_sum += log($proportion_bp_in_peaks) - log $i;
        }
        $sum += exp $log_sum;
        $prob = 1 - exp( -$proportion_bp_in_peaks ) * $sum;
    }

    # Sort widths and get median
    my $total_sig_peaks = scalar @sig_peak_widths;
    @sig_peak_widths = sort { $a <=> $b } @sig_peak_widths;
    my $median_sig_peak_width = $sig_peak_widths[ int( $total_sig_peaks / 2 ) ];

    ## no critic (ProhibitMagicNumbers)
    my %summary = (
        total_read_count_per_mb     => $total_read_count / 1_000_000,
        total_sig_read_count_per_mb => $total_sig_read_count / 1_000_000,
        total_sig_peak_width_in_mb  => $total_sig_peak_width / 1_000_000,
        median_sig_peak_width       => $median_sig_peak_width || 0,
        total_sig_peaks             => $total_sig_peaks,
        peak_buffer_width           => $arg_ref->{peak_buffer_width},
        read_threshold              => $read_threshold,
    );
    ## use critic

    return \%summary;
}

# Calculate log sum
sub _calc_log_sum {
    my ( $i, $log ) = @_;

    if ( $i == 0 ) {
        return 0;
    }
    else {
        return $log - log($i) + _calc_log_sum( $i - 1, $log );
    }
}

=func run_peak_hmm

  Usage       : my $hmm_ref = DETCT::Misc::PeakHMM::run_peak_hmm( {
                    dir           => '.',
                    hmm_sig_level => 0.001,
                    seq_name      => '1',
                    seq_bp        => 10_000_000,
                    bin_size      => 100,
                    read_bins     => $read_bins_hash_ref,
                    summary       => $summary_hash_ref,
                    hmm_binary    => 'chiphmmnew',
                } );
  Purpose     : Run peak HMM
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
                        Arrayref [
                            Int (bin),
                            Int (read count),
                            Float (log probability),
                        ],
                        ... (peaks)
                }
  Parameters  : Hashref {
                    dir           => String (the working directory),
                    hmm_sig_level => Float (the HMM significance level),
                    seq_name      => String (the sequence name),
                    seq_bp        => Int (the sequence bp),
                    bin_size      => Int (the bin size),
                    read_bins     => Hashref (of read bins),
                    summary       => Hashref (of summary),
                    hmm_binary    => String (the HMM binary)
                }
  Throws      : If directory is missing
                If HMM significance level is missing
                If sequence name is missing
                If sequence bp is missing
                If bin size is missing
                If read bins are missing
                If summary is missing
                If HMM binary is missing
                If command line can't be run
  Comments    : None

=cut

sub run_peak_hmm {
    my ($arg_ref) = @_;

    confess 'No directory specified' if !defined $arg_ref->{dir};
    confess 'No HMM significance level specified'
      if !defined $arg_ref->{hmm_sig_level};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No sequence bp specified'   if !defined $arg_ref->{seq_bp};
    confess 'No bin size specified'      if !defined $arg_ref->{bin_size};
    confess 'No read bins specified'     if !defined $arg_ref->{read_bins};
    confess 'No summary specified'       if !defined $arg_ref->{summary};
    confess 'No HMM binary specified'    if !defined $arg_ref->{hmm_binary};

    if ( !scalar keys %{ $arg_ref->{summary} } ) {

        # No summary (i.e. no peaks), so won't run HMM
        return { $arg_ref->{seq_name} => [] };
    }

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Store original directory and change to working directory
    my $original_dir = getcwd();
    chdir $arg_ref->{dir};

    # Sanitise sequence name for using in filenames
    my $safe_seq_name = $arg_ref->{seq_name};
    $safe_seq_name =~ s/\W+//xmsg;

    # Write read bins to file
    my $bin_file = $safe_seq_name . '.bins';
    open my $bin_fh, '>', $bin_file;
    foreach my $bin ( sort { $a <=> $b } keys %{ $arg_ref->{read_bins} } ) {
        write_or_die( $bin_fh, $bin, "\t", $arg_ref->{read_bins}->{$bin},
            "\n" );
    }
    close $bin_fh;

    # Write summary to file
    my $sum_file = $safe_seq_name . '.params';
    ## no critic (RequireBriefOpen)
    open my $sum_fh, '>', $sum_file;
    write_or_die( $sum_fh, $arg_ref->{summary}->{total_read_count_per_mb},
        "\n" );
    write_or_die( $sum_fh, $arg_ref->{summary}->{total_sig_read_count_per_mb},
        "\n" );
    write_or_die( $sum_fh, $arg_ref->{summary}->{total_sig_peak_width_in_mb},
        "\n" );
    write_or_die( $sum_fh, $arg_ref->{summary}->{median_sig_peak_width}, "\n" );
    write_or_die( $sum_fh, $arg_ref->{summary}->{total_sig_peaks},       "\n" );
    write_or_die( $sum_fh, $arg_ref->{summary}->{peak_buffer_width},     "\n" );
    write_or_die( $sum_fh, $arg_ref->{summary}->{read_threshold},        "\n" );
    close $sum_fh;
    ## use critic

    my $hmm_file    = $safe_seq_name . '.hmm';
    my $stdout_file = $safe_seq_name . '.o';
    my $stderr_file = $safe_seq_name . '.e';

    my $num_bins = ceil( $arg_ref->{seq_bp} / $arg_ref->{bin_size} );

    my $cmd = join q{ }, $arg_ref->{hmm_binary}, $bin_file, $sum_file,
      $hmm_file, $arg_ref->{bin_size}, $num_bins;
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    # Reformat output into array of arrayrefs
    my $log_hmm_sig_level = log $arg_ref->{hmm_sig_level};
    my @hmm_output        = ();
    if ( -r $hmm_file ) {    # Peak HMM can fail
        foreach my $line ( read_file($hmm_file) ) {
            chomp $line;
            my ( $bin, undef, undef, $read_count, $log_prob ) = split /\t/xms,
              $line;
            next if $log_prob >= $log_hmm_sig_level;
            push @hmm_output, [ $bin, $read_count, $log_prob ];
        }
    }

    # Change back to original directory
    chdir $original_dir;

    return { $arg_ref->{seq_name} => \@hmm_output };
}

=func join_hmm_bins

  Usage       : my $regions_ref = DETCT::Misc::PeakHMM::join_hmm_bins( {
                    bin_size => 100,
                    seq_name => '1',
                    hmm_bins => $hmm_bins_ary_ref,
                } );
  Purpose     : Join reads bins output by peak HMM into regions
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
                        Arrayref [
                            Int (region start),
                            Int (region end),
                            Int (region maximum read count),
                            Float (region log probability sum),
                        ],
                        ... (regions)
                }
  Parameters  : Hashref {
                    bin_size => Int (the bin size),
                    seq_name => String (the sequence name),
                    hmm_bins => Arrayref (of HMM bins),
                }
  Throws      : If bin size is missing
                If sequence name is missing
                If HMM bins are missing
  Comments    : None

=cut

sub join_hmm_bins {
    my ($arg_ref) = @_;

    confess 'No bin size specified'      if !defined $arg_ref->{bin_size};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No HMM bins specified'      if !defined $arg_ref->{hmm_bins};

    my @regions;

    # Region variables (where a region is a set of merged bins)
    my $region_bin_start;
    my $region_bin_end;
    my $region_max_read_count;
    my $region_log_prob_sum;

    foreach my $hmm_bin ( @{ $arg_ref->{hmm_bins} } ) {
        my ( $bin, $read_count, $log_prob ) = @{$hmm_bin};

        # We're starting the first region
        if ( !defined $region_bin_start ) {
            $region_bin_start      = $bin;
            $region_bin_end        = $bin;
            $region_max_read_count = $read_count;
            $region_log_prob_sum   = $log_prob;
            next;
        }

        # Extend or finish current region?
        if ( $bin == $region_bin_end + 1 ) {

            # Next bin, so extend current region
            $region_bin_end = $bin;
            if ( $read_count > $region_max_read_count ) {
                $region_max_read_count = $read_count;
            }
            $region_log_prob_sum += $log_prob;
        }
        else {
            # Finish current region and convert to genomic coordinates
            push @regions,
              [
                $region_bin_start * $arg_ref->{bin_size} + 1,
                ( $region_bin_end + 1 ) * $arg_ref->{bin_size},
                $region_max_read_count,
                $region_log_prob_sum,
              ];

            # Start new region
            $region_bin_start      = $bin;
            $region_bin_end        = $bin;
            $region_max_read_count = $read_count;
            $region_log_prob_sum   = $log_prob;
        }
    }

    # Finish last region
    if ( defined $region_bin_start ) {
        push @regions,
          [
            $region_bin_start * $arg_ref->{bin_size} + 1,
            ( $region_bin_end + 1 ) * $arg_ref->{bin_size},
            $region_max_read_count,
            $region_log_prob_sum,
          ];
    }

    return { $arg_ref->{seq_name} => \@regions };
}

1;
