use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 103;

use DETCT::Misc::PeakHMM qw(
  merge_read_peaks
  summarise_read_peaks
  run_peak_hmm
  join_hmm_bins
);

use File::Temp qw( tempdir );
use File::Path qw( make_path );
use POSIX qw( WIFEXITED);

# Compile quince_chiphmmnew if necessary
if ( !-r 'bin/quince_chiphmmnew' ) {
    make_path('bin');
    my $cmd = 'g++ -o bin/quince_chiphmmnew src/quince_chiphmmnew.cpp';
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd";
}

my $input_peaks;
my $output_peaks;

# Check merging read peaks required parameters
throws_ok {
    merge_read_peaks(
        {
            seq_name => 1,
            peaks    => [ [ 1, 2, 1 ] ],
        }
    );
}
qr/No peak buffer width specified/ms, 'No peak buffer width';
throws_ok {
    merge_read_peaks(
        {
            peak_buffer_width => 100,
            peaks             => [ [ 1, 2, 1 ] ],
        }
    );
}
qr/No sequence name specified/ms, 'No sequence name';
throws_ok {
    merge_read_peaks(
        {
            peak_buffer_width => 100,
            seq_name          => 1,
        }
    );
}
qr/No peaks specified/ms, 'No peaks';

# Two peaks but no merging
$input_peaks = [ [ 100, 200, 1 ], [ 500, 600, 1 ], ];
$output_peaks = merge_read_peaks(
    {
        peak_buffer_width => 100,
        seq_name          => 1,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$output_peaks},     1,   '1 sequence' );
is( scalar @{ $output_peaks->{'1'} }, 2,   '2 peaks' );
is( $output_peaks->{'1'}->[0]->[0],   100, 'Start of first peak' );
is( $output_peaks->{'1'}->[0]->[1],   200, 'End of first peak' );
is( $output_peaks->{'1'}->[0]->[2],   1,   'First peak read count' );
is( $output_peaks->{'1'}->[-1]->[0],  500, 'Start of last peak' );
is( $output_peaks->{'1'}->[-1]->[1],  600, 'End of last peak' );
is( $output_peaks->{'1'}->[-1]->[2],  1,   'Last peak read count' );

# Two peaks merged into one
$input_peaks = [ [ 100, 200, 1 ], [ 250, 350, 1 ], ];
$output_peaks = merge_read_peaks(
    {
        peak_buffer_width => 100,
        seq_name          => 1,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$output_peaks},     1,   '1 sequence' );
is( scalar @{ $output_peaks->{'1'} }, 1,   '1 peak' );
is( $output_peaks->{'1'}->[0]->[0],   100, 'Start of first peak' );
is( $output_peaks->{'1'}->[0]->[1],   350, 'End of first peak' );
is( $output_peaks->{'1'}->[0]->[2],   2,   'First peak read count' );

# Three peaks with first two merged
$input_peaks = [ [ 100, 200, 1 ], [ 250, 350, 1 ], [ 500, 600, 1 ], ];
$output_peaks = merge_read_peaks(
    {
        peak_buffer_width => 100,
        seq_name          => 1,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$output_peaks},     1,   '1 sequence' );
is( scalar @{ $output_peaks->{'1'} }, 2,   '2 peaks' );
is( $output_peaks->{'1'}->[0]->[0],   100, 'Start of first peak' );
is( $output_peaks->{'1'}->[0]->[1],   350, 'End of first peak' );
is( $output_peaks->{'1'}->[0]->[2],   2,   'First peak read count' );
is( $output_peaks->{'1'}->[-1]->[0],  500, 'Start of last peak' );
is( $output_peaks->{'1'}->[-1]->[1],  600, 'End of last peak' );
is( $output_peaks->{'1'}->[-1]->[2],  1,   'Last peak read count' );

# Three peaks with second two merged
$input_peaks = [ [ 100, 200, 1 ], [ 500, 600, 1 ], [ 550, 650, 1 ], ];
$output_peaks = merge_read_peaks(
    {
        peak_buffer_width => 100,
        seq_name          => 1,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$output_peaks},     1,   '1 sequence' );
is( scalar @{ $output_peaks->{'1'} }, 2,   '2 peaks' );
is( $output_peaks->{'1'}->[0]->[0],   100, 'Start of first peak' );
is( $output_peaks->{'1'}->[0]->[1],   200, 'End of first peak' );
is( $output_peaks->{'1'}->[0]->[2],   1,   'First peak read count' );
is( $output_peaks->{'1'}->[-1]->[0],  500, 'Start of last peak' );
is( $output_peaks->{'1'}->[-1]->[1],  650, 'End of last peak' );
is( $output_peaks->{'1'}->[-1]->[2],  2,   'Last peak read count' );

# Two peaks separated by buffer width
$input_peaks = [ [ 100, 200, 1 ], [ 300, 400, 1 ], ];
$output_peaks = merge_read_peaks(
    {
        peak_buffer_width => 100,
        seq_name          => 1,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$output_peaks},     1,   '1 sequence' );
is( scalar @{ $output_peaks->{'1'} }, 2,   '2 peaks' );
is( $output_peaks->{'1'}->[0]->[0],   100, 'Start of first peak' );
is( $output_peaks->{'1'}->[0]->[1],   200, 'End of first peak' );
is( $output_peaks->{'1'}->[0]->[2],   1,   'First peak read count' );
is( $output_peaks->{'1'}->[-1]->[0],  300, 'Start of last peak' );
is( $output_peaks->{'1'}->[-1]->[1],  400, 'End of last peak' );
is( $output_peaks->{'1'}->[-1]->[2],  1,   'Last peak read count' );

# Two peaks separated by just under buffer width
$input_peaks = [ [ 100, 200, 1 ], [ 299, 400, 1 ], ];
$output_peaks = merge_read_peaks(
    {
        peak_buffer_width => 100,
        seq_name          => 1,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$output_peaks},     1,   '1 sequence' );
is( scalar @{ $output_peaks->{'1'} }, 1,   '1 peak' );
is( $output_peaks->{'1'}->[0]->[0],   100, 'Start of first peak' );
is( $output_peaks->{'1'}->[0]->[1],   400, 'End of first peak' );
is( $output_peaks->{'1'}->[0]->[2],   2,   'First peak read count' );

# Two peaks with same start
$input_peaks = [ [ 100, 200, 1 ], [ 100, 300, 1 ], ];
$output_peaks = merge_read_peaks(
    {
        peak_buffer_width => 100,
        seq_name          => 1,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$output_peaks},     1,   '1 sequence' );
is( scalar @{ $output_peaks->{'1'} }, 1,   '1 peak' );
is( $output_peaks->{'1'}->[0]->[0],   100, 'Start of first peak' );
is( $output_peaks->{'1'}->[0]->[1],   300, 'End of first peak' );
is( $output_peaks->{'1'}->[0]->[2],   2,   'First peak read count' );

# No peaks
$input_peaks  = [];
$output_peaks = merge_read_peaks(
    {
        peak_buffer_width => 100,
        seq_name          => 1,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$output_peaks},     1, '1 sequence' );
is( scalar @{ $output_peaks->{'1'} }, 0, '0 peaks' );

my $summary;

# Check summarising read peaks required parameters
throws_ok {
    summarise_read_peaks(
        {
            peak_buffer_width => 100,
            hmm_sig_level     => 0.001,
            seq_name          => '1',
            seq_bp            => 1000,
            read_length       => 54,
            peaks             => [ [ 1, 2, 1 ] ],
        }
    );
}
qr/No bin size specified/ms, 'No bin size';
throws_ok {
    summarise_read_peaks(
        {
            bin_size      => 100,
            hmm_sig_level => 0.001,
            seq_name      => '1',
            seq_bp        => 1000,
            read_length   => 54,
            peaks         => [ [ 1, 2, 1 ] ],
        }
    );
}
qr/No peak buffer width specified/ms, 'No peak buffer width';
throws_ok {
    summarise_read_peaks(
        {
            bin_size          => 100,
            peak_buffer_width => 100,
            seq_name          => '1',
            seq_bp            => 1000,
            read_length       => 54,
            peaks             => [ [ 1, 2, 1 ] ],
        }
    );
}
qr/No HMM significance level specified/ms, 'No HMM significance level';
throws_ok {
    summarise_read_peaks(
        {
            bin_size          => 100,
            peak_buffer_width => 100,
            hmm_sig_level     => 0.001,
            seq_bp            => 1000,
            read_length       => 54,
            peaks             => [ [ 1, 2, 1 ] ],
        }
    );
}
qr/No sequence name specified/ms, 'No sequence name';
throws_ok {
    summarise_read_peaks(
        {
            bin_size          => 100,
            peak_buffer_width => 100,
            hmm_sig_level     => 0.001,
            seq_name          => '1',
            read_length       => 54,
            peaks             => [ [ 1, 2, 1 ] ],
        }
    );
}
qr/No sequence bp specified/ms, 'No sequence bp';
throws_ok {
    summarise_read_peaks(
        {
            bin_size          => 100,
            peak_buffer_width => 100,
            hmm_sig_level     => 0.001,
            seq_name          => '1',
            seq_bp            => 1000,
            peaks             => [ [ 1, 2, 1 ] ],
        }
    );
}
qr/No read length specified/ms, 'No read length';
throws_ok {
    summarise_read_peaks(
        {
            bin_size          => 100,
            peak_buffer_width => 100,
            hmm_sig_level     => 0.001,
            seq_name          => '1',
            seq_bp            => 1000,
            read_length       => 54,
        }
    );
}
qr/No peaks specified/ms, 'No peaks';

# Two peaks, one significant
$input_peaks = [ [ 100, 199, 5 ], [ 300, 399, 1 ], ];
$summary = summarise_read_peaks(
    {
        bin_size          => 100,
        peak_buffer_width => 100,
        hmm_sig_level     => 0.001,
        seq_name          => '1',
        seq_bp            => 1000,
        read_length       => 54,
        peaks             => $input_peaks,
    }
);
is( scalar keys %{$summary},          1, '1 sequence' );
is( scalar keys %{ $summary->{'1'} }, 9, '9 keys' );
is(
    $summary->{'1'}->{total_read_count_per_mb},
    6 / 1_000_000,
    'Total read count per Mb'
);
is(
    $summary->{'1'}->{total_sig_read_count_per_mb},
    5 / 1_000_000,
    'Total significant read count per Mb'
);
is(
    $summary->{'1'}->{total_sig_peak_width_in_mb},
    100 / 1_000_000,
    'Total significant peak width in Mb'
);
is( $summary->{'1'}->{median_sig_peak_width},
    100, 'Median significant peak width' );
is( $summary->{'1'}->{total_sig_peaks},   1,   'Total significant peaks' );
is( $summary->{'1'}->{peak_buffer_width}, 100, 'Peak buffer width' );
ok( $summary->{'1'}->{read_threshold} < 5, 'Read threshold' );
is( $summary->{'1'}->{bin_size}, 100, 'Bin size' );
is( $summary->{'1'}->{num_bins}, 10,  'Number of bins' );

# No significant peaks
$input_peaks = [ [ 300, 399, 1 ], ];
$summary = summarise_read_peaks(
    {
        bin_size          => 100,
        peak_buffer_width => 100,
        hmm_sig_level     => 0.001,
        seq_name          => '1',
        seq_bp            => 1000,
        read_length       => 54,
        peaks             => $input_peaks,
    }
);
is( $summary->{'1'}->{median_sig_peak_width},
    0, 'Median significant peak width' );

# Three significant peaks
$input_peaks = [ [ 100, 149, 500 ], [ 300, 399, 500 ], [ 600, 759, 500 ], ];
$summary = summarise_read_peaks(
    {
        bin_size          => 100,
        peak_buffer_width => 100,
        hmm_sig_level     => 0.001,
        seq_name          => '1',
        seq_bp            => 1000,
        read_length       => 54,
        peaks             => $input_peaks,
    }
);
is( $summary->{'1'}->{median_sig_peak_width},
    100, 'Median significant peak width' );

# No peaks
$summary = summarise_read_peaks(
    {
        bin_size          => 100,
        peak_buffer_width => 100,
        hmm_sig_level     => 0.001,
        seq_name          => '1',
        seq_bp            => 1000,
        read_length       => 54,
        peaks             => [],
    }
);
is( scalar keys %{ $summary->{'1'} }, 0, 'No summary' );

my $tmp_dir = tempdir( CLEANUP => 1 );

# Check running peak HMM required parameters
my $read_bins = {
    1 => 500,
    3 => 1,
};
$summary = {
    total_read_count_per_mb     => 501 / 1_000_000,
    total_sig_read_count_per_mb => 500 / 1_000_000,
    total_sig_peak_width_in_mb  => 100 / 1_000_000,
    median_sig_peak_width       => 100,
    total_sig_peaks             => 1,
    peak_buffer_width           => 100,
    read_threshold              => 3,
    bin_size                    => 100,
    num_bins                    => 10,
};
throws_ok {
    run_peak_hmm(
        {
            hmm_sig_level => 0.001,
            seq_name      => '1',
            read_bins     => $read_bins,
            summary       => $summary,
            hmm_binary    => 'bin/quince_chiphmmnew',
        }
    );
}
qr/No directory specified/ms, 'No directory';
throws_ok {
    run_peak_hmm(
        {
            dir        => $tmp_dir,
            seq_name   => '1',
            read_bins  => $read_bins,
            summary    => $summary,
            hmm_binary => 'bin/quince_chiphmmnew',
        }
    );
}
qr/No HMM significance level specified/ms, 'No HMM significance level';
throws_ok {
    run_peak_hmm(
        {
            dir           => $tmp_dir,
            hmm_sig_level => 0.001,
            read_bins     => $read_bins,
            summary       => $summary,
            hmm_binary    => 'bin/quince_chiphmmnew',
        }
    );
}
qr/No sequence name specified/ms, 'No sequence name';
throws_ok {
    run_peak_hmm(
        {
            dir           => $tmp_dir,
            hmm_sig_level => 0.001,
            seq_name      => '1',
            summary       => $summary,
            hmm_binary    => 'bin/quince_chiphmmnew',
        }
    );
}
qr/No read bins specified/ms, 'No read bins';
throws_ok {
    run_peak_hmm(
        {
            dir           => $tmp_dir,
            hmm_sig_level => 0.001,
            seq_name      => '1',
            read_bins     => $read_bins,
            hmm_binary    => 'bin/quince_chiphmmnew',
        }
    );
}
qr/No summary specified/ms, 'No summary';
throws_ok {
    run_peak_hmm(
        {
            dir           => $tmp_dir,
            hmm_sig_level => 0.001,
            seq_name      => '1',
            read_bins     => $read_bins,
            summary       => $summary,
        }
    );
}
qr/No HMM binary specified/ms, 'No HMM binary';

my $hmm;

# Run peak HMM
$hmm = run_peak_hmm(
    {
        dir           => $tmp_dir,
        hmm_sig_level => 0.001,
        seq_name      => '1',
        read_bins     => $read_bins,
        summary       => $summary,
        hmm_binary    => 'bin/quince_chiphmmnew',
    }
);
is( scalar keys %{$hmm},     1,   '1 sequence' );
is( scalar @{ $hmm->{'1'} }, 1,   '1 peak' );
is( $hmm->{'1'}->[0]->[0],   1,   'Bin 1' );
is( $hmm->{'1'}->[0]->[1],   500, '500 reads' );
ok( $hmm->{'1'}->[0]->[2] < 0, 'Log probability negative' );

# Run peak HMM with no summary
$hmm = run_peak_hmm(
    {
        dir           => $tmp_dir,
        hmm_sig_level => 0.001,
        seq_name      => '1',
        read_bins     => $read_bins,
        summary       => {},
        hmm_binary    => 'bin/quince_chiphmmnew',
    }
);
is( scalar keys %{$hmm},     1, '1 sequence' );
is( scalar @{ $hmm->{'1'} }, 0, '0 peaks' );

# Run peak HMM with non-existent working directory
$hmm = run_peak_hmm(
    {
        dir           => $tmp_dir . '/test',
        hmm_sig_level => 0.001,
        seq_name      => '1',
        read_bins     => $read_bins,
        summary       => $summary,
        hmm_binary    => 'bin/quince_chiphmmnew',
    }
);
is( scalar keys %{$hmm},     1, '1 sequence' );
is( scalar @{ $hmm->{'1'} }, 1, '1 peak' );

# Check running peak HMM required parameters
my $hmm_bins = [
    [ 1, 10, -2.3 ],
    [ 2, 20, -2.3 ],
    [ 4, 10, -2.3 ],
    [ 5, 30, -2.3 ],
    [ 6, 20, -2.3 ],
];
throws_ok {
    join_hmm_bins(
        {
            seq_name => '1',
            hmm_bins => $hmm_bins,
        }
    );
}
qr/No bin size specified/ms, 'No bin size';
throws_ok {
    join_hmm_bins(
        {
            bin_size => 100,
            hmm_bins => $hmm_bins,
        }
    );
}
qr/No sequence name specified/ms, 'No sequence name';
throws_ok {
    join_hmm_bins(
        {
            bin_size => 100,
            seq_name => '1',
        }
    );
}
qr/No HMM bins specified/ms, 'No HMM bins';

my $regions;

# Five peaks joined to two regions
$regions = join_hmm_bins(
    {
        bin_size => 100,
        seq_name => '1',
        hmm_bins => $hmm_bins,
    }
);
is( scalar keys %{$regions},     1,    '1 sequence' );
is( scalar @{ $regions->{'1'} }, 2,    '2 peaks' );
is( $regions->{'1'}->[0]->[0],   101,  'Region 1 start' );
is( $regions->{'1'}->[0]->[1],   300,  'Region 1 end' );
is( $regions->{'1'}->[0]->[2],   20,   'Region 1 max read count' );
is( $regions->{'1'}->[0]->[3],   -4.6, 'Region 1 log probability sum' );
is( $regions->{'1'}->[1]->[0],   401,  'Region 2 start' );
is( $regions->{'1'}->[1]->[1],   700,  'Region 2 end' );
is( $regions->{'1'}->[1]->[2],   30,   'Region 2 max read count' );
is( $regions->{'1'}->[1]->[3],   -6.9, 'Region 2 log probability sum' );

# No peaks
$regions = join_hmm_bins(
    {
        bin_size => 100,
        seq_name => '1',
        hmm_bins => [],
    }
);
is( scalar keys %{$regions},     1, '1 sequence' );
is( scalar @{ $regions->{'1'} }, 0, '0 peaks' );
