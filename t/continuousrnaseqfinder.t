use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 22;

use DETCT::ContinuousRNASeqFinder;

my $continuousrnaseq_finder = DETCT::ContinuousRNASeqFinder->new(
    {
        table_file => 't/data/test_transcripts.tsv',
    }
);

isa_ok( $continuousrnaseq_finder, 'DETCT::ContinuousRNASeqFinder' );

# Test table file attribute
is( $continuousrnaseq_finder->table_file,
    't/data/test_transcripts.tsv', 'Get table file' );
throws_ok { $continuousrnaseq_finder->set_table_file('nonexistent') }
qr/cannot be read/ms, 'Unreadable table file';

# Test table format attribute
$continuousrnaseq_finder->set_table_file('t/data/test_transcripts.tsv');
is( $continuousrnaseq_finder->table_format, 'tsv', 'Get guessed table format' );
is( $continuousrnaseq_finder->set_table_format('csv'),
    undef, 'Set table format' );
is( $continuousrnaseq_finder->table_format, 'csv', 'Get new table format' );
throws_ok { $continuousrnaseq_finder->set_table_format('invalid') }
qr/Invalid table format/ms, 'Invalid table format';

# Test getting containing continuous RNA-Seq
$continuousrnaseq_finder = DETCT::ContinuousRNASeqFinder->new(
    {
        table_file => 't/data/test_transcripts.tsv',
    }
);
my @transcripts;

(@transcripts) =
  $continuousrnaseq_finder->get_containing_continuous_rnaseq( '1', 15241, 1 );
is( scalar @transcripts, 0, 'No RNA-Seq' );
(@transcripts) =
  $continuousrnaseq_finder->get_containing_continuous_rnaseq( '1', 15139, 1 );
is( scalar @transcripts,  2,                    '2 transcripts' );
is( $transcripts[0]->[0], 'ENSDART00000109479', '1st transcript' );
is( $transcripts[1]->[0], 'ENSDART00000152276', '2nd transcript' );
(@transcripts) =
  $continuousrnaseq_finder->get_containing_continuous_rnaseq( '1', 15139, -1 );
is( scalar @transcripts, 0, 'No RNA-Seq on other strand' );
(@transcripts) =
  $continuousrnaseq_finder->get_containing_continuous_rnaseq( '1', 47195, -1 );
is( scalar @transcripts,  1,                    '1 transcript' );
is( $transcripts[0]->[0], 'ENSDART00000018741', 'Only transcript' );
(@transcripts) =
  $continuousrnaseq_finder->get_containing_continuous_rnaseq( '1', 47195, 1 );
is( scalar @transcripts, 0, 'No RNA-Seq on other strand' );
(@transcripts) =
  $continuousrnaseq_finder->get_containing_continuous_rnaseq( 'nonexistent', 1,
    1 );
is( scalar @transcripts, 0, 'Non-existent sequence' );

# Adding continuous RNA-Seq annotation
my $regions = [
    [
        1, 10, 10, -10, 1,
        [
            [ '1', 15241, 1,  20, 1, 'A', 'A', 1, 'A', 1, 1 ],
            [ '1', 15139, 1,  10, 1, 'T', 'T', 0, 'C', 1, 1 ],
            [ '1', 47195, -1, 10, 1, 'A', 'A', 1, 'A', 1, 1 ],
        ]
    ]
];
my $annotated_regions =
  $continuousrnaseq_finder->add_continuous_rnaseq($regions);
is( scalar @{ $annotated_regions->[0]->[-1]->[0]->[-1] }, 0, 'No RNA-Seq' );
is( scalar @{ $annotated_regions->[0]->[-1]->[1]->[-1] }, 2, '2 transcripts' );
is( $annotated_regions->[0]->[-1]->[1]->[-1]->[0]->[0],
    'ENSDART00000109479', '1st transcript' );
is( $annotated_regions->[0]->[-1]->[1]->[-1]->[1]->[0],
    'ENSDART00000152276', '2nd transcript' );
is( scalar @{ $annotated_regions->[0]->[-1]->[2]->[-1] }, 1, '1 transcript' );
is( $annotated_regions->[0]->[-1]->[2]->[-1]->[0]->[0],
    'ENSDART00000018741', 'Only transcript' );
