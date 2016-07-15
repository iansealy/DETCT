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

use DETCT::TransposonFinder;

# Mock repeat features
my @repeat_features;
my @repeat_parameters =
  ( [ 10, 20, 'Transposon' ], [ 30, 40, 'Not' ], [ 50, 60, 'Transposon' ] );
foreach my $repeat_parameters (@repeat_parameters) {
    my ( $start, $end, $type ) = @{$repeat_parameters};

    # Create reepat consensus
    my $consensus = Test::MockObject->new();
    $consensus->set_always( 'repeat_type', $type );

    # Create repeat feature
    my $feature = Test::MockObject->new();
    $feature->set_always( 'seq_region_start', $start );
    $feature->set_always( 'seq_region_end',   $end );
    $feature->set_always( 'repeat_consensus', $consensus );
    push @repeat_features, $feature;
}

# Mock slice
my $slice = Test::MockObject->new();
$slice->set_always( 'get_all_RepeatFeatures', \@repeat_features );

# Mock slice adaptor
my $slice_adaptor = Test::MockObject->new();
$slice_adaptor->set_isa('Bio::EnsEMBL::DBSQL::SliceAdaptor');
$slice_adaptor->set_always( 'fetch_by_region', $slice );

my $transposon_finder =
  DETCT::TransposonFinder->new( { slice_adaptor => $slice_adaptor, } );

isa_ok( $transposon_finder, 'DETCT::TransposonFinder' );

# Test Ensembl slice adaptor attribute
isa_ok( $transposon_finder->slice_adaptor,
    'Bio::EnsEMBL::DBSQL::SliceAdaptor' );
throws_ok { $transposon_finder->set_slice_adaptor() }
qr/No Ensembl slice adaptor specified/ms, 'No Ensembl slice adaptor';
throws_ok { $transposon_finder->set_slice_adaptor('invalid') }
qr/Class of Ensembl slice adaptor/ms, 'Invalid Ensembl slice adaptor';

my $distance;
my $nearest_transposon_pos;

# Nearest to single transposon
( $distance, $nearest_transposon_pos ) =
  $transposon_finder->get_nearest_transposon( '1', 1, 1 );
is( $distance,               9,  q{Transposon is 9 bp downstream} );
is( $nearest_transposon_pos, 10, q{Transposon starts at 10 bp} );
( $distance, $nearest_transposon_pos ) =
  $transposon_finder->get_nearest_transposon( '1', 1, -1 );
is( $distance,               -9, q{Transposon is 9 bp upstream} );
is( $nearest_transposon_pos, 10, q{Transposon starts at 10 bp} );

# Nearest to non-transposon repeat
( $distance, $nearest_transposon_pos ) =
  $transposon_finder->get_nearest_transposon( '1', 29, 1 );
is( $distance,               -9, q{Transposon is 9 bp upstream} );
is( $nearest_transposon_pos, 20, q{Transposon ends at 20 bp} );
( $distance, $nearest_transposon_pos ) =
  $transposon_finder->get_nearest_transposon( '1', 29, -1 );
is( $distance,               9,  q{Transposon is 9 bp downstream} );
is( $nearest_transposon_pos, 20, q{Transposon ends at 20 bp} );
( $distance, $nearest_transposon_pos ) =
  $transposon_finder->get_nearest_transposon( '1', 41, 1 );
is( $distance,               9,  q{Transposon is 9 bp downstream} );
is( $nearest_transposon_pos, 50, q{Transposon starts at 50 bp} );
( $distance, $nearest_transposon_pos ) =
  $transposon_finder->get_nearest_transposon( '1', 41, -1 );
is( $distance,               -9, q{Transposon is 9 bp upstream} );
is( $nearest_transposon_pos, 50, q{Transposon starts at 50 bp} );

# Adding transposon annotation
my $regions = [
    [
        1, 10, 10, -10, 1,
        [
            [ '1', 10, 1,  20, 1, 'A', 'A', 1, 'A' ],
            [ '1', 20, -1, 10, 1, 'T', 'T', 0, 'C' ],
            [ '1', 30, 1,  10, 1, 'A', 'A', 1, 'A' ],
        ]
    ]
];
my $annotated_regions = $transposon_finder->add_transposon_annotation($regions);
is( $annotated_regions->[0]->[-1]->[0]->[-2], 0,   'Distance to transposon' );
is( $annotated_regions->[0]->[-1]->[0]->[-1], 10,  'Position of transposon' );
is( $annotated_regions->[0]->[-1]->[1]->[-2], 0,   'Distance to transposon' );
is( $annotated_regions->[0]->[-1]->[1]->[-1], 20,  'Position of transposon' );
is( $annotated_regions->[0]->[-1]->[2]->[-2], -10, 'Distance to transposon' );
is( $annotated_regions->[0]->[-1]->[2]->[-1], 20,  'Position of transposon' );
