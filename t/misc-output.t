use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 13;

use DETCT::Misc::Output qw(
  dump_as_table
);

use File::Temp qw( tempdir );
use File::Spec;

my $tmp_dir = tempdir( CLEANUP => 1 );

# Mock sample objects
my $sample1 = Test::MockObject->new();
$sample1->set_isa('DETCT::Sample');
$sample1->set_always( 'name',      'wt1' );
$sample1->set_always( 'condition', 'sibling' );
$sample1->set_always( 'group',     '1' );
my $sample2 = Test::MockObject->new();
$sample2->set_isa('DETCT::Sample');
$sample2->set_always( 'name',      'wt2' );
$sample2->set_always( 'condition', 'sibling' );
$sample2->set_always( 'group',     '2' );
my $sample3 = Test::MockObject->new();
$sample3->set_isa('DETCT::Sample');
$sample3->set_always( 'name',      'mut1' );
$sample3->set_always( 'condition', 'mutant' );
$sample3->set_always( 'group',     '1' );
my $sample4 = Test::MockObject->new();
$sample4->set_isa('DETCT::Sample');
$sample4->set_always( 'name',      'mut2' );
$sample4->set_always( 'condition', 'mutant' );
$sample4->set_always( 'group',     '2' );
my $samples = [ $sample1, $sample2, $sample3, $sample4 ];

# Mock analysis object
my $analysis = Test::MockObject->new();
$analysis->set_isa('DETCT::Analysis');
$analysis->set_always( 'get_all_samples', $samples );
$analysis->set_always( 'ensembl_species', 'danio_rerio' );

my $regions = [
    [
        '1', 1, 110, 10, -10, '1', 110,
        1,   10,
        [ 4,   1,   2,   7 ],
        [ 4.6, 1.1, 2.1, 4.6 ],
        undef, undef,
        [ 1.18, 0.233 ],
        [ [ 0.46, -1.13 ], [ 4.18, 2.06 ] ],
        {
            e61 => [
                [
                    'ENSDARG00000095747',
                    'cxc64',
                    'CXC chemokine 64',
                    'protein_coding',
                    5,
                    [ [ 'ENSDART00000133571', 'protein_coding' ] ]
                ]
            ]
        }
    ],
    [
        '1',   1,
        1000,  10,
        -10,   undef,
        undef, undef,
        undef, [ 4, 1, 2, 7 ],
        [ 4.6, 1.1, 2.1, 4.6 ], undef,
        undef, [ 1.18, 0.233 ],
        [ [ 0.46, -1.13 ], [ 4.18, 2.06 ] ], {}
    ],
];

is(
    dump_as_table(
        { analysis => $analysis, dir => $tmp_dir, regions => $regions, }
    ),
    undef, 'Dump'
);

foreach my $format ( 'csv', 'tsv', 'html' ) {
    foreach my $level ( 'all', 'sig' ) {
        my $file = $level . q{.} . $format;
        my $filepath = File::Spec->catfile( $tmp_dir, $file );
        ok( -e $filepath,  $file . ' exists' );
        ok( !-z $filepath, $file . ' is not empty' );
    }
}

# TODO: Actually test output
