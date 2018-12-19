use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 233;

use DETCT::Misc::Output qw(
  dump_as_table
  parse_table
);

use File::Temp qw( tempdir );
use File::Spec;

my $tmp_dir = tempdir( CLEANUP => 1 );

# Mock sample objects
my $sample1 = Test::MockObject->new();
$sample1->set_isa('DETCT::Sample');
$sample1->set_always( 'name',      'wt1' );
$sample1->set_always( 'condition', 'sibling' );
$sample1->set_always( 'group',  ['1'] );
$sample1->set_always( 'groups', ['1'] );
my $sample2 = Test::MockObject->new();
$sample2->set_isa('DETCT::Sample');
$sample2->set_always( 'name',      'wt2' );
$sample2->set_always( 'condition', 'sibling' );
$sample2->set_always( 'group',  ['2'] );
$sample2->set_always( 'groups', ['2'] );
my $sample3 = Test::MockObject->new();
$sample3->set_isa('DETCT::Sample');
$sample3->set_always( 'name',      'mut1' );
$sample3->set_always( 'condition', 'mutant' );
$sample3->set_always( 'group',  ['1'] );
$sample3->set_always( 'groups', ['1'] );
my $sample4 = Test::MockObject->new();
$sample4->set_isa('DETCT::Sample');
$sample4->set_always( 'name',      'mut2' );
$sample4->set_always( 'condition', 'mutant' );
$sample4->set_always( 'group',  ['2'] );
$sample4->set_always( 'groups', ['2'] );
my $samples = [ $sample1, $sample2, $sample3, $sample4 ];

# Mock analysis object
my $analysis = Test::MockObject->new();
$analysis->set_isa('DETCT::Analysis');
$analysis->set_always( 'get_all_samples', $samples );
$analysis->set_list( 'list_all_conditions', 'mutant', 'sibling' );
$analysis->set_list( 'list_all_groups',     '1',      '2' );
$analysis->set_always( 'ensembl_species',        'danio_rerio' );
$analysis->set_always( 'control_condition',      undef );
$analysis->set_always( 'experimental_condition', undef );

my $regions;

$regions = [
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
        my $file     = $level . q{.} . $format;
        my $filepath = File::Spec->catfile( $tmp_dir, $file );
        ok( -e $filepath,  $file . ' exists' );
        ok( !-z $filepath, $file . ' is not empty' );
    }
}

# Parse output

my $new_regions;

$new_regions = parse_table(
    {
        analysis     => $analysis,
        table_file   => File::Spec->catfile( $tmp_dir, 'all.tsv' ),
        table_format => 'tsv',
    }
);

is( $new_regions->[0]->[0],       '1',   'Region sequence name' );
is( $new_regions->[0]->[1],       1,     'Region start' );
is( $new_regions->[0]->[2],       110,   'Region end' );
is( $new_regions->[0]->[3],       undef, 'Region maximum read count' );
is( $new_regions->[0]->[4],       undef, 'Region log probability sum' );
is( $new_regions->[0]->[5],       '1',   q{3' end sequence name} );
is( $new_regions->[0]->[6],       110,   q{3' end position} );
is( $new_regions->[0]->[7],       1,     q{3' end strand} );
is( $new_regions->[0]->[8],       10,    q{3' end read count} );
is( $new_regions->[0]->[9]->[0],  4,     'Count' );
is( $new_regions->[0]->[9]->[1],  1,     'Count' );
is( $new_regions->[0]->[9]->[2],  2,     'Count' );
is( $new_regions->[0]->[9]->[3],  7,     'Count' );
is( $new_regions->[0]->[10]->[0], 4.6,   'Normalised count' );
is( $new_regions->[0]->[10]->[1], 1.1,   'Normalised count' );
is( $new_regions->[0]->[10]->[2], 2.1,   'Normalised count' );
is( $new_regions->[0]->[10]->[3], 4.6,   'Normalised count' );
is( $new_regions->[0]->[11],      undef, 'p value' );
is( $new_regions->[0]->[12],      undef, 'Adjusted p value' );
is( ( sprintf '%.2f', $new_regions->[0]->[13]->[0] ),
    1.18, 'Condition fold change' );
is( ( sprintf '%.3f', $new_regions->[0]->[13]->[1] ),
    0.233, 'Condition fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[0]->[0] ),
    0.46, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[0]->[1] ),
    -1.13, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[1]->[0] ),
    4.18, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[1]->[1] ),
    2.06, 'Group fold change' );
is( scalar keys %{ $new_regions->[0]->[15] }, 1, 'Gene annotation' );
is( $new_regions->[0]->[15]->{e61}->[0]->[0], 'ENSDARG00000095747', 'Gene ID' );
is( $new_regions->[0]->[15]->{e61}->[0]->[1], 'cxc64', 'Gene name' );
is(
    $new_regions->[0]->[15]->{e61}->[0]->[2],
    'CXC chemokine 64',
    'Gene description'
);
is( $new_regions->[0]->[15]->{e61}->[0]->[3], 'protein_coding', 'Gene type' );
is( $new_regions->[0]->[15]->{e61}->[0]->[4], 5, q{3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript ID' );
is( $new_regions->[0]->[15]->{e61}->[0]->[5]->[0]->[1],
    'protein_coding', 'Transcript type' );

is( $new_regions->[1]->[0],       '1',   'Region sequence name' );
is( $new_regions->[1]->[1],       1,     'Region start' );
is( $new_regions->[1]->[2],       1000,  'Region end' );
is( $new_regions->[1]->[3],       undef, 'Region maximum read count' );
is( $new_regions->[1]->[4],       undef, 'Region log probability sum' );
is( $new_regions->[1]->[5],       undef, q{3' end sequence name} );
is( $new_regions->[1]->[6],       undef, q{3' end position} );
is( $new_regions->[1]->[7],       undef, q{3' end strand} );
is( $new_regions->[1]->[8],       undef, q{3' end read count} );
is( $new_regions->[1]->[9]->[0],  4,     'Count' );
is( $new_regions->[1]->[9]->[1],  1,     'Count' );
is( $new_regions->[1]->[9]->[2],  2,     'Count' );
is( $new_regions->[1]->[9]->[3],  7,     'Count' );
is( $new_regions->[1]->[10]->[0], 4.6,   'Normalised count' );
is( $new_regions->[1]->[10]->[1], 1.1,   'Normalised count' );
is( $new_regions->[1]->[10]->[2], 2.1,   'Normalised count' );
is( $new_regions->[1]->[10]->[3], 4.6,   'Normalised count' );
is( $new_regions->[1]->[11],      undef, 'p value' );
is( $new_regions->[1]->[12],      undef, 'Adjusted p value' );
is( ( sprintf '%.2f', $new_regions->[1]->[13]->[0] ),
    1.18, 'Condition fold change' );
is( ( sprintf '%.3f', $new_regions->[1]->[13]->[1] ),
    0.233, 'Condition fold change' );
is( ( sprintf '%.2f', $new_regions->[1]->[14]->[0]->[0] ),
    0.46, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[1]->[14]->[0]->[1] ),
    -1.13, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[1]->[14]->[1]->[0] ),
    4.18, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[1]->[14]->[1]->[1] ),
    2.06, 'Group fold change' );
is( scalar keys %{ $new_regions->[1]->[15] }, 0, 'No gene annotation' );

# Dump as table if multiple 3' end positions for single region

$regions = [
    [
        '1', 1, 110, 10, -10, '1',
        [ 110, 120 ],
        1,
        [ 10, 20 ],
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
                    [ 5, 15 ],
                    [ [ 'ENSDART00000133571', 'protein_coding' ] ]
                ]
            ]
        }
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
        my $file     = $level . q{.} . $format;
        my $filepath = File::Spec->catfile( $tmp_dir, $file );
        ok( -e $filepath,  $file . ' exists' );
        ok( !-z $filepath, $file . ' is not empty' );
    }
}

# Parse output if multiple 3' end positions for single region

$new_regions = parse_table(
    {
        analysis     => $analysis,
        table_file   => File::Spec->catfile( $tmp_dir, 'all.tsv' ),
        table_format => 'tsv',
    }
);

is( $new_regions->[0]->[0],       '1',   'Region sequence name' );
is( $new_regions->[0]->[1],       1,     'Region start' );
is( $new_regions->[0]->[2],       110,   'Region end' );
is( $new_regions->[0]->[3],       undef, 'Region maximum read count' );
is( $new_regions->[0]->[4],       undef, 'Region log probability sum' );
is( $new_regions->[0]->[5],       '1',   q{3' end sequence name} );
is( $new_regions->[0]->[6]->[0],  110,   q{Nearest 3' end position} );
is( $new_regions->[0]->[6]->[1],  120,   q{Furthest 3' end position} );
is( $new_regions->[0]->[7],       1,     q{3' end strand} );
is( $new_regions->[0]->[8]->[0],  10,    q{Nearest 3' end read count} );
is( $new_regions->[0]->[8]->[1],  20,    q{Furthest 3' end read count} );
is( $new_regions->[0]->[9]->[0],  4,     'Count' );
is( $new_regions->[0]->[9]->[1],  1,     'Count' );
is( $new_regions->[0]->[9]->[2],  2,     'Count' );
is( $new_regions->[0]->[9]->[3],  7,     'Count' );
is( $new_regions->[0]->[10]->[0], 4.6,   'Normalised count' );
is( $new_regions->[0]->[10]->[1], 1.1,   'Normalised count' );
is( $new_regions->[0]->[10]->[2], 2.1,   'Normalised count' );
is( $new_regions->[0]->[10]->[3], 4.6,   'Normalised count' );
is( $new_regions->[0]->[11],      undef, 'p value' );
is( $new_regions->[0]->[12],      undef, 'Adjusted p value' );
is( ( sprintf '%.2f', $new_regions->[0]->[13]->[0] ),
    1.18, 'Condition fold change' );
is( ( sprintf '%.3f', $new_regions->[0]->[13]->[1] ),
    0.233, 'Condition fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[0]->[0] ),
    0.46, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[0]->[1] ),
    -1.13, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[1]->[0] ),
    4.18, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[1]->[1] ),
    2.06, 'Group fold change' );
is( scalar keys %{ $new_regions->[0]->[15] }, 1, 'Gene annotation' );
is( $new_regions->[0]->[15]->{e61}->[0]->[0], 'ENSDARG00000095747', 'Gene ID' );
is( $new_regions->[0]->[15]->{e61}->[0]->[1], 'cxc64', 'Gene name' );
is(
    $new_regions->[0]->[15]->{e61}->[0]->[2],
    'CXC chemokine 64',
    'Gene description'
);
is( $new_regions->[0]->[15]->{e61}->[0]->[3], 'protein_coding', 'Gene type' );
is( $new_regions->[0]->[15]->{e61}->[0]->[4]->[0],
    5, q{Nearest 3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[0]->[4]->[1],
    15, q{Furthest 3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript ID' );
is( $new_regions->[0]->[15]->{e61}->[0]->[5]->[0]->[1],
    'protein_coding', 'Transcript type' );

# Dump as table if multiple 3' end positions and genes with single distance

$regions = [
    [
        '1', 1, 110, 10, -10, '1',
        [ 110, 120 ],
        1,
        [ 10, 20 ],
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
                ],
                [
                    'ENSDARG00000024771',
                    'slc24a5',
                    'solute carrier family 24',
                    'protein_coding',
                    15,
                    [ [ 'ENSDART00000033574', 'protein_coding' ] ]
                ]
            ]
        }
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
        my $file     = $level . q{.} . $format;
        my $filepath = File::Spec->catfile( $tmp_dir, $file );
        ok( -e $filepath,  $file . ' exists' );
        ok( !-z $filepath, $file . ' is not empty' );
    }
}

# Parse output if multiple 3' end positions and genes with single distance

$new_regions = parse_table(
    {
        analysis     => $analysis,
        table_file   => File::Spec->catfile( $tmp_dir, 'all.tsv' ),
        table_format => 'tsv',
    }
);

is( $new_regions->[0]->[0],       '1',   'Region sequence name' );
is( $new_regions->[0]->[1],       1,     'Region start' );
is( $new_regions->[0]->[2],       110,   'Region end' );
is( $new_regions->[0]->[3],       undef, 'Region maximum read count' );
is( $new_regions->[0]->[4],       undef, 'Region log probability sum' );
is( $new_regions->[0]->[5],       '1',   q{3' end sequence name} );
is( $new_regions->[0]->[6]->[0],  110,   q{Nearest 3' end position} );
is( $new_regions->[0]->[6]->[1],  120,   q{Furthest 3' end position} );
is( $new_regions->[0]->[7],       1,     q{3' end strand} );
is( $new_regions->[0]->[8]->[0],  10,    q{Nearest 3' end read count} );
is( $new_regions->[0]->[8]->[1],  20,    q{Furthest 3' end read count} );
is( $new_regions->[0]->[9]->[0],  4,     'Count' );
is( $new_regions->[0]->[9]->[1],  1,     'Count' );
is( $new_regions->[0]->[9]->[2],  2,     'Count' );
is( $new_regions->[0]->[9]->[3],  7,     'Count' );
is( $new_regions->[0]->[10]->[0], 4.6,   'Normalised count' );
is( $new_regions->[0]->[10]->[1], 1.1,   'Normalised count' );
is( $new_regions->[0]->[10]->[2], 2.1,   'Normalised count' );
is( $new_regions->[0]->[10]->[3], 4.6,   'Normalised count' );
is( $new_regions->[0]->[11],      undef, 'p value' );
is( $new_regions->[0]->[12],      undef, 'Adjusted p value' );
is( ( sprintf '%.2f', $new_regions->[0]->[13]->[0] ),
    1.18, 'Condition fold change' );
is( ( sprintf '%.3f', $new_regions->[0]->[13]->[1] ),
    0.233, 'Condition fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[0]->[0] ),
    0.46, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[0]->[1] ),
    -1.13, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[1]->[0] ),
    4.18, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[1]->[1] ),
    2.06, 'Group fold change' );
is( scalar keys %{ $new_regions->[0]->[15] }, 1, 'Gene annotation' );
is( $new_regions->[0]->[15]->{e61}->[0]->[0], 'ENSDARG00000095747', 'Gene ID' );
is( $new_regions->[0]->[15]->{e61}->[0]->[1], 'cxc64', 'Gene name' );
is(
    $new_regions->[0]->[15]->{e61}->[0]->[2],
    'CXC chemokine 64',
    'Gene description'
);
is( $new_regions->[0]->[15]->{e61}->[0]->[3], 'protein_coding', 'Gene type' );
is( $new_regions->[0]->[15]->{e61}->[0]->[4], 5, q{3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript ID' );
is( $new_regions->[0]->[15]->{e61}->[0]->[5]->[0]->[1],
    'protein_coding', 'Transcript type' );
is( $new_regions->[0]->[15]->{e61}->[1]->[0], 'ENSDARG00000024771', 'Gene ID' );
is( $new_regions->[0]->[15]->{e61}->[1]->[1], 'slc24a5', 'Gene name' );
is(
    $new_regions->[0]->[15]->{e61}->[1]->[2],
    'solute carrier family 24',
    'Gene description'
);
is( $new_regions->[0]->[15]->{e61}->[1]->[3], 'protein_coding', 'Gene type' );
is( $new_regions->[0]->[15]->{e61}->[1]->[4], 15, q{3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[1]->[5]->[0]->[0],
    'ENSDART00000033574', 'Transcript ID' );
is( $new_regions->[0]->[15]->{e61}->[1]->[5]->[0]->[1],
    'protein_coding', 'Transcript type' );

# Dump as table if multiple 3' end positions and genes with multiple distances

$regions = [
    [
        '1', 1, 110, 10, -10, '1',
        [ 110, 120 ],
        1,
        [ 10, 20 ],
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
                    [ 5, 15 ],
                    [ [ 'ENSDART00000133571', 'protein_coding' ] ]
                ],
                [
                    'ENSDARG00000024771',
                    'slc24a5',
                    'solute carrier family 24',
                    'protein_coding',
                    [ 5, 15 ],
                    [ [ 'ENSDART00000033574', 'protein_coding' ] ]
                ]
            ]
        }
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
        my $file     = $level . q{.} . $format;
        my $filepath = File::Spec->catfile( $tmp_dir, $file );
        ok( -e $filepath,  $file . ' exists' );
        ok( !-z $filepath, $file . ' is not empty' );
    }
}

# Parse output if multiple 3' end positions and genes with multiple distances

$new_regions = parse_table(
    {
        analysis     => $analysis,
        table_file   => File::Spec->catfile( $tmp_dir, 'all.tsv' ),
        table_format => 'tsv',
    }
);

is( $new_regions->[0]->[0],       '1',   'Region sequence name' );
is( $new_regions->[0]->[1],       1,     'Region start' );
is( $new_regions->[0]->[2],       110,   'Region end' );
is( $new_regions->[0]->[3],       undef, 'Region maximum read count' );
is( $new_regions->[0]->[4],       undef, 'Region log probability sum' );
is( $new_regions->[0]->[5],       '1',   q{3' end sequence name} );
is( $new_regions->[0]->[6]->[0],  110,   q{Nearest 3' end position} );
is( $new_regions->[0]->[6]->[1],  120,   q{Furthest 3' end position} );
is( $new_regions->[0]->[7],       1,     q{3' end strand} );
is( $new_regions->[0]->[8]->[0],  10,    q{Nearest 3' end read count} );
is( $new_regions->[0]->[8]->[1],  20,    q{Furthest 3' end read count} );
is( $new_regions->[0]->[9]->[0],  4,     'Count' );
is( $new_regions->[0]->[9]->[1],  1,     'Count' );
is( $new_regions->[0]->[9]->[2],  2,     'Count' );
is( $new_regions->[0]->[9]->[3],  7,     'Count' );
is( $new_regions->[0]->[10]->[0], 4.6,   'Normalised count' );
is( $new_regions->[0]->[10]->[1], 1.1,   'Normalised count' );
is( $new_regions->[0]->[10]->[2], 2.1,   'Normalised count' );
is( $new_regions->[0]->[10]->[3], 4.6,   'Normalised count' );
is( $new_regions->[0]->[11],      undef, 'p value' );
is( $new_regions->[0]->[12],      undef, 'Adjusted p value' );
is( ( sprintf '%.2f', $new_regions->[0]->[13]->[0] ),
    1.18, 'Condition fold change' );
is( ( sprintf '%.3f', $new_regions->[0]->[13]->[1] ),
    0.233, 'Condition fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[0]->[0] ),
    0.46, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[0]->[1] ),
    -1.13, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[1]->[0] ),
    4.18, 'Group fold change' );
is( ( sprintf '%.2f', $new_regions->[0]->[14]->[1]->[1] ),
    2.06, 'Group fold change' );
is( scalar keys %{ $new_regions->[0]->[15] }, 1, 'Gene annotation' );
is( $new_regions->[0]->[15]->{e61}->[0]->[0], 'ENSDARG00000095747', 'Gene ID' );
is( $new_regions->[0]->[15]->{e61}->[0]->[1], 'cxc64', 'Gene name' );
is(
    $new_regions->[0]->[15]->{e61}->[0]->[2],
    'CXC chemokine 64',
    'Gene description'
);
is( $new_regions->[0]->[15]->{e61}->[0]->[3], 'protein_coding', 'Gene type' );
is( $new_regions->[0]->[15]->{e61}->[0]->[4]->[0],
    5, q{Nearest 3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[0]->[4]->[1],
    15, q{Furthest 3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript ID' );
is( $new_regions->[0]->[15]->{e61}->[0]->[5]->[0]->[1],
    'protein_coding', 'Transcript type' );
is( $new_regions->[0]->[15]->{e61}->[1]->[0], 'ENSDARG00000024771', 'Gene ID' );
is( $new_regions->[0]->[15]->{e61}->[1]->[1], 'slc24a5', 'Gene name' );
is(
    $new_regions->[0]->[15]->{e61}->[1]->[2],
    'solute carrier family 24',
    'Gene description'
);
is( $new_regions->[0]->[15]->{e61}->[1]->[3], 'protein_coding', 'Gene type' );
is( $new_regions->[0]->[15]->{e61}->[1]->[4]->[0],
    5, q{Nearest 3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[1]->[4]->[1],
    15, q{Furthest 3' end distance} );
is( $new_regions->[0]->[15]->{e61}->[1]->[5]->[0]->[0],
    'ENSDART00000033574', 'Transcript ID' );
is( $new_regions->[0]->[15]->{e61}->[1]->[5]->[0]->[1],
    'protein_coding', 'Transcript type' );

# TODO: Actually test output
