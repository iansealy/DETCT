use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 117;

use DETCT::GeneFinder;

# Mock genes
my @genes;
my @transcript_three_prime_ends =
  ( [ 100, 1 ], [ 200, 1 ], [ 300, -1 ], [ 400, -1 ], );
foreach my $transcript_three_prime_end (@transcript_three_prime_ends) {
    my ( $pos, $strand ) = @{$transcript_three_prime_end};

    # Construct start and end so $pos is always 3' end
    my $start = $strand == 1 ? $pos - 50 : $pos;
    my $end   = $strand == 1 ? $pos      : $pos + 50;

    # Create genes named by 3' end position
    my $gene = Test::MockObject->new();
    $gene->set_always( 'stable_id',         'ENSDARG00000095747' );
    $gene->set_always( 'external_name',     q{g} . $pos . q{:} . $strand );
    $gene->set_always( 'description',       undef );
    $gene->set_always( 'biotype',           'protein_coding' );
    $gene->set_always( 'seq_region_start',  $start );
    $gene->set_always( 'seq_region_end',    $end );
    $gene->set_always( 'seq_region_strand', $strand );
    my $attribute1 = Test::MockObject->new();
    $attribute1->set_always( 'code', 'appris_pi' );
    my $attribute2 = Test::MockObject->new();
    $attribute2->set_always( 'code', 'gencode_basic' );
    my $attribute3 = Test::MockObject->new();
    $attribute3->set_always( 'code',  'TSL' );
    $attribute3->set_always( 'value', 'tsl3 (assigned to version 2)' );
    my $attribute4 = Test::MockObject->new();
    $attribute4->set_always( 'code', 'synonym' );
    my $transcript = Test::MockObject->new();
    $transcript->set_always( 'get_all_Attributes',
        [ $attribute1, $attribute2, $attribute3, $attribute4 ] );
    $transcript->set_always( 'stable_id',        'ENSDART00000133571' );
    $transcript->set_always( 'external_name',    q{t} . $pos . q{:} . $strand );
    $transcript->set_always( 'description',      undef );
    $transcript->set_always( 'biotype',          'protein_coding' );
    $transcript->set_always( 'seq_region_start', $start );
    $transcript->set_always( 'seq_region_end',   $end );
    $transcript->set_always( 'seq_region_strand', $strand );
    my $transcript_far = Test::MockObject->new();
    $transcript_far->set_always( 'get_all_Attributes', undef );
    $transcript_far->set_always( 'stable_id',          'ENSDART00000133572' );
    $transcript_far->set_always( 'external_name',      'cxc64-001' );
    $transcript_far->set_always( 'description',        undef );
    $transcript_far->set_always( 'biotype',            'protein_coding' );
    $transcript_far->set_always( 'seq_region_start',   100_000 );
    $transcript_far->set_always( 'seq_region_end',     100_100 );
    $transcript_far->set_always( 'seq_region_strand',  $strand );
    $gene->set_always( 'get_all_Transcripts',
        [ $transcript, $transcript_far ] );
    push @genes, $gene;
}

# Mock slice
my $slice = Test::MockObject->new();
$slice->set_always( 'get_all_Genes', \@genes );

# Mock slice adaptor
my $slice_adaptor = Test::MockObject->new();
$slice_adaptor->set_isa('Bio::EnsEMBL::DBSQL::SliceAdaptor');
$slice_adaptor->set_always( 'fetch_by_region', $slice );

my $gene_finder =
  DETCT::GeneFinder->new( { slice_adaptor => $slice_adaptor, } );

isa_ok( $gene_finder, 'DETCT::GeneFinder' );

# Test Ensembl slice adaptor attribute
isa_ok( $gene_finder->slice_adaptor, 'Bio::EnsEMBL::DBSQL::SliceAdaptor' );
throws_ok { $gene_finder->set_slice_adaptor() }
qr/No Ensembl slice adaptor specified/ms, 'No Ensembl slice adaptor';
throws_ok { $gene_finder->set_slice_adaptor('invalid') }
qr/Class of Ensembl slice adaptor/ms, 'Invalid Ensembl slice adaptor';

# Test skip transcripts
is( $gene_finder->is_skip_transcript('ENSDART00000135768'),
    0, 'ENSDART00000135768 is not skip transcript' );
is( $gene_finder->set_skip_transcripts( ['ENSDART00000135768'] ),
    undef, 'Add skip transcripts' );
is( $gene_finder->is_skip_transcript('ENSDART00000135768'),
    1, 'ENSDART00000135768 is skip transcript' );

my $genes;
my $transcripts;
my $distance;
my $nearest_end_pos;

# Near to one gene on forward strand
( $genes, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_genes( '1', 110, 1 );
is( $genes->[0]->name, 'g100:1',
    q{Gene with 3' end at 100 bp on forward strand} );
is( $distance,        10,  q{3' end is 10 bp downstream} );
is( $nearest_end_pos, 100, q{3' end at 100 bp} );
( $genes, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_genes( '1', 90, 1 );
is( $genes->[0]->name, 'g100:1',
    q{Gene with 3' end at 100 bp on forward strand} );
is( $distance,        -10, q{3' end is 10 bp upstream} );
is( $nearest_end_pos, 100, q{3' end at 100 bp} );

# Near to one gene on reverse strand
( $genes, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_genes( '1', 290, -1 );
is( $genes->[0]->name, 'g300:-1',
    q{Gene with 3' end at 300 bp on reverse strand} );
is( $distance,        10,  q{3' end is 10 bp downstream} );
is( $nearest_end_pos, 300, q{3' end at 300 bp} );
( $genes, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_genes( '1', 310, -1 );
is( $genes->[0]->name, 'g300:-1',
    q{Gene with 3' end at 300 bp on reverse strand} );
is( $distance,        -10, q{3' end is 10 bp upstream} );
is( $nearest_end_pos, 300, q{3' end at 300 bp} );

# Between two genes on forward strand
( $genes, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_genes( '1', 150, 1 );
is( $genes->[0]->name, 'g100:1',
    q{Gene with 3' end at 100 bp on forward strand} );
is( $distance,        50,  q{3' end is 50 bp upstream} );
is( $nearest_end_pos, 100, q{3' end at 100 bp} );

# Between two genes on reverse strand
( $genes, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_genes( '1', 350, -1 );
is( $genes->[0]->name, 'g400:-1',
    q{Gene with 3' end at 400 bp on reverse strand} );
is( $distance,        50,  q{3' end is 50 bp upstream} );
is( $nearest_end_pos, 400, q{3' end at 400 bp} );

# Near to one transcript on forward strand
( $transcripts, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_transcripts( '1', 110, 1 );
is( $transcripts->[0]->name,
    't100:1', q{Transcript with 3' end at 100 bp on forward strand} );
is( $distance,        10,  q{3' end is 10 bp downstream} );
is( $nearest_end_pos, 100, q{3' end at 100 bp} );
( $transcripts, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_transcripts( '1', 90, 1 );
is( $transcripts->[0]->name,
    't100:1', q{Transcript with 3' end at 100 bp on forward strand} );
is( $distance,        -10, q{3' end is 10 bp upstream} );
is( $nearest_end_pos, 100, q{3' end at 100 bp} );

# Near to one transcript on reverse strand
( $transcripts, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_transcripts( '1', 290, -1 );
is( $transcripts->[0]->name,
    't300:-1', q{Transcript with 3' end at 300 bp on reverse strand} );
is( $distance,        10,  q{3' end is 10 bp downstream} );
is( $nearest_end_pos, 300, q{3' end at 300 bp} );
( $transcripts, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_transcripts( '1', 310, -1 );
is( $transcripts->[0]->name,
    't300:-1', q{Transcript with 3' end at 300 bp on reverse strand} );
is( $distance,        -10, q{3' end is 10 bp upstream} );
is( $nearest_end_pos, 300, q{3' end at 300 bp} );

# Between two transcripts on forward strand
( $transcripts, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_transcripts( '1', 150, 1 );
is( $transcripts->[0]->name,
    't100:1', q{Transcript with 3' end at 100 bp on forward strand} );
is( $distance,        50,  q{3' end is 50 bp upstream} );
is( $nearest_end_pos, 100, q{3' end at 100 bp} );

# Between two transcripts on reverse strand
( $transcripts, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_transcripts( '1', 350, -1 );
is( $transcripts->[0]->name,
    't400:-1', q{Transcript with 3' end at 400 bp on reverse strand} );
is( $distance,        50,  q{3' end is 50 bp upstream} );
is( $nearest_end_pos, 400, q{3' end at 400 bp} );

# Check adding gene annotation required parameters
throws_ok { $gene_finder->add_gene_annotation() } qr/No regions specified/ms,
  'No regions';

my $regions;
my $annotated_regions;
my $gv;

# Adding gene annotation
$regions = [
    [ '1', 1, 1000, 10, -10, '1', 110, 1,  10, [], undef, undef, [], [] ],
    [ '1', 1, 1000, 10, -10, '1', 290, -1, 10, [], undef, undef, [], [] ],
    [ '1', 1, 1000, 10, -10, '1', 100, 1,  10, [], undef, undef, [], [] ],
    [ '1', 1, 1000, 10, -10, '1', 300, -1, 10, [], undef, undef, [], [] ],
];
$annotated_regions = $gene_finder->add_gene_annotation($regions);
($gv) = keys %{ $annotated_regions->[0]->[-1] };    # Genebuild version varies
is( scalar keys %{ $annotated_regions->[0]->[-1] },   1, '1 genebuild' );
is( scalar @{ $annotated_regions->[0]->[-1]->{$gv} }, 1, '1 gene' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[0],
    'ENSDARG00000095747', 'Stable id' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[1],
    'g100:1', q{3' end as name} );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[2], undef, 'Description' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[3],
    'protein_coding', 'Biotype' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[4], 10, 'Distance downstream' );
is( scalar @{ $annotated_regions->[0]->[-1]->{$gv}->[0]->[5] },
    1, '1 transcript' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript stable id' );
is(
    $annotated_regions->[0]->[-1]->{$gv}->[0]->[5]->[0]->[1],
    'protein_coding:appris_pi:gencode_basic:tsl3',
    'Transcript biotype'
);
is( scalar keys %{ $annotated_regions->[1]->[-1] },   1, '1 genebuild' );
is( scalar @{ $annotated_regions->[1]->[-1]->{$gv} }, 1, '1 gene' );
is( $annotated_regions->[1]->[-1]->{$gv}->[0]->[0],
    'ENSDARG00000095747', 'Stable id' );
is( $annotated_regions->[1]->[-1]->{$gv}->[0]->[1],
    'g300:-1', q{3' end as name} );
is( $annotated_regions->[1]->[-1]->{$gv}->[0]->[2], undef, 'Description' );
is( $annotated_regions->[1]->[-1]->{$gv}->[0]->[3],
    'protein_coding', 'Biotype' );
is( $annotated_regions->[1]->[-1]->{$gv}->[0]->[4], 10, 'Distance downstream' );
is( scalar @{ $annotated_regions->[1]->[-1]->{$gv}->[0]->[5] },
    1, '1 transcript' );
is( $annotated_regions->[1]->[-1]->{$gv}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript stable id' );
is(
    $annotated_regions->[1]->[-1]->{$gv}->[0]->[5]->[0]->[1],
    'protein_coding:appris_pi:gencode_basic:tsl3',
    'Transcript biotype'
);
is( scalar keys %{ $annotated_regions->[2]->[-1] },   1, '1 genebuild' );
is( scalar @{ $annotated_regions->[2]->[-1]->{$gv} }, 1, '1 gene' );
is( $annotated_regions->[2]->[-1]->{$gv}->[0]->[0],
    'ENSDARG00000095747', 'Stable id' );
is( $annotated_regions->[2]->[-1]->{$gv}->[0]->[1],
    'g100:1', q{3' end as name} );
is( $annotated_regions->[2]->[-1]->{$gv}->[0]->[2], undef, 'Description' );
is( $annotated_regions->[2]->[-1]->{$gv}->[0]->[3],
    'protein_coding', 'Biotype' );
is( $annotated_regions->[2]->[-1]->{$gv}->[0]->[4], 0, 'Distance downstream' );
is( scalar @{ $annotated_regions->[2]->[-1]->{$gv}->[0]->[5] },
    1, '1 transcript' );
is( $annotated_regions->[2]->[-1]->{$gv}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript stable id' );
is(
    $annotated_regions->[2]->[-1]->{$gv}->[0]->[5]->[0]->[1],
    'protein_coding:appris_pi:gencode_basic:tsl3',
    'Transcript biotype'
);
is( scalar keys %{ $annotated_regions->[3]->[-1] },   1, '1 genebuild' );
is( scalar @{ $annotated_regions->[3]->[-1]->{$gv} }, 1, '1 gene' );
is( $annotated_regions->[3]->[-1]->{$gv}->[0]->[0],
    'ENSDARG00000095747', 'Stable id' );
is( $annotated_regions->[3]->[-1]->{$gv}->[0]->[1],
    'g300:-1', q{3' end as name} );
is( $annotated_regions->[3]->[-1]->{$gv}->[0]->[2], undef, 'Description' );
is( $annotated_regions->[3]->[-1]->{$gv}->[0]->[3],
    'protein_coding', 'Biotype' );
is( $annotated_regions->[3]->[-1]->{$gv}->[0]->[4], 0, 'Distance downstream' );
is( scalar @{ $annotated_regions->[3]->[-1]->{$gv}->[0]->[5] },
    1, '1 transcript' );
is( $annotated_regions->[3]->[-1]->{$gv}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript stable id' );
is(
    $annotated_regions->[3]->[-1]->{$gv}->[0]->[5]->[0]->[1],
    'protein_coding:appris_pi:gencode_basic:tsl3',
    'Transcript biotype'
);

# Adding gene annotation if no 3' end
$regions = [
    [
        '1', 1, 1000, 10, -10, undef, undef, undef,
        undef, [], undef, undef, [], []
    ],
];
$annotated_regions = $gene_finder->add_gene_annotation($regions);
is( scalar keys %{ $annotated_regions->[0]->[-1] }, 0, '0 genebuilds' );

# Adding gene annotation if multiple 3' end positions for single region
$regions = [
    [
        '1', 1, 1000, 10, -10, '1', [ 110, 120 ],
        1,                          [ 10,  20 ],
        [], undef, undef, [], []
    ],
];
$annotated_regions = $gene_finder->add_gene_annotation($regions);
($gv) = keys %{ $annotated_regions->[0]->[-1] };    # Genebuild version varies
is( scalar keys %{ $annotated_regions->[0]->[-1] },   1, '1 genebuild' );
is( scalar @{ $annotated_regions->[0]->[-1]->{$gv} }, 1, '1 gene' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[0],
    'ENSDARG00000095747', 'Stable id' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[1],
    'g100:1', q{3' end as name} );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[2], undef, 'Description' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[3],
    'protein_coding', 'Biotype' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[4]->[0],
    10, 'Nearest distance downstream' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[4]->[1],
    20, 'Furthest distance downstream' );
is( scalar @{ $annotated_regions->[0]->[-1]->{$gv}->[0]->[5] },
    1, '1 transcript' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript stable id' );
is(
    $annotated_regions->[0]->[-1]->{$gv}->[0]->[5]->[0]->[1],
    'protein_coding:appris_pi:gencode_basic:tsl3',
    'Transcript biotype'
);

# Check skipping transcripts
$gene_finder = DETCT::GeneFinder->new(
    {
        slice_adaptor    => $slice_adaptor,
        skip_transcripts => ['ENSDART00000133571'],
    }
);
isa_ok( $gene_finder, 'DETCT::GeneFinder' );
( $genes, $distance, $nearest_end_pos ) =
  $gene_finder->get_nearest_genes( '1', 100, 1 );
is( $distance,        -100_000, q{3' end is 100,000 bp upstream} );
is( $nearest_end_pos, 100_100,  q{3' end at 100,100 bp} );

# Mock genes all on one strand
@genes = ();
@transcript_three_prime_ends =
  ( [ 100, 1 ], [ 200, 1 ], [ 300, 1 ], [ 400, 1 ], );
foreach my $transcript_three_prime_end (@transcript_three_prime_ends) {
    my ( $pos, $strand ) = @{$transcript_three_prime_end};

    # Create genes named by 3' end position
    my $gene = Test::MockObject->new();
    $gene->set_always( 'stable_id',         'ENSDARG00000095747' );
    $gene->set_always( 'external_name',     q{g} . $pos . q{:} . $strand );
    $gene->set_always( 'description',       undef );
    $gene->set_always( 'biotype',           'protein_coding' );
    $gene->set_always( 'seq_region_start',  1 );
    $gene->set_always( 'seq_region_end',    $pos );
    $gene->set_always( 'seq_region_strand', $strand );
    my $transcript = Test::MockObject->new();
    $transcript->set_always( 'get_all_Attributes', undef );
    $transcript->set_always( 'stable_id',          'ENSDART00000133571' );
    $transcript->set_always( 'external_name',      'cxc64-001' );
    $transcript->set_always( 'description',        undef );
    $transcript->set_always( 'biotype',            'protein_coding' );
    $transcript->set_always( 'seq_region_start',   $pos - 50 );
    $transcript->set_always( 'seq_region_end',     $pos );
    $transcript->set_always( 'seq_region_strand',  $strand );
    $gene->set_always( 'get_all_Transcripts', [$transcript] );
    push @genes, $gene;
}

# Mock slice
$slice = Test::MockObject->new();
$slice->set_always( 'get_all_Genes', \@genes );

# Mock slice adaptor
$slice_adaptor = Test::MockObject->new();
$slice_adaptor->set_isa('Bio::EnsEMBL::DBSQL::SliceAdaptor');
$slice_adaptor->set_always( 'fetch_by_region', $slice );

$gene_finder = DETCT::GeneFinder->new( { slice_adaptor => $slice_adaptor, } );

isa_ok( $gene_finder, 'DETCT::GeneFinder' );

# Adding gene annotation with genes only on one strand
$regions = [
    [ '1', 1, 1000, 10, -10, '1', 110, 1,  10, [], undef, undef, [], [] ],
    [ '1', 1, 1000, 10, -10, '1', 290, -1, 10, [], undef, undef, [], [] ],
];
$annotated_regions = $gene_finder->add_gene_annotation($regions);
($gv) = keys %{ $annotated_regions->[0]->[-1] };    # Genebuild version varies
is( scalar keys %{ $annotated_regions->[0]->[-1] },   1, '1 genebuild' );
is( scalar @{ $annotated_regions->[0]->[-1]->{$gv} }, 1, '1 gene' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[0],
    'ENSDARG00000095747', 'Stable id' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[1],
    'g100:1', q{3' end as name} );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[2], undef, 'Description' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[3],
    'protein_coding', 'Biotype' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[4], 10, 'Distance downstream' );
is( scalar @{ $annotated_regions->[0]->[-1]->{$gv}->[0]->[5] },
    1, '1 transcript' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[5]->[0]->[0],
    'ENSDART00000133571', 'Transcript stable id' );
is( $annotated_regions->[0]->[-1]->{$gv}->[0]->[5]->[0]->[1],
    'protein_coding:::', 'Transcript biotype' );
is( scalar keys %{ $annotated_regions->[1]->[-1] },
    0, 'No genes on reverse strand' );

# Mock slice adaptor with non-existent slice
$slice_adaptor = Test::MockObject->new();
$slice_adaptor->set_isa('Bio::EnsEMBL::DBSQL::SliceAdaptor');
$slice_adaptor->set_always( 'fetch_by_region', undef );

$gene_finder = DETCT::GeneFinder->new( { slice_adaptor => $slice_adaptor, } );

isa_ok( $gene_finder, 'DETCT::GeneFinder' );

$annotated_regions = $gene_finder->add_gene_annotation($regions);

is( scalar keys %{ $annotated_regions->[0]->[-1] },
    0, 'No genes on reverse strand' );

# Check Ensembl database types

# Mock slice
$slice = Test::MockObject->new();
$slice->set_always( 'get_all_Genes', \@genes );

# Mock slice adaptor
$slice_adaptor = Test::MockObject->new();
$slice_adaptor->set_isa('Bio::EnsEMBL::DBSQL::SliceAdaptor');
$slice_adaptor->set_always( 'fetch_by_region', $slice );

my $args;

$gene_finder = DETCT::GeneFinder->new( { slice_adaptor => $slice_adaptor, } );
$annotated_regions = $gene_finder->add_gene_annotation($regions);
( undef, $args ) = $slice->next_call;
is( $args->[2], undef, 'No Ensembl database type' );

$gene_finder = DETCT::GeneFinder->new(
    {
        slice_adaptor    => $slice_adaptor,
        ensembl_db_types => ['core'],
    }
);
$annotated_regions = $gene_finder->add_gene_annotation($regions);
( undef, $args ) = $slice->next_call;
is( $args->[2], 'core', 'Core Ensembl database type' );

$gene_finder = DETCT::GeneFinder->new(
    {
        slice_adaptor    => $slice_adaptor,
        ensembl_db_types => [ 'core', 'otherfeatures' ],
    }
);
$annotated_regions = $gene_finder->add_gene_annotation($regions);
( undef, $args ) = $slice->next_call;
is( $args->[2], 'core', 'Core Ensembl database type' );
( undef, $args ) = $slice->next_call;
is( $args->[2], 'otherfeatures', 'OtherFeatures Ensembl database type' );
