use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 33;

use DETCT::GeneOntologyTerm;

my $term = DETCT::GeneOntologyTerm->new(
    {
        accession  => 'GO:0005622',
        namespace  => 'cellular_component',
        name       => 'intracellular',
        definition => 'The living contents of a cell',
    }
);

isa_ok( $term, 'DETCT::GeneOntologyTerm' );

# Test accession attribute
is( $term->accession,                   'GO:0005622', 'Get accession' );
is( $term->set_accession('GO:0016874'), undef,        'Set accession' );
is( $term->accession,                   'GO:0016874', 'Get new accession' );
throws_ok { $term->set_accession() } qr/No accession specified/ms,
  'No accession';
throws_ok { $term->set_accession('invalid') } qr/Invalid accession/ms,
  'Invalid accession';

# Test namespace attribute
is( $term->namespace, 'cellular_component', 'Get namespace' );
is( $term->set_namespace('molecular_function'), undef, 'Set namespace' );
is( $term->namespace, 'molecular_function', 'Get new namespace' );
throws_ok { $term->set_namespace() } qr/No namespace specified/ms,
  'No namespace';
throws_ok { $term->set_namespace('invalid') } qr/Invalid namespace/ms,
  'Invalid namespace';

# Test name attribute
is( $term->name,                        'intracellular',   'Get name' );
is( $term->set_name('ligase activity'), undef,             'Set name' );
is( $term->name,                        'ligase activity', 'Get new name' );
throws_ok { $term->set_name() } qr/No name specified/ms, 'No name';

# Test definition attribute
is( $term->definition, 'The living contents of a cell', 'Get definition' );
is( $term->set_definition('Catalysis of the joining of two substances'),
    undef, 'Set definition' );
is(
    $term->definition,
    'Catalysis of the joining of two substances',
    'Get new definition'
);
throws_ok { $term->set_definition() } qr/No definition specified/ms,
  'No definition';

# Mock gene objects
my $gene1 = Test::MockObject->new();
$gene1->set_isa('DETCT::Gene');
$gene1->set_always( 'stable_id', 'ENSDARG00000095747' );
my $gene2 = Test::MockObject->new();
$gene2->set_isa('DETCT::Gene');
$gene2->set_always( 'stable_id', 'ENSDARG00000024771' );

# Test adding and retrieving genes
my $genes;
$genes = $term->get_all_genes();
is( scalar @{$genes},        0,     'No genes' );
is( $term->add_gene($gene1), undef, 'Add gene' );
$genes = $term->get_all_genes();
is( scalar @{$genes}, 1, 'Get one gene' );
$term->add_gene($gene2);
is( scalar @{$genes}, 2, 'Get two genes' );
throws_ok { $term->add_gene() } qr/No gene specified/ms, 'No gene specified';
throws_ok { $term->add_gene('invalid') } qr/Class of gene/ms, 'Invalid gene';

# Test adding and retrieving evidence codes
is( $term->get_evidence_code($gene1),          undef, 'Get evidence code' );
is( $term->add_evidence_code( 'EXP', $gene1 ), undef, 'Set evidence code' );
is( $term->get_evidence_code($gene1), 'EXP', 'Get new evidence code by term' );
is( $term->get_evidence_code( $gene1->stable_id ),
    'EXP', 'Get new evidence code by stable_id' );
throws_ok { $term->add_evidence_code( 'XXX', $gene1 ) }
qr/Invalid evidence code/ms, 'Invalid evidence code';
throws_ok { $term->get_evidence_code() } qr/No gene specified/ms, 'No gene';
throws_ok { $term->get_evidence_code( [] ) } qr/Class of gene/ms,
  'Invalid gene';
throws_ok { $term->get_evidence_code('invalid') } qr/Invalid stable id/ms,
  'Invalid stable id';
