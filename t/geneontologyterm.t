use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 19;

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
