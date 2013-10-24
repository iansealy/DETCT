use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 48;

use DETCT::Transcript;

my $transcript = DETCT::Transcript->new(
    {
        stable_id => 'ENSDART00000133571',
        biotype   => 'protein_coding',
        seq_name  => '5',
        start     => 40352744,
        end       => 40354399,
        strand    => 1,
    }
);

isa_ok( $transcript, 'DETCT::Transcript' );

# Test stable id attribute
is( $transcript->stable_id, 'ENSDART00000133571', 'Get stable id' );
is( $transcript->set_stable_id('ENSDART00000033574'), undef, 'Set stable id' );
is( $transcript->stable_id, 'ENSDART00000033574', 'Get new stable id' );
throws_ok { $transcript->set_stable_id() } qr/No stable id specified/ms,
  'No stable id';
throws_ok { $transcript->set_stable_id('#invalid#') } qr/Invalid stable id/ms,
  'Invalid stable id';

# Test name attribute
is( $transcript->name,                  undef,       'Get name' );
is( $transcript->set_name('cxc64-001'), undef,       'Set name' );
is( $transcript->name,                  'cxc64-001', 'Get new name' );
is( $transcript->set_name(),            undef,       'Set undef name' );
is( $transcript->name,                  undef,       'Get undef name' );
my $long_name = 'X' x ( $DETCT::Transcript::MAX_NAME_LENGTH + 1 );
throws_ok { $transcript->set_name('') } qr/Name is empty/ms, 'Empty name';
throws_ok { $transcript->set_name($long_name) }
qr/longer than \d+ characters/ms, 'Invalid name';

# Test description attribute
is( $transcript->description, undef, 'Get description' );
is( $transcript->set_description('CXC chemokine 64'), undef,
    'Set description' );
is( $transcript->description, 'CXC chemokine 64', 'Get new description' );
is( $transcript->set_description(), undef, 'Set undef description' );
is( $transcript->description,       undef, 'Get undef description' );

# Test biotype attribute
is( $transcript->biotype, 'protein_coding', 'Get biotype' );
is( $transcript->set_biotype('nonsense_mediated_decay'), undef, 'Set biotype' );
is( $transcript->biotype, 'nonsense_mediated_decay', 'Get new biotype' );
throws_ok { $transcript->set_biotype() } qr/No biotype specified/ms,
  'No biotype';
throws_ok { $transcript->set_biotype('#invalid#') } qr/Invalid biotype/ms,
  'Invalid biotype';

# Test sequence name attribute
is( $transcript->seq_name,          '5',   'Get sequence name' );
is( $transcript->set_seq_name('6'), undef, 'Set sequence name' );
is( $transcript->seq_name,          '6',   'Get new sequence name' );
is( $transcript->set_seq_name('GL456210.1'),
    undef, 'Set sequence name with full stop' );
throws_ok { $transcript->set_seq_name() } qr/No sequence name specified/ms,
  'No sequence name';
throws_ok { $transcript->set_seq_name('#invalid#') }
qr/Invalid sequence name/ms, 'Invalid sequence name';

# Test start attribute
is( $transcript->start,               40352744, 'Get start' );
is( $transcript->set_start(30352744), undef,    'Set start' );
is( $transcript->start,               30352744, 'Get new start' );
throws_ok { $transcript->set_start() } qr/No start specified/ms, 'No start';
throws_ok { $transcript->set_start(-1) } qr/Invalid start/ms, 'Invalid start';

# Test end attribute
is( $transcript->end,               40354399, 'Get end' );
is( $transcript->set_end(30354399), undef,    'Set end' );
is( $transcript->end,               30354399, 'Get new end' );
throws_ok { $transcript->set_end() } qr/No end specified/ms, 'No end';
throws_ok { $transcript->set_end(-2) } qr/Invalid end/ms, 'Invalid end';

# Test strand attribute
is( $transcript->strand,         1,     'Get strand' );
is( $transcript->set_strand(-1), undef, 'Set strand' );
is( $transcript->strand,         -1,    'Get new strand' );
throws_ok { $transcript->set_strand() } qr/No strand specified/ms, 'No strand';
throws_ok { $transcript->set_strand(0) } qr/Invalid strand/ms, 'Invalid strand';

# Mock gene object
my $gene = Test::MockObject->new();
$gene->set_isa('DETCT::Gene');
$gene->set_always( 'name', 'cxc64' );

# Test gene attribute
is( $transcript->gene,            undef,   'Get gene' );
is( $transcript->set_gene($gene), undef,   'Set gene' );
is( $transcript->gene->name,      'cxc64', 'Get new gene' );
throws_ok { $transcript->set_gene('invalid') } qr/Class of gene/ms,
  'Invalid gene';
