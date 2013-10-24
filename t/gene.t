use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 55;

use DETCT::Gene;

my $gene = DETCT::Gene->new(
    {
        genebuild_version => 'e61',
        stable_id         => 'ENSDARG00000095747',
        biotype           => 'protein_coding',
        seq_name          => '5',
        start             => 40352744,
        end               => 40354399,
        strand            => 1,
    }
);

isa_ok( $gene, 'DETCT::Gene' );

# Test genebuild version attribute
is( $gene->genebuild_version,            'e61', 'Get genebuild version' );
is( $gene->set_genebuild_version('e62'), undef, 'Set genebuild version' );
is( $gene->genebuild_version,            'e62', 'Get new genebuild version' );
throws_ok { $gene->set_genebuild_version() }
qr/No genebuild version specified/ms, 'No genebuild version';
throws_ok { $gene->set_genebuild_version('#invalid#') }
qr/Invalid genebuild version/ms, 'Invalid genebuild version';

# Test stable id attribute
is( $gene->stable_id, 'ENSDARG00000095747', 'Get stable id' );
is( $gene->set_stable_id('ENSDARG00000024771'), undef, 'Set stable id' );
is( $gene->stable_id, 'ENSDARG00000024771', 'Get new stable id' );
throws_ok { $gene->set_stable_id() } qr/No stable id specified/ms,
  'No stable id';
throws_ok { $gene->set_stable_id('#invalid#') } qr/Invalid stable id/ms,
  'Invalid stable id';

# Test name attribute
is( $gene->name,              undef,   'Get name' );
is( $gene->set_name('cxc64'), undef,   'Set name' );
is( $gene->name,              'cxc64', 'Get new name' );
is( $gene->set_name(),        undef,   'Set undef name' );
is( $gene->name,              undef,   'Get undef name' );
my $long_name = 'X' x ( $DETCT::Gene::MAX_NAME_LENGTH + 1 );
throws_ok { $gene->set_name('') } qr/Name is empty/ms, 'Empty name';
throws_ok { $gene->set_name($long_name) } qr/longer than \d+ characters/ms,
  'Invalid name';

# Test description attribute
is( $gene->description,                         undef, 'Get description' );
is( $gene->set_description('CXC chemokine 64'), undef, 'Set description' );
is( $gene->description,       'CXC chemokine 64', 'Get new description' );
is( $gene->set_description(), undef,              'Set undef description' );
is( $gene->description,       undef,              'Get undef description' );

# Test biotype attribute
is( $gene->biotype, 'protein_coding', 'Get biotype' );
is( $gene->set_biotype('nonsense_mediated_decay'), undef, 'Set biotype' );
is( $gene->biotype, 'nonsense_mediated_decay', 'Get new biotype' );
throws_ok { $gene->set_biotype() } qr/No biotype specified/ms, 'No biotype';
throws_ok { $gene->set_biotype('#invalid#') } qr/Invalid biotype/ms,
  'Invalid biotype';

# Test sequence name attribute
is( $gene->seq_name,          '5',   'Get sequence name' );
is( $gene->set_seq_name('6'), undef, 'Set sequence name' );
is( $gene->seq_name,          '6',   'Get new sequence name' );
is( $gene->set_seq_name('GL456210.1'),
    undef, 'Set sequence name with full stop' );
throws_ok { $gene->set_seq_name() } qr/No sequence name specified/ms,
  'No sequence name';
throws_ok { $gene->set_seq_name('#invalid#') } qr/Invalid sequence name/ms,
  'Invalid sequence name';

# Test start attribute
is( $gene->start,               40352744, 'Get start' );
is( $gene->set_start(30352744), undef,    'Set start' );
is( $gene->start,               30352744, 'Get new start' );
throws_ok { $gene->set_start() } qr/No start specified/ms, 'No start';
throws_ok { $gene->set_start(-1) } qr/Invalid start/ms, 'Invalid start';

# Test end attribute
is( $gene->end,               40354399, 'Get end' );
is( $gene->set_end(30354399), undef,    'Set end' );
is( $gene->end,               30354399, 'Get new end' );
throws_ok { $gene->set_end() } qr/No end specified/ms, 'No end';
throws_ok { $gene->set_end(-2) } qr/Invalid end/ms, 'Invalid end';

# Test strand attribute
is( $gene->strand,         1,     'Get strand' );
is( $gene->set_strand(-1), undef, 'Set strand' );
is( $gene->strand,         -1,    'Get new strand' );
throws_ok { $gene->set_strand() } qr/No strand specified/ms, 'No strand';
throws_ok { $gene->set_strand(0) } qr/Invalid strand/ms, 'Invalid strand';

# Mock transcript objects
my $transcript1 = Test::MockObject->new();
$transcript1->set_isa('DETCT::Transcript');
my $transcript2 = Test::MockObject->new();
$transcript2->set_isa('DETCT::Transcript');

# Test adding and retrieving transcripts
my $transcripts;
$transcripts = $gene->get_all_transcripts();
is( scalar @{$transcripts},              0,     'No transcripts' );
is( $gene->add_transcript($transcript1), undef, 'Add transcript' );
$transcripts = $gene->get_all_transcripts();
is( scalar @{$transcripts}, 1, 'Get one transcript' );
$gene->add_transcript($transcript2);
is( scalar @{$transcripts}, 2, 'Get two transcripts' );
throws_ok { $gene->add_transcript() } qr/No transcript specified/ms,
  'No transcript specified';
throws_ok { $gene->add_transcript('invalid') } qr/Class of transcript/ms,
  'Invalid transcript';

