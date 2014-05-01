use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 36;

use DETCT::Sample;

my $sample = DETCT::Sample->new(
    {
        name      => 'zmp_ph1_1m',
        condition => 'mutant',
        group     => '1',
        tag       => 'NNNNBGAGGC',
        bam_file  => 't/data/test1.bam',
    }
);

isa_ok( $sample, 'DETCT::Sample' );

# Test name attribute
is( $sample->name,                   'zmp_ph1_1m', 'Get name' );
is( $sample->set_name('zmp_ph1_1s'), undef,        'Set name' );
is( $sample->name,                   'zmp_ph1_1s', 'Get new name' );
throws_ok { $sample->set_name() } qr/No name specified/ms, 'No name';
my $long_name = 'X' x ( $DETCT::Sample::MAX_NAME_LENGTH + 1 );
throws_ok { $sample->set_name(' ') } qr/Invalid name specified/ms,
  'Invalid name';
throws_ok { $sample->set_name('') } qr/Empty name specified/ms, 'Empty name';
throws_ok { $sample->set_name($long_name) } qr/longer than \d+ characters/ms,
  'Long name';

# Test description attribute
is( $sample->description, undef, 'Get description' );
is( $sample->set_description('ZMP phenotype 1.1 mutant'),
    undef, 'Set description' );
is( $sample->description, 'ZMP phenotype 1.1 mutant', 'Get new description' );
is( $sample->set_description(), undef, 'Set undef description' );
is( $sample->description,       undef, 'Get undef description' );

# Test condition attribute
is( $sample->condition,                'mutant',  'Get condition' );
is( $sample->set_condition('sibling'), undef,     'Set condition' );
is( $sample->condition,                'sibling', 'Get new condition' );
throws_ok { $sample->set_condition() } qr/No condition specified/ms,
  'No condition';
my $long_condition = 'X' x ( $DETCT::Sample::MAX_CONDITION_LENGTH + 1 );
throws_ok { $sample->set_condition('') } qr/Empty condition specified/ms,
  'Empty condition';
throws_ok { $sample->set_condition($long_condition) }
qr/longer than \d+ characters/ms, 'Long condition';

# Test group attribute
is( $sample->group,          '1',   'Get group' );
is( $sample->set_group('2'), undef, 'Set group' );
is( $sample->group,          '2',   'Get new group' );
is( $sample->set_group(),    undef, 'Set undefined group' );
is( $sample->group,          undef, 'Get undefined group' );
my $long_group = 'X' x ( $DETCT::Sample::MAX_GROUP_LENGTH + 1 );
throws_ok { $sample->set_group('') } qr/Empty group specified/ms, 'Empty group';
throws_ok { $sample->set_group($long_group) } qr/longer than \d+ characters/ms,
  'Long group';

# Test tag attribute
is( $sample->tag,                   'NNNNBGAGGC', 'Get tag' );
is( $sample->set_tag('NNNNBCAGAG'), undef,        'Set tag' );
is( $sample->tag,                   'NNNNBCAGAG', 'Get new tag' );
throws_ok { $sample->set_tag() } qr/No tag specified/ms, 'No tag';
throws_ok { $sample->set_tag('NNNNBCAGAN') } qr/Invalid tag/ms, 'Invalid tag';

# Test bam file attribute
is( $sample->bam_file, 't/data/test1.bam', 'Get BAM file' );
is( $sample->set_bam_file('t/data/test2.bam'), undef, 'Set BAM file' );
is( $sample->bam_file, 't/data/test2.bam', 'Get new BAM file' );
throws_ok { $sample->set_bam_file() } qr/No BAM file specified/ms,
  'No BAM file';
throws_ok { $sample->set_bam_file('nonexistent.bam') }
qr/does not exist or cannot be read/ms, 'Missing BAM file';
