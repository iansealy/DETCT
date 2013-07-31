use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 12;

use DETCT::Sequence;

my $sequence = DETCT::Sequence->new(
    {
        name => '1',
        bp   => 60_348_388,
    }
);

isa_ok( $sequence, 'DETCT::Sequence' );

# Test name attribute
is( $sequence->name,          '1',   'Get name' );
is( $sequence->set_name('2'), undef, 'Set name' );
is( $sequence->name,          '2',   'Get new name' );
throws_ok { $sequence->set_name() } qr/No name specified/ms, 'No name';
my $long_name = 'X' x ( $DETCT::Sequence::MAX_NAME_LENGTH + 1 );
throws_ok { $sequence->set_name('') } qr/Empty name specified/ms, 'Empty name';
throws_ok { $sequence->set_name($long_name) } qr/longer than \d+ characters/ms,
  'Long name';

# Test bp attribute
is( $sequence->bp,                 60_348_388, 'Get bp' );
is( $sequence->set_bp(60_300_536), undef,      'Set bp' );
is( $sequence->bp,                 60_300_536, 'Get new bp' );
throws_ok { $sequence->set_bp() } qr/No bp specified/ms, 'No bp';
throws_ok { $sequence->set_bp(-1) } qr/Invalid bp/ms, 'Invalid bp';
