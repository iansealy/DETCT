use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 26;

use DETCT::Pipeline::Stage;

my $stage = DETCT::Pipeline::Stage->new(
    {
        name           => 'count_tags',
        default_memory => 3000,
    }
);

isa_ok( $stage, 'DETCT::Pipeline::Stage' );

# Test name attribute
is( $stage->name,                  'count_tags', 'Get name' );
is( $stage->set_name('bin_reads'), undef,        'Set name' );
is( $stage->name,                  'bin_reads',  'Get new name' );
throws_ok { $stage->set_name() } qr/No name specified/ms, 'No name';
throws_ok { $stage->set_name('/') } qr/Invalid name/ms, 'Invalid name';

# Test default memory attribute
is( $stage->default_memory,           3000,  'Get default memory' );
is( $stage->set_default_memory(2000), undef, 'Set default memory' );
is( $stage->default_memory,           2000,  'Get new default memory' );
throws_ok { $stage->set_default_memory() } qr/No default memory specified/ms,
  'No default memory';
throws_ok { $stage->set_default_memory(-1) } qr/Invalid default memory/ms,
  'Invalid default memory';

# Test threads attribute
is( $stage->threads,        undef, 'Get threads' );
is( $stage->set_threads(2), undef, 'Set threads' );
is( $stage->threads,        2,     'Get new threads' );
throws_ok { $stage->set_threads(-1) } qr/Invalid threads/ms, 'Invalid threads';

# Test all jobs run attribute
is( $stage->all_jobs_run,         0,     'Get all jobs run' );
is( $stage->set_all_jobs_run(10), undef, 'Set all jobs run to true' );
is( $stage->all_jobs_run,         1,     'Get new true all jobs run' );
is( $stage->set_all_jobs_run(),   undef, 'Set all jobs run to false' );
is( $stage->all_jobs_run,         0,     'Get new false all jobs run' );

# Mock stage object with different name
my $stage1 = Test::MockObject->new();
$stage1->set_isa('DETCT::Pipeline::Stage');
$stage1->set_always( 'name', 'get_read_peaks' );
my $stage2 = Test::MockObject->new();
$stage2->set_isa('DETCT::Pipeline::Stage');
$stage2->set_always( 'name', 'merge_read_peaks' );

# Test adding and retrieving prerequisites
my $prerequisites;
$prerequisites = $stage->get_all_prerequisites();
is( scalar @{$prerequisites},          0,     'No prerequisites' );
is( $stage->add_prerequisite($stage1), undef, 'Add prerequisite' );
$prerequisites = $stage->get_all_prerequisites();
is( scalar @{$prerequisites}, 1, 'Get one prerequisite' );
$stage->add_prerequisite($stage2);
is( scalar @{$prerequisites}, 2, 'Get two prerequisites' );
throws_ok { $stage->add_prerequisite() } qr/No prerequisite specified/ms,
  'No prerequisite specified';
throws_ok { $stage->add_prerequisite('invalid') } qr/Class of prerequisite/ms,
  'Invalid prerequisite';
