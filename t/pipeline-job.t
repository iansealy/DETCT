use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 61;

use DETCT::Pipeline::Job;

# Mock stage objects with different names
my $stage1 = Test::MockObject->new();
$stage1->set_isa('DETCT::Pipeline::Stage');
$stage1->set_always( 'name', 'get_read_peaks' );
my $stage2 = Test::MockObject->new();
$stage2->set_isa('DETCT::Pipeline::Stage');
$stage2->set_always( 'name', 'merge_read_peaks' );

my $job = DETCT::Pipeline::Job->new(
    {
        stage         => $stage1,
        component     => 2,
        scheduler     => 'local',
        base_filename => './run_deseq/1',
    }
);

isa_ok( $job, 'DETCT::Pipeline::Job' );

# Test stage attribute
is( $job->stage->name,        'get_read_peaks',   'Get stage' );
is( $job->set_stage($stage2), undef,              'Set stage' );
is( $job->stage->name,        'merge_read_peaks', 'Get new stage' );
throws_ok { $job->set_stage() } qr/No stage specified/ms, 'No stage';
throws_ok { $job->set_stage('invalid') } qr/Class of stage/ms, 'Invalid stage';

# Test component attribute
is( $job->component,        2,     'Get component' );
is( $job->set_component(3), undef, 'Set component' );
is( $job->component,        3,     'Get new component' );
throws_ok { $job->set_component() } qr/No component specified/ms,
  'No component';
throws_ok { $job->set_component(-1) } qr/Invalid component/ms,
  'Invalid component';

# Test scheduler attribute
is( $job->scheduler,            'local', 'Get scheduler' );
is( $job->set_scheduler('lsf'), undef,   'Set scheduler' );
is( $job->scheduler,            'lsf',   'Get new scheduler' );
throws_ok { $job->set_scheduler() } qr/Invalid scheduler specified/ms,
  'No scheduler';
throws_ok { $job->set_scheduler('invalid') } qr/Invalid scheduler specified/ms,
  'Invalid scheduler';

# Test base_filename attribute
is( $job->base_filename, './run_deseq/1', 'Get base filename' );
is( $job->set_base_filename('./count_reads/2'), undef, 'Set base filename' );
is( $job->base_filename, './count_reads/2', 'Get new base filename' );
throws_ok { $job->set_base_filename() } qr/No base filename specified/ms,
  'No base filename';
throws_ok { $job->set_base_filename('') } qr/No base filename specified/ms,
  'Empty base filename';

# Test parameters attribute
is( $job->parameters,             undef,  'Get parameters' );
is( $job->set_parameters('test'), undef,  'Set parameters' );
is( $job->parameters,             'test', 'Get new parameters' );
is( $job->set_parameters(),       undef,  'Set undef parameters' );
is( $job->parameters,             undef,  'Get undef parameters' );

# Test memory attribute
is( $job->memory,           undef, 'Get memory' );
is( $job->set_memory(1000), undef, 'Set memory' );
is( $job->memory,           1000,  'Get new memory' );
throws_ok { $job->set_memory(-1) } qr/Invalid memory/ms, 'Invalid memory';

# Test queue attribute
is( $job->queue,                undef,     'Get queue' );
is( $job->set_queue('normal'),  undef,     'Set normal queue' );
is( $job->queue,                'normal',  'Get normal queue' );
is( $job->set_queue('long'),    undef,     'Set long queue' );
is( $job->queue,                'long',    'Get long queue' );
is( $job->set_queue('hugemem'), undef,     'Set hugemem queue' );
is( $job->queue,                'hugemem', 'Get hugemem queue' );
throws_ok { $job->set_queue('invalid') } qr/Invalid queue/ms, 'Invalid queue';

# Test retries attribute
is( $job->retries,        undef, 'Get retries' );
is( $job->set_retries(5), undef, 'Set retries' );
is( $job->retries,        5,     'Get new retries' );
throws_ok { $job->set_retries(-1) } qr/Invalid retries/ms, 'Invalid retries';

# Test LSF job id attribute
is( $job->lsf_job_id,          undef, 'Get LSF job id' );
is( $job->set_lsf_job_id(123), undef, 'Set LSF job id' );
is( $job->lsf_job_id,          123,   'Get new LSF job id' );
throws_ok { $job->set_lsf_job_id(-1) } qr/Invalid LSF job id/ms,
  'Invalid LSF job id';
is( $job->print_lsf_job_id, ' (123)', 'Get bracketed LSF job id' );

# Test status code attribute
is( $job->status_code,                'NOT_RUN', 'Get not run status code' );
is( $job->set_status_code('DONE'),    undef,     'Set done status code' );
is( $job->status_code,                'DONE',    'Get done status code' );
is( $job->set_status_code('RUNNING'), undef,     'Set running status code' );
is( $job->status_code,               'RUNNING', 'Get new running status code' );
is( $job->set_status_code('FAILED'), undef,     'Set failed status code' );
is( $job->status_code,               'FAILED',  'Get new failed status code' );
throws_ok { $job->set_status_code() } qr/No status code specified/ms,
  'No status code';
throws_ok { $job->set_status_code('invalid') } qr/Invalid status code/ms,
  'Invalid status code';

# Test status text attribute
is( $job->status_text,                            undef, 'Get status text' );
is( $job->set_status_text('Job killed by owner'), undef, 'Set status text' );
is( $job->status_text,       'Job killed by owner', 'Get new status text' );
is( $job->set_status_text(), undef,                 'Set undef status text' );
is( $job->status_text,       undef,                 'Get undef status text' );
