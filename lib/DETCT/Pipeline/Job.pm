## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Pipeline::Job;
## use critic

# ABSTRACT: Object representing a pipeline job

## Author         : is1
## Maintainer     : is1
## Created        : 2013-01-17
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Readonly;
use Class::InsideOut qw( private register id );
use List::MoreUtils qw( any );
use English qw( -no_match_vars );
use File::ReadBackwards;
use YAML::Tiny qw( LoadFile );

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private stage         => my %stage;            # DETCT::Pipeline::Stage object
private component     => my %component;        # e.g. 2
private scheduler     => my %scheduler;        # e.g. lsf
private base_filename => my %base_filename;    # e.g. ./run_deseq/1
private parameters    => my %parameters;       # e.g. arrayref or scalar
private retries       => my %retries;          # e.g. 5
private lsf_job_id    => my %lsf_job_id;       # e.g. 123
private memory        => my %memory;           # e.g. 3000
private queue         => my %queue;            # e.g. normal
private status_code   => my %status_code;      # e.g. DONE
private status_text   => my %status_text;      # e.g. Job killed by owner

# Constants
Readonly our %STATUS_FOR => (
    PEND  => 'RUNNING',
    PSUSP => 'RUNNING',
    RUN   => 'RUNNING',
    USUSP => 'RUNNING',
    SSUSP => 'RUNNING',
    WAIT  => 'RUNNING',
    UNKWN => 'RUNNING',
    EXIT  => 'FAILED',
    ZOMBI => 'FAILED',
    DONE  => 'DONE',
);

=method new

  Usage       : my $job = DETCT::Pipeline::Job->new( {
                    stage         => $stage,
                    component     => 2,
                    scheduler     => 'lsf',
                    base_filename => './run_deseq/1',
                    parameters    => $parameters,
                } );
  Purpose     : Constructor for job objects
  Returns     : DETCT::Pipeline::Job
  Parameters  : Hashref {
                    stage         => DETCT::Pipeline::Stage,
                    component     => Int,
                    scheduler     => String,
                    base_filename => String,
                    parameters    => Any (probably arrayref or scalar)
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_stage( $arg_ref->{stage} );
    $self->set_component( $arg_ref->{component} );
    $self->set_scheduler( $arg_ref->{scheduler} );
    $self->set_base_filename( $arg_ref->{base_filename} );
    $self->set_parameters( $arg_ref->{parameters} );
    $self->set_state_from_filesystem();
    return $self;
}

=method stage

  Usage       : my $stage = $job->stage;
  Purpose     : Getter for stage attribute
  Returns     : DETCT::Pipeline::Stage
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub stage {
    my ($self) = @_;
    return $stage{ id $self};
}

=method set_stage

  Usage       : $job->set_stage($stage);
  Purpose     : Setter for stage attribute
  Returns     : undef
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub set_stage {
    my ( $self, $arg ) = @_;
    $stage{ id $self} = _check_stage($arg);
    return;
}

# Usage       : $stage = _check_stage($stage);
# Purpose     : Check for valid stage object
# Returns     : DETCT::Pipeline::Stage
# Parameters  : DETCT::Pipeline::Stage
# Throws      : If stage object is missing or invalid (i.e. not a
#               DETCT::Pipeline::Stage object)
# Comments    : None
sub _check_stage {
    my ($stage) = @_;
    return $stage if defined $stage && $stage->isa('DETCT::Pipeline::Stage');
    confess 'No stage specified' if !defined $stage;
    confess 'Class of stage (', ref $stage, ') not DETCT::Pipeline::Stage';
}

=method component

  Usage       : my $component = $job->component;
  Purpose     : Getter for component attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub component {
    my ($self) = @_;
    return $component{ id $self};
}

=method set_component

  Usage       : $job->set_component(2);
  Purpose     : Setter for component attribute
  Returns     : undef
  Parameters  : +ve Int (the component)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_component {
    my ( $self, $arg ) = @_;
    $component{ id $self} = _check_component($arg);
    return;
}

# Usage       : $component = _check_component($component);
# Purpose     : Check for valid component
# Returns     : +ve Int (the valid component)
# Parameters  : +ve Int (the component)
# Throws      : If component is missing or not a positive integer
# Comments    : None
sub _check_component {
    my ($component) = @_;
    return $component if defined $component && $component =~ m/\A \d+ \z/xms;
    confess 'No component specified' if !defined $component;
    confess "Invalid component ($component) specified";
}

=method scheduler

  Usage       : my $scheduler = $job->scheduler;
  Purpose     : Getter for scheduler attribute
  Returns     : String (e.g. "lsf")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub scheduler {
    my ($self) = @_;
    return $scheduler{ id $self};
}

=method set_scheduler

  Usage       : $job->set_scheduler('lsf');
  Purpose     : Setter for scheduler attribute
  Returns     : undef
  Parameters  : String (the scheduler)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_scheduler {
    my ( $self, $arg ) = @_;
    $scheduler{ id $self} = _check_scheduler($arg);
    return;
}

# Usage       : $scheduler = _check_scheduler($scheduler);
# Purpose     : Check for valid scheduler
# Returns     : String (the valid scheduler)
# Parameters  : String (the scheduler)
# Throws      : If scheduler is not lsf or local
# Comments    : None
sub _check_scheduler {
    my ($scheduler) = @_;

    confess 'Invalid scheduler specified'
      if !defined $scheduler
      || ( $scheduler ne 'lsf' && $scheduler ne 'local' );

    return $scheduler;
}

=method base_filename

  Usage       : my $base_filename = $job->base_filename;
  Purpose     : Getter for the base filename attribute
  Returns     : String (e.g. "./run_deseq/1")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub base_filename {
    my ($self) = @_;
    return $base_filename{ id $self};
}

=method set_base_filename

  Usage       : $job->set_base_filename('./run_deseq/1');
  Purpose     : Setter for the base filename attribute
  Returns     : undef
  Parameters  : String (the base filename)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_base_filename {
    my ( $self, $arg ) = @_;
    $base_filename{ id $self} = _check_base_filename($arg);
    return;
}

# Usage       : $base_filename = _check_base_filename($base_filename);
# Purpose     : Check for valid base filename
# Returns     : String (the valid base filename) or undef
# Parameters  : String (the base filename)
# Throws      : If base filename is missing
# Comments    : None
sub _check_base_filename {
    my ($base_filename) = @_;

    confess 'No base filename specified'
      if !defined $base_filename || !$base_filename;

    return $base_filename;
}

=method parameters

  Usage       : my $parameters = $job->parameters;
  Purpose     : Getter for parameters attribute
  Returns     : Any (usually arrayref or scalar)
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub parameters {
    my ($self) = @_;
    return $parameters{ id $self};
}

=method set_parameters

  Usage       : $job->set_parameters($parameters);
  Purpose     : Setter for parameters attribute
  Returns     : undef
  Parameters  : Any (the parameters; usually arrayref or scalar)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_parameters {
    my ( $self, $arg ) = @_;
    $parameters{ id $self} = $arg;
    return;
}

=method retries

  Usage       : my $retries = $job->retries;
  Purpose     : Getter for retries attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub retries {
    my ($self) = @_;
    return $retries{ id $self};
}

=method set_retries

  Usage       : $job->set_retries(5);
  Purpose     : Setter for retries attribute
  Returns     : undef
  Parameters  : +ve Int (the retries)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_retries {
    my ( $self, $arg ) = @_;
    $retries{ id $self} = _check_retries($arg);
    return;
}

# Usage       : $retries = _check_retries($retries);
# Purpose     : Check for valid retries
# Returns     : +ve Int (the valid retries)
# Parameters  : +ve Int (the retries)
# Throws      : If retries is not a positive integer
# Comments    : None
sub _check_retries {
    my ($retries) = @_;

    confess "Invalid retries ($retries) specified"
      if defined $retries && $retries !~ m/\A \d+ \z/xms;

    return $retries;
}

=method lsf_job_id

  Usage       : my $lsf_job_id = $job->lsf_job_id;
  Purpose     : Getter for LSF job id attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub lsf_job_id {
    my ($self) = @_;
    return $lsf_job_id{ id $self};
}

=method set_lsf_job_id

  Usage       : $job->set_lsf_job_id(123);
  Purpose     : Setter for LSF job id attribute
  Returns     : undef
  Parameters  : +ve Int (the LSF job id)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_lsf_job_id {
    my ( $self, $arg ) = @_;
    $lsf_job_id{ id $self} = _check_lsf_job_id($arg);
    return;
}

# Usage       : $lsf_job_id = _check_lsf_job_id($lsf_job_id);
# Purpose     : Check for valid LSF job id
# Returns     : +ve Int (the valid LSF job id)
# Parameters  : +ve Int (the LSF job id)
# Throws      : If LSF job id is not a positive integer
# Comments    : None
sub _check_lsf_job_id {
    my ($lsf_job_id) = @_;

    confess "Invalid LSF job id ($lsf_job_id) specified"
      if defined $lsf_job_id && $lsf_job_id !~ m/\A \d+ \z/xms;

    return $lsf_job_id;
}

=method print_lsf_job_id

  Usage       : print $job->print_lsf_job_id;
  Purpose     : Getter for LSF job id attribute in brackets
  Returns     : String (e.g. ' (123)')
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub print_lsf_job_id {
    my ($self) = @_;
    return $lsf_job_id{ id $self}
      ? sprintf ' (%s)', $lsf_job_id{ id $self}
      : q{};
}

=method memory

  Usage       : my $memory = $job->memory;
  Purpose     : Getter for memory attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub memory {
    my ($self) = @_;
    return $memory{ id $self};
}

=method set_memory

  Usage       : $job->set_memory(1000);
  Purpose     : Setter for memory attribute
  Returns     : undef
  Parameters  : +ve Int (the memory)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_memory {
    my ( $self, $arg ) = @_;
    $memory{ id $self} = _check_memory($arg);
    return;
}

# Usage       : $memory = _check_memory($memory);
# Purpose     : Check for valid memory
# Returns     : +ve Int (the valid memory)
# Parameters  : +ve Int (the memory)
# Throws      : If memory is not a positive integer
# Comments    : None
sub _check_memory {
    my ($memory) = @_;

    confess "Invalid memory ($memory) specified"
      if defined $memory && $memory !~ m/\A \d+ \z/xms;

    return $memory;
}

=method queue

  Usage       : my $queue = $job->queue;
  Purpose     : Getter for the queue attribute
  Returns     : String (e.g. "normal")
  Parameters  : None
  Throws      : No exceptions
  Comments    : Queue can be normal, long, or hugemem

=cut

sub queue {
    my ($self) = @_;
    return $queue{ id $self};
}

=method set_queue

  Usage       : $job->set_queue('long');
  Purpose     : Setter for the queue attribute
  Returns     : undef
  Parameters  : String (the queue)
  Throws      : No exceptions
  Comments    : Defaults to normal

=cut

sub set_queue {
    my ( $self, $arg ) = @_;
    $queue{ id $self} = _check_queue( $arg || 'normal' );
    return;
}

# Usage       : $queue = _check_queue($queue);
# Purpose     : Check for valid queue
# Returns     : String (the valid queue)
# Parameters  : String (the queue)
# Throws      : If queue is not valid
# Comments    : None
sub _check_queue {
    my ($queue) = @_;

    return $queue
      if any { $_ eq $queue } qw(normal long hugemem);
    confess "Invalid queue ($queue) specified";
}

=method status_code

  Usage       : my $status_code = $job->status_code;
  Purpose     : Getter for the status code attribute
  Returns     : String (e.g. "DONE")
  Parameters  : None
  Throws      : No exceptions
  Comments    : Status code can be RUNNING, FAILED, DONE or NOT_RUN

=cut

sub status_code {
    my ($self) = @_;
    return $status_code{ id $self};
}

=method set_status_code

  Usage       : $job->set_status_code('DONE');
  Purpose     : Setter for the status code attribute
  Returns     : undef
  Parameters  : String (the status code)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_status_code {
    my ( $self, $arg ) = @_;
    $status_code{ id $self} = _check_status_code($arg);
    return;
}

# Usage       : $status_code = _check_status_code($status_code);
# Purpose     : Check for valid status code
# Returns     : String (the valid status code)
# Parameters  : String (the status code)
# Throws      : If status code is not valid
# Comments    : None
sub _check_status_code {
    my ($status_code) = @_;

    return $status_code
      if defined $status_code
      && ( $status_code eq 'RUNNING'
        || $status_code eq 'FAILED'
        || $status_code eq 'DONE'
        || $status_code eq 'NOT_RUN' );
    confess 'No status code specified' if !defined $status_code;
    confess "Invalid status code ($status_code) specified";
}

=method status_text

  Usage       : my $status_text = $job->status_text;
  Purpose     : Getter for status text attribute
  Returns     : String (e.g. "Job killed by owner")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub status_text {
    my ($self) = @_;
    return $status_text{ id $self};
}

=method set_status_text

  Usage       : $job->set_status_text('Job killed by owner');
  Purpose     : Setter for status text attribute
  Returns     : undef
  Parameters  : String (the status text)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_status_text {
    my ( $self, $arg ) = @_;
    $status_text{ id $self} = $arg;
    return;
}

=method set_state_from_filesystem

  Usage       : $job->set_state_from_filesystem();
  Purpose     : Set state-related attributes of a job from filesystem
  Returns     : undef
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub set_state_from_filesystem {
    my ($self) = @_;

    if ( $self->scheduler eq 'local' ) {
        $self->_set_state_from_filesystem_for_local();
    }
    elsif ( $self->scheduler eq 'lsf' ) {
        $self->_set_state_from_filesystem_for_lsf();
    }

    return;
}

# Usage       : $self->_set_state_from_filesystem_for_local();
# Purpose     : Set state-related attributes of a job run locally (probably for
#               testing)
# Returns     : None
# Parameters  : None
# Throws      : No exceptions
# Comments    : None
sub _set_state_from_filesystem_for_local {
    my ($self) = @_;

    my $output_file = $self->base_filename . '.out';

    my $job_file = $self->base_filename . '.job';

    # Check if job has even been run yet
    if ( !-e $job_file ) {
        $self->set_status_code('NOT_RUN');
        return;
    }

    # Get number of retries
    my $yaml    = LoadFile($job_file);
    my $retries = $yaml->{retries};
    $self->set_retries($retries);

    if ( -e $output_file && !-z $output_file ) {
        $self->set_status_code('DONE');
    }
    elsif ( !-e $output_file || -z $output_file ) {
        $self->set_status_code('FAILED');
        $self->set_status_text( 'Empty output file: ' . $output_file );
    }

    return;
}

# Usage       : $self->_set_state_from_filesystem_for_lsf();
# Purpose     : Set state-related attributes of a job submitted to LSF
# Returns     : None
# Parameters  : None
# Throws      : If the status returned by bjobs is not recognised
# Comments    : None
sub _set_state_from_filesystem_for_lsf {
    my ($self) = @_;

    my $output_file = $self->base_filename . '.out';

    my $job_file = $self->base_filename . '.job';

    # Check if job has even been run yet
    if ( !-e $job_file ) {
        $self->set_status_code('NOT_RUN');
        return;
    }

    my $yaml = LoadFile($job_file);

    # Get LSF job id
    my $lsf_job_id = $yaml->{id};
    $self->set_lsf_job_id($lsf_job_id);

    # Get number of retries
    my $retries = $yaml->{retries};
    $self->set_retries($retries);

    # Get memory requested
    my $memory = $yaml->{memory};
    $self->set_memory($memory);

    # Get queue requested
    my $queue = $yaml->{queue};
    $self->set_queue($queue);

    my ( $status_code, $status_text );

    # Get job status for LSF job id from bjobs command
    my $lsf_status;
    open my $pipe, q{-|}, 'bjobs ' . $lsf_job_id . ' 2>/dev/null'; # Hide STDERR
    while ( my $job_line = <$pipe> ) {
        if ( $job_line =~ m/\A $lsf_job_id \s+ \S+ \s+ (\S+)/xms ) {
            $lsf_status = $1;
        }
    }
    close $pipe;
    if ($lsf_status) {

        # Got job status from bjobs
        if ( !exists $STATUS_FOR{$lsf_status} ) {
            confess "Unknown LSF status ($lsf_status)";
        }
        $status_code = $STATUS_FOR{$lsf_status};
        $status_text = 'LSF status: ' . $lsf_status;
    }
    if ( !$status_code || $status_code eq 'FAILED' ) {

        # If bjobs doesn't return status or failed then check job's STDOUT
        ( $status_code, $status_text ) = $self->_parse_lsf_stdout();
    }

    $self->set_status_code($status_code);
    $self->set_status_text($status_text);

    return;
}

# Usage       : ($status_code, $status_text) = $pipeline->_parse_lsf_stdout();
# Purpose     : Parses LSF's STDOUT to get a job's status
# Returns     : String (status code: DONE or FAILED)
#               String (status info) or undef
# Parameters  : None
# Throws      : If STDOUT file can't be read
# Comments    : Based on
#               https://github.com/VertebrateResequencing/vr-pipe/blob/master/modules/VRPipe/Parser/lsf.pm
sub _parse_lsf_stdout {
    my ($self) = @_;

    # Check STDOUT file exists at all (in case job was killed whilst pending)
    my $stdout_file = $self->base_filename . '.o';
    if ( !-e $stdout_file ) {
        return 'FAILED', 'Job did not run';
    }

    # STDOUT file is overwritten so no need to read backwards to get last job
    my ( $status_code, $status_text, $stdout_job_id );
    my $found_start = 0;
    my $found_end   = 0;
    my $bw          = File::ReadBackwards->new($stdout_file)
      or confess "Can't read $stdout_file: $OS_ERROR";
    while ( defined( my $line = $bw->readline ) ) {
        if ( $line =~ m/\A Resource \s usage \s summary: /xms ) {
            $found_end = 1;
            next;
        }
        elsif ( $line =~ m/\A Sender: \s LSF \s System/xms ) {
            $found_start = 1;
            last;    # Will find start after end
        }
        elsif ($found_end) {

            # Get job id
            if ( $line =~ m/\A Subject: \s Job \s (\d+):/xms ) {
                $stdout_job_id = $1;
            }

            # Get job's status code
            if ( $line =~ m/\A Successfully \s completed[.] /xms ) {
                $status_code = 'DONE';
            }
            elsif ( !$status_code
                && $line =~
                m/\A Exited \s with \s exit \s code \s (\d+) [.] /xms )
            {
                $status_code = 'FAILED';
                $status_text = "Exit code: $1";
            }
            elsif ( $line =~ m/\A TERM_ (\w+: .*) [.] /xms ) {
                $status_code = 'FAILED';
                $status_text = $1;
            }
        }
    }

    # Ensure correct job
    if ( defined $stdout_job_id && $self->lsf_job_id != $stdout_job_id ) {
        $status_code = 'FAILED';
        $status_text = sprintf 'Wrong LSF job id (expecting %s, got %s)',
          $self->lsf_job_id, $stdout_job_id;
    }

    # If no status then STDOUT could not be parsed
    if ( !defined $status_code ) {
        $status_code = 'FAILED';
        $status_text = "Could not parse job's STDOUT: $stdout_file";
    }

    return $status_code, $status_text;
}

1;
