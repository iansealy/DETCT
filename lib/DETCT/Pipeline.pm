## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Pipeline;
## VERSION
## use critic

# ABSTRACT: Object representing a pipeline

## Author         : is1
## Maintainer     : is1
## Created        : 2013-01-09
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
use Scalar::Util qw( refaddr );
use English qw( -no_match_vars );
use IPC::System::Simple qw( capture );
use POSIX qw( WIFEXITED WIFSIGNALED WTERMSIG );
use Path::Tiny;
use File::Spec;
use File::Path qw( make_path );
use Hash::Merge;
use File::ReadBackwards;
use YAML::Tiny qw( DumpFile );
use YAML;
use JSON;
use Sys::Hostname;
use File::Basename;
use File::Find;
use DETCT;
use DETCT::Pipeline::Job;
use DETCT::Pipeline::Stage;

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private scheduler     => my %scheduler;        # e.g. lsf
private analysis_dir  => my %analysis_dir;     # e.g. .
private analysis      => my %analysis;         # DETCT::Analysis object
private cmd_line      => my %cmd_line;         # e.g. run_pipeline.pl
private avoid_node    => my %avoid_node;       # arrayref of nodes to be avoided
private max_jobs      => my %max_jobs;         # e.g. 1000
private max_retries   => my %max_retries;      # e.g. 10
private sleep_time    => my %sleep_time;       # e.g. 600
private skip_clean_up => my %skip_clean_up;    # e.g. 1
private serialiser_format       => my %serialiser_format;          # e.g. json
private memory_limit_multiplier => my %memory_limit_multiplier;    # e.g. 1000
private stage_to_run => my %stage_to_run;    # DETCT::Pipeline::Stage object
private component_to_run => my %component_to_run;    # e.g. 5
private all_stages_run   => my %all_stages_run;      # e.g. 1
private jobs_running     => my %jobs_running;        # e.g. 1
private verbose          => my %verbose;             # e.g. 1
private hash_merge       => my %hash_merge;          # Hash::Merge object
private stage            => my %stage;               # arrayref of stages

# Constants
Readonly our $MAX_MEMORY_BEFORE_HUGEMEM => 28_000;
Readonly our %EXTENSION_TO_KEEP => map { $_ => 1 } qw(
  csv html pdf tsv txt
);
Readonly our %EXTENSION_TO_DELETE => map { $_ => 1 } qw(
  bam bai
);

=method new

  Usage       : my $pipeline = DETCT::Pipeline->new( {
                    scheduler    => 'lsf',
                    analysis_dir => '.',
                    analysis     => $analysis,
                    cmd_line     => 'run_pipeline.pl',
                    max_jobs     => 1000,
                    max_retries  => 10,
                    sleep_time   => 600,
                    verbose      => 1,
                } );
  Purpose     : Constructor for pipeline objects
  Returns     : DETCT::Pipeline
  Parameters  : Hashref {
                    scheduler         => String,
                    analysis_dir      => String,
                    analysis          => DETCT::Analysis,
                    cmd_line          => String,
                    max_jobs          => Int,
                    max_retries       => Int,
                    sleep_time        => Int,
                    skip_clean_up     => Boolean or undef,
                    serialiser_format => String or undef,
                    verbose           => Boolean or undef,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_scheduler( $arg_ref->{scheduler} );
    $self->set_analysis_dir( $arg_ref->{analysis_dir} );
    $self->set_analysis( $arg_ref->{analysis} );
    $self->set_cmd_line( $arg_ref->{cmd_line} );
    $self->set_max_jobs( $arg_ref->{max_jobs} );
    $self->set_max_retries( $arg_ref->{max_retries} );
    $self->set_sleep_time( $arg_ref->{sleep_time} );
    $self->set_skip_clean_up( $arg_ref->{skip_clean_up} );
    $self->set_serialiser_format( $arg_ref->{serialiser_format} );
    $self->set_verbose( $arg_ref->{verbose} );
    return $self;
}

=method scheduler

  Usage       : my $scheduler = $pipeline->scheduler;
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

  Usage       : $pipeline->set_scheduler('lsf');
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

=method analysis_dir

  Usage       : my $analysis_dir = $pipeline->analysis_dir;
  Purpose     : Getter for analysis directory attribute
  Returns     : String (e.g. ".")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub analysis_dir {
    my ($self) = @_;
    return $analysis_dir{ id $self};
}

=method set_analysis_dir

  Usage       : $pipeline->set_analysis_dir('.');
  Purpose     : Setter for analysis directory attribute
  Returns     : undef
  Parameters  : String (the analysis directory)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_analysis_dir {
    my ( $self, $arg ) = @_;
    $analysis_dir{ id $self} = _check_analysis_dir($arg);
    return;
}

# Usage       : $analysis_dir = _check_analysis_dir($analysis_dir);
# Purpose     : Check for valid analysis directory
# Returns     : String (the valid analysis directory)
# Parameters  : String (the analysis directory)
# Throws      : If analysis directory is missing or invalid
# Comments    : None
sub _check_analysis_dir {
    my ($analysis_dir) = @_;

    # Make sure analysis directory exists
    if ( defined $analysis_dir && !-d $analysis_dir ) {
        make_path($analysis_dir);
    }

    return $analysis_dir if defined $analysis_dir && -d $analysis_dir;
    confess 'No analysis_dir specified' if !defined $analysis_dir;
    confess "Invalid analysis_dir ($analysis_dir) specified";
}

=method analysis

  Usage       : my $analysis = $pipeline->analysis;
  Purpose     : Getter for analysis attribute
  Returns     : DETCT::Analysis
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub analysis {
    my ($self) = @_;
    return $analysis{ id $self};
}

=method set_analysis

  Usage       : $pipeline->set_analysis($analysis);
  Purpose     : Setter for analysis attribute
  Returns     : undef
  Parameters  : DETCT::Analysis
  Throws      : No exceptions
  Comments    : None

=cut

sub set_analysis {
    my ( $self, $arg ) = @_;
    $analysis{ id $self} = _check_analysis($arg);
    return;
}

# Usage       : $analysis = _check_analysis($analysis);
# Purpose     : Check for valid analysis object
# Returns     : DETCT::Analysis
# Parameters  : DETCT::Analysis
# Throws      : If analysis object is missing or invalid (i.e. not a
#               DETCT::Analysis object)
# Comments    : None
sub _check_analysis {
    my ($analysis) = @_;
    return $analysis if defined $analysis && $analysis->isa('DETCT::Analysis');
    confess 'No analysis specified' if !defined $analysis;
    confess 'Class of analysis (', ref $analysis, ') not DETCT::Analysis';
}

=method cmd_line

  Usage       : my $cmd_line = $pipeline->cmd_line;
  Purpose     : Getter for command line attribute
  Returns     : String (e.g. "run_pipeline.pl")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub cmd_line {
    my ($self) = @_;
    return $cmd_line{ id $self};
}

=method set_cmd_line

  Usage       : $pipeline->set_cmd_line('run_pipeline.pl');
  Purpose     : Setter for command line attribute
  Returns     : undef
  Parameters  : String (the command line)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_cmd_line {
    my ( $self, $arg ) = @_;
    $cmd_line{ id $self} = _check_cmd_line($arg);
    return;
}

# Usage       : $cmd_line = _check_cmd_line($cmd_line);
# Purpose     : Check for valid command line
# Returns     : String (the valid command line)
# Parameters  : String (the command line)
# Throws      : If command line is missing
# Comments    : None
sub _check_cmd_line {
    my ($cmd_line) = @_;

    confess 'No command line specified' if !defined $cmd_line || !$cmd_line;

    return $cmd_line;
}

=method add_avoid_node

  Usage       : $pipeline->add_avoid_node('bc-25-2-02');
  Purpose     : Add a node to be avoided to a pipeline
  Returns     : undef
  Parameters  : String (the node to be avoided)
  Throws      : If node to be avoided is missing
  Comments    : None

=cut

sub add_avoid_node {
    my ( $self, $avoid_node ) = @_;

    confess 'No node to be avoided specified' if !defined $avoid_node;

    if ( !exists $avoid_node{ id $self} ) {
        $avoid_node{ id $self} = [$avoid_node];
    }
    else {
        push @{ $avoid_node{ id $self} }, $avoid_node;
    }

    return;
}

=method get_all_avoid_nodes

  Usage       : $avoid_nodes = $pipeline->get_all_avoid_nodes();
  Purpose     : Get all nodes to be avoided of a pipeline
  Returns     : Arrayref of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_avoid_nodes {
    my ($self) = @_;

    return $avoid_node{ id $self} || [];
}

=method max_jobs

  Usage       : my $max_jobs = $pipeline->max_jobs;
  Purpose     : Getter for max jobs attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub max_jobs {
    my ($self) = @_;
    return $max_jobs{ id $self};
}

=method set_max_jobs

  Usage       : $pipeline->set_max_jobs(1000);
  Purpose     : Setter for max jobs attribute
  Returns     : undef
  Parameters  : +ve Int (the max jobs)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_max_jobs {
    my ( $self, $arg ) = @_;
    $max_jobs{ id $self} = _check_max_jobs($arg);
    return;
}

# Usage       : $max_jobs = _check_max_jobs($max_jobs);
# Purpose     : Check for valid max jobs
# Returns     : +ve Int (the valid max jobs)
# Parameters  : +ve Int (the max jobs)
# Throws      : If max jobs is missing or not a positive integer
# Comments    : None
sub _check_max_jobs {
    my ($max_jobs) = @_;
    return $max_jobs
      if defined $max_jobs && $max_jobs =~ m/\A \d+ \z/xms;
    confess 'No max jobs specified' if !defined $max_jobs;
    confess "Invalid max jobs ($max_jobs) specified";
}

=method max_retries

  Usage       : my $max_retries = $pipeline->max_retries;
  Purpose     : Getter for max retries attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub max_retries {
    my ($self) = @_;
    return $max_retries{ id $self};
}

=method set_max_retries

  Usage       : $pipeline->set_max_retries(10);
  Purpose     : Setter for max retries attribute
  Returns     : undef
  Parameters  : +ve Int (the max retries)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_max_retries {
    my ( $self, $arg ) = @_;
    $max_retries{ id $self} = _check_max_retries($arg);
    return;
}

# Usage       : $max_retries = _check_max_retries($max_retries);
# Purpose     : Check for valid max retries
# Returns     : +ve Int (the valid max retries)
# Parameters  : +ve Int (the max retries)
# Throws      : If max retries is missing or not a positive integer
# Comments    : None
sub _check_max_retries {
    my ($max_retries) = @_;
    return $max_retries
      if defined $max_retries && $max_retries =~ m/\A \d+ \z/xms;
    confess 'No max retries specified' if !defined $max_retries;
    confess "Invalid max retries ($max_retries) specified";
}

=method sleep_time

  Usage       : my $sleep_time = $pipeline->sleep_time;
  Purpose     : Getter for sleep time attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub sleep_time {
    my ($self) = @_;
    return $sleep_time{ id $self};
}

=method set_sleep_time

  Usage       : $pipeline->set_sleep_time(600);
  Purpose     : Setter for sleep time attribute
  Returns     : undef
  Parameters  : +ve Int (the sleep time)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_sleep_time {
    my ( $self, $arg ) = @_;
    $sleep_time{ id $self} = _check_sleep_time($arg);
    return;
}

# Usage       : $sleep_time = _check_sleep_time($sleep_time);
# Purpose     : Check for valid sleep time
# Returns     : +ve Int (the valid sleep time)
# Parameters  : +ve Int (the sleep time)
# Throws      : If sleep time is missing or not a positive integer
# Comments    : None
sub _check_sleep_time {
    my ($sleep_time) = @_;
    return $sleep_time
      if defined $sleep_time && $sleep_time =~ m/\A \d+ \z/xms;
    confess 'No sleep time specified' if !defined $sleep_time;
    confess "Invalid sleep time ($sleep_time) specified";
}

=method memory_limit_multiplier

  Usage       : my $multiplier = $pipeline->memory_limit_multiplier;
  Purpose     : Getter for memory limit multiplier attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub memory_limit_multiplier {
    my ($self) = @_;

    # For LSF, units for memory limit (-M option) can vary by compute farm
    if ( !defined $memory_limit_multiplier{ id $self}
        && $self->scheduler eq 'lsf' )
    {
        ## no critic (ProhibitMagicNumbers)
        my $memory_limit_multiplier = 1000;    # Default
        ## use critic
        my $output = capture('lsadmin showconf lim');
        if ( $output =~ m/LSF_UNIT_FOR_LIMITS\s=\sMB/xms ) {
            $memory_limit_multiplier = 1;      # Only currently handle MB
        }
        $self->set_memory_limit_multiplier($memory_limit_multiplier);
    }

    return $memory_limit_multiplier{ id $self};
}

=method set_memory_limit_multiplier

  Usage       : $pipeline->set_memory_limit_multiplier(1000);
  Purpose     : Setter for memory limit multiplier attribute
  Returns     : undef
  Parameters  : +ve Int (the memory limit multiplier)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_memory_limit_multiplier {
    my ( $self, $arg ) = @_;
    $memory_limit_multiplier{ id $self} = _check_memory_limit_multiplier($arg);
    return;
}

# Usage       : $multiplier
#                   = _check_memory_limit_multiplier($memory_limit_multiplier);
# Purpose     : Check for valid memory limit multiplier
# Returns     : +ve Int (the valid memory limit multiplier)
# Parameters  : +ve Int (the memory limit multiplier)
# Throws      : If memory limit multiplier is missing or not a positive integer
# Comments    : None
sub _check_memory_limit_multiplier {
    my ($memory_limit_multiplier) = @_;
    return $memory_limit_multiplier
      if defined $memory_limit_multiplier
      && $memory_limit_multiplier =~ m/\A \d+ \z/xms;
    confess 'No memory limit multiplier specified'
      if !defined $memory_limit_multiplier;
    confess
      "Invalid memory limit multiplier ($memory_limit_multiplier) specified";
}

=method stage_to_run

  Usage       : my $stage = $pipeline->stage_to_run;
  Purpose     : Getter for stage to be run attribute
  Returns     : DETCT::Pipeline::Stage
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub stage_to_run {
    my ($self) = @_;
    return $stage_to_run{ id $self};
}

=method set_stage_to_run

  Usage       : $pipeline->set_stage_to_run($stage);
  Purpose     : Setter for stage to be run attribute
  Returns     : undef
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub set_stage_to_run {
    my ( $self, $arg ) = @_;
    $stage_to_run{ id $self} = _check_stage_to_run($arg);
    return;
}

# Usage       : $stage = _check_stage_to_run($stage);
# Purpose     : Check for valid stage to be run object
# Returns     : DETCT::Pipeline::Stage
# Parameters  : DETCT::Pipeline::Stage
# Throws      : If stage to be run object is missing or invalid (i.e. not a
#               DETCT::Pipeline::Stage object)
# Comments    : None
sub _check_stage_to_run {
    my ($stage_to_run) = @_;
    return $stage_to_run
      if defined $stage_to_run && $stage_to_run->isa('DETCT::Pipeline::Stage');
    confess 'No stage to be run specified' if !defined $stage_to_run;
    confess 'Class of stage to be run (', ref $stage_to_run,
      ') not DETCT::Pipeline::Stage';
}

=method component_to_run

  Usage       : my $component = $pipeline->component_to_run;
  Purpose     : Getter for component to be run attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub component_to_run {
    my ($self) = @_;
    return $component_to_run{ id $self};
}

=method set_component_to_run

  Usage       : $pipeline->set_component_to_run(5);
  Purpose     : Setter for component to be run attribute
  Returns     : undef
  Parameters  : +ve Int (the component to be run)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_component_to_run {
    my ( $self, $arg ) = @_;
    $component_to_run{ id $self} = _check_component_to_run($arg);
    return;
}

# Usage       : $component = _check_component_to_run($component);
# Purpose     : Check for valid component to be run
# Returns     : +ve Int (the valid component to be run)
# Parameters  : +ve Int (the component to be run)
# Throws      : If component to be run is missing or not a positive integer
# Comments    : None
sub _check_component_to_run {
    my ($component_to_run) = @_;
    return $component_to_run
      if defined $component_to_run && $component_to_run =~ m/\A \d+ \z/xms;
    confess 'No component to be run specified' if !defined $component_to_run;
    confess "Invalid component to be run ($component_to_run) specified";
}

=method all_stages_run

  Usage       : my $all_stages_run = $stage->all_stages_run;
  Purpose     : Getter for all stages run flag
  Returns     : Boolean
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_stages_run {
    my ($self) = @_;
    return $all_stages_run{ id $self} || 0;
}

=method set_all_stages_run

  Usage       : $stage->set_all_stages_run(1);
  Purpose     : Setter for all stages run flag
  Returns     : undef
  Parameters  : Boolean
  Throws      : No exceptions
  Comments    : None

=cut

sub set_all_stages_run {
    my ( $self, $arg ) = @_;
    $all_stages_run{ id $self} = $arg ? 1 : 0;
    return;
}

=method jobs_running

  Usage       : my $jobs_running = $stage->jobs_running;
  Purpose     : Getter for jobs running flag
  Returns     : Boolean
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub jobs_running {
    my ($self) = @_;
    return $jobs_running{ id $self} || 0;
}

=method set_jobs_running

  Usage       : $stage->set_jobs_running(1);
  Purpose     : Setter for jobs running flag
  Returns     : undef
  Parameters  : Boolean
  Throws      : No exceptions
  Comments    : None

=cut

sub set_jobs_running {
    my ( $self, $arg ) = @_;
    $jobs_running{ id $self} = $arg ? 1 : 0;
    return;
}

=method verbose

  Usage       : my $verbose = $pipeline->verbose;
  Purpose     : Getter for verbose flag
  Returns     : Boolean
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub verbose {
    my ($self) = @_;
    return $verbose{ id $self} || 0;
}

=method set_verbose

  Usage       : $pipeline->set_verbose(1);
  Purpose     : Setter for verbose flag
  Returns     : undef
  Parameters  : Boolean
  Throws      : No exceptions
  Comments    : None

=cut

sub set_verbose {
    my ( $self, $arg ) = @_;
    $verbose{ id $self} = $arg ? 1 : 0;
    return;
}

=method skip_clean_up

  Usage       : my $skip_clean_up = $pipeline->skip_clean_up;
  Purpose     : Getter for skip clean up flag
  Returns     : Boolean
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub skip_clean_up {
    my ($self) = @_;
    return $skip_clean_up{ id $self} || 0;
}

=method set_skip_clean_up

  Usage       : $pipeline->set_skip_clean_up(1);
  Purpose     : Setter for skip clean up flag
  Returns     : undef
  Parameters  : Boolean
  Throws      : No exceptions
  Comments    : None

=cut

sub set_skip_clean_up {
    my ( $self, $arg ) = @_;
    $skip_clean_up{ id $self} = $arg ? 1 : 0;
    return;
}

=method serialiser_format

  Usage       : my $serialiser_format = $pipeline->serialiser_format;
  Purpose     : Getter for serialiser format attribute
  Returns     : String (e.g. "json")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub serialiser_format {
    my ($self) = @_;
    return $serialiser_format{ id $self};
}

=method set_serialiser_format

  Usage       : $pipeline->set_serialiser_format('json');
  Purpose     : Setter for serialiser format attribute
  Returns     : undef
  Parameters  : String (the serialiser format)
  Throws      : No exceptions
  Comments    : Defaults to JSON

=cut

sub set_serialiser_format {
    my ( $self, $arg ) = @_;
    $serialiser_format{ id $self} = _check_serialiser_format( $arg || 'json' );
    return;
}

# Usage       : $serialiser_format = _check_serialiser_format($serialiser_format);
# Purpose     : Check for valid serialiser format
# Returns     : String (the valid serialiser format)
# Parameters  : String (the serialiser format)
# Throws      : If serialiser format is not json or yaml
# Comments    : None
sub _check_serialiser_format {
    my ($serialiser_format) = @_;

    confess 'Invalid serialiser format specified'
      if !defined $serialiser_format
      || ( $serialiser_format ne 'json' && $serialiser_format ne 'yaml' );

    return $serialiser_format;
}

=method hash_merge

  Usage       : %chunk_hmm
                    = %{ $pipeline->hash_merge->merge(\%chunk_hmm, $seq_hmm) };
  Purpose     : Return a Hash::Merge object for merging job output
  Returns     : Hash::Merge
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub hash_merge {
    my ($self) = @_;

    if ( !exists $hash_merge{ id $self} ) {
        ## no critic (ProtectPrivateSubs)
        Hash::Merge::specify_behavior(
            {
                SCALAR => {
                    SCALAR => sub { $_[0] + $_[1] },    # Add scalars
                    ARRAY  => sub { undef },
                    HASH   => sub { undef },
                },
                ARRAY => {
                    SCALAR => sub { undef },
                    ARRAY  => sub { [ @{ $_[0] }, @{ $_[1] } ] },  # Join arrays
                    HASH   => sub { undef },
                },
                HASH => {
                    SCALAR => sub { undef },
                    ARRAY  => sub { undef },
                    HASH => sub { Hash::Merge::_merge_hashes( $_[0], $_[1] ) },
                },
            },
            'detct',
        );
        ## use critic
        $hash_merge{ id $self} = Hash::Merge->new('detct');
    }

    return $hash_merge{ id $self};
}

=method add_stages_from_yaml

  Usage       : $pipeline->add_stages_from_yaml( 'detct.yaml' );
  Purpose     : Add stages from a YAML file
  Returns     : undef
  Parameters  : String (the YAML file)
  Throws      : If YAML file is missing or not readable or invalid
  Comments    : None

=cut

sub add_stages_from_yaml {
    my ( $self, $yaml_file ) = @_;

    confess "YAML file ($yaml_file) does not exist or cannot be read"
      if !-r $yaml_file;

    my $yaml = YAML::Tiny->read($yaml_file);

    if ( !$yaml ) {
        confess sprintf 'YAML file (%s) is invalid: %s', $yaml_file,
          YAML::Tiny->errstr;
    }

    my %tmp_cache;    # Temporarily store stages by name

    foreach my $stage_hash ( @{ $yaml->[0] } ) {
        my $stage = DETCT::Pipeline::Stage->new(
            {
                name           => $stage_hash->{name},
                default_memory => $stage_hash->{default_memory},
                threads        => $stage_hash->{threads},
            }
        );
        foreach my $prerequisite_name ( @{ $stage_hash->{prerequisites} } ) {
            $stage->add_prerequisite( $tmp_cache{$prerequisite_name} );
        }
        $self->add_stage($stage);

        $tmp_cache{ $stage_hash->{name} } = $stage;
    }

    return;
}

=method add_stage

  Usage       : $pipeline->add_stage($stage);
  Purpose     : Add a stage to a pipeline
  Returns     : undef
  Parameters  : DETCT::Pipeline::Stage
  Throws      : If stage is missing or invalid (i.e. not a
                DETCT::Pipeline::Stage object)
  Comments    : None

=cut

sub add_stage {
    my ( $self, $stage ) = @_;

    confess 'No stage specified' if !defined $stage;
    confess 'Class of stage (', ref $stage, ') not DETCT::Pipeline::Stage'
      if !$stage->isa('DETCT::Pipeline::Stage');

    if ( !exists $stage{ id $self} ) {
        $stage{ id $self} = [$stage];
    }
    else {
        push @{ $stage{ id $self} }, $stage;
    }

    return;
}

=method get_all_stages

  Usage       : $stages = $pipeline->get_all_stages();
  Purpose     : Get all stages of a pipeline
  Returns     : Arrayref of DETCT::Pipeline::Stage objects
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_stages {
    my ($self) = @_;

    return $stage{ id $self} || [];
}

=method get_stage_by_name

  Usage       : $stage = $pipeline->get_stage_by_name('run_deseq');
  Purpose     : Get a named stage of a pipeline
  Returns     : DETCT::Pipeline::Stage
  Parameters  : String (the stage name)
  Throws      : If stage with specified name does not exist
  Comments    : None

=cut

sub get_stage_by_name {
    my ( $self, $name ) = @_;

    foreach my $stage ( @{ $stage{ id $self} } ) {
        return $stage if $stage->name eq $name;
    }

    confess "Invalid stage name ($name)";
}

=method run

  Usage       : $pipeline->run();
  Purpose     : Run whole pipeline
  Returns     : undef
  Parameters  : None
  Throws      : If jobs can't be run but all stages have not finished
  Comments    : None

=cut

sub run {
    my ($self) = @_;

    $self->init_run();

    while ( !$self->all_stages_run ) {

        # Assume all run until find otherwise
        $self->set_all_stages_run(1);

        # Assume no jobs running until find otherwise
        $self->set_jobs_running(0);

        # Run a single cycle of the pipeline if number of pending jobs not high
        if ( $self->num_pending_jobs < $self->max_jobs ) {
            $self->_run_cycle();
        }
        else {
            $self->say_if_verbose(
                sprintf 'Not running jobs because >%d jobs pending',
                $self->max_jobs );
            $self->set_all_stages_run(0);
            $self->set_jobs_running(1);
        }

        if ( !$self->all_stages_run && !$self->jobs_running ) {
            $self->_delete_lock();
            die 'Stopping pipeline - no jobs to run' . "\n";
        }

        if ( !$self->all_stages_run ) {
            $self->say_if_verbose( sprintf 'Sleeping for %d seconds',
                $self->sleep_time );
            sleep $self->sleep_time;
        }
    }

    print 'Pipeline finished - all jobs run' . "\n";

    $self->clean_up();

    $self->end_run();

    return;
}

=method run_once

  Usage       : $pipeline->run_once();
  Purpose     : Run a single cycle of the pipeline
  Returns     : undef
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub run_once {
    my ($self) = @_;

    $self->init_run();

    # Assume all run until find otherwise
    $self->set_all_stages_run(1);

    # Assume no jobs running until find otherwise
    $self->set_jobs_running(0);

    # Run a single cycle of the pipeline if number of pending jobs not high
    if ( $self->num_pending_jobs < $self->max_jobs ) {
        $self->_run_cycle();
    }
    else {
        $self->say_if_verbose(
            sprintf 'Not running jobs because >%d jobs pending',
            $self->max_jobs );
        $self->set_all_stages_run(0);
        $self->set_jobs_running(1);
    }

    $self->end_run();

    return;
}

# Usage       : $self->_run_cycle();
# Purpose     : Run a single cycle of the pipeline
# Returns     : undef
# Parameters  : None
# Throws      : No exceptions
# Comments    : None
sub _run_cycle {
    my ($self) = @_;

    # Iterate over all stages of pipeline
  STAGE: foreach my $stage ( @{ $self->get_all_stages() } ) {

        # Check prerequisites have already run and skip this stage if not
        foreach my $prereq_stage ( @{ $stage->get_all_prerequisites() } ) {
            if ( !$prereq_stage->all_jobs_run ) {
                $self->say_if_verbose( sprintf 'Skipping %s because %s not run',
                    $stage->name, $prereq_stage->name );
                next STAGE;
            }
        }

        # Create directory for current stage of analysis
        my $dir = $self->get_and_create_stage_dir($stage);

        # Assume all jobs have run OK until we know otherwise
        $stage->set_all_jobs_run(1);

        # All jobs marked as having run OK?
        my $done_marker_file = $dir . '.done';
        if ( -e $done_marker_file ) {
            $self->say_if_verbose( sprintf 'Stage %s has finished',
                $stage->name );
            next STAGE;
        }

        # Running a specific stage, but not this one
        if ( $self->stage_to_run
            && refaddr( $self->stage_to_run ) != refaddr($stage) )
        {
            next STAGE;
        }

        # Get all parameters for all components of current stage
        my @all_parameters = $self->all_parameters($stage);

        $self->say_if_verbose( sprintf 'Stage %s has %d components',
            $stage->name, scalar @all_parameters );

        my $component = 0;    # Index for current component of current stage
        foreach my $parameters (@all_parameters) {
            $component++;

            # Running a specific component, but not this one
            if (
                   $self->stage_to_run
                && $self->component_to_run
                && ( refaddr( $self->stage_to_run ) != refaddr($stage)
                    || $self->component_to_run != $component )
              )
            {
                next;
            }

            my $job = DETCT::Pipeline::Job->new(
                {
                    stage         => $stage,
                    component     => $component,
                    scheduler     => $self->scheduler,
                    base_filename => File::Spec->catfile( $dir, $component ),
                    parameters    => $parameters,
                }
            );

            # Run job if running a specific component of a specific stage
            if ( $self->stage_to_run && $self->component_to_run ) {
                $self->run_job($job);
                exit;
            }

            $self->process_job($job);
        }

        if ( $stage->all_jobs_run ) {
            path($done_marker_file)->spew('1');
            $self->summarise_memory_usage( $stage, $component );
        }
        else {
            $self->set_all_stages_run(0);
        }
    }

    return;
}

=method init_run

  Usage       : $self->init_run();
  Purpose     : Initialise a pipeline run
  Returns     : undef
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub init_run {
    my ($self) = @_;

    if ( !$self->stage_to_run && !$self->component_to_run ) {
        ## no critic (RequireLocalizedPunctuationVars)
        $SIG{INT} = sub {
            $self->_delete_lock();
            die "\n" . 'Interrupted' . "\n";
        };
        ## use critic
        $self->_create_lock();
    }

    return;
}

=method end_run

  Usage       : $self->end_run();
  Purpose     : End a pipeline run
  Returns     : undef
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub end_run {
    my ($self) = @_;

    $SIG{'INT'} = 'DEFAULT';    ## no critic (RequireLocalizedPunctuationVars)
    $self->_delete_lock();

    return;
}

=method get_and_create_stage_dir

  Usage       : my $dir = $pipeline->get_and_create_stage_dir( $stage );
  Purpose     : Get (and create if necessary) a directory for the current stage
  Returns     : String (the directory)
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub get_and_create_stage_dir {
    my ( $self, $stage ) = @_;

    my $stage_dir = File::Spec->catdir( $self->analysis_dir, $stage->name );
    if ( !-d $stage_dir ) {
        make_path($stage_dir);
    }

    return $stage_dir;
}

=method get_and_check_output_file

  Usage       : my $file = $pipeline->get_and_check_output_file('run_deseq', 1);
  Purpose     : Get an output file for a particular component of a stage
  Returns     : String (the file)
  Parameters  : String (the stage)
                Int (the component)
  Throws      : If output file doesn't exist
  Comments    : None

=cut

sub get_and_check_output_file {
    my ( $self, $stage_name, $component ) = @_;

    my $output_file = File::Spec->catfile( $self->analysis_dir, $stage_name,
        $component . '.out' );
    if ( !-e $output_file ) {
        confess "$output_file doesn't exist, but should";
    }

    return $output_file;
}

=method process_job

  Usage       : $pipeline->process_job($job);
  Purpose     : Process a job to see if needs to be submitted
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub process_job {
    my ( $self, $job ) = @_;

    if ( $job->status_code eq 'NOT_RUN' ) {

        # Job not yet run so submit it
        $job->stage->set_all_jobs_run(0);
        $self->set_jobs_running(1);
        $self->submit_job($job);
        $self->say_if_verbose( sprintf '  Running component %d of %s%s',
            $job->component, $job->stage->name, $job->print_lsf_job_id );
    }
    elsif ( $job->status_code eq 'RUNNING' ) {

        # Job is running
        $job->stage->set_all_jobs_run(0);
        $self->set_jobs_running(1);
        $self->say_if_verbose(
            sprintf '  Component %d of %s is still running or pending%s',
            $job->component, $job->stage->name, $job->print_lsf_job_id );
    }
    elsif ( $job->status_code eq 'FAILED' ) {

        # Job has failed, so submit again
        $job->stage->set_all_jobs_run(0);
        $self->say_if_verbose(
            sprintf '  Component %d of %s has FAILED: %s%s',
            $job->component,   $job->stage->name,
            $job->status_text, $job->print_lsf_job_id
        );
        if ( $job->retries < $self->max_retries ) {
            $self->set_jobs_running(1);
            $self->submit_job($job);
            $self->say_if_verbose( sprintf '  Running component %d of %s%s',
                $job->component, $job->stage->name, $job->print_lsf_job_id );
        }
        else {
            $self->say_if_verbose(
                sprintf
                  '  Not running component %d of %s because retried %d times',
                $job->component, $job->stage->name, $job->retries );
        }
    }

    return;
}

=method submit_job

  Usage       : $pipeline->submit_job($job);
  Purpose     : Submit a job
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : If job or bsub can't be run
                If LSF job id can't be extracted from bsub output
  Comments    : None

=cut

sub submit_job {
    my ( $self, $job ) = @_;

    my $stdout_file = $job->base_filename . '.o';
    my $stderr_file = $job->base_filename . '.e';

    my $cmd =
        $self->cmd_line
      . ' --stage '
      . $job->stage->name
      . ' --component '
      . $job->component;

    if ( $job->scheduler eq 'local' ) {

        # Just run job
        $cmd .= ' 1>' . $stdout_file;
        $cmd .= ' 2>' . $stderr_file;
        my $cmd_status = system $cmd;

        # Die if the command was interrupted
        if ( WIFSIGNALED($cmd_status) && WTERMSIG($cmd_status) == 2 ) {
            $self->_delete_lock();
            die "\n" . 'Interrupted' . "\n";
        }

        # Die if the command couldn't be run
        confess "Couldn't run $cmd ($OS_ERROR)" if !WIFEXITED($cmd_status);

        if ( defined $job->retries ) {
            $job->set_retries( $job->retries + 1 );
        }
        else {
            $job->set_retries(0);
        }
        my $dump = { retries => $job->retries, };
        my $job_file = $job->base_filename . '.job';
        DumpFile( $job_file, $dump );
    }
    elsif ( $job->scheduler eq 'lsf' ) {

        # Either use default memory or increase by 50% (if retrying failed job)
        if ( !$job->memory ) {
            $job->set_memory( $job->stage->default_memory );
        }
        elsif ( $job->status_text =~ m/\A MEMLIMIT /xms ) {
            ## no critic (ProhibitMagicNumbers)
            $job->set_memory( int( $job->memory * 1.5 ) );
            ## use critic
        }

        # Switch queue if required
        if ( $job->memory >= $MAX_MEMORY_BEFORE_HUGEMEM ) {
            $job->set_queue('hugemem');
        }
        elsif ( $job->status_text && $job->status_text =~ m/\A RUNLIMIT /xms ) {
            $job->set_queue('long');
        }

        # bsub job
        my $bsub_stdout_file = $job->base_filename . '.bsub.o';
        my $bsub_stderr_file = $job->base_filename . '.bsub.e';
        my $queue_clause     = sprintf q{ -q %s }, $job->queue;
        my $memory_clause =
          sprintf q{ -R'span[hosts=1] select[mem>%d] rusage[mem=%d]' -M%d },
          $job->memory, $job->memory,
          $job->memory * $self->memory_limit_multiplier;
        my $threads_clause =
          $job->stage->threads == 1 ? q{} : sprintf q{ -n%d },
          $job->stage->threads;
        my $avoid_nodes_clause = q{};
        foreach my $avoid_node ( @{ $self->get_all_avoid_nodes() } ) {
            $avoid_nodes_clause .= sprintf q{ -R"select[hname!='%s']" },
              $avoid_node;
        }
        $cmd =
            'bsub' . ' -oo '
          . $stdout_file . ' -eo '
          . $stderr_file
          . $queue_clause
          . $memory_clause
          . $threads_clause
          . $avoid_nodes_clause
          . $cmd . ' 1>'
          . $bsub_stdout_file . ' 2>'
          . $bsub_stderr_file;
        WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

        # Extract LSF job id from bsub output and store with other parameters
        my $bsub_stdout = path($bsub_stdout_file)->slurp;
        if ( $bsub_stdout =~ m/Job \s <(\d+)> \s is \s submitted/xms ) {
            my $id = $1;
            $job->set_lsf_job_id($id);
            if ( defined $job->retries ) {
                $job->set_retries( $job->retries + 1 );
            }
            else {
                $job->set_retries(0);
            }
            my $dump = {
                id      => $id,
                retries => $job->retries,
                memory  => $job->memory,
                queue   => $job->queue,
            };
            my $job_file = $job->base_filename . '.job';
            DumpFile( $job_file, $dump );
        }
        else {
            confess "Couldn't get LSF job id from $bsub_stdout_file";
        }
    }

    return;
}

=method all_parameters

  Usage       : @all_parameters = $pipeline->all_parameters( $stage );
  Purpose     : Get all the parameters for a stage
  Returns     : Array
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : This function calls the all_parameters_for_ method associated
                with the current stage and gets all the parameters for that
                stage as an array of arbitrary data (e.g. arrayref or scalar)

=cut

sub all_parameters {
    my ( $self, $stage ) = @_;

    my $sub_name = 'all_parameters_for_' . $stage->name;

    return $self->$sub_name($stage);
}

=method run_job

  Usage       : $pipeline->run_job($job);
  Purpose     : Run a job
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : This function calls the run_ method associated with the current
                stage and passes along the parameters for the current component,
                which are arbitrary (e.g. arrayref or scalar)

=cut

sub run_job {
    my ( $self, $job ) = @_;

    my $sub_name = 'run_' . $job->stage->name;

    $self->$sub_name($job);

    return;
}

=method dump_serialised

  Usage       : $pipeline->dump_serialised($file, $data);
  Purpose     : Store serialised data in file
  Returns     : undef
  Parameters  : String (the file)
                Any reference (the data; usually hashref)
  Throws      : No exceptions
  Comments    : None

=cut

sub dump_serialised {
    my ( $self, $file, $data ) = @_;

    if ( $self->serialiser_format eq 'json' ) {
        path($file)->spew( JSON::to_json( $data, { pretty => 1 } ) );
    }
    elsif ( $self->serialiser_format eq 'yaml' ) {
        YAML::DumpFile( $file, $data );
    }

    return;
}

=method load_serialised

  Usage       : $data = $pipeline->load_serialised($file);
  Purpose     : Load serialised data from a file
  Returns     : Any reference (the data; usually hashref)
  Parameters  : String (the file)
  Throws      : No exceptions
  Comments    : None

=cut

sub load_serialised {
    my ( $self, $file ) = @_;

    if ( $self->serialiser_format eq 'json' ) {
        return JSON::from_json( path($file)->slurp );
    }
    elsif ( $self->serialiser_format eq 'yaml' ) {
        return YAML::LoadFile($file);
    }
}

=method input_overview

  Usage       : $pipeline->say_if_verbose($pipeline->input_overview);
  Purpose     : Return textual overview of pipeline's input
  Returns     : Array of Strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub input_overview {
    my ($self) = @_;

    my @output;

    push @output, 'Command line:', $self->cmd_line;
    if ( defined $DETCT::VERSION ) {
        push @output, 'DETCT version: ' . $DETCT::VERSION;
    }
    push @output, 'Working directory: ' . $self->analysis_dir;

    return @output;
}

=method say_if_verbose

  Usage       : $pipeline->say_if_verbose( 'Command line:', $cmd_line );
  Purpose     : Print output if pipeline is set to verbose
  Returns     : undef
  Parameters  : Array of Strings
  Throws      : No exceptions
  Comments    : Each string is a line without carriage returns or newlines

=cut

sub say_if_verbose {
    my ( $self, @output ) = @_;
    if ( $self->verbose ) {
        print join "\n", @output;
        print "\n";
    }
    return;
}

=method write_log_file

  Usage       : $pipeline->write_log_file( @output );
  Purpose     : Write data to a specified log file
  Returns     : undef
  Parameters  : String (the filename)
                Array of Strings
  Throws      : No exceptions
  Comments    : None

=cut

sub write_log_file {
    my ( $self, $filename, @output ) = @_;

    my $log_file = File::Spec->catfile( $self->analysis_dir, $filename );
    path($log_file)->spew(@output);

    return;
}

# Usage       : $self->_create_lock();
# Purpose     : Create lock file
# Returns     : undef
# Parameters  : None
# Throws      : If lock file already exists
# Comments    : None
sub _create_lock {
    my ($self) = @_;

    my $lock_file = File::Spec->catfile( $self->analysis_dir, 'pipeline.lock' );

    if ( -e $lock_file ) {
        my $message =
            "\nERROR: Is another pipeline running?\n"
          . "Make sure before deleting $lock_file and restarting.\n"
          . "Lock file contains:\n\n";
        $message .= path($lock_file)->slurp;
        die $message . "\n";
    }

    my $hostname  = hostname();
    my $timestamp = localtime;
    path($lock_file)->spew( $hostname . "\n" . $timestamp . "\n" );

    return;
}

# Usage       : $self->_delete_lock();
# Purpose     : Delete lock file
# Returns     : undef
# Parameters  : None
# Throws      : No exceptions
# Comments    : None
sub _delete_lock {
    my ($self) = @_;

    my $lock_file = File::Spec->catfile( $self->analysis_dir, 'pipeline.lock' );

    unlink $lock_file;

    return;
}

=method clean_up

  Usage       : $self->clean_up();
  Purpose     : Move results, archive stages and delete data
  Returns     : undef
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub clean_up {
    my ($self) = @_;

    return if $self->skip_clean_up;

    # Already cleaned up?
    my $done_marker_file =
      File::Spec->catfile( $self->analysis_dir, 'cleanup.done' );
    return if -e $done_marker_file;

    print 'Cleaning up...' . "\n";

    my @stage_dirs =
      map { $self->get_and_create_stage_dir($_) } @{ $self->get_all_stages() };

    my $wanted;

    # Delete files, so not archived
    $wanted = \&_delete;
    find(
        {
            wanted   => sub { $wanted->( $self->analysis_dir ) },
            no_chdir => 1,
        },
        @stage_dirs
    );

    # Tar all stages
    my $tarball_file =
      File::Spec->catfile( $self->analysis_dir, 'archive.tar.gz' );
    my $cmd = join q{ }, 'tar', 'cf', q{-}, @stage_dirs, q{|}, 'gzip', '-9',
      '-c', q{>}, $tarball_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    # Delete or move files
    $wanted = \&_move_or_delete;
    find(
        {
            wanted      => sub { $wanted->( $self->analysis_dir ) },
            postprocess => sub { rmdir $File::Find::dir },
            no_chdir    => 1,
        },
        @stage_dirs
    );

    path($done_marker_file)->spew('1');

    print 'Done' . "\n";

    return;
}

# Usage       : find(\&_move_or_delete, $dir);
# Purpose     : Move results files and delete other files
# Returns     : undef
# Parameters  : None
# Throws      : No exceptions
# Comments    : None
sub _move_or_delete {
    my ($archive_dir) = @_;

    return if -d;    # Ignore directories

    my ( $filename, undef, $extension ) =
      fileparse( $File::Find::name, qr/[.][^.]+/xms );
    $filename .= $extension;
    $extension =~ s/\A [.]//xms;

    # Move or delete?
    if ( $EXTENSION_TO_KEEP{$extension} ) {
        rename $File::Find::name,
          File::Spec->catfile( $archive_dir, $filename );
    }
    else {
        unlink $File::Find::name;
    }

    return;
}

# Usage       : find(\&_delete, $dir);
# Purpose     : Delete intermediate files with specific extensions
# Returns     : undef
# Parameters  : None
# Throws      : No exceptions
# Comments    : None
sub _delete {
    my ($archive_dir) = @_;

    return if -d;    # Ignore directories

    my ( $filename, undef, $extension ) =
      fileparse( $File::Find::name, qr/[.][^.]+/xms );
    $filename .= $extension;
    $extension =~ s/\A [.]//xms;

    # Delete?
    if ( $EXTENSION_TO_DELETE{$extension} ) {
        unlink $File::Find::name;
    }

    return;
}

=method num_pending_jobs

  Usage       : my $num_pending_jobs = $pipeline->num_pending_jobs;
  Purpose     : Get number of pending jobs
  Returns     : +ve Int
  Parameters  : None
  Throws      : If don't capture a positive integer
  Comments    : None

=cut

sub num_pending_jobs {
    my ($self) = @_;

    return 0 if $self->scheduler ne 'lsf';

    my @output = capture('busers');
    my @fields = split /\s+/xms, $output[1];    # Ignore header line
    my $num_pending_jobs = $fields[4];    ## no critic (ProhibitMagicNumbers)

    if ( $num_pending_jobs !~ m/\A \d+ \z/xms ) {
        confess 'Number of pending jobs not numeric';
    }

    return $num_pending_jobs;
}

=method summarise_memory_usage

  Usage       : $self->summarise_memory_usage( $stage, $last_component );
  Purpose     : Summarise memory usage of all components of a stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Stage
                Int (the last component)
  Throws      : If STDOUT file can't be read
  Comments    : None

=cut

sub summarise_memory_usage {
    my ( $self, $stage, $last_component ) = @_;

    return if $self->scheduler ne 'lsf';

    # Get directory for this stage
    my $dir = $self->get_and_create_stage_dir($stage);

    my @max_memory;
    foreach my $component ( 1 .. $last_component ) {
        my $stdout_file = File::Spec->catfile( $dir, $component . '.o' );
        my $bw = File::ReadBackwards->new($stdout_file)
          or confess "Can't read $stdout_file: $OS_ERROR";
        while ( defined( my $line = $bw->readline ) ) {
            if ( $line =~ m/\A \s+ Max \s Memory \s : \s+ (\d+) \s MB/xms ) {
                push @max_memory, $1;
                last;
            }
        }
    }

    # Output summary
    my $memory_summary_file = $dir . '.memory.txt';
    DumpFile( $memory_summary_file, \@max_memory );

    return;
}

1;
