#!/usr/bin/env perl

# PODNAME: run_multiple_pipelines.pl
# ABSTRACT: Run multiple DETCT pipelines simultaneously

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-06
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Probe::Perl;
use Getopt::Long;
use Pod::Usage;
use English qw( -no_match_vars );
use File::Spec;
use Path::Tiny;
use DETCT::Analysis::DiffExpr;
use DETCT::Analysis::Downsample;
use DETCT::Pipeline::DiffExpr;
use DETCT::Pipeline::Downsample;

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $pipeline_type;
my $scheduler = 'lsf';
my @analysis_yamls;
my $stages_yaml = 'stages.yaml';
my @avoid_nodes;
## no critic (ProhibitMagicNumbers)
my $max_jobs    = 1000;
my $max_retries = 10;
my $sleep_time  = 600;    # 10 minutes
## use critic
my $skip_clean_up;
my $serialiser_format = 'json';
my $stage_to_run;
my $component_to_run;
my $verbose;
my ( $help, $man );

# Get base command line (including interpreter and options)
my $base_cmd_line = get_base_cmd_line();

# Get and check command line options
get_and_check_options();

# Modules
my $module          = $pipeline_type eq 'de' ? 'DiffExpr' : 'Downsample';
my $analysis_module = 'DETCT::Analysis::' . $module;
my $pipeline_module = 'DETCT::Pipeline::' . $module;

# Extend base command line with common options
$base_cmd_line .=
    ' --pipeline '
  . $pipeline_type
  . ' --scheduler '
  . $scheduler
  . ' --stages_yaml '
  . $stages_yaml
  . ' --max_jobs '
  . $max_jobs
  . ' --max_retries '
  . $max_retries
  . ' --serialiser_format '
  . $serialiser_format;
if (@avoid_nodes) {
    $base_cmd_line .= ' --avoid_nodes ' . join q{ }, @avoid_nodes;
}

# Create Pipeline object for each pipeline
my @pipelines;
foreach my $analysis_yaml (@analysis_yamls) {

    # Create analysis
    my $analysis = $analysis_module->new_from_yaml($analysis_yaml);

    # Working directory is based on YAML filename
    my $analysis_dir = $analysis_yaml;
    $analysis_dir =~ s/[.]yaml \z//xms;

    # Skip already completed pipelines
    next if -e File::Spec->catfile( $analysis_dir, 'cleanup.done' );

    # Finish command line
    my $cmd_line = $base_cmd_line;
    $cmd_line .=
      ' --dir ' . $analysis_dir . ' --analysis_yaml ' . $analysis_yaml;

    # Create pipeline
    my $pipeline = $pipeline_module->new(
        {
            scheduler         => $scheduler,
            analysis_dir      => $analysis_dir,
            analysis          => $analysis,
            cmd_line          => $cmd_line,
            max_jobs          => $max_jobs,
            max_retries       => $max_retries,
            sleep_time        => $sleep_time,
            skip_clean_up     => $skip_clean_up,
            serialiser_format => $serialiser_format,
            verbose           => $verbose,
        }
    );

    # Add nodes to be avoided
    foreach my $avoid_node (@avoid_nodes) {
        $pipeline->add_avoid_node($avoid_node);
    }

    # Add stages to pipeline
    $pipeline->add_stages_from_yaml($stages_yaml);

    # Write overview of pipeline input and config files to log file
    my @log = map { "$_\n" } $pipeline->input_overview;
    push @log, "\nYAML analysis config file:\n\n", path($analysis_yaml)->slurp;
    push @log, "\nYAML stages config file:\n\n",   path($stages_yaml)->slurp;
    $pipeline->write_log_file( $pipeline_type . '.log', @log );

    push @pipelines, $pipeline;
}

# Run pipelines consecutively
my $all_pipelines_run = 0;
my $all_pipelines_ok  = 1;
while ( !$all_pipelines_run ) {

    # Assume all run until find otherwise
    $all_pipelines_run = 1;
    foreach my $pipeline (@pipelines) {
        printf "* Running %s pipeline *\n", $pipeline->analysis_dir;
        $pipeline->run_once();
        if ( !$pipeline->all_stages_run && !$pipeline->jobs_running ) {
            printf "* No jobs to run for %s pipeline *\n",
              $pipeline->analysis_dir;
            $all_pipelines_ok = 0;
        }
        elsif ( $pipeline->all_stages_run ) {
            printf "* %s pipeline finished - all jobs run *\n",
              $pipeline->analysis_dir;
            $pipeline->clean_up();
        }
        else {
            $all_pipelines_run = 0;
        }
    }
    if ( !$all_pipelines_run ) {
        printf "* Sleeping for %d seconds *\n", $sleep_time;
        sleep $sleep_time;
    }
}

if ($all_pipelines_ok) {
    print '* All pipelines finished OK *', "\n";
}
else {
    print '* All pipelines not finished, but no jobs to run *', "\n";
}

# Get base command line
sub get_base_cmd_line {

    # Get all lib directories
    my %lib = map { $_ => 1 } @INC;

    # Remove default lib directories
    foreach my $lib ( Probe::Perl->perl_inc() ) {
        delete $lib{$lib};
    }

    # Remove PERL5LIB lib directories
    foreach my $lib ( split /:/xms, $ENV{PERL5LIB} ) {
        delete $lib{$lib};
    }

    # Reconstruct -I lib directories
    my @libs;
    foreach my $lib ( keys %lib ) {
        push @libs, '-I' . $lib;
    }

    # Alter script name
    my $program_name = $PROGRAM_NAME;
    $program_name =~ s/run_multiple_(de_)*pipelines/run_pipeline/xms;

    return join q{ }, Probe::Perl->find_perl_interpreter(), @libs,
      $program_name;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'pipeline=s'           => \$pipeline_type,
        'scheduler=s'          => \$scheduler,
        'analysis_yamls=s{1,}' => \@analysis_yamls,
        'stages_yaml=s'        => \$stages_yaml,
        'avoid_nodes=s@{,}'    => \@avoid_nodes,
        'max_jobs=i'           => \$max_jobs,
        'max_retries=i'        => \$max_retries,
        'sleep_time=i'         => \$sleep_time,
        'skip_clean_up'        => \$skip_clean_up,
        'serialiser_format=s'  => \$serialiser_format,
        'stage=s'              => \$stage_to_run,
        'component=i'          => \$component_to_run,
        'verbose'              => \$verbose,
        'help'                 => \$help,
        'man'                  => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # If necessary, assign pipeline based on script name
    if ( !$pipeline_type ) {
        if ( $PROGRAM_NAME =~ m/run_multiple_de_pipelines[.]pl \z/xms ) {
            $pipeline_type = 'de';
        }
        elsif (
            $PROGRAM_NAME =~ m/run_multiple_downsample_pipelines[.]pl \z/xms )
        {
            $pipeline_type = 'downsample';
        }
    }

    # Check options
    if ( !$pipeline_type
        || ( $pipeline_type ne 'de' && $pipeline_type ne 'downsample' ) )
    {
        pod2usage("--pipeline must be 'de' or 'downsample'\n");
    }
    if ( $scheduler ne 'lsf' && $scheduler ne 'local' ) {
        pod2usage("--scheduler must be 'lsf' or 'local'\n");
    }
    if ( $stage_to_run && !$component_to_run
        || !$stage_to_run && $component_to_run )
    {
        pod2usage("--stage and --component must be specified together\n");
    }

    return;
}

=head1 USAGE

    run_multiple_pipelines.pl
        [--pipeline de|downsample]
        [--scheduler lsf|local]
        [--analysis_yamls file...]
        [--stages_yaml file]
        [--avoid_nodes node...]
        [--max_jobs int]
        [--max_retries int]
        [--sleep_time int]
        [--skip_clean_up]
        [--serialiser_format json|yaml]
        [--stage stage]
        [--component int]
        [--verbose]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--pipeline de|downsample>

Pipeline to run - de or downsample.

=item B<--scheduler lsf|local>

Job scheduler - lsf (default) or local (for testing).

=item B<--analysis_yamls FILE...>

YAML analysis configuration files.

=item B<--stages_yaml FILE>

YAML stages configuration file.

=item B<--avoid_nodes NODE...>

Nodes to be avoided when submitting LSF jobs.

=item B<--max_jobs INT>

Approximate maximum number of jobs user can have pending.

=item B<--max_retries INT>

Maximum number of times to retry a failing job.

=item B<--sleep_time INT>

Time to sleep, in seconds, between each cycle of the pipeline.

=item B<--skip_clean_up>

Skip final clean up stage.

=item B<--serialiser_format json|yaml>

Internal serialisation format - json (default) or yaml.

=item B<--stage STAGE>

The specific stage of the pipeline to be run.

=item B<--component INT>

The index of the component of the specified stage of the pipeline to be run.

=item B<--verbose>

Print information about the pipeline as it runs.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
