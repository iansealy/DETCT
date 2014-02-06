#!/usr/bin/env perl

# PODNAME: run_multiple_de_pipelines.pl
# ABSTRACT: Run multiple DETCT differential expression pipeline simultaneously

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
use File::Slurp;
use DETCT::Pipeline::WithDiffExprStages;
use DETCT::Analysis;

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $scheduler = 'lsf';
my @analysis_yamls;
my $stages_yaml = 'stages.yaml';
my @avoid_nodes;
## no critic (ProhibitMagicNumbers)
my $max_retries = 10;
my $sleep_time  = 600;    # 10 minutes
## use critic
my $stage_to_run;
my $component_to_run;
my $verbose;
my ( $help, $man );

# Get base command line (including interpreter and options)
my $base_cmd_line = get_base_cmd_line();

# Get and check command line options
get_and_check_options();

# Extend base command line with common options
$base_cmd_line .=
    ' --scheduler '
  . $scheduler
  . ' --stages_yaml '
  . $stages_yaml
  . ' --max_retries '
  . $max_retries;
if (@avoid_nodes) {
    $base_cmd_line .= ' --avoid_nodes ' . join q{ }, @avoid_nodes;
}

# Create Pipeline object for each pipeline
my @pipelines;
foreach my $analysis_yaml (@analysis_yamls) {

    # Create analysis
    my $analysis = DETCT::Analysis->new_from_yaml($analysis_yaml);

    # Working directory is based on YAML filename
    my $analysis_dir = $analysis_yaml;
    $analysis_dir =~ s/[.]yaml \z//xms;

    # Finish command line
    my $cmd_line = $base_cmd_line;
    $cmd_line .=
      ' --dir ' . $analysis_dir . ' --analysis_yaml ' . $analysis_yaml;

    # Create pipeline
    my $pipeline = DETCT::Pipeline::WithDiffExprStages->new(
        {
            scheduler    => $scheduler,
            analysis_dir => $analysis_dir,
            analysis     => $analysis,
            cmd_line     => $cmd_line,
            max_retries  => $max_retries,
            sleep_time   => $sleep_time,
            verbose      => $verbose,
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
    push @log, "\nYAML analysis config file:\n\n", read_file($analysis_yaml);
    push @log, "\nYAML stages config file:\n\n",   read_file($stages_yaml);
    $pipeline->write_log_file( 'de.log', @log );

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
    $program_name =~ s/run_multiple_de_pipelines/run_de_pipeline/xms;

    return join q{ }, Probe::Perl->find_perl_interpreter(), @libs,
      $program_name;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'scheduler=s'          => \$scheduler,
        'analysis_yamls=s{1,}' => \@analysis_yamls,
        'stages_yaml=s'        => \$stages_yaml,
        'avoid_nodes=s@{,}'    => \@avoid_nodes,
        'max_retries=i'        => \$max_retries,
        'sleep_time=i'         => \$sleep_time,
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

    # Check options
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

    run_multiple_de_pipelines.pl
        [--scheduler lsf|local]
        [--analysis_yamls file...]
        [--stages_yaml file]
        [--avoid_nodes node...]
        [--max_retries int]
        [--sleep_time int]
        [--stage stage]
        [--component int]
        [--verbose]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--scheduler lsf|local>

Job scheduler - lsf (default) or local (for testing).

=item B<--analysis_yamls FILE...>

YAML analysis configuration files.

=item B<--stages_yaml FILE>

YAML stages configuration file.

=item B<--avoid_nodes NODE...>

Nodes to be avoided when submitting LSF jobs.

=item B<--max_retries INT>

Maximum number of times to retry a failing job.

=item B<--sleep_time INT>

Time to sleep, in seconds, between each cycle of the pipeline.

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
