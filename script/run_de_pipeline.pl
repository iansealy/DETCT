#!/usr/bin/env perl

# PODNAME: run_de_pipeline.pl
# ABSTRACT: Run DETCT differential expression pipeline

## Author         : is1
## Maintainer     : is1
## Created        : 2012-09-26
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
my $scheduler     = 'lsf';
my $analysis_dir  = q{.};
my $analysis_yaml = File::Spec->catfile( $analysis_dir, 'analysis.yaml' );
my $stages_yaml   = File::Spec->catfile( $analysis_dir, 'stages.yaml' );
my @avoid_nodes;
## no critic (ProhibitMagicNumbers)
my $max_retries = 10;
my $sleep_time  = 600;    # 10 minutes
## use critic
my $stage_to_run;
my $component_to_run;
my $verbose;
my ( $help, $man );

# Get command line (including interpreter and options)
my $cmd_line = get_cmd_line();

# Get and check command line options
get_and_check_options();

# Create analysis
my $analysis = DETCT::Analysis->new_from_yaml($analysis_yaml);

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

# Are we running the main pipeline or running a specific component of a specific
# stage (i.e. a job to be run under LSF or locally)?
if ($stage_to_run) {
    $pipeline->set_stage_to_run( $pipeline->get_stage_by_name($stage_to_run) );
}
if ($component_to_run) {
    $pipeline->set_component_to_run($component_to_run);
}

# Turn off verbose output when running specific components
if ( $pipeline->stage_to_run && $pipeline->component_to_run ) {
    $pipeline->set_verbose(0);
}

# Write overview of pipeline input and config files to log file
if ( !$pipeline->stage_to_run && !$pipeline->component_to_run ) {
    my @log = map { "$_\n" } $pipeline->input_overview;
    push @log, "\nYAML analysis config file:\n\n", read_file($analysis_yaml);
    push @log, "\nYAML stages config file:\n\n",   read_file($stages_yaml);
    $pipeline->write_log_file( 'de.log', @log );
}

# Print overview of pipeline input
$pipeline->say_if_verbose( $pipeline->input_overview );

# Run pipeline
$pipeline->run();

# Get entire command line
sub get_cmd_line {

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

    return join q{ }, Probe::Perl->find_perl_interpreter(), @libs,
      $PROGRAM_NAME, @ARGV;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'scheduler=s'       => \$scheduler,
        'dir=s'             => \$analysis_dir,
        'analysis_yaml=s'   => \$analysis_yaml,
        'stages_yaml=s'     => \$stages_yaml,
        'avoid_nodes=s@{,}' => \@avoid_nodes,
        'max_retries=i'     => \$max_retries,
        'sleep_time=i'      => \$sleep_time,
        'stage=s'           => \$stage_to_run,
        'component=i'       => \$component_to_run,
        'verbose'           => \$verbose,
        'help'              => \$help,
        'man'               => \$man,
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

    run_de_pipeline.pl
        [--scheduler lsf|local]
        [--dir directory]
        [--analysis_yaml file]
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

=item B<--dir DIRECTORY>

Working directory for analysis.

=item B<--analysis_yaml FILE>

YAML analysis configuration file.

=item B<--stages_yaml FILE>

YAML stages configuration file.

=item B<--avoid_nodes NODE...>

Nodes to be avoided when submitting LSF jobs.

=item B<--max_retries INT>

Maximum number of times to retry a failing job.

=item B<--sleep_time INT>

Time to sleep, in seconds, between each iteration of the pipeline.

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
