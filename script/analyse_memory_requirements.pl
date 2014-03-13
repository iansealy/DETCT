#!/usr/bin/env perl

# PODNAME: analyse_memory_requirements.pl
# ABSTRACT: Analyse memory requirements of pipeline based on previous runs

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-11
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Getopt::Long;
use Pod::Usage;
use File::Spec;
use YAML::Tiny;
use Statistics::Descriptive;

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $stages_yaml = 'stages.yaml';
my @analysis_dirs;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Read stages config
confess "YAML file ($stages_yaml) does not exist or cannot be read"
  if !-r $stages_yaml;
my $yaml = YAML::Tiny->read($stages_yaml);
if ( !$yaml ) {
    confess sprintf 'YAML file (%s) is invalid: %s', $stages_yaml,
      YAML::Tiny->errstr;
}

# Iterate over each stage
foreach my $stage_hash ( @{ $yaml->[0] } ) {
    my $stage_name  = $stage_hash->{name};
    my $default_mem = $stage_hash->{default_memory};
    printf "%s:\n  Current default memory: %d MB\n", $stage_name, $default_mem;

    # Get all previous memory usage
    my @mem_usage;
    foreach my $dir (@analysis_dirs) {
        next if !-d $dir;

        my $mem_file = File::Spec->catfile( $dir, $stage_name . '.memory.txt' );
        next if !-r $mem_file;
        my $mem_yaml = YAML::Tiny->read($mem_file);
        next if !$mem_yaml;
        push @mem_usage, @{ $mem_yaml->[0] };
    }

    # Show stats
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@mem_usage);
    next if !$stat->count;
    printf "  Job count:\t%d\n",          $stat->count;
    printf "  Min memory:\t%d\n",         $stat->min;
    printf "  Median memory:\t%d\n",      $stat->median;
    printf "  Mean memory:\t%d\n",        $stat->mean;
    printf "  Max memory:\t%d\n",         $stat->max;
    printf "  Standard deviation:\t%d\n", $stat->standard_deviation;

    # Calculate minimum memory wasted for each different max retries
    my %optimal_mem;
    foreach my $test_mem ( reverse $stat->min .. $stat->max ) {
        next if $test_mem == 1;    # Won't increase
        my ( $mem_wasted, $retries ) = simulate_stage( $test_mem, \@mem_usage );
        $optimal_mem{$retries}{$mem_wasted} = $test_mem;
    }

    foreach my $max_retries ( sort { $a <=> $b } keys %optimal_mem ) {
        my $min_mem_wasted =
          ( sort { $a <=> $b } keys %{ $optimal_mem{$max_retries} } )[0];
        my $optimal_mem = $optimal_mem{$max_retries}{$min_mem_wasted};
        printf "    Max retries: %d\tOptimal memory: %d\tMemory wasted: %d\n",
          $max_retries, $optimal_mem, $min_mem_wasted;
    }
}

# For a particular default memory, simulate number of retries and calculate
# memory wasted
sub simulate_stage {
    my ( $initial_mem, $mem_usage ) = @_;

    my $mem_wasted  = 0;
    my $max_retries = 0;

    foreach my $actual_mem ( @{$mem_usage} ) {
        my $mem     = $initial_mem;
        my $retries = 0;
        while ( $mem < $actual_mem ) {
            $mem_wasted += $mem;
            $retries++;
            $mem = int( $mem * 1.5 );    ## no critic (ProhibitMagicNumbers)
        }
        $mem_wasted += $mem - $actual_mem;
        if ( $retries > $max_retries ) {
            $max_retries = $retries;
        }
    }

    return $mem_wasted, $max_retries;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'stages_yaml=s' => \$stages_yaml,
        'dirs=s@{,}'    => \@analysis_dirs,
        'help'          => \$help,
        'man'           => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

=head1 USAGE

    analyse_memory_requirements.pl
        [--stages_yaml file]
        [--dirs directory...]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--stages_yaml FILE>

YAML stages configuration file.

=item B<--dirs DIRECTORY...>

Working directories for analysis.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
