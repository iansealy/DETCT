#!/usr/bin/env perl

# PODNAME: make_all_analysis_yaml_conditions_pairs.pl
# ABSTRACT: Create analysis YAML files with all pairs of conditions

## Author         : is1
## Maintainer     : is1
## Created        : 2015-11-25
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
use YAML::Tiny;
use File::Slurp;
use DETCT::Misc qw( write_or_die );

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $base_yaml = 'analysis.yaml';
my $table_file;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Read analysis config
confess "YAML file ($base_yaml) does not exist or cannot be read"
  if !-r $base_yaml;
my $yaml      = read_file($base_yaml);
my $yaml_tiny = YAML::Tiny->read_string($yaml);
if ( !$yaml_tiny ) {
    confess sprintf 'YAML file (%s) is invalid: %s', $base_yaml,
      YAML::Tiny->errstr;
}

# Get all conditions
my %is_condition;
foreach my $sample ( @{ $yaml_tiny->[0]->{'samples'} } ) {
    $is_condition{ $sample->{'condition'} } = 1;
}
my @conditions = sort keys %is_condition;

# Get all pairs
my @pairs;
while ( scalar @conditions > 1 ) {
    my $condition1 = shift @conditions;
    foreach my $condition2 (@conditions) {
        push @pairs, [ $condition1, $condition2 ];
        push @pairs, [ $condition2, $condition1 ];
    }
}

# Write YAML files
foreach my $pair (@pairs) {
    my $out_file = $pair->[0] . '_vs_' . $pair->[1] . '.yaml';
    my $output   = $yaml;
    if ($table_file) {
        $output .= sprintf "table_file: %s\n", $table_file;
    }
    $output .= sprintf "control_condition: %s\n",      $pair->[0];
    $output .= sprintf "experimental_condition: %s\n", $pair->[1];
    open my $fh, '>', $out_file;
    write_or_die( $fh, $output );
    close $fh;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'base_yaml=s'  => \$base_yaml,
        'table_file=s' => \$table_file,
        'help'         => \$help,
        'man'          => \$man,
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

    make_all_analysis_yaml_conditions_pairs.pl
        [--base_yaml file]
        [--table_file file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--base_yaml FILE>

YAML file to use as base for new YAML files (and to extract conditions from).

=item B<--table_file FILE>

Differential expression pipeline CSV or TSV output file.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
