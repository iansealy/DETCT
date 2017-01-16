#!/usr/bin/env perl

# PODNAME: subset_analysis_yaml.pl
# ABSTRACT: Subset an analysis YAML file

## Author         : is1
## Maintainer     : is1
## Created        : 2017-01-16
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
use Path::Tiny;
use DETCT::Misc qw( print_or_die );

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $base_yaml = 'analysis.yaml';
my $name;
my $control_condition;
my $experimental_condition;
my @conditions;
my @groups;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Read base analysis config
confess "YAML file ($base_yaml) does not exist or cannot be read"
  if !-r $base_yaml;
my $yaml      = path($base_yaml)->slurp;
my $yaml_tiny = YAML::Tiny->read_string($yaml);
if ( !$yaml_tiny ) {
    confess sprintf 'YAML file (%s) is invalid: %s', $base_yaml,
      YAML::Tiny->errstr;
}

# Add analysis name and conditions if specified
if ($name) {
    $yaml_tiny->[0]->{'name'} = $name;
}
if ($control_condition) {
    $yaml_tiny->[0]->{'control_condition'}      = $control_condition;
    $yaml_tiny->[0]->{'experimental_condition'} = $experimental_condition;
}

# Get all required conditions and groups
my %is_required_condition = map { $_ => 1 } @conditions;
my %is_required_group     = map { $_ => 1 } @groups;
if ($control_condition) {
    $is_required_condition{$control_condition} = 1;
}
if ($experimental_condition) {
    $is_required_condition{$experimental_condition} = 1;
}

# Subset samples
my $subset_on_condition = scalar keys %is_required_condition ? 1 : 0;
my $subset_on_group     = scalar keys %is_required_group     ? 1 : 0;
my @subset;
foreach my $sample ( @{ $yaml_tiny->[0]->{'samples'} } ) {
    next
      if $subset_on_condition
      && !exists $is_required_condition{ $sample->{'condition'} };
    next
      if $subset_on_group && !exists $is_required_group{ $sample->{'group'} };
    push @subset, $sample;
}
$yaml_tiny->[0]->{'samples'} = \@subset;

# Dump YAML file
$yaml = $yaml_tiny->write_string;
$yaml =~ s/\A---\n//xms;
$yaml =~ s/: \s+ '([\w.-]+)'/: $1/xmsg;    # Remove quotes round numbers
$yaml =~ s/: \s+ ~/:/xmsg;
print_or_die($yaml);

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'base_yaml=s'              => \$base_yaml,
        'name=s'                   => \$name,
        'control_condition=s'      => \$control_condition,
        'experimental_condition=s' => \$experimental_condition,
        'conditions=s{1,}'         => \@conditions,
        'groups=s{1,}'             => \@groups,
        'help'                     => \$help,
        'man'                      => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( $control_condition && !$experimental_condition ) {
        pod2usage(
            sprintf
              "--%s_condition must be specified if %s condition specified\n",
            'experimental', 'control'
        );
    }
    if ( $experimental_condition && !$control_condition ) {
        pod2usage(
            sprintf
              "--%s_condition must be specified if %s condition specified\n",
            'control', 'experimental'
        );
    }

    return;
}

=head1 USAGE

    subset_analysis_yaml.pl
        [--base_yaml file]
        [--name name]
        [--control_condition condition]
        [--experimental_condition condition]
        [--conditions condition...]
        [--groups group...]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--base_yaml FILE>

YAML file to use as base for new YAML file.

=item B<--name NAME>

The name for the analysis.

=item B<--control_condition CONDITION>

The control condition.

=item B<--experimental_condition CONDITION>

The experimental condition.

=item B<--conditions CONDITION...>

The conditions to retain.

=item B<--groups GROUP...>

The groups to retain.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
