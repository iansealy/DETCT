#!/usr/bin/env perl

# PODNAME: remove_sample_from_analysis_yaml.pl
# ABSTRACT: Remove a sample from an analysis YAML file

## Author         : is1
## Maintainer     : is1
## Created        : 2017-03-02
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
my $analysis_yaml = 'analysis.yaml';
my $sample_name;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Read analysis config
confess "YAML file ($analysis_yaml) does not exist or cannot be read"
  if !-r $analysis_yaml;
my $yaml      = path($analysis_yaml)->slurp;
my $yaml_tiny = YAML::Tiny->read_string($yaml);
if ( !$yaml_tiny ) {
    confess sprintf 'YAML file (%s) is invalid: %s', $analysis_yaml,
      YAML::Tiny->errstr;
}

# Subset samples
my @subset;
foreach my $sample ( @{ $yaml_tiny->[0]->{'samples'} } ) {
    next if $sample->{'name'} eq $sample_name;
    push @subset, $sample;
}
$yaml_tiny->[0]->{'samples'} = \@subset;

# Dump YAML file
$yaml = $yaml_tiny->write_string;
$yaml =~ s/\A---\n//xms;
$yaml =~ s/: \s+ '([^']+)'/: $1/xmsg;    # Remove quotes round values
$yaml =~ s/: \s+ ~/:/xmsg;
print_or_die($yaml);

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'analysis_yaml=s' => \$analysis_yaml,
        'sample_name=s'   => \$sample_name,
        'help'            => \$help,
        'man'             => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$sample_name ) {
        pod2usage("--sample_name must be specified\n");
    }

    return;
}

=head1 USAGE

    remove_sample_from_analysis_yaml.pl
        [--analysis_yaml file]
        [--sample_name name]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--analysis_yaml FILE>

Analysis YAML file.

=item B<--sample_name NAME>

The sample name to remove.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
