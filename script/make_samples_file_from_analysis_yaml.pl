#!/usr/bin/env perl

# PODNAME: make_samples_file_from_analysis_yaml.pl
# ABSTRACT: Make DESeq2 samples.txt file from analysis YAML file

## Author         : is1
## Maintainer     : is1
## Created        : 2017-08-30
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
use DETCT::Sample;
use DETCT::Misc::R;
use DETCT::Misc qw( print_or_die );

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $analysis_yaml = 'analysis.yaml';
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

# Get samples
my @samples;
foreach my $sample ( @{ $yaml_tiny->[0]->{'samples'} } ) {
    if ( !exists $sample->{'group'} ) {
        $sample->{'group'} = undef;
    }
    push @samples,
      DETCT::Sample->new(
        {
            name        => $sample->{'name'},
            description => $sample->{'description'},
            condition   => $sample->{'condition'},
            group       => $sample->{'group'},
            tag         => $sample->{'tag'},
            bam_file    => $sample->{'bam_file'},
        }
      );
}

# Print samples.txt
my $condition_prefix = DETCT::Misc::R::condition_prefix( \@samples );
my $group_prefix     = DETCT::Misc::R::group_prefix( \@samples );
my $samples_text;
foreach my $sample (@samples) {
    my @row = (
        $sample->name,
        $condition_prefix . $sample->condition,
        map { $group_prefix . $_ } @{ $sample->groups },
    );
    $samples_text .= ( join "\t", @row ) . "\n";
}
my @header     = ( q{}, 'condition' );
my $num_groups = scalar @{ $samples[0]->groups };
if ( $num_groups == 1 ) {
    push @header, 'group';
}
elsif ( $num_groups > 1 ) {
    push @header, map { 'group' . $_ } ( 1 .. $num_groups );
}
print_or_die( ( join "\t", @header ), "\n" );
print_or_die($samples_text);

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'analysis_yaml=s' => \$analysis_yaml,
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

    return;
}

=head1 USAGE

    make_samples_file_from_analysis_yaml.pl
        [--analysis_yaml file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--analysis_yaml FILE>

Analysis YAML file.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
