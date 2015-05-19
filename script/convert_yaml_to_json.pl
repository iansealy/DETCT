#!/usr/bin/env perl

# PODNAME: convert_yaml_to_json.pl
# ABSTRACT: Convert YAML output files to JSON

## Author         : is1
## Maintainer     : is1
## Created        : 2015-05-19
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
use File::Slurp;
use YAML;
use JSON;

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my @files;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

foreach my $file (@files) {
    my $data = YAML::LoadFile($file);
    rename $file, $file . '.orig';
    write_file( $file, JSON::to_json( $data, { pretty => 1 } ) );
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'files=s@{,}' => \@files,
        'help'        => \$help,
        'man'         => \$man,
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

    convert_yaml_to_json.pl
        [--files file...]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--files FILE...>

YAML files to convert to JSON.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
