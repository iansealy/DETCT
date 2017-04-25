#!/usr/bin/env perl

# PODNAME: redump.pl
# ABSTRACT: Dump data from output file

## Author         : is1
## Maintainer     : is1
## Created        : 2017-04-18
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
use Readonly;
use Path::Tiny;
use File::Spec;
use File::Path qw( make_path );
use DETCT::Analysis::DiffExpr;
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES

    perl redump.pl \
        --dir fixed \
        --analysis_yaml analysis.yaml \
        --all_file all.tsv

=cut

# Default options
my $analysis_dir = q{.};
my $analysis_yaml = File::Spec->catfile( $analysis_dir, 'analysis.yaml' );
my $all_file;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

if ( !-d $analysis_dir ) {
    make_path($analysis_dir);
}

# Create analysis
my $analysis = DETCT::Analysis::DiffExpr->new_from_yaml($analysis_yaml);

my $regions = get_regions( $all_file, $analysis );

DETCT::Misc::Output::dump_as_table(
    {
        analysis => $analysis,
        dir      => $analysis_dir,
        regions  => $regions,
    }
);

sub get_regions {
    ## no critic (ProhibitReusedNames)
    my ( $file, $analysis ) = @_;
    ## use critic

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    ## no critic (ProhibitReusedNames)
    my $regions = DETCT::Misc::Output::parse_table(
        ## use critic
        {
            analysis     => $analysis,
            table_file   => $file,
            table_format => $extension,
        }
    );

    return $regions;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'dir=s'           => \$analysis_dir,
        'analysis_yaml=s' => \$analysis_yaml,
        'all_file=s'      => \$all_file,
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
    if ( !$all_file ) {
        pod2usage("--all_file must be specified\n");
    }

    return;
}

=head1 USAGE

    redump.pl
        [--dir directory]
        [--analysis_yaml file]
        [--all_file file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--dir DIRECTORY>

Output directory.

=item B<--analysis_yaml FILE>

YAML analysis configuration file.

=item B<--all_file FILE>

DETCT output file (containing all regions).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
