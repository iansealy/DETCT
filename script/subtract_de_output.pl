#!/usr/bin/env perl

# PODNAME: subtract_de_output.pl
# ABSTRACT: Subtract differential expression pipeline output file from another

## Author         : is1
## Maintainer     : is1
## Created        : 2015-10-22
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
use DETCT::Misc qw( print_or_die );
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Default options
my $full_file;
my $minus_file;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my $is_minus_region = get_regions($minus_file);
output_regions( $full_file, $is_minus_region );

# Get regions of interest
sub get_regions {
    my ($file) = @_;

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %is_region;

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID by joining chr, start, end, 3' end position and strand
        ## no critic (ProhibitMagicNumbers)
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        ## use critic
        $is_region{$region} = 1;
    }
    close $fh;

    return \%is_region;
}

# Output regions of interest with new counts
sub output_regions {
    my ( $file, $is_minus_region ) = @_;    ## no critic (ProhibitReusedNames)

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    open my $fh, '<', $file;                ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    print_or_die($header);
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID by joining chr, start, end, 3' end position and strand
        ## no critic (ProhibitMagicNumbers)
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        ## use critic
        next if exists $is_minus_region->{$region};
        print_or_die($line);
    }
    close $fh;

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'full_file=s'  => \$full_file,
        'minus_file=s' => \$minus_file,
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

    # Check options
    if ( !$full_file ) {
        pod2usage("--full_file must be specified\n");
    }
    if ( !$minus_file ) {
        pod2usage("--minus_file must be specified\n");
    }

    return;
}

=head1 USAGE

    subtract_de_output.pl
        [--full_file file]
        [--minus_file file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--full_file FILE>

Differential expression pipeline output file to have regions removed from (e.g.
all.tsv).

=item B<--minus_file FILE>

Differential expression pipeline output file to be subtracted (by region) from
full file (e.g. sig.tsv).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
