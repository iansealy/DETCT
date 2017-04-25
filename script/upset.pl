#!/usr/bin/env perl

# PODNAME: upset.pl
# ABSTRACT: Convert output files to UpSet format

## Author         : is1
## Maintainer     : is1
## Created        : 2017-04-25
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
use Sort::Versions;
use DETCT::Misc qw( printf_or_die );
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Default options
my @files;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my @sets = remove_common_prefix_suffix(@files);

my %upset;

my $i = 0;
foreach my $file (@files) {
    my $set = $sets[$i];    ## no critic (ProhibitAmbiguousNames)
    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }
    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID by joining chr, start, end, strand and optional gene ID
        ## no critic (ProhibitMagicNumbers)
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        if ( $fields[9] ) {
            $region .= q{:} . $fields[9];
        }
        ## use critic
        $upset{$region}{$set} = 1;
    }
    close $fh;
    $i++;
}

printf_or_die( "%s\r\n", join q{,}, 'Region', @sets );
foreach my $region ( sort { versioncmp( $a, $b ) } keys %upset ) {
    printf_or_die( "%s\r\n", join q{,}, $region,
        map { $upset{$region}{$_} || 0 } @sets );
}

sub remove_common_prefix_suffix {
    my @strings = @_;

    foreach ( 0 .. 1 ) {
        @strings = map { scalar reverse } @strings;
        my $common_prefix_idx = 0;
      CHAR: while (1) {
            my $char = substr $strings[0], $common_prefix_idx, 1;
            foreach my $string (@strings) {
                last CHAR if substr( $string, $common_prefix_idx, 1 ) ne $char;
            }
            $common_prefix_idx++;
        }
        if ($common_prefix_idx) {
            @strings = map { substr $_, $common_prefix_idx } @strings;
        }
    }

    return @strings;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'files=s@{2,}' => \@files,
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
    if ( !@files ) {
        pod2usage("--files must be specified\n");
    }

    return;
}

=head1 USAGE

    upset.pl
        [--files file...]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--files FILE...>

Differential expression pipeline output files (e.g. sig.tsv).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
