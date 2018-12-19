#!/usr/bin/env perl

# PODNAME: subset_de_output.pl
# ABSTRACT: Subset output file from differential expression pipeline

## Author         : is1
## Maintainer     : is1
## Created        : 2015-06-16
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
my $regions_file;
my $counts_file;
my $ignore_missing_regions;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my $is_region = get_regions($regions_file);
my ( $counts_for, $counts_headings ) =
  get_counts_by_region( $counts_file, $is_region );
output_regions( $regions_file, $counts_for, $counts_headings );

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

        # Get region ID by joining chr, start, end and strand
        ## no critic (ProhibitMagicNumbers)
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        ## use critic
        $is_region{$region} = 1;
    }
    close $fh;

    return \%is_region;
}

# Get counts (raw and normalised) by region
sub get_counts_by_region {
    my ( $file, $is_region ) = @_;    ## no critic (ProhibitReusedNames)

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %counts_for;

    open my $fh, '<', $file;          ## no critic (RequireBriefOpen)
                                      # Get range for count columns
    my $header   = <$fh>;
    my @headings = DETCT::Misc::Output::parse_line( $header, $extension );
    my $counts_start = 15;            ## no critic (ProhibitMagicNumbers)
    my $counts_end;
    my $ordinal = $counts_start;
    while ( defined $headings[$ordinal] ) {

        if ( $headings[$ordinal] =~ m/\s count \z/xms ) {
            $counts_end = $ordinal;
        }
        $ordinal++;
    }
    @headings = @headings[ $counts_start .. $counts_end ];
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID by joining chr, start, end and strand
        ## no critic (ProhibitMagicNumbers)
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        ## use critic
        @fields = @fields[ $counts_start .. $counts_end ];
        $counts_for{$region} = \@fields;
    }
    close $fh;

    return \%counts_for, \@headings;
}

# Output regions of interest with new counts
sub output_regions {
    ## no critic (ProhibitReusedNames)
    my ( $file, $counts_for, $counts_headings ) = @_;
    ## use critic

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header   = <$fh>;
    my @headings = DETCT::Misc::Output::parse_line( $header, $extension );
    ## no critic (ProhibitMagicNumbers)
    @headings = ( @headings[ 0 .. 14 ], @{$counts_headings} );
    ## use critic
    write_line( \@headings, $extension );
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID by joining chr, start, end and strand
        ## no critic (ProhibitMagicNumbers)
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        if ( !exists $counts_for->{$region} ) {
            next if $ignore_missing_regions;
            confess "No counts for region $region";
        }
        @fields = ( @fields[ 0 .. 14 ], @{ $counts_for->{$region} } );
        ## use critic
        write_line( \@fields, $extension );
    }
    close $fh;

    return;
}

# Write line as CSV or TSV
sub write_line {
    my ( $fields, $extension ) = @_;

    if ( $extension eq 'csv' ) {
        @{$fields} = map { defined $_ ? $_ : q{} } @{$fields};
        @{$fields} = map { quote_for_csv($_) } @{$fields};
        print_or_die( ( join q{,}, @{$fields} ), "\r\n" );
    }
    elsif ( $extension eq 'tsv' ) {
        @{$fields} = map { defined $_ ? $_ : q{-} } @{$fields};
        print_or_die( ( join "\t", @{$fields} ), "\n" );
    }

    return;
}

# Escape quotes and add quotes for CSV output
sub quote_for_csv {
    my ($field) = @_;

    $field =~ s/"/""/xmsg;

    return q{"} . $field . q{"};
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'regions_file=s'         => \$regions_file,
        'counts_file=s'          => \$counts_file,
        'ignore_missing_regions' => \$ignore_missing_regions,
        'help'                   => \$help,
        'man'                    => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$regions_file ) {
        pod2usage("--regions_file must be specified\n");
    }
    if ( !$counts_file ) {
        pod2usage("--counts_file must be specified\n");
    }

    return;
}

=head1 USAGE

    subset_de_output.pl
        [--regions_file file]
        [--counts_file file]
        [--ignore_missing_regions]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--regions_file FILE>

Differential expression pipeline output file to use to select regions for
subsetting (e.g. sig.tsv).

=item B<--counts_file FILE>

Differential expression pipeline output file to use to extract counts from
(e.g. all.tsv).

=item B<--ignore_missing_regions>

Don't die if region in regions file is missing from counts file.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
