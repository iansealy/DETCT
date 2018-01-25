#!/usr/bin/env perl

# PODNAME: extract_regions_with_multiple_ends.pl
# ABSTRACT: Extract regions that have multiple 3' ends

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
my $input_file;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my $is_required_region = check_regions($input_file);
output_regions( $input_file, $is_required_region );

# Check regions to identify those with multiple ends
sub check_regions {
    my ($file) = @_;

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %gene_count;
    my %end_count;

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID by joining chr, start, end and strand
        ## no critic (ProhibitMagicNumbers)
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        $gene_count{ $fields[9] }{$region} = 1;
        foreach my $three_prime_end_pos ( split /,/xms, $fields[3] ) {
            ## use critic
            my $three_prime_end = join q{:}, $fields[0], $three_prime_end_pos;
            $end_count{$three_prime_end}{$region} = 1;
        }
    }
    close $fh;

    my %is_required_region;
    foreach my $gene ( keys %gene_count ) {
        if ( scalar keys %{ $gene_count{$gene} } > 1 ) {
            foreach my $region ( keys %{ $gene_count{$gene} } ) {
                $is_required_region{$region} = 1;
            }
        }
    }
    foreach my $end ( keys %end_count ) {
        if ( scalar keys %{ $end_count{$end} } > 1 ) {
            foreach my $region ( keys %{ $end_count{$end} } ) {
                $is_required_region{$region} = 1;
            }
        }
    }

    return \%is_required_region;
}

# Output regions of interest
sub output_regions {
    my ( $file, $is_required_region ) = @_;   ## no critic (ProhibitReusedNames)

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    open my $fh, '<', $file;                  ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    print_or_die($header);
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID by joining chr, start, end and strand
        ## no critic (ProhibitMagicNumbers)
        my $region = join q{:}, @fields[ 0, 1, 2, 4 ];
        ## use critic
        next if !exists $is_required_region->{$region};
        print_or_die($line);
    }
    close $fh;

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'input_file=s' => \$input_file,
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
    if ( !$input_file ) {
        pod2usage("--input_file must be specified\n");
    }

    return;
}

=head1 USAGE

    extract_regions_with_multiple_ends.pl
        [--input_file file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--input_file FILE>

Differential expression pipeline output file (e.g. all.tsv).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
