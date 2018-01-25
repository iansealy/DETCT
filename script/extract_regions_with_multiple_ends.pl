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

my $is_required_gene = check_regions($input_file);
output_regions( $input_file, $is_required_gene );

# Check regions to identify those with multiple ends
sub check_regions {
    my ($file) = @_;

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %region_count_for;
    my %end_count_for;

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        ## no critic (ProhibitMagicNumbers)
        # Get region ID by joining chr, start, end and strand
        my $region              = join q{:}, @fields[ 0, 1, 2, 4 ];
        my $genes               = $fields[9];
        my @three_prime_end_pos = split /,/xms, $fields[3];
        ## use critic
        next if !$genes;
        foreach my $gene ( split /,/xms, $genes ) {
            $region_count_for{$gene}++;
            foreach my $three_prime_end_pos (@three_prime_end_pos) {

                # Get 3' end ID by join chr and end position
                my $three_prime_end = join q{:}, $fields[0],
                  $three_prime_end_pos;
                $end_count_for{$gene}{$region}{$three_prime_end} = 1;
            }
        }
    }
    close $fh;

    my %is_required_gene;
  GENE: foreach my $gene ( keys %region_count_for ) {
        next if $region_count_for{$gene} == 1;
        foreach my $region ( keys %{ $end_count_for{$gene} } ) {
            next GENE
              if scalar keys %{ $end_count_for{$gene}{$region} } ==
              $region_count_for{$gene};
        }
        $is_required_gene{$gene} = 1;
    }

    return \%is_required_gene;
}

# Output regions of interest
sub output_regions {
    my ( $file, $is_required_gene ) = @_;    ## no critic (ProhibitReusedNames)

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    open my $fh, '<', $file;                 ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    print_or_die($header);
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        my $genes = $fields[9];              ## no critic (ProhibitMagicNumbers)
        next if !$genes;
        foreach my $gene ( split /,/xms, $genes ) {
            if ( exists $is_required_gene->{$gene} ) {
                print_or_die($line);
                last;
            }
        }
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
