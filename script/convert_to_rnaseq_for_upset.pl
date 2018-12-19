#!/usr/bin/env perl

# PODNAME: convert_to_rnaseq_for_upset.pl
# ABSTRACT: Convert output files to RNA-Seq format for UpSet

## Author         : is1
## Maintainer     : is1
## Created        : 2017-05-04
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
my $file;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
    confess sprintf '%s is not .csv or .tsv file', $file;
}
open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
my $line     = <$fh>;
my @headings = DETCT::Misc::Output::parse_line( $line, $extension );
my ($log2fc_col) =
  grep { $headings[$_] =~ m/\A Log2 \s fold \s change \s [(] .* [)] \z/xms }
  ( 0 .. scalar @headings - 1 );
splice @headings, 0, 15;    ## no critic (ProhibitMagicNumbers)
@headings = grep { !m/\A Log2 \s fold \s change/xms } @headings;
@headings = (
    'Gene ID', 'pval',  'adjp',        'log2fc',
    'Chr',     'Start', 'End',         'Strand',
    'Biotype', 'Name',  'Description', @headings,
);
print_or_die( ( join "\t", @headings ), "\n" );

while ( $line = <$fh> ) {
    my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

    # Ignore line if no gene annotation
    next if !defined $fields[9];    ## no critic (ProhibitMagicNumbers)

    # Note:
    # Output can contain multiple rows per gene (unlike RNA-Seq output)
    # Start and end correspond to region, not gene
    # Regions associated with multiple genes are duplicated on multiple rows
    ## no critic (ProhibitMagicNumbers)
    my @output = (
        @fields[ 6 .. 7 ],
        $fields[$log2fc_col],
        @fields[ 0 .. 2 ],
        $fields[4],
        $fields[10],
        @fields[ 13 .. 14 ],
        @fields[ 15 .. ( $log2fc_col - 1 ) ],
    );
    foreach my $gene_id ( split /,/xms, $fields[9] ) {
        ## use critic
        print_or_die( ( join "\t", $gene_id, @output ), "\n" );
    }
}
close $fh;

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'file=s' => \$file,
        'help'   => \$help,
        'man'    => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$file ) {
        pod2usage("--file must be specified\n");
    }

    return;
}

=head1 USAGE

    convert_to_rnaseq_for_upset.pl
        [--file file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--file FILE>

Differential expression pipeline output file (e.g. sig.tsv).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
