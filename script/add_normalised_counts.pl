#!/usr/bin/env perl

# PODNAME: add_normalised_counts.pl
# ABSTRACT: Add normalised counts to DESeq input file

## Author         : is1
## Maintainer     : is1
## Created        : 2014-07-04
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
use Path::Tiny;
use DETCT::Misc qw( print_or_die );

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Default options
my $input_file;
my $size_factors_file;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Read in size factors
my @size_factors = path($size_factors_file)->lines( { chomp => 1 } );

open my $fh, '<', $input_file;
print_header($fh);
print_rows( $fh, \@size_factors );
close $fh;

# Print header with additional columns for normalised counts
sub print_header {
    my ($fh) = @_;

    my $header = <$fh>;
    chomp $header;
    my @columns = split /\t/xms, $header;
    shift @columns;

    my @normalised_columns = map { $_ . ' normalised' } @columns;

    print_or_die( ( join "\t", q{}, @columns, @normalised_columns ), "\n" );

    return;
}

# Print rows with the addition of normalised counts
sub print_rows {
    my ( $fh, $size_factors ) = @_;

    while ( my $line = <$fh> ) {
        chomp $line;
        my @counts = split /\t/xms, $line;
        my $id = shift @counts;
        my @normalised_counts;
        foreach my $i ( 0 .. ( scalar @{$size_factors} ) - 1 ) {
            push @normalised_counts, $counts[$i] / $size_factors->[$i];
        }
        print_or_die( ( join "\t", $id, @counts, @normalised_counts ), "\n" );
    }

    return;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'input_file=s'        => \$input_file,
        'size_factors_file=s' => \$size_factors_file,
        'help'                => \$help,
        'man'                 => \$man,
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
    if ( !$size_factors_file ) {
        pod2usage("--size_factors_file must be specified\n");
    }

    return;
}

=head1 USAGE

    add_normalised_counts.pl
        [--input_file file]
        [--size_factors_file file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--input_file FILE>

DESeq input file.

=item B<--size_factors_file FILE>

Size factors file.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
