#!/usr/bin/env perl

# PODNAME: gather_npg_qc_data.pl
# ABSTRACT: Gather NPG data for DETCT experiments for QC

## Author         : is1
## Maintainer     : is1
## Created        : 2016-05-04
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
use DBI;
use Sort::Naturally;
use DETCT::Misc qw( printf_or_die );

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Constants
Readonly our $HOST => 'seqw-db';
Readonly our $PORT => '3379';
Readonly our $USER => 'warehouse_ro';
Readonly our $PASS => q{};
Readonly our $NAME => 'sequencescape_warehouse';

Readonly our @COMMON_NAMES =>
  qw( tag_index tag_sequence insert_size_quartile1 insert_size_median insert_size_quartile3 gc_percent_forward_read gc_percent_reverse_read sequence_mismatch_percent_forward_read sequence_mismatch_percent_reverse_read adapters_percent_forward_read adapters_percent_reverse_read tag_decode_percent tag_decode_count bam_num_reads bam_percent_mapped bam_percent_duplicate accession_number );
Readonly our $OUTPUT_NAMES => join "\t",
  ( qw( sample run lane ), @COMMON_NAMES );
Readonly our $FIELD_NAMES => join q{, },
  ( qw( id_run position ), @COMMON_NAMES );

# Default options
my @expts;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Connect to database
my $dsn = "dbi:mysql:database=$NAME;host=$HOST;port=$PORT";
my $dbh = DBI->connect( $dsn, $USER, $PASS );

printf_or_die( "%s\n", $OUTPUT_NAMES );

foreach my $expt (@expts) {
    my $expt_q = $dbh->quote($expt);
    $expt_q =~ s/'\z/%'/xms;
    my $ary_ref = $dbh->selectall_arrayref(
        <<"SQL"
        SELECT name, supplier_name
        FROM   current_samples
        WHERE  (name LIKE $expt_q OR supplier_name LIKE $expt_q)
        AND    description LIKE '%polyT%'
SQL
    );
    my @samples;
    foreach ( @{$ary_ref} ) {
        my ( $name, $supplier_name ) = @{$_};
        if ( $name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $name;
        }
        elsif ( $supplier_name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $supplier_name;
        }
    }

    foreach my $sample ( sort { ncmp( $a, $b ) } @samples ) {
        my $sample_q = $dbh->quote($sample);
        $ary_ref = $dbh->selectall_arrayref(
            <<"SQL"
            SELECT   $FIELD_NAMES
            FROM     npg_plex_information npi, current_samples cs
            WHERE    npi.sample_id = cs.internal_id
            AND      (supplier_name = $sample_q OR name = $sample_q)
            ORDER BY id_run, position
SQL
        );
        foreach ( @{$ary_ref} ) {
            my @values = map { !defined $_ ? q{} : $_ } @{$_};
            printf_or_die( "%s\n", join "\t", $sample, @values );
        }
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'expts=s@{1,}' => \@expts,
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
    if ( !@expts ) {
        pod2usage("--expts must be specified\n");
    }

    return;
}

=head1 USAGE

    gather_npg_qc_data.pl
        [--expts prefixes]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--expts PREFIXES>

Experiment prefixes (e.g. zmp_ph80).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
