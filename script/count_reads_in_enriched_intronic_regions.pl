#!/usr/bin/env perl

# PODNAME: count_reads_in_enriched_intronic_regions.pl
# ABSTRACT: Count reads in base-enriched intronic regions.

## Author         : is1
## Maintainer     : is1
## Created        : 2016-03-14
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
use DETCT::Misc qw( printf_or_die );
use Bio::DB::HTS;

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $bam_file;
my $region_file;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my $hts = Bio::DB::HTS->new( -bam => $bam_file );

open my $fh, '<', $region_file;    ## no critic (RequireBriefOpen)
while ( my $line = <$fh> ) {
    chomp $line;
    my ( undef, $chromosome, $strand, $start, $end ) = split /\t/xms, $line;
    my $read2_alignments = $hts->features(
        -seq_id   => $chromosome,
        -start    => $start,
        -end      => $end,
        -flags    => { SECOND_MATE => 1 },
        -iterator => 1,
    );
    my $forward_strand_count = 0;
    my $reverse_strand_count = 0;

    while ( my $alignment = $read2_alignments->next_seq ) {
        next if $alignment->get_tag_values('FLAGS') =~ m/\bDUPLICATE\b/xms;
        if ( $alignment->strand == $strand ) {
            $forward_strand_count++;
        }
        else {
            $reverse_strand_count++;
        }
    }

    printf_or_die( "%s\t%d\t%d\n", $line, $forward_strand_count,
        $reverse_strand_count );
}

close $fh;

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'bam_file=s'    => \$bam_file,
        'region_file=s' => \$region_file,
        'help'          => \$help,
        'man'           => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$bam_file ) {
        pod2usage("--bam_file must be specified\n");
    }
    if ( !$region_file ) {
        pod2usage("--region_file must be specified\n");
    }

    return;
}

=head1 USAGE

    count_reads_in_enriched_intronic_regions.pl
        [--bam_file file]
        [--region_file file]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--bam_file FILE>

BAM file.

=item B<--region_file FILE>

Region file (from get_base_enriched_intronic_regions.pl).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
