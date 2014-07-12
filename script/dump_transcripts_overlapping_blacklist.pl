#!/usr/bin/env perl

# PODNAME: dump_transcripts_overlapping_blacklist.pl
# ABSTRACT: Dump info for Ensembl transcripts that overlap a blacklist

## Author         : is1
## Maintainer     : is1
## Created        : 2014-07-11
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
use YAML::Tiny;
use Bio::EnsEMBL::Registry;
use DETCT::Misc qw( printf_or_die );

=head1 DESCRIPTION



=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-74/ensembl/modules \
        dump_transcripts_overlapping_blacklist.pl \
        --blacklist_yaml data/e74-transcript-blacklist.txt \
        > tmp/e74-dump.txt

    perl \
        -Ibranch-ensembl-75/ensembl/modules \
        dump_transcripts_overlapping_blacklist.pl \
        --blacklist_yaml data/e74-transcript-blacklist.txt \
        > tmp/e75-dump.txt

    diff tmp/e74-dump.txt tmp/e75-dump.txt

=cut

# Default options
my $blacklist_yaml;
my $species        = 'Danio rerio';
my $ensembl_dbhost = 'ensembldb.ensembl.org';
my $ensembl_dbport;
my $ensembl_dbuser = 'anonymous';
my $ensembl_dbpass;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Connnect to Ensembl database
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => $ensembl_dbhost,
    -port => $ensembl_dbport,
    -user => $ensembl_dbuser,
    -pass => $ensembl_dbpass,
);

# Get genebuild version
my $genebuild_version = 'e' . Bio::EnsEMBL::ApiVersion::software_version();
printf_or_die( "Genebuild version: %s\n", $genebuild_version );

# Get Ensembl adaptor
my $ta = Bio::EnsEMBL::Registry->get_adaptor( $species, 'core', 'Transcript' );

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Read blacklist
confess "YAML file ($blacklist_yaml) does not exist or cannot be read"
  if !-r $blacklist_yaml;
my $yaml = YAML::Tiny->read($blacklist_yaml);
if ( !$yaml ) {
    confess sprintf 'YAML file (%s) is invalid: %s', $blacklist_yaml,
      YAML::Tiny->errstr;
}

# Iterate over all blacklisted transcripts
foreach my $transcript_id ( sort @{ $yaml->[0]->{skip_transcripts} } ) {
    my $transcript = $ta->fetch_by_stable_id($transcript_id);

    # Check transcript exists
    if ( !$transcript ) {
        printf_or_die( "Blacklisted transcript %s not found\n",
            $transcript_id );
        next;
    }

    printf_or_die(
        "Blacklisted transcript %s / %s %s:%d-%d\n",
        $transcript_id,
        $transcript->get_Gene->stable_id,
        $transcript->seq_region_name,
        $transcript->seq_region_start,
        $transcript->seq_region_end
    );

    # Get overlapping transcripts
    my $slice                   = $transcript->feature_Slice;
    my $overlapping_transcripts = $slice->get_all_Transcripts();
    foreach
      my $overlapping_transcript ( sort { $a->stable_id cmp $b->stable_id }
        @{$overlapping_transcripts} )
    {
        next if $overlapping_transcript->stable_id eq $transcript_id;
        printf_or_die(
            "    %s overlaps %s / %s %s:%d-%d\n",
            $transcript->stable_id,
            $overlapping_transcript->stable_id,
            $overlapping_transcript->get_Gene->stable_id,
            $overlapping_transcript->seq_region_name,
            $overlapping_transcript->seq_region_start,
            $overlapping_transcript->seq_region_end
        );
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'blacklist_yaml=s' => \$blacklist_yaml,
        'species=s'        => \$species,
        'ensembl_dbhost=s' => \$ensembl_dbhost,
        'ensembl_dbport=i' => \$ensembl_dbport,
        'ensembl_dbuser=s' => \$ensembl_dbuser,
        'ensembl_dbpass=s' => \$ensembl_dbpass,
        'help'             => \$help,
        'man'              => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$blacklist_yaml ) {
        pod2usage("--blacklist_yaml must be specified\n");
    }

    return;
}

=head1 USAGE

    dump_transcripts_overlapping_blacklist.pl
        [--blacklist_yaml file]
        [--species species]
        [--ensembl_dbhost host]
        [--ensembl_dbport port]
        [--ensembl_dbuser username]
        [--ensembl_dbpass password]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--blacklist_yaml FILE>

Transcript blacklist in YAML format.

=item B<--species SPECIES>

Species (defaults to Danio rerio).

=item B<--ensembl_dbhost HOST>

Ensembl MySQL database host.

=item B<--ensembl_dbport PORT>

Ensembl MySQL database port.

=item B<--ensembl_dbuser USERNAME>

Ensembl MySQL database username.

=item B<--ensembl_dbpass PASSWORD>

Ensembl MySQL database password.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
