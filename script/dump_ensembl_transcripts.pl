#!/usr/bin/env perl

# PODNAME: dump_ensembl_transcripts.pl
# ABSTRACT: Dump info for Ensembl transcripts

## Author         : is1
## Maintainer     : is1
## Created        : 2015-07-10
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
use Bio::EnsEMBL::Registry;
use Sort::Naturally;
use Set::IntervalTree;
use DETCT::Misc qw( printf_or_die );
use DETCT::Misc::BAM;

=head1 DESCRIPTION



=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-74/ensembl/modules \
        dump_ensembl_transcripts.pl \
        > data/e74-transcripts.txt

=cut

# Default options
my $species        = 'Danio rerio';
my $ensembl_dbhost = 'ensembldb.ensembl.org';
my $ensembl_dbport;
my $ensembl_dbuser = 'anonymous';
my $ensembl_dbpass;
my $slice_regexp;
my $wgs_prefix = 'CABZ';
my $rnaseq_bedgraph;
## no critic (ProhibitMagicNumbers)
my $rnaseq_threshold = 3;
my $repeat_padding   = 20;
## use critic
my ( $debug, $help, $man );

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
warn 'Genebuild version: ', $genebuild_version, "\n" if $debug;

# Get Ensembl adaptors
my $sa = Bio::EnsEMBL::Registry->get_adaptor( $species, 'core', 'Slice' );
my $asia =
  Bio::EnsEMBL::Registry->get_adaptor( $species, 'core', 'ArchiveStableId' );

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Cache RNA-Seq locations
my $rnaseq_cache = get_rnaseq_locations();
warn "Loaded RNA-Seq cache\n" if $debug;

# Print header
my @HEADERS = (
    'Chr',
    'Start',
    'End',
    'Strand',
    'Transcript ID',
    'Gene ID',
    'Transcript Biotype',
    'Gene Biotype',
    q{3' UTR Length},
    'APPRIS',
    'Transcript ID Initial Version',
    'Gene ID Initial Version',
    'Overlapping Genes',
    'Distance Hexamer Upstream (up to 40 bp)',
    'Hexamer',
    'polyA?',
    '14 bp Downstream',
    'Sequences',
    'Distance to Gap',
    'Distance to WGS',
    'Distance to Transcript CDS End',
    'Continuous RNA-Seq End',
    'Distance to Continuous RNA-Seq End',
    'Genes Extended Over',
    'Distance Hexamer Upstream of RNA-Seq End (up to 40 bp)',
    'RNA-Seq End Hexamer',
    'RNA-Seq End polyA?',
    'RNA-Seq End 14 bp Downstream',
);
printf_or_die( "#%s\n", join "\t", @HEADERS );

# Get all transcripts
my $slices = $sa->fetch_all('toplevel');
warn scalar @{$slices}, " slices\n" if $debug;
if ($slice_regexp) {
    @{$slices} =
      grep { $_->seq_region_name =~ m/\A $slice_regexp \z/xms } @{$slices};
    warn scalar @{$slices}, " slices after filtering\n" if $debug;
}
foreach my $slice ( sort { ncmp( $a->seq_region_name, $b->seq_region_name ) }
    @{$slices} )
{
    warn 'Slice: ', $slice->name, "\n" if $debug;

    my %gene_id_version_cache;

    # Cache gap and WGS locations
    my ( $gap_cache, $wgs_cache ) = get_gap_and_wgs_locations($slice);
    warn "Loaded gap and WGS caches\n" if $debug;

    # Cache simple repeat locations
    my $repeat_cache = get_simple_repeat_locations($slice);
    warn "Loaded simple repeat cache\n" if $debug;

    my $transcripts = $slice->get_all_Transcripts();
    warn scalar @{$transcripts}, " transcripts\n" if $debug;
    my $n = 0;
    foreach my $transcript (
        sort {
                 $a->seq_region_start <=> $b->seq_region_start
              || $a->seq_region_end <=> $b->seq_region_end
        } @{$transcripts}
      )
    {
        dump_transcript( $transcript, $slice, \%gene_id_version_cache,
            $gap_cache, $wgs_cache, $rnaseq_cache->{ $slice->seq_region_name },
            $repeat_cache );
    }
}

sub dump_transcript {    ## no critic (ProhibitManyArgs)
    ## no critic (ProhibitReusedNames)
    my ( $transcript, $slice, $gene_id_version_cache, $gap_cache,
        $wgs_cache, $rnaseq_cache, $repeat_cache )
      = @_;
    ## use critic

    warn 'Transcript: ', $transcript->stable_id, "\n" if $debug;

    my @fields;

    push @fields, $transcript->seq_region_name;
    push @fields, $transcript->seq_region_start;
    push @fields, $transcript->seq_region_end;
    push @fields, $transcript->seq_region_strand;
    push @fields, $transcript->stable_id;
    push @fields, $transcript->get_Gene->stable_id;
    push @fields, $transcript->biotype;
    push @fields, $transcript->get_Gene->biotype;
    push @fields, get_three_prime_utr_length($transcript);
    push @fields, get_appris($transcript);
    push @fields, get_stable_id_initial_version( $transcript->stable_id );

    if ( !exists $gene_id_version_cache->{ $transcript->get_Gene->stable_id } )
    {
        $gene_id_version_cache->{ $transcript->get_Gene->stable_id } =
          get_stable_id_initial_version( $transcript->get_Gene->stable_id );
    }
    push @fields, $gene_id_version_cache->{ $transcript->get_Gene->stable_id };
    push @fields, get_overlapping_genes($transcript);
    push @fields,
      get_hexamer(
        $transcript->seq_region_name, $transcript->seq_region_start,
        $transcript->seq_region_end,  $transcript->seq_region_strand
      );
    push @fields,
      get_polya(
        $transcript->seq_region_name, $transcript->seq_region_start,
        $transcript->seq_region_end,  $transcript->seq_region_strand
      );
    push @fields, get_sequence_type($transcript);
    push @fields,
      get_downstream_gap_and_wgs( $transcript, $slice->length, $gap_cache,
        $wgs_cache );
    my ( $cds_upstream, $extension_pos, $extension_dist, $genes_extended_over )
      = get_rnaseq_extension( $transcript, $gap_cache, $rnaseq_cache,
        $repeat_cache );
    push @fields, $cds_upstream, $extension_pos, $extension_dist,
      $genes_extended_over;

    if ($extension_dist) {
        push @fields,
          get_hexamer( $transcript->seq_region_name,
            $extension_pos, $extension_pos, $transcript->seq_region_strand );
        push @fields,
          get_polya( $transcript->seq_region_name,
            $extension_pos, $extension_pos, $transcript->seq_region_strand );
    }
    else {
        push @fields, undef, undef, undef, undef;
    }

    @fields = map { ( !defined $_ || $_ eq q{} ) ? q{-} : $_ } @fields;
    printf_or_die( "%s\n", join "\t", @fields );

    return;
}

sub get_three_prime_utr_length {
    my ($transcript) = @_;

    if ( defined $transcript->three_prime_utr ) {
        return $transcript->three_prime_utr->length;
    }

    return undef;    ## no critic (ProhibitExplicitReturnUndef)
}

sub get_appris {
    my ($transcript) = @_;

    foreach my $attr ( @{ $transcript->get_all_Attributes || [] } ) {
        return $attr->code if ( $attr->code =~ /\A appris/xms );
    }

    return undef;    ## no critic (ProhibitExplicitReturnUndef)
}

sub get_stable_id_initial_version {
    my ($stable_id) = @_;

    my $history = $asia->fetch_history_tree_by_stable_id($stable_id);
    my $initial_version;
    foreach my $archive ( @{ $history->get_all_ArchiveStableIds } ) {
        next if $archive->stable_id() ne $stable_id;
        if ( !defined $initial_version || $archive->release < $initial_version )
        {
            $initial_version = $archive->release;
        }
    }

    return $initial_version;
}

sub get_overlapping_genes {
    my ($transcript) = @_;

    my $overlapping_genes = $transcript->feature_Slice->get_all_Genes();
    my @genes;
    foreach my $overlapping_gene ( @{$overlapping_genes} ) {
        next
          if $overlapping_gene->stable_id eq $transcript->get_Gene->stable_id;
        next
          if $overlapping_gene->seq_region_strand !=
          $transcript->seq_region_strand;
        push @genes, $overlapping_gene->stable_id;
    }

    return join q{,}, @genes;
}

sub get_hexamer {
    my ( $chr, $transcript_start, $transcript_end, $transcript_strand ) = @_;

    my $SEQ_TO_CHECK_FOR_HEXAMER = 40;

    # From http://www.ncbi.nlm.nih.gov/pmc/articles/PMC310884/
    my @PRIMARY_HEXAMERS = qw( AATAAA ATTAAA );
    my @SECONDARY_HEXAMERS =
      qw( AGTAAA TATAAA CATAAA GATAAA AATATA AATACA AATAGA ACTAAA AAGAAA AATGAA );
    my @primary_hexamers   = map { scalar reverse } @PRIMARY_HEXAMERS;
    my @secondary_hexamers = map { scalar reverse } @SECONDARY_HEXAMERS;

    my $start;
    my $end;
    if ( $transcript_strand > 0 ) {
        $start = $transcript_end - $SEQ_TO_CHECK_FOR_HEXAMER + 1;
        $end   = $transcript_end;
    }
    else {
        $start = $transcript_start;
        $end   = $transcript_start + $SEQ_TO_CHECK_FOR_HEXAMER - 1;
    }
    my $upstream_seq =
      reverse $sa->fetch_by_region( 'toplevel', $chr, $start, $end,
        $transcript_strand )->seq;

    foreach my $hexamer (@primary_hexamers) {
        my $pos = index $upstream_seq, $hexamer;
        if ( $pos != -1 ) {
            return $pos + 5, scalar reverse $hexamer;
        }
    }

    my $nearest_pos;
    my $nearest_hexamer;
    foreach my $hexamer (@secondary_hexamers) {
        my $pos = index $upstream_seq, $hexamer;
        if ( $pos != -1
            && ( !defined $nearest_pos || $pos + 5 < $nearest_pos ) )
        {
            $nearest_pos     = $pos + 5;
            $nearest_hexamer = scalar reverse $hexamer;
        }
    }

    return $nearest_pos, $nearest_hexamer;
}

sub get_polya {
    my ( $chr, $transcript_start, $transcript_end, $transcript_strand ) = @_;

    my $start;
    my $end;
    if ( $transcript_strand > 0 ) {
        $start = $transcript_end + 1;
        $end   = $transcript_end + 14;
    }
    else {
        $start = $transcript_start - 14;
        $end   = $transcript_start - 1;
    }
    my $seq =
      $sa->fetch_by_region( 'toplevel', $chr, $start, $end, $transcript_strand )
      ->seq;

    return DETCT::Misc::BAM::is_polya( substr $seq, 0, 10 ) ? q{y} : q{n}, $seq;
}

sub get_sequence_type {
    my ($transcript) = @_;

    my @sequences;

    my $projection = $transcript->feature_Slice->project('seqlevel');
    foreach my $segment ( @{$projection} ) {
        push @sequences, $segment->to_Slice->seq_region_name;
    }

    return join q{,}, @sequences;
}

sub get_downstream_gap_and_wgs {
    my ( $transcript, $slice_length, $gap_cache, $wgs_cache ) = @_;

    my $gap_distance;
    my $wgs_distance;
    if ( $transcript->seq_region_strand > 0 ) {
        $gap_distance =
          get_next_forward_gap_or_wgs( $transcript, $gap_cache, $slice_length );
        $wgs_distance =
          get_next_forward_gap_or_wgs( $transcript, $wgs_cache, $slice_length );
    }
    else {
        $gap_distance = get_next_reverse_gap_or_wgs( $transcript, $gap_cache );
        $wgs_distance = get_next_reverse_gap_or_wgs( $transcript, $wgs_cache );
    }

    return $gap_distance, $wgs_distance;
}

sub get_next_forward_gap_or_wgs {
    my ( $transcript, $cache, $slice_length ) = @_;

    my $distance;
    my $intervals =
      $cache->fetch( $transcript->seq_region_end, $slice_length + 1 );
    @{$intervals} = sort { $a->{start} <=> $b->{start} } @{$intervals};
    if ( @{$intervals} ) {
        $distance = ${$intervals}[0]->{start} - $transcript->seq_region_end;
        $distance = $distance < 0 ? 0 : $distance;
    }

    return $distance;
}

sub get_next_reverse_gap_or_wgs {
    my ( $transcript, $cache ) = @_;

    my $distance;
    my $intervals = $cache->fetch( 0, $transcript->seq_region_start );
    @{$intervals} = sort { $a->{end} <=> $b->{end} } @{$intervals};
    if ( @{$intervals} ) {
        $distance = $transcript->seq_region_start - ${$intervals}[-1]->{end};
        $distance = $distance < 0 ? 0 : $distance;
    }

    return $distance;
}

sub get_rnaseq_extension {
    ## no critic (ProhibitReusedNames)
    my ( $transcript, $gap_cache, $rnaseq_cache, $repeat_cache ) = @_;
    ## use critic

    my $cds_upstream;
    my $extension_pos;
    if ( $transcript->seq_region_strand > 0 ) {
        $extension_pos = $transcript->seq_region_end;
        if ( defined $transcript->translation ) {
            $cds_upstream =
              $transcript->seq_region_end - $transcript->coding_region_end;
            $extension_pos = $transcript->coding_region_end;
        }
    }
    else {
        $extension_pos = $transcript->seq_region_start;
        if ( defined $transcript->translation ) {
            $cds_upstream =
              $transcript->coding_region_start - $transcript->seq_region_start;
            $extension_pos = $transcript->coding_region_start;
        }
    }

    my $still_extending = 1;
    while ($still_extending) {
        $extension_pos += $transcript->seq_region_strand;
        my $gap_intervals =
          $gap_cache->fetch( $extension_pos, $extension_pos + 1 );
        if ( @{$gap_intervals} ) {
            $still_extending = 0;
        }
        my $rnaseq_intervals =
          $rnaseq_cache->fetch( $extension_pos, $extension_pos + 1 );
        my $repeat_intervals =
          $repeat_cache->fetch( $extension_pos, $extension_pos + 1 );
        if ( !@{$rnaseq_intervals} && !@{$repeat_intervals} ) {
            $still_extending = 0;
        }
    }
    $extension_pos -= $transcript->seq_region_strand;

    my $extension_diff;
    if ( $transcript->seq_region_strand > 0 ) {
        $extension_diff = $extension_pos - $transcript->seq_region_end;
    }
    else {
        $extension_diff = $transcript->seq_region_start - $extension_pos;
    }

    my $genes_extended_over_count = 0;
    if ( $extension_diff > 0 ) {
        my ( $slice_start, $slice_end );
        if ( $transcript->seq_region_strand > 0 ) {
            $slice_start = $transcript->seq_region_end + 1;
            $slice_end   = $extension_pos;
        }
        else {
            $slice_start = $extension_pos;
            $slice_end   = $transcript->seq_region_start - 1;
        }
        my $extension_slice =
          $sa->fetch_by_region( 'toplevel', $transcript->seq_region_name,
            $slice_start, $slice_end, $transcript->seq_region_strand );
        my $genes_extended_over = $extension_slice->get_all_Genes();
        foreach my $gene_extended_over ( @{$genes_extended_over} ) {
            next
              if $gene_extended_over->stable_id eq
              $transcript->get_Gene->stable_id;
            next
              if $gene_extended_over->seq_region_strand !=
              $transcript->seq_region_strand;
            $genes_extended_over_count++;
        }
    }

    # Don't extend at all if extending over other genes
    if ($genes_extended_over_count) {
        $extension_pos  = undef;
        $extension_diff = undef;
    }

    return $cds_upstream, $extension_pos, $extension_diff,
      $genes_extended_over_count;
}

sub get_gap_and_wgs_locations {
    my ($slice) = @_;

    my $gap_tree   = Set::IntervalTree->new;
    my $wgs_tree   = Set::IntervalTree->new;
    my $projection = $slice->project('seqlevel');
    my $prev_end   = 0;
    foreach my $segment ( @{$projection} ) {
        my $start = $segment->from_start;
        my $end   = $segment->from_end;
        if ( $prev_end + 1 < $start ) {

            # Got gap
            $gap_tree->insert(
                {
                    start => $prev_end + 1,
                    end   => $start - 1,
                },
                $prev_end + 1,
                $start
            );
        }
        if ( $segment->to_Slice->seq_region_name =~ m/\A $wgs_prefix/xms ) {

            # Got WGS
            $wgs_tree->insert(
                {
                    start => $start,
                    end   => $end,
                },
                $start,
                $end + 1
            );
        }
        $prev_end = $end;
    }

    return $gap_tree, $wgs_tree;
}

sub get_rnaseq_locations {
    my ($slice_regexp) = @_;    ## no critic (ProhibitReusedNames)

    my %rnaseq_cache;

    my $chr;
    my $cur_chr;
    my $cur_start;
    my $cur_end;
    open my $fh, '<', $rnaseq_bedgraph;    ## no critic (RequireBriefOpen)
    while ( my $line = <$fh> ) {
        chomp $line;
        ## no critic (ProhibitReusedNames)
        my ( $chr, $start, $end, $reads ) = split /\t/xms, $line;
        ## use critic
        next if $slice_regexp && $chr !~ m/\A $slice_regexp \z/xms;
        $start++;
        if ( !exists $rnaseq_cache{$chr} ) {
            $rnaseq_cache{$chr} = Set::IntervalTree->new;
        }
        next if $reads < $rnaseq_threshold;
        if ( !defined $cur_chr ) {
            $cur_chr   = $chr;
            $cur_start = $start;
            $cur_end   = $end;
        }
        if ( $cur_chr ne $chr || $start - $cur_end > 1 ) {

            # Starting new region so store current one
            $rnaseq_cache{$cur_chr}->insert(
                {
                    start => $cur_start,
                    end   => $cur_end,
                },
                $cur_start,
                $cur_end + 1
            );
            $cur_chr   = $chr;
            $cur_start = $start;
            $cur_end   = $end;
        }
        else {
            # Extend current region
            $cur_end = $end;
        }
    }

    # Store final region
    if ( $cur_start < $cur_end ) {
        $rnaseq_cache{$cur_chr}->insert(
            {
                start => $cur_start,
                end   => $cur_end,
            },
            $cur_start,
            $cur_end + 1
        );
    }

    return \%rnaseq_cache;
}

sub get_simple_repeat_locations {
    my ($slice) = @_;

    my $repeat_tree = Set::IntervalTree->new;

    my $dust_features = $slice->get_all_RepeatFeatures('dust');
    my $trf_features  = $slice->get_all_RepeatFeatures('trf');
    foreach my $repeat ( @{$dust_features}, @{$trf_features} ) {
        $repeat_tree->insert(
            {
                start => $repeat->seq_region_start - $repeat_padding,
                end   => $repeat->seq_region_end + $repeat_padding,
            },
            $repeat->seq_region_start - $repeat_padding,
            $repeat->seq_region_end + $repeat_padding + 1
        );
    }

    return $repeat_tree;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'species=s'          => \$species,
        'ensembl_dbhost=s'   => \$ensembl_dbhost,
        'ensembl_dbport=i'   => \$ensembl_dbport,
        'ensembl_dbuser=s'   => \$ensembl_dbuser,
        'ensembl_dbpass=s'   => \$ensembl_dbpass,
        'slice_regexp=s'     => \$slice_regexp,
        'wgs_prefix=s'       => \$wgs_prefix,
        'rnaseq_bedgraph=s'  => \$rnaseq_bedgraph,
        'rnaseq_threshold=i' => \$rnaseq_threshold,
        'repeat_padding=i'   => \$repeat_padding,
        'debug'              => \$debug,
        'help'               => \$help,
        'man'                => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

=head1 USAGE

    dump_ensembl_transcripts.pl
        [--species species]
        [--ensembl_dbhost host]
        [--ensembl_dbport port]
        [--ensembl_dbuser username]
        [--ensembl_dbpass password]
        [--slice_regexp regexp]
        [--wgs_prefix prefix]
        [--rnaseq_bedgraph file]
        [--rnaseq_threshold int]
        [--repeat_padding int]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

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

=item B<--slice_regexp REGEXP>

Regular expression for matching slice names.

=item B<--wgs_prefix PREFIX>

Prefix identifying WGS contig names.

=item B<--rnaseq_bedgraph FILE>

BedGraph file of RNA-Seq data.

=item B<--rnaseq_threshold INT>

RNA-Seq read threshold.

=item B<--repeat_padding INT>

Padding around repeats.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
