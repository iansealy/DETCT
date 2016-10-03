#!/usr/bin/env perl

# PODNAME: filter_output.pl
# ABSTRACT: Filter DETCT output using 3' ends

## Author         : is1
## Maintainer     : is1
## Created        : 2016-09-21
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
use Path::Tiny;
use File::Spec;
use Sort::Naturally;
use Set::IntervalTree;
use DETCT::GeneFinder;
use DETCT::Analysis::DiffExpr;
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-74/ensembl/modules \
        filter_output.pl \
            --dir filter \
            --analysis_yaml analysis.yaml \
            --all_file all.tsv \
            --ends_file ends.tsv

=cut

# Constants
Readonly our $THREE_PRIME_END_POS_FIELD           => 0;
Readonly our $THREE_PRIME_END_READ_COUNT_FIELD    => 1;
Readonly our $IS_POLYA_FIELD                      => 2;
Readonly our $UPSTREAM_14_BP_FIELD                => 3;
Readonly our $DOWNSTREAM_14_BP_FIELD              => 4;
Readonly our $DISTANCE_HEXAMER_UPSTREAM_FIELD     => 5;
Readonly our $HEXAMER_FIELD                       => 6;
Readonly our $TRANSPOSON_DISTANCE_FIELD           => 7;
Readonly our $TRANSPOSON_POSITION_FIELD           => 8;
Readonly our $CONTINUOUS_RNASEQ_TRANSCRIPTS_FIELD => 9;

Readonly our $DISTANCE_TO_POLY_A => 10;
Readonly our $DISTANCE_THRESHOLD => 100;

Readonly our %IS_PRIMARY_HEXAMER   => map { $_ => 1 } qw( AATAAA ATTAAA );
Readonly our %IS_SECONDARY_HEXAMER => map { $_ => 1 }
  qw( AGTAAA TATAAA CATAAA GATAAA AATATA AATACA AATAGA ACTAAA AAGAAA AATGAA );

# Default options
my $analysis_dir = q{.};
my $analysis_yaml = File::Spec->catfile( $analysis_dir, 'analysis.yaml' );
my $all_file;
my $ends_file;
my $slice_regexp;
my @biotype;
my $keep_regions_without_ends;
my $keep_ends_in_cds;
my $keep_polya_ends;
my $keep_ends_without_hexamer;
my $keep_ends_without_rnaseq;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Create analysis
my $analysis = DETCT::Analysis::DiffExpr->new_from_yaml($analysis_yaml);

my $gene_finder = DETCT::GeneFinder->new(
    {
        slice_adaptor     => $analysis->slice_adaptor,
        skip_transcripts  => $analysis->get_all_skip_transcripts,
        ensembl_db_types  => $analysis->get_all_ensembl_db_types,
        required_biotypes => \@biotype,
    }
);

my $cds_cache = get_cds_cache( $analysis->slice_adaptor, $slice_regexp );

my $ends_for = get_ends( $ends_file, $slice_regexp );

my $regions = get_regions( $all_file, $analysis, $slice_regexp );

if (@biotype) {
    $regions = remove_ends_by_biotype( $regions, $gene_finder );
}

if ( !$keep_ends_in_cds ) {
    $regions = remove_ends_in_cds( $regions, $cds_cache );
}

if ( !$keep_polya_ends ) {
    $regions = remove_polya( $regions, $ends_for );
}

if ( !$keep_ends_without_hexamer ) {
    $regions = remove_ends_without_hexamer( $regions, $ends_for, $gene_finder );
}

if ( !$keep_ends_without_rnaseq ) {
    $regions = remove_ends_without_rnaseq( $regions, $ends_for, $gene_finder );
}

if ( !$keep_regions_without_ends ) {
    $regions = remove_regions_without_ends($regions);
}

# Reannotate (perhaps only with required biotypes)
$regions = $gene_finder->add_gene_annotation($regions);

DETCT::Misc::Output::dump_as_table(
    {
        analysis => $analysis,
        dir      => $analysis_dir,
        regions  => $regions,
    }
);

sub get_cds_cache {
    my ( $sa, $slice_regexp ) = @_;    ## no critic (ProhibitReusedNames)

    my %cds_cache;

    # Get all transcripts
    my $slices = $sa->fetch_all('toplevel');
    if ($slice_regexp) {
        @{$slices} =
          grep { $_->seq_region_name =~ m/\A $slice_regexp \z/xms } @{$slices};
    }
    foreach
      my $slice ( sort { ncmp( $a->seq_region_name, $b->seq_region_name ) }
        @{$slices} )
    {
        my $chr = $slice->seq_region_name;
        if ( !exists $cds_cache{$chr} ) {
            $cds_cache{$chr}{-1} = Set::IntervalTree->new;
            $cds_cache{$chr}{1}  = Set::IntervalTree->new;
        }

        my $transcripts = $slice->get_all_Transcripts();
        foreach my $transcript ( @{$transcripts} ) {
            if ( defined $transcript->translation ) {
                my $strand = $transcript->seq_region_strand;
                my $exons  = $transcript->get_all_Exons();
                foreach my $exon ( @{$exons} ) {
                    next if !defined $exon->coding_region_start($transcript);
                    $cds_cache{$chr}{$strand}->insert(
                        {
                            start => $exon->coding_region_start($transcript),
                            end   => $exon->coding_region_end($transcript),
                        },
                        $exon->coding_region_start($transcript),
                        $exon->coding_region_end($transcript) + 1
                    );
                }
            }
        }
    }

    return \%cds_cache;
}

sub get_ends {
    ## no critic (ProhibitReusedNames)
    my ( $file, $slice_regexp ) = @_;
    ## use critic

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %ends_for;

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );
        next if $slice_regexp && $fields[0] !~ m/\A $slice_regexp \z/xms;

        # Get region ID from chr, start or end (depending on strand) and strand
        ## no critic (ProhibitMagicNumbers)
        my $strand       = $fields[3];
        my $start_or_end = $strand > 0 ? 1 : 2;
        my @region       = @fields[ 0 .. 3 ];
        my $ends         = $fields[5];
        my $region       = join q{:}, @region[ 0, $start_or_end, 3 ];
        ## use critic
        $ends_for{$region} = $ends;
    }
    close $fh;

    return \%ends_for;
}

sub get_regions {
    ## no critic (ProhibitReusedNames)
    my ( $file, $analysis, $slice_regexp ) = @_;
    ## use critic

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    ## no critic (ProhibitReusedNames)
    my $regions = DETCT::Misc::Output::parse_table(
        ## use critic
        {
            analysis     => $analysis,
            table_file   => $file,
            table_format => $extension,
        }
    );

    return $regions if !$slice_regexp;

    my @new_regions;
    foreach my $region ( @{$regions} ) {
        next if $region->[0] !~ m/\A $slice_regexp \z/xms;
        push @new_regions, $region;
    }

    return \@new_regions;
}

sub remove_ends_by_biotype {
    my ( $regions, $gene_finder ) = @_;    ## no critic (ProhibitReusedNames)

    foreach my $region ( @{$regions} ) {
        my $seq_name = $region->[0];
        ## no critic (ProhibitMagicNumbers)
        my $three_prime_end_pos = $region->[6];
        my $strand              = $region->[7];
        ## use critic
        next if !defined $three_prime_end_pos;
        my @three_prime_end_pos =
          ref $three_prime_end_pos eq 'ARRAY'
          ? @{$three_prime_end_pos}
          : ($three_prime_end_pos);
        foreach my $pos (@three_prime_end_pos) {
            my ( undef, $distance ) =
              $gene_finder->get_nearest_transcripts( $seq_name, $pos, $strand );
            if ( !defined $distance ) {
                $region = remove_end_from_region( $region, $pos );
            }
        }
    }

    return $regions;
}

sub remove_ends_in_cds {
    my ( $regions, $cds_cache ) = @_;    ## no critic (ProhibitReusedNames)

    foreach my $region ( @{$regions} ) {
        my $chr = $region->[0];
        ## no critic (ProhibitMagicNumbers)
        my $strand              = $region->[7];
        my $three_prime_end_pos = $region->[6];
        ## use critic
        next if !defined $three_prime_end_pos;
        my @three_prime_end_pos =
          ref $three_prime_end_pos eq 'ARRAY'
          ? @{$three_prime_end_pos}
          : ($three_prime_end_pos);
        foreach my $pos (@three_prime_end_pos) {
            next if !exists $cds_cache->{$chr};    # e.g. spike sequences
            my $cds_intervals =
              $cds_cache->{$chr}{$strand}->fetch( $pos, $pos + 1 );
            if ( @{$cds_intervals} ) {
                $region = remove_end_from_region( $region, $pos );
            }
        }
    }

    return $regions;
}

sub remove_polya {
    my ( $regions, $ends_for ) = @_;    ## no critic (ProhibitReusedNames)

    foreach my $region ( @{$regions} ) {
        my @all_ends = parse_ends( $region, $ends_for );
        ## no critic (ProhibitMagicNumbers)
        my $three_prime_end_pos = $region->[6];
        ## use critic
        next if !defined $three_prime_end_pos;
        my @three_prime_end_pos =
          ref $three_prime_end_pos eq 'ARRAY'
          ? @{$three_prime_end_pos}
          : ($three_prime_end_pos);
        foreach my $pos (@three_prime_end_pos) {
            foreach my $end (@all_ends) {
                if ( $end->[$IS_POLYA_FIELD] eq 'y'
                    && abs( $end->[$THREE_PRIME_END_POS_FIELD] - $pos ) <
                    $DISTANCE_TO_POLY_A )
                {
                    $region = remove_end_from_region( $region, $pos );
                    last;
                }
            }
        }
    }

    return $regions;
}

sub remove_ends_without_hexamer {
    ## no critic (ProhibitReusedNames)
    my ( $regions, $ends_for, $gene_finder ) = @_;
    ## use critic

    foreach my $region ( @{$regions} ) {
        my @all_ends = parse_ends( $region, $ends_for );
        my $seq_name = $region->[0];
        ## no critic (ProhibitMagicNumbers)
        my $three_prime_end_pos = $region->[6];
        my $strand              = $region->[7];
        ## use critic
        next if !defined $three_prime_end_pos;
        my @three_prime_end_pos =
          ref $three_prime_end_pos eq 'ARRAY'
          ? @{$three_prime_end_pos}
          : ($three_prime_end_pos);
        foreach my $pos (@three_prime_end_pos) {
            my ($end) =
              grep { $_->[$THREE_PRIME_END_POS_FIELD] == $pos } @all_ends;

            # Keep end if primary hexamer
            next
              if defined $end->[$HEXAMER_FIELD]
              && exists $IS_PRIMARY_HEXAMER{ $end->[$HEXAMER_FIELD] };

            # Or keep end if secondary hexamer and near annotated 3' end
            my $remove = 1;
            if ( defined $end->[$HEXAMER_FIELD]
                && exists $IS_SECONDARY_HEXAMER{ $end->[$HEXAMER_FIELD] } )
            {
                my ( undef, $distance ) =
                  $gene_finder->get_nearest_transcripts( $seq_name, $pos,
                    $strand );
                if ( defined $distance && abs $distance <= $DISTANCE_THRESHOLD )
                {
                    $remove = 0;
                }
            }
            if ($remove) {
                $region = remove_end_from_region( $region, $pos );
            }
        }
    }

    return $regions;
}

sub remove_ends_without_rnaseq {
    ## no critic (ProhibitReusedNames)
    my ( $regions, $ends_for, $gene_finder ) = @_;
    ## use critic

    foreach my $region ( @{$regions} ) {
        my @all_ends = parse_ends( $region, $ends_for );
        my $seq_name = $region->[0];
        ## no critic (ProhibitMagicNumbers)
        my $three_prime_end_pos = $region->[6];
        my $strand              = $region->[7];
        ## use critic
        next if !defined $three_prime_end_pos;
        my @three_prime_end_pos =
          ref $three_prime_end_pos eq 'ARRAY'
          ? @{$three_prime_end_pos}
          : ($three_prime_end_pos);
        foreach my $pos (@three_prime_end_pos) {
            my ($end) =
              grep { $_->[$THREE_PRIME_END_POS_FIELD] == $pos } @all_ends;

            # Keep end if continuous RNA-Seq
            next if $end->[$CONTINUOUS_RNASEQ_TRANSCRIPTS_FIELD];

            # Or keep end if near annotated 3' end
            my ( undef, $distance ) =
              $gene_finder->get_nearest_transcripts( $seq_name, $pos, $strand );
            if ( !defined $distance || abs $distance > $DISTANCE_THRESHOLD ) {
                $region = remove_end_from_region( $region, $pos );
            }
        }
    }

    return $regions;
}

sub remove_regions_without_ends {
    my ($regions) = @_;    ## no critic (ProhibitReusedNames)

    my @new_regions;
    foreach my $region ( @{$regions} ) {
        ## no critic (ProhibitMagicNumbers)
        my $three_prime_end_pos = $region->[6];
        ## use critic
        next if !defined $three_prime_end_pos;
        push @new_regions, $region;
    }

    return \@new_regions;
}

sub remove_end_from_region {
    my ( $region, $pos_to_remove ) = @_;

    ## no critic (ProhibitMagicNumbers)
    my $three_prime_end_pos = $region->[6];
    ## use critic

    if ( ref $three_prime_end_pos ne 'ARRAY'
        && $three_prime_end_pos != $pos_to_remove )
    {
        confess sprintf q{Trying to remove unknown 3' end (%d)}, $pos_to_remove;
    }

    if ( ref $three_prime_end_pos ne 'ARRAY' ) {
        ## no critic (ProhibitMagicNumbers)
        $region->[6]  = undef;    # 3' end position
        $region->[8]  = undef;    # 3' end read count
        $region->[15] = {};       # Gene annotation
        ## use critic
    }
    else {
        my $end_count = scalar @{$three_prime_end_pos};
        my $index     = 0;
        while ( $three_prime_end_pos->[$index] != $pos_to_remove ) {
            $index++;
            confess sprintf q{Trying to remove unknown 3' end (%d)},
              $pos_to_remove
              if !defined $three_prime_end_pos->[$index];
        }
        ## no critic (ProhibitMagicNumbers)
        splice @{ $region->[6] }, $index, 1;
        splice @{ $region->[8] }, $index, 1;
        if ( scalar @{ $region->[6] } == 1 ) {
            $region->[6] = $region->[6]->[0];
        }
        if ( scalar @{ $region->[8] } == 1 ) {
            $region->[8] = $region->[8]->[0];
        }

        # Attempt to remove 3' end distance, but really should reannotate
        my $genebuild_version = ( keys %{ $region->[15] } )[0];
        if ($genebuild_version) {
            my $genes = $region->[15]->{$genebuild_version};
            my @new_genes;
            foreach my $gene ( @{$genes} ) {
                next if ref $gene->[4] ne 'ARRAY';
                if ( scalar @{ $gene->[4] } == $end_count ) {
                    splice @{ $gene->[4] }, $index, 1;
                    if ( scalar @{ $gene->[4] } == 1 ) {
                        $gene->[4] = $gene->[4]->[0];
                    }
                }
                push @new_genes, $gene;
            }
            if (@new_genes) {
                $region->[15]->{$genebuild_version} = \@new_genes;
            }
            else {
                $region->[15] = {};
            }
        }
        ## use critic
    }

    return $region;
}

sub parse_ends {
    my ( $region, $ends_for ) = @_;    ## no critic (ProhibitReusedNames)

    # Get region ID from chr, start or end (depending on strand) and strand
    ## no critic (ProhibitMagicNumbers)
    my $strand = $region->[7];
    my $start_or_end = $strand > 0 ? 1 : 2;
    my $ends = $ends_for->{ join q{:}, @{$region}[ 0, $start_or_end, 7 ] };
    ## use critic
    my @ends;
    if ($ends) {
        foreach my $end ( split /,/xms, $ends ) {
            my @end_data = split /[:\/]/xms, $end;
            splice @end_data, 1, 1;    # Remove redundant strand
            if ( scalar @end_data == 9 ) {   ## no critic (ProhibitMagicNumbers)
                push @end_data, q{};         # Add missing last field
            }
            push @ends, \@end_data;
        }
        @ends = sort { $a->[0] <=> $b->[0] } @ends;
        if ( $strand < 0 ) {
            @ends = reverse @ends;
        }
    }

    return @ends;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'dir=s'                     => \$analysis_dir,
        'analysis_yaml=s'           => \$analysis_yaml,
        'all_file=s'                => \$all_file,
        'ends_file=s'               => \$ends_file,
        'slice_regexp=s'            => \$slice_regexp,
        'biotype=s@{,}'             => \@biotype,
        'keep_regions_without_ends' => \$keep_regions_without_ends,
        'keep_ends_in_cds'          => \$keep_ends_in_cds,
        'keep_polya_ends'           => \$keep_polya_ends,
        'keep_ends_without_hexamer' => \$keep_ends_without_hexamer,
        'keep_ends_without_rnaseq'  => \$keep_ends_without_rnaseq,
        'help'                      => \$help,
        'man'                       => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$all_file ) {
        pod2usage("--all_file must be specified\n");
    }
    if ( !$ends_file ) {
        pod2usage("--ends_file must be specified\n");
    }

    return;
}

=head1 USAGE

    filter_output.pl
        [--dir directory]
        [--analysis_yaml file]
        [--all_file file]
        [--ends_file file]
        [--slice_regexp regexp]
        [--biotype biotype...]
        [--keep_regions_without_ends]
        [--keep_ends_in_cds]
        [--keep_polya_ends]
        [--keep_ends_without_hexamer]
        [--keep_ends_without_rnaseq]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--dir DIRECTORY>

Output directory.

=item B<--analysis_yaml FILE>

YAML analysis configuration file.

=item B<--all_file FILE>

DETCT output file (containing all regions).

=item B<--ends_file FILE>

DETCT 3' ends file.

=item B<--slice_regexp REGEXP>

Regular expression for matching slice names.

=item B<--biotype BIOTYPES>

Required biotype(s).

=item B<--keep_regions_without_ends>

Don't filter out regions lacking a 3' end.

=item B<--keep_ends_in_cds>

Don't filter 3' ends in a transcript's CDS.

=item B<--keep_polya_ends>

Don't filter out polyA 3' ends labelled.

=item B<--keep_ends_without_hexamer>

Don't filter out 3' ends lacking a primary hexamer.

=item B<--keep_ends_without_rnaseq>

Don't filter out 3' ends lacking continuous RNA-Seq transcripts.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
