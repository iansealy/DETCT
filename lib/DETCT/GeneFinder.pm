## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::GeneFinder;
## use critic

# ABSTRACT: Object for finding genes (and transcripts) by location

## Author         : is1
## Maintainer     : is1
## Created        : 2012-11-24
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Class::InsideOut qw( private register id );
use DETCT::Gene;
use DETCT::Transcript;

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private slice_adaptor => my %slice_adaptor;  # Bio::EnsEMBL::DBSQL::SliceAdaptor
private cache         => my %cache;          # Hashref

=method new

  Usage       : my $gene_finder = DETCT::GeneFinder->new( {
                    slice_adaptor => $slice_adaptor,
                } );
  Purpose     : Constructor for gene finder objects
  Returns     : DETCT::GeneFinder
  Parameters  : Hashref {
                    slice_adaptor => Bio::EnsEMBL::DBSQL::SliceAdaptor,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_slice_adaptor( $arg_ref->{slice_adaptor} );
    return $self;
}

=method slice_adaptor

  Usage       : my $slice_adaptor = $analysis->slice_adaptor;
  Purpose     : Getter for Ensembl slice adaptor attribute
  Returns     : Bio::EnsEMBL::DBSQL::SliceAdaptor
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub slice_adaptor {
    my ($self) = @_;
    return $slice_adaptor{ id $self};
}

=method set_slice_adaptor

  Usage       : $analysis->set_slice_adaptor($slice_adaptor);
  Purpose     : Setter for Ensembl slice adaptor attribute
  Returns     : undef
  Parameters  : Bio::EnsEMBL::DBSQL::SliceAdaptor
  Throws      : No exceptions
  Comments    : None

=cut

sub set_slice_adaptor {
    my ( $self, $arg ) = @_;
    $slice_adaptor{ id $self} = _check_slice_adaptor($arg);
    return;
}

# Usage       : $slice_adaptor = _check_slice_adaptor($slice_adaptor);
# Purpose     : Check for valid Ensembl slice adaptor
# Returns     : Bio::EnsEMBL::DBSQL::SliceAdaptor
# Parameters  : Bio::EnsEMBL::DBSQL::SliceAdaptor
# Throws      : If slice adaptor is missing or invalid (i.e. not a
#               Bio::EnsEMBL::DBSQL::SliceAdaptor object)
# Comments    : None
sub _check_slice_adaptor {
    my ($slice_adaptor) = @_;
    return $slice_adaptor
      if defined $slice_adaptor
      && $slice_adaptor->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor');
    confess 'No Ensembl slice adaptor specified' if !defined $slice_adaptor;
    confess 'Class of Ensembl slice adaptor (', ref $slice_adaptor,
      ') not Bio::EnsEMBL::DBSQL::SliceAdaptor';
}

=method get_nearest_transcripts

  Usage       : $gene_finder->get_nearest_transcripts($seq_name, $pos, $strand);
  Purpose     : Retrieve the nearest transcripts to a 3' end
  Returns     : Arrayref (of DETCT::Transcript objects)
                Int (distance)
                Int (nearest 3' end position)
  Parameters  : String (the 3' end sequence name)
                Int (the 3' end position)
                Int (the 3' end strand)
  Throws      : No exceptions
  Comments    : Distance is positive if downstream of 3' end and negative if
                upstream

=cut

sub get_nearest_transcripts {
    my ( $self, $seq_name, $pos, $strand ) = @_;

    # Ensure cache is filled
    $self->_fill_cache_from_ensembl($seq_name);

    my $nearest_distance;
    my $nearest_transcripts = [];
    my $nearest_end_pos;

    # Iterate over all 3' end transcript positions in relevant portion of cache
    my @transcript_positions = keys %{ $cache{ id $self}->{$seq_name} };

    # Favour upstream if get two transcripts same distance upstream and
    # downstream (and strand is known)
    @transcript_positions = sort { $a <=> $b } @transcript_positions;
    ## no critic (ProhibitMagicNumbers)
    if ( defined $strand && $strand == -1 ) {
        ## use critic
        @transcript_positions = reverse @transcript_positions;
    }

    foreach my $transcript_position (@transcript_positions) {
        my @transcripts =
          @{ $cache{ id $self}->{$seq_name}->{$transcript_position} };

        # Only consider transcripts matching strand (if specified)
        if ( defined $strand ) {
            @transcripts = grep { $_->strand == $strand } @transcripts;
        }
        next if !@transcripts;

        my $distance = $pos - $transcript_position;
        ## no critic (ProhibitMagicNumbers)
        if ( $transcripts[0]->strand == -1 ) {
            ### use critic
            $distance = -$distance;
        }

        # Keep transcripts if nearer than seen before
        if (  !defined $nearest_distance
            || abs $distance < abs $nearest_distance )
        {
            $nearest_transcripts = \@transcripts;
            $nearest_distance    = $distance;
            $nearest_end_pos     = $transcript_position;
        }
    }

    # Sort by stable id
    @{$nearest_transcripts} =
      sort { $a->stable_id cmp $b->stable_id } @{$nearest_transcripts};

    return $nearest_transcripts, $nearest_distance, $nearest_end_pos;
}

=method get_nearest_genes

  Usage       : $gene_finder->get_nearest_genes($seq_name, $pos, $strand);
  Purpose     : Retrieve the nearest genes to a 3' end
  Returns     : Arrayref (of DETCT::Gene objects)
                Int (distance)
                Int (nearest 3' end position)
  Parameters  : String (the 3' end sequence name)
                Int (the 3' end position)
                Int (the 3' end strand)
  Throws      : No exceptions
  Comments    : Distance is positive if downstream of 3' end and negative if
                upstream

=cut

sub get_nearest_genes {
    my ( $self, $seq_name, $pos, $strand ) = @_;

    my ( $transcripts, $distance, $nearest_end_pos ) =
      $self->get_nearest_transcripts( $seq_name, $pos, $strand );

    my %tmp_cache;    # Temporarily store genes by stable id

    # Get all genes corresponding to these transcripts
    foreach my $transcript ( @{$transcripts} ) {
        $tmp_cache{ $transcript->gene->stable_id } = $transcript->gene;
    }

    my $nearest_genes = [ values %tmp_cache ];

    # Sort by stable id
    @{$nearest_genes} =
      sort { $a->stable_id cmp $b->stable_id } @{$nearest_genes};

    return $nearest_genes, $distance, $nearest_end_pos;
}

# Usage       : $self->_fill_cache_from_ensembl( $seq_name );
# Purpose     : Fill the cache from Ensembl for a particular sequence
# Returns     : undef
# Parameters  : String (the sequence name)
# Throws      : No exceptions
# Comments    : Cache is a hashref (keyed by sequence name) of hashrefs (keyed
#               by 3' end position) of arrayrefs of transcripts

sub _fill_cache_from_ensembl {
    my ( $self, $seq_name ) = @_;

    # Skip if cache already filled
    return if exists $cache{ id $self}->{$seq_name};

    # Make sure default key exists (in case there are no genes)
    $cache{ id $self}->{$seq_name} = {};

    my $slice = $self->slice_adaptor->fetch_by_region( 'toplevel', $seq_name );

    require Bio::EnsEMBL::ApiVersion;
    my $genebuild_version = 'e' . Bio::EnsEMBL::ApiVersion::software_version();

    my $ens_genes = $slice->get_all_Genes( undef, undef, 1 ); # Plus transcripts
    foreach my $ens_gene ( @{$ens_genes} ) {
        my $gene = DETCT::Gene->new(
            {
                genebuild_version => $genebuild_version,
                stable_id         => $ens_gene->stable_id,
                name              => $ens_gene->external_name,
                description       => $ens_gene->description,
                biotype           => $ens_gene->biotype,
                seq_name          => $seq_name,
                start             => $ens_gene->seq_region_start,
                end               => $ens_gene->seq_region_end,
                strand            => $ens_gene->seq_region_strand,
            }
        );

        # Get 3' end position for each transcript
        my $ens_transcripts = $ens_gene->get_all_Transcripts();
        foreach my $ens_transcript ( @{$ens_transcripts} ) {
            my $transcript = DETCT::Transcript->new(
                {
                    stable_id   => $ens_transcript->stable_id,
                    name        => $ens_transcript->external_name,
                    description => $ens_transcript->description,
                    biotype     => $ens_transcript->biotype,
                    seq_name    => $seq_name,
                    start       => $ens_transcript->seq_region_start,
                    end         => $ens_transcript->seq_region_end,
                    strand      => $ens_transcript->seq_region_strand,
                    gene        => $gene,
                }
            );
            $gene->add_transcript($transcript);

            my $pos =
                $ens_transcript->seq_region_strand == 1
              ? $ens_transcript->seq_region_end
              : $ens_transcript->seq_region_start;

            push @{ $cache{ id $self}->{$seq_name}->{$pos} }, $transcript;
        }
    }

    return;
}

=method add_gene_annotation

  Usage       : my $regions_ref
                    = $gene_finder->add_gene_annotation($regions_ary_ref);
  Purpose     : Add gene annotation to regions with 3' ends
  Returns     : Arrayref [
                    Arrayref [
                        String (region sequence name),
                        Int (region start),
                        Int (region end),
                        Int (region maximum read count),
                        Float (region log probability sum),
                        String (3' end sequence name) or undef,
                        Int (3' end position) or undef,
                        Int (3' end strand) or undef,
                        Int (3' end read count) or undef,
                        Arrayref [
                            Int (count)
                            ...
                        ],
                        Arrayref [
                            Int (normalised count)
                            ...
                        ],
                        Int (p value) or undef,
                        Int (adjusted p value) or undef,
                        Arrayref [
                            Int (condition fold change) or undef,
                            Int (log2 condition fold change) or undef,
                        ],
                        Arrayref [
                            Arrayref [
                                Int (group fold change) or undef,
                                Int (log2 group fold change) or undef,
                            ],
                            ... (groups)
                        ],
                        Hashref {
                            String (genebuild version) => Arrayref [
                                Arrayref [
                                    String (gene stable id),
                                    String (gene name) or undef,
                                    String (gene description) or undef,
                                    String (gene biotype),
                                    Int (distance to 3' end),
                                    Arrayref [
                                        Arrayref [
                                            String (transcript stable id),
                                            String (transcript biotype),
                                        ],
                                        ... (transcripts)
                                    ],
                                ],
                                ... (genes)
                            ],
                        }
                    ],
                    ... (regions)
                }
  Parameters  : Arrayref (of regions)
  Throws      : If regions are missing
  Comments    : None

=cut

sub add_gene_annotation {
    my ( $self, $regions ) = @_;

    confess 'No regions specified' if !defined $regions;

    my @output;

    foreach my $region ( @{$regions} ) {

        # Get details for region and 3' end
        my $region_seq_name = $region->[0];
        my $region_start    = $region->[1];
        my $region_end      = $region->[2];
        ## no critic (ProhibitMagicNumbers)
        my $three_prime_seq_name = $region->[5];
        my $three_prime_pos      = $region->[6];
        my $three_prime_strand   = $region->[7];
        ## use critic

        my %gene_annotation = ();
        my $genes;
        my $distance;
        my $nearest_end_pos;

        if ( defined $three_prime_seq_name ) {

            # Find nearest genes to 3' end (taking strand into account)
            ( $genes, $distance, $nearest_end_pos ) =
              $self->get_nearest_genes( $three_prime_seq_name, $three_prime_pos,
                $three_prime_strand );
        }

        # Add annotation if got genes
        foreach my $gene ( @{$genes} ) {
            my @transcripts;
            foreach my $transcript ( @{ $gene->get_all_transcripts() } ) {

                # Only add those transcripts nearest to 3' end
                ## no critic (ProhibitMagicNumbers)
                if (
                    (
                           $transcript->strand == 1
                        && $transcript->end == $nearest_end_pos
                    )
                    || (   $transcript->strand == -1
                        && $transcript->start == $nearest_end_pos )
                  )
                {
                    ## use critic
                    push @transcripts,
                      [ $transcript->stable_id, $transcript->biotype, ];
                }
            }
            push @{ $gene_annotation{ $gene->genebuild_version } },
              [
                $gene->stable_id, $gene->name, $gene->description,
                $gene->biotype,   $distance,   \@transcripts,
              ];
        }

        push @{$region}, \%gene_annotation;
        push @output, $region;
    }

    return \@output;
}

1;
