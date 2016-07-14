## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::TransposonFinder;
## VERSION
## use critic

# ABSTRACT: Object for finding transposons by location

## Author         : is1
## Maintainer     : is1
## Created        : 2016-07-14
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
use Number::Closest::XS qw( find_closest_numbers );

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private slice_adaptor => my %slice_adaptor;  # Bio::EnsEMBL::DBSQL::SliceAdaptor
private cache         => my %cache;          # hashref

=method new

  Usage       : my $transposon_finder = DETCT::TransposonFinder->new( {
                    slice_adaptor => $slice_adaptor,
                } );
  Purpose     : Constructor for transposon finder objects
  Returns     : DETCT::TransposonFinder
  Parameters  : Hashref {
                    slice_adaptor    => Bio::EnsEMBL::DBSQL::SliceAdaptor,
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

  Usage       : my $slice_adaptor = $transposon_finder->slice_adaptor;
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

  Usage       : $transposon_finder->set_slice_adaptor($slice_adaptor);
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

=method get_nearest_transposon

  Usage       : $finder->get_nearest_transposon($seq_name, $pos, $strand);
  Purpose     : Retrieve the nearest transposon to a 3' end
  Returns     : Int (distance)
                Int (nearest transposon position)
  Parameters  : String (the 3' end sequence name)
                Int (the 3' end position)
                Int (the 3' end strand)
  Throws      : No exceptions
  Comments    : Distance is positive if downstream of 3' end and negative if
                upstream

=cut

sub get_nearest_transposon {
    my ( $self, $seq_name, $pos, $strand ) = @_;

    # Ensure cache is filled
    $self->_fill_cache_from_ensembl($seq_name);

    my $transposon_positions = $cache{ id $self}->{$seq_name};

    my $distance;
    my $nearest_transposon_pos;

    my $nearest = find_closest_numbers( $pos, $transposon_positions );
    if ( @{$nearest} ) {
        $nearest_transposon_pos = $nearest->[0];
    }
    $distance = $pos - $nearest_transposon_pos;
    if ( $strand > 0 ) {
        $distance *= -1;    ## no critic (ProhibitMagicNumbers)
    }

    return $distance, $nearest_transposon_pos;
}

sub _fill_cache_from_ensembl {
    my ( $self, $seq_name ) = @_;

    # Skip if cache already filled
    return if exists $cache{ id $self}->{$seq_name};

    # Make sure default key exists (in case there are no transposons)
    $cache{ id $self}->{$seq_name} = [];

    my $slice = $self->slice_adaptor->fetch_by_region( 'toplevel', $seq_name );
    return if !defined $slice;    # No transposons if non-existent sequence name

    my %ends;
    my $repeat_features = $slice->get_all_RepeatFeatures();
    foreach my $repeat_feature ( @{$repeat_features} ) {
        next
          if $repeat_feature->repeat_consensus->repeat_type !~ m/Transposon/xms;
        $ends{ $repeat_feature->seq_region_start } = 1;
        $ends{ $repeat_feature->seq_region_end }   = 1;
    }
    $cache{ id $self}->{$seq_name} = [ sort { $a <=> $b } keys %ends ];

    return;
}

1;
