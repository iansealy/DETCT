## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Analysis::Downsample;
## use critic

# ABSTRACT: Object representing downsampling analysis

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-17
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use parent qw(DETCT::Analysis);

use Readonly;
use Class::InsideOut qw( private register id );
use List::MoreUtils qw( any );
use YAML::Tiny;

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private read_count_type => my %read_count_type;    # e.g. paired
private round_down_to   => my %round_down_to;      # e.g. 1000000

=method new

  Usage       : my $analysis = DETCT::Analysis::Downsample->new( {
                    name            => 'zmp_ph1',
                    read_count_type => 'paired',
                    round_down_to   => 1000000,
                    chunk_total     => 20,
                } );
  Purpose     : Constructor for analysis objects
  Returns     : DETCT::Analysis::Downsample
  Parameters  : Hashref {
                    name            => String,
                    read_count_type => String ('paired', 'mapped' or 'proper'),
                    round_down_to   => Int,
                    ensembl_host    => String or undef,
                    ensembl_port    => Int or undef,
                    ensembl_user    => String or undef,
                    ensembl_pass    => String or undef,
                    ensembl_name    => String or undef,
                    ensembl_species => String or undef,
                    chunk_total     => Int,
                    test_chunk      => Int or undef,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = $class->SUPER::new($arg_ref);
    $self->set_read_count_type( $arg_ref->{read_count_type} );
    $self->set_round_down_to( $arg_ref->{round_down_to} );
    return $self;
}

=method new_from_yaml

  Usage       : my $analysis
                    = DETCT::Analysis::DiffExpr->new_from_yaml('zmp_ph1.yaml');
  Purpose     : Constructor for creating analysis objects from a YAML file
  Returns     : DETCT::Analysis::Downsample
  Parameters  : String (the YAML file)
  Throws      : If YAML file is missing or not readable
  Comments    : None

=cut

sub new_from_yaml {
    my ( $class, $yaml_file ) = @_;
    my $self = $class->SUPER::new_from_yaml($yaml_file);

    confess "YAML file ($yaml_file) does not exist or cannot be read"
      if !-r $yaml_file;

    my $yaml = YAML::Tiny->read($yaml_file);

    if ( !$yaml ) {
        confess sprintf 'YAML file (%s) is invalid: %s', $yaml_file,
          YAML::Tiny->errstr;
    }

    $self->set_read_count_type( $yaml->[0]->{read_count_type} );
    $self->set_round_down_to( $yaml->[0]->{round_down_to} );

    return $self;
}

=method read_count_type

  Usage       : my $read_count_type = $analysis->read_count_type;
  Purpose     : Getter for read count type attribute
  Returns     : String ('paired', 'mapped', or 'proper')
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub read_count_type {
    my ($self) = @_;
    return $read_count_type{ id $self};
}

=method set_read_count_type

  Usage       : $analysis->set_read_count_type(20);
  Purpose     : Setter for read count type attribute
  Returns     : undef
  Parameters  : String (the read count type)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_read_count_type {
    my ( $self, $arg ) = @_;
    $read_count_type{ id $self} = _check_read_count_type($arg);
    return;
}

# Usage       : $read_count_type = _check_read_count_type($read_count_type);
# Purpose     : Check for valid read count type
# Returns     : String (the valid read count type)
# Parameters  : string (the read count type)
# Throws      : If read count type is missing or invalid
# Comments    : None
sub _check_read_count_type {
    my ($read_count_type) = @_;
    return $read_count_type
      if defined $read_count_type && any { $_ eq $read_count_type }
    qw(paired mapped proper);
    confess 'No read count type specified' if !defined $read_count_type;
    confess "Invalid read count type ($read_count_type) specified";
}

=method round_down_to

  Usage       : my $round_down_to = $analysis->round_down_to;
  Purpose     : Getter for round down to attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub round_down_to {
    my ($self) = @_;
    return $round_down_to{ id $self};
}

=method set_round_down_to

  Usage       : $analysis->set_round_down_to(1000000);
  Purpose     : Setter for round down to attribute
  Returns     : undef
  Parameters  : +ve Int (the read count to round down to)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_round_down_to {
    my ( $self, $arg ) = @_;
    $round_down_to{ id $self} = _check_round_down_to($arg);
    return;
}

# Usage       : $round_down_to = _check_round_down_to($round_down_to);
# Purpose     : Check for valid round down to
# Returns     : +ve Int (the valid round down to)
# Parameters  : +ve Int (the read count to round down to)
# Throws      : If round down to is missing or not a positive integer
# Comments    : None
sub _check_round_down_to {
    my ($round_down_to) = @_;
    return $round_down_to
      if defined $round_down_to && $round_down_to =~ m/\A \d+ \z/xms;
    confess 'No round down to specified' if !defined $round_down_to;
    confess "Invalid round down to ($round_down_to) specified";
}

1;
