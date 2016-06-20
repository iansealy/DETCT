## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Sequence;
# VERSION
## use critic

# ABSTRACT: Object representing a sequence (a component of a reference sequence)

## Author         : is1
## Maintainer     : is1
## Created        : 2012-09-21
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Readonly;
use Class::InsideOut qw( private register id );

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private name => my %name;    # e.g. 1
private bp   => my %bp;      # e.g. 60348388

# Constants
Readonly our $MAX_NAME_LENGTH => 128;

=method new

  Usage       : my $sequence = DETCT::Sequence->new( {
                    name => '1',
                    bp   => 60_348_388,
                } );
  Purpose     : Constructor for sequence objects
  Returns     : DETCT::Sequence
  Parameters  : Hashref {
                    name => String,
                    bp   => Int,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_name( $arg_ref->{name} );
    $self->set_bp( $arg_ref->{bp} );
    return $self;
}

=method name

  Usage       : my $name = $sequence->name;
  Purpose     : Getter for name attribute
  Returns     : String (e.g. "1")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub name {
    my ($self) = @_;
    return $name{ id $self};
}

=method set_name

  Usage       : $sequence->set_name('1');
  Purpose     : Setter for name attribute
  Returns     : undef
  Parameters  : String (the name)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_name {
    my ( $self, $arg ) = @_;
    $name{ id $self} = _check_name($arg);
    return;
}

# Usage       : $name = _check_name($name);
# Purpose     : Check for valid name
# Returns     : String (the valid name)
# Parameters  : String (the name)
# Throws      : If name is missing
#               If name is empty
#               If name > $MAX_NAME_LENGTH characters
# Comments    : None
sub _check_name {
    my ($name) = @_;

    confess 'No name specified' if !defined $name;
    confess 'Empty name specified' if !length $name;
    confess "Name ($name) longer than $MAX_NAME_LENGTH characters"
      if length $name > $MAX_NAME_LENGTH;

    return $name;
}

=method bp

  Usage       : my $bp = $sequence->bp;
  Purpose     : Getter for bp attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub bp {
    my ($self) = @_;
    return $bp{ id $self};
}

=method set_bp

  Usage       : $sequence->set_bp(40352744);
  Purpose     : Setter for bp attribute
  Returns     : undef
  Parameters  : +ve Int (bp)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_bp {
    my ( $self, $arg ) = @_;
    $bp{ id $self} = _check_bp($arg);
    return;
}

# Usage       : $bp = _check_bp($bp);
# Purpose     : Check for valid bp
# Returns     : +ve Int (valid bp)
# Parameters  : +ve Int (bp)
# Throws      : If bp is missing or not a positive integer
# Comments    : None
sub _check_bp {
    my ($bp) = @_;
    return $bp
      if defined $bp && $bp =~ m/\A \d+ \z/xms;
    confess 'No bp specified' if !defined $bp;
    confess "Invalid bp ($bp) specified";
}

1;
