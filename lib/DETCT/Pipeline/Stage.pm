## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Pipeline::Stage;
## VERSION
## use critic

# ABSTRACT: Object representing a pipeline stage

## Author         : is1
## Maintainer     : is1
## Created        : 2013-01-09
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

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private name           => my %name;              # e.g. count_tags
private default_memory => my %default_memory;    # e.g. 3000
private threads        => my %threads;           # e.g. 1
private all_jobs_run   => my %all_jobs_run;      # e.g. 1
private prerequisite   => my %prerequisite;      # arrayref of stages

=method new

  Usage       : my $stage = DETCT::Pipeline::Stage->new( {
                    name           => 'count_tags',
                    default_memory => 3000,
                } );
  Purpose     : Constructor for stage objects
  Returns     : DETCT::Pipeline::Stage
  Parameters  : Hashref {
                    name           => String,
                    default_memory => Int,
                    threads        => Int,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_name( $arg_ref->{name} );
    $self->set_default_memory( $arg_ref->{default_memory} );
    $self->set_threads( $arg_ref->{threads} );
    return $self;
}

=method name

  Usage       : my $name = $sample->name;
  Purpose     : Getter for name attribute
  Returns     : String (e.g. "count_tags")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub name {
    my ($self) = @_;
    return $name{ id $self};
}

=method set_name

  Usage       : $sample->set_name('count_tags');
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
# Throws      : If name is missing or invalid (i.e. not alphanumeric)
# Comments    : None
sub _check_name {
    my ($name) = @_;

    return $name if defined $name && $name =~ m/\A \w+ \z/xms;
    confess 'No name specified' if !defined $name;
    confess "Invalid name ($name) specified";
}

=method default_memory

  Usage       : my $default_memory = $stage->default_memory;
  Purpose     : Getter for default memory attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub default_memory {
    my ($self) = @_;
    return $default_memory{ id $self};
}

=method set_default_memory

  Usage       : $stage->set_default_memory(3000);
  Purpose     : Setter for default memory attribute
  Returns     : undef
  Parameters  : +ve Int (the default memory)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_default_memory {
    my ( $self, $arg ) = @_;
    $default_memory{ id $self} = _check_default_memory($arg);
    return;
}

# Usage       : $default_memory = _check_default_memory($default_memory);
# Purpose     : Check for valid default memory
# Returns     : +ve Int (the valid default memory)
# Parameters  : +ve Int (the default memory)
# Throws      : If default memory is missing or not a positive integer
# Comments    : None
sub _check_default_memory {
    my ($default_memory) = @_;
    return $default_memory
      if defined $default_memory && $default_memory =~ m/\A \d+ \z/xms;
    confess 'No default memory specified' if !defined $default_memory;
    confess "Invalid default memory ($default_memory) specified";
}

=method threads

  Usage       : my $threads = $stage->threads;
  Purpose     : Getter for threads attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub threads {
    my ($self) = @_;
    return $threads{ id $self} || 1;
}

=method set_threads

  Usage       : $stage->set_threads(2);
  Purpose     : Setter for threads attribute
  Returns     : undef
  Parameters  : +ve Int (the threads)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_threads {
    my ( $self, $arg ) = @_;
    $threads{ id $self} = _check_threads($arg);
    return;
}

# Usage       : $threads = _check_threads($threads);
# Purpose     : Check for valid threads
# Returns     : +ve Int (the valid threads)
# Parameters  : +ve Int (the threads)
# Throws      : If threads is not a positive integer
# Comments    : None
sub _check_threads {
    my ($threads) = @_;

    confess "Invalid threads ($threads) specified"
      if defined $threads && $threads !~ m/\A \d+ \z/xms;

    return $threads;
}

=method all_jobs_run

  Usage       : my $all_jobs_run = $stage->all_jobs_run;
  Purpose     : Getter for all jobs run flag
  Returns     : Boolean
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_jobs_run {
    my ($self) = @_;
    return $all_jobs_run{ id $self} || 0;
}

=method set_all_jobs_run

  Usage       : $stage->set_all_jobs_run(1);
  Purpose     : Setter for all jobs run flag
  Returns     : undef
  Parameters  : Boolean
  Throws      : No exceptions
  Comments    : None

=cut

sub set_all_jobs_run {
    my ( $self, $arg ) = @_;
    $all_jobs_run{ id $self} = $arg ? 1 : 0;
    return;
}

=method add_prerequisite

  Usage       : $stage->add_prerequisite($prerequisite);
  Purpose     : Add a prerequisite to a stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Stage
  Throws      : If prerequisite is missing or invalid (i.e. not a
                DETCT::Pipeline::Stage object)
  Comments    : None

=cut

sub add_prerequisite {
    my ( $self, $prerequisite ) = @_;

    confess 'No prerequisite specified' if !defined $prerequisite;
    confess 'Class of prerequisite (', ref $prerequisite,
      ') not DETCT::Pipeline::Stage'
      if !$prerequisite->isa('DETCT::Pipeline::Stage');

    if ( !exists $prerequisite{ id $self} ) {
        $prerequisite{ id $self} = [$prerequisite];
    }
    else {
        push @{ $prerequisite{ id $self} }, $prerequisite;
    }

    return;
}

=method get_all_prerequisites

  Usage       : $prerequisites = $stage->get_all_prerequisites();
  Purpose     : Get all prerequisites of a stage
  Returns     : Arrayref of DETCT::Pipeline::Stage objects
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_prerequisites {
    my ($self) = @_;

    return $prerequisite{ id $self} || [];
}

1;
