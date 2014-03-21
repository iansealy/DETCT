## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::GeneOntologyTerm;
## use critic

# ABSTRACT: Object representing a GO term

## Author         : is1
## Maintainer     : is1
## Created        : 2014-03-21
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
use Scalar::Util qw( weaken );

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private accession  => my %accession;     # e.g. GO:0005622
private namespace  => my %namespace;     # e.g. cellular_component
private name       => my %name;          # e.g. intracellular
private definition => my %definition;    # e.g. The living contents of a cell

=method new

  Usage       : my $term = DETCT::GeneOntologyTerm->new( {
                    accession  => 'GO:0005622',
                    namespace  => 'cellular_component',
                    name       => 'intracellular',
                    definition => 'The living contents of a cell',
                } );
  Purpose     : Constructor for GO term objects
  Returns     : DETCT::GeneOntologyTerm
  Parameters  : Hashref {
                    accession  => String,
                    namespace  => String,
                    name       => String,
                    definition => String,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_accession( $arg_ref->{accession} );
    $self->set_namespace( $arg_ref->{namespace} );
    $self->set_name( $arg_ref->{name} );
    $self->set_definition( $arg_ref->{definition} );
    return $self;
}

=method accession

  Usage       : my $accession = $term->accession;
  Purpose     : Getter for accession attribute
  Returns     : String (e.g. "GO:0005622")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub accession {
    my ($self) = @_;
    return $accession{ id $self};
}

=method set_accession

  Usage       : $term->set_accession('GO:0005622');
  Purpose     : Setter for accession attribute
  Returns     : undef
  Parameters  : String (the accession)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_accession {
    my ( $self, $arg ) = @_;
    $accession{ id $self} = check_accession($arg);
    return;
}

=method check_accession

  Usage       : $accession = check_accession($accession);
  Purpose     : Check for valid accession
  Returns     : String (the valid accession)
  Parameters  : String (the accession)
  Throws      : If accession is missing or invalid
  Comments    : None

=cut

sub check_accession {
    my ($accession) = @_;
    return $accession
      if defined $accession && $accession =~ m/\A GO: \d{7} \z/xms;
    confess 'No accession specified' if !defined $accession;
    confess "Invalid accession ($accession) specified";
}

=method namespace

  Usage       : my $namespace = $term->namespace;
  Purpose     : Getter for namespace attribute
  Returns     : String (e.g. "cellular_component")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub namespace {
    my ($self) = @_;
    return $namespace{ id $self};
}

=method set_namespace

  Usage       : $term->set_namespace('cellular_component');
  Purpose     : Setter for namespace attribute
  Returns     : undef
  Parameters  : String (the namespace)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_namespace {
    my ( $self, $arg ) = @_;
    $namespace{ id $self} = check_namespace($arg);
    return;
}

=method check_namespace

  Usage       : $namespace = check_namespace($namespace);
  Purpose     : Check for valid namespace
  Returns     : String (the valid namespace)
  Parameters  : String (the namespace)
  Throws      : If namespace is missing or invalid
  Comments    : None

=cut

sub check_namespace {
    my ($namespace) = @_;
    return $namespace
      if defined $namespace
      && ( $namespace eq 'cellular_component'
        || $namespace eq 'biological_process'
        || $namespace eq 'molecular_function' );
    confess 'No namespace specified' if !defined $namespace;
    confess "Invalid namespace ($namespace) specified";
}

=method name

  Usage       : my $name = $term->name;
  Purpose     : Getter for name attribute
  Returns     : String (e.g. "intracellular")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub name {
    my ($self) = @_;
    return $name{ id $self};
}

=method set_name

  Usage       : $term->set_name('intracellular');
  Purpose     : Setter for name attribute
  Returns     : undef
  Parameters  : String (the name)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_name {
    my ( $self, $arg ) = @_;
    $name{ id $self} = check_name($arg);
    return;
}

=method check_name

  Usage       : $name = check_name($name);
  Purpose     : Check for valid name
  Returns     : String (the valid name)
  Parameters  : String (the name)
  Throws      : If name is missing
  Comments    : None

=cut

sub check_name {
    my ($name) = @_;
    return $name if defined $name;
    confess 'No name specified' if !defined $name;
}

=method definition

  Usage       : my $definition = $term->definition;
  Purpose     : Getter for definition attribute
  Returns     : String (e.g. "The living contents of a cell")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub definition {
    my ($self) = @_;
    return $definition{ id $self};
}

=method set_definition

  Usage       : $term->set_definition('The living contents of a cell');
  Purpose     : Setter for definition attribute
  Returns     : undef
  Parameters  : String (the definition)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_definition {
    my ( $self, $arg ) = @_;
    $definition{ id $self} = check_definition($arg);
    return;
}

=method check_definition

  Usage       : $definition = check_definition($definition);
  Purpose     : Check for valid definition
  Returns     : String (the valid definition)
  Parameters  : String (the definition)
  Throws      : If definition is missing
  Comments    : None

=cut

sub check_definition {
    my ($definition) = @_;
    return $definition if defined $definition;
    confess 'No definition specified' if !defined $definition;
}

1;
