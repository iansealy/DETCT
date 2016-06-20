## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::GeneOntologyTerm;
## VERSION
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
use Scalar::Util qw( blessed );
use DETCT::GeneOntologyEvidenceCode;

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private accession     => my %accession;     # e.g. GO:0005622
private namespace     => my %namespace;     # e.g. cellular_component
private name          => my %name;          # e.g. intracellular
private definition    => my %definition;    # e.g. The living contents of a cell
private gene          => my %gene;          # DETCT::Gene
private evidence_code => my %evidence_code; # e.g. IEA

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
    confess 'No name specified' if !defined $name;
    return $name;
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
    confess 'No definition specified' if !defined $definition;
    return $definition;
}

=method add_gene

  Usage       : $term->add_gene($gene);
  Purpose     : Add a gene to a Gene Ontology term
  Returns     : undef
  Parameters  : DETCT::Gene
  Throws      : If gene is missing or invalid (i.e. not a DETCT::Gene object)
  Comments    : None

=cut

sub add_gene {
    my ( $self, $gene ) = @_;

    confess 'No gene specified' if !defined $gene;
    confess 'Class of gene (', ref $gene, ') not DETCT::Gene'
      if !$gene->isa('DETCT::Gene');

    if ( !exists $gene{ id $self} ) {
        $gene{ id $self} = [$gene];
    }
    else {
        push @{ $gene{ id $self} }, $gene;
    }

    return;
}

=method get_all_genes

  Usage       : $genes = $term->get_all_genes();
  Purpose     : Get all genes for a Gene Ontology term
  Returns     : Arrayref of DETCT::Gene objects
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_genes {
    my ($self) = @_;

    return $gene{ id $self} || [];
}

=method add_evidence_code

  Usage       : $term->add_evidence_code('IEA', $gene);
  Purpose     : Add evidence code for a gene to a Gene Ontology term
  Returns     : undef
  Parameters  : String (the evidence code)
                DETCT::Gene or String (the stable id)
  Throws      : If evidence code is missing or invalid
  Comments    : None

=cut

sub add_evidence_code {
    my ( $self, $evidence_code, $gene_or_stable_id ) = @_;

    confess "Invalid evidence code ($evidence_code) specified"
      if !
      exists $DETCT::GeneOntologyEvidenceCode::DESCRIPTION_FOR_EVIDENCE_CODE{
        $evidence_code};

    my $stable_id = _check_gene_or_stable_id($gene_or_stable_id);

    $evidence_code{ id $self}{$stable_id} = $evidence_code;

    return;
}

=method get_evidence_code

  Usage       : $evidence_code = $term->get_evidence_code($gene);
  Purpose     : Get evidence code for a particular gene
  Returns     : String
  Parameters  : DETCT::Gene or String (the stable id)
  Throws      : No exceptions
  Comments    : None

=cut

sub get_evidence_code {
    my ( $self, $gene_or_stable_id ) = @_;

    my $stable_id = _check_gene_or_stable_id($gene_or_stable_id);

    return $evidence_code{ id $self}{$stable_id};
}

# Usage       : my $stable_id = _check_gene_or_stable_id($gene);
# Purpose     : Check for valid gene or stable id
# Returns     : String (the valid stable id)
# Parameters  : DETCT::Gene or String (the stable id)
# Throws      : If gene is missing or invalid (i.e. not a DETCT::Gene object or
#               not a stable id)
# Comments    : None
sub _check_gene_or_stable_id {
    my ($gene_or_stable_id) = @_;

    # Transform to stable id?
    if ( blessed($gene_or_stable_id) && $gene_or_stable_id->isa('DETCT::Gene') )
    {
        $gene_or_stable_id = $gene_or_stable_id->stable_id;
    }

    confess 'No gene specified' if !defined $gene_or_stable_id;
    confess 'Class of gene (', ref $gene_or_stable_id, ') not DETCT::Gene'
      if ref $gene_or_stable_id;
    confess "Invalid stable id ($gene_or_stable_id) specified"
      if $gene_or_stable_id !~ m/\A [[:upper:]]+ \d{11} \z/xms;

    return $gene_or_stable_id;
}

1;
