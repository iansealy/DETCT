## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Gene;
## use critic

# ABSTRACT: Object representing a gene

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

use Readonly;
use Class::InsideOut qw( private register id );
use Scalar::Util qw( weaken );

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private genebuild_version => my %genebuild_version;   # e.g. e69
private stable_id         => my %stable_id;           # e.g. ENSDARG00000095747
private name              => my %name;                # e.g. cxc64
private description       => my %description;         # e.g. CXC chemokine 64...
private biotype           => my %biotype;             # e.g. protein_coding
private seq_name          => my %seq_name;            # e.g. 5
private start             => my %start;               # e.g. 40352744
private end               => my %end;                 # e.g. 40354399
private strand            => my %strand;              # e.g. 1
private transcript        => my %transcript;          # DETCT::Transcript

# Constants
Readonly our $MAX_NAME_LENGTH => 128;

=method new

  Usage       : my $gene = DETCT::Gene->new( {
                    genebuild_version => 'e61',
                    stable_id         => 'ENSDARG00000095747',
                    biotype           => 'protein_coding',
                    seq_name          => '5',
                    start             => 40352744,
                    end               => 40354399,
                    strand            => 1,
                } );
  Purpose     : Constructor for gene objects
  Returns     : DETCT::Gene
  Parameters  : Hashref {
                    genebuild_version => String,
                    stable_id         => String,
                    name              => String or undef,
                    description       => String or undef,
                    biotype           => String,
                    seq_name          => String,
                    start             => +ve Int,
                    end               => +ve Int,
                    strand            => Int (1 or -1),
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_genebuild_version( $arg_ref->{genebuild_version} );
    $self->set_stable_id( $arg_ref->{stable_id} );
    $self->set_name( $arg_ref->{name} );
    $self->set_description( $arg_ref->{description} );
    $self->set_biotype( $arg_ref->{biotype} );
    $self->set_seq_name( $arg_ref->{seq_name} );
    $self->set_start( $arg_ref->{start} );
    $self->set_end( $arg_ref->{end} );
    $self->set_strand( $arg_ref->{strand} );
    return $self;
}

=method genebuild_version

  Usage       : my $gv = $gene->genebuild_version;
  Purpose     : Getter for genebuild version attribute
  Returns     : String (e.g. "e61" for Ensembl 61)
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub genebuild_version {
    my ($self) = @_;
    return $genebuild_version{ id $self};
}

=method set_genebuild_version

  Usage       : $gene->set_genebuild_version('e61');
  Purpose     : Setter for genebuild version attribute
  Returns     : undef
  Parameters  : String (the genebuild version)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_genebuild_version {
    my ( $self, $arg ) = @_;
    $genebuild_version{ id $self} = check_genebuild_version($arg);
    return;
}

=method check_genebuild_version

  Usage       : $gv = check_genebuild_version($gv);
  Purpose     : Check for valid genebuild version
  Returns     : String (the valid genebuild version)
  Parameters  : String (the genebuild version)
  Throws      : If genebuild version is missing or invalid (i.e. not
                alphanumeric)
  Comments    : None

=cut

sub check_genebuild_version {
    my ($genebuild_version) = @_;
    return $genebuild_version
      if defined $genebuild_version && $genebuild_version =~ m/\A \w+ \z/xms;
    confess 'No genebuild version specified' if !defined $genebuild_version;
    confess "Invalid genebuild version ($genebuild_version) specified";
}

=method stable_id

  Usage       : my $stable_id = $gene->stable_id;
  Purpose     : Getter for stable id attribute
  Returns     : String (e.g. "ENSDARG00000095747")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub stable_id {
    my ($self) = @_;
    return $stable_id{ id $self};
}

=method set_stable_id

  Usage       : $gene->set_stable_id('ENSDARG00000095747');
  Purpose     : Setter for stable id attribute
  Returns     : undef
  Parameters  : String (the stable id)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_stable_id {
    my ( $self, $arg ) = @_;
    $stable_id{ id $self} = check_stable_id($arg);
    return;
}

=method check_stable_id

  Usage       : $stable_id = check_stable_id($stable_id);
  Purpose     : Check for valid stable id
  Returns     : String (the valid stable id)
  Parameters  : String (the stable id)
  Throws      : If stable id is missing or invalid
  Comments    : None

=cut

sub check_stable_id {
    my ($stable_id) = @_;
    return $stable_id
      if defined $stable_id && $stable_id =~ m/\A [[:upper:]]+ \d{11} \z/xms;
    confess 'No stable id specified' if !defined $stable_id;
    confess "Invalid stable id ($stable_id) specified";
}

=method name

  Usage       : my $name = $gene->name;
  Purpose     : Getter for name attribute
  Returns     : String (e.g. "cxc64")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub name {
    my ($self) = @_;
    return $name{ id $self};
}

=method set_name

  Usage       : $gene->set_name('cxc64');
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
# Throws      : If name > $MAX_NAME_LENGTH characters
# Comments    : None
sub _check_name {
    my ($name) = @_;
    return $name
      if !defined $name
      || ( length $name > 0 && length $name <= $MAX_NAME_LENGTH );
    confess 'Name is empty' if !length $name;
    confess "Name ($name) longer than $MAX_NAME_LENGTH characters";
}

=method description

  Usage       : my $description = $gene->description;
  Purpose     : Getter for description attribute
  Returns     : String (e.g. "CXC chemokine 64")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub description {
    my ($self) = @_;
    return $description{ id $self};
}

=method set_description

  Usage       : $gene->set_description('CXC chemokine 64');
  Purpose     : Setter for description attribute
  Returns     : undef
  Parameters  : String (the description)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_description {
    my ( $self, $arg ) = @_;
    $description{ id $self} = $arg;
    return;
}

=method biotype

  Usage       : my $biotype = $gene->biotype;
  Purpose     : Getter for biotype attribute
  Returns     : String (e.g. "protein_coding")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub biotype {
    my ($self) = @_;
    return $biotype{ id $self};
}

=method set_biotype

  Usage       : $gene->set_biotype('protein_coding');
  Purpose     : Setter for biotype attribute
  Returns     : undef
  Parameters  : String (the biotype)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_biotype {
    my ( $self, $arg ) = @_;
    $biotype{ id $self} = check_biotype($arg);
    return;
}

=method check_biotype

  Usage       : $biotype = check_biotype($biotype);
  Purpose     : Check for valid biotype
  Returns     : String (the valid biotype)
  Parameters  : String (the biotype)
  Throws      : If biotype is missing or invalid (i.e. not alphanumeric)
  Comments    : None

=cut

sub check_biotype {
    my ($biotype) = @_;
    return $biotype if defined $biotype && $biotype =~ m/\A \w+ \z/xms;
    confess 'No biotype specified' if !defined $biotype;
    confess "Invalid biotype ($biotype) specified";
}

=method seq_name

  Usage       : my $seq_name = $gene->seq_name;
  Purpose     : Getter for sequence name attribute
  Returns     : String (e.g. "5")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub seq_name {
    my ($self) = @_;
    return $seq_name{ id $self};
}

=method set_seq_name

  Usage       : $gene->set_seq_name('5');
  Purpose     : Setter for sequence name attribute
  Returns     : undef
  Parameters  : String (the sequence name)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_seq_name {
    my ( $self, $arg ) = @_;
    $seq_name{ id $self} = check_seq_name($arg);
    return;
}

=method check_seq_name

  Usage       : $seq_name = check_seq_name($seq_name);
  Purpose     : Check for valid sequence name
  Returns     : String (the valid sequence name)
  Parameters  : String (the sequence name)
  Throws      : If sequence name is missing or invalid (i.e. not alphanumeric)
  Comments    : None

=cut

sub check_seq_name {
    my ($seq_name) = @_;
    return $seq_name if defined $seq_name && $seq_name =~ m/\A \w+ \z/xms;
    confess 'No sequence name specified' if !defined $seq_name;
    confess "Invalid sequence name ($seq_name) specified";
}

=method start

  Usage       : my $start = $gene->start;
  Purpose     : Getter for start attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub start {
    my ($self) = @_;
    return $start{ id $self};
}

=method set_start

  Usage       : $gene->set_start(40352744);
  Purpose     : Setter for start attribute
  Returns     : undef
  Parameters  : +ve Int (the start)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_start {
    my ( $self, $arg ) = @_;
    $start{ id $self} = check_start($arg);
    return;
}

=method check_start

  Usage       : $start = check_start($start);
  Purpose     : Check for valid start
  Returns     : +ve Int (the valid start)
  Parameters  : +ve Int (the start)
  Throws      : If start is missing or not a positive integer
  Comments    : None

=cut

sub check_start {
    my ($start) = @_;
    return $start if defined $start && $start =~ m/\A \d+ \z/xms;
    confess 'No start specified' if !defined $start;
    confess "Invalid start ($start) specified";
}

=method end

  Usage       : my $end = $gene->end;
  Purpose     : Getter for end attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub end {
    my ($self) = @_;
    return $end{ id $self};
}

=method set_end

  Usage       : $gene->set_end(40352744);
  Purpose     : Setter for end attribute
  Returns     : undef
  Parameters  : +ve Int (the end)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_end {
    my ( $self, $arg ) = @_;
    $end{ id $self} = check_end($arg);
    return;
}

=method check_end

  Usage       : $end = check_end($end);
  Purpose     : Check for valid end
  Returns     : +ve Int (the valid end)
  Parameters  : +ve Int (the end)
  Throws      : If end is missing or not a positive integer
  Comments    : None

=cut

sub check_end {
    my ($end) = @_;
    return $end if defined $end && $end =~ m/\A \d+ \z/xms;
    confess 'No end specified' if !defined $end;
    confess "Invalid end ($end) specified";
}

=method strand

  Usage       : my $strand = $gene->strand;
  Purpose     : Getter for strand attribute
  Returns     : Int (1 or -1)
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub strand {
    my ($self) = @_;
    return $strand{ id $self};
}

=method set_strand

  Usage       : $gene->set_strand(1);
  Purpose     : Setter for strand attribute
  Returns     : undef
  Parameters  : Int (the strand)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_strand {
    my ( $self, $arg ) = @_;
    $strand{ id $self} = _check_strand($arg);
    return;
}

# Usage       : $strand = _check_strand($strand);
# Purpose     : Check for valid strand
# Returns     : Int (1 or -1) (the valid strand)
# Parameters  : Int (1 or -1) (the strand)
# Throws      : If strand is missing or not 1 or -1
# Comments    : None
sub _check_strand {
    my ($strand) = @_;
    return $strand if defined $strand && $strand =~ m/\A \-? 1 \z/xms;
    confess 'No strand specified' if !defined $strand;
    confess "Invalid strand ($strand) specified";
}

=method add_transcript

  Usage       : $gene->add_transcript($transcript);
  Purpose     : Add a transcript to a gene
  Returns     : undef
  Parameters  : DETCT::Transcript
  Throws      : If transcript is missing or invalid (i.e. not a
                DETCT::Transcript object)
  Comments    : None

=cut

sub add_transcript {
    my ( $self, $transcript ) = @_;

    confess 'No transcript specified' if !defined $transcript;
    confess 'Class of transcript (', ref $transcript, ') not DETCT::Transcript'
      if !$transcript->isa('DETCT::Transcript');

    weaken($transcript);    # Avoid circular references

    if ( !exists $transcript{ id $self} ) {
        $transcript{ id $self} = [$transcript];
    }
    else {
        push @{ $transcript{ id $self} }, $transcript;
    }

    return;
}

=method get_all_transcripts

  Usage       : $transcripts = $gene->get_all_transcripts();
  Purpose     : Get all transcripts of a gene
  Returns     : Arrayref of DETCT::Transcript objects
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_transcripts {
    my ($self) = @_;

    return $transcript{ id $self} || [];
}

1;
