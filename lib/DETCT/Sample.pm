## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Sample;
## use critic

# ABSTRACT: Object representing a sample

## Author         : is1
## Maintainer     : is1
## Created        : 2012-09-19
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
private name        => my %name;           # e.g. zmp_ph1_1m
private description => my %description;    # e.g. ZMP phenotype 1.1 mutant
private condition   => my %condition;      # e.g. mutant
private group       => my %group;          # e.g. 1
private tag         => my %tag;            # e.g. NNNNBGAGGC
private bam_file    => my %bam_file;       # e.g. 8295_6#1.bam

# Constants
Readonly our $MAX_NAME_LENGTH      => 128;
Readonly our $MAX_CONDITION_LENGTH => 128;
Readonly our $MAX_GROUP_LENGTH     => 128;

=method new

  Usage       : my $sample = DETCT::Sample->new( {
                    name      => 'zmp_ph1_1m',
                    condition => 'mutant',
                    group     => '1',
                    tag       => 'NNNNBGAGGC',
                    bam_file  => '8295_6#1.bam',
                } );
  Purpose     : Constructor for sample objects
  Returns     : DETCT::Sample
  Parameters  : Hashref {
                    name        => String,
                    description => String or undef,
                    condition   => String,
                    group       => String or undef,
                    tag         => String,
                    bam_file    => String,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_name( $arg_ref->{name} );
    $self->set_description( $arg_ref->{description} );
    $self->set_condition( $arg_ref->{condition} );
    $self->set_group( $arg_ref->{group} );
    $self->set_tag( $arg_ref->{tag} );
    $self->set_bam_file( $arg_ref->{bam_file} );
    return $self;
}

=method name

  Usage       : my $name = $sample->name;
  Purpose     : Getter for name attribute
  Returns     : String (e.g. "zmp_ph1_1m")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub name {
    my ($self) = @_;
    return $name{ id $self};
}

=method set_name

  Usage       : $sample->set_name('zmp_ph1_1m');
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
#               If name is invalid (i.e. not alphanumeric)
#               If name is empty
#               If name > $MAX_NAME_LENGTH characters
# Comments    : None
sub _check_name {
    my ($name) = @_;

    confess 'No name specified'      if !defined $name;
    confess 'Empty name specified'   if !length $name;
    confess 'Invalid name specified' if $name !~ m/\A [\w.-]+ \z/xms;
    confess "Name ($name) longer than $MAX_NAME_LENGTH characters"
      if length $name > $MAX_NAME_LENGTH;

    return $name;
}

=method description

  Usage       : my $description = $sample->description;
  Purpose     : Getter for description attribute
  Returns     : String (e.g. "ZMP phenotype 1.1 mutant")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub description {
    my ($self) = @_;
    return $description{ id $self};
}

=method set_description

  Usage       : $sample->set_description('ZMP phenotype 1.1 mutant');
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

=method condition

  Usage       : my $condition = $sample->condition;
  Purpose     : Getter for condition attribute
  Returns     : String (e.g. "mutant")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub condition {
    my ($self) = @_;
    return $condition{ id $self};
}

=method set_condition

  Usage       : $sample->set_condition('mutant');
  Purpose     : Setter for condition attribute
  Returns     : undef
  Parameters  : String (the condition)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_condition {
    my ( $self, $arg ) = @_;
    $condition{ id $self} = _check_condition($arg);
    return;
}

# Usage       : $condition = _check_condition($condition);
# Purpose     : Check for valid condition
# Returns     : String (the valid condition)
# Parameters  : String (the condition)
# Throws      : If condition is missing
#               If condition is empty
#               If condition > $MAX_CONDITION_LENGTH characters
# Comments    : None
sub _check_condition {
    my ($condition) = @_;

    confess 'No condition specified' if !defined $condition;
    confess 'Empty condition specified' if !length $condition;
    confess
      "Condition ($condition) longer than $MAX_CONDITION_LENGTH characters"
      if length $condition > $MAX_CONDITION_LENGTH;

    return $condition;
}

=method group

  Usage       : my $group = $sample->group;
  Purpose     : Getter for group attribute
  Returns     : String (e.g. "1")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub group {
    my ($self) = @_;
    return $group{ id $self};
}

=method set_group

  Usage       : $sample->set_group('1');
  Purpose     : Setter for group attribute
  Returns     : undef
  Parameters  : String (the group)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_group {
    my ( $self, $arg ) = @_;
    $group{ id $self} = _check_group($arg);
    return;
}

# Usage       : $group = _check_group($group);
# Purpose     : Check for valid group
# Returns     : String (the valid group)
# Parameters  : String (the group)
# Throws      : If group is empty
#               If group > $MAX_GROUP_LENGTH characters
# Comments    : None
sub _check_group {
    my ($group) = @_;

    confess 'Empty group specified' if defined $group && !length $group;
    confess "Group ($group) longer than $MAX_GROUP_LENGTH characters"
      if defined $group && length $group > $MAX_GROUP_LENGTH;

    return $group;
}

=method tag

  Usage       : my $tag = $sample->tag;
  Purpose     : Getter for tag attribute
  Returns     : String (e.g. "NNNNBGAGGC")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub tag {
    my ($self) = @_;
    return $tag{ id $self};
}

=method set_tag

  Usage       : $sample->set_tag('NNNNBGAGGC');
  Purpose     : Setter for tag attribute
  Returns     : undef
  Parameters  : String (the tag)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_tag {
    my ( $self, $arg ) = @_;
    $tag{ id $self} = _check_tag($arg);
    return;
}

# Usage       : $tag = _check_tag($tag);
# Purpose     : Check for valid tag
# Returns     : String (the valid tag)
# Parameters  : String (the tag)
# Throws      : If tag is missing or invalid
# Comments    : None
sub _check_tag {
    my ($tag) = @_;
    return $tag
      if defined $tag && $tag =~ m/\A [NRYKMSWBDHV]+ [AGCT]+ \z/xms;
    confess 'No tag specified' if !defined $tag;
    confess "Invalid tag ($tag) specified";
}

=method bam_file

  Usage       : my $bam_file = $sample->bam_file;
  Purpose     : Getter for BAM file attribute
  Returns     : String (e.g. "8295_6#1.bam")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub bam_file {
    my ($self) = @_;
    return $bam_file{ id $self};
}

=method set_bam_file

  Usage       : $sample->set_bam('8295_6#1.bam');
  Purpose     : Setter for BAM file attribute
  Returns     : undef
  Parameters  : String (the BAM file)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_bam_file {
    my ( $self, $arg ) = @_;
    $bam_file{ id $self} = check_bam_file($arg);
    return;
}

=method check_bam_file

  Usage       : $bam_file = check_bam_file($bam_file);
  Purpose     : Check for valid BAM file
  Returns     : String (the valid BAM file)
  Parameters  : String (the BAM file)
  Throws      : If BAM file is missing or not readable
  Comments    : None

=cut

sub check_bam_file {
    my ($bam_file) = @_;
    return $bam_file if defined $bam_file && -r $bam_file;
    confess 'No BAM file specified' if !defined $bam_file;
    confess "BAM file ($bam_file) does not exist or cannot be read";
}

1;
