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
private target_read_count => my %target_read_count;    # e.g. 15000000
private read_count_type   => my %read_count_type;      # e.g. paired
private round_down_to     => my %round_down_to;        # e.g. 1000000
private samtools_binary   => my %samtools_binary;      # e.g. samtools
private java_binary       => my %java_binary;          # e.g. java
private mark_duplicates_jar =>
  my %mark_duplicates_jar;                             # e.g. MarkDuplicates.jar
private merge_sam_files_jar => my %merge_sam_files_jar; # e.g. MergeSamFiles.jar
private mark_duplicates_method => my %mark_duplicates_method;    # e.g. picard

=method new

  Usage       : my $analysis = DETCT::Analysis::Downsample->new( {
                    name                => 'zmp_ph1',
                    target_read_count   => 15000000,
                    read_count_type     => 'paired',
                    round_down_to       => 1000000,
                    samtools_binary     => 'samtools',
                    java_binary         => 'java',
                    merge_sam_files_jar => 'picard/MergeSamFiles.jar',
                    chunk_total         => 20,
                } );
  Purpose     : Constructor for analysis objects
  Returns     : DETCT::Analysis::Downsample
  Parameters  : Hashref {
                    name                   => String,
                    target_read_count      => Int or undef,
                    read_count_type        => String (paired, mapped or proper),
                    round_down_to          => Int or undef,
                    samtools_binary        => String,
                    java_binary            => String,
                    mark_duplicates_jar    => String,
                    merge_sam_files_jar    => String,
                    mark_duplicates_method => String (native or picard),
                    ensembl_host           => String or undef,
                    ensembl_port           => Int or undef,
                    ensembl_user           => String or undef,
                    ensembl_pass           => String or undef,
                    ensembl_name           => String or undef,
                    ensembl_species        => String or undef,
                    chunk_total            => Int,
                    test_chunk             => Int or undef,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = $class->SUPER::new($arg_ref);
    $self->set_target_read_count( $arg_ref->{target_read_count} );
    $self->set_read_count_type( $arg_ref->{read_count_type} );
    $self->set_round_down_to( $arg_ref->{round_down_to} );
    $self->set_samtools_binary( $arg_ref->{samtools_binary} );
    $self->set_java_binary( $arg_ref->{java_binary} );
    $self->set_mark_duplicates_jar( $arg_ref->{mark_duplicates_jar} );
    $self->set_merge_sam_files_jar( $arg_ref->{merge_sam_files_jar} );
    $self->set_mark_duplicates_method( $arg_ref->{mark_duplicates_method} );
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

    $self->set_target_read_count( $yaml->[0]->{target_read_count} );
    $self->set_read_count_type( $yaml->[0]->{read_count_type} );
    $self->set_round_down_to( $yaml->[0]->{round_down_to} );
    $self->set_samtools_binary( $yaml->[0]->{samtools_binary} );
    $self->set_java_binary( $yaml->[0]->{java_binary} );
    $self->set_mark_duplicates_jar( $yaml->[0]->{mark_duplicates_jar} );
    $self->set_merge_sam_files_jar( $yaml->[0]->{merge_sam_files_jar} );
    $self->set_mark_duplicates_method( $yaml->[0]->{mark_duplicates_method} );

    return $self;
}

=method target_read_count

  Usage       : my $target_read_count = $analysis->target_read_count;
  Purpose     : Getter for target read count attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub target_read_count {
    my ($self) = @_;
    return $target_read_count{ id $self};
}

=method set_target_read_count

  Usage       : $analysis->set_target_read_count(15000000);
  Purpose     : Setter for target read count attribute
  Returns     : undef
  Parameters  : +ve Int (the target read count)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_target_read_count {
    my ( $self, $arg ) = @_;
    $target_read_count{ id $self} = _check_target_read_count($arg);
    return;
}

# Usage       : $target_read_count
#                   = _check_target_read_count($target_read_count);
# Purpose     : Check for valid target read count
# Returns     : +ve Int (the valid target read count)
# Parameters  : +ve Int (the target read count)
# Throws      : If target read count is defined but not a positive integer
# Comments    : None
sub _check_target_read_count {
    my ($target_read_count) = @_;
    return $target_read_count
      if !defined $target_read_count || $target_read_count =~ m/\A \d+ \z/xms;
    confess "Invalid target read count ($target_read_count) specified";
}

=method read_count_type

  Usage       : my $read_count_type = $analysis->read_count_type;
  Purpose     : Getter for read count type attribute
  Returns     : String ('paired', 'mapped' or 'proper')
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub read_count_type {
    my ($self) = @_;
    return $read_count_type{ id $self};
}

=method set_read_count_type

  Usage       : $analysis->set_read_count_type('paired');
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
# Parameters  : String (the read count type)
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
# Throws      : If round down to is defined but not a positive integer
# Comments    : None
sub _check_round_down_to {
    my ($round_down_to) = @_;
    return $round_down_to
      if !defined $round_down_to || $round_down_to =~ m/\A \d+ \z/xms;
    confess "Invalid round down to ($round_down_to) specified";
}

=method samtools_binary

  Usage       : my $samtools_binary = $analysis->samtools_binary;
  Purpose     : Getter for SAMtools binary attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub samtools_binary {
    my ($self) = @_;
    return $samtools_binary{ id $self};
}

=method set_samtools_binary

  Usage       : $analysis->set_samtools_binary('samtools');
  Purpose     : Setter for SAMtools binary attribute
  Returns     : undef
  Parameters  : String (the SAMtools binary)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_samtools_binary {
    my ( $self, $arg ) = @_;
    $samtools_binary{ id $self} = _check_samtools_binary($arg);
    return;
}

# Usage       : $samtools_binary = _check_samtools_binary($samtools_binary);
# Purpose     : Check for valid SAMtools binary
# Returns     : String (the valid SAMtools binary)
# Parameters  : String (the SAMtools binary)
# Throws      : If SAMtools binary is missing
# Comments    : None
sub _check_samtools_binary {
    my ($samtools_binary) = @_;
    return $samtools_binary if defined $samtools_binary;
    confess 'No SAMtools binary specified';
}

=method java_binary

  Usage       : my $java_binary = $analysis->java_binary;
  Purpose     : Getter for Java binary attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub java_binary {
    my ($self) = @_;
    return $java_binary{ id $self};
}

=method set_java_binary

  Usage       : $analysis->set_java_binary('java');
  Purpose     : Setter for Java binary attribute
  Returns     : undef
  Parameters  : String (the Java binary)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_java_binary {
    my ( $self, $arg ) = @_;
    $java_binary{ id $self} = _check_java_binary($arg);
    return;
}

# Usage       : $java_binary = _check_java_binary($java_binary);
# Purpose     : Check for valid Java binary
# Returns     : String (the valid Java binary)
# Parameters  : String (the Java binary)
# Throws      : If Java binary is missing
# Comments    : None
sub _check_java_binary {
    my ($java_binary) = @_;
    return $java_binary if defined $java_binary;
    confess 'No Java binary specified';
}

=method mark_duplicates_jar

  Usage       : my $mark_duplicates_jar = $analysis->mark_duplicates_jar;
  Purpose     : Getter for MarkDuplicates JAR attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub mark_duplicates_jar {
    my ($self) = @_;
    return $mark_duplicates_jar{ id $self};
}

=method set_mark_duplicates_jar

  Usage       : $analysis->set_mark_duplicates_jar('picard/MarkDuplicates.jar');
  Purpose     : Setter for MarkDuplicates JAR attribute
  Returns     : undef
  Parameters  : String (the MarkDuplicates JAR)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_mark_duplicates_jar {
    my ( $self, $arg ) = @_;
    $mark_duplicates_jar{ id $self} = _check_mark_duplicates_jar($arg);
    return;
}

# Usage       : $mark_duplicates_jar
#                   = _check_mark_duplicates_jar($mark_duplicates_jar);
# Purpose     : Check for valid MarkDuplicates JAR
# Returns     : String (the valid MarkDuplicates JAR)
# Parameters  : String (the MarkDuplicates JAR)
# Throws      : If MarkDuplicates JAR is defined but not readable
# Comments    : None
sub _check_mark_duplicates_jar {
    my ($mark_duplicates_jar) = @_;
    return $mark_duplicates_jar
      if !defined $mark_duplicates_jar || -r $mark_duplicates_jar;
    confess sprintf 'MarkDuplicates JAR (%s) cannot be read',
      $mark_duplicates_jar;
}

=method merge_sam_files_jar

  Usage       : my $merge_sam_files_jar = $analysis->merge_sam_files_jar;
  Purpose     : Getter for MergeSamFiles JAR attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub merge_sam_files_jar {
    my ($self) = @_;
    return $merge_sam_files_jar{ id $self};
}

=method set_merge_sam_files_jar

  Usage       : $analysis->set_merge_sam_files_jar('picard/MergeSamFiles.jar');
  Purpose     : Setter for MergeSamFiles JAR attribute
  Returns     : undef
  Parameters  : String (the MergeSamFiles JAR)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_merge_sam_files_jar {
    my ( $self, $arg ) = @_;
    $merge_sam_files_jar{ id $self} = _check_merge_sam_files_jar($arg);
    return;
}

# Usage       : $merge_sam_files_jar
#                   = _check_merge_sam_files_jar($merge_sam_files_jar);
# Purpose     : Check for valid MergeSamFiles JAR
# Returns     : String (the valid MergeSamFiles JAR)
# Parameters  : String (the MergeSamFiles JAR)
# Throws      : If MergeSamFiles JAR is missing or not readable
# Comments    : None
sub _check_merge_sam_files_jar {
    my ($merge_sam_files_jar) = @_;
    return $merge_sam_files_jar
      if defined $merge_sam_files_jar && -r $merge_sam_files_jar;
    confess 'No MergeSamFiles JAR specified' if !defined $merge_sam_files_jar;
    confess sprintf 'MergeSamFiles JAR (%s) does not exist or cannot be read',
      $merge_sam_files_jar;
}

=method mark_duplicates_method

  Usage       : my $mark_duplicates_method = $analysis->mark_duplicates_method;
  Purpose     : Getter for mark duplicates method attribute
  Returns     : String ('native' or 'picard')
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub mark_duplicates_method {
    my ($self) = @_;
    return $mark_duplicates_method{ id $self};
}

=method set_mark_duplicates_method

  Usage       : $analysis->set_mark_duplicates_method('picard');
  Purpose     : Setter for mark duplicates method attribute
  Returns     : undef
  Parameters  : String (the mark duplicates method)
  Throws      : No exceptions
  Comments    : Defaults to native

=cut

sub set_mark_duplicates_method {
    my ( $self, $arg ) = @_;
    $mark_duplicates_method{ id $self} =
      _check_mark_duplicates_method( $arg || 'native' );
    return;
}

# Usage       : $mark_duplicates_method =
#                   _check_mark_duplicates_method($mark_duplicates_method);
# Purpose     : Check for valid mark duplicates method
# Returns     : String (the valid mark duplicates method)
# Parameters  : String (the mark duplicates method)
# Throws      : If mark duplicates method is invalid
# Comments    : None
sub _check_mark_duplicates_method {
    my ($mark_duplicates_method) = @_;
    return $mark_duplicates_method
      if any { $_ eq $mark_duplicates_method } qw(native picard);
    confess
      "Invalid mark duplicates method ($mark_duplicates_method) specified";
}

1;
