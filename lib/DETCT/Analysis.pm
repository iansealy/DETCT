## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Analysis;
## VERSION
## use critic

# ABSTRACT: Object representing an analysis of a collection of samples

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
use List::MoreUtils qw( uniq );
use Sort::Naturally;
use YAML::Tiny;
use Data::Compare;
use DETCT::Sample;
use DETCT::Sequence;
use DETCT::Misc::BAM;

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private name          => my %name;           # e.g. zmp_ph1
private skip_sequence => my %skip_sequence;  # hashref of skipped sequence names
private sample        => my %sample;         # arrayref of samples
private sequence      => my %sequence;       # arrayref of sequences
private total_bp      => my %total_bp;       # e.g. 1412464843
private ref_fasta     => my %ref_fasta;      # e.g. zv9.fa
private fasta_index   => my %fasta_index;    # Bio::DB::HTS::Fai
private ensembl_host  => my %ensembl_host;   # e.g. ensembldb.ensembl.org
private ensembl_port  => my %ensembl_port;   # e.g. 3306
private ensembl_user  => my %ensembl_user;   # e.g. anonymous
private ensembl_pass  => my %ensembl_pass;   # e.g. secret
private ensembl_name  => my %ensembl_name;   # e.g. zv9_core
private ensembl_species => my %ensembl_species;    # e.g. danio_rerio
private ensembl_db_type => my %ensembl_db_type;    # hashref of database types
private slice_adaptor => my %slice_adaptor; # Bio::EnsEMBL::DBSQL::SliceAdaptor
private chunk_total   => my %chunk_total;   # e.g. 20
private chunk         => my %chunk;         # arrayref of arrayrefs of sequences
private test_chunk    => my %test_chunk;    # e.g. 1

# Constants
Readonly our $MAX_NAME_LENGTH         => 128;
Readonly our $DEFAULT_ENSEMBL_HOST    => 'ensembldb.ensembl.org';
Readonly our $DEFAULT_ENSEMBL_USER    => 'anonymous';
Readonly our $DEFAULT_ENSEMBL_DB_TYPE => 'core';

=method new

  Usage       : my $analysis = DETCT::Analysis->new( {
                    name        => 'zmp_ph1',
                    chunk_total => 20,
                } );
  Purpose     : Constructor for analysis objects
  Returns     : DETCT::Analysis
  Parameters  : Hashref {
                    name            => String,
                    ref_fasta       => String or undef,
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
    my $self = register($class);
    $self->set_name( $arg_ref->{name} );
    $self->set_ref_fasta( $arg_ref->{ref_fasta} );
    $self->set_ensembl_host( $arg_ref->{ensembl_host} );
    $self->set_ensembl_port( $arg_ref->{ensembl_port} );
    $self->set_ensembl_user( $arg_ref->{ensembl_user} );
    $self->set_ensembl_pass( $arg_ref->{ensembl_pass} );
    $self->set_ensembl_name( $arg_ref->{ensembl_name} );
    $self->set_ensembl_species( $arg_ref->{ensembl_species} );
    $self->set_chunk_total( $arg_ref->{chunk_total} );
    $self->set_test_chunk( $arg_ref->{test_chunk} );
    return $self;
}

=method new_from_yaml

  Usage       : my $analysis = DETCT::Analysis->new_from_yaml( 'zmp_ph1.yaml' );
  Purpose     : Constructor for creating analysis objects from a YAML file
  Returns     : DETCT::Analysis
  Parameters  : String (the YAML file)
  Throws      : If YAML file is missing or not readable
  Comments    : None

=cut

sub new_from_yaml {
    my ( $class, $yaml_file ) = @_;
    my $self = register($class);

    confess "YAML file ($yaml_file) does not exist or cannot be read"
      if !-r $yaml_file;

    my $yaml = YAML::Tiny->read($yaml_file);

    if ( !$yaml ) {
        confess sprintf 'YAML file (%s) is invalid: %s', $yaml_file,
          YAML::Tiny->errstr;
    }

    $self->set_name( $yaml->[0]->{name} );
    $self->set_ref_fasta( $yaml->[0]->{ref_fasta} );
    $self->set_ensembl_host( $yaml->[0]->{ensembl_host} );
    $self->set_ensembl_port( $yaml->[0]->{ensembl_port} );
    $self->set_ensembl_user( $yaml->[0]->{ensembl_user} );
    $self->set_ensembl_pass( $yaml->[0]->{ensembl_pass} );
    $self->set_ensembl_name( $yaml->[0]->{ensembl_name} );
    $self->set_ensembl_species( $yaml->[0]->{ensembl_species} );
    $self->set_chunk_total( $yaml->[0]->{chunk_total} );
    $self->set_test_chunk( $yaml->[0]->{test_chunk} );

    $self->add_all_ensembl_db_types( $yaml->[0]->{ensembl_db_types} );

    # Need to add skip sequences before adding samples because add_sample()
    # calls add_all_sequences()
    $self->add_all_skip_sequences( $yaml->[0]->{skip_sequences} );

    foreach my $sample_hash ( @{ $yaml->[0]->{samples} } ) {
        my $sample = DETCT::Sample->new(
            {
                name        => $sample_hash->{name},
                description => $sample_hash->{description},
                condition   => $sample_hash->{condition},
                group       => $sample_hash->{group},
                tag         => $sample_hash->{tag},
                bam_file    => $sample_hash->{bam_file},
            }
        );
        $self->add_sample( $sample, 1 );    # 1 = do not validate
    }

    $self->validate();

    return $self;
}

=method name

  Usage       : my $name = $analysis->name;
  Purpose     : Getter for name attribute
  Returns     : String (e.g. "zmp_ph1")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub name {
    my ($self) = @_;
    return $name{ id $self};
}

=method set_name

  Usage       : $analysis->set_name('zmp_ph1');
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

=method add_all_skip_sequences

  Usage       : $analysis->add_all_skip_sequences(['MT']);
  Purpose     : Add all skip sequences to an analysis
  Returns     : undef
  Parameters  : Arrayref of strings (the skip sequences) or undef
  Throws      : No exceptions
  Comments    : None

=cut

sub add_all_skip_sequences {
    my ( $self, $skip_sequences ) = @_;

    $skip_sequence{ id $self} = {};

    foreach my $seq_name ( @{ $skip_sequences || [] } ) {
        $skip_sequence{ id $self}->{$seq_name} = 1;
    }

    return;
}

=method get_all_skip_sequences

  Usage       : $skip_sequences = $analysis->get_all_skip_sequences();
  Purpose     : Get all skip sequences of an analysis
  Returns     : Arrayref of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_skip_sequences {
    my ($self) = @_;

    return [ sort keys %{ $skip_sequence{ id $self} || {} } ];
}

=method is_skip_sequence

  Usage       : next if $analysis->is_skip_sequence('MT');
  Purpose     : Check if sequence should be skipped
  Returns     : 1 or 0
  Parameters  : String (the sequence name)
  Throws      : No exceptions
  Comments    : None

=cut

sub is_skip_sequence {
    my ( $self, $seq_name ) = @_;

    return exists $skip_sequence{ id $self}->{$seq_name} ? 1 : 0;
}

=method add_sample

  Usage       : $analysis->add_sample($sample);
  Purpose     : Add a sample to an analysis
  Returns     : undef
  Parameters  : DETCT::Sample
                Defined or undef (indicating if validation is needed)
  Throws      : If sample is missing or invalid (i.e. not a DETCT::Sample
                object)
  Comments    : None

=cut

sub add_sample {
    my ( $self, $sample, $no_validaton ) = @_;

    confess 'No sample specified' if !defined $sample;
    confess 'Class of sample (', ref $sample, ') not DETCT::Sample'
      if !$sample->isa('DETCT::Sample');

    if ( !exists $sample{ id $self} ) {
        $sample{ id $self} = [$sample];
        $self->add_all_sequences( $sample->bam_file );    # Because first sample
    }
    else {
        push @{ $sample{ id $self} }, $sample;
    }

    if ( !defined $no_validaton ) {
        $self->validate();
    }

    return;
}

=method get_all_samples

  Usage       : $samples = $analysis->get_all_samples();
  Purpose     : Get all samples of an analysis
  Returns     : Arrayref of DETCT::Sample objects
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_samples {
    my ($self) = @_;

    return $sample{ id $self} || [];
}

=method add_all_sequences

  Usage       : $analysis->add_all_sequences($bam_file);
  Purpose     : Add all sequences (sorted by decreasing length) to an analysis
  Returns     : undef
  Parameters  : String (the BAM file)
  Throws      : No exceptions
  Comments    : None

=cut

sub add_all_sequences {
    my ( $self, $bam_file ) = @_;

    undef $total_bp{ id $self};

    $bam_file = DETCT::Sample::check_bam_file($bam_file);

    $sequence{ id $self} = [];

    my %len =
      DETCT::Misc::BAM::get_reference_sequence_lengths( $bam_file,
        $self->get_all_skip_sequences() );

    foreach my $name ( reverse sort { $len{$a} <=> $len{$b} } keys %len ) {

        my $sequence = DETCT::Sequence->new(
            {
                name => $name,
                bp   => $len{$name},
            }
        );

        push @{ $sequence{ id $self} }, $sequence;
    }

    # Group sequences into chunks
    $self->add_all_chunks();

    return;
}

=method get_all_sequences

  Usage       : $sequences = $analysis->get_all_sequences();
  Purpose     : Get all sequences (sorted by decreasing length) of an analysis
  Returns     : Arrayref of DETCT::Sequence objects
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_sequences {
    my ($self) = @_;

    return $sequence{ id $self} || [];
}

=method validate

  Usage       : $analysis->validate();
  Purpose     : Check analysis
  Returns     : 1
  Parameters  : None
  Throws      : If reference sequences don't match
                If BAM files has no index
                If sample names are duplicated
                If multiple samples have the same tag and BAM file
                If samples do not all have same number of groups
                If group name is duplicated between groups
                If BAM file doesn't contain tag in read group
  Comments    : None

=cut

sub validate {
    my ($self) = @_;

    my @bam_files = $self->list_all_bam_files();

    # Compare reference sequence from first BAM file to all other BAM files
    my $first_bam_file = shift @bam_files;
    my %first_bam_length =
      DETCT::Misc::BAM::get_reference_sequence_lengths( $first_bam_file,
        $self->get_all_skip_sequences() );
    foreach my $bam_file (@bam_files) {
        my %bam_length =
          DETCT::Misc::BAM::get_reference_sequence_lengths( $bam_file,
            $self->get_all_skip_sequences() );
        if ( !Compare( \%first_bam_length, \%bam_length ) ) {
            confess "$first_bam_file and $bam_file use different reference";
        }
    }
    unshift @bam_files, $first_bam_file;

    # Check for missing BAM file index
    foreach my $bam_file (@bam_files) {
        if (
            !Bio::DB::HTSfile->index_load( Bio::DB::HTSfile->open($bam_file) ) )
        {
            confess "$bam_file has no index";
        }
    }

    # Check samples
    my ( %seen_name, %sample_bam_tag, %groups_per_sample, %group_ordinal );
    foreach my $sample ( @{ $self->get_all_samples } ) {

        # Check for duplicated sample names
        if ( exists $seen_name{ $sample->name } ) {
            confess 'Sample name (', $sample->name, ') is duplicated';
        }
        else {
            $seen_name{ $sample->name } = 1;
        }

        # Check for samples which have the same tag and BAM file
        if ( exists $sample_bam_tag{ $sample->bam_file }{ $sample->tag } ) {
            confess 'Multiple samples have the same tag (', $sample->tag,
              ') and BAM file (', $sample->bam_file, q{)};
        }
        else {
            $sample_bam_tag{ $sample->bam_file }{ $sample->tag } = 1;
        }

        # Check all samples have same number of groups
        $groups_per_sample{ scalar @{ $sample->groups } }++;

        # Check same group name not repeated in multiple groups
        my $ordinal = 0;
        foreach my $group ( @{ $sample->groups } ) {
            $ordinal++;
            $group_ordinal{$group}{$ordinal}++;
        }
    }

    # Check all samples have same number of groups
    if ( scalar keys %groups_per_sample > 1 ) {
        confess 'Samples do not all have same number of groups';
    }

    # Check same group name not repeated in multiple groups
    foreach my $group ( sort keys %group_ordinal ) {
        if ( scalar keys %{ $group_ordinal{$group} } > 1 ) {
            confess 'Group name (', $group, ') is duplicated between groups';
        }
    }

    # If BAM file contains read groups then check sample tags are present
    foreach my $bam_file (@bam_files) {
        my $bam        = Bio::DB::HTSfile->open( $bam_file, q{r} );
        my $header     = $bam->header_read->text;
        my %tag_in_bam = map { $_ => 1 }
          $header =~ m/\@RG .*? LB: [NRYKMSWBDHV]* ([AGCT]+)/xmsg;
        if ( keys %tag_in_bam ) {
            foreach my $tag ( keys %{ $sample_bam_tag{$bam_file} } ) {
                my ($index) = $tag =~ m/([AGCT]+) \z/xms;
                if ( !exists $tag_in_bam{$index} ) {
                    confess $bam_file, ' does not contain tag ', $tag;
                }
            }
        }
    }

    return 1;
}

=method total_bp

  Usage       : my $total_bp = $analysis->total_bp;
  Purpose     : Getter for total bp
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub total_bp {
    my ($self) = @_;

    if ( !defined $total_bp{ id $self} ) {
        my @seqs = @{ $self->get_all_sequences() };
        $total_bp{ id $self} = 0;
        foreach my $seq (@seqs) {
            $total_bp{ id $self} += $seq->bp;
        }
    }

    return $total_bp{ id $self};
}

=method ref_fasta

  Usage       : my $ref_fasta = $analysis->ref_fasta;
  Purpose     : Getter for reference FASTA attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ref_fasta {
    my ($self) = @_;
    return $ref_fasta{ id $self};
}

=method set_ref_fasta

  Usage       : $analysis->set_ref_fasta('zv9.fa');
  Purpose     : Setter for reference FASTA attribute
  Returns     : undef
  Parameters  : String (the reference FASTA)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ref_fasta {
    my ( $self, $arg ) = @_;
    $ref_fasta{ id $self} = _check_ref_fasta($arg);
    return;
}

# Usage       : $ref_fasta = _check_ref_fasta($ref_fasta);
# Purpose     : Check for valid reference FASTA
# Returns     : String (the valid reference FASTA)
# Parameters  : String (the reference FASTA)
# Throws      : If reference FASTA is defined but not readable
# Comments    : None
sub _check_ref_fasta {
    my ($ref_fasta) = @_;
    return $ref_fasta if !defined $ref_fasta || -r $ref_fasta;
    confess "Reference FASTA ($ref_fasta) cannot be read";
}

=method fasta_index

  Usage       : my $fai = $analysis->fasta_index;
  Purpose     : Getter for FASTA index attribute
  Returns     : Bio::DB::HTS::Fai
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub fasta_index {
    my ($self) = @_;

    if ( !defined $fasta_index{ id $self} && $self->ref_fasta ) {

        # We can create a FASTA index object
        $self->set_fasta_index( Bio::DB::HTS::Fai->load( $self->ref_fasta ) );
    }

    return $fasta_index{ id $self};
}

=method set_fasta_index

  Usage       : $analysis->set_fasta_index($fai);
  Purpose     : Setter for FASTA index attribute
  Returns     : undef
  Parameters  : Bio::DB::HTS::Fai
  Throws      : No exceptions
  Comments    : None

=cut

sub set_fasta_index {
    my ( $self, $arg ) = @_;
    $fasta_index{ id $self} = _check_fasta_index($arg);
    return;
}

# Usage       : $fai = _check_fasta_index($fai);
# Purpose     : Check for valid FASTA index
# Returns     : Bio::DB::HTS::Fai
# Parameters  : Bio::DB::HTS::Fai
# Throws      : If FASTA index is missing or invalid (i.e. not a
#               Bio::DB::HTS::Fai object)
# Comments    : None
sub _check_fasta_index {
    my ($fasta_index) = @_;
    return $fasta_index
      if defined $fasta_index
      && $fasta_index->isa('Bio::DB::HTS::Fai');
    confess 'No FASTA index specified' if !defined $fasta_index;
    confess 'Class of FASTA index (', ref $fasta_index,
      ') not Bio::DB::HTS::Fai';
}

=method ensembl_host

  Usage       : my $ensembl_host = $analysis->ensembl_host;
  Purpose     : Getter for Ensembl host attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ensembl_host {
    my ($self) = @_;
    return $ensembl_host{ id $self};
}

=method set_ensembl_host

  Usage       : $analysis->set_ensembl_host('ensembldb.ensembl.org');
  Purpose     : Setter for Ensembl host attribute
  Returns     : undef
  Parameters  : String (the Ensembl host)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ensembl_host {
    my ( $self, $arg ) = @_;
    $ensembl_host{ id $self} = $arg;
    return;
}

=method ensembl_port

  Usage       : my $ensembl_port = $analysis->ensembl_port;
  Purpose     : Getter for Ensembl port attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ensembl_port {
    my ($self) = @_;
    return $ensembl_port{ id $self};
}

=method set_ensembl_port

  Usage       : $analysis->set_ensembl_port(3306);
  Purpose     : Setter for Ensembl port attribute
  Returns     : undef
  Parameters  : +ve Int (the Ensembl port)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ensembl_port {
    my ( $self, $arg ) = @_;
    $ensembl_port{ id $self} = _check_ensembl_port($arg);
    return;
}

# Usage       : $ensembl_port = _check_ensembl_port($ensembl_port);
# Purpose     : Check for valid Ensembl port
# Returns     : +ve Int (the valid Ensembl port)
# Parameters  : +ve Int (the Ensembl port)
# Throws      : If Ensembl port is defined but not a positive integer
# Comments    : None
sub _check_ensembl_port {
    my ($ensembl_port) = @_;
    return $ensembl_port
      if !defined $ensembl_port || $ensembl_port =~ m/\A \d+ \z/xms;
    confess "Invalid Ensembl port ($ensembl_port) specified";
}

=method ensembl_user

  Usage       : my $ensembl_user = $analysis->ensembl_user;
  Purpose     : Getter for Ensembl username attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ensembl_user {
    my ($self) = @_;
    return $ensembl_user{ id $self};
}

=method set_ensembl_user

  Usage       : $analysis->set_ensembl_user('anonymous');
  Purpose     : Setter for Ensembl username attribute
  Returns     : undef
  Parameters  : String (the Ensembl username)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ensembl_user {
    my ( $self, $arg ) = @_;
    $ensembl_user{ id $self} = $arg;
    return;
}

=method ensembl_pass

  Usage       : my $ensembl_pass = $analysis->ensembl_pass;
  Purpose     : Getter for Ensembl password attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ensembl_pass {
    my ($self) = @_;
    return $ensembl_pass{ id $self};
}

=method set_ensembl_pass

  Usage       : $analysis->set_ensembl_pass('secret');
  Purpose     : Setter for Ensembl password attribute
  Returns     : undef
  Parameters  : String (the Ensembl password)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ensembl_pass {
    my ( $self, $arg ) = @_;
    $ensembl_pass{ id $self} = $arg;
    return;
}

=method ensembl_name

  Usage       : my $ensembl_name = $analysis->ensembl_name;
  Purpose     : Getter for Ensembl database name attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ensembl_name {
    my ($self) = @_;
    return $ensembl_name{ id $self};
}

=method set_ensembl_name

  Usage       : $analysis->set_ensembl_name('zv9_core');
  Purpose     : Setter for Ensembl database name attribute
  Returns     : undef
  Parameters  : String (the Ensembl database name)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ensembl_name {
    my ( $self, $arg ) = @_;
    $ensembl_name{ id $self} = $arg;
    return;
}

=method ensembl_species

  Usage       : my $ensembl_species = $analysis->ensembl_species;
  Purpose     : Getter for Ensembl species attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ensembl_species {
    my ($self) = @_;
    return $ensembl_species{ id $self};
}

=method set_ensembl_species

  Usage       : $analysis->set_ensembl_species('danio_rerio');
  Purpose     : Setter for Ensembl species attribute
  Returns     : undef
  Parameters  : String (the Ensembl species)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ensembl_species {
    my ( $self, $arg ) = @_;
    $ensembl_species{ id $self} = $arg;
    return;
}

=method add_all_ensembl_db_types

  Usage       : $analysis->add_all_ensembl_db_types(['core', 'otherfeatures']);
  Purpose     : Add all Ensembl database types to an analysis
  Returns     : undef
  Parameters  : Arrayref of strings (the Ensembl database types) or undef
  Throws      : No exceptions
  Comments    : None

=cut

sub add_all_ensembl_db_types {
    my ( $self, $ensembl_db_types ) = @_;

    $ensembl_db_type{ id $self} = {};

    foreach my $ensembl_db_type ( @{ $ensembl_db_types || [] } ) {
        $ensembl_db_type{ id $self}->{$ensembl_db_type} = 1;
    }

    return;
}

=method get_all_ensembl_db_types

  Usage       : $ensembl_db_types = $analysis->get_all_ensembl_db_types();
  Purpose     : Get all Ensembl database types of an analysis
  Returns     : Arrayref of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_ensembl_db_types {
    my ($self) = @_;

    return [
        sort keys %{ $ensembl_db_type{ id $self}
              || { $DEFAULT_ENSEMBL_DB_TYPE => 1 } }
    ];
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

    if ( !defined $slice_adaptor{ id $self}
        && ( $self->ensembl_species || $self->ensembl_name ) )
    {
        # We can create an Ensembl slice adaptor
        $self->_create_slice_adaptor();
    }

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
    confess 'No Ensembl slice adaptor specified'
      if !defined $slice_adaptor;
    confess 'Class of Ensembl slice adaptor (', ref $slice_adaptor,
      ') not Bio::EnsEMBL::DBSQL::SliceAdaptor';
}

=method chunk_total

  Usage       : my $chunk_total = $analysis->chunk_total;
  Purpose     : Getter for chunk total attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub chunk_total {
    my ($self) = @_;
    return $chunk_total{ id $self};
}

=method set_chunk_total

  Usage       : $analysis->set_chunk_total(20);
  Purpose     : Setter for chunk total attribute
  Returns     : undef
  Parameters  : +ve Int (the chunk total)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_chunk_total {
    my ( $self, $arg ) = @_;
    $chunk_total{ id $self} = _check_chunk_total($arg);

    # Recalculate chunks if necessary
    if ( scalar @{ $self->get_all_samples() } ) {
        $self->add_all_chunks();
    }

    return;
}

# Usage       : $chunk_total = _check_chunk_total($chunk_total);
# Purpose     : Check for valid chunk total
# Returns     : +ve Int (the valid chunk total)
# Parameters  : +ve Int (the chunk total)
# Throws      : If chunk total is missing or not a positive integer
# Comments    : None
sub _check_chunk_total {
    my ($chunk_total) = @_;
    return $chunk_total
      if defined $chunk_total && $chunk_total =~ m/\A \d+ \z/xms;
    confess 'No chunk total specified' if !defined $chunk_total;
    confess "Invalid chunk total ($chunk_total) specified";
}

=method test_chunk

  Usage       : my $test_chunk = $analysis->test_chunk;
  Purpose     : Getter for test chunk attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub test_chunk {
    my ($self) = @_;
    return $test_chunk{ id $self};
}

=method set_test_chunk

  Usage       : $analysis->set_test_chunk(1);
  Purpose     : Setter for test chunk attribute
  Returns     : undef
  Parameters  : +ve Int (the test chunk)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_test_chunk {
    my ( $self, $arg ) = @_;
    $test_chunk{ id $self} = $arg;
    return;
}

=method add_all_chunks

  Usage       : $analysis->add_all_chunks();
  Purpose     : Add all chunks (groups of sequences) to an analysis
  Returns     : undef
  Parameters  : None
  Throws      : No exceptions
  Comments    : Groups all sequences into a specific number of (roughly equally
                sized) chunks

=cut

sub add_all_chunks {
    my ($self) = @_;

    my $total_bp = $self->total_bp;

    # Get chunk target size (+ 1 to ensure slight overestimate)
    my $target_chunk_size = int( $total_bp / $self->chunk_total + 1 );

    my @chunks;
    my @chunk_size = map { 0 } 1 .. $self->chunk_total;

    # Iterate over sequences
    my @seqs = @{ $self->get_all_sequences() };
  SEQ: foreach my $seq (@seqs) {

        # Iterate over each chunk
        foreach my $chunk_index ( 0 .. $self->chunk_total - 1 ) {

            # Add sequence to chunk if there's room or if the chunk is empty
            if (   $chunk_size[$chunk_index] + $seq->bp <= $target_chunk_size
                || $chunk_size[$chunk_index] == 0 )
            {
                push @{ $chunks[$chunk_index] }, $seq;
                $chunk_size[$chunk_index] += $seq->bp;
                next SEQ;    # Next sequence
            }
        }

        # Sequence hasn't been added to a chunk, so add to chunk with most room
        my $roomy_chunk_index = 0;
        foreach my $chunk_index ( 0 .. $self->chunk_total - 1 ) {
            if ( $chunk_size[$chunk_index] < $chunk_size[$roomy_chunk_index] ) {
                $roomy_chunk_index = $chunk_index;
            }
        }
        push @{ $chunks[$roomy_chunk_index] }, $seq;
        $chunk_size[$roomy_chunk_index] += $seq->bp;
    }

    # Iterate over empty chunks in order to attempt to add sequences to them
    foreach my $empty_chunk_index ( 0 .. $self->chunk_total - 1 ) {
        next
          if defined $chunks[$empty_chunk_index];    # Only want empty chunks

        # Find chunk with highest number of sequences (but more than one)
        my $max_seqs_chunk_index;
        my $max_seqs;
        foreach my $chunk_index ( 0 .. $self->chunk_total - 1 ) {
            next
              if !defined $chunks[$chunk_index];    # Only want non-empty chunks
            my $seqs = scalar @{ $chunks[$chunk_index] };
            if ( $seqs > 1
                && ( !defined $max_seqs || $seqs > $max_seqs ) )
            {
                $max_seqs_chunk_index = $chunk_index;
                $max_seqs             = $seqs;
            }
        }

        last if !defined $max_seqs;                 # No splittable chunks

        # Split chosen chunk into empty chunk
        my $split_index = int( $max_seqs / 2 );
        @{ $chunks[$empty_chunk_index] } =
          splice @{ $chunks[$max_seqs_chunk_index] }, 0, $split_index;
    }

    $chunk{ id $self} = \@chunks;

    # Number of chunks may be smaller than requested chunk total, so adjust
    $chunk_total{ id $self} = scalar @chunks;

    return;
}

=method get_all_chunks

  Usage       : $chunks = $analysis->get_all_chunks();
  Purpose     : Get all chunks (groups of sequences) of an analysis
  Returns     : Arrayref of arrayrefs of DETCT::Sequence objects
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_chunks {
    my ($self) = @_;

    my $chunks = $chunk{ id $self} || [];

    # If a test chunk is specified then only return that chunk not all chunks
    if ( $self->test_chunk
        && exists $chunks->[ $self->test_chunk - 1 ] )
    {
        $chunks = [ $chunks->[ $self->test_chunk - 1 ] ];
    }

    return $chunks;
}

=method list_all_bam_files

  Usage       : @bam_files = $analysis->list_all_bam_files();
  Purpose     : Get all BAM files used in an analysis
  Returns     : Arrayref of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub list_all_bam_files {
    my ($self) = @_;

    my $samples = $self->get_all_samples();

    my @bam_files = map { $_->bam_file } @{$samples};

    return uniq( sort @bam_files );
}

=method list_all_tags_by_bam_file

  Usage       : @tags = $analysis->list_all_tags_by_bam_file($bam_file);
  Purpose     : Get all tags used in an analysis in a particular BAM file
  Returns     : Arrayref of strings
  Parameters  : String (the BAM file)
  Throws      : No exceptions
  Comments    : None

=cut

sub list_all_tags_by_bam_file {
    my ( $self, $bam_file ) = @_;

    my $samples = $self->get_all_samples();

    my @tags =
      map { $_->tag } grep { $_->bam_file eq $bam_file } @{$samples};

    return uniq( sort @tags );
}

=method get_sample_name_by_bam_file_and_tag

  Usage       : $name = $analysis->get_sample_name_by_bam_file_and_tag(
                    $bam_file, $tag);
  Purpose     : Get sample name corresponding to a particular BAM file and tag
  Returns     : String (sample name)
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_sample_name_by_bam_file_and_tag {
    my ( $self, $bam_file, $tag ) = @_;

    my $samples = $self->get_all_samples();

    my ($name) =
      map { $_->name }
      grep { $_->bam_file eq $bam_file && $_->tag eq $tag } @{$samples};

    return $name;
}

=method get_sample_names_by_bam_file

  Usage       : @names = $analysis->get_sample_names_by_bam_file( $bam_file);
  Purpose     : Get sample names corresponding to a particular BAM file
  Returns     : Arrayref of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_sample_names_by_bam_file {
    my ( $self, $bam_file, $tag ) = @_;

    my $samples = $self->get_all_samples();

    my @names =
      map { $_->name } grep { $_->bam_file eq $bam_file } @{$samples};

    return uniq( sort @names );
}

=method list_all_conditions

  Usage       : @conditions = $analysis->list_all_conditions();
  Purpose     : Get all conditions used in an analysis
  Returns     : Array of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub list_all_conditions {
    my ($self) = @_;

    my $samples = $self->get_all_samples();

    my @conditions = map { $_->condition } @{$samples};

    return uniq( nsort(@conditions) );
}

=method list_all_groups

  Usage       : @groups = $analysis->list_all_groups();
  Purpose     : Get all groups used in an analysis
  Returns     : Array of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub list_all_groups {
    my ($self) = @_;

    my $samples = $self->get_all_samples();

    my @groups = map { @{ $_->groups } } @{$samples};

    return uniq( nsort(@groups) );
}

=method get_subsequence

  Usage       : $seq = $analysis->get_subsequence('1', 1, 10);
  Purpose     : Get subsequence from reference
  Returns     : String (sequence)
  Parameters  : String (the sequence name)
                Int (the sequence start)
                Int (the sequence end)
                Int (the sequence strand)
  Throws      : If sequence name is missing
                If sequence start is missing
                If sequence end is missing
                If sequence strand is missing
  Comments    : None

=cut

sub get_subsequence {
    my ( $self, $seq_name, $start, $end, $strand ) = @_;

    confess 'No sequence name specified'   if !defined $seq_name;
    confess 'No sequence start specified'  if !defined $start;
    confess 'No sequence end specified'    if !defined $end;
    confess 'No sequence strand specified' if !defined $strand;

    # Avoid negative positions (but don't worry if end is larger than sequence)
    if ( $start < 1 ) {
        $start = 1;
    }
    if ( $end < 1 ) {
        $end = 1;
    }

    my $subseq;

    if ( $self->fasta_index ) {
        $subseq = DETCT::Misc::BAM::get_sequence(
            {
                fasta_index => $self->fasta_index,
                seq_name    => $seq_name,
                start       => $start,
                end         => $end,
                strand      => $strand,
            }
        );
    }
    elsif ( $self->slice_adaptor ) {
        $subseq =
          $self->slice_adaptor->fetch_by_region( 'toplevel',
            $seq_name, $start, $end, $strand )->seq;
    }
    else {
        confess 'No reference FASTA or Ensembl database';
    }

    return uc $subseq;
}

# Usage       : $self->_create_slice_adaptor();
# Purpose     : Create an Ensembl slice adaptor
# Returns     : undef
# Parameters  : None
# Throws      : No exceptions
# Comments    : None
sub _create_slice_adaptor {
    my ($self) = @_;

    my $host =
        $self->ensembl_host
      ? $self->ensembl_host
      : $DEFAULT_ENSEMBL_HOST;
    my $port = $self->ensembl_port;
    my $user =
        $self->ensembl_user
      ? $self->ensembl_user
      : $DEFAULT_ENSEMBL_USER;
    my $pass = $self->ensembl_pass;
    my $slice_adaptor;
    if ( !$self->ensembl_name ) {

        # Get slice adaptor via registry
        require Bio::EnsEMBL::Registry;
        Bio::EnsEMBL::Registry->load_registry_from_db(
            -host    => $host,
            -port    => $port,
            -user    => $user,
            -pass    => $pass,
            -species => $self->ensembl_species,
        );
        $slice_adaptor =
          Bio::EnsEMBL::Registry->get_adaptor( $self->ensembl_species,
            'core', 'slice' );
    }
    else {
        # Get slice adaptor from specific database
        require Bio::EnsEMBL::DBSQL::DBAdaptor;
        my $ensembl_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
            -host   => $host,
            -port   => $port,
            -user   => $user,
            -pass   => $pass,
            -dbname => $self->ensembl_name,
        );
        $slice_adaptor = $ensembl_db->get_SliceAdaptor();
    }

    $self->set_slice_adaptor($slice_adaptor);

    return;
}

1;
