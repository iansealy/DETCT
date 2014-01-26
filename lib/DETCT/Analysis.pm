## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Analysis;
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
use YAML::Tiny;
use Data::Compare;
use DETCT::Sample;
use DETCT::Sequence;
use DETCT::Misc::BAM;

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private name               => my %name;               # e.g. zmp_ph1
private sample             => my %sample;             # arrayref of samples
private sequence           => my %sequence;           # arrayref of sequences
private read1_length       => my %read1_length;       # e.g. 30
private read2_length       => my %read2_length;       # e.g. 54
private mismatch_threshold => my %mismatch_threshold; # e.g. 2
private bin_size           => my %bin_size;           # e.g. 100
private peak_buffer_width  => my %peak_buffer_width;  # e.g. 100
private hmm_sig_level      => my %hmm_sig_level;      # e.g. 0.001
private hmm_binary         => my %hmm_binary;         # e.g. chiphmmnew
private r_binary           => my %r_binary;           # e.g. R
private deseq_script       => my %deseq_script;       # e.g. ~/run_deseq.R
private output_sig_level   => my %output_sig_level;   # e.g. 0.05
private ref_fasta          => my %ref_fasta;          # e.g. zv9.fa
private fasta_index        => my %fasta_index;        # Bio::DB::Sam::Fai
private ensembl_host    => my %ensembl_host;       # e.g. ensembldb.ensembl.org
private ensembl_port    => my %ensembl_port;       # e.g. 3306
private ensembl_user    => my %ensembl_user;       # e.g. anonymous
private ensembl_pass    => my %ensembl_pass;       # e.g. secret
private ensembl_name    => my %ensembl_name;       # e.g. zv9_core
private ensembl_species => my %ensembl_species;    # e.g. danio_rerio
private slice_adaptor => my %slice_adaptor; # Bio::EnsEMBL::DBSQL::SliceAdaptor
private chunk_total   => my %chunk_total;   # e.g. 20
private chunk         => my %chunk;         # arrayref of arrayrefs of sequences
private test_chunk    => my %test_chunk;    # e.g. 1

# Constants
Readonly our $MAX_NAME_LENGTH      => 128;
Readonly our $DEFAULT_ENSEMBL_HOST => 'ensembldb.ensembl.org';
Readonly our $DEFAULT_ENSEMBL_USER => 'anonymous';

=method new

  Usage       : my $analysis = DETCT::Analysis->new( {
                    name               => 'zmp_ph1',
                    read1_length       => 30,
                    read2_length       => 54,
                    mismatch_threshold => 2,
                    bin_size           => 100,
                    peak_buffer_width  => 100,
                    hmm_sig_level      => 0.001,
                    hmm_binary         => 'chiphmmnew',
                    r_binary           => 'R',
                    deseq_script       => 'script/run_deseq.R',
                    output_sig_level   => 0.05,
                    chunk_total        => 20,
                } );
  Purpose     : Constructor for analysis objects
  Returns     : DETCT::Analysis
  Parameters  : Hashref {
                    name               => String,
                    read1_length       => Int,
                    read2_length       => Int,
                    mismatch_threshold => Int,
                    bin_size           => Int,
                    peak_buffer_width  => Int,
                    hmm_sig_level      => Float,
                    hmm_binary         => String,
                    r_binary           => String,
                    deseq_script       => String,
                    output_sig_level   => Float,
                    ref_fasta          => String or undef,
                    ensembl_host       => String or undef,
                    ensembl_port       => Int or undef,
                    ensembl_user       => String or undef,
                    ensembl_pass       => String or undef,
                    ensembl_name       => String or undef,
                    ensembl_species    => String or undef,
                    chunk_total        => Int,
                    test_chunk         => Int or undef,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_name( $arg_ref->{name} );
    $self->set_read1_length( $arg_ref->{read1_length} );
    $self->set_read2_length( $arg_ref->{read2_length} );
    $self->set_mismatch_threshold( $arg_ref->{mismatch_threshold} );
    $self->set_bin_size( $arg_ref->{bin_size} );
    $self->set_peak_buffer_width( $arg_ref->{peak_buffer_width} );
    $self->set_hmm_sig_level( $arg_ref->{hmm_sig_level} );
    $self->set_hmm_binary( $arg_ref->{hmm_binary} );
    $self->set_r_binary( $arg_ref->{r_binary} );
    $self->set_deseq_script( $arg_ref->{deseq_script} );
    $self->set_output_sig_level( $arg_ref->{output_sig_level} );
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

    $self->set_name( $yaml->[0]->{name} );
    $self->set_read1_length( $yaml->[0]->{read1_length} );
    $self->set_read2_length( $yaml->[0]->{read2_length} );
    $self->set_mismatch_threshold( $yaml->[0]->{mismatch_threshold} );
    $self->set_bin_size( $yaml->[0]->{bin_size} );
    $self->set_peak_buffer_width( $yaml->[0]->{peak_buffer_width} );
    $self->set_hmm_sig_level( $yaml->[0]->{hmm_sig_level} );
    $self->set_hmm_binary( $yaml->[0]->{hmm_binary} );
    $self->set_r_binary( $yaml->[0]->{r_binary} );
    $self->set_deseq_script( $yaml->[0]->{deseq_script} );
    $self->set_output_sig_level( $yaml->[0]->{output_sig_level} );
    $self->set_ref_fasta( $yaml->[0]->{ref_fasta} );
    $self->set_ensembl_host( $yaml->[0]->{ensembl_host} );
    $self->set_ensembl_port( $yaml->[0]->{ensembl_port} );
    $self->set_ensembl_user( $yaml->[0]->{ensembl_user} );
    $self->set_ensembl_pass( $yaml->[0]->{ensembl_pass} );
    $self->set_ensembl_name( $yaml->[0]->{ensembl_name} );
    $self->set_ensembl_species( $yaml->[0]->{ensembl_species} );
    $self->set_chunk_total( $yaml->[0]->{chunk_total} );
    $self->set_test_chunk( $yaml->[0]->{test_chunk} );

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

    $self->_validate();

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
        $self->_validate();
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

    $bam_file = DETCT::Sample::check_bam_file($bam_file);

    $sequence{ id $self} = [];

    my %len = DETCT::Misc::BAM::get_reference_sequence_lengths($bam_file);

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

# Usage       : $analysis->_validate();
# Purpose     : Check analysis
# Returns     : 1
# Parameters  : None
# Throws      : If reference sequences don't match
# Comments    : None
sub _validate {
    my ($self) = @_;

    my @bam_files = $self->list_all_bam_files();

    # Compare reference sequence from first BAM file to all other BAM files
    my $first_bam_file = shift @bam_files;
    my %first_bam_length =
      DETCT::Misc::BAM::get_reference_sequence_lengths($first_bam_file);
    foreach my $bam_file (@bam_files) {
        my %bam_length =
          DETCT::Misc::BAM::get_reference_sequence_lengths($bam_file);
        if ( !Compare( \%first_bam_length, \%bam_length ) ) {
            confess "$first_bam_file and $bam_file use different reference";
        }
    }

    return 1;
}

=method read1_length

  Usage       : my $read1_length = $analysis->read1_length;
  Purpose     : Getter for read 1 length attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub read1_length {
    my ($self) = @_;
    return $read1_length{ id $self};
}

=method set_read1_length

  Usage       : $analysis->set_read1_length(20);
  Purpose     : Setter for read 1 length attribute
  Returns     : undef
  Parameters  : +ve Int (the read 1 length)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_read1_length {
    my ( $self, $arg ) = @_;
    $read1_length{ id $self} = _check_read1_length($arg);
    return;
}

# Usage       : $read1_length = _check_read1_length($read1_length);
# Purpose     : Check for valid read 1 length
# Returns     : +ve Int (the valid read 1 length)
# Parameters  : +ve Int (the read 1 length)
# Throws      : If read 1 length is missing or not a positive integer
# Comments    : None
sub _check_read1_length {
    my ($read1_length) = @_;
    return $read1_length
      if defined $read1_length && $read1_length =~ m/\A \d+ \z/xms;
    confess 'No read 1 length specified' if !defined $read1_length;
    confess "Invalid read 1 length ($read1_length) specified";
}

=method read2_length

  Usage       : my $read2_length = $analysis->read2_length;
  Purpose     : Getter for read 2 length attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub read2_length {
    my ($self) = @_;
    return $read2_length{ id $self};
}

=method set_read2_length

  Usage       : $analysis->set_read2_length(20);
  Purpose     : Setter for read 2 length attribute
  Returns     : undef
  Parameters  : +ve Int (the read 2 length)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_read2_length {
    my ( $self, $arg ) = @_;
    $read2_length{ id $self} = _check_read2_length($arg);
    return;
}

# Usage       : $read2_length = _check_read2_length($read2_length);
# Purpose     : Check for valid read 2 length
# Returns     : +ve Int (the valid read 2 length)
# Parameters  : +ve Int (the read 2 length)
# Throws      : If read 2 length is missing or not a positive integer
# Comments    : None
sub _check_read2_length {
    my ($read2_length) = @_;
    return $read2_length
      if defined $read2_length && $read2_length =~ m/\A \d+ \z/xms;
    confess 'No read 2 length specified' if !defined $read2_length;
    confess "Invalid read 2 length ($read2_length) specified";
}

=method mismatch_threshold

  Usage       : my $mismatch_threshold = $analysis->mismatch_threshold;
  Purpose     : Getter for mismatch threshold attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub mismatch_threshold {
    my ($self) = @_;
    return $mismatch_threshold{ id $self};
}

=method set_mismatch_threshold

  Usage       : $analysis->set_mismatch_threshold(20);
  Purpose     : Setter for mismatch threshold attribute
  Returns     : undef
  Parameters  : +ve Int (the mismatch threshold)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_mismatch_threshold {
    my ( $self, $arg ) = @_;
    $mismatch_threshold{ id $self} = _check_mismatch_threshold($arg);
    return;
}

# Usage       : $mismatch_threshold
#                   = _check_mismatch_threshold($mismatch_threshold);
# Purpose     : Check for valid mismatch threshold
# Returns     : +ve Int (the valid mismatch threshold)
# Parameters  : +ve Int (the mismatch threshold)
# Throws      : If mismatch threshold is missing or not a positive integer
# Comments    : None
sub _check_mismatch_threshold {
    my ($mismatch_threshold) = @_;
    return $mismatch_threshold
      if defined $mismatch_threshold && $mismatch_threshold =~ m/\A \d+ \z/xms;
    confess 'No mismatch threshold specified' if !defined $mismatch_threshold;
    confess "Invalid mismatch threshold ($mismatch_threshold) specified";
}

=method bin_size

  Usage       : my $bin_size = $analysis->bin_size;
  Purpose     : Getter for bin size attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub bin_size {
    my ($self) = @_;
    return $bin_size{ id $self};
}

=method set_bin_size

  Usage       : $analysis->set_bin_size(100);
  Purpose     : Setter for bin size attribute
  Returns     : undef
  Parameters  : +ve Int (the bin size)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_bin_size {
    my ( $self, $arg ) = @_;
    $bin_size{ id $self} = _check_bin_size($arg);
    return;
}

# Usage       : $bin_size = _check_bin_size($bin_size);
# Purpose     : Check for valid bin size
# Returns     : +ve Int (the valid bin size)
# Parameters  : +ve Int (the bin size)
# Throws      : If bin size is missing or not a positive integer
# Comments    : None
sub _check_bin_size {
    my ($bin_size) = @_;
    return $bin_size
      if defined $bin_size && $bin_size =~ m/\A \d+ \z/xms;
    confess 'No bin size specified' if !defined $bin_size;
    confess "Invalid bin size ($bin_size) specified";
}

=method peak_buffer_width

  Usage       : my $peak_buffer_width = $analysis->peak_buffer_width;
  Purpose     : Getter for peak buffer width attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub peak_buffer_width {
    my ($self) = @_;
    return $peak_buffer_width{ id $self};
}

=method set_peak_buffer_width

  Usage       : $analysis->set_peak_buffer_width(100);
  Purpose     : Setter for peak buffer width attribute
  Returns     : undef
  Parameters  : +ve Int (the peak buffer width)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_peak_buffer_width {
    my ( $self, $arg ) = @_;
    $peak_buffer_width{ id $self} = _check_peak_buffer_width($arg);
    return;
}

# Usage       : $peak_buffer_width = _check_peak_buffer_width($peak_buffer_width);
# Purpose     : Check for valid peak buffer width
# Returns     : +ve Int (the valid peak buffer width)
# Parameters  : +ve Int (the peak buffer width)
# Throws      : If peak buffer width is missing or not a positive integer
# Comments    : None
sub _check_peak_buffer_width {
    my ($peak_buffer_width) = @_;
    return $peak_buffer_width
      if defined $peak_buffer_width && $peak_buffer_width =~ m/\A \d+ \z/xms;
    confess 'No peak buffer width specified' if !defined $peak_buffer_width;
    confess "Invalid peak buffer width ($peak_buffer_width) specified";
}

=method hmm_sig_level

  Usage       : my $hmm_sig_level = $analysis->hmm_sig_level;
  Purpose     : Getter for HMM significance level attribute
  Returns     : +ve Float
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub hmm_sig_level {
    my ($self) = @_;
    return $hmm_sig_level{ id $self};
}

=method set_hmm_sig_level

  Usage       : $analysis->set_hmm_sig_level(0.001);
  Purpose     : Setter for HMM significance level attribute
  Returns     : undef
  Parameters  : +ve Float (the HMM significance level)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_hmm_sig_level {
    my ( $self, $arg ) = @_;
    $hmm_sig_level{ id $self} = _check_hmm_sig_level($arg);
    return;
}

# Usage       : $hmm_sig_level = _check_hmm_sig_level($hmm_sig_level);
# Purpose     : Check for valid HMM significance level
# Returns     : +ve Float (the valid HMM significance level)
# Parameters  : +ve Float (the HMM significance level)
# Throws      : If HMM significance level is missing or not a positive float
# Comments    : None
sub _check_hmm_sig_level {
    my ($hmm_sig_level) = @_;
    return $hmm_sig_level
      if defined $hmm_sig_level && $hmm_sig_level =~ m/\A \d* [.] \d+ \z/xms;
    confess 'No HMM significance level specified' if !defined $hmm_sig_level;
    confess "Invalid HMM significance level ($hmm_sig_level) specified";
}

=method hmm_binary

  Usage       : my $hmm_binary = $analysis->hmm_binary;
  Purpose     : Getter for HMM binary attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub hmm_binary {
    my ($self) = @_;
    return $hmm_binary{ id $self};
}

=method set_hmm_binary

  Usage       : $analysis->set_hmm_binary('chiphmmnew');
  Purpose     : Setter for HMM binary attribute
  Returns     : undef
  Parameters  : String (the HMM binary)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_hmm_binary {
    my ( $self, $arg ) = @_;
    $hmm_binary{ id $self} = _check_hmm_binary($arg);
    return;
}

# Usage       : $hmm_binary = _check_hmm_binary($hmm_binary);
# Purpose     : Check for valid HMM binary
# Returns     : String (the valid HMM binary)
# Parameters  : String (the HMM binary)
# Throws      : If HMM binary is missing
# Comments    : None
sub _check_hmm_binary {
    my ($hmm_binary) = @_;
    return $hmm_binary if defined $hmm_binary;
    confess 'No HMM binary specified';
}

=method r_binary

  Usage       : my $r_binary = $analysis->r_binary;
  Purpose     : Getter for R binary attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub r_binary {
    my ($self) = @_;
    return $r_binary{ id $self};
}

=method set_r_binary

  Usage       : $analysis->set_r_binary('R');
  Purpose     : Setter for R binary attribute
  Returns     : undef
  Parameters  : String (the R binary)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_r_binary {
    my ( $self, $arg ) = @_;
    $r_binary{ id $self} = _check_r_binary($arg);
    return;
}

# Usage       : $r_binary = _check_r_binary($r_binary);
# Purpose     : Check for valid R binary
# Returns     : String (the valid R binary)
# Parameters  : String (the R binary)
# Throws      : If R binary is missing
# Comments    : None
sub _check_r_binary {
    my ($r_binary) = @_;
    return $r_binary if defined $r_binary;
    confess 'No R binary specified';
}

=method deseq_script

  Usage       : my $deseq_script = $analysis->deseq_script;
  Purpose     : Getter for DESeq script attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub deseq_script {
    my ($self) = @_;
    return $deseq_script{ id $self};
}

=method set_deseq_script

  Usage       : $analysis->set_deseq_script('script/run_deseq.R');
  Purpose     : Setter for DESeq script attribute
  Returns     : undef
  Parameters  : String (the DESeq script)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_deseq_script {
    my ( $self, $arg ) = @_;
    $deseq_script{ id $self} = _check_deseq_script($arg);
    return;
}

# Usage       : $deseq_script = _check_deseq_script($deseq_script);
# Purpose     : Check for valid DESeq script
# Returns     : String (the valid DESeq script)
# Parameters  : String (the DESeq script)
# Throws      : If DESeq script is missing or not readable
# Comments    : None
sub _check_deseq_script {
    my ($deseq_script) = @_;
    return $deseq_script if defined $deseq_script && -r $deseq_script;
    confess 'No DESeq script specified' if !defined $deseq_script;
    confess "DESeq script ($deseq_script) does not exist or cannot be read";
}

=method output_sig_level

  Usage       : my $output_sig_level = $analysis->output_sig_level;
  Purpose     : Getter for output significance level attribute
  Returns     : +ve Float
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub output_sig_level {
    my ($self) = @_;
    return $output_sig_level{ id $self};
}

=method set_output_sig_level

  Usage       : $analysis->set_output_sig_level(0.001);
  Purpose     : Setter for output significance level attribute
  Returns     : undef
  Parameters  : +ve Float (the output significance level)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_output_sig_level {
    my ( $self, $arg ) = @_;
    $output_sig_level{ id $self} = _check_output_sig_level($arg);
    return;
}

# Usage       : $output_sig_level = _check_output_sig_level($output_sig_level);
# Purpose     : Check for valid output significance level
# Returns     : +ve Float (the valid output significance level)
# Parameters  : +ve Float (the output significance level)
# Throws      : If output significance level is missing or not a positive float
# Comments    : None
sub _check_output_sig_level {
    my ($output_sig_level) = @_;
    return $output_sig_level
      if defined $output_sig_level
      && $output_sig_level =~ m/\A \d* [.] \d+ \z/xms;
    confess 'No output significance level specified'
      if !defined $output_sig_level;
    confess "Invalid output significance level ($output_sig_level) specified";
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

=method ensembl_host

  Usage       : my $ensembl_host = $analysis->ensembl_host;
  Purpose     : Getter for Ensembl host attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

=method fasta_index

  Usage       : my $fai = $analysis->fasta_index;
  Purpose     : Getter for FASTA index attribute
  Returns     : Bio::DB::Sam::Fai
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub fasta_index {
    my ($self) = @_;

    if ( !defined $fasta_index{ id $self} && $self->ref_fasta ) {

        # We can create a FASTA index object
        $self->set_fasta_index( Bio::DB::Sam::Fai->load( $self->ref_fasta ) );
    }

    return $fasta_index{ id $self};
}

=method set_fasta_index

  Usage       : $analysis->set_fasta_index($fai);
  Purpose     : Setter for FASTA index attribute
  Returns     : undef
  Parameters  : Bio::DB::Sam::Fai
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
# Returns     : Bio::DB::Sam::Fai
# Parameters  : Bio::DB::Sam::Fai
# Throws      : If FASTA index is missing or invalid (i.e. not a
#               Bio::DB::Sam::Fai object)
# Comments    : None
sub _check_fasta_index {
    my ($fasta_index) = @_;
    return $fasta_index
      if defined $fasta_index && $fasta_index->isa('Bio::DB::Sam::Fai');
    confess 'No FASTA index specified' if !defined $fasta_index;
    confess 'Class of FASTA index (', ref $fasta_index,
      ') not Bio::DB::Sam::Fai';
}

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
    confess 'No Ensembl slice adaptor specified' if !defined $slice_adaptor;
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

    my @seqs = @{ $self->get_all_sequences() };

    # Get total sequence length
    my $total_bp = 0;
    foreach my $seq (@seqs) {
        $total_bp += $seq->bp;
    }

    # Get chunk target size (+ 1 to ensure slight overestimate)
    my $target_chunk_size = int( $total_bp / $self->chunk_total + 1 );

    my @chunks;
    my @chunk_size = map { 0 } 1 .. $self->chunk_total;

    # Iterate over sequences
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
        next if defined $chunks[$empty_chunk_index];    # Only want empty chunks

        # Find chunk with highest number of sequences (but more than one)
        my $max_seqs_chunk_index;
        my $max_seqs;
        foreach my $chunk_index ( 0 .. $self->chunk_total - 1 ) {
            next if !defined $chunks[$chunk_index]; # Only want non-empty chunks
            my $seqs = scalar @{ $chunks[$chunk_index] };
            if ( $seqs > 1 && ( !defined $max_seqs || $seqs > $max_seqs ) ) {
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
    if ( $self->test_chunk && exists $chunks->[ $self->test_chunk - 1 ] ) {
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

  Usage       : @tags = $analysis->list_all_tags_by_bam_file();
  Purpose     : Get all tags used in an analysis in a particular BAM file
  Returns     : Arrayref of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub list_all_tags_by_bam_file {
    my ( $self, $bam_file ) = @_;

    my $samples = $self->get_all_samples();

    my @tags = map { $_->tag } grep { $_->bam_file eq $bam_file } @{$samples};

    return uniq( sort @tags );
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
          $self->slice_adaptor->fetch_by_region( 'toplevel', $seq_name, $start,
            $end, $strand )->seq;
    }
    else {
        confess 'No reference FASTA or Ensembl database';
    }

    return uc $subseq;
}

# Usage       : $self->_create_slice_adaptor();
# Purpose     : Create an Ensembl slice adaptor
# Returns     : Undef
# Parameters  : None
# Throws      : No exceptions
# Comments    : None
sub _create_slice_adaptor {
    my ($self) = @_;

    my $host =
      $self->ensembl_host ? $self->ensembl_host : $DEFAULT_ENSEMBL_HOST;
    my $port = $self->ensembl_port;
    my $user =
      $self->ensembl_user ? $self->ensembl_user : $DEFAULT_ENSEMBL_USER;
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
          Bio::EnsEMBL::Registry->get_adaptor( $self->ensembl_species, 'core',
            'slice' );
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
