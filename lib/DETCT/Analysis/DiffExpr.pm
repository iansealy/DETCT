## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Analysis::DiffExpr;
## VERSION
## use critic

# ABSTRACT: Object representing differential expression analysis

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-14
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
private read1_length           => my %read1_length;           # e.g. 30
private read2_length           => my %read2_length;           # e.g. 54
private mismatch_threshold     => my %mismatch_threshold;     # e.g. 2
private mapq_threshold         => my %mapq_threshold;         # e.g. 10
private bin_size               => my %bin_size;               # e.g. 100
private peak_buffer_width      => my %peak_buffer_width;      # e.g. 100
private hmm_sig_level          => my %hmm_sig_level;          # e.g. 0.001
private hmm_binary             => my %hmm_binary;             # e.g. chiphmmnew
private exonerate_binary       => my %exonerate_binary;       # e.g. exonerate
private r_binary               => my %r_binary;               # e.g. R
private deseq_script           => my %deseq_script;           # e.g. run_deseq.R
private filter_percentile      => my %filter_percentile;      # e.g. 40
private spike_prefix           => my %spike_prefix;           # e.g. ERCC
private normalisation_method   => my %normalisation_method;   # e.g. spike
private deseq_model            => my %deseq_model;            # e.g. additive
private output_sig_level       => my %output_sig_level;       # e.g. 0.05
private table_file             => my %table_file;             # e.g. all.tsv
private table_format           => my %table_format;           # e.g. tsv
private control_condition      => my %control_condition;      # e.g. sibling
private experimental_condition => my %experimental_condition; # e.g. mutant
private skip_transcript => my %skip_transcript; # hashref of skipped transcripts
private ensembl_transcript_table_file =>
  my %ensembl_transcript_table_file;            # e.g. transcripts.tsv
private ensembl_transcript_table_format =>
  my %ensembl_transcript_table_format;          # e.g. tsv

# Constants
Readonly our $MAX_CONDITION_LENGTH => 128;

=method new

  Usage       : my $analysis = DETCT::Analysis::DiffExpr->new( {
                    name               => 'zmp_ph1',
                    read1_length       => 30,
                    read2_length       => 54,
                    mismatch_threshold => 2,
                    mapq_threshold     => 10,
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
  Returns     : DETCT::Analysis::DiffExpr
  Parameters  : Hashref {
                    name                            => String,
                    read1_length                    => Int,
                    read2_length                    => Int,
                    mismatch_threshold              => Int,
                    mapq_threshold                  => Int,
                    bin_size                        => Int,
                    peak_buffer_width               => Int,
                    hmm_sig_level                   => Float,
                    hmm_binary                      => String,
                    exonerate_binary                => String,
                    r_binary                        => String,
                    deseq_script                    => String,
                    filter_percentile               => Int,
                    spike_prefix                    => String,
                    normalisation_method            => String,
                    deseq_model                     => String,
                    output_sig_level                => Float,
                    table_file                      => String,
                    table_format                    => String,
                    control_condition               => String,
                    experimental_condition          => String,
                    ensembl_transcript_table_file   => String,
                    ensembl_transcript_table_format => String,
                    ref_fasta                       => String or undef,
                    ensembl_host                    => String or undef,
                    ensembl_port                    => Int or undef,
                    ensembl_user                    => String or undef,
                    ensembl_pass                    => String or undef,
                    ensembl_name                    => String or undef,
                    ensembl_species                 => String or undef,
                    chunk_total                     => Int,
                    test_chunk                      => Int or undef,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = $class->SUPER::new($arg_ref);
    $self->set_read1_length( $arg_ref->{read1_length} );
    $self->set_read2_length( $arg_ref->{read2_length} );
    $self->set_mismatch_threshold( $arg_ref->{mismatch_threshold} );
    $self->set_mapq_threshold( $arg_ref->{mapq_threshold} );
    $self->set_bin_size( $arg_ref->{bin_size} );
    $self->set_peak_buffer_width( $arg_ref->{peak_buffer_width} );
    $self->set_hmm_sig_level( $arg_ref->{hmm_sig_level} );
    $self->set_hmm_binary( $arg_ref->{hmm_binary} );
    $self->set_exonerate_binary( $arg_ref->{exonerate_binary} );
    $self->set_r_binary( $arg_ref->{r_binary} );
    $self->set_deseq_script( $arg_ref->{deseq_script} );
    $self->set_filter_percentile( $arg_ref->{filter_percentile} );
    $self->set_spike_prefix( $arg_ref->{spike_prefix} );
    $self->set_normalisation_method( $arg_ref->{normalisation_method} );
    $self->set_deseq_model( $arg_ref->{deseq_model} );
    $self->set_output_sig_level( $arg_ref->{output_sig_level} );
    $self->set_table_file( $arg_ref->{table_file} );
    $self->set_table_format( $arg_ref->{table_format} );
    $self->set_control_condition( $arg_ref->{control_condition} );
    $self->set_experimental_condition( $arg_ref->{experimental_condition} );
    $self->set_ensembl_transcript_table_file(
        $arg_ref->{ensembl_transcript_table_file} );
    $self->set_ensembl_transcript_table_format(
        $arg_ref->{ensembl_transcript_table_format} );
    return $self;
}

=method new_from_yaml

  Usage       : my $analysis
                    = DETCT::Analysis::DiffExpr->new_from_yaml('zmp_ph1.yaml');
  Purpose     : Constructor for creating analysis objects from a YAML file
  Returns     : DETCT::Analysis::DiffExpr
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

    $self->set_read1_length( $yaml->[0]->{read1_length} );
    $self->set_read2_length( $yaml->[0]->{read2_length} );
    $self->set_mismatch_threshold( $yaml->[0]->{mismatch_threshold} );
    $self->set_mapq_threshold( $yaml->[0]->{mapq_threshold} );
    $self->set_bin_size( $yaml->[0]->{bin_size} );
    $self->set_peak_buffer_width( $yaml->[0]->{peak_buffer_width} );
    $self->set_hmm_sig_level( $yaml->[0]->{hmm_sig_level} );
    $self->set_hmm_binary( $yaml->[0]->{hmm_binary} );
    $self->set_exonerate_binary( $yaml->[0]->{exonerate_binary} );
    $self->set_r_binary( $yaml->[0]->{r_binary} );
    $self->set_deseq_script( $yaml->[0]->{deseq_script} );
    $self->set_filter_percentile( $yaml->[0]->{filter_percentile} );
    $self->set_spike_prefix( $yaml->[0]->{spike_prefix} );
    $self->set_normalisation_method( $yaml->[0]->{normalisation_method} );
    $self->set_deseq_model( $yaml->[0]->{deseq_model} );
    $self->set_output_sig_level( $yaml->[0]->{output_sig_level} );
    $self->set_table_file( $yaml->[0]->{table_file} );
    $self->set_table_format( $yaml->[0]->{table_format} );
    $self->set_control_condition( $yaml->[0]->{control_condition} );
    $self->set_experimental_condition( $yaml->[0]->{experimental_condition} );
    $self->add_all_skip_transcripts( $yaml->[0]->{skip_transcripts} );
    $self->set_ensembl_transcript_table_file(
        $yaml->[0]->{ensembl_transcript_table_file} );
    $self->set_ensembl_transcript_table_format(
        $yaml->[0]->{ensembl_transcript_table_format} );

    return $self;
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

=method mapq_threshold

  Usage       : my $mapq_threshold = $analysis->mapq_threshold;
  Purpose     : Getter for MAPQ threshold attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : Defaults to 0

=cut

sub mapq_threshold {
    my ($self) = @_;
    return $mapq_threshold{ id $self} || 0;
}

=method set_mapq_threshold

  Usage       : $analysis->set_mapq_threshold(20);
  Purpose     : Setter for MAPQ threshold attribute
  Returns     : undef
  Parameters  : +ve Int (the MAPQ threshold)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_mapq_threshold {
    my ( $self, $arg ) = @_;
    $mapq_threshold{ id $self} = _check_mapq_threshold($arg);
    return;
}

# Usage       : $mapq_threshold = _check_mapq_threshold($mapq_threshold);
# Purpose     : Check for valid MAPQ threshold
# Returns     : +ve Int (the valid MAPQ threshold)
# Parameters  : +ve Int (the MAPQ threshold)
# Throws      : If MAPQ threshold is not a positive integer
# Comments    : None
sub _check_mapq_threshold {
    my ($mapq_threshold) = @_;
    return $mapq_threshold
      if !defined $mapq_threshold || $mapq_threshold =~ m/\A \d+ \z/xms;
    confess "Invalid MAPQ threshold ($mapq_threshold) specified";
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

=method exonerate_binary

  Usage       : my $exonerate_binary = $analysis->exonerate_binary;
  Purpose     : Getter for Exonerate binary attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub exonerate_binary {
    my ($self) = @_;
    return $exonerate_binary{ id $self};
}

=method set_exonerate_binary

  Usage       : $analysis->set_exonerate_binary('exonerate');
  Purpose     : Setter for Exonerate binary attribute
  Returns     : undef
  Parameters  : String (the Exonerate binary)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_exonerate_binary {
    my ( $self, $arg ) = @_;
    $exonerate_binary{ id $self} = $arg;
    return;
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

=method filter_percentile

  Usage       : my $filter_percentile = $analysis->filter_percentile;
  Purpose     : Getter for filter percentile attribute
  Returns     : +ve Int
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub filter_percentile {
    my ($self) = @_;
    return $filter_percentile{ id $self};
}

=method set_filter_percentile

  Usage       : $analysis->set_filter_percentile(40);
  Purpose     : Setter for filter percentile attribute
  Returns     : undef
  Parameters  : +ve Int (the filter percentile)
  Throws      : No exceptions
  Comments    : Defaults to 0 (i.e. no filtering)

=cut

sub set_filter_percentile {
    my ( $self, $arg ) = @_;
    $filter_percentile{ id $self} = _check_filter_percentile( $arg || 0 );
    return;
}

# Usage       : $filter_percentile
#                   = _check_filter_percentile($filter_percentile);
# Purpose     : Check for valid filter percentile
# Returns     : +ve Int (the valid filter percentile)
# Parameters  : +ve Int (the filter percentile)
# Throws      : If filter percentile is defined but not a positive integer
# Comments    : None
sub _check_filter_percentile {
    my ($filter_percentile) = @_;
    return $filter_percentile
      if !defined $filter_percentile || $filter_percentile =~ m/\A \d+ \z/xms;
    confess "Invalid filter percentile ($filter_percentile) specified";
}

=method spike_prefix

  Usage       : my $spike_prefix = $analysis->spike_prefix;
  Purpose     : Getter for spike prefix attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub spike_prefix {
    my ($self) = @_;
    return $spike_prefix{ id $self};
}

=method set_spike_prefix

  Usage       : $analysis->set_spike_prefix('ERCC');
  Purpose     : Setter for spike prefix attribute
  Returns     : undef
  Parameters  : String (the spike prefix)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_spike_prefix {
    my ( $self, $arg ) = @_;
    $spike_prefix{ id $self} = $arg;
    return;
}

=method normalisation_method

  Usage       : my $normalisation_method = $analysis->normalisation_method;
  Purpose     : Getter for normalisation method attribute
  Returns     : String ('deseq', 'spike' or 'none')
  Parameters  : None
  Throws      : No exceptions
  Comments    : Defaults to 'deseq'

=cut

sub normalisation_method {
    my ($self) = @_;
    return $normalisation_method{ id $self} || 'deseq';
}

=method set_normalisation_method

  Usage       : $analysis->set_normalisation_method('spike');
  Purpose     : Setter for normalisation method attribute
  Returns     : undef
  Parameters  : String (the normalisation method)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_normalisation_method {
    my ( $self, $arg ) = @_;
    $normalisation_method{ id $self} = _check_normalisation_method($arg);
    return;
}

# Usage       : $normalisation_method = _check_normalisation_method($normalisation_method);
# Purpose     : Check for valid normalisation method
# Returns     : String (the valid normalisation method)
# Parameters  : String (the normalisation method)
# Throws      : If normalisation method is defined but invalid
# Comments    : None
sub _check_normalisation_method {
    my ($normalisation_method) = @_;
    return $normalisation_method
      if !defined $normalisation_method
      || any { $_ eq $normalisation_method } qw(deseq spike none);
    confess "Invalid normalisation method ($normalisation_method) specified";
}

=method deseq_model

  Usage       : my $deseq_model = $analysis->deseq_model;
  Purpose     : Getter for DESeq model attribute
  Returns     : String ('additive' or 'interaction')
  Parameters  : None
  Throws      : No exceptions
  Comments    : Defaults to 'additive'

=cut

sub deseq_model {
    my ($self) = @_;
    return $deseq_model{ id $self} || 'additive';
}

=method set_deseq_model

  Usage       : $analysis->set_deseq_model('interaction');
  Purpose     : Setter for DESeq model attribute
  Returns     : undef
  Parameters  : String (the DESeq model)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_deseq_model {
    my ( $self, $arg ) = @_;
    $deseq_model{ id $self} = _check_deseq_model($arg);
    return;
}

# Usage       : $deseq_model = _check_deseq_model($deseq_model);
# Purpose     : Check for valid DESeq model
# Returns     : String (the valid DESeq model)
# Parameters  : String (the DESeq model)
# Throws      : If DESeq model is defined but invalid
# Comments    : None
sub _check_deseq_model {
    my ($deseq_model) = @_;
    return $deseq_model
      if !defined $deseq_model
      || any { $_ eq $deseq_model } qw(additive interaction);
    confess "Invalid DESeq model ($deseq_model) specified";
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

=method table_file

  Usage       : my $table_file = $analysis->table_file;
  Purpose     : Getter for table file attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub table_file {
    my ($self) = @_;
    return $table_file{ id $self};
}

=method set_table_file

  Usage       : $analysis->set_table_file('all.tsv');
  Purpose     : Setter for table file attribute
  Returns     : undef
  Parameters  : String (the table file)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_table_file {
    my ( $self, $arg ) = @_;
    $table_file{ id $self} = _check_table_file($arg);
    return;
}

# Usage       : $table_file = _check_table_file($table_file);
# Purpose     : Check for valid table file
# Returns     : String (the valid table file)
# Parameters  : String (the table file)
# Throws      : If table file is defined but not readable
# Comments    : None
sub _check_table_file {
    my ($table_file) = @_;
    return $table_file if !defined $table_file || -r $table_file;
    confess "Table file ($table_file) cannot be read";
}

=method table_format

  Usage       : my $table_format = $analysis->table_format;
  Purpose     : Getter for table format attribute
  Returns     : String ('csv' or 'tsv')
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub table_format {
    my ($self) = @_;

    if ( !$table_format{ id $self} && $table_file{ id $self} ) {

        # Attempt to guess format from filename
        my ($extension) =
          $table_file{ id $self} =~ m/[.] ([[:lower:]]{3}) \z/xms;

        try {
            $self->set_table_format($extension);
        };
    }

    return $table_format{ id $self};
}

=method set_table_format

  Usage       : $analysis->set_table_format('tsv');
  Purpose     : Setter for table format attribute
  Returns     : undef
  Parameters  : String (the table format)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_table_format {
    my ( $self, $arg ) = @_;
    $table_format{ id $self} = _check_table_format($arg);
    return;
}

# Usage       : $table_format = _check_table_format($table_format);
# Purpose     : Check for valid table format
# Returns     : String (the valid table format)
# Parameters  : String (the table format)
# Throws      : If table format is defined but invalid
# Comments    : None
sub _check_table_format {
    my ($table_format) = @_;
    return $table_format
      if !defined $table_format || any { $_ eq $table_format } qw(csv tsv);
    confess "Invalid table format ($table_format) specified";
}

=method control_condition

  Usage       : my $control_condition = $analysis->control_condition;
  Purpose     : Getter for control condition attribute
  Returns     : String (e.g. "sibling")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub control_condition {
    my ($self) = @_;
    return $control_condition{ id $self};
}

=method set_control_condition

  Usage       : $analysis->set_control_condition('sibling');
  Purpose     : Setter for control condition attribute
  Returns     : undef
  Parameters  : String (the control condition)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_control_condition {
    my ( $self, $arg ) = @_;
    $control_condition{ id $self} = _check_control_condition($arg);
    return;
}

# Usage       : $control_condition = _check_control_condition($control_condition);
# Purpose     : Check for valid control condition
# Returns     : String (the valid control condition)
# Parameters  : String (the control condition)
# Throws      : If control condition is empty
#               If control condition > $MAX_CONDITION_LENGTH characters
# Comments    : None
sub _check_control_condition {
    my ($control_condition) = @_;

    confess 'Empty control condition specified'
      if defined $control_condition && !length $control_condition;
    confess( sprintf 'Control condition (%s) longer than %d characters',
        $control_condition, $MAX_CONDITION_LENGTH )
      if defined $control_condition
      && length $control_condition > $MAX_CONDITION_LENGTH;

    return $control_condition;
}

=method experimental_condition

  Usage       : my $experimental_condition = $analysis->experimental_condition;
  Purpose     : Getter for experimental condition attribute
  Returns     : String (e.g. "mutant")
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub experimental_condition {
    my ($self) = @_;
    return $experimental_condition{ id $self};
}

=method set_experimental_condition

  Usage       : $analysis->set_experimental_condition('mutant');
  Purpose     : Setter for experimental condition attribute
  Returns     : undef
  Parameters  : String (the experimental condition)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_experimental_condition {
    my ( $self, $arg ) = @_;
    $experimental_condition{ id $self} = _check_experimental_condition($arg);
    return;
}

# Usage       : $experimental_condition
#                   = _check_experimental_condition($experimental_condition);
# Purpose     : Check for valid experimental condition
# Returns     : String (the valid experimental condition)
# Parameters  : String (the experimental condition)
# Throws      : If experimental condition is empty
#               If experimental condition > $MAX_CONDITION_LENGTH characters
# Comments    : None
sub _check_experimental_condition {
    my ($experimental_condition) = @_;

    confess 'Empty experimental condition specified'
      if defined $experimental_condition && !length $experimental_condition;
    confess( sprintf 'Control condition (%s) longer than %d characters',
        $experimental_condition, $MAX_CONDITION_LENGTH )
      if defined $experimental_condition
      && length $experimental_condition > $MAX_CONDITION_LENGTH;

    return $experimental_condition;
}

=method add_all_skip_transcripts

  Usage       : $analysis->add_all_skip_transcripts(['ENSDART00000135768']);
  Purpose     : Add all skip transcripts to an analysis
  Returns     : undef
  Parameters  : Arrayref of strings (the skip transcripts) or undef
  Throws      : No exceptions
  Comments    : None

=cut

sub add_all_skip_transcripts {
    my ( $self, $skip_transcripts ) = @_;

    $skip_transcript{ id $self} = {};

    foreach my $id ( @{ $skip_transcripts || [] } ) {
        $skip_transcript{ id $self}->{$id} = 1;
    }

    return;
}

=method get_all_skip_transcripts

  Usage       : $skip_transcripts = $analysis->get_all_skip_transcripts();
  Purpose     : Get all skip transcripts of an analysis
  Returns     : Arrayref of strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub get_all_skip_transcripts {
    my ($self) = @_;

    return [ sort keys %{ $skip_transcript{ id $self} || {} } ];
}

=method ensembl_transcript_table_file

  Usage       : my $ens_trans_file = $analysis->ensembl_transcript_table_file;
  Purpose     : Getter for Ensembl transcript table file attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ensembl_transcript_table_file {
    my ($self) = @_;
    return $ensembl_transcript_table_file{ id $self};
}

=method set_ensembl_transcript_table_file

  Usage       : $analysis->set_ensembl_transcript_table_file('transcripts.tsv');
  Purpose     : Setter for Ensembl transcript table file attribute
  Returns     : undef
  Parameters  : String (the Ensembl transcript table file)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ensembl_transcript_table_file {
    my ( $self, $arg ) = @_;
    $ensembl_transcript_table_file{ id $self} =
      _check_ensembl_transcript_table_file($arg);
    return;
}

# Usage       : $ens_trans_file
#                   = _check_ensembl_transcript_table_file($ens_trans_file);
# Purpose     : Check for valid Ensembl transcript table file
# Returns     : String (the valid Ensembl transcript table file)
# Parameters  : String (the Ensembl transcript table file)
# Throws      : If Ensembl transcript table file is defined but not readable
# Comments    : None
sub _check_ensembl_transcript_table_file {
    my ($ensembl_transcript_table_file) = @_;
    return $ensembl_transcript_table_file
      if !defined $ensembl_transcript_table_file
      || -r $ensembl_transcript_table_file;
    confess sprintf 'Ensembl transcript table file (%s) cannot be read',
      $ensembl_transcript_table_file;
}

=method ensembl_transcript_table_format

  Usage       : my $ens_trans_format = $analysis->ensembl_transcript_table_format;
  Purpose     : Getter for Ensembl transcript table format attribute
  Returns     : String ('csv' or 'tsv')
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub ensembl_transcript_table_format {
    my ($self) = @_;

    if (  !$ensembl_transcript_table_format{ id $self}
        && $ensembl_transcript_table_file{ id $self} )
    {

        # Attempt to guess format from filename
        my ($extension) = $ensembl_transcript_table_file{ id $self} =~
          m/[.] ([[:lower:]]{3}) \z/xms;

        try {
            $self->set_ensembl_transcript_table_format($extension);
        };
    }

    return $ensembl_transcript_table_format{ id $self};
}

=method set_ensembl_transcript_table_format

  Usage       : $analysis->set_ensembl_transcript_table_format('tsv');
  Purpose     : Setter for Ensembl transcript table format attribute
  Returns     : undef
  Parameters  : String (the Ensembl transcript table format)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_ensembl_transcript_table_format {
    my ( $self, $arg ) = @_;
    $ensembl_transcript_table_format{ id $self} =
      _check_ensembl_transcript_table_format($arg);
    return;
}

# Usage       : $ens_trans_format
#                   = _check_ensembl_transcript_table_format($ens_trans_format);
# Purpose     : Check for valid Ensembl transcript table format
# Returns     : String (the valid Ensembl transcript table format)
# Parameters  : String (the Ensembl transcript table format)
# Throws      : If Ensembl transcript table format is defined but invalid
# Comments    : None
sub _check_ensembl_transcript_table_format {
    my ($ensembl_transcript_table_format) = @_;
    return $ensembl_transcript_table_format
      if !defined $ensembl_transcript_table_format
      || any { $_ eq $ensembl_transcript_table_format } qw(csv tsv);
    confess sprintf 'Invalid Ensembl transcript table format (%s) specified',
      $ensembl_transcript_table_format;
}

=method validate

  Usage       : $analysis->validate();
  Purpose     : Check analysis
  Returns     : 1
  Parameters  : None
  Throws      : If spike prefix not present in BAM files
                If both control and experimental conditions not specified
                If control or experimental conditions not specified for a sample
  Comments    : None

=cut

sub validate {
    my ($self) = @_;

    $self->SUPER::validate();

    if ( $self->spike_prefix ) {
        my $spike      = $self->spike_prefix;
        my $seen_spike = 0;
        foreach my $seq ( @{ $self->get_all_sequences() } ) {
            if ( $seq->name =~ m/^$spike/xms ) {
                $seen_spike = 1;
            }
        }
        if ( !$seen_spike ) {
            confess "Spike prefix ($spike) not seen in BAM files";
        }
    }

    if ( $self->control_condition && !$self->experimental_condition ) {
        confess 'Control condition specified, but no experimental condition';
    }
    elsif ( $self->experimental_condition && !$self->control_condition ) {
        confess 'Control condition specified, but no experimental condition';
    }
    elsif ( $self->control_condition && $self->experimental_condition ) {
        my %is_condition = map { $_ => 1 } $self->list_all_conditions();
        if ( !exists $is_condition{ $self->control_condition } ) {
            confess 'Control condition (', $self->control_condition,
              ') not specified for any sample';
        }
        if ( !exists $is_condition{ $self->experimental_condition } ) {
            confess 'Experimental condition (', $self->experimental_condition,
              ') not specified for any sample';
        }
    }

    return 1;
}

1;
