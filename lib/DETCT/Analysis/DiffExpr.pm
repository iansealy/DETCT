## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Analysis::DiffExpr;
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
use YAML::Tiny;

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private read1_length       => my %read1_length;          # e.g. 30
private read2_length       => my %read2_length;          # e.g. 54
private mismatch_threshold => my %mismatch_threshold;    # e.g. 2
private bin_size           => my %bin_size;              # e.g. 100
private peak_buffer_width  => my %peak_buffer_width;     # e.g. 100
private hmm_sig_level      => my %hmm_sig_level;         # e.g. 0.001
private hmm_binary         => my %hmm_binary;            # e.g. chiphmmnew
private r_binary           => my %r_binary;              # e.g. R
private deseq_script       => my %deseq_script;          # e.g. ~/run_deseq.R
private output_sig_level   => my %output_sig_level;      # e.g. 0.05
private input_tsv_file     => my %input_tsv_file;        # e.g. all.tsv

=method new

  Usage       : my $analysis = DETCT::Analysis::DiffExpr->new( {
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
  Returns     : DETCT::Analysis::DiffExpr
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
                    input_tsv_file     => String,
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
    my $self = $class->SUPER::new($arg_ref);
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
    $self->set_input_tsv_file( $arg_ref->{input_tsv_file} );
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
    $self->set_bin_size( $yaml->[0]->{bin_size} );
    $self->set_peak_buffer_width( $yaml->[0]->{peak_buffer_width} );
    $self->set_hmm_sig_level( $yaml->[0]->{hmm_sig_level} );
    $self->set_hmm_binary( $yaml->[0]->{hmm_binary} );
    $self->set_r_binary( $yaml->[0]->{r_binary} );
    $self->set_deseq_script( $yaml->[0]->{deseq_script} );
    $self->set_output_sig_level( $yaml->[0]->{output_sig_level} );
    $self->set_input_tsv_file( $yaml->[0]->{input_tsv_file} );

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

=method input_tsv_file

  Usage       : my $input_tsv_file = $analysis->input_tsv_file;
  Purpose     : Getter for input TSV file attribute
  Returns     : String
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub input_tsv_file {
    my ($self) = @_;
    return $input_tsv_file{ id $self};
}

=method set_input_tsv_file

  Usage       : $analysis->set_input_tsv_file('all.tsv');
  Purpose     : Setter for input TSV file attribute
  Returns     : undef
  Parameters  : String (the input TSV file)
  Throws      : No exceptions
  Comments    : None

=cut

sub set_input_tsv_file {
    my ( $self, $arg ) = @_;
    $input_tsv_file{ id $self} = _check_input_tsv_file($arg);
    return;
}

# Usage       : $input_tsv_file = _check_input_tsv_file($input_tsv_file);
# Purpose     : Check for valid input TSV file
# Returns     : String (the valid input TSV file)
# Parameters  : String (the input TSV file)
# Throws      : If input TSV file is defined but not readable
# Comments    : None
sub _check_input_tsv_file {
    my ($input_tsv_file) = @_;
    return $input_tsv_file if !defined $input_tsv_file || -r $input_tsv_file;
    confess "Input TSV file ($input_tsv_file) cannot be read";
}

1;
