## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::Picard;
## use critic

# ABSTRACT: Miscellaneous functions wrapping Picard

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-19
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use English qw( -no_match_vars );
use POSIX qw( WIFEXITED);
use File::Spec;
use File::Path qw( make_path );

use base qw( Exporter );
our @EXPORT_OK = qw(
  mark_duplicates
  merge
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func mark_duplicates

  Usage       : DETCT::Misc::Picard::mark_duplicates( {
                    dir                 => '.',
                    input_bam_file      => $input_bam_file,
                    output_bam_file     => $output_bam_file,
                    metrics_file        => $metrics_file,
                    java_binary         => 'java',
                    mark_duplicates_jar => 'MarkDuplicates.jar',
                    memory              => 4000,
                    consider_tags       => 1,
                } );
  Purpose     : Run MarkDuplicates
  Returns     : undef
  Parameters  : Hashref {
                    dir                 => String (the working directory),
                    input_bam_file      => String (the input BAM file),
                    output_bam_file     => String (the output BAM file),
                    metrics_file        => String (the metrics file),
                    java_binary         => String (the Java binary),
                    mark_duplicates_jar => String (the MarkDuplicates JAR),
                    memory              => Int (the memory allocated),
                    consider_tags       => Boolean (whether to consider tags),
                }
  Throws      : If directory is missing
                If input BAM file is missing
                If output BAM file is missing
                If metrics file is missing
                If Java binary is missing
                If MarkDuplicates JAR is missing
                If command line can't be run
  Comments    : None

=cut

sub mark_duplicates {
    my ($arg_ref) = @_;

    confess 'No directory specified' if !defined $arg_ref->{dir};
    confess 'No input BAM file specified'
      if !defined $arg_ref->{input_bam_file};
    confess 'No output BAM file specified'
      if !defined $arg_ref->{output_bam_file};
    confess 'No metrics file specified' if !defined $arg_ref->{metrics_file};
    confess 'No Java binary specified'  if !defined $arg_ref->{java_binary};
    confess 'No MarkDuplicates JAR specified'
      if !defined $arg_ref->{mark_duplicates_jar};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Options
    my %option = (
        INPUT                 => $arg_ref->{input_bam_file},
        OUTPUT                => $arg_ref->{output_bam_file},
        METRICS_FILE          => $arg_ref->{metrics_file},
        REMOVE_DUPLICATES     => 'false',
        TMP_DIR               => $arg_ref->{dir},
        VERBOSITY             => 'WARNING',
        QUIET                 => 'true',
        VALIDATION_STRINGENCY => 'SILENT',
        CREATE_INDEX          => 'false',
    );
    if ( $arg_ref->{consider_tags} ) {
        $option{CONSIDER_TAGS} = 'true';
    }
    my @options = map { $_ . q{=} . $option{$_} } sort keys %option;

    my $stdout_file =
      File::Spec->catfile( $arg_ref->{dir}, 'markduplicates.o' );
    my $stderr_file =
      File::Spec->catfile( $arg_ref->{dir}, 'markduplicates.e' );

    my $memory =
      $arg_ref->{memory}
      ? sprintf '-Xmx%dm', $arg_ref->{memory}
      : q{};

    my $cmd = join q{ }, $arg_ref->{java_binary}, '-XX:ParallelGCThreads=1',
      $memory, '-jar', $arg_ref->{mark_duplicates_jar}, @options;
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

=func merge

  Usage       : DETCT::Misc::Picard::merge( {
                    dir                 => '.',
                    input_bam_files     => \@input_bam_files,
                    output_bam_file     => $output_bam_file,
                    java_binary         => 'java',
                    merge_sam_files_jar => 'MergeSamFiles.jar',
                    memory              => 4000,
                } );
  Purpose     : Run MergeSamFiles
  Returns     : undef
  Parameters  : Hashref {
                    dir                 => String (the working directory),
                    input_bam_files     => Arrayref of strings (the BAM files),
                    output_bam_file     => String (the output BAM file),
                    java_binary         => String (the Java binary),
                    merge_sam_files_jar => String (the MergeSamFiles JAR),
                    memory              => Int (the memory allocated),
                }
  Throws      : If directory is missing
                If input BAM files are missing
                If output BAM file is missing
                If Java binary is missing
                If MergeSamFiles JAR is missing
                If command line can't be run
  Comments    : None

=cut

sub merge {
    my ($arg_ref) = @_;

    confess 'No directory specified' if !defined $arg_ref->{dir};
    confess 'No input BAM files specified'
      if !defined $arg_ref->{input_bam_files};
    confess 'No output BAM file specified'
      if !defined $arg_ref->{output_bam_file};
    confess 'No Java binary specified' if !defined $arg_ref->{java_binary};
    confess 'No MergeSamFiles JAR specified'
      if !defined $arg_ref->{merge_sam_files_jar};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Options
    my %option = (
        OUTPUT                      => $arg_ref->{output_bam_file},
        ASSUME_SORTED               => 'false',
        MERGE_SEQUENCE_DICTIONARIES => 'true',
        TMP_DIR                     => $arg_ref->{dir},
        VERBOSITY                   => 'WARNING',
        QUIET                       => 'true',
        VALIDATION_STRINGENCY       => 'SILENT',
        CREATE_INDEX                => 'false',
    );
    my @options = map { $_ . q{=} . $option{$_} } sort keys %option;
    foreach my $input_bam_file ( @{ $arg_ref->{input_bam_files} } ) {
        unshift @options, 'INPUT=' . $input_bam_file;
    }

    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'mergesamfiles.o' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'mergesamfiles.e' );

    my $memory =
      $arg_ref->{memory}
      ? sprintf '-Xmx%dm', $arg_ref->{memory}
      : q{};

    my $cmd = join q{ }, $arg_ref->{java_binary}, '-XX:ParallelGCThreads=1',
      $memory, '-jar', $arg_ref->{merge_sam_files_jar}, @options;
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

=func bam_to_fastq

  Usage       : DETCT::Misc::Picard::bam_to_fastq( {
                    dir              => '.',
                    bam_file         => $bam_file,
                    read1_fastq_file => $read1_fastq_file,
                    read2_fastq_file => $read2_fastq_file,
                    java_binary      => 'java',
                    bam_to_fastq_jar => 'SamToFastq.jar',
                    memory           => 2000,
                } );
  Purpose     : Run SamToFastq
  Returns     : undef
  Parameters  : Hashref {
                    dir              => String (the working directory),
                    bam_file         => String (the BAM file),
                    read1_fastq_file => String (the read 1 FASTQ file)
                    read2_fastq_file => String (the read 2 FASTQ file)
                    java_binary      => String (the Java binary),
                    bam_to_fastq_jar => String (the SamToFastq JAR),
                    memory           => Int (the memory allocated),
                }
  Throws      : If directory is missing
                If BAM file is missing
                If read 1 FASTQ file is missing
                If read 2 FASTQ file is missing
                If Java binary is missing
                If SamToFastq JAR is missing
                If command line can't be run
  Comments    : None

=cut

sub bam_to_fastq {
    my ($arg_ref) = @_;

    confess 'No directory specified' if !defined $arg_ref->{dir};
    confess 'No BAM file specified'  if !defined $arg_ref->{bam_file};
    confess 'No read 1 FASTQ file specified'
      if !defined $arg_ref->{read1_fastq_file};
    confess 'No read 2 FASTQ file specified'
      if !defined $arg_ref->{read2_fastq_file};
    confess 'No Java binary specified' if !defined $arg_ref->{java_binary};
    confess 'No SamToFastq JAR specified'
      if !defined $arg_ref->{bam_to_fastq_jar};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Options
    my %option = (
        INPUT                 => $arg_ref->{bam_file},
        FASTQ                 => $arg_ref->{read1_fastq_file},
        SECOND_END_FASTQ      => $arg_ref->{read2_fastq_file},
        TMP_DIR               => $arg_ref->{dir},
        VERBOSITY             => 'WARNING',
        QUIET                 => 'true',
        VALIDATION_STRINGENCY => 'SILENT',
        CREATE_INDEX          => 'false',
    );
    my @options = map { $_ . q{=} . $option{$_} } sort keys %option;

    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'bamtofastq.o' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'bamtofastq.e' );

    my $memory =
      $arg_ref->{memory}
      ? sprintf '-Xmx%dm', $arg_ref->{memory}
      : q{};

    my $cmd = join q{ }, $arg_ref->{java_binary}, '-XX:ParallelGCThreads=1',
      $memory, '-jar', $arg_ref->{bam_to_fastq_jar}, @options;
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

1;
