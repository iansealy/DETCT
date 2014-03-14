## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::SAMtools;
## use critic

# ABSTRACT: Miscellaneous functions wrapping SAMtools

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-21
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
  make_index
  flagstats
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func make_index

  Usage       : DETCT::Misc::SAMtools::make_index( {
                    dir             => '.',
                    bam_file        => $bam_file,
                    samtools_binary => 'samtools',
                } );
  Purpose     : Run SAMtools index
  Returns     : undef
  Parameters  : Hashref {
                    dir             => String (the working directory),
                    bam_file        => String (the BAM file),
                    samtools_binary => String (the SAMtools binary),
                }
  Throws      : If directory is missing
                If BAM file is missing
                If SAMtools binary is missing
                If command line can't be run
  Comments    : None

=cut

sub make_index {
    my ($arg_ref) = @_;

    confess 'No directory specified' if !defined $arg_ref->{dir};
    confess 'No BAM file specified'  if !defined $arg_ref->{bam_file};
    confess 'No SAMtools binary specified'
      if !defined $arg_ref->{samtools_binary};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'index.o' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'index.e' );

    my $cmd = join q{ }, $arg_ref->{samtools_binary}, 'index',
      $arg_ref->{bam_file};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

=func flagstats

  Usage       : DETCT::Misc::SAMtools::flagstats( {
                    dir             => '.',
                    bam_file        => $bam_file,
                    samtools_binary => 'samtools',
                } );
  Purpose     : Run SAMtools flagstat
  Returns     : undef
  Parameters  : Hashref {
                    dir             => String (the working directory),
                    bam_file        => String (the BAM file),
                    samtools_binary => String (the SAMtools binary),
                }
  Throws      : If directory is missing
                If BAM file is missing
                If SAMtools binary is missing
                If command line can't be run
  Comments    : None

=cut

sub flagstats {
    my ($arg_ref) = @_;

    confess 'No directory specified' if !defined $arg_ref->{dir};
    confess 'No BAM file specified'  if !defined $arg_ref->{bam_file};
    confess 'No SAMtools binary specified'
      if !defined $arg_ref->{samtools_binary};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'flagstat.txt' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'flagstat.e' );

    my $cmd = join q{ }, $arg_ref->{samtools_binary}, 'flagstat',
      $arg_ref->{bam_file};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

=func sort

  Usage       : DETCT::Misc::SAMtools::sort( {
                    dir             => '.',
                    input_bam_file  => $input_bam_file,
                    output_bam_file => $output_bam_file,
                    samtools_binary => 'samtools',
                } );
  Purpose     : Run SAMtools sort
  Returns     : undef
  Parameters  : Hashref {
                    dir             => String (the working directory),
                    input_bam_file  => String (the input BAM file),
                    output_bam_file => String (the output BAM file),
                    samtools_binary => String (the SAMtools binary),
                    sort_order      => String (the sort order),
                }
  Throws      : If directory is missing
                If input BAM file is missing
                If output BAM file is missing
                If SAMtools binary is missing
                If command line can't be run
                If sort order is invalid (not coordinate, queryname or undef)
  Comments    : Defaults to coordinate sorting

=cut

sub sort {
    my ($arg_ref) = @_;

    confess 'No directory specified' if !defined $arg_ref->{dir};
    confess 'No input BAM file specified'
      if !defined $arg_ref->{input_bam_file};
    confess 'No output BAM file specified'
      if !defined $arg_ref->{output_bam_file};
    confess 'No SAMtools binary specified'
      if !defined $arg_ref->{samtools_binary};

    my $sort_order =
      defined $arg_ref->{sort_order} ? $arg_ref->{sort_order} : 'coordinate';
    confess 'Invalid sort order specified'
      if $sort_order ne 'coordinate' && $sort_order ne 'queryname';

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'sort.o' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'sort.e' );

    my @options = ('-f');
    if ( $sort_order eq 'queryname' ) {
        push @options, '-n';
    }

    my $cmd = join q{ }, $arg_ref->{samtools_binary}, 'sort', @options,
      $arg_ref->{input_bam_file}, $arg_ref->{output_bam_file};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

1;
