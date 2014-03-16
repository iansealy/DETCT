## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::BWA;
## use critic

# ABSTRACT: Miscellaneous functions wrapping BWA

## Author         : is1
## Maintainer     : is1
## Created        : 2014-03-16
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
  aln
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func aln

  Usage       : DETCT::Misc::BWA::aln( {
                    dir        => '.',
                    sai_file   => 'out.sai',
                    ref_fasta  => 'zv9.fa',
                    fastq_file => 'in.fastq',
                    bwa_binary => 'bwa',
                } );
  Purpose     : Run BWA aln
  Returns     : undef
  Parameters  : Hashref {
                    dir        => String (the working directory),
                    sai_file   => String (the sai file),
                    ref_fasta  => String (the reference FASTA file),
                    fastq_file => String (the FASTQ file)
                    bwa_binary => String (the BWA binary),
                }
  Throws      : If directory is missing
                If sai file is missing
                If reference FASTA is missing
                If FASTQ file is missing
                If BWA binary is missing
                If command line can't be run
  Comments    : None

=cut

sub aln {
    my ($arg_ref) = @_;

    confess 'No directory specified'       if !defined $arg_ref->{dir};
    confess 'No sai file specified'        if !defined $arg_ref->{sai_file};
    confess 'No reference FASTA specified' if !defined $arg_ref->{ref_fasta};
    confess 'No FASTQ file specified'      if !defined $arg_ref->{fastq_file};
    confess 'No BWA binary specified'      if !defined $arg_ref->{bwa_binary};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'aln.o' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'aln.e' );

    my $cmd = join q{ }, $arg_ref->{bwa_binary}, 'aln', '-f',
      $arg_ref->{sai_file}, $arg_ref->{ref_fasta}, $arg_ref->{fastq_file};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

1;
