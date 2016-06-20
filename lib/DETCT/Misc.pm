## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc;
## VERSION
## use critic

# ABSTRACT: Miscellaneous functions

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-07
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

use base qw( Exporter );
our @EXPORT_OK = qw(
  write_or_die
  print_or_die
  printf_or_die
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func write_or_die

  Usage       : write_or_die( $fh, @fastq );
  Purpose     : Print to a filehandle or confess
  Returns     : undef
  Parameters  : Filehandle
                Array
  Throws      : If printing fails
  Comments    : None

=cut

sub write_or_die {
    my ( $fh, @data ) = @_;

    print {$fh} @data or confess "Can't write: $OS_ERROR";

    return;
}

=func print_or_die

  Usage       : print_or_die( $fasta );
  Purpose     : Print or confess
  Returns     : undef
  Parameters  : Array
  Throws      : If printing fails
  Comments    : None

=cut

sub print_or_die {
    my (@data) = @_;

    print @data or confess "Can't print: $OS_ERROR";

    return;
}

=func printf_or_die

  Usage       : printf_or_die( ">%s\n", $seq_region );
  Purpose     : Print formatted string or confess
  Returns     : undef
  Parameters  : String (the format)
                Array
  Throws      : If printing fails
  Comments    : None

=cut

sub printf_or_die {
    my ( $format, @data ) = @_;

    printf $format, @data or confess "Can't printf: $OS_ERROR";

    return;
}

1;
