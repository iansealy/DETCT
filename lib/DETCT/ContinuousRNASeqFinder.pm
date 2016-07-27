## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::ContinuousRNASeqFinder;
## VERSION
## use critic

# ABSTRACT: Object for finding continuous RNA-Seq by location

## Author         : is1
## Maintainer     : is1
## Created        : 2016-07-22
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Class::InsideOut qw( private register id );
use List::MoreUtils qw( any );
use DETCT::Misc::Output qw(
  parse_ensembl_transcripts_table
);

=head1 SYNOPSIS

    # Brief code examples

=cut

# Attributes:
private cache        => my %cache;           # hashref
private table_file   => my %table_file;      # e.g. transcripts.tsv
private table_format => my %table_format;    # e.g. tsv

=method new

  Usage       : my $rnaseq_finder = DETCT::ContinuousRNASeqFinder->new( {
                    table_file   => $table_file,
                    table_format => $table_format,
                } );
  Purpose     : Constructor for continuous RNA-Seq finder objects
  Returns     : DETCT::ContinuousRNASeqFinder
  Parameters  : Hashref {
                    table_file   => String (the table file),
                    table_format => String (the table format),
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub new {
    my ( $class, $arg_ref ) = @_;
    my $self = register($class);
    $self->set_table_file( $arg_ref->{table_file} );
    $self->set_table_format( $arg_ref->{table_format} );
    return $self;
}

=method table_file

  Usage       : my $table_file = $rnaseq_finder->table_file;
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

  Usage       : $rnaseq_finder->set_table_file('transcripts.tsv');
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

  Usage       : my $table_format = $rnaseq_finder->table_format;
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

  Usage       : $rnaseq_finder->set_table_format('tsv');
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

=method get_containing_continuous_rnaseq

  Usage       : $finder->get_containing_continuous_rnaseq($seq_name, $pos,
                    $strand);
  Purpose     : Retrieve continuous RNA-Seq containing a 3' end
  Returns     : Arrayref of transcript stable ids
  Parameters  : String (the 3' end sequence name)
                Int (the 3' end position)
                Int (the 3' end strand)
  Throws      : No exceptions
  Comments    : None

=cut

sub get_containing_continuous_rnaseq {
    my ( $self, $seq_name, $pos, $strand ) = @_;

    # Ensure cache is filled
    $self->_fill_cache_from_file();

    return () if !defined $cache{ id $self}->{$seq_name}->{$strand};

    my @rnaseq = @{ $cache{ id $self}->{$seq_name}->{$strand} };
    if ( $strand > 0 ) {
        @rnaseq = grep { $pos >= $_->[1] && $pos <= $_->[2] } @rnaseq;
    }
    else {
        @rnaseq = grep { $pos >= $_->[2] && $pos <= $_->[1] } @rnaseq;
    }
    my @stable_ids = sort map { $_->[0] } @rnaseq;

    return @stable_ids;
}

# Usage       : $self->_fill_cache_from_file();
# Purpose     : Fill the cache from a file
# Returns     : undef
# Parameters  : None
# Throws      : No exceptions
# Comments    : Cache is a hashref (keyed by sequence name) of hashrefs (keyed
#               by strand) of arrayrefs of arrayrefs of continuous RNA-Seq
#               transcript IDs CDS end position and RNA-Seq end position

sub _fill_cache_from_file {
    my ($self) = @_;

    # Skip if cache already filled
    return if exists $cache{ id $self};

    my $transcripts = parse_ensembl_transcripts_table(
        {
            table_file   => $self->table_file,
            table_format => $self->table_format,
        }
    );

    foreach my $transcript ( @{$transcripts} ) {
        my ( $seq_name, $stable_id, $strand, $cds_end, $rnaseq_end ) =
          @{$transcript};
        push @{ $cache{ id $self}->{$seq_name}->{$strand} },
          [ $stable_id, $cds_end, $rnaseq_end ];
    }

    return;
}

=method add_continuous_rnaseq

  Usage       : my $regions_ref
                    = $finder->add_continuous_rnaseq($regions_ary_ref);
  Purpose     : Add continuous RNA-Seq to 3' ends
  Returns     : Arrayref [
                    Arrayref [
                        Int (region start),
                        Int (region end),
                        Int (region maximum read count),
                        Float (region log probability sum),
                        Int ( 1 or -1 ) (3' end strand)
                        Arrayref [
                            Arrayref [
                                String (3' end sequence name),
                                Int (3' end position),
                                Int (3' end strand),
                                Int (3' end read count),
                                Boolean (whether polyA),
                                String (upstream 14 bp),
                                String (downstream 14 bp),
                                Int (distance hexamer upstream) or undef,
                                String (hexamer sequence),
                                Int (distance to nearest transposon),
                                Int (position of nearest transposon),
                                Arrayref [
                                    String (continuous RNA-Seq transcript id),
                                    ... (continuous RNA-Seq transcript ids)
                                ]
                            ],
                            ... (3' ends)
                        ],
                    ],
                    ... (regions)
                ]
  Parameters  : Arrayref (of regions)
  Throws      : If regions are missing
  Comments    : None

=cut

sub add_continuous_rnaseq {
    my ( $self, $regions ) = @_;

    confess 'No regions specified' if !defined $regions;

    my @output;

    foreach my $region ( @{$regions} ) {
        my ( undef, undef, undef, undef, undef, $three_prime_ends ) =
          @{$region};

        foreach my $three_prime_end ( @{$three_prime_ends} ) {
            my ( $seq_name, $pos, $strand ) = @{$three_prime_end};
            my @transcripts =
              $self->get_containing_continuous_rnaseq( $seq_name, $pos,
                $strand );
            push @{$three_prime_end}, \@transcripts;
        }

        push @output, $region;
    }

    return \@output;
}

1;
