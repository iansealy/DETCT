## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::Tag;
## use critic

# ABSTRACT: Miscellaneous functions for interacting with DETCT read tags

## Author         : is1
## Maintainer     : is1
## Created        : 2013-01-07
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use base qw( Exporter );
our @EXPORT_OK = qw(
  detag_trim_fastq
  convert_tag_to_regexp
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func detag_trim_fastq

  Usage       : DETCT::Misc::Tag::detag_trim_fastq( {
                    fastq_read1_input     => $fastq_read1_input,
                    fastq_read2_input     => $fastq_read2_input,
                    fastq_output_prefix   => $fastq_output_prefix,
                    read1_required_length => $read1_required_length,
                    read2_required_length => $read2_required_length,
                    polyt_trim_length     => $polyt_trim_length,
                    polyt_min_length      => $polyt_min_length,
                    read_tags             => \@read_tags,
                } );
  Purpose     : Detag and trim FASTQ files
  Returns     : undef
  Parameters  : Hashref {
                    fastq_read1_input     => String (read 1 FASTQ file),
                    fastq_read2_input     => String (read 2 FASTQ file),
                    fastq_output_prefix   => String (prefix for output FASTQs),
                    read1_required_length => Int (required length of read 1),
                    read2_required_length => Int (required length of read 2),
                    polyt_trim_length     => Int (polyT length to be trimmed),
                    polyt_min_length      => Int (min Ts to define polyT),
                    read_tags             => Arrayref (of read tags),
                    no_pair_suffix        => Boolean or undef,
                    treat_n_in_polyt_as_t => Boolean or undef,
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub detag_trim_fastq {
    my ($arg_ref) = @_;

    # Assume all tags are same length
    my $tag_length = length $arg_ref->{read_tags}[0];

    my $min_polyt = q{T} x $arg_ref->{polyt_min_length};
    my $polyt_re = qr/$min_polyt/xms;    # Regexp for polyT matching

    my $read1_required_length = $arg_ref->{read1_required_length};
    my $read2_required_length = $arg_ref->{read2_required_length};
    my $polyt_trim_length     = $arg_ref->{polyt_trim_length};

    # Convert tags to regular expressions
    my @read_tags  = @{ $arg_ref->{read_tags} };
    my %re_tag_for = convert_tag_to_regexp(@read_tags);

    ## no critic (RequireBriefOpen)
    open my $fh1_in, '<', $arg_ref->{fastq_read1_input};
    open my $fh2_in, '<', $arg_ref->{fastq_read2_input};
    ## use critic
    my $fh_out_for = _open_output_fhs( $arg_ref->{fastq_output_prefix},
        $tag_length, @read_tags );

    # Check last tag seen first because likely to be seen next
    my $prev_tag = $read_tags[0];    # Arbitrarily choose first tag

    while ( my $read1_id = <$fh1_in> ) {
        my $read2_id   = <$fh2_in>;
        my $read1_seq  = <$fh1_in>;
        my $read2_seq  = <$fh2_in>;
        my $read1_plus = <$fh1_in>;
        my $read2_plus = <$fh2_in>;
        my $read1_qual = <$fh1_in>;
        my $read2_qual = <$fh2_in>;

        chomp $read1_id;
        chomp $read2_id;
        chomp $read1_seq;
        chomp $read2_seq;
        chomp $read1_plus;
        chomp $read2_plus;
        chomp $read1_qual;
        chomp $read2_qual;

        # Do we need to add pair suffix to read IDs?
        if ( $arg_ref->{no_pair_suffix} ) {
            $read1_id .= '/1';
            $read2_id .= '/2';
        }

        # Remove /1 or /2 from read ids and then check they match
        my $read1_id_no_suffix = $read1_id;
        my $read2_id_no_suffix = $read2_id;
        ## no critic (ProhibitMagicNumbers)
        substr $read1_id_no_suffix, -2, 2, q{};
        substr $read2_id_no_suffix, -2, 2, q{};
        ## use critic
        if ( $read1_id_no_suffix ne $read2_id_no_suffix ) {
            confess 'Read order does not match in input '
              . "($read1_id_no_suffix does not match $read2_id_no_suffix)";
        }

        # Trim reads to specified length if necessary
        my $read1_pre_detag_length =
          $read1_required_length + $tag_length + $polyt_trim_length;
        if ( length $read1_seq > $read1_pre_detag_length ) {
            $read1_seq  = substr $read1_seq,  0, $read1_pre_detag_length;
            $read1_qual = substr $read1_qual, 0, $read1_pre_detag_length;
        }
        if ( length $read2_seq > $read2_required_length ) {
            $read2_seq  = substr $read2_seq,  0, $read2_required_length;
            $read2_qual = substr $read2_qual, 0, $read2_required_length;
        }

        # Get tag and putative polyT from read 1
        my $tag_in_read = substr $read1_seq, 0, $tag_length;
        my $polyt_seq = substr $read1_seq, $tag_length, $polyt_trim_length;
        if ( $arg_ref->{treat_n_in_polyt_as_t} ) {
            $polyt_seq =~ s/N/T/xmsg;
        }

        # Default tag to add to id if no match
        my $tag_for_id = q{X} x $tag_length;
        my $tag_found  = q{X} x $tag_length;

        # Make sure a tag matches and polyT is present
      TAG: foreach my $tag ( $prev_tag, sort keys %re_tag_for ) {
            my $regexps = $re_tag_for{$tag};
            foreach my $re ( @{$regexps} ) {
                if ( $tag_in_read =~ $re && $polyt_seq =~ $polyt_re ) {
                    $tag_for_id = $tag_in_read;
                    $tag_found  = $tag;
                    substr $read1_seq, 0, $tag_length + $polyt_trim_length, q{};
                    substr $read1_qual, 0, $tag_length + $polyt_trim_length,
                      q{};
                    $prev_tag = $tag;
                    last TAG;    # Skip rest if got a match
                }
            }
        }

        # Add tag to id
        $read1_id =~ s{ /1 \z}{#$tag_for_id/1}xms;
        $read2_id =~ s{ /2 \z}{#$tag_for_id/2}xms;

        print { $fh_out_for->{$tag_found}->{1} } $read1_id,   "\n";
        print { $fh_out_for->{$tag_found}->{1} } $read1_seq,  "\n";
        print { $fh_out_for->{$tag_found}->{1} } $read1_plus, "\n";
        print { $fh_out_for->{$tag_found}->{1} } $read1_qual, "\n";
        print { $fh_out_for->{$tag_found}->{2} } $read2_id,   "\n";
        print { $fh_out_for->{$tag_found}->{2} } $read2_seq,  "\n";
        print { $fh_out_for->{$tag_found}->{2} } $read2_plus, "\n";
        print { $fh_out_for->{$tag_found}->{2} } $read2_qual, "\n";
    }

    close $fh1_in;
    close $fh2_in;
    _close_output_fhs($fh_out_for);

    return;
}

# Usage       : my $fh_out_for = _open_output_fhs(
#                   $fastq_output_prefix, $tag_length, @read_tags
#               );
# Purpose     : Open filehandles for all output FASTQ files
# Returns     : Hashref of hashref of filehandles
# Parameters  : String (prefix for output FASTQs)
#               Int (tag length)
#               Array of strings (the tags)
# Throws      : No exceptions
# Comments    : None
sub _open_output_fhs {
    my ( $fastq_output_prefix, $tag_length, @tags ) = @_;

    push @tags, q{X} x $tag_length;    # Default tag if no match

    my %fh_for;
    foreach my $tag (@tags) {
        foreach my $read ( 1, 2 ) {
            my $filename = join q{_}, $fastq_output_prefix, $tag, $read;
            $filename .= '.fastq';
            open my $fh, '>', $filename;    ## no critic (RequireBriefOpen)
            $fh_for{$tag}->{$read} = $fh;
        }
    }

    return \%fh_for;
}

# Usage       : _close_output_fhs($fh_out_for);
# Purpose     : Close filehandles for all output FASTQ files
# Returns     : undef
# Parameters  : Hashref of hashref of filehandles
# Throws      : No exceptions
# Comments    : None
sub _close_output_fhs {
    my ($fh_for) = @_;

    foreach my $tag ( keys %{$fh_for} ) {
        foreach my $read ( keys %{ $fh_for->{$tag} } ) {
            close $fh_for->{$tag}->{$read};
        }
    }

    return;
}

=func convert_tag_to_regexp

  Usage       : %re_for = convert_tag_to_regexp( 'NNNNBGAGGC', 'NNNNBAGAAG' );
  Purpose     : Convert tags to regular expressions for matching
  Returns     : Hash (
                    String (tag) => Arrayref (of Regexps)
                )
  Parameters  : Array of strings (the tags)
  Throws      : No exceptions
  Comments    : None

=cut

sub convert_tag_to_regexp {
    my @tags = @_;

    my %re_for;
    foreach my $tag (@tags) {
        my @mismatch_tags = ($tag);    # Start with tag without mismatches

        # Add tag with each possible mismatch
        foreach my $i ( 0 .. length($tag) - 1 ) {
            my $mismatch_tag = $tag;
            my $base = substr $mismatch_tag, $i, 1, q{N};    # Replace with N
            if ( $base ne q{N} ) {

                # Not completely random base already
                push @mismatch_tags, $mismatch_tag;
            }
        }

        # Convert IUPAC codes to AGCT (or N)
        foreach my $re (@mismatch_tags) {
            $re =~ s/N/[NAGCT]/xmsg;    # Random bases can be called as N
            $re =~ s/B/[GCT]/xmsg;
            $re =~ s/D/[AGT]/xmsg;
            $re =~ s/H/[ACT]/xmsg;
            $re =~ s/V/[AGC]/xmsg;
            $re =~ s/R/[AG]/xmsg;
            $re =~ s/Y/[CT]/xmsg;
            $re =~ s/K/[GT]/xmsg;
            $re =~ s/M/[AC]/xmsg;
            $re =~ s/S/[GC]/xmsg;
            $re =~ s/W/[AT]/xmsg;
            push @{ $re_for{$tag} }, qr/\A $re \Z/xms;
        }
    }

    return %re_for;
}

1;
