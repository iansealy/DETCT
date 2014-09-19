## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::BAM;
## use critic

# ABSTRACT: Miscellaneous functions for interacting with BAM files

## Author         : is1
## Maintainer     : is1
## Created        : 2012-09-20
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
use Bio::DB::Sam;
use List::Util qw( min sum );
use List::MoreUtils qw( all );
use Data::Compare;
use DETCT::Misc::Tag;
use DETCT::Misc::BAM::Flag;

use base qw( Exporter );
our @EXPORT_OK = qw(
  get_reference_sequence_lengths
  get_sequence
  count_tags
  bin_reads
  get_read_peaks
  get_three_prime_ends
  merge_three_prime_ends
  filter_three_prime_ends
  choose_three_prime_end
  count_reads
  merge_read_counts
  stats_by_tag
  stats_all_reads
  downsample_by_tag
  downsample_all_reads
  mark_duplicates
  filter_by_tag
);

=head1 SYNOPSIS

    # Brief code examples

=cut

# Constants

# Regexps for checking for polyA
Readonly our @POLYA_REGEXP => (
    qr/\A AAA.AAA... \z/xms,
    qr/\A AAA.AA.A.. \z/xms,
    qr/\A AAA.A.AA.. \z/xms,
    qr/\A AA.AAAA... \z/xms,
    qr/\A AA.AAA.A.. \z/xms,
    qr/\A AA.A.AAA.. \z/xms,
    qr/\A A.AAAAA... \z/xms,
    qr/\A A.AAAA.A.. \z/xms,
    qr/\A A.AAA.AA.. \z/xms,
    qr/\A A.AA.AAA.. \z/xms,
    qr/\A A.A.AAAA.. \z/xms,
    qr/\A AA.AA.AA.. \z/xms,
);

=func get_reference_sequence_lengths

  Usage       : my %length_of
                    = DETCT::Misc::BAM::get_reference_sequence_lengths($bam_file);
  Purpose     : Get length of each reference sequence from a BAM file
  Returns     : Hash (
                    seq_region => length
                )
  Parameters  : String (the BAM file)
  Throws      : If BAM file is missing
  Comments    : None

=cut

sub get_reference_sequence_lengths {
    my ($bam_file) = @_;

    confess 'No BAM file specified' if !defined $bam_file;

    my $sam = Bio::DB::Sam->new( -bam => $bam_file );

    my %length_of;

    foreach my $seq_id ( $sam->seq_ids ) {
        $length_of{$seq_id} = $sam->length($seq_id);
    }

    return %length_of;
}

=func get_sequence

  Usage       : my $seq = DETCT::Misc::BAM::get_sequence( {
                    fasta_index => $fai,
                    seq_name    => '1',
                    start       => 1,
                    end         => 1000,
                    strand      => 1,
                } );
  Purpose     : Get sequence from FASTA file
  Returns     : String (sequence)
  Parameters  : Hashref {
                    fasta_index => Bio::DB::Sam::Fai
                    ref_fasta   => String (the FASTA file)
                    seq_name    => String (the sequence name)
                    start       => Int (the sequence start)
                    end         => Int (the sequence end)
                    strand      => Int (the sequence strand)
                }
  Throws      : If FASTA index and file are both missing
                If sequence name is missing
                If sequence start is missing
                If sequence end is missing
                If sequence strand is missing
  Comments    : None

=cut

sub get_sequence {
    my ($arg_ref) = @_;

    confess 'No FASTA index or FASTA file specified'
      if !defined $arg_ref->{fasta_index} && !defined $arg_ref->{ref_fasta};
    confess 'No sequence name specified'   if !defined $arg_ref->{seq_name};
    confess 'No sequence start specified'  if !defined $arg_ref->{start};
    confess 'No sequence end specified'    if !defined $arg_ref->{end};
    confess 'No sequence strand specified' if !defined $arg_ref->{strand};

    my $fai =
        $arg_ref->{fasta_index}
      ? $arg_ref->{fasta_index}
      : Bio::DB::Sam::Fai->load( $arg_ref->{ref_fasta} );

    my $query = sprintf '%s:%d-%d', $arg_ref->{seq_name}, $arg_ref->{start},
      $arg_ref->{end};

    my $seq = uc $fai->fetch($query);

    if ( $arg_ref->{strand} == -1 ) {    ## no critic (ProhibitMagicNumbers)
        $seq = reverse $seq;
        $seq =~ tr/ACGT/TGCA/;
    }

    return $seq;
}

=func count_tags

  Usage       : my $count_ref = DETCT::Misc::BAM::count_tags( {
                    bam_file           => $bam_file,
                    mismatch_threshold => 2,
                    seq_name           => '1',
                    start              => 1,
                    end                => 1000,
                    tags               => ['NNNNBGAGGC', 'NNNNBAGAAG'],
                } );
  Purpose     : Count tags and random bases in a BAM file
  Returns     : Hashref {
                    String (tag) => Hashref {
                        String (random bases) => Int (count)
                    }
                }
  Parameters  : Hashref {
                    bam_file           => String (the BAM file)
                    mismatch_threshold => Int (the mismatch threshold)
                    seq_name           => String (the sequence name)
                    start              => Int (the start) or undef
                    end                => Int (the end) or undef
                    tags               => Arrayref of strings (the tags)
                }
  Throws      : If BAM file is missing
                If mismatch threshold is missing
                If tags are missing
  Comments    : None

=cut

sub count_tags {
    my ($arg_ref) = @_;

    confess 'No BAM file specified' if !defined $arg_ref->{bam_file};
    confess 'No mismatch threshold specified'
      if !defined $arg_ref->{mismatch_threshold};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No tags specified'          if !defined $arg_ref->{tags};

    my @tags = @{ $arg_ref->{tags} };

    # Convert tags to regular expressions
    my %re_for = DETCT::Misc::Tag::convert_tag_to_regexp(@tags);

    # Count random bases per tag
    my %random_count_for;
    foreach my $tag (@tags) {
        my @random = $tag =~ m/[NRYKMSWBDHV]/xmsg;
        $random_count_for{$tag} = scalar @random;
    }

    my $sam = Bio::DB::Sam->new( -bam => $arg_ref->{bam_file} );

    my %count;

    # Callback for filtering
    my $callback = sub {
        my ($alignment) = @_;
        return if !is_read2($alignment);
        return if is_duplicate($alignment);
        return if $alignment->unmapped;
        return
          if above_mismatch_threshold( $alignment,
            $arg_ref->{mismatch_threshold} );

        # Match tag
        my ($tag_in_read) = $alignment->query->name =~ m/[#] ([NAGCT]+) \z/xmsg;
        return if !$tag_in_read;
      TAG: foreach my $tag ( sort keys %re_for ) {
            my $regexps = $re_for{$tag};
            foreach my $re ( @{$regexps} ) {
                if ( $tag_in_read =~ $re ) {
                    my $random = substr $tag_in_read, 0,
                      $random_count_for{$tag};
                    $count{$tag}{$random}++;
                    last TAG;
                }
            }
        }

        return;
    };

    # Construct region
    my $region = $arg_ref->{seq_name};
    if ( exists $arg_ref->{start} ) {
        $region .= q{:} . $arg_ref->{start};
        if ( exists $arg_ref->{end} ) {
            $region .= q{-} . $arg_ref->{end};
        }
    }

    $sam->fetch( $region, $callback );

    return \%count;
}

=func bin_reads

  Usage       : my $bin_ref = DETCT::Misc::BAM::bin_reads( {
                    bam_file           => $bam_file,
                    mismatch_threshold => 2,
                    bin_size           => 100,
                    seq_name           => '1',
                    tags               => ['NNNNBGAGGC', 'NNNNBAGAAG'],
                } );
  Purpose     : Bin read 2s in a BAM file
  Returns     : Hashref {
                    Int (1 or -1) (strand) => Hashref {
                        Int (bin) => Int (count)
                    }
                }
  Parameters  : Hashref {
                    bam_file           => String (the BAM file)
                    mismatch_threshold => Int (the mismatch threshold)
                    bin_size           => Int (the bin size)
                    seq_name           => String (the sequence name)
                    tags               => Arrayref of strings (the tags)
                }
  Throws      : If BAM file is missing
                If mismatch threshold is missing
                If bin size is missing
                If sequence name is missing
                If tags are missing
  Comments    : None

=cut

sub bin_reads {
    my ($arg_ref) = @_;

    confess 'No BAM file specified' if !defined $arg_ref->{bam_file};
    confess 'No mismatch threshold specified'
      if !defined $arg_ref->{mismatch_threshold};
    confess 'No bin size specified'      if !defined $arg_ref->{bin_size};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No tags specified'          if !defined $arg_ref->{tags};

    my @tags = @{ $arg_ref->{tags} };

    # Convert tags to regular expressions
    my %re_for = DETCT::Misc::Tag::convert_tag_to_regexp(@tags);

    my $sam = Bio::DB::Sam->new( -bam => $arg_ref->{bam_file} );

    my %read_count_for = (
        '1'  => {},
        '-1' => {},
    );

    # Callback for filtering
    my $callback = sub {
        my ($alignment) = @_;
        return if !is_read2($alignment);
        return if is_duplicate($alignment);
        return if $alignment->unmapped;
        return
          if above_mismatch_threshold( $alignment,
            $arg_ref->{mismatch_threshold} );
        return if !matched_tag( $alignment, \%re_for );

        # Read can span multiple bins
        my $start_bin = int( ( $alignment->start - 1 ) / $arg_ref->{bin_size} );
        my $end_bin   = int( ( $alignment->end - 1 ) / $arg_ref->{bin_size} );

        foreach my $bin ( $start_bin .. $end_bin ) {
            $read_count_for{ $alignment->strand }{$bin}++;
        }

        return;
    };

    $sam->fetch( $arg_ref->{seq_name}, $callback );

    return { $arg_ref->{seq_name} => \%read_count_for };
}

=func get_read_peaks

  Usage       : my $peaks_ref = DETCT::Misc::BAM::get_read_peaks( {
                    bam_file           => $bam_file,
                    mismatch_threshold => 2,
                    peak_buffer_width  => 100,
                    seq_name           => '1',
                    tags               => ['NNNNBGAGGC', 'NNNNBAGAAG'],
                } );
  Purpose     : Get read 2 peaks (overlapping reads) for a BAM file
  Returns     : Hashref {
                    String (sequence name) => Hashref {
                        Int (1 or -1) (strand) => Arrayref [
                            Arrayref [
                                Int (peak start),
                                Int (peak end),
                                Int (peak read count),
                            ],
                            ... (peaks)
                        ]
                    }
                }
  Parameters  : Hashref {
                    bam_file           => String (the BAM file)
                    mismatch_threshold => Int (the mismatch threshold)
                    peak_buffer_width  => Int (the peak buffer size),
                    seq_name           => String (the sequence name)
                    tags               => Arrayref of strings (the tags)
                }
  Throws      : If BAM file is missing
                If mismatch threshold is missing
                If peak buffer width is missing
                If sequence name is missing
                If tags are missing
  Comments    : BAM file must be sorted by coordinate

=cut

sub get_read_peaks {
    my ($arg_ref) = @_;

    confess 'No BAM file specified' if !defined $arg_ref->{bam_file};
    confess 'No mismatch threshold specified'
      if !defined $arg_ref->{mismatch_threshold};
    confess 'No peak buffer width specified'
      if !defined $arg_ref->{peak_buffer_width};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No tags specified'          if !defined $arg_ref->{tags};

    my @tags = @{ $arg_ref->{tags} };

    # Convert tags to regular expressions
    my %re_for = DETCT::Misc::Tag::convert_tag_to_regexp(@tags);

    my $sam = Bio::DB::Sam->new( -bam => $arg_ref->{bam_file} );

    # Peak variables (all keyed by strand)
    my %peaks = (
        '1'  => [],
        '-1' => [],
    );
    my %current_peak_read_count;
    my %current_peak_start;
    my %current_peak_end;

    # Read variables (all keyed by strand)
    my %current_read_start;
    my %current_read_end;

    # Callback for filtering
    my $callback = sub {
        my ($alignment) = @_;
        return if !is_read2($alignment);
        return if is_duplicate($alignment);
        return if $alignment->unmapped;
        return
          if above_mismatch_threshold( $alignment,
            $arg_ref->{mismatch_threshold} );
        return if !matched_tag( $alignment, \%re_for );

        my $strand = $alignment->strand;

        $current_read_start{$strand} = $alignment->start;
        $current_read_end{$strand}   = $alignment->end;

        # We're starting the first peak
        if ( !exists $current_peak_start{$strand} ) {
            $current_peak_start{$strand}      = $current_read_start{$strand};
            $current_peak_end{$strand}        = $current_read_end{$strand};
            $current_peak_read_count{$strand} = 1;
            return;
        }

        # Extend or finish current peak?
        if ( $current_read_start{$strand} - $current_peak_end{$strand} <
            $arg_ref->{peak_buffer_width} )
        {
            # Extend current peak
            $current_peak_end{$strand} = $current_read_end{$strand};
            $current_peak_read_count{$strand}++;
        }
        else {
            # Finish current peak
            push @{ $peaks{$strand} },
              [
                $current_peak_start{$strand}, $current_peak_end{$strand},
                $current_peak_read_count{$strand}
              ];

            # Start new peak
            $current_peak_start{$strand}      = $current_read_start{$strand};
            $current_peak_end{$strand}        = $current_read_end{$strand};
            $current_peak_read_count{$strand} = 1;
        }

        return;
    };

    # Identify peaks (where peaks are read 2s separated by a buffer of specific
    # size)
    $sam->fetch( $arg_ref->{seq_name}, $callback );

    # Finish last peaks
    foreach my $strand ( 1, -1 ) {    ## no critic (ProhibitMagicNumbers)
        if ( $current_peak_read_count{$strand} ) {
            push @{ $peaks{$strand} },
              [
                $current_peak_start{$strand}, $current_peak_end{$strand},
                $current_peak_read_count{$strand}
              ];
        }
    }

    return { $arg_ref->{seq_name} => \%peaks };
}

=func get_three_prime_ends

  Usage       : my $three_prime_ref = DETCT::Misc::BAM::get_three_prime_ends( {
                    bam_file           => $bam_file,
                    mismatch_threshold => 2,
                    seq_name           => '1',
                    strand             => 1,
                    tags               => ['NNNNBGAGGC', 'NNNNBAGAAG'],
                    regions            => $regions_ary_ref,
                } );
  Purpose     : Get all 3' ends for a list of regions
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
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
                                ],
                                ... (3' ends)
                            ],
                        ],
                        ... (regions)
                }
  Parameters  : Hashref {
                    bam_file           => String (the BAM file)
                    mismatch_threshold => Int (the mismatch threshold)
                    seq_name           => String (the sequence name)
                    strand             => Int ( 1 or -1 ) (the 3' end strand)
                    tags               => Arrayref of strings (the tags)
                    regions            => Arrayref (of regions)
                }
  Throws      : If BAM file is missing
                If mismatch threshold is missing
                If sequence name is missing
                If strand is missing
                If tags are missing
                If regions are missing
  Comments    : regions parameter is a list of regions, unlike the regions
                parameter for merge_three_prime_ends where it is a list of lists
                of regions
                strand parameter is strand of 3' end, which is the same as the
                strand of read 2

=cut

sub get_three_prime_ends {
    my ($arg_ref) = @_;

    confess 'No BAM file specified' if !defined $arg_ref->{bam_file};
    confess 'No mismatch threshold specified'
      if !defined $arg_ref->{mismatch_threshold};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No strand specified'        if !defined $arg_ref->{strand};
    confess 'No tags specified'          if !defined $arg_ref->{tags};
    confess 'No regions specified'       if !defined $arg_ref->{regions};

    my @tags = @{ $arg_ref->{tags} };

    my $three_prime_strand = $arg_ref->{strand};

    # Convert tags to regular expressions
    my %re_for = DETCT::Misc::Tag::convert_tag_to_regexp(@tags);

    my $sam = Bio::DB::Sam->new( -bam => $arg_ref->{bam_file} );

    my @regions_with_three_prime_ends;

    foreach my $region ( @{ $arg_ref->{regions} } ) {
        my ( $start, $end, $max_read_count, $log_prob_sum ) = @{$region};

        my %count_for;

        # Get all second reads in region
        my $read2_alignments = $sam->features(
            -seq_id   => $arg_ref->{seq_name},
            -start    => $start,
            -end      => $end,
            -flags    => { SECOND_MATE => 1 },
            -iterator => 1,
        );

        # Get all 3' ends
        while ( my $alignment = $read2_alignments->next_seq ) {
            next if is_duplicate($alignment);

            # next if $alignment->unmapped; # Not needed; always mapped
            next if $alignment->munmapped;    # Want read 1 mapped too
            next
              if above_mismatch_threshold( $alignment,
                $arg_ref->{mismatch_threshold} );
            next if !matched_tag( $alignment, \%re_for );

            # Strand of read 1 is opposite to 3' end strand
            next if $alignment->mstrand == $three_prime_strand;

            # Strand of read 2 is same as 3' end strand
            next if $alignment->strand != $three_prime_strand;

            # Skip if 3' end is on a different chromosome
            # Hopefully not significant number of real 3' ends on different
            # chromosomes because are hard to deal with
            # If reads are on different chromosomes then TLEN will be 0 and
            # mate_end will return undefined (i.e. can't get 3' end position
            # without querying by read name, which is slow for a BAM file
            # sorted by coordinate)
            next if $alignment->mate_seq_id ne $arg_ref->{seq_name};

            # Identify 3' end position based on alignment of read 1
            my $three_prime_seq = $alignment->mate_seq_id;
            my $three_prime_pos =
                $three_prime_strand == 1
              ? $alignment->mate_end
              : $alignment->mate_start;

            # Count number of reads supporting each 3' end
            my $three_prime = join q{:}, $three_prime_seq, $three_prime_pos;
            $count_for{$three_prime}++;
        }

        # Turn counts into an array
        my @three_prime_ends;
        foreach my $three_prime (
            reverse sort { $count_for{$a} <=> $count_for{$b} }
            keys %count_for
          )
        {
            my ( $seq, $pos ) = split /:/xms, $three_prime;
            push @three_prime_ends,
              [ $seq, $pos, $three_prime_strand, $count_for{$three_prime} ];
        }

        # Add three prime ends to regions
        push @regions_with_three_prime_ends,
          [
            $start,        $end,                $max_read_count,
            $log_prob_sum, $three_prime_strand, \@three_prime_ends,
          ];
    }

    return { $arg_ref->{seq_name} => \@regions_with_three_prime_ends };
}

=func merge_three_prime_ends

  Usage       : my $three_prime_ref
                    = DETCT::Misc::BAM::merge_three_prime_ends( {
                    seq_name => '1',
                    regions  => $regions_ary_ref,
                } );
  Purpose     : Merge multiple lists of regions with 3' ends
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
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
                                ],
                                ... (3' ends)
                            ],
                        ],
                        ... (regions)
                }
  Parameters  : Hashref {
                    seq_name => String (the sequence name)
                    regions  => Arrayref (of arrayrefs of regions)
                }
  Throws      : If sequence name is missing
                If regions are missing
                If each list of regions doesn't have same number of regions
                If regions are not in the same order or not the same in each
                list
  Comments    : regions parameter is a list of lists of regions, unlike
                the regions parameter for get_three_prime_ends where it is a
                list of regions

=cut

sub merge_three_prime_ends {
    my ($arg_ref) = @_;

    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No regions specified'       if !defined $arg_ref->{regions};

    my @list_of_lists_of_regions = @{ $arg_ref->{regions} };

    # No need to merge if only one list of regions
    my $num_lists = scalar @list_of_lists_of_regions;
    if ( $num_lists == 1 ) {
        return { $arg_ref->{seq_name} => $list_of_lists_of_regions[0] };
    }

    # Ensure each list has same number of regions as first list
    my $num_regions1 = scalar @{ $list_of_lists_of_regions[0] };
    foreach my $list_index ( 1 .. $num_lists - 1 ) {
        my $num_regions2 = scalar @{ $list_of_lists_of_regions[$list_index] };
        if ( $num_regions1 != $num_regions2 ) {
            confess 'Number of regions does not match in all lists';
        }
    }

    my @regions_with_three_prime_ends;

    # Merge all lists
    foreach my $region_index ( 0 .. $num_regions1 - 1 ) {

        # Ensure region from first list is same in each list
        my $region1 = $list_of_lists_of_regions[0]->[$region_index];
        my ( $start1, $end1, $max_read_count1, $log_prob_sum1, $strand1 ) =
          @{$region1};
        foreach my $list_index ( 1 .. $num_lists - 1 ) {
            my $region2 =
              $list_of_lists_of_regions[$list_index]->[$region_index];
            my ( $start2, $end2, $max_read_count2, $log_prob_sum2, $strand2 ) =
              @{$region2};
            if (   $start1 != $start2
                || $end1 != $end2
                || $max_read_count1 != $max_read_count2
                || $log_prob_sum1 != $log_prob_sum2
                || $strand1 != $strand2 )
            {
                confess
                  'Regions not in the same order or not the same in each list';
            }
        }

        # Get all the 3' ends
        my @unmerged_three_prime_ends;
        foreach my $list_index ( 0 .. $num_lists - 1 ) {
            my $list   = $list_of_lists_of_regions[$list_index];
            my $region = $list->[$region_index];
            my ( undef, undef, undef, undef, undef, $three_prime_ends ) =
              @{$region};
            push @unmerged_three_prime_ends, @{$three_prime_ends};
        }

        # Add up counts for identical 3' ends
        my %count_for;
        foreach my $three_prime_end (@unmerged_three_prime_ends) {
            my ( $seq, $pos, $strand, $read_count ) = @{$three_prime_end};
            my $three_prime = join q{:}, $seq, $pos, $strand;
            $count_for{$three_prime} += $read_count;
        }

        # Turn counts into an array
        my @three_prime_ends;
        foreach my $three_prime (
            reverse sort { $count_for{$a} <=> $count_for{$b} }
            keys %count_for
          )
        {
            my ( $seq, $pos, $strand ) = split /:/xms, $three_prime;
            push @three_prime_ends,
              [ $seq, $pos, $strand, $count_for{$three_prime} ];
        }

        # Add three prime ends to regions
        push @regions_with_three_prime_ends,
          [
            $start1,        $end1,    $max_read_count1,
            $log_prob_sum1, $strand1, \@three_prime_ends,
          ];

        $region_index++;
    }

    return { $arg_ref->{seq_name} => \@regions_with_three_prime_ends };
}

=func filter_three_prime_ends

  Usage       : my $three_prime_ref
                    = DETCT::Misc::BAM::filter_three_prime_ends( {
                    analysis => $analysis,
                    seq_name => '1',
                    regions  => $regions_ary_ref,
                } );
  Purpose     : Filter list of regions with 3' ends
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
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
                                ],
                                ... (3' ends)
                            ],
                        ],
                        ... (regions)
                }
  Parameters  : Hashref {
                    analysis => DETCT::Analysis
                    seq_name => String (the sequence name)
                    regions  => Arrayref (of regions)
                }
  Throws      : If analysis is missing
                If sequence name is missing
                If regions are missing
  Comments    : regions parameter is a list of regions, unlike the regions
                parameter for merge_three_prime_ends where it is a list of lists
                of regions

=cut

sub filter_three_prime_ends {
    my ($arg_ref) = @_;

    confess 'No analysis specified'      if !defined $arg_ref->{analysis};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No regions specified'       if !defined $arg_ref->{regions};

    my @regions_with_three_prime_ends;

    # Iterate over regions
    foreach my $region ( @{ $arg_ref->{regions} } ) {
        my ( $region_start, $region_end, $region_max_read_count,
            $region_log_prob_sum, $three_prime_strand,
            $unfiltered_three_prime_ends )
          = @{$region};

        # Filter 3' ends
        my @three_prime_ends;
        foreach my $three_prime_end ( @{$unfiltered_three_prime_ends} ) {
            my ( $seq_name, $pos, $strand, $read_count ) = @{$three_prime_end};

            # Must be supported by more than 3 reads
            next if $read_count <= 3;    ## no critic (ProhibitMagicNumbers)

            # Check 10 bp downstream of 3' end for polyA
            my $ten_bp_start;
            my $ten_bp_end;
            if ( $strand == 1 ) {
                $ten_bp_start = $pos + 1;
                $ten_bp_end   = $pos + 10;   ## no critic (ProhibitMagicNumbers)
            }
            else {
                $ten_bp_start = $pos - 10;   ## no critic (ProhibitMagicNumbers)
                $ten_bp_end   = $pos - 1;
            }
            my $ten_bp_seq =
              $arg_ref->{analysis}
              ->get_subsequence( $seq_name, $ten_bp_start, $ten_bp_end,
                $strand );

            # Check if 10 bp downstream is polyA
            next if is_polya($ten_bp_seq);

            push @three_prime_ends, $three_prime_end;
        }

        # Add three prime ends to regions
        push @regions_with_three_prime_ends,
          [
            $region_start,          $region_end,
            $region_max_read_count, $region_log_prob_sum,
            $three_prime_strand,    \@three_prime_ends,
          ];
    }

    return { $arg_ref->{seq_name} => \@regions_with_three_prime_ends };
}

=func choose_three_prime_end

  Usage       : my $three_prime_ref
                    = DETCT::Misc::BAM::choose_three_prime_end( {
                    seq_name => '1',
                    regions  => $regions_ary_ref,
                } );
  Purpose     : Filter and adjust list of regions and choose best 3' end
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
                        Arrayref [
                            Int (region start),
                            Int (region end),
                            Int (region maximum read count),
                            Float (region log probability sum),
                            String (3' end sequence name) or undef,
                            Int (3' end position) or undef,
                            Int (3' end strand),
                            Int (3' end read count) or undef,
                        ],
                        ... (regions)
                }
  Parameters  : Hashref {
                    seq_name => String (the sequence name)
                    regions  => Arrayref (of regions)
                }
  Throws      : If sequence name is missing
                If regions are missing
  Comments    : regions parameter is a list of regions, unlike the regions
                parameter for merge_three_prime_ends where it is a list of lists
                of regions

=cut

sub choose_three_prime_end {
    my ($arg_ref) = @_;

    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No regions specified'       if !defined $arg_ref->{regions};

    my @regions_with_three_prime_ends;

    # Iterate over regions
    foreach my $region ( @{ $arg_ref->{regions} } ) {
        my ( $region_start, $region_end, $region_max_read_count,
            $region_log_prob_sum, $strand, $three_prime_ends )
          = @{$region};

        my (
            $three_prime_seq_name, $three_prime_pos,
            $three_prime_strand,   $three_prime_read_count
        );

        @{$three_prime_ends} = reverse sort {
            _sort_three_prime_end( $a, $b, $arg_ref->{seq_name}, $region_start,
                $region_end )
        } @{$three_prime_ends};

        # Get best 3' end (highest read count)
        if ( @{$three_prime_ends} ) {
            (
                $three_prime_seq_name, $three_prime_pos,
                $three_prime_strand,   $three_prime_read_count
            ) = @{ $three_prime_ends->[0] };
        }

        # Reduce size of region if appropriate
        ## no critic (ProhibitMagicNumbers)
        if ( defined $three_prime_seq_name
            && $three_prime_seq_name eq $arg_ref->{seq_name} )
        {
            if (   $three_prime_strand == 1
                && $three_prime_pos < $region_end
                && $three_prime_pos > $region_start )
            {
                $region_end = $three_prime_pos;
            }
            elsif ($three_prime_strand == -1
                && $three_prime_pos > $region_start
                && $three_prime_pos < $region_end )
            {
                $region_start = $three_prime_pos;
            }
        }
        ## use critic

        # Use strand for region (i.e. based on read 2 alignments) if necessary
        $three_prime_strand = $three_prime_strand || $strand;

        # Add three prime ends to regions
        push @regions_with_three_prime_ends,
          [
            $region_start,          $region_end,
            $region_max_read_count, $region_log_prob_sum,
            $three_prime_seq_name,  $three_prime_pos,
            $three_prime_strand,    $three_prime_read_count,
          ];
    }

    return { $arg_ref->{seq_name} => \@regions_with_three_prime_ends };
}

# Sort by read count then distance to region
sub _sort_three_prime_end {
    my ( $a, $b, $seq_name, $region_start, $region_end ) = @_;

    my $seq_name_a = $a->[0];
    my $seq_name_b = $b->[0];
    my $pos_a      = $a->[1];
    my $pos_b      = $b->[1];
    ## no critic (ProhibitMagicNumbers)
    my $read_count_a = $a->[3];
    my $read_count_b = $b->[3];
    ## use critic

    # Get minimum distance to region
    my $dist_a = min( abs $region_start - $pos_a, abs $region_end - $pos_a );
    my $dist_b = min( abs $region_start - $pos_b, abs $region_end - $pos_b );

    # Make sure 3' end is on same chromosome as region
    # (1e+100 is bigger than any chromosome can be to ensure sorting last)
    if ( $seq_name_a ne $seq_name ) {
        $dist_a = 1e+100;    ## no critic (ProhibitMagicNumbers)
    }
    if ( $seq_name_b ne $seq_name ) {
        $dist_b = 1e+100;    ## no critic (ProhibitMagicNumbers)
    }

    return $read_count_a <=> $read_count_b || $dist_b <=> $dist_a;
}

=func count_reads

  Usage       : my $count_ref = DETCT::Misc::BAM::count_reads( {
                    bam_file           => $bam_file,
                    mismatch_threshold => 2,
                    seq_name           => '1',
                    regions            => $regions_ary_ref,
                    tags               => ['NNNNBGAGGC', 'NNNNBAGAAG'],
                } );
  Purpose     : Count reads in regions of a BAM file
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
                        Arrayref [
                            Int (region start),
                            Int (region end),
                            Int (region maximum read count),
                            Float (region log probability sum),
                            String (3' end sequence name) or undef,
                            Int (3' end position) or undef,
                            Int (3' end strand) or undef,
                            Int (3' end read count) or undef,
                            Hashref {
                                String (tag) => Int (count)
                            }
                        ],
                        ... (regions)
                }
  Parameters  : Hashref {
                    bam_file           => String (the BAM file)
                    mismatch_threshold => Int (the mismatch threshold)
                    seq_name           => String (the sequence name) or undef
                    regions            => Arrayref (of regions)
                    tags               => Arrayref of strings (the tags)
                }
  Throws      : If BAM file is missing
                If mismatch threshold is missing
                If sequence name is missing
                If regions are missing
                If tags are missing
  Comments    : regions parameter is a list of regions, unlike the regions
                parameter for merge_read_counts where it is a hash keyed by BAM
                file with values being lists of regions

=cut

sub count_reads {
    my ($arg_ref) = @_;

    confess 'No BAM file specified' if !defined $arg_ref->{bam_file};
    confess 'No mismatch threshold specified'
      if !defined $arg_ref->{mismatch_threshold};
    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No regions specified'       if !defined $arg_ref->{regions};
    confess 'No tags specified'          if !defined $arg_ref->{tags};

    my @tags = @{ $arg_ref->{tags} };

    # Convert tags to regular expressions
    my %re_for = DETCT::Misc::Tag::convert_tag_to_regexp(@tags);

    my $sam = Bio::DB::Sam->new( -bam => $arg_ref->{bam_file} );

    my @regions_with_three_prime_ends;

    # Iterate over regions
    foreach my $region ( @{ $arg_ref->{regions} } ) {
        my (
            $region_start,          $region_end,
            $region_max_read_count, $region_log_prob_sum,
            $three_prime_seq_name,  $three_prime_pos,
            $three_prime_strand,    $three_prime_read_count
        ) = @{$region};

        my %count = map { $_ => 0 } @tags;

        # Get first read from each pair
        my $read2_alignments = $sam->features(
            -seq_id   => $arg_ref->{seq_name},
            -start    => $region_start,
            -end      => $region_end,
            -flags    => { SECOND_MATE => 1 },
            -iterator => 1,
        );
        while ( my $alignment = $read2_alignments->next_seq ) {
            next if is_duplicate($alignment);

            #next if $alignment->unmapped; # Not needed; always mapped
            next
              if above_mismatch_threshold( $alignment,
                $arg_ref->{mismatch_threshold} );

            # Only count reads on 3' end strand
            # (which is the same as read 2 strand)
            next if $alignment->strand != $three_prime_strand;

            # Match tag
            my $tag = matched_tag( $alignment, \%re_for );
            next if !$tag;
            $count{$tag}++;
        }

        # Add read counts to regions
        push @regions_with_three_prime_ends,
          [
            $region_start,          $region_end,
            $region_max_read_count, $region_log_prob_sum,
            $three_prime_seq_name,  $three_prime_pos,
            $three_prime_strand,    $three_prime_read_count,
            \%count,
          ];
    }

    return { $arg_ref->{seq_name} => \@regions_with_three_prime_ends };
}

=func merge_read_counts

  Usage       : my $count_ref
                    = DETCT::Misc::BAM::merge_read_counts( {
                    seq_name => '1',
                    regions  => $regions_hash_ref,
                    samples  => $samples_ary_ref,
                } );
  Purpose     : Merge multiple lists of regions with read counts
  Returns     : Hashref {
                    String (sequence name) => Arrayref [
                        Arrayref [
                            Int (region start),
                            Int (region end),
                            Int (region maximum read count),
                            Float (region log probability sum),
                            String (3' end sequence name) or undef,
                            Int (3' end position) or undef,
                            Int (3' end strand) or undef,
                            Int (3' end read count) or undef,
                            Arrayref [
                                Int (count)
                                ...
                            ]
                        ],
                        ... (regions)
                }
  Parameters  : Hashref {
                    seq_name => String (the sequence name)
                    regions  => Arrayref (of arrayrefs of regions)
                    samples  => Arrayref (of samples)
                }
  Throws      : If sequence name is missing
                If regions are missing
                If samples are missing
                If each list of regions doesn't have same number of regions
                If regions are not in the same order or not the same in each
                list
  Comments    : regions parameter is a hash keyed by BAM file with values being
                lists of regions, unlike the regions parameter for count_reads
                where it is a list of regions

=cut

sub merge_read_counts {
    my ($arg_ref) = @_;

    confess 'No sequence name specified' if !defined $arg_ref->{seq_name};
    confess 'No regions specified'       if !defined $arg_ref->{regions};
    confess 'No samples specified'       if !defined $arg_ref->{samples};

    my %hash_of_lists_of_regions = %{ $arg_ref->{regions} };

    # Ensure each list has same number of regions
    my @bam_files    = keys %hash_of_lists_of_regions;
    my $num_regions1 = scalar @{ $hash_of_lists_of_regions{ $bam_files[0] } };
    foreach my $list_index ( 1 .. scalar @bam_files - 1 ) {
        my $num_regions2 =
          scalar @{ $hash_of_lists_of_regions{ $bam_files[$list_index] } };
        if ( $num_regions1 != $num_regions2 ) {
            confess 'Number of regions does not match in all lists';
        }
    }

    # Get index for each sample
    my %sample_index_for;
    my $index = 0;
    foreach my $sample ( @{ $arg_ref->{samples} } ) {
        my $bam_file = $sample->bam_file;
        my $tag      = $sample->tag;
        $sample_index_for{$bam_file}{$tag} = $index;
        $index++;
    }

    my @regions_with_three_prime_ends;

    # Merge all lists
    foreach my $region_index ( 0 .. $num_regions1 - 1 ) {

        # Ensure regions are the same in each list
        my $region1 =
          $hash_of_lists_of_regions{ $bam_files[0] }->[$region_index];
        my @region1 = @{$region1}[ 0 .. 7 ]; ## no critic (ProhibitMagicNumbers)
        foreach my $list_index ( 1 .. scalar @bam_files - 1 ) {
            my $region2 =
              $hash_of_lists_of_regions{ $bam_files[$list_index] }
              ->[$region_index];

            # Check first 8 fields of each region are identical
            my @region2 =
              @{$region2}[ 0 .. 7 ];         ## no critic (ProhibitMagicNumbers)
            if ( !Compare( \@region1, \@region2 ) ) {
                confess
                  'Regions not in the same order or not the same in each list';
            }
        }

        my @read_counts;

        # Get read count for each BAM file / tag
        foreach my $bam_file (@bam_files) {
            my $region = $hash_of_lists_of_regions{$bam_file}->[$region_index];
            my $read_counts_ref = $region->[-1];    # Read counts are last field
            foreach my $tag ( keys %{$read_counts_ref} ) {
                my $read_count = $read_counts_ref->{$tag};
                if ( !exists $sample_index_for{$bam_file}{$tag} ) {
                    confess "Unknown BAM file ($bam_file) / tag ($tag) pair";
                }
                my $sample_index = $sample_index_for{$bam_file}{$tag};
                $read_counts[$sample_index] = $read_count;
            }
        }

        push @regions_with_three_prime_ends, [ @region1, \@read_counts ];

        $region_index++;
    }

    return { $arg_ref->{seq_name} => \@regions_with_three_prime_ends };
}

=func stats_by_tag

  Usage       : my $stats = DETCT::Misc::BAM::stats_by_tag( {
                    bam_file => $bam_file,
                    tags     => ['NNNNBGAGGC', 'NNNNBAGAAG'],
                } );
  Purpose     : Generate stats per tag in a BAM file
  Returns     : Hashref {
                    String (tag) => Hashref {
                        paired => Int (paired read count),
                        mapped => Int (mapped paired read count),
                        proper => Int (properly paired read count),
                    }
                }
  Parameters  : Hashref {
                    bam_file       => String (the BAM file)
                    tags           => Arrayref of strings (the tags)
                    skip_sequences => Arrayref of strings (the skip sequences)
                }
  Throws      : If BAM file is missing
                If tags are missing
  Comments    : None

=cut

sub stats_by_tag {
    my ($arg_ref) = @_;

    confess 'No BAM file specified' if !defined $arg_ref->{bam_file};
    confess 'No tags specified'     if !defined $arg_ref->{tags};

    my @tags = @{ $arg_ref->{tags} };

    # Convert tags to regular expressions
    my %re_for = DETCT::Misc::Tag::convert_tag_to_regexp(@tags);

    my %is_skip_sequence =
      map { $_ => 1 } @{ $arg_ref->{skip_sequences} || [] };

    my $sam = Bio::DB::Sam->new( -bam => $arg_ref->{bam_file} );

    my %stats = map { $_ => { paired => 0, mapped => 0, proper => 0, } } @tags;

    # Get all reads
    my $alignments = $sam->features( -iterator => 1, );
    while ( my $alignment = $alignments->next_seq ) {
        next
          if ( $alignment->seq_id
            && exists $is_skip_sequence{ $alignment->seq_id } )
          || ( $alignment->mate_seq_id
            && exists $is_skip_sequence{ $alignment->mate_seq_id } );

        # Match tag
        my $tag_found = matched_tag( $alignment, \%re_for );
        next if !$tag_found;

        # Counts
        if ( is_paired($alignment) ) {
            $stats{$tag_found}{paired}++;
        }
        if ( is_mapped_pair($alignment) ) {
            $stats{$tag_found}{mapped}++;
        }
        if ( is_properly_paired($alignment) ) {
            $stats{$tag_found}{proper}++;
        }
    }

    return \%stats;
}

=func stats_all_reads

  Usage       : my $stats = DETCT::Misc::BAM::stats_all_reads( {
                    bam_file => $bam_file,
                } );
  Purpose     : Generate stats for all reads in a BAM file
  Returns     : Hashref {
                    paired => Int (paired read count),
                    mapped => Int (mapped paired read count),
                    proper => Int (properly paired read count),
                }
  Parameters  : Hashref {
                    bam_file       => String (the BAM file)
                    skip_sequences => Arrayref of strings (the skip sequences)
                }
  Throws      : If BAM file is missing
  Comments    : None

=cut

sub stats_all_reads {
    my ($arg_ref) = @_;

    confess 'No BAM file specified' if !defined $arg_ref->{bam_file};

    my %is_skip_sequence =
      map { $_ => 1 } @{ $arg_ref->{skip_sequences} || [] };

    my $sam = Bio::DB::Sam->new( -bam => $arg_ref->{bam_file} );

    my %stats = ( paired => 0, mapped => 0, proper => 0, );

    # Get all reads
    my $alignments = $sam->features( -iterator => 1, );
    while ( my $alignment = $alignments->next_seq ) {
        next
          if ( $alignment->seq_id
            && exists $is_skip_sequence{ $alignment->seq_id } )
          || ( $alignment->mate_seq_id
            && exists $is_skip_sequence{ $alignment->mate_seq_id } );

        # Counts
        if ( is_paired($alignment) ) {
            $stats{paired}++;
        }
        if ( is_mapped_pair($alignment) ) {
            $stats{mapped}++;
        }
        if ( is_properly_paired($alignment) ) {
            $stats{proper}++;
        }
    }

    return \%stats;
}

=func downsample_by_tag

  Usage       : DETCT::Misc::BAM::downsample_by_tag( {
                    source_bam_file   => $bam_file,
                    source_read_count => 1_234_234,
                    tag               => 'NNNNBGAGGC',
                    target_bam_file   => $new_bam_file,
                    target_read_count => 1_000_000,
                    read_count_type   => 'proper',
                } );
  Purpose     : Downsample one tag in a BAM file to a target read count of
                chosen type
  Returns     : +ve Int (the number of reads in the target BAM file)
  Parameters  : Hashref {
                    source_bam_file   => String (the original BAM file)
                    source_read_count => Int (the original read count)
                    tag               => String (the tag)
                    target_bam_file   => String (the downsampled BAM file)
                    target_read_count => Int (the target read count)
                    read_count_type   => String ('paired', 'mapped' or 'proper')
                    skip_sequences    => Arrayref of strings (skip sequences)
                }
  Throws      : If tag is missing
  Comments    : None

=cut

sub downsample_by_tag {
    my ($arg_ref) = @_;

    confess 'No tag specified' if !defined $arg_ref->{tag};

    return downsample($arg_ref);
}

=func downsample_all_reads

  Usage       : DETCT::Misc::BAM::downsample_all_reads( {
                    source_bam_file   => $bam_file,
                    source_read_count => 1_234_234,
                    target_bam_file   => $new_bam_file,
                    target_read_count => 1_000_000,
                    read_count_type   => 'proper',
                } );
  Purpose     : Downsample all reads in a BAM file to a target read count of
                chosen type
  Returns     : +ve Int (the number of reads in the target BAM file)
  Parameters  : Hashref {
                    source_bam_file   => String (the original BAM file)
                    source_read_count => Int (the original read count)
                    target_bam_file   => String (the downsampled BAM file)
                    target_read_count => Int (the target read count)
                    read_count_type   => String ('paired', 'mapped' or 'proper')
                    skip_sequences    => Arrayref of strings (skip sequences)
                }
  Throws      : No exceptions
  Comments    : None

=cut

sub downsample_all_reads {
    my ($arg_ref) = @_;
    return downsample($arg_ref);
}

=func downsample

  Usage       : DETCT::Misc::BAM::downsample( {
                    source_bam_file   => $bam_file,
                    source_read_count => 1_234_234,
                    target_bam_file   => $new_bam_file,
                    target_read_count => 1_000_000,
                    read_count_type   => 'proper',
                } );
  Purpose     : Downsample a BAM file to a target read count of chosen type
  Returns     : +ve Int (the number of reads in the target BAM file)
  Parameters  : Hashref {
                    source_bam_file   => String (the original BAM file)
                    source_read_count => Int (the original read count)
                    tag               => String (the tag)
                    target_bam_file   => String (the downsampled BAM file)
                    target_read_count => Int (the target read count)
                    read_count_type   => String ('paired', 'mapped' or 'proper')
                    skip_sequences    => Arrayref of strings (skip sequences)
                }
  Throws      : If source BAM file is missing
                If source read count is missing
                If target BAM file is missing
                If target read count is missing
                If read count type is missing or invalid
  Comments    : None

=cut

sub downsample {    ## no critic (ProhibitExcessComplexity)
    my ($arg_ref) = @_;

    confess 'No source BAM file specified'
      if !defined $arg_ref->{source_bam_file};
    confess 'No source read count specified'
      if !defined $arg_ref->{source_read_count};
    confess 'No target BAM file specified'
      if !defined $arg_ref->{target_bam_file};
    confess 'No target read count specified'
      if !defined $arg_ref->{target_read_count};
    confess 'No read count type specified'
      if !defined $arg_ref->{read_count_type};
    if ( all { $_ ne $arg_ref->{read_count_type} } qw(paired mapped proper) ) {
        confess sprintf 'Invalid read count type (%s) specified',
          $arg_ref->{read_count_type};
    }

    # Convert tag to regular expressions
    my %re_for =
      exists $arg_ref->{tag}
      ? DETCT::Misc::Tag::convert_tag_to_regexp( $arg_ref->{tag} )
      : ();

    # Open source and target and write header to target
    my $bam_in  = Bio::DB::Bam->open( $arg_ref->{source_bam_file}, q{r} );
    my $bam_out = Bio::DB::Bam->open( $arg_ref->{target_bam_file}, q{w} );
    $bam_out->header_write( $bam_in->header );

    # Convert skip sequences to BAM refID
    my $i = 0;
    my %sam_to_bam = map { $_ => $i++ } @{ $bam_in->header->target_name };
    my %is_skip_sequence =
      map { $sam_to_bam{$_} => 1 } @{ $arg_ref->{skip_sequences} || [] };

    # Calculate probability of keeping read
    my $keep_chance =
      $arg_ref->{target_read_count} / $arg_ref->{source_read_count};

    # Cache read names whilst waiting for other of pair
    my %keep;
    my %discard;

    my $total_kept      = 0;
    my $total_half_kept = 0;

    # Get all reads
    while ( my $alignment = $bam_in->read1 ) {
        last
          if $total_kept == $arg_ref->{target_read_count} && !$total_half_kept;

        next
          if (
            $alignment->tid != -1    ## no critic (ProhibitMagicNumbers)
            && exists $is_skip_sequence{ $alignment->tid }
          )
          || (
            $alignment->mtid != -1    ## no critic (ProhibitMagicNumbers)
            && exists $is_skip_sequence{ $alignment->mtid }
          );

        # Write read if kept mate
        if ( exists $keep{ $alignment->qname } ) {
            delete $keep{ $alignment->qname };
            $bam_out->write1($alignment);
            $total_kept++;
            $total_half_kept--;
            next;
        }

        # Discard read if discarded mate
        if ( exists $discard{ $alignment->qname } ) {
            delete $discard{ $alignment->qname };
            next;
        }

        # Check if more pairs not needed (i.e. just waiting for mates)
        next
          if ( $total_kept + $total_half_kept ) >=
          $arg_ref->{target_read_count};

        # Skip reads of wrong type
        next if !is_paired($alignment);
        next
          if $arg_ref->{read_count_type} eq 'mapped'
          && !is_mapped_pair($alignment);
        next
          if $arg_ref->{read_count_type} eq 'proper'
          && !is_properly_paired($alignment);

        # Match tag (if required)
        next if exists $arg_ref->{tag} && !matched_tag( $alignment, \%re_for );

        if ( rand() < $keep_chance ) {
            $keep{ $alignment->qname } = 1;
            $bam_out->write1($alignment);
            $total_kept++;
            $total_half_kept++;
        }
        else {
            $discard{ $alignment->qname } = 1;
        }
    }

    return $total_kept;
}

=func mark_duplicates

  Usage       : my $metrics = DETCT::Misc::BAM::mark_duplicates( {
                    input_bam_file  => $bam_file,
                    output_bam_file => $new_bam_file,
                    consider_tags   => 1,
                    tags            => ['NNNNBGAGGC', 'NNNNBAGAAG'],
                } );
  Purpose     : Mark duplicates in a name-sorted BAM file
  Returns     : Hashref {
                    String (tag or pseudo-tag) => Hashref {
                        mapped_reads_without_mapped_mate           => Int,
                        mapped_read_pairs                          => Int,
                        mapped_reads                               => Int,
                        unmapped_reads                             => Int,
                        duplicate_mapped_reads_without_mapped_mate => Int,
                        duplicate_mapped_read_pairs                => Int,
                        optical_duplicate_mapped_read_pairs        => Int,
                        duplicate_reads                            => Int,
                        duplication_rate                           => Float,
                        estimated_library_size                     => Int,
                    }
                }
  Parameters  : Hashref {
                    input_bam_file  => String (the original BAM file),
                    output_bam_file => String (the output BAM file),
                    consider_tags   => Boolean (whether to consider tags),
                    tags            => Arrayref of strings (the tags),
                }
  Throws      : If input BAM file is missing
                If output BAM file is missing
                If input BAM file is not paired end, sorted by queryname
  Comments    : Input BAM file must be sorted by queryname
                If consider tags parameter is true, then tags are still optional
                (i.e. consider tags, but don't collect tag-specific metrics)

=cut

sub mark_duplicates {    ## no critic (ProhibitExcessComplexity)
    my ($arg_ref) = @_;

    confess 'No input BAM file specified'
      if !defined $arg_ref->{input_bam_file};
    confess 'No output BAM file specified'
      if !defined $arg_ref->{output_bam_file};

    # Ignore any specified tags if not considering tags
    my @tags = $arg_ref->{consider_tags}
      && $arg_ref->{tags} ? @{ $arg_ref->{tags} } : ();

    # Convert tag to regular expressions
    my %re_for = @tags ? DETCT::Misc::Tag::convert_tag_to_regexp(@tags) : ();

    # Open input (for first pass)
    my $sam_in = Bio::DB::Sam->new( -bam => $arg_ref->{input_bam_file} );

    # Track signatures of each read pair and separately of each read
    # Value is best total base quality or, in the case of reads from mapped
    # pairs, undef
    my %is_pe_dupe;
    my %is_se_dupe;

    # Get all reads in pairs where both reads are mapped
    my $alignments_pe = $sam_in->features(
        ## no critic (RequireInterpolationOfMetachars)
        -filter => 'return if $a->unmapped || $a->munmapped',
        ## use critic
        -iterator => 1,
    );
    while ( my $alignment1 = $alignments_pe->next_seq ) {
        my $alignment2 = $alignments_pe->next_seq;

        # Ensure paired end reads, sorted by read name
        confess sprintf 'Read names do not match (%s and %s) in %s',
          $alignment1->qname, $alignment2->qname, $arg_ref->{input_bam_file}
          if $alignment1->qname ne $alignment2->qname;

        # Get reference sequence, 5' end position and strand for each read
        my @signature_components;
        foreach my $alignment ( $alignment1, $alignment2 ) {
            my $pos = get_five_prime_pos_plus_soft_clip($alignment);
            push @signature_components,
              [ $alignment->tid, $pos, $alignment->strand ];
        }

        # Sort signature components, so order is consistent
        # (i.e. pair can be duplicates even if read 1 and read 2 are in
        # opposite orientations, so long as positions are same)
        @signature_components = sort _sort_sig @signature_components;

        # Make signatures
        my $signature_pe = join q{:},
          ( @{ $signature_components[0] }, @{ $signature_components[1] } );
        my $signature_se1 = join q{:}, @{ $signature_components[0] };
        my $signature_se2 = join q{:}, @{ $signature_components[1] };

        # Add tag to signatures if necessary
        if ( $arg_ref->{consider_tags} ) {
            my $tag_bytes = get_tag_for_signature($alignment1);
            $signature_pe  .= q{:} . $tag_bytes;
            $signature_se1 .= q{:} . $tag_bytes;
            $signature_se2 .= q{:} . $tag_bytes;
        }

        # Store total base quality score for read pair if best score
        my $score =
          get_base_qual_sum($alignment1) + get_base_qual_sum($alignment2);
        if ( !exists $is_pe_dupe{$signature_pe}
            || $score > $is_pe_dupe{$signature_pe} )
        {
            $is_pe_dupe{$signature_pe} = $score;
        }

        # Mark each read as coming from a read pair where both reads are mapped
        $is_se_dupe{$signature_se1} = undef;
        $is_se_dupe{$signature_se2} = undef;
    }

    # Get all reads where at least one of pair isn't mapped
    my $alignments_se = $sam_in->features(
        ## no critic (RequireInterpolationOfMetachars)
        -filter => 'return if !$a->unmapped && !$a->munmapped',
        ## use critic
        -iterator => 1,
    );
    while ( my $alignment = $alignments_se->next_seq ) {
        next if $alignment->unmapped;

        # Get 5' end position
        my $pos = get_five_prime_pos_plus_soft_clip($alignment);

        # Make signature
        my $signature_se = join q{:}, $alignment->tid, $pos, $alignment->strand;

        # Add tag to signature if necessary
        if ( $arg_ref->{consider_tags} ) {
            my $tag_bytes = get_tag_for_signature($alignment);
            $signature_se .= q{:} . $tag_bytes;
        }

        # Store total base quality score for read pair if best score and not
        # seen as part of mapped pair
        my $score = get_base_qual_sum($alignment);
        if (
            !exists $is_se_dupe{$signature_se}
            || ( defined $is_se_dupe{$signature_se}
                && $score > $is_se_dupe{$signature_se} )
          )
        {
            $is_se_dupe{$signature_se} = $score;
        }
    }

    # Make bitmask for marking non-duplicate reads
    my $not_dupe_mask =
      $DETCT::Misc::BAM::Flag::READ_PAIRED +
      $DETCT::Misc::BAM::Flag::PROPER_PAIR +
      $DETCT::Misc::BAM::Flag::READ_UNMAPPED +
      $DETCT::Misc::BAM::Flag::MATE_UNMAPPED +
      $DETCT::Misc::BAM::Flag::READ_REVERSE_STRAND +
      $DETCT::Misc::BAM::Flag::MATE_REVERSE_STRAND +
      $DETCT::Misc::BAM::Flag::FIRST_IN_PAIR +
      $DETCT::Misc::BAM::Flag::SECOND_IN_PAIR +
      $DETCT::Misc::BAM::Flag::SECONDARY +
      $DETCT::Misc::BAM::Flag::QC_FAILED +
      $DETCT::Misc::BAM::Flag::SUPPLEMENTARY;

    # Initialise metrics
    my $metrics = _init_duplication_metrics(@tags);

    # Open input (for second pass) and output and write header to output
    $sam_in = undef;
    my $bam_in  = Bio::DB::Bam->open( $arg_ref->{input_bam_file},  q{r} );
    my $bam_out = Bio::DB::Bam->open( $arg_ref->{output_bam_file}, q{w} );
    $bam_out->header_write( $bam_in->header );

    # Get all reads in pairs
    while ( my $alignment1 = $bam_in->read1 ) {
        my $alignment2 = $bam_in->read1;

        my $is_dupe          = 0;
        my $num_reads_mapped = 0;
        my $tag;

        if ( !$alignment1->unmapped && !$alignment2->unmapped ) {

            # Both reads mapped
            $num_reads_mapped = 2;

            # Get reference sequence, 5' end position and strand for each read
            my @signature_components;
            foreach my $alignment ( $alignment1, $alignment2 ) {
                my $pos = get_five_prime_pos_plus_soft_clip($alignment);
                push @signature_components,
                  [ $alignment->tid, $pos, $alignment->strand ];
            }

            # Sort signature components, so order is consistent
            # (i.e. pair can be duplicates even if read 1 and read 2 are in
            # opposite orientations, so long as positions are same)
            @signature_components = sort _sort_sig @signature_components;

            # Make signature
            my $signature_pe = join q{:},
              ( @{ $signature_components[0] }, @{ $signature_components[1] } );

            # Add tag to signature if necessary
            if ( $arg_ref->{consider_tags} ) {
                my $tag_bytes = get_tag_for_signature($alignment1);
                $signature_pe .= q{:} . $tag_bytes;
                $tag = matched_tag( $alignment1, \%re_for );
            }

            # Not a duplicate if best score
            my $score =
              get_base_qual_sum($alignment1) + get_base_qual_sum($alignment2);
            if ( defined $is_pe_dupe{$signature_pe}
                && $is_pe_dupe{$signature_pe} == $score )
            {
                # Not a duplicate, so change flags
                $alignment1->flag( $alignment1->flag & $not_dupe_mask );
                $alignment2->flag( $alignment2->flag & $not_dupe_mask );

                # All others will be dupes
                $is_pe_dupe{$signature_pe} = undef;
            }
            else {
                # Is a duplicate, so change flags
                $is_dupe = 1;
                $alignment1->flag(
                    $alignment1->flag | $DETCT::Misc::BAM::Flag::DUPLICATE );
                $alignment2->flag(
                    $alignment2->flag | $DETCT::Misc::BAM::Flag::DUPLICATE );
            }
        }
        elsif ( !( $alignment1->unmapped && $alignment2->unmapped ) ) {

            # One read mapped
            $num_reads_mapped = 1;

            foreach my $alignment ( $alignment1, $alignment2 ) {

                # If unmapped then not a duplicate
                if ( $alignment->unmapped ) {
                    $alignment->flag( $alignment->flag & $not_dupe_mask );
                    next;
                }

                # Get 5' end position
                my $pos = get_five_prime_pos_plus_soft_clip($alignment);

                # Make signature
                my $signature_se = join q{:}, $alignment->tid, $pos,
                  $alignment->strand;

                # Add tag to signature if necessary
                if ( $arg_ref->{consider_tags} ) {
                    my $tag_bytes = get_tag_for_signature($alignment);
                    $signature_se .= q{:} . $tag_bytes;
                    $tag = matched_tag( $alignment, \%re_for );
                }

                # Not a duplicate if best score
                my $score = get_base_qual_sum($alignment);
                if ( defined $is_se_dupe{$signature_se}
                    && $is_se_dupe{$signature_se} == $score )
                {
                    # Not a duplicate, so change flags
                    $alignment->flag( $alignment->flag & $not_dupe_mask );

                    # All others will be dupes
                    $is_se_dupe{$signature_se} = undef;
                }
                else {
                    # Is a duplicate, so change flags
                    $is_dupe = 1;
                    $alignment->flag(
                        $alignment->flag | $DETCT::Misc::BAM::Flag::DUPLICATE );
                }
            }
        }
        else {
            # Neither read mapped, but get tag if necessary
            if ( $arg_ref->{consider_tags} ) {
                $tag = matched_tag( $alignment1, \%re_for );
            }
            $alignment1->flag( $alignment1->flag & $not_dupe_mask );
            $alignment2->flag( $alignment2->flag & $not_dupe_mask );
        }

        # Update metrics
        $metrics = _update_duplication_metrics(
            {
                metrics          => $metrics,
                tags             => \@tags,
                tag              => $tag,
                num_reads_mapped => $num_reads_mapped,
                is_dupe          => $is_dupe,
            }
        );

        # Write reads
        $bam_out->write1($alignment1);
        $bam_out->write1($alignment2);
    }

    # Calculate derived metrics
    $metrics = _calc_derived_duplication_metrics($metrics);

    return $metrics;
}

# Sort signature components
sub _sort_sig {
    return $a->[0] cmp $b->[0]    # Reference sequence
      || $a->[1] cmp $b->[1]      # 5' end position
      || $a->[2] cmp $b->[2]      # Strand
}

=func get_five_prime_pos_plus_soft_clip

  Usage       : my $pos = get_five_prime_pos_plus_soft_clip( $alignment );
  Purpose     : Get 5' end position adjusted for soft-clipped bases
  Returns     : Int (the adjusted 5' position)
  Parameters  : Bio::DB::Bam::Alignment or Bio::DB::Bam::AlignWrapper
  Throws      : No exceptions
  Comments    : None

=cut

sub get_five_prime_pos_plus_soft_clip {
    my ($alignment) = @_;

    my $pos;

    my $cigar_ref = $alignment->cigar_array;
    if ( $alignment->strand == 1 ) {
        $pos = $alignment->start;
        my $pair_ref = shift @{$cigar_ref};
        my ( $op, $count ) = @{$pair_ref};
        if ( $op eq q{S} ) {
            $pos -= $count;
        }
    }
    else {
        $pos = $alignment->end;
        my $pair_ref = pop @{$cigar_ref};
        my ( $op, $count ) = @{$pair_ref};
        if ( $op eq q{S} ) {
            $pos += $count;
        }
    }

    return $pos;
}

=func get_base_qual_sum

  Usage       : my $score = get_base_qual_sum( $alignment );
  Purpose     : Get total base quality (15+) for a read
  Returns     : Int (the total base quality)
  Parameters  : Bio::DB::Bam::Alignment or Bio::DB::Bam::AlignWrapper
  Throws      : No exceptions
  Comments    : None

=cut

sub get_base_qual_sum {
    my ($alignment) = @_;

    ## no critic (ProhibitMagicNumbers)
    return sum( grep { $_ >= 15 } $alignment->qscore ) || 0;
    ## use critic
}

=func get_tag_for_signature

  Usage       : my $tag_bytes = get_tag_for_signature( $alignment );
  Purpose     : Get tag as bytes
  Returns     : String (the tag as bytes)
  Parameters  : Bio::DB::Bam::Alignment or Bio::DB::Bam::AlignWrapper
  Throws      : No exceptions
  Comments    : None

=cut

sub get_tag_for_signature {
    my ($alignment) = @_;

    my ($tag) = $alignment->qname =~ m/[#] ([NAGCTX]+) \z/xmsg;

    # Convert tag to binary
    # N = 000; A = 001; G = 010; C = 011; T = 100; X = 101
    $tag =~ s/N/000/xmsg;
    $tag =~ s/A/001/xmsg;
    $tag =~ s/G/010/xmsg;
    $tag =~ s/C/011/xmsg;
    $tag =~ s/T/100/xmsg;
    $tag =~ s/X/101/xmsg;
    $tag = pack 'b*', $tag;

    return $tag;
}

# Usage       : my $metrics = _init_duplication_metrics( @tags );
# Purpose     : Initial duplication metrics
# Returns     : Hashref {
#                   String (tag or pseudo-tag) => Hashref {
#                       mapped_reads_without_mapped_mate           => 0,
#                       mapped_read_pairs                          => 0,
#                       unmapped_reads                             => 0,
#                       duplicate_mapped_reads_without_mapped_mate => 0,
#                       duplicate_mapped_read_pairs                => 0,
#                       optical_duplicate_mapped_read_pairs        => 0,
#                   }
#               }
# Parameters  : Array of strings (the tags)
# Throws      : No exceptions
# Comments    : None
sub _init_duplication_metrics {
    my (@tags) = @_;

    my $metrics;

    my @tags_and_pseudotags = ('_all');
    if (@tags) {
        push @tags_and_pseudotags, @tags, '_other';
    }
    foreach my $tag (@tags_and_pseudotags) {
        $metrics->{$tag}{mapped_reads_without_mapped_mate}           = 0;
        $metrics->{$tag}{mapped_read_pairs}                          = 0;
        $metrics->{$tag}{unmapped_reads}                             = 0;
        $metrics->{$tag}{duplicate_mapped_reads_without_mapped_mate} = 0;
        $metrics->{$tag}{duplicate_mapped_read_pairs}                = 0;
        $metrics->{$tag}{optical_duplicate_mapped_read_pairs}        = 0;
    }

    return $metrics;
}

# Usage       : $metrics = _update_duplication_metrics( {
#                   metrics          => $metrics,
#                   tags             => ['NNNNBGAGGC', 'NNNNBAGAAG'],
#                   tag              => 'NNNNBGAGGC',
#                   num_reads_mapped => 1,
#                   is_dupe          => 1,
#               } );
# Purpose     : Update duplication metrics
# Returns     : Hashref {
#                   String (tag or pseudo-tag) => Hashref {
#                       mapped_reads_without_mapped_mate           => Int,
#                       mapped_read_pairs                          => Int,
#                       unmapped_reads                             => Int,
#                       duplicate_mapped_reads_without_mapped_mate => Int,
#                       duplicate_mapped_read_pairs                => Int,
#                       optical_duplicate_mapped_sread_pairs       => Int,
#                   }
#               }
# Parameters  : Hashref {
#                   metrics          => Hashref (the metrics),
#                   tags             => Arrayref of strings (the tags),
#                   tag              => String (the tag) or undef,
#                   num_reads_mapped => Int (the number of reads mapped),
#                   is_dupe          => Boolean (whether reads are duplicates),
#               }
# Throws      : No exceptions
# Comments    : Metrics are updated for a single read pair
sub _update_duplication_metrics {
    my ($arg_ref) = @_;

    my $metrics = $arg_ref->{metrics};

    my @tags_and_pseudotags = ('_all');
    if ( @{ $arg_ref->{tags} } && $arg_ref->{tag} ) {
        push @tags_and_pseudotags, $arg_ref->{tag};
    }
    elsif ( @{ $arg_ref->{tags} } ) {
        push @tags_and_pseudotags, '_other';
    }
    foreach my $tag (@tags_and_pseudotags) {
        if ( $arg_ref->{num_reads_mapped} == 0 ) {

            # Neither read mapped
            $metrics->{$tag}{unmapped_reads} += 2;
        }
        elsif ( $arg_ref->{num_reads_mapped} == 1 ) {

            # One read of pair mapped
            $metrics->{$tag}{mapped_reads_without_mapped_mate}++;
            $metrics->{$tag}{unmapped_reads}++;
            if ( $arg_ref->{is_dupe} ) {
                $metrics->{$tag}{duplicate_mapped_reads_without_mapped_mate}++;
            }
        }
        else {
            # Both reads of pair mapped
            $metrics->{$tag}{mapped_read_pairs}++;
            if ( $arg_ref->{is_dupe} ) {
                $metrics->{$tag}{duplicate_mapped_read_pairs}++;
            }
        }
    }

    return $metrics;
}

# Usage       : $metrics = _calc_derived_duplication_metrics( $metrics );
# Purpose     : Calculate metrics derived from final metrics
# Returns     : Hashref {
#                   String (tag or pseudo-tag) => Hashref {
#                       mapped_reads_without_mapped_mate           => Int,
#                       mapped_read_pairs                          => Int,
#                       mapped_reads                               => Int,
#                       unmapped_reads                             => Int,
#                       duplicate_mapped_reads_without_mapped_mate => Int,
#                       duplicate_mapped_read_pairs                => Int,
#                       optical_duplicate_mapped_read_pairs        => Int,
#                       duplicate_reads                            => Int,
#                       duplication_rate                           => Float,
#                       estimated_library_size                     => Int,
#                   }
#               }
# Parameters  : Hashref {
#                   String (tag or pseudo-tag) => Hashref {
#                       mapped_reads_without_mapped_mate           => Int,
#                       mapped_read_pairs                          => Int,
#                       unmapped_reads                             => Int,
#                       duplicate_mapped_reads_without_mapped_mate => Int,
#                       duplicate_mapped_read_pairs                => Int,
#                       optical_duplicate_mapped_read_pairs        => Int,
#                   }
#               }
# Throws      : No exceptions
# Comments    : None
sub _calc_derived_duplication_metrics {
    my ($metrics) = @_;

    # Calculate derived metrics for each tag or pseudotag (_all or _other)
    foreach my $tag ( keys %{$metrics} ) {

        # Mapped reads is sum of read pairs and reads without mapped mates
        $metrics->{$tag}{mapped_reads} =
          $metrics->{$tag}{mapped_reads_without_mapped_mate} +
          $metrics->{$tag}{mapped_read_pairs} * 2;

        # Duplicate reads is sum of duplicate read pairs and duplicate reads
        # without mapped mates
        $metrics->{$tag}{duplicate_reads} =
          $metrics->{$tag}{duplicate_mapped_reads_without_mapped_mate} +
          $metrics->{$tag}{duplicate_mapped_read_pairs} * 2;

        # Duplication rate is duplicated reads divided by mapped reads (or 0 if
        # no mapped reads)
        if ( $metrics->{$tag}{mapped_reads} ) {
            $metrics->{$tag}{duplication_rate} =
              $metrics->{$tag}{duplicate_reads} /
              $metrics->{$tag}{mapped_reads};
        }
        else {
            $metrics->{$tag}{duplication_rate} = 0;
        }

        # Estimated library size (0 if no mapped read pairs)
        if ( $metrics->{$tag}{mapped_read_pairs} ) {
            $metrics->{$tag}{estimated_library_size} = _estimate_library_size(
                {
                    num_read_pairs => $metrics->{$tag}{mapped_read_pairs},
                    num_duplicate_read_pairs =>
                      $metrics->{$tag}{duplicate_mapped_read_pairs},
                }
            ) || 0;
        }
        else {
            $metrics->{$tag}{estimated_library_size} = 0;
        }
    }

    return $metrics;
}

# Usage       : my $size = _estimate_library_size( {
#                   num_read_pairs           => 100000,
#                   num_duplicate_read_pairs => 1000,
#               } );
# Purpose     : Estimate the size of a library
# Returns     : Int (the estimated library size) or undef
# Parameters  : Hashref {
#                   num_read_pairs           => Int (read pair count),
#                   num_duplicate_read_pairs => Int (duplicate read pair count),
#               }
# Throws      : If read pair count or duplicate read pair count are not positive
#               integers or duplicate read pair count is higher than read pair
#               count
# Comments    : Based on code from Picard's DuplicationMetrics.java
#               Uses the Lander-Waterman equation:
#                   C/X = 1 - exp( -N/X )
#               Where:
#                   X = number of distinct molecules in library
#                   N = number of read pairs
#                   C = number of unique read pairs
sub _estimate_library_size {
    my ($arg_ref) = @_;

    # n = number of read pairs
    my $n = $arg_ref->{num_read_pairs};

    # c = number of unique read pairs
    my $c = $arg_ref->{num_read_pairs} - $arg_ref->{num_duplicate_read_pairs};

    if ( $n < 1 || $c < 1 || $c > $n ) {
        confess sprintf 'Invalid read pairs (%d) or unique read pairs (%d)',
          $n, $c;
    }

    # Can't estimate library size if all unique
    if ( $c == $n ) {
        return;
    }

    my $lo = 1;
    my $hi = 10;    ## no critic (ProhibitMagicNumbers)

    # Make upper limit high enough
    while ( ( 1 / $hi - 1 + exp -$n / $hi / $c ) >= 0 ) {
        $hi *= 10;    ## no critic (ProhibitMagicNumbers)
    }

    # Converge on estimate of library size
    foreach ( 1 .. 40 ) {    ## no critic (ProhibitMagicNumbers)
        my $avg  = ( $lo + $hi ) / 2;
        my $diff = 1 / $avg - 1 + exp -$n / $avg / $c;
        if ( $diff > 0 ) {
            $lo = $avg;
        }
        elsif ( $diff < 0 ) {
            $hi = $avg;
        }
    }

    return int( $c * ( $lo + $hi ) / 2 );
}

=func filter_by_tag

  Usage       : DETCT::Misc::BAM::filter_by_tag( {
                    source_bam_file => $bam_file,
                    target_bam_file => $new_bam_file,
                    tags            => ['NNNNBGAGGC', 'NNNNBAGAAG'],
                } );
  Purpose     : Filter a BAM file to only include specific tags
  Returns     : undef
  Parameters  : Hashref {
                    source_bam_file => String (the original BAM file)
                    target_bam_file => String (the filtered BAM file)
                    tags            => Arrayref of strings (the tags)
                }
  Throws      : If source BAM file is missing
                If target BAM file is missing
                If tags are missing
  Comments    : None

=cut

sub filter_by_tag {
    my ($arg_ref) = @_;

    confess 'No source BAM file specified'
      if !defined $arg_ref->{source_bam_file};
    confess 'No target BAM file specified'
      if !defined $arg_ref->{target_bam_file};
    confess 'No tags specified' if !defined $arg_ref->{tags};

    my @tags = @{ $arg_ref->{tags} };

    # Convert tags to regular expressions
    my %re_for = DETCT::Misc::Tag::convert_tag_to_regexp(@tags);

    # Open source and target and write header to target
    my $bam_in  = Bio::DB::Bam->open( $arg_ref->{source_bam_file}, q{r} );
    my $bam_out = Bio::DB::Bam->open( $arg_ref->{target_bam_file}, q{w} );
    $bam_out->header_write( $bam_in->header );

    # Get all reads
    while ( my $alignment = $bam_in->read1 ) {
        next if !matched_tag( $alignment, \%re_for );
        $bam_out->write1($alignment);
    }

    return;
}

=func matched_tag

  Usage       : next if !matched_tag($alignment, \%re_for);
  Purpose     : Get tag matching alignment
  Returns     : String (the matched tag) or undef
  Parameters  : Bio::DB::Bam::Alignment or Bio::DB::Bam::AlignWrapper
              : Hashref of regular expressions
  Throws      : No exceptions
  Comments    : None

=cut

sub matched_tag {
    my ( $alignment, $re_for ) = @_;

    # Match tag
    my ($tag_in_read) = $alignment->query->name =~ m/[#] ([NAGCT]+) \z/xmsg;
    if ($tag_in_read) {
        foreach my $tag ( sort keys %{$re_for} ) {
            my $regexps = $re_for->{$tag};
            foreach my $re ( @{$regexps} ) {
                if ( $tag_in_read =~ $re ) {
                    return $tag;
                }
            }
        }
    }

    return;
}

=func is_read2

  Usage       : next if is_read2($alignment);
  Purpose     : Check if alignment is from read 2 (not read 1)
  Returns     : 1 or 0
  Parameters  : Bio::DB::Bam::AlignWrapper
  Throws      : No exceptions
  Comments    : None

=cut

sub is_read2 {
    my ($alignment) = @_;

    return ( $alignment->get_tag_values('FLAGS') =~ m/\bSECOND_MATE\b/xms )
      ? 1
      : 0;
}

=func is_duplicate

  Usage       : next if is_duplicate($alignment);
  Purpose     : Check if alignment is marked as a duplicate
  Returns     : 1 or 0
  Parameters  : Bio::DB::Bam::AlignWrapper
  Throws      : No exceptions
  Comments    : None

=cut

sub is_duplicate {
    my ($alignment) = @_;

    return ( $alignment->get_tag_values('FLAGS') =~ m/\bDUPLICATE\b/xms )
      ? 1
      : 0;
}

=func is_paired

  Usage       : next if is_paired($alignment);
  Purpose     : Check if alignment represents a read paired in sequencing
  Returns     : 1 or 0
  Parameters  : Bio::DB::Bam::AlignWrapper
  Throws      : No exceptions
  Comments    : None

=cut

sub is_paired {
    my ($alignment) = @_;

    return ( $alignment->get_tag_values('FLAGS') =~ m/\bPAIRED\b/xms ) ? 1 : 0;
}

=func is_mapped_pair

  Usage       : next if is_mapped_pair($alignment);
  Purpose     : Check if alignment represents a pair of mapped reads
  Returns     : 1 or 0
  Parameters  : Bio::DB::Bam::AlignWrapper
  Throws      : No exceptions
  Comments    : None

=cut

sub is_mapped_pair {
    my ($alignment) = @_;

    return ( $alignment->unmapped || $alignment->munmapped ) ? 0 : 1;
}

=func is_properly_paired

  Usage       : next if is_properly_paired($alignment);
  Purpose     : Check if alignment represents a properly paired read
  Returns     : 1 or 0
  Parameters  : Bio::DB::Bam::AlignWrapper
  Throws      : No exceptions
  Comments    : None

=cut

sub is_properly_paired {
    my ($alignment) = @_;

    return ( $alignment->get_tag_values('FLAGS') =~ m/\bMAP_PAIR\b/xms )
      ? 1
      : 0;
}

=func above_mismatch_threshold

  Usage       : next if above_mismatch_threshold($alignment, 2);
  Purpose     : Check if alignment has too many mismatches
  Returns     : 1 or 0
  Parameters  : Bio::DB::Bam::Alignment or Bio::DB::Bam::AlignWrapper
              : Int (mismatch threshold)
  Throws      : No exceptions
  Comments    : None

=cut

sub above_mismatch_threshold {
    my ( $alignment, $threshold ) = @_;

    # Count soft clipped bases
    my $cigar_ref          = $alignment->cigar_array;
    my $soft_clipped_bases = 0;
    foreach my $pair_ref ( @{$cigar_ref} ) {
        my ( $op, $count ) = @{$pair_ref};
        if ( $op eq q{S} ) {
            $soft_clipped_bases += $count;
        }
    }

    # Get edit distance / number of mismatches
    my $nm = $alignment->aux_get('NM');

    # Check if above mismatch threshold
    return ( $nm + $soft_clipped_bases > $threshold ) ? 1 : 0;
}

=func is_polya

  Usage       : next if is_polya($seq);
  Purpose     : Check if sequence contains polyA
  Returns     : 1 or 0
  Parameters  : String (sequence)
  Throws      : No exceptions
  Comments    : None

=cut

sub is_polya {
    my ($seq) = @_;

    my $is_polya = 0;

    # Check for more than 3 As at start
    if ( $seq =~ m/\A AAAA /xms ) {
        $is_polya = 1;
    }

    # Check for more than 6 As in total
    if ( !$is_polya ) {
        my $a = $seq =~ tr/A/A/;
        if ( $a > 6 ) {    ## no critic (ProhibitMagicNumbers)
            $is_polya = 1;
        }
    }

    # Check specific patterns for polyA
    if ( !$is_polya ) {
        foreach my $regexp (@POLYA_REGEXP) {
            if ( $seq =~ $regexp ) {
                $is_polya = 1;
                last;
            }
        }
    }

    return $is_polya;
}

1;
