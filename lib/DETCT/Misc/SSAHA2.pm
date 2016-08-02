## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::SSAHA2;
## VERSION
## use critic

# ABSTRACT: Miscellaneous functions wrapping SSAHA2

## Author         : is1
## Maintainer     : is1
## Created        : 2016-08-02
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
use Path::Tiny;
use DETCT::Misc qw( write_or_die );

use base qw( Exporter );
our @EXPORT_OK = qw(
  add_region_alignments
  run_ssaha2
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func add_region_alignments

  Usage       : my $regions_ref = add_region_alignments( {
                    regions       => $regions_hash_ref,
                    dir           => '.',
                    ref_fasta     => 'zv9.fa',
                    analysis      => $analysis,
                    ssaha2_binary => 'ssaha2',
                } );
  Purpose     : Add region alignments using SSAHA2
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
                                    Boolean (whether polyA),
                                    String (upstream 14 bp),
                                    String (downstream 14 bp),
                                    Int (distance hexamer upstream) or undef,
                                    String (hexamer sequence),
                                    Int (distance to nearest transposon),
                                    Int (position of nearest transposon),
                                    Arrayref [
                                        String (transcript stable id),
                                        ... (continuous RNA-Seq transcript ids)
                                    ]
                                ],
                                ... (3' ends)
                            ],
                            Arrayref [
                                Arrayref [
                                    String (target sequence name),
                                    Int (target start),
                                    Int (target end),
                                    String (target strand),
                                    Int (percent ID),
                                    Int (alignment length),
                                    Int (query length)
                                ],
                                ... (alignments)
                            ],
                        ],
                        ... (regions)
                    ]
                }
  Parameters  : Hashref {
                    regions       => Hashref (of arrayrefs of regions),
                    dir           => String (the working directory),
                    ref_fasta     => String (the reference FASTA file),
                    analysis      => DETCT::Analysis,
                    ssaha2_binary => String (the SSAHA2 binary)
                }
  Throws      : If regions are missing
                If directory is missing
                If reference FASTA is missing
                If analysis is missing
                If SSAHA2 binary is missing
  Comments    : None

=cut

sub add_region_alignments {
    my ($arg_ref) = @_;

    confess 'No regions specified'         if !defined $arg_ref->{regions};
    confess 'No directory specified'       if !defined $arg_ref->{dir};
    confess 'No reference FASTA specified' if !defined $arg_ref->{ref_fasta};
    confess 'No analysis specified'        if !defined $arg_ref->{analysis};
    confess 'No SSAHA2 binary specified' if !defined $arg_ref->{ssaha2_binary};

    my $regions = $arg_ref->{regions};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Write regions as FASTA
    my $fasta_file = File::Spec->catfile( $arg_ref->{dir}, 'ssaha2.fa' );
    open my $fh, '>', $fasta_file;    ## no critic (RequireBriefOpen)
    foreach my $seq_name ( sort keys %{$regions} ) {
        my $i = 0;
        foreach my $region ( @{ $regions->{$seq_name} } ) {
            my ( $region_start, $region_end, undef, undef, $strand ) =
              @{$region};
            my $seq =
              $arg_ref->{analysis}
              ->get_subsequence( $seq_name, $region_start, $region_end,
                $strand );
            write_or_die( $fh, sprintf ">%s:%d\n%s\n", $seq_name, $i, $seq );
            $i++;
        }
    }
    close $fh;

    run_ssaha2(
        {
            dir           => $arg_ref->{dir},
            ref_fasta     => $arg_ref->{ref_fasta},
            fasta_file    => $fasta_file,
            ssaha2_binary => $arg_ref->{ssaha2_binary},
        }
    );

    my %alignments = parse_ssaha2( $arg_ref->{dir} );
    foreach my $seq_name ( sort keys %{$regions} ) {
        my $i = 0;
        foreach my $region ( @{ $regions->{$seq_name} } ) {
            if ( exists $alignments{$seq_name}{$i} ) {
                push @{ $regions->{$seq_name}->[$i] },
                  $alignments{$seq_name}{$i};
            }
            $i++;
        }
    }

    return $regions;
}

=func run_ssaha2

  Usage       : DETCT::Misc::SSAHA2::run_ssaha2( {
                    dir           => '.',
                    ref_fasta     => 'zv9.fa',
                    fasta_file    => 'input.fa',
                    ssaha2_binary => 'ssaha2',
                } );
  Purpose     : Run SSAHA2
  Returns     : undef
  Parameters  : Hashref {
                    dir           => String (the working directory),
                    ref_fasta     => String (the reference FASTA file),
                    fasta_file    => String (the FASTA file),
                    ssaha2_binary => String (the SSAHA2 binary),
                }
  Throws      : If directory is missing
                If reference FASTA is missing
                If FASTA file is missing
                If SSAHA2 binary is missing
                If command line can't be run
  Comments    : None

=cut

sub run_ssaha2 {
    my ($arg_ref) = @_;

    confess 'No directory specified'       if !defined $arg_ref->{dir};
    confess 'No reference FASTA specified' if !defined $arg_ref->{ref_fasta};
    confess 'No FASTA file specified'      if !defined $arg_ref->{fasta_file};
    confess 'No SSAHA2 binary specified' if !defined $arg_ref->{ssaha2_binary};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'ssaha2.o' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'ssaha2.e' );

    # If no regions then don't run SSAHA2
    if ( -z $arg_ref->{fasta_file} ) {
        path($stdout_file)->spew();
        return;
    }

    my $cmd = join q{ }, $arg_ref->{ssaha2_binary}, '-seeds', '1', '-kmer',
      '13', '-score', '12', '-skip', '1', '-cut', '50000', '-save',
      $arg_ref->{ref_fasta}, $arg_ref->{fasta_file};
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    if ( -s $stderr_file ) {
        confess "Couldn't run $cmd";
    }

    return;
}

=func parse_ssaha2

  Usage       : DETCT::Misc::SSAHA2::parse_ssaha2($dir);
  Purpose     : Parse SSAHA2 output
  Returns     : Arrayref [
                    Arrayref [
                        String (target sequence name),
                        Int (target start),
                        Int (target end),
                        String (target strand),
                        Int (percent ID),
                        Int (alignment length),
                        Int (query length)
                    ],
                    ... (alignments)
                ]
  Parameters  : String (the working directory)
  Throws      : If directory is missing
  Comments    : None

=cut

sub parse_ssaha2 {
    my ($dir) = @_;

    confess 'No directory specified' if !defined $dir;

    my %alignments;

    my $stdout_file = File::Spec->catfile( $dir, 'ssaha2.o' );
    my $query_length;
    open my $fh, '<', $stdout_file;    ## no critic (RequireBriefOpen)
    while ( my $line = <$fh> ) {
        if ( $line =~ m/\A Matches \s For \s Query \s \d+ \s [(] (\d+) \s /xms )
        {
            $query_length = $1;
        }
        next if $line !~ m/\A ALIGNMENT /xms;
        chomp $line;
        my (
            undef,             $score,       $query_id,
            $target_id,        $query_start, $query_end,
            $target_start,     $target_end,  $strand,
            $alignment_length, $percent_id
        ) = split /\s+/xms, $line;
        my ( $seq_name, $index ) = split /:/xms, $query_id;
        $strand = $strand eq q{F} ? q{+} : q{-};
        $percent_id = int $percent_id + 0.5; ## no critic (ProhibitMagicNumbers)
        push @{ $alignments{$seq_name}{$index} },
          [
            $target_id, $target_start, $target_end,
            $strand,    $percent_id,   $alignment_length,
            $query_length
          ];
    }
    close $fh;

    return %alignments;
}

1;
