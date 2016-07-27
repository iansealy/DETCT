## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::Exonerate;
## VERSION
## use critic

# ABSTRACT: Miscellaneous functions wrapping Exonerate

## Author         : is1
## Maintainer     : is1
## Created        : 2016-07-15
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
  run_exonerate
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=func add_region_alignments

  Usage       : my $regions_ref = add_region_alignments( {
                    regions          => $regions_hash_ref,
                    dir              => '.',
                    ref_fasta        => 'zv9.fa',
                    analysis         => $analysis,
                    exonerate_binary => 'exonerate',
                } );
  Purpose     : Add region alignments using Exonerate
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
                    regions          => Hashref (of arrayrefs of regions),
                    dir              => String (the working directory),
                    ref_fasta        => String (the reference FASTA file),
                    analysis         => DETCT::Analysis,
                    exonerate_binary => String (the Exonerate binary)
                }
  Throws      : If regions are missing
                If directory is missing
                If reference FASTA is missing
                If analysis is missing
                If Exonerate binary is missing
  Comments    : None

=cut

sub add_region_alignments {
    my ($arg_ref) = @_;

    confess 'No regions specified'         if !defined $arg_ref->{regions};
    confess 'No directory specified'       if !defined $arg_ref->{dir};
    confess 'No reference FASTA specified' if !defined $arg_ref->{ref_fasta};
    confess 'No analysis specified'        if !defined $arg_ref->{analysis};
    confess 'No Exonerate binary specified'
      if !defined $arg_ref->{exonerate_binary};

    my $regions = $arg_ref->{regions};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Write regions as FASTA
    my $fasta_file = File::Spec->catfile( $arg_ref->{dir}, 'exonerate.fa' );
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

    run_exonerate(
        {
            dir              => $arg_ref->{dir},
            ref_fasta        => $arg_ref->{ref_fasta},
            fasta_file       => $fasta_file,
            exonerate_binary => $arg_ref->{exonerate_binary},
        }
    );

    my %alignments = parse_exonerate( $arg_ref->{dir} );
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

=func run_exonerate

  Usage       : DETCT::Misc::Exonerate::run_exonerate( {
                    dir              => '.',
                    ref_fasta        => 'zv9.fa',
                    fasta_file       => 'input.fa',
                    exonerate_binary => 'exonerate',
                } );
  Purpose     : Run Exonerate
  Returns     : undef
  Parameters  : Hashref {
                    dir              => String (the working directory),
                    ref_fasta        => String (the reference FASTA file),
                    fasta_file       => String (the FASTA file),
                    exonerate_binary => String (the Exonerate binary),
                }
  Throws      : If directory is missing
                If reference FASTA is missing
                If FASTA file is missing
                If Exonerate binary is missing
                If command line can't be run
  Comments    : None

=cut

sub run_exonerate {
    my ($arg_ref) = @_;

    confess 'No directory specified'       if !defined $arg_ref->{dir};
    confess 'No reference FASTA specified' if !defined $arg_ref->{ref_fasta};
    confess 'No FASTA file specified'      if !defined $arg_ref->{fasta_file};
    confess 'No Exonerate binary specified'
      if !defined $arg_ref->{exonerate_binary};

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    my $stdout_file = File::Spec->catfile( $arg_ref->{dir}, 'exonerate.o' );
    my $stderr_file = File::Spec->catfile( $arg_ref->{dir}, 'exonerate.e' );

    ## no critic (RequireInterpolationOfMetachars)
    my $cmd = join q{ }, $arg_ref->{exonerate_binary}, '--query',
      $arg_ref->{fasta_file}, '--target', $arg_ref->{ref_fasta}, '--querytype',
      'dna', '--targettype', 'dna', '--model', 'affine:local', '--percent',
      '50', '--showalignment', 'FALSE', '--showvulgar', 'FALSE', '--ryo',
      q{'RESULT: %S %pi %ql %tl\n'};
    ## use critic
    $cmd .= ' 1>' . $stdout_file;
    $cmd .= ' 2>' . $stderr_file;
    WIFEXITED( system $cmd) or confess "Couldn't run $cmd ($OS_ERROR)";

    return;
}

=func parse_exonerate

  Usage       : DETCT::Misc::Exonerate::parse_exonerate($dir);
  Purpose     : Parse Exonerate output
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

sub parse_exonerate {
    my ($dir) = @_;

    confess 'No directory specified' if !defined $dir;

    my %alignments;

    my $stdout_file = File::Spec->catfile( $dir, 'exonerate.o' );
    open my $fh, '<', $stdout_file;    ## no critic (RequireBriefOpen)
    while ( my $line = <$fh> ) {
        next if $line !~ m/\A RESULT: /xms;
        chomp $line;
        my (
            undef,          $query_id,  $query_start,  $query_end,
            $query_strand,  $target_id, $target_start, $target_end,
            $target_strand, $score,     $percent_id,   $query_length,
            $target_length
        ) = split /\s/xms, $line;
        my ( $seq_name, $index ) = split /:/xms, $query_id;
        push @{ $alignments{$seq_name}{$index} },
          [
            $target_id,  $target_start,
            $target_end, $target_strand,
            $percent_id, abs( $query_end - $query_start ),
            $query_length
          ];
    }
    close $fh;

    return %alignments;
}

1;
