#!/usr/bin/env perl

# PODNAME: make_test_sam.pl
# ABSTRACT: Make transcript counting test file in SAM format

## Author         : is1
## Maintainer     : is1
## Created        : 2012-09-14
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Getopt::Long;
use Pod::Usage;
use Readonly;
use DETCT::Misc qw( print_or_die );
use DETCT::Misc::BAM::Flag;

=head1 DESCRIPTION

This script generates test transcript counting SAM files. The number and maximum
length of chromosomes can be varied along with the number and length of reads.
Read tags must be specified.

=head1 EXAMPLES

    # Generate random BAM file using default values
    perl -Ilib script/make_test_sam.pl --read_tags NNNNCTACCA \
        | samtools view -bS - | samtools sort - test

    # Generate BAM file with reproducible chromosomes using default values
    perl -Ilib script/make_test_sam.pl --seed 1 --read_tags NNNNCTACCA \
        | samtools view -bS - | samtools sort - test

    # Generate BAM file with 25 chromosomes (each up to 50 Mbp long), 1000
    # alignments per chromosome and four 10mer tags
    perl -Ilib script/make_test_sam.pl \
        --seq_region_count 25 \
        --seq_region_max_length 50_000_000 \
        --read_pair_count 1000 \
        --read_tags NNNNCTACCA NNNNAAGTTA NNNNTTAATC NNNNTAGACA \
        | samtools view -bS - | samtools sort - test

=cut

# Constants from http://samtools.sourceforge.net/SAM1.pdf

# Regexps for checking alignment line mandatory fields
Readonly our %ALIGNMENT_REGEXP_MANDATORY => (
    qname => qr/\A [!-?A-~]{1,255} \z/xms,
    rname => qr/\A [*] | [!-()+-<>-~][!-~]* \z/xms,
    cigar => qr/\A [*] | (\d+[MIDNSHPX=])+ \z/xms,
    rnext => qr/\A [*] | = | [!-()+-<>-~][!-~]* \z/xms,
    seq   => qr/\A [*] | [[:alpha:]=.]+ \z/xms,
    qual  => qr/\A [!-~]+ \z/xms,
);

# Ranges for checking alignment line mandatory fields
Readonly our %ALIGNMENT_RANGE_MANDATORY => (
    flag  => [ 0,          2**16 - 1 ],
    pos   => [ 0,          2**29 - 1 ],
    mapq  => [ 0,          2**8 - 1 ],
    pnext => [ 0,          2**29 - 1 ],
    tlen  => [ -2**29 + 1, 2**29 - 1 ],
);

# Regexps for checking alignment line optional fields
Readonly our %ALIGNMENT_REGEXP_OPTIONAL => (
    A => qr/\A [!-~] \z/xms,
    i => qr/\A [-+]?\d+ \z/xms,
    f => qr/\A [-+]?\d*[.]?\d+([eE][-+]?\d+)? \z/xms,
    Z => qr/\A [ !-~]+ \z/xms,
    H => qr/\A [\dA-F]+ \z/xms,
    B => qr/\A [cCsSiIf](,[-+]?\d*[.]?\d+([eE][-+]?\d+)?)+ \z/xms,
);

# Chance one read of a pair is unmapped
Readonly our $CHANCE_UNMAPPED => 0.1;

# Default options
## no critic (ProhibitMagicNumbers)
my $seed;
my $seq_region_count      = 1;
my $seq_region_max_length = 1_000_000;
my $read_pair_count       = 100;
my @read_tags;
my $read1_length = 30;
my $read2_length = 54;
my ( $help, $man );
## use critic

# Get and check command line options
get_and_check_options();

# Ensure reproducible chromosome lengths if seed set
if ( defined $seed ) {
    srand $seed;
}

# Construct command line
my @cl = ('make_test_sam.pl');
if ($seed) {
    push @cl, '--seed', $seed;
}
push @cl, '--seq_region_count',      $seq_region_count;
push @cl, '--seq_region_max_length', $seq_region_max_length;
push @cl, '--read_pair_count',       $read_pair_count;
push @cl, '--read_tags',             @read_tags;
push @cl, '--read1_length',          $read1_length;
push @cl, '--read2_length',          $read2_length;
my $cl = join q{ }, @cl;

# Print HD and RG SAM header
print_or_die( header_line( 'HD', [ 'VN', '1.4' ], [ 'SO', 'unsorted' ] ) );
print_or_die( header_line( 'RG', [ 'ID', q{1} ],  [ 'SM', 'TC' ] ) );
print_or_die(
    header_line(
        'PG',
        [ 'ID', q{1} ],
        [ 'PN', 'make_test_sam.pl' ],
        [ 'CL', $cl ]
    )
);

# Make each chromosome of random length and print SQ SAM headers
my %length_of;
foreach my $seq_region ( 1 .. $seq_region_count ) {
    my $length = int rand( $seq_region_max_length + 1 );
    $length_of{$seq_region} = $length;
    print_or_die(
        header_line( 'SQ', [ 'SN', $seq_region ], [ 'LN', $length ] ) );
}

# Ensure alignments are always random
srand;

# Generate start of each read name
## no critic (ProhibitMagicNumbers)
my $qname_base = 'HS';
$qname_base .= ( int rand 50 ) + 1;        # Instrument name
$qname_base .= q{_};
$qname_base .= ( int rand 20_000 ) + 1;    # Run
$qname_base .= q{:};
$qname_base .= ( int rand 8 ) + 1;         # Flowcell lane
$qname_base .= q{:};
## use critic

# Generate alignments for each chromosome one by one
foreach my $seq_region ( 1 .. $seq_region_count ) {
    foreach ( 1 .. $read_pair_count ) {
        my $read1_qname = get_qname( $qname_base, get_read_tag() );
        my $read2_qname = $read1_qname;    # Always the same
        my ( $read1_pos, $read2_pos ) =
          get_pos( $length_of{$seq_region}, $read1_length, $read2_length );
        my ( $read1_flag, $read2_flag ) = get_flag( $read1_pos, $read2_pos );
        my ( $read1_tlen, $read2_tlen ) =
          get_tlen( $read1_pos, $read2_pos, $read1_length, $read2_length );
        ( $read1_flag, $read2_flag, $read1_pos, $read2_pos ) =
          get_unmapped( $read1_flag, $read2_flag, $read1_pos, $read2_pos );
        my ($read1_cigar) = get_cigar($read1_length);
        my ($read2_cigar) = get_cigar($read2_length);
        my ($read1_nm)    = get_nm();
        my ($read2_nm)    = get_nm();

        # Rarely generate 50 to 99 real duplicates to simulate peaks
        ## no critic (ProhibitMagicNumbers)
        my $num_real_duplicates = 0;
        if ( rand $read_pair_count < 2 ) {
            $num_real_duplicates = int( rand 50 ) + 50;
        }
        ## use critic

        # Generate PCR duplicates (i.e. marked as duplicates)
        ## no critic (ProhibitMagicNumbers)
        my $num_pcr_duplicates = poisson_number(0.6);
        ## use critic

        my $num_duplicates = $num_real_duplicates + $num_pcr_duplicates;
        foreach my $read_pair_count ( 1 .. $num_duplicates + 1 ) {

            # First read
            print_or_die(
                alignment_line(
                    qname => $read1_qname,
                    flag  => $read1_flag,
                    rname => $seq_region,
                    pos   => $read1_pos,
                    mapq  => 255,
                    cigar => $read1_cigar,
                    rnext => q{=},
                    pnext => $read2_pos,
                    tlen  => $read1_tlen,
                    seq   => get_seq($read1_length),
                    qual  => get_qual($read1_length),
                    opt   => {
                        'NM:i' => $read1_nm,
                        'RG:Z' => q{1},
                    },
                )
            );

            # Second read
            print_or_die(
                alignment_line(
                    qname => $read2_qname,
                    flag  => $read2_flag,
                    rname => $seq_region,
                    pos   => $read2_pos,
                    mapq  => 255,
                    cigar => $read2_cigar,
                    rnext => q{=},
                    pnext => $read1_pos,
                    tlen  => $read2_tlen,
                    seq   => get_seq($read2_length),
                    qual  => get_qual($read2_length),
                    opt   => {
                        'NM:i' => $read2_nm,
                        'RG:Z' => q{1},
                    },
                )
            );

            # New read name for duplicates
            $read1_qname = get_qname( $qname_base, get_read_tag() );
            $read2_qname = $read1_qname;    # Always the same

            if ( $read_pair_count == $num_real_duplicates + 1 ) {

                # Mark rest of reads as duplicates
                $read1_flag = $read1_flag | $DETCT::Misc::BAM::Flag::DUPLICATE;
                $read2_flag = $read2_flag | $DETCT::Misc::BAM::Flag::DUPLICATE;
            }
        }
    }
}

# Generate SAM header line
sub header_line {
    my ( $record_type, @data ) = @_;

    my $header_line = q{};

    if ( $record_type !~ m/\A [[:alpha:]][[:alpha:]] \z/xms ) {
        confess 'Invalid record type (', $record_type, q{)};
    }

    $header_line .= q{@} . $record_type;

    foreach my $datum (@data) {
        if ( ref $datum ne 'ARRAY' ) {
            confess 'Arrayref of tag / value pairs is required (not ',
              ref $datum, q{)};
        }
        my ( $tag, $value ) = @{$datum};
        if ( $tag !~ m/\A [[:alpha:]][[:alpha:]\d] \z/xms ) {
            confess 'Invalid tag (', $tag, q{)};
        }
        if ( $value !~ m/\A [ -~]+ \z/xms ) {
            confess 'Invalid value (', $value, q{)};
        }

        $header_line .= "\t" . $tag . q{:} . $value;
    }

    $header_line .= "\n";

    return $header_line;
}

# Generate SAM alignment line
sub alignment_line {
    my (%data) = @_;

    # Check string fields
    foreach my $field ( sort keys %ALIGNMENT_REGEXP_MANDATORY ) {
        if ( $data{$field} !~ $ALIGNMENT_REGEXP_MANDATORY{$field} ) {
            confess 'Invalid ', uc $field, ' (', $data{$field}, q{)};
        }
    }

    # Check int fields
    foreach my $field ( sort keys %ALIGNMENT_RANGE_MANDATORY ) {
        if (   $data{$field} < $ALIGNMENT_RANGE_MANDATORY{$field}->[0]
            || $data{$field} > $ALIGNMENT_RANGE_MANDATORY{$field}->[1] )
        {
            confess 'Invalid ', uc $field, ' (', $data{$field}, q{)};
        }
    }

    # Mandatory fields
    my $alignment_line = join "\t", $data{qname}, $data{flag}, $data{rname},
      $data{pos}, $data{mapq}, $data{cigar}, $data{rnext}, $data{pnext},
      $data{tlen}, $data{seq}, $data{qual};

    # Optional fields
    if ( exists $data{opt} ) {
        foreach my $tag_type ( keys %{ $data{opt} } ) {
            my $value = $data{opt}->{$tag_type};
            my ( $tag, $type ) = split /:/xms, $tag_type;

            # Validate tag
            if ( $tag !~ /\A [[:alpha:]][[:alpha:]\d] \z/xms ) {
                confess 'Invalid tag (', $tag, q{)};
            }

            # Validate type
            if ( !exists $ALIGNMENT_REGEXP_OPTIONAL{$type} ) {
                confess 'Invalid type (', $type, q{)};
            }

            # Validate value
            if ( $value !~ $ALIGNMENT_REGEXP_OPTIONAL{$type} ) {
                confess 'Invalid ', $tag, ' (', $value, q{)};
            }

            $alignment_line .= "\t";
            $alignment_line .= join q{:}, $tag, $type, $value;
        }
    }

    $alignment_line .= "\n";

    return $alignment_line;
}

# Get a random read tag and substitute random bases
sub get_read_tag {
    my $tag = $read_tags[ int rand $#read_tags + 1 ];

    # Replace IUPAC code with random bases
    $tag =~ s/ N / qw( A G C T )[ int rand 4 ] /xmsge;
    $tag =~ s/ B / qw( G C T   )[ int rand 3 ] /xmsge;
    $tag =~ s/ D / qw( A G T   )[ int rand 3 ] /xmsge;
    $tag =~ s/ H / qw( A C T   )[ int rand 3 ] /xmsge;
    $tag =~ s/ V / qw( A G C   )[ int rand 3 ] /xmsge;
    $tag =~ s/ R / qw( A G     )[ int rand 2 ] /xmsge;
    $tag =~ s/ Y / qw( C T     )[ int rand 2 ] /xmsge;
    $tag =~ s/ K / qw( G T     )[ int rand 2 ] /xmsge;
    $tag =~ s/ M / qw( A C     )[ int rand 2 ] /xmsge;
    $tag =~ s/ S / qw( G C     )[ int rand 2 ] /xmsge;
    $tag =~ s/ W / qw( A T     )[ int rand 2 ] /xmsge;

    return $tag;
}

# Construct read name
sub get_qname {
    my ( $qname, $read_tag ) = @_;

    ## no critic (ProhibitMagicNumbers)
    $qname .= ( int rand 3_000 ) + 1;      # Tile number
    $qname .= q{:};
    $qname .= ( int rand 20_000 ) + 1;     # Cluster x coordinate
    $qname .= q{:};
    $qname .= ( int rand 200_000 ) + 1;    # Cluster y coordinate
    $qname .= q{#};
    $qname .= $read_tag;
    ## use critic

    return $qname;
}

# Get position for both reads
sub get_pos {
    my ( $seq_region_len, $read1_len, $read2_len ) = @_;

    my ( $read1_pos, $read2_pos );

    my $pair_ok = 0;

    while ( !$pair_ok ) {
        $read1_pos = ( int rand $seq_region_len ) + 1;
        $read2_pos = ( int rand $seq_region_len ) + 1;
        $pair_ok   = 1;

        my $read1_end = $read1_pos + $read1_len - 1;
        my $read2_end = $read2_pos + $read2_len - 1;

        # Check reads are within seq region
        if ( $read1_end > $seq_region_len ) {
            $pair_ok = 0;
        }
        if ( $read2_end > $seq_region_len ) {
            $pair_ok = 0;
        }

        # Check reads don't overlap
        if ( $read1_pos <= $read2_end && $read1_end >= $read2_pos ) {
            $pair_ok = 0;
        }
    }

    return $read1_pos, $read2_pos;
}

# Get flags for both reads (http://picard.sourceforge.net/explain-flags.html)
sub get_flag {
    my ( $read1_pos, $read2_pos ) = @_;

    my $read1_flag =
      $DETCT::Misc::BAM::Flag::READ_PAIRED |
      $DETCT::Misc::BAM::Flag::PROPER_PAIR;
    my $read2_flag =
      $DETCT::Misc::BAM::Flag::READ_PAIRED |
      $DETCT::Misc::BAM::Flag::PROPER_PAIR;

    if ( $read1_pos < $read2_pos ) {
        $read1_flag =
          $read1_flag | $DETCT::Misc::BAM::Flag::MATE_REVERSE_STRAND;
        $read2_flag =
          $read2_flag | $DETCT::Misc::BAM::Flag::READ_REVERSE_STRAND;
    }
    else {
        $read1_flag =
          $read1_flag | $DETCT::Misc::BAM::Flag::READ_REVERSE_STRAND;
        $read2_flag =
          $read2_flag | $DETCT::Misc::BAM::Flag::MATE_REVERSE_STRAND;
    }

    $read1_flag = $read1_flag | $DETCT::Misc::BAM::Flag::FIRST_IN_PAIR;
    $read2_flag = $read2_flag | $DETCT::Misc::BAM::Flag::SECOND_IN_PAIR;

    return $read1_flag, $read2_flag;
}

# Get template length for both reads
sub get_tlen {
    my ( $read1_pos, $read2_pos, $read1_len, $read2_len ) = @_;

    my ( $read1_tlen, $read2_tlen );

    if ( $read1_pos < $read2_pos ) {
        $read1_tlen = $read2_pos - $read1_pos + $read2_len;
        $read2_tlen = -$read1_tlen;
    }
    else {
        $read2_tlen = $read1_pos - $read2_pos + $read1_len;
        $read1_tlen = -$read2_tlen;
    }

    return $read1_tlen, $read2_tlen;
}

# Adjust flags and positions if a read is unmapped
sub get_unmapped {
    my ( $read1_flag, $read2_flag, $read1_pos, $read2_pos ) = @_;

    if ( rand() < $CHANCE_UNMAPPED ) {
        $read1_flag = $read1_flag ^ $DETCT::Misc::BAM::Flag::PROPER_PAIR;
        $read2_flag = $read2_flag ^ $DETCT::Misc::BAM::Flag::PROPER_PAIR;
        if ( rand() < 0.5 ) {    ## no critic (ProhibitMagicNumbers)
                                 # Read 1 unmapped
            $read1_flag = $read1_flag | $DETCT::Misc::BAM::Flag::READ_UNMAPPED;
            $read2_flag = $read2_flag | $DETCT::Misc::BAM::Flag::MATE_UNMAPPED;
            $read1_pos  = $read2_pos;
        }
        else {
            # Read 2 unmapped
            $read2_flag = $read2_flag | $DETCT::Misc::BAM::Flag::READ_UNMAPPED;
            $read1_flag = $read1_flag | $DETCT::Misc::BAM::Flag::MATE_UNMAPPED;
            $read2_pos  = $read1_pos;
        }
    }

    return $read1_flag, $read2_flag, $read1_pos, $read2_pos;
}

# Get sequence (just random)
sub get_seq {
    my ($read_len) = @_;

    ## no critic (ProhibitMagicNumbers)
    return join q{}, map { qw( A G C T ) [ int rand 4 ] } 1 .. $read_len;
    ## use critic
}

# Get CIGAR string containing random soft clipping
sub get_cigar {
    my ($read_len) = @_;

    my $m = $read_len;    # Length of alignment match

    ## no critic (ProhibitMagicNumbers)
    my $s1 = poisson_number(0.7);    # Soft clipping at start of alignment
    my $s2 = poisson_number(0.7);    # Soft clipping at end of alignment
    ## use critic

    $m = $m - $s1 - $s2;

    # Construct CIGAR

    my $cigar = $m . q{M};

    if ($s1) {
        $cigar = $s1 . q{S} . $cigar;
    }
    if ($s2) {
        $cigar = $cigar . $s2 . q{S};
    }

    return $cigar;
}

# Get quality
sub get_qual {
    my ($read_len) = @_;

    return q{~} x $read_len;
}

# Get random number of mismatches for a read
sub get_nm {
    ## no critic (ProhibitMagicNumbers)
    return poisson_number(0.6);    # ~ e^-0.5, so skewed towards 0 and 1
    ## use critic
}

# Generate random Poisson-distributed number using Knuth's algorithm
sub poisson_number {
    my ($l) = @_;                  # e^-lambda

    my $k = 0;
    my $p = 1;

    while ( $p > $l ) {
        $k++;
        $p = $p * rand;
    }

    return $k - 1;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'seed=i'                  => \$seed,
        'seq_region_count=i'      => \$seq_region_count,
        'seq_region_max_length=i' => \$seq_region_max_length,
        'read_pair_count=i'       => \$read_pair_count,
        'read_tags=s@{1,}'        => \@read_tags,
        'read1_length=i'          => \$read1_length,
        'read2_length=i'          => \$read2_length,
        'help'                    => \$help,
        'man'                     => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$seq_region_count ) {
        pod2usage("--seq_region_count must be a positive integer\n");
    }
    if ( !$seq_region_max_length ) {
        pod2usage("--seq_region_max_length must be a positive integer\n");
    }
    if ( !$read_pair_count ) {
        pod2usage("--read_pair_count must be a positive integer\n");
    }
    if ( !$read1_length ) {
        pod2usage("--read1_length must be a positive integer\n");
    }
    if ( !$read2_length ) {
        pod2usage("--read2_length must be a positive integer\n");
    }
    if ( !@read_tags ) {
        pod2usage("--read_tags must be specified\n");
    }

    return;
}

=head1 USAGE

    make_test_sam.pl
        [--seed seed]
        [--seq_region_count int]
        [--seq_region_max_length int]
        [--read_pair_count int]
        [--read_tags tags...]
        [--read1_length int]
        [--read2_length int]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--seed INT>

Random seed (to get reproducible chromosome lengths).

=item B<--seq_region_count INT>

Number of seq regions (default to 1).

=item B<--seq_region_max_length INT>

Maximum length of each seq region (defaults to 1,000,000 bp).

=item B<--read_pair_count INT>

Number of read pairs aligned to each seq region (defaults to 100).

=item B<--read_tags TAGS>

Read tags.

=item B<--read1_length INT>

Length of read 1 after trimming (defaults to 30 bp).

=item B<--read2_length INT>

Length of read 2 (defaults to 54 bp).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
