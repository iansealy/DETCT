use warnings;
use strict;
use autodie;
use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 117;

use DETCT::Misc::Tag qw(
  detag_trim_fastq
  convert_tag_to_regexp
);

use File::Temp qw( tempdir );
use Path::Tiny;

=for comment

Test random FASTQ files can be regenerated using:

perl -Ilib script/make_test_fastq.pl --seed 1 --output_prefix test1 \
--read_tags NNNNBGAGGC NNNNBAGAAG
perl -Ilib script/make_test_fastq.pl --seed 2 --output_prefix test2 \
--read_tags NNNNBGAGGC NNNNBAGAAG --read_length 54
perl -Ilib script/make_test_fastq.pl --seed 3 --output_prefix test3 \
--read_tags NNNNBBBBNNNNATCACGTT NNNNBBBBNNNNCGATGTTT
mv test* t/data/

Some numbers in tests below will then need updating.

test1	NNNNBGAGGC:	25
test1	NNNNBAGAAG:	24
test1	XXXXXXXXXX:	51

test2	NNNNBGAGGC:	24
test2	NNNNBAGAAG:	35
test2	XXXXXXXXXX:	41

test3	NNNNBBBBNNNNATCACGTT:	27
test3	NNNNBBBBNNNNCGATGTTT:	21
test3	XXXXXXXXXXXXXXXXXXXX:	52

=cut

my $tmp_dir = tempdir( CLEANUP => 1 );

# Check detagging and trimming FASTQ files
is(
    detag_trim_fastq(
        {
            fastq_read1_input     => 't/data/test1_1.fastq',
            fastq_read2_input     => 't/data/test1_2.fastq',
            fastq_output_prefix   => $tmp_dir . '/test1',
            read1_required_length => 30,
            read2_required_length => 54,
            polyt_trim_length     => 14,
            polyt_min_length      => 10,
            read_tags             => [ 'NNNNBGAGGC', 'NNNNBAGAAG' ],
        }
    ),
    undef,
    'Detag and trim FASTQ'
);

my @fastq;

@fastq = path( $tmp_dir . '/test1_NNNNBGAGGC_1.fastq' )->lines;
is( scalar @fastq / 4, 25, '25 read 1s' );
@fastq = path( $tmp_dir . '/test1_NNNNBGAGGC_2.fastq' )->lines;
is( scalar @fastq / 4, 25, '25 read 2s' );
@fastq = path( $tmp_dir . '/test1_NNNNBAGAAG_1.fastq' )->lines;
is( scalar @fastq / 4, 24, '24 read 1s' );
@fastq = path( $tmp_dir . '/test1_NNNNBAGAAG_2.fastq' )->lines;
is( scalar @fastq / 4, 24, '24 read 2s' );
@fastq = path( $tmp_dir . '/test1_XXXXXXXXXX_1.fastq' )->lines;
is( scalar @fastq / 4, 51, '51 read 1s' );
@fastq = path( $tmp_dir . '/test1_XXXXXXXXXX_2.fastq' )->lines;
is( scalar @fastq / 4, 51, '51 read 2s' );

my $read_name;
my $read_seq;
my $read_qual;

@fastq     = path( $tmp_dir . '/test1_NNNNBGAGGC_1.fastq' )->lines;
$read_name = $fastq[0];
chomp $read_name;
is( substr( $read_name, -7 ), 'GAGGC/1', 'Tag added to read name' );
$read_seq = $fastq[1];
chomp $read_seq;
is( length $read_seq, 30, 'Sequence trimmed to 30 bp' );
$read_qual = $fastq[3];
chomp $read_qual;
is( length $read_qual, 30, 'Quality trimmed to 30 bp' );

@fastq     = path( $tmp_dir . '/test1_NNNNBGAGGC_2.fastq' )->lines;
$read_name = $fastq[0];
chomp $read_name;
is( substr( $read_name, -7 ), 'GAGGC/2', 'Tag added to read name' );
$read_seq = $fastq[1];
chomp $read_seq;
is( length $read_seq, 54, 'Sequence trimmed to 54 bp' );
$read_qual = $fastq[3];
chomp $read_qual;
is( length $read_qual, 54, 'Quality trimmed to 54 bp' );

@fastq     = path( $tmp_dir . '/test1_XXXXXXXXXX_1.fastq' )->lines;
$read_name = $fastq[0];
chomp $read_name;
is( substr( $read_name, -13 ), '#XXXXXXXXXX/1', 'Tag added to read name' );
$read_seq = $fastq[1];
chomp $read_seq;
is( length $read_seq, 54, 'Sequence trimmed to 54 bp' );
$read_qual = $fastq[3];
chomp $read_qual;
is( length $read_qual, 54, 'Quality trimmed to 54 bp' );

@fastq     = path( $tmp_dir . '/test1_XXXXXXXXXX_2.fastq' )->lines;
$read_name = $fastq[0];
chomp $read_name;
is( substr( $read_name, -13 ), '#XXXXXXXXXX/2', 'Tag added to read name' );
$read_seq = $fastq[1];
chomp $read_seq;
is( length $read_seq, 54, 'Sequence trimmed to 54 bp' );
$read_qual = $fastq[3];
chomp $read_qual;
is( length $read_qual, 54, 'Quality trimmed to 54 bp' );

# Check detagging and trimming FASTQ files
is(
    detag_trim_fastq(
        {
            fastq_read1_input     => 't/data/test2_1.fastq',
            fastq_read2_input     => 't/data/test2_2.fastq',
            fastq_output_prefix   => $tmp_dir . '/test2',
            read1_required_length => 30,
            read2_required_length => 54,
            polyt_trim_length     => 14,
            polyt_min_length      => 10,
            read_tags             => [ 'NNNNBGAGGC', 'NNNNBAGAAG' ],
        }
    ),
    undef,
    'Detag and trim FASTQ'
);

@fastq = path( $tmp_dir . '/test2_NNNNBGAGGC_1.fastq' )->lines;
is( scalar @fastq / 4, 24, '30 read 1s' );
@fastq = path( $tmp_dir . '/test2_NNNNBGAGGC_2.fastq' )->lines;
is( scalar @fastq / 4, 24, '30 read 2s' );
@fastq = path( $tmp_dir . '/test2_NNNNBAGAAG_1.fastq' )->lines;
is( scalar @fastq / 4, 35, '24 read 1s' );
@fastq = path( $tmp_dir . '/test2_NNNNBAGAAG_2.fastq' )->lines;
is( scalar @fastq / 4, 35, '24 read 2s' );
@fastq = path( $tmp_dir . '/test2_XXXXXXXXXX_1.fastq' )->lines;
is( scalar @fastq / 4, 41, '46 read 1s' );
@fastq = path( $tmp_dir . '/test2_XXXXXXXXXX_2.fastq' )->lines;
is( scalar @fastq / 4, 41, '46 read 2s' );

throws_ok {
    detag_trim_fastq(
        {
            fastq_read1_input     => 't/data/test1_1.fastq',
            fastq_read2_input     => 't/data/test2_2.fastq',
            fastq_output_prefix   => $tmp_dir . '/test',
            read1_required_length => 30,
            read2_required_length => 54,
            polyt_trim_length     => 14,
            polyt_min_length      => 10,
            read_tags             => [ 'NNNNBGAGGC', 'NNNNBAGAAG' ],
        }
    );
}
qr/Read order does not match in input/ms, 'FASTQ files not matched';

# Check detagging and trimming FASTQ files
is(
    detag_trim_fastq(
        {
            fastq_read1_input     => 't/data/test3_1.fastq',
            fastq_read2_input     => 't/data/test3_2.fastq',
            fastq_output_prefix   => $tmp_dir . '/test3',
            read1_required_length => 30,
            read2_required_length => 54,
            polyt_trim_length     => 14,
            polyt_min_length      => 10,
            read_tags => [ 'NNNNBBBBNNNNATCACGTT', 'NNNNBBBBNNNNCGATGTTT' ],
        }
    ),
    undef,
    'Detag and trim FASTQ'
);

@fastq     = path( $tmp_dir . '/test3_NNNNBBBBNNNNATCACGTT_1.fastq' )->lines;
$read_name = $fastq[0];
chomp $read_name;
is( substr( $read_name, -10 ), 'ATCACGTT/1', 'Tag added to read name' );
$read_seq = $fastq[1];
chomp $read_seq;
is( length $read_seq, 30, 'Sequence trimmed to 30 bp' );
$read_qual = $fastq[3];
chomp $read_qual;
is( length $read_qual, 30, 'Quality trimmed to 30 bp' );

@fastq     = path( $tmp_dir . '/test3_NNNNBBBBNNNNATCACGTT_2.fastq' )->lines;
$read_name = $fastq[0];
chomp $read_name;
is( substr( $read_name, -10 ), 'ATCACGTT/2', 'Tag added to read name' );
$read_seq = $fastq[1];
chomp $read_seq;
is( length $read_seq, 54, 'Sequence trimmed to 54 bp' );
$read_qual = $fastq[3];
chomp $read_qual;
is( length $read_qual, 54, 'Quality trimmed to 54 bp' );

# Check detagging and trimming FASTQ files
is(
    detag_trim_fastq(
        {
            fastq_read1_input        => 't/data/test1_1.fastq',
            fastq_read2_input        => 't/data/test1_2.fastq',
            fastq_output_prefix      => $tmp_dir . '/test4',
            read1_required_length    => 30,
            read2_required_length    => 54,
            read1_5prime_trim_length => 2,
            polyt_trim_length        => 14,
            polyt_min_length         => 10,
            read_tags                => [ 'NNNNBGAGGC', 'NNNNBAGAAG' ],
        }
    ),
    undef,
    'Detag and trim FASTQ'
);
is(
    detag_trim_fastq(
        {
            fastq_read1_input        => 't/data/test1_1.fastq',
            fastq_read2_input        => 't/data/test1_2.fastq',
            fastq_output_prefix      => $tmp_dir . '/test5',
            read1_required_length    => 30,
            read2_required_length    => 54,
            read2_5prime_trim_length => 2,
            polyt_trim_length        => 14,
            polyt_min_length         => 10,
            read_tags                => [ 'NNNNBGAGGC', 'NNNNBAGAAG' ],
        }
    ),
    undef,
    'Detag and trim FASTQ'
);

my @fastq_read1        = path( $tmp_dir . '/test5_NNNNBGAGGC_1.fastq' )->lines;
my @fastq_5prime_read1 = path( $tmp_dir . '/test4_NNNNBGAGGC_1.fastq' )->lines;
my @fastq_read2        = path( $tmp_dir . '/test4_NNNNBGAGGC_2.fastq' )->lines;
my @fastq_5prime_read2 = path( $tmp_dir . '/test5_NNNNBGAGGC_2.fastq' )->lines;

chomp @fastq_read1;
chomp @fastq_5prime_read1;
chomp @fastq_read2;
chomp @fastq_5prime_read2;

my $seq_qual_match_read1 = 1;
my $names_match_read1    = 1;
my $seq_qual_match_read2 = 1;
my $names_match_read2    = 1;
foreach my $i ( 0 .. ( ( scalar @fastq_read1 ) / 2 ) - 1 ) {
    if ( $fastq_read1[ $i * 2 ] ne $fastq_5prime_read1[ $i * 2 ] ) {
        $names_match_read1 = 0;
    }
    if ( $fastq_read2[ $i * 2 ] ne $fastq_5prime_read2[ $i * 2 ] ) {
        $names_match_read2 = 0;
    }
    my $untrimmed1 = $fastq_read1[ $i * 2 + 1 ];
    my $trimmed1   = $fastq_5prime_read1[ $i * 2 + 1 ];
    if ( ( substr $untrimmed1, 2 ) ne substr $trimmed1,
        0, ( length $trimmed1 ) - 2 )
    {
        $seq_qual_match_read1 = 0;
    }
    my $untrimmed2 = $fastq_read2[ $i * 2 + 1 ];
    my $trimmed2   = $fastq_5prime_read2[ $i * 2 + 1 ];
    if ( ( substr $untrimmed2, 2 ) ne substr $trimmed2,
        0, ( length $trimmed2 ) - 2 )
    {
        $seq_qual_match_read2 = 0;
    }
}
is( $seq_qual_match_read1, 1, q{Read 1 5' end trimming sequence or quality} );
is( $seq_qual_match_read2, 1, q{Read 2 5' end trimming sequence or quality} );
is( $names_match_read1,    1, q{Read 1 5' end trimming names} );
is( $names_match_read2,    1, q{Read 2 5' end trimming names} );

# Check converting tags to regular expressions

my @tags = qw( N B D H V R Y K M S W A G C T AA );

my %re_for = convert_tag_to_regexp(@tags);

ok( q{A} =~ $re_for{A}->[0], 'A matches A' );
ok( q{G} !~ $re_for{A}->[0], 'A does not match G' );
ok( q{C} !~ $re_for{A}->[0], 'A does not match C' );
ok( q{T} !~ $re_for{A}->[0], 'A does not match T' );
ok( q{N} !~ $re_for{A}->[0], 'A does not match N' );

ok( q{A} !~ $re_for{G}->[0], 'G does not match A' );
ok( q{G} =~ $re_for{G}->[0], 'G matches G' );
ok( q{C} !~ $re_for{G}->[0], 'G does not match C' );
ok( q{T} !~ $re_for{G}->[0], 'G does not match T' );
ok( q{N} !~ $re_for{G}->[0], 'G does not match N' );

ok( q{A} !~ $re_for{C}->[0], 'C does not match A' );
ok( q{G} !~ $re_for{C}->[0], 'C does not match G' );
ok( q{C} =~ $re_for{C}->[0], 'C matches C' );
ok( q{T} !~ $re_for{C}->[0], 'C does not match T' );
ok( q{N} !~ $re_for{C}->[0], 'C does not match N' );

ok( q{A} !~ $re_for{T}->[0], 'T does not match A' );
ok( q{G} !~ $re_for{T}->[0], 'T does not match G' );
ok( q{C} !~ $re_for{T}->[0], 'T does not match C' );
ok( q{T} =~ $re_for{T}->[0], 'T matches T' );
ok( q{N} !~ $re_for{T}->[0], 'T does not match N' );

ok( q{A} =~ $re_for{R}->[0], 'R matches A' );
ok( q{G} =~ $re_for{R}->[0], 'R matches G' );
ok( q{C} !~ $re_for{R}->[0], 'R does not match C' );
ok( q{T} !~ $re_for{R}->[0], 'R does not match T' );
ok( q{N} !~ $re_for{R}->[0], 'R does not match N' );

ok( q{A} !~ $re_for{Y}->[0], 'Y does not match A' );
ok( q{G} !~ $re_for{Y}->[0], 'Y does not match G' );
ok( q{C} =~ $re_for{Y}->[0], 'Y matches C' );
ok( q{T} =~ $re_for{Y}->[0], 'Y matches T' );
ok( q{N} !~ $re_for{Y}->[0], 'Y does not match N' );

ok( q{A} !~ $re_for{S}->[0], 'S does not match A' );
ok( q{G} =~ $re_for{S}->[0], 'S matches G' );
ok( q{C} =~ $re_for{S}->[0], 'S matches C' );
ok( q{T} !~ $re_for{S}->[0], 'S does not match T' );
ok( q{N} !~ $re_for{S}->[0], 'S does not match N' );

ok( q{A} =~ $re_for{W}->[0], 'W matches A' );
ok( q{G} !~ $re_for{W}->[0], 'W does not match G' );
ok( q{C} !~ $re_for{W}->[0], 'W does not match C' );
ok( q{T} =~ $re_for{W}->[0], 'W matches T' );
ok( q{N} !~ $re_for{W}->[0], 'W does not match N' );

ok( q{A} !~ $re_for{K}->[0], 'K does not match A' );
ok( q{G} =~ $re_for{K}->[0], 'K matches G' );
ok( q{C} !~ $re_for{K}->[0], 'K does not match C' );
ok( q{T} =~ $re_for{K}->[0], 'K matches T' );
ok( q{N} !~ $re_for{K}->[0], 'K does not match N' );

ok( q{A} =~ $re_for{M}->[0], 'M matches A' );
ok( q{G} !~ $re_for{M}->[0], 'M does not match G' );
ok( q{C} =~ $re_for{M}->[0], 'M matches C' );
ok( q{T} !~ $re_for{M}->[0], 'M does not match T' );
ok( q{N} !~ $re_for{M}->[0], 'M does not match N' );

ok( q{A} !~ $re_for{B}->[0], 'B does not match A' );
ok( q{G} =~ $re_for{B}->[0], 'B matches G' );
ok( q{C} =~ $re_for{B}->[0], 'B matches C' );
ok( q{T} =~ $re_for{B}->[0], 'B matches T' );
ok( q{N} !~ $re_for{B}->[0], 'B does not match N' );

ok( q{A} =~ $re_for{D}->[0], 'D matches A' );
ok( q{G} =~ $re_for{D}->[0], 'D matches G' );
ok( q{C} !~ $re_for{D}->[0], 'D does not match C' );
ok( q{T} =~ $re_for{D}->[0], 'D matches T' );
ok( q{N} !~ $re_for{D}->[0], 'D does not match N' );

ok( q{A} =~ $re_for{H}->[0], 'M matches A' );
ok( q{G} !~ $re_for{H}->[0], 'M does not match G' );
ok( q{C} =~ $re_for{H}->[0], 'M matches C' );
ok( q{T} =~ $re_for{H}->[0], 'M matches T' );
ok( q{N} !~ $re_for{H}->[0], 'M does not match N' );

ok( q{A} =~ $re_for{V}->[0], 'V matches A' );
ok( q{G} =~ $re_for{V}->[0], 'V matches G' );
ok( q{C} =~ $re_for{V}->[0], 'V matches C' );
ok( q{T} !~ $re_for{V}->[0], 'V does not match T' );
ok( q{N} !~ $re_for{V}->[0], 'V does not match N' );

ok( q{A} =~ $re_for{N}->[0], 'N matches A' );
ok( q{G} =~ $re_for{N}->[0], 'N matches G' );
ok( q{C} =~ $re_for{N}->[0], 'N matches C' );
ok( q{T} =~ $re_for{N}->[0], 'N matches T' );
ok( q{N} =~ $re_for{N}->[0], 'N matches N' );

ok( q{A} !~ $re_for{AA}->[0],  'AA does not match A' );
ok( q{AA} =~ $re_for{AA}->[0], 'AA matches AA' );
