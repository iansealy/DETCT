use Test::More;
use Test::Exception;
use Test::Warn;
use Test::DatabaseRow;
use Test::MockObject;
use Carp;

plan tests => 145;

use DETCT::Analysis;

use File::Path qw( make_path );
use IO::Socket::INET;

my $is_ensembl_reachable = is_ensembl_reachable();

my $analysis = DETCT::Analysis->new(
    {
        name               => 'zmp_ph1',
        read1_length       => 30,
        read2_length       => 54,
        mismatch_threshold => 2,
        bin_size           => 100,
        peak_buffer_width  => 100,
        hmm_sig_level      => 0.001,
        hmm_binary         => 'chiphmmnew',
        r_binary           => 'R',
        deseq_script       => 'script/run_deseq.R',
        output_sig_level   => 0.05,
        chunk_total        => 20,
    }
);

isa_ok( $analysis, 'DETCT::Analysis' );

# Test name attribute
is( $analysis->name,                'zmp_ph1', 'Get name' );
is( $analysis->set_name('zmp_ph2'), undef,     'Set name' );
is( $analysis->name,                'zmp_ph2', 'Get new name' );
throws_ok { $analysis->set_name() } qr/No name specified/ms, 'No name';
my $long_name = 'X' x ( $DETCT::Analysis::MAX_NAME_LENGTH + 1 );
throws_ok { $analysis->set_name('') } qr/Empty name specified/ms, 'Empty name';
throws_ok { $analysis->set_name($long_name) } qr/longer than \d+ characters/ms,
  'Long name';

# Test read 1 length attribute
is( $analysis->read1_length,         30,    'Get read 1 length' );
is( $analysis->set_read1_length(40), undef, 'Set read 1 length' );
is( $analysis->read1_length,         40,    'Get new read 1 length' );
throws_ok { $analysis->set_read1_length() } qr/No read 1 length specified/ms,
  'No read 1 length';
throws_ok { $analysis->set_read1_length(-1) } qr/Invalid read 1 length/ms,
  'Invalid read 1 length';

# Test read 2 length attribute
is( $analysis->read2_length,         54,    'Get read 2 length' );
is( $analysis->set_read2_length(64), undef, 'Set read 2 length' );
is( $analysis->read2_length,         64,    'Get new read 2 length' );
throws_ok { $analysis->set_read2_length() } qr/No read 2 length specified/ms,
  'No read 2 length';
throws_ok { $analysis->set_read2_length(-2) } qr/Invalid read 2 length/ms,
  'Invalid read 2 length';

# Test mismatch threshold attribute
is( $analysis->mismatch_threshold,        2,     'Get mismatch threshold' );
is( $analysis->set_mismatch_threshold(3), undef, 'Set mismatch threshold' );
is( $analysis->mismatch_threshold,        3,     'Get new mismatch threshold' );
throws_ok { $analysis->set_mismatch_threshold() }
qr/No mismatch threshold specified/ms, 'No mismatch threshold';
throws_ok { $analysis->set_mismatch_threshold(-1) }
qr/Invalid mismatch threshold/ms, 'Invalid mismatch threshold';

# Test bin size attribute
is( $analysis->bin_size,          100,   'Get bin size' );
is( $analysis->set_bin_size(200), undef, 'Set bin size' );
is( $analysis->bin_size,          200,   'Get new bin size' );
throws_ok { $analysis->set_bin_size() } qr/No bin size specified/ms,
  'No bin size';
throws_ok { $analysis->set_bin_size(-1) } qr/Invalid bin size/ms,
  'Invalid bin size';

# Test peak buffer width attribute
is( $analysis->peak_buffer_width,          100,   'Get peak buffer width' );
is( $analysis->set_peak_buffer_width(200), undef, 'Set peak buffer width' );
is( $analysis->peak_buffer_width,          200,   'Get new peak buffer width' );
throws_ok { $analysis->set_peak_buffer_width() }
qr/No peak buffer width specified/ms, 'No peak buffer width';
throws_ok { $analysis->set_peak_buffer_width(-1) }
qr/Invalid peak buffer width/ms, 'Invalid peak buffer width';

# Test HMM significance level attribute
is( $analysis->hmm_sig_level,          0.001, 'Get HMM significance level' );
is( $analysis->set_hmm_sig_level(0.1), undef, 'Set HMM significance level' );
is( $analysis->hmm_sig_level, 0.1, 'Get new HMM significance level' );
throws_ok { $analysis->set_hmm_sig_level() }
qr/No HMM significance level specified/ms, 'No HMM significance level';
throws_ok { $analysis->set_hmm_sig_level(1) }
qr/Invalid HMM significance level/ms, 'Invalid HMM significance level';

# Test HMM binary attribute
is( $analysis->hmm_binary, 'chiphmmnew', 'Get HMM binary' );
is( $analysis->set_hmm_binary('chiphmmnew'), undef, 'Set HMM binary' );
is( $analysis->hmm_binary, 'chiphmmnew', 'Get new HMM binary' );
throws_ok { $analysis->set_hmm_binary() } qr/No HMM binary specified/ms,
  'No HMM binary';

# Test R binary attribute
is( $analysis->r_binary,          'R',   'Get R binary' );
is( $analysis->set_r_binary('S'), undef, 'Set R binary' );
is( $analysis->r_binary,          'S',   'Get new R binary' );
throws_ok { $analysis->set_r_binary() } qr/No R binary specified/ms,
  'No R binary';

# Test DESeq script attribute
is( $analysis->deseq_script, 'script/run_deseq.R', 'Get DESeq script' );
is( $analysis->set_deseq_script('script'), undef,    'Set DESeq script' );
is( $analysis->deseq_script,               'script', 'Get new DESeq script' );
throws_ok { $analysis->set_deseq_script() } qr/No DESeq script specified/ms,
  'No DESeq script';
throws_ok { $analysis->set_deseq_script('nonexistent') }
qr/does not exist or cannot be read/ms, 'Missing DESeq script';

# Test output significance level attribute
is( $analysis->output_sig_level, 0.05, 'Get output significance level' );
is( $analysis->set_output_sig_level(0.01),
    undef, 'Set output significance level' );
is( $analysis->output_sig_level, 0.01, 'Get new output significance level' );
throws_ok { $analysis->set_output_sig_level() }
qr/No output significance level specified/ms, 'No output significance level';
throws_ok { $analysis->set_output_sig_level(1) }
qr/Invalid output significance level/ms, 'Invalid output significance level';

# Test reference FASTA attribute
is( $analysis->ref_fasta, undef, 'Get reference FASTA' );
is( $analysis->set_ref_fasta('t/data/test12.fa'), undef,
    'Set reference FASTA' );
is( $analysis->ref_fasta, 't/data/test12.fa', 'Get new reference FASTA' );
throws_ok { $analysis->set_ref_fasta('nonexistent') } qr/cannot be read/ms,
  'Missing reference FASTA';

# Test Ensembl host attribute
is( $analysis->ensembl_host, undef, 'Get Ensembl host' );
is( $analysis->set_ensembl_host('ensembldb.ensembl.org'),
    undef, 'Set Ensembl host' );
is( $analysis->ensembl_host, 'ensembldb.ensembl.org', 'Get new Ensembl host' );

# Test Ensembl port attribute
is( $analysis->ensembl_port,           undef, 'Get Ensembl port' );
is( $analysis->set_ensembl_port(3306), undef, 'Set Ensembl port' );
is( $analysis->ensembl_port,           3306,  'Get new Ensembl port' );
throws_ok { $analysis->set_ensembl_port(-1) } qr/Invalid Ensembl port/ms,
  'Invalid Ensembl port';

# Test Ensembl username attribute
is( $analysis->ensembl_user,                  undef, 'Get Ensembl username' );
is( $analysis->set_ensembl_user('anonymous'), undef, 'Set Ensembl username' );
is( $analysis->ensembl_user, 'anonymous', 'Get new Ensembl username' );

# Test Ensembl password attribute
is( $analysis->ensembl_pass,               undef, 'Get Ensembl password' );
is( $analysis->set_ensembl_pass('secret'), undef, 'Set Ensembl password' );
is( $analysis->ensembl_pass, 'secret', 'Get new Ensembl password' );

# Test Ensembl database name attribute
is( $analysis->ensembl_name, undef, 'Get Ensembl database name' );
is( $analysis->set_ensembl_name('zv9_core'),
    undef, 'Set Ensembl database name' );
is( $analysis->ensembl_name, 'zv9_core', 'Get new Ensembl database name' );

# Test Ensembl species attribute
is( $analysis->ensembl_species, undef, 'Get Ensembl species' );
is( $analysis->set_ensembl_species('danio_rerio'),
    undef, 'Set Ensembl species' );
is( $analysis->ensembl_species, 'danio_rerio', 'Get new Ensembl species' );

# Test chunk total attribute
is( $analysis->chunk_total,         20,    'Get chunk total' );
is( $analysis->set_chunk_total(30), undef, 'Set chunk total' );
is( $analysis->chunk_total,         30,    'Get new chunk total' );
throws_ok { $analysis->set_chunk_total() } qr/No chunk total specified/ms,
  'No chunk total';
throws_ok { $analysis->set_chunk_total(-1) } qr/Invalid chunk total/ms,
  'Invalid chunk total';

# Test sequences and chunks before adding samples
my $sequences = $analysis->get_all_sequences();
is( scalar @{$sequences}, 0, 'No sequences' );
my $chunks = $analysis->get_all_chunks();
is( scalar @{$chunks}, 0, 'No chunks' );

# Mock sample object
my $sample = Test::MockObject->new();
$sample->set_isa('DETCT::Sample');
$sample->set_always( 'bam_file', 't/data/test1.bam' );

# Mock sample object with different reference sequence
my $sample_diff = Test::MockObject->new();
$sample_diff->set_isa('DETCT::Sample');
$sample_diff->set_always( 'bam_file', 't/data/test3.bam' );

# Test adding and retrieving samples
my $samples;
$samples = $analysis->get_all_samples();
is( scalar @{$samples},             0,     'No samples' );
is( $analysis->add_sample($sample), undef, 'Add sample' );
$samples = $analysis->get_all_samples();
is( scalar @{$samples}, 1, 'Get one sample' );
$analysis->add_sample($sample);
is( scalar @{$samples}, 2, 'Get two samples' );
throws_ok { $analysis->add_sample($sample_diff) } qr/use different reference/ms,
  'Different reference for sample';
throws_ok { $analysis->add_sample() } qr/No sample specified/ms,
  'No sample specified';
throws_ok { $analysis->add_sample('invalid') } qr/Class of sample/ms,
  'Invalid sample';

# Test sequences and chunks after adding samples
$sequences = $analysis->get_all_sequences();
is( scalar @{$sequences}, 5, '5 sequences' );
$chunks = $analysis->get_all_chunks();
ok( scalar @{$chunks} > 0, 'Chunks' );

# Count sequence in chunks
my $sequence_total = 0;
foreach my $chunk ( @{$chunks} ) {
    $sequence_total += scalar @{$chunk};
}
is( $sequence_total, 5, '5 sequences in chunks' );

# Recalculate chunks so one sequence per chunk
$analysis->set_chunk_total(10000);
$chunks = $analysis->get_all_chunks();
is( scalar @{$chunks}, 5, '5 chunks' );

# Recalculate chunks so 5/3 sequences per chunk on average
$analysis->set_chunk_total(3);
$chunks = $analysis->get_all_chunks();
is( scalar @{$chunks}, 3, '3 chunks' );

# Count sequence in chunks
$sequence_total = 0;
foreach my $chunk ( @{$chunks} ) {
    $sequence_total += scalar @{$chunk};
}
is( $sequence_total, 5, '5 sequences in chunks' );

# Test test chunk attribute
is( $analysis->test_chunk,        undef, 'Get test chunk' );
is( $analysis->set_test_chunk(1), undef, 'Set test chunk' );
is( $analysis->test_chunk,        1,     'Get new test chunk' );
$chunks = $analysis->get_all_chunks();
is( scalar @{$chunks},            1,     '1 chunk' );
is( $analysis->set_test_chunk(4), undef, 'Set test chunk' );
$chunks = $analysis->get_all_chunks();
is( scalar @{$chunks}, 3, '3 chunks' );

# Test constructing from YAML
$analysis = DETCT::Analysis->new_from_yaml('t/data/test_analysis12.yaml');
isa_ok( $analysis, 'DETCT::Analysis' );
$samples = $analysis->get_all_samples();
is( scalar @{$samples}, 2, 'Get two YAML samples' );
throws_ok { $analysis = DETCT::Analysis->new_from_yaml('nonexistent.yaml') }
qr/does not exist or cannot be read/ms, 'Missing YAML file';

# Test validating analysis
throws_ok {
    $analysis = DETCT::Analysis->new_from_yaml('t/data/test_analysis13.yaml');
}
qr/use different reference/ms, 'Different reference';

# Test summary info
$analysis = DETCT::Analysis->new_from_yaml('t/data/test_analysis1122.yaml');
my @bam_files = $analysis->list_all_bam_files();
is( scalar @bam_files, 2, '2 BAM files' );
is( $bam_files[0], 't/data/test1.bam', 'Got BAM file' );
my @tags = $analysis->list_all_tags_by_bam_file('t/data/test1.bam');
is( scalar @tags, 2,            '2 tags' );
is( $tags[0],     'NNNNBAGAAG', 'Got tag' );

my $seq;

# Set FASTA index
$analysis = DETCT::Analysis->new_from_yaml('t/data/test_analysis12.yaml');
throws_ok { $analysis->set_fasta_index(); } qr/No FASTA index specified/ms,
  'No FASTA index';
throws_ok { $analysis->set_fasta_index('invalid'); } qr/Class of FASTA index/ms,
  'Invalid FASTA index';

# Set Ensembl slice adaptor
$analysis = DETCT::Analysis->new_from_yaml('t/data/test_analysis12.yaml');
throws_ok { $analysis->set_slice_adaptor(); }
qr/No Ensembl slice adaptor specified/ms, 'No slice adaptor';
throws_ok { $analysis->set_slice_adaptor('invalid'); }
qr/Class of Ensembl slice adaptor/ms, 'Invalid slice adaptor';

# Get subsequence with missing parameters
$analysis = DETCT::Analysis->new_from_yaml('t/data/test_analysis12.yaml');
throws_ok { $analysis->get_subsequence(); } qr/No sequence name specified/ms,
  'No sequence name';
throws_ok { $analysis->get_subsequence('1'); }
qr/No sequence start specified/ms, 'No sequence start';
throws_ok { $analysis->get_subsequence( '1', 1 ); }
qr/No sequence end specified/ms, 'No sequence end';
throws_ok { $analysis->get_subsequence( '1', 1, 10 ); }
qr/No sequence strand specified/ms, 'No sequence strand';

# Check getting sequence from test FASTA file
# First 10 bp of chromosome 1 should be CCAGGCGCGG according to:

=for comment
head -2 t/data/test12.fa
=cut

$analysis = DETCT::Analysis->new_from_yaml('t/data/test_analysis12.yaml');
$seq = $analysis->get_subsequence( '1', 1, 10, 1 );
is( length $seq, 10,           'FASTA subsequence length' );
is( $seq,        'CCAGGCGCGG', 'FASTA subsequence' );
$seq = $analysis->get_subsequence( '1', 1, 10, -1 );
is( length $seq, 10,           'FASTA reverse complement subsequence length' );
is( $seq,        'CCGCGCCTGG', 'FASTA reverse complement subsequence' );

# Check getting subsequence outside size of sequence
$seq = $analysis->get_subsequence( '1', -1, 10, 1 );
is( length $seq, 10,           'Negative start FASTA subsequence length' );
is( $seq,        'CCAGGCGCGG', 'Negative start FASTA subsequence' );
$seq = $analysis->get_subsequence( '1', -1, -1, 1 );
is( length $seq, 1,   'Negative start and end FASTA subsequence length' );
is( $seq,        'C', 'Negative start and end FASTA subsequence' );
$seq = $analysis->get_subsequence( '1', 1_000_000_001, 1_000_000_010, 1 );
is( length $seq, 0,  'Large start and end FASTA subsequence length' );
is( $seq,        '', 'Large start and end FASTA subsequence' );

# Check getting sequence from Ensembl database
# First 10 bp of chromosome 1 should be TTCTTCTGGG according to:
# http://www.ensembl.org/Danio_rerio/Location/View?r=1%3A1-10
SKIP: {
    skip 'Ensembl not reachable', 4 if !$is_ensembl_reachable;

    $analysis = DETCT::Analysis->new(
        {
            name               => 'zmp_ph1',
            read1_length       => 30,
            read2_length       => 54,
            mismatch_threshold => 2,
            bin_size           => 100,
            peak_buffer_width  => 100,
            hmm_sig_level      => 0.001,
            hmm_binary         => 'chiphmmnew',
            r_binary           => 'R',
            deseq_script       => 'script/run_deseq.R',
            output_sig_level   => 0.05,
            chunk_total        => 20,
            ensembl_species    => 'danio_rerio',
        }
    );
    $seq = $analysis->get_subsequence( '1', 1, 10, 1 );
    is( length $seq, 10,           'Ensembl subsequence length' );
    is( $seq,        'TTCTTCTGGG', 'Ensembl subsequence' );
    $seq = $analysis->get_subsequence( '1', 1, 10, -1 );
    is( length $seq, 10, 'Ensembl reverse complement subsequence length' );
    is( $seq, 'CCCAGAAGAA', 'Ensembl reverse complement subsequence' );
}

# Check getting sequence without FASTA file or Ensembl database
$analysis = DETCT::Analysis->new(
    {
        name               => 'zmp_ph1',
        read1_length       => 30,
        read2_length       => 54,
        mismatch_threshold => 2,
        bin_size           => 100,
        peak_buffer_width  => 100,
        hmm_sig_level      => 0.001,
        hmm_binary         => 'chiphmmnew',
        r_binary           => 'R',
        deseq_script       => 'script/run_deseq.R',
        output_sig_level   => 0.05,
        chunk_total        => 20,
    }
);
throws_ok { $analysis->get_subsequence( '1', 1, 10, 1 ); }
qr/No reference FASTA or Ensembl database/ms, 'No FASTA or Ensembl';

# Check getting sequence from Ensembl database with explicit connection
SKIP: {
    skip 'Ensembl not reachable', 2 if !$is_ensembl_reachable;

    # Avoid warnings about loading the same databases again
    Bio::EnsEMBL::Registry::clear();

    $analysis = DETCT::Analysis->new(
        {
            name               => 'zmp_ph1',
            read1_length       => 30,
            read2_length       => 54,
            mismatch_threshold => 2,
            bin_size           => 100,
            peak_buffer_width  => 100,
            hmm_sig_level      => 0.001,
            hmm_binary         => 'chiphmmnew',
            r_binary           => 'R',
            deseq_script       => 'script/run_deseq.R',
            output_sig_level   => 0.05,
            chunk_total        => 20,
            ensembl_host       => 'ensembldb.ensembl.org',
            ensembl_port       => 5306,
            ensembl_user       => 'anonymous',
            ensembl_pass       => '',
            ensembl_species    => 'danio_rerio',
        }
    );
    $seq = $analysis->get_subsequence( '1', 1, 10, 1 );
    is( length $seq, 10,           'Ensembl subsequence length' );
    is( $seq,        'TTCTTCTGGG', 'Ensembl subsequence' );
}

# Check getting sequence from specific Ensembl database
# Get database name via:

=for comment
mysql -u anonymous -h ensembldb.ensembl.org -P 5306 -Bse \
"SHOW DATABASES LIKE 'danio_rerio_core\_%'" | sort | tail -1
=cut

SKIP: {
    skip 'Ensembl not reachable', 2 if !$is_ensembl_reachable;

    $analysis = DETCT::Analysis->new(
        {
            name               => 'zmp_ph1',
            read1_length       => 30,
            read2_length       => 54,
            mismatch_threshold => 2,
            bin_size           => 100,
            peak_buffer_width  => 100,
            hmm_sig_level      => 0.001,
            hmm_binary         => 'chiphmmnew',
            r_binary           => 'R',
            deseq_script       => 'script/run_deseq.R',
            output_sig_level   => 0.05,
            chunk_total        => 20,
            ensembl_host       => 'ensembldb.ensembl.org',
            ensembl_port       => 5306,
            ensembl_user       => 'anonymous',
            ensembl_pass       => '',
            ensembl_name       => 'danio_rerio_core_73_9',
        }
    );
    $seq = $analysis->get_subsequence( '1', 1, 10, 1 );
    is( length $seq, 10,           'Ensembl subsequence length' );
    is( $seq,        'TTCTTCTGGG', 'Ensembl subsequence' );
}

# Check getting subsequence outside size of sequence
SKIP: {
    skip 'Ensembl not reachable', 6 if !$is_ensembl_reachable;

    $seq = $analysis->get_subsequence( '1', -1, 10, 1 );
    is( length $seq, 10, 'Negative start Ensembl subsequence length' );
    is( $seq, 'TTCTTCTGGG', 'Negative start Ensembl subsequence' );
    $seq = $analysis->get_subsequence( '1', -1, -1, 1 );
    is( length $seq, 1,   'Negative start and end Ensembl subsequence length' );
    is( $seq,        'T', 'Negative start and end Ensembl subsequence' );
    $seq = $analysis->get_subsequence( '1', 1_000_000_001, 1_000_000_010, 1 );
    is( length $seq, 10, 'Large start and end Ensembl subsequence length' );
    is( $seq, 'NNNNNNNNNN', 'Large start and end Ensembl subsequence' );
}

# Check if Ensembl is reachable
sub is_ensembl_reachable {
    my $handle = IO::Socket::INET->new(
        PeerAddr => 'ensembldb.ensembl.org:5306',
        Timeout  => 1,
        Proto    => 'tcp',
    );

    if ( defined $handle && $handle ) {
        $handle->close();
        return 1;
    }
    else {
        return 0;
    }
}
