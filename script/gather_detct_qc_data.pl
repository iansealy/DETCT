#!/usr/bin/env perl

# PODNAME: gather_detct_qc_data.pl
# ABSTRACT: Gather DETCT data for QC

## Author         : is1
## Maintainer     : is1
## Created        : 2016-05-05
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
use DBI;
use Sort::Naturally;
use File::Spec;
use File::Find;
use DETCT::Misc qw( printf_or_die );
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Constants
Readonly our $HOST => 'seqw-db';
Readonly our $PORT => '3379';
Readonly our $USER => 'warehouse_ro';
Readonly our $PASS => q{};
Readonly our $NAME => 'sequencescape_warehouse';

Readonly our $SPIKE_PREFIX => 'ERCC';

Readonly our $OUTPUT_NAMES => join "\t",
  qw( sample transcriptome_insert_size genome_insert_size %dup estimated_library_size unmapped_reads genomic_region_counts filtered_genomic_region_counts spike_counts );

# Default options
my $input_file;
my $filtered_input_file;
my @expts;
my $bam_dir = '/lustre/scratch110/sanger/is1/detct-bam';
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my ( $genomic_count, $spike_count ) = get_counts($input_file);
my ($filtered_genomic_count) = get_counts($filtered_input_file);

# Connect to database
my $dsn = "dbi:mysql:database=$NAME;host=$HOST;port=$PORT";
my $dbh = DBI->connect( $dsn, $USER, $PASS );

printf_or_die( "%s\n", $OUTPUT_NAMES );

foreach my $expt (@expts) {
    my $expt_q = $dbh->quote($expt);
    $expt_q =~ s/'\z/%'/xms;
    my $ary_ref = $dbh->selectall_arrayref(
        <<"SQL"
        SELECT DISTINCT name, supplier_name, tag_sequence
        FROM   npg_plex_information npi, current_samples cs
        WHERE  npi.sample_id = cs.internal_id
        AND    (name LIKE $expt_q OR supplier_name LIKE $expt_q)
        AND    description LIKE '%polyT%'
        AND    tag_sequence IS NOT NULL
SQL
    );
    my %sample_for;
    my @samples;
    foreach ( @{$ary_ref} ) {
        my ( $name, $supplier_name, $tag_sequence ) = @{$_};
        if ( $name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $name;
            $sample_for{$tag_sequence} = $name;
        }
        elsif ( $supplier_name =~ m/\A $expt [^[:alpha:]\d] /xms ) {
            push @samples, $supplier_name;
            $sample_for{$tag_sequence} = $supplier_name;
        }
    }

    my @stats_files = get_stats_files( $bam_dir, $expt );
    ## no critic (ProhibitMagicNumbers)
    my %transcriptome_insert_size =
      get_file_data( \@stats_files, 'insertsize-transcriptome-stats.txt',
        4, \%sample_for );
    my %genome_insert_size =
      get_file_data( \@stats_files, 'insertsize-stats.txt', 4, \%sample_for );
    my %dup =
      get_file_data( \@stats_files, 'markdup-stats.txt', 7, \%sample_for );
    my %estimated_library_size =
      get_file_data( \@stats_files, 'markdup-stats.txt', 8, \%sample_for );
    my %unmapped_reads =
      get_file_data( \@stats_files, 'markdup-stats.txt', 3, \%sample_for );
    ## use critic

    foreach my $sample ( sort { ncmp( $a, $b ) } @samples ) {
        my @values = ($sample);
        push @values, $transcriptome_insert_size{$sample};
        push @values, $genome_insert_size{$sample};
        push @values, $dup{$sample};
        push @values, $estimated_library_size{$sample};
        push @values, $unmapped_reads{$sample};
        push @values, $genomic_count->{$sample};
        push @values, $filtered_genomic_count->{$sample};
        push @values, $spike_count->{$sample};
        @values = map { !defined $_ ? q{} : $_ } @values;
        printf_or_die( "%s\n", join "\t", @values );
    }
}

# Get counts from input file
sub get_counts {
    my ($file) = @_;

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }
    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)

    # Get input headings
    my $header = <$fh>;
    $header =~ s/\A [#]//xms;
    my @headings = DETCT::Misc::Output::parse_line( $header, $extension );

    # Get columns
    my @sample_cols;
    my %col_to_sample;
    my $col = -1;               ## no critic (ProhibitMagicNumbers)
    foreach my $heading (@headings) {
        $col++;
        if ( $heading =~ m/\A (\S+) \s count \z/xms ) {
            $col_to_sample{$col} = $1;
            push @sample_cols, $col;
        }
    }

    # Get counts
    my %genomic_count = map { $_ => 0 } values %col_to_sample;
    my %spike_count   = map { $_ => 0 } values %col_to_sample;
    while ( my $line = <$fh> ) {
        my @values = DETCT::Misc::Output::parse_line( $line, $extension );
        foreach my $col (@sample_cols) {
            if ( $values[0] =~ m/\A $SPIKE_PREFIX/xms ) {
                $spike_count{ $col_to_sample{$col} } += $values[$col];
            }
            else {
                $genomic_count{ $col_to_sample{$col} } += $values[$col];
            }
        }
    }

    close $fh;

    return \%genomic_count, \%spike_count;
}

# Get all possible stats files
sub get_stats_files {
    my ( $root_dir, $expt ) = @_;

    my @files;
    opendir my ($dh), $root_dir;
    while ( my $dir = readdir $dh ) {
        next if $dir !~ m/ $expt (?:[^[:alpha:]\d]|\z) /xms;
        find(
            sub { -f && m/stats.txt \z/xms && push @files, $File::Find::name },
            File::Spec->catfile( $root_dir, $dir )
        );
    }
    closedir $dh;

    return @files;
}

# Get data from file
sub get_file_data {
    my ( $files, $file_name, $col, $sample_for ) = @_;

    my %data_for = map { $_ => q{} } values %{$sample_for};

    my ($file) = grep { m/$file_name \z/xms } @{$files};
    return %data_for if !$file;

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        chomp $line;
        my @fields = split /\s+/xms, $line;
        my $sample;
        foreach my $tag_sequence ( keys %{$sample_for} ) {
            if ( $fields[0] =~ m/$tag_sequence/xms ) {
                $sample = $sample_for->{$tag_sequence};
                last;
            }
        }
        if ($sample) {
            $fields[$col] =~ s/, \z//xms;
            $data_for{$sample} = $fields[$col];
        }
    }

    close $fh;

    return %data_for;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'input_file=s'          => \$input_file,
        'filtered_input_file=s' => \$filtered_input_file,
        'expts=s@{1,}'          => \@expts,
        'bam_dir=s'             => \$bam_dir,
        'help'                  => \$help,
        'man'                   => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$input_file ) {
        pod2usage("--input_file must be specified\n");
    }
    if ( !$filtered_input_file ) {
        pod2usage("--filtered_input_file must be specified\n");
    }
    if ( !@expts ) {
        pod2usage("--expts must be specified\n");
    }

    return;
}

=head1 USAGE

    gather_detct_qc_data.pl
        [--input_file file]
        [--filtered_input_file file]
        [--expts prefixes]
        [--bam_dir dir]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--input_file FILE>

Differential expression pipeline output file (e.g. all.tsv).

=item B<--filtered_input_file FILE>

Filtered differential expression pipeline output file (e.g. all.filtered.tsv).

=item B<--expts PREFIXES>

Experiment prefixes (e.g. zmp_ph80).

=item B<--bam_dir DIR>

Top level directory containing all BAM files and QC files.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
