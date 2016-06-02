#!/usr/bin/env perl

# PODNAME: convert_to_biolayout.pl
# ABSTRACT: Convert DETCT output into BioLayout Express3D format

## Author         : is1
## Maintainer     : is1
## Created        : 2016-04-22
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
use DETCT::Misc qw( printf_or_die );
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Constants
Readonly our @FIELDS_BEFORE_LFC => ( 0 .. 8 );
Readonly our @FIELDS_AFTER_LFC  => ( 9 .. 14 );

# Default options
my $input_file;
my $samples_file;
my $readable_samples;
my $config_file;
my @metadata_files;
my $sample_regexp;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my ( $sample_cols, $lfc_cols ) = output_header(
    $input_file,  $samples_file,  $readable_samples,
    $config_file, $sample_regexp, @metadata_files
);
output_regions( $input_file, $sample_cols, $lfc_cols );

# Output header
sub output_header {    ## no critic (ProhibitManyArgs)
    ## no critic (ProhibitReusedNames)
    my (
        $input_file,  $samples_file,  $readable_samples,
        $config_file, $sample_regexp, @metadata_files
    ) = @_;
    ## use critic

    my ($extension) = $input_file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $input_file;
    }

    # Get input headings
    open my $input_fh, '<', $input_file;
    my $header = <$input_fh>;
    $header =~ s/\A [#]//xms;
    my @headings = DETCT::Misc::Output::parse_line( $header, $extension );
    close $input_fh;

    my @sample_cols;
    my @lfc_cols;
    my %sample_to_col;
    my $col = 0;
    my $i   = -1;    ## no critic (ProhibitMagicNumbers)
    foreach my $heading (@headings) {
        $i++;
        if ( $heading =~ m/normalised \s count \z/xms ) {
            $heading =~ s/\s normalised \s count \z//xms;
            if ( !$sample_regexp || $heading =~ m/$sample_regexp/xms ) {
                push @sample_cols, $i;
                $sample_to_col{$heading} = ++$col;
            }
        }
        if ( $heading =~ m/\A Log2 \s fold \s change/xms ) {
            push @lfc_cols, $i;
        }
    }

    my @output_headings = ('ID');
    push @output_headings, @headings[@FIELDS_BEFORE_LFC];
    push @output_headings, @headings[@lfc_cols];
    push @output_headings, @headings[@FIELDS_AFTER_LFC];
    push @output_headings, @headings[@sample_cols];

    my @tmp_headings;
    foreach my $heading (@output_headings) {
        $heading =~ s/\s normalised \s count \z//xms;
        push @tmp_headings, $heading;
    }
    @output_headings = @tmp_headings;

    my @all_output_headings = ( \@output_headings );

    my $cols_before_samples =
      scalar @FIELDS_BEFORE_LFC + scalar @lfc_cols + scalar @FIELDS_AFTER_LFC;
    if ($samples_file) {
        push @all_output_headings,
          output_samples_header( $samples_file, $cols_before_samples,
            \%sample_to_col );
    }
    if (@metadata_files) {
        push @all_output_headings,
          output_metadata_header(
            $cols_before_samples, $config_file,
            \%sample_to_col,      @metadata_files
          );
    }

    if ($readable_samples) {

        # Swap first two sample headings
        ( $all_output_headings[0], $all_output_headings[1] ) =
          ( $all_output_headings[1], $all_output_headings[0] );
        @tmp_headings =
          @{ $all_output_headings[0] }[ 0 .. $cols_before_samples ];
        @{ $all_output_headings[0] }[ 0 .. $cols_before_samples ] =
          @{ $all_output_headings[1] }[ 0 .. $cols_before_samples ];
        @{ $all_output_headings[1] }[ 0 .. $cols_before_samples ] =
          @tmp_headings;
    }

    foreach my $output_headings (@all_output_headings) {
        printf "%s\n", join "\t", tidy_output( @{$output_headings} );
    }

    return \@sample_cols, \@lfc_cols;
}

# Output samples header
sub output_samples_header {
    ## no critic (ProhibitReusedNames)
    my ( $samples_file, $cols_before_samples, $sample_to_col ) = @_;
    ## use critic

    open my $samples_fh, '<', $samples_file;    ## no critic (RequireBriefOpen)
    my $header = <$samples_fh>;
    chomp $header;
    if ( $header =~ m/\A \s/xms ) {
        $header =~ s/\A \s+//xms;
    }
    else {
        $header =~ s/\A \S+ \s+//xms;
    }
    my @columns = split /\s+/xms, $header;
    close $samples_fh;

    # Get all sample headings
    my @all_sample_headings;
    foreach my $col ( 1 .. scalar @columns ) {
        my @sample_headings =
          ( $columns[ $col - 1 ], (q{}) x $cols_before_samples );
        open my $samples_fh, '<', $samples_file;
        $header = <$samples_fh>;
        while ( my $line = <$samples_fh> ) {
            chomp $line;
            my @fields = split /\s+/xms, $line;
            next if !exists $sample_to_col->{ $fields[0] };
            push @sample_headings, $fields[$col];
        }
        close $samples_fh;
        push @all_sample_headings, \@sample_headings;
    }

    # Concatenate all factors if more than one
    if ( scalar @all_sample_headings > 1 ) {
        my @concat_headings = ( join q{_}, @columns );
        push @concat_headings, (q{}) x scalar $cols_before_samples;
        foreach my $sample_idx (
            scalar $cols_before_samples +
            1 .. scalar @{ $all_sample_headings[0] } -
            1 )
        {
            my @values;
            foreach my $factor_idx ( 0 .. scalar @columns - 1 ) {
                push @values, $all_sample_headings[$factor_idx][$sample_idx];
            }
            push @concat_headings, join q{_}, @values;
        }
        unshift @all_sample_headings, \@concat_headings;
    }

    # Convert to unique readable sample names
    my @headings = @{ $all_sample_headings[0] };
    my @readable = ('Sample');
    shift @headings;
    push @readable, (q{}) x scalar $cols_before_samples;
    my %heading_ordinal;
    foreach my $heading (@headings) {
        next if !$heading;
        push @readable, $heading . q{_} . ++$heading_ordinal{$heading};
    }
    unshift @all_sample_headings, \@readable;

    return @all_sample_headings;
}

# Output metadata header
sub output_metadata_header {
    ## no critic (ProhibitReusedNames)
    my ( $cols_before_samples, $config_file, $sample_to_col, @metadata_files )
      = @_;
    ## use critic

    my %round;
    my %skip;
    if ($config_file) {
        open my $fh, '<', $config_file;    ## no critic (RequireBriefOpen)
        while ( my $line = <$fh> ) {
            chomp $line;
            my ( $type, $value ) = split /\s+/xms, $line;
            next if !$value;
            if ( $value eq q{X} ) {
                $skip{$type} = 1;
            }
            elsif ( $value =~ m/\A [\d.]+ \z/xms ) {
                $round{$type} = $value;
            }
        }
        close $fh;
    }

    my @all_metadata_headings;
    foreach my $file (@metadata_files) {
        open my $fh, '<', $file;
        my $header = <$fh>;
        chomp $header;
        my @columns = split /\t/xms, $header;
        shift @columns;    # Ignore sample name column
        close $fh;

        # Get all metadata headings
        foreach my $col ( 1 .. scalar @columns ) {
            next if exists $skip{ $columns[ $col - 1 ] };
            my @headings =
              ( $columns[ $col - 1 ], (q{}) x $cols_before_samples );
            open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
            $header = <$fh>;
            while ( my $line = <$fh> ) {
                chomp $line;
                my @fields = split /\t/xms, $line;
                my $sample = $fields[0];
                next if !exists $sample_to_col->{$sample};
                if ( $fields[$col] && exists $round{ $columns[ $col - 1 ] } ) {
                    my $round_to = $round{ $columns[ $col - 1 ] };
                    $fields[$col] =
                      int( ( $fields[$col] + $round_to / 2 ) / $round_to ) *
                      $round_to;
                }
                $headings[ $sample_to_col->{$sample} + $cols_before_samples ] =
                  $fields[$col];
            }
            close $fh;
            push @all_metadata_headings, \@headings;
        }
    }

    return @all_metadata_headings;
}

# Output regions of interest
sub output_regions {
    ## no critic (ProhibitReusedNames)
    my ( $file, $sample_cols, $lfc_cols ) = @_;
    ## use critic

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get ID by joining gene name, chr, start, end, 3' end position & strand
        ## no critic (ProhibitMagicNumbers)
        $fields[13] = $fields[13] || q{};
        my $id = join q{:}, @fields[ 13, 0, 1, 2, 4 ];
        ## use critic

        my @output_fields = ($id);
        push @output_fields, @fields[@FIELDS_BEFORE_LFC];
        push @output_fields, @fields[ @{$lfc_cols} ];
        push @output_fields, @fields[@FIELDS_AFTER_LFC];
        push @output_fields, @fields[ @{$sample_cols} ];

        printf_or_die( "%s\n", join "\t", tidy_output(@output_fields) );
    }
    close $fh;

    return;
}

sub tidy_output {
    my (@fields) = @_;

    @fields = map { defined $_ ? $_       : q{-} } @fields;
    @fields = map { /\s/xms    ? qq{"$_"} : $_ } @fields;
    @fields = map {
        /\A ([-]?[\d.]+)(e[-]?\d+) \z/xms
          ? ( sprintf '%.8f', $1 ) . uc $2
          : $_
    } @fields;
    @fields =
      map { /\A [-]?\d+[.]\d+ \z/xms ? sprintf '%.8f', $_ : $_ } @fields;

    return @fields;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'input_file=s'         => \$input_file,
        'samples_file=s'       => \$samples_file,
        'readable_samples'     => \$readable_samples,
        'config_file=s'        => \$config_file,
        'metadata_files=s@{,}' => \@metadata_files,
        'sample_regexp=s'      => \$sample_regexp,
        'help'                 => \$help,
        'man'                  => \$man,
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
    if ( $readable_samples && !$samples_file ) {
        pod2usage("--readable_samples specified without --samples_file\n");
    }

    return;
}

=head1 USAGE

    convert_to_biolayout.pl
        [--input_file file]
        [--samples_file file]
        [--readable_samples]
        [--config_file file]
        [--metadata_files files]
        [--sample_regexp regexp]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--input_file FILE>

Differential expression pipeline output file (e.g. all.tsv).

=item B<--samples_file FILE>

Differential expression pipeline DESeq2 samples file (e.g. samples.txt).

=item B<--readable_samples>

Force sample names to be readable.

=item B<--config_file FILE>

Metadata config file. One data type per row with name in first column and either
an integer (specifying how to round values) or X (indicating types to ignore) in
second column.

=item B<--metadata_files FILES>

QC metadata files (e.g. QC_all_stages_lab.tsv).

=item B<--sample_regexp REGEXP>

Regular expression for matching sample names.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
