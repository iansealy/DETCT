## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Misc::Output;
## use critic

# ABSTRACT: Miscellaneous functions for outputting data

## Author         : is1
## Maintainer     : is1
## Created        : 2012-11-25
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
use File::Slurp;
use File::Spec;
use File::Path qw( make_path );
use Sort::Naturally;
use Text::CSV;
use List::MoreUtils qw( uniq all );
use DETCT::Misc qw( write_or_die );
use DETCT::Misc::R;

use base qw( Exporter );
our @EXPORT_OK = qw(
  dump_as_table
  dump_duplication_metrics
  parse_table
  convert_table_for_deseq
);

=head1 SYNOPSIS

    # Brief code examples

=cut

# Constants

# Types
Readonly our $STRING => 1;
Readonly our $INT    => 2;
Readonly our $FLOAT  => 3;

# Output formats
Readonly our @FORMATS => qw( csv tsv html );

=func dump_as_table

  Usage       : DETCT::Misc::Output::dump_as_table( {
                    analysis => $analysis,
                    dir      => '.',
                    regions  => $regions_hash_ref,
                } );
  Purpose     : Dump regions in tabular format
  Returns     : undef
  Parameters  : Hashref {
                    analysis => DETCT::Analysis,
                    dir      => String (the working directory),
                    regions  => Arrayref (of regions),
                }
  Throws      : If analysis is missing
                If regions are missing
                If directory is missing
  Comments    : None

=cut

sub dump_as_table {
    my ($arg_ref) = @_;

    confess 'No analysis specified'  if !defined $arg_ref->{analysis};
    confess 'No directory specified' if !defined $arg_ref->{dir};
    confess 'No regions specified'   if !defined $arg_ref->{regions};

    # Get conditions and groups
    my @samples = @{ $arg_ref->{analysis}->get_all_samples() };
    my @conditions = uniq( nsort( map { $_->condition } @samples ) );
    my @groups =
      uniq( nsort( grep { defined $_ } map { $_->group } @samples ) );

    # Get regions sorted by p value then location
    my $regions = sort_regions( $arg_ref->{regions} );

    # Get genebuild version
    my $genebuild_version;
    foreach my $region ( @{$regions} ) {
        ## no critic (ProhibitMagicNumbers)
        ($genebuild_version) = ( sort keys %{ $region->[15] } )[-1];   # Highest
        ## use critic
        last if $genebuild_version;
    }

    # Get definition for all columns (which determines formatting)
    my $definition = get_definition( $arg_ref->{analysis}->ensembl_species,
        $genebuild_version, \@samples, \@conditions, \@groups );

    # Make sure working directory exists
    if ( !-d $arg_ref->{dir} ) {
        make_path( $arg_ref->{dir} );
    }

    # Open filehandles and begin all output files
    my $fh = begin_output( $arg_ref->{dir}, $definition );

    foreach my $region ( @{$regions} ) {
        my @row;

        # Region
        my $seq_name = $region->[0];
        my $start    = $region->[1];
        my $end      = $region->[2];
        push @row, [ $seq_name, [ $seq_name, $start, $end, $seq_name ] ];
        push @row, [ $start,    [ $seq_name, $start, $end, $start ] ];
        push @row, [ $end,      [ $seq_name, $start, $end, $end ] ];

        # 3' end
        ## no critic (ProhibitMagicNumbers)
        my $tpe_seq_name   = $region->[5];
        my $tpe_pos        = $region->[6];
        my $tpe_strand     = $region->[7];
        my $tpe_read_count = $region->[8];
        ## use critic
        push @row,
          [ $tpe_pos, [ $tpe_seq_name, $tpe_pos, $tpe_pos, $tpe_pos ] ];
        push @row, [$tpe_strand];
        push @row, [$tpe_read_count];

        # p values
        ## no critic (ProhibitMagicNumbers)
        my $pval = $region->[11];
        my $padj = $region->[12];
        ## use critic
        push @row, [$pval];
        push @row, [$padj];

        # Gene details
        ## no critic (ProhibitMagicNumbers)
        my %gene = %{ $region->[15] };
        my ($genebuild) = ( sort keys %gene )[-1];    # Highest
        ## use critic
        my @distance;
        my ( @gene_stable_id, @gene_stable_id_to_link );
        my @gene_biotype;
        my ( @transcript_stable_id, @transcript_stable_id_to_link );
        my @transcript_biotype;
        my ( @name, @name_to_link );
        my @description;

        if ($genebuild) {
            foreach my $gene ( @{ $gene{$genebuild} } ) {
                my ( $gene_stable_id, $name, $description, $gene_biotype,
                    $distance, $transcripts )
                  = @{$gene};
                push @distance,       $distance;
                push @gene_stable_id, $gene_stable_id;
                push @gene_stable_id_to_link,
                  [ $gene_stable_id, $gene_stable_id ];
                push @gene_biotype, $gene_biotype;
                foreach my $transcript ( @{$transcripts} ) {
                    my ( $transcript_stable_id, $transcript_biotype ) =
                      @{$transcript};
                    push @transcript_stable_id, $transcript_stable_id;
                    push @transcript_stable_id_to_link,
                      [ $transcript_stable_id, $transcript_stable_id ];
                    push @transcript_biotype, $transcript_biotype;
                }
                push @name, $name;
                push @name_to_link, [ $gene_stable_id, $name ];
                push @description, $description;
            }
        }
        push @row, [ \@distance ];
        push @row, [ \@gene_stable_id, [@gene_stable_id_to_link] ];
        push @row, [ \@gene_biotype ];
        push @row, [ \@transcript_stable_id, [@transcript_stable_id_to_link] ];
        push @row, [ \@transcript_biotype ];
        push @row, [ \@name, [@name_to_link] ];
        push @row, [ \@description ];

        # Counts and normalised counts
        ## no critic (ProhibitMagicNumbers)
        my @counts            = map { [$_] } @{ $region->[9] };
        my @normalised_counts = map { [$_] } @{ $region->[10] };
        ## use critic
        foreach my $count (@counts) {
            push @row, [$count];
        }
        foreach my $normalised_count (@normalised_counts) {
            push @row, [$normalised_count];
        }

        # Condition fold changes
        if ( scalar @conditions == 2 ) {
            ## no critic (ProhibitMagicNumbers)
            my ( $condition_fold_change, $log2_condition_fold_change ) =
              @{ $region->[13] };
            ## use critic
            push @row, [$log2_condition_fold_change];
        }

        # Group fold changes
        if ( scalar @conditions == 2 && scalar @groups > 1 ) {
            ## no critic (ProhibitMagicNumbers)
            my @group_fold_changes = @{ $region->[14] };
            ## use critic
            foreach my $group_fold_change (@group_fold_changes) {
                my (
                    $condition_group_fold_change,
                    $log2_condition_group_fold_change
                ) = @{$group_fold_change};
                push @row, [$log2_condition_group_fold_change];
            }
        }

        # Dump row in all required formats
        my @levels = qw( all );
        if (   defined $padj
            && $padj ne 'NA'
            && $padj < $arg_ref->{analysis}->output_sig_level )
        {
            push @levels, 'sig';
        }
        dump_output( \@levels, $fh, $definition, \@row );
    }

    # End all output files and close filehandles
    end_output($fh);

    return;
}

=func sort_regions

  Usage       : $regions = sort_regions( $regions );
  Purpose     : Sort regions by p value then location
  Returns     : Arrayref (of regions)
  Parameters  : Arrayref of regions
  Throws      : No exceptions
  Comments    : None

=cut

sub sort_regions {
    my ($regions) = @_;

    # Separate regions with no p value from rest
    my @regions_with_pval;
    my @regions_no_pval;
    foreach my $region ( @{$regions} ) {
        ## no critic (ProhibitMagicNumbers)
        if ( defined $region->[11] && $region->[11] ne 'NA' ) {
            ## use critic
            push @regions_with_pval, $region;
        }
        else {
            push @regions_no_pval, $region;
        }
    }

    # Sort by adjusted p value and p value then regions without p value
    ## no critic (ProhibitMagicNumbers)
    my @regions =
      sort { $a->[12] <=> $b->[12] || $a->[11] <=> $b->[11] }
      @regions_with_pval;
    ## use critic
    push @regions, @regions_no_pval;    # Sorted by location already

    return \@regions;
}

=func get_definition

  Usage       : $definition = get_definition($genebuild_version, $samples,
                    $conditions, $groups);
  Purpose     : Return the definitions for all columns of the table
  Returns     : Arrayref (of column definitions)
  Parameters  : String (Ensembl species)
                String (genebuild version)
                Arrayref of samples
                Arrayref of conditions
                Arrayref of groups
  Throws      : No exceptions
  Comments    : None

=cut

sub get_definition {
    my ( $species, $genebuild_version, $samples, $conditions, $groups ) = @_;

    # Ensembl links
    my $loc_link =
      $species
      ? qq{<a href="http://www.ensembl.org/$species/psychic?q=%s:%d-%d" target="_blank">%s</a>}
      : undef;
    my $gene_link =
      q{<a href="http://www.ensembl.org/id/%s" target="_blank">%s</a>};

    my @def;

    push @def, [ 'Chr',              $STRING, $loc_link, ];
    push @def, [ 'Region start',     $INT,    $loc_link, ];
    push @def, [ 'Region end',       $INT,    $loc_link, ];
    push @def, [ q{3' end position}, $INT,    $loc_link, ];
    push @def, [ q{3' end strand},   $INT, ];
    push @def, [ q{3' end read count},  $INT, ];
    push @def, [ 'p value',             $FLOAT, ];
    push @def, [ 'Adjusted p value',    $FLOAT, ];
    push @def, [ q{Distance to 3' end}, $INT, ];
    push @def,
      [ $genebuild_version . ' Ensembl Gene ID', $STRING, $gene_link, ];
    push @def, [ 'Gene type', $STRING, ];
    push @def,
      [ $genebuild_version . ' Ensembl Transcript ID', $STRING, $gene_link, ];
    push @def, [ 'Transcript type',  $STRING, ];
    push @def, [ 'Gene name',        $STRING, $gene_link, ];
    push @def, [ 'Gene description', $STRING, ];

    foreach my $sample ( @{$samples} ) {
        push @def, [ $sample->name . ' count', $INT ];
    }
    foreach my $sample ( @{$samples} ) {
        push @def, [ $sample->name . ' normalised count', $FLOAT ];
    }

    if ( scalar @{$conditions} == 2 ) {
        my $heading = sprintf 'Log2 fold change (%s/%s)', $conditions->[0],
          $conditions->[1];
        push @def, [ $heading, $FLOAT ];
    }

    if ( scalar @{$conditions} == 2 && scalar @{$groups} > 1 ) {
        foreach my $group ( @{$groups} ) {
            my $heading = sprintf 'Log2 fold change (%s/%s) for group %s',
              $conditions->[0], $conditions->[1], $group;
            push @def, [ $heading, $FLOAT ];
        }
    }

    return \@def;
}

=func begin_output

  Usage       : my $fh = begin_output( $dir, $defintion );
  Purpose     : Open filehandles and begin all output files
  Returns     : undef
  Parameters  : String (the directory)
                Arrayref (the definition)
  Throws      : No exceptions
  Comments    : None

=cut

sub begin_output {
    my ( $dir, $definition ) = @_;

    my %fh;
    foreach my $format (@FORMATS) {

        # Level determines whether output all regions or just significant ones
        foreach my $level (qw( all sig )) {
            my $file = File::Spec->catfile( $dir, $level . q{.} . $format );
            open my $fh, '>', $file;    ## no critic (RequireBriefOpen)
            $fh{$format}{$level} = $fh;
            my $begin_sub_name = 'begin_' . $format;
            my $sub_ref        = \&{$begin_sub_name};
            &{$sub_ref}( $fh, $definition );
        }
    }

    return \%fh;
}

=func end_output

  Usage       : end_output( $fh );
  Purpose     : End all output files and close filehandles
  Returns     : undef
  Parameters  : Hashref (of filehandles)
  Throws      : No exceptions
  Comments    : None

=cut

sub end_output {
    my ($fh) = @_;

    foreach my $format (@FORMATS) {
        foreach my $level (qw( all sig )) {
            my $end_sub_name = 'end_' . $format;
            my $sub_ref      = \&{$end_sub_name};
            &{$sub_ref}( $fh->{$format}{$level} );
            close $fh->{$format}{$level};
        }
    }

    return;
}

=func dump_output

  Usage       : dump_output( $levels, $fh, $definition, $row );
  Purpose     : Dump row in all required formats at all required levels
  Returns     : undef
  Parameters  : Arrayref (of levels)
                Hashref (of filehandles)
                Arrayref (the definition)
                Arrayref (of row data)
  Throws      : No exceptions
  Comments    : None

=cut

sub dump_output {
    my ( $levels, $fh, $definition, $row ) = @_;

    foreach my $format (@FORMATS) {
        foreach my $level ( @{$levels} ) {
            my $dump_sub_name = 'dump_' . $format;
            my $sub_ref       = \&{$dump_sub_name};
            &{$sub_ref}( $fh->{$format}{$level}, $definition, $row );
        }
    }

    return;
}

=func begin_csv

  Usage       : begin_csv( $fh, $defintion );
  Purpose     : Begin CSV table
  Returns     : undef
  Parameters  : Filehandle
                Arrayref (the definition)
  Throws      : No exceptions
  Comments    : None

=cut

sub begin_csv {
    my ( $fh, $definition ) = @_;

    my @headings;
    foreach my $column ( @{$definition} ) {
        my ($heading) = @{$column};
        $heading =~ s/"/""/xmsg;
        push @headings, q{"} . $heading . q{"};
    }
    write_or_die( $fh, ( join q{,}, @headings ), "\r\n" );

    return;
}

=func end_csv

  Usage       : end_csv( $fh );
  Purpose     : End CSV table
  Returns     : undef
  Parameters  : Filehandle
  Throws      : No exceptions
  Comments    : None

=cut

sub end_csv {
    return;
}

=func dump_csv

  Usage       : dump_csv( $fh, $definition, $row );
  Purpose     : Dump the data in a CSV table
  Returns     : undef
  Parameters  : Filehandle
                Arrayref (the defintion)
                Arrayref (the row data)
  Throws      : No exceptions
  Comments    : None

=cut

sub dump_csv {
    my ( $fh, $definition, $row ) = @_;

    my @output_row;
    my $i = 0;    # Index to definition
    foreach my $cell ( @{$row} ) {
        my $type = $definition->[$i]->[1];
        my ($data) = @{$cell};

        # Turn into a list of data, even if just one
        if ( ref $data ne 'ARRAY' ) {
            $data = [$data];
        }

        # Substitute default if undefined
        my @output_cell;
        foreach my $datum ( @{$data} ) {
            $datum = defined $datum ? $datum : q{};
            push @output_cell, $datum;
        }

        # Add default if necessary
        if ( !@output_cell ) {
            push @output_cell, q{};
        }

        my $output_cell = join q{,}, @output_cell;

        # Strings and lists need quoting
        if ( $type == $STRING || scalar @output_cell > 1 ) {
            $output_cell =~ s/"/""/xmsg;
            $output_cell = q{"} . $output_cell . q{"};
        }

        push @output_row, $output_cell;

        $i++;
    }
    write_or_die( $fh, ( join q{,}, @output_row ), "\r\n" );

    return;
}

=func begin_tsv

  Usage       : begin_tsv( $fh, $defintion );
  Purpose     : Begin TSV table
  Returns     : undef
  Parameters  : Filehandle
                Arrayref (the definition)
  Throws      : No exceptions
  Comments    : None

=cut

sub begin_tsv {
    my ( $fh, $definition ) = @_;

    my @headings = map { $_->[0] } @{$definition};
    write_or_die( $fh, q{#}, ( join "\t", @headings ), "\n" );

    return;
}

=func end_tsv

  Usage       : end_tsv( $fh );
  Purpose     : End TSV table
  Returns     : undef
  Parameters  : Filehandle
  Throws      : No exceptions
  Comments    : None

=cut

sub end_tsv {
    return;
}

=func dump_tsv

  Usage       : dump_tsv( $fh, $definition, $row );
  Purpose     : Dump the data in a TSV table
  Returns     : undef
  Parameters  : Filehandle
                Arrayref (the defintion)
                Arrayref (the row data)
  Throws      : No exceptions
  Comments    : None

=cut

sub dump_tsv {
    my ( $fh, $definition, $row ) = @_;

    my @output_row;
    my $i = 0;    # Index to definition
    foreach my $cell ( @{$row} ) {
        my $type = $definition->[$i]->[1];
        my ($data) = @{$cell};

        # Turn into a list of data, even if just one
        if ( ref $data ne 'ARRAY' ) {
            $data = [$data];
        }

        # Substitute default if undefined
        my @output_cell;
        foreach my $datum ( @{$data} ) {
            $datum = defined $datum && length $datum > 0 ? $datum : q{-};
            push @output_cell, $datum;
        }

        # Add default if necessary
        if ( !@output_cell ) {
            push @output_cell, q{-};
        }

        push @output_row, ( join q{,}, @output_cell );

        $i++;
    }
    write_or_die( $fh, ( join "\t", @output_row ), "\n" );

    return;
}

=func begin_html

  Usage       : begin_html( $fh, $defintion );
  Purpose     : Begin HTML table
  Returns     : undef
  Parameters  : Filehandle
                Arrayref (the definition)
  Throws      : No exceptions
  Comments    : None

=cut

sub begin_html {
    my ( $fh, $definition ) = @_;

    write_or_die( $fh, <<'HTML');
<!DOCTYPE html>
<html>
    <head>
        <title>DETCT</title>
        <link href="http://netdna.bootstrapcdn.com/twitter-bootstrap/2.2.1/css/bootstrap-combined.min.css" rel="stylesheet">
    </head>
    <body>
        <table class="table table-bordered table-hover table-condensed">
            <thead>
                <tr>
HTML

    foreach my $column ( @{$definition} ) {
        my ($heading) = @{$column};
        write_or_die( $fh, '<th>', $heading, '</th>', "\n" );
    }

    write_or_die( $fh, <<'HTML');
                </tr>
            </thead>
            <tbody>
HTML

    return;
}

=func end_html

  Usage       : end_html( $fh );
  Purpose     : End HTML table
  Returns     : undef
  Parameters  : Filehandle
  Throws      : No exceptions
  Comments    : None

=cut

sub end_html {
    my ($fh) = @_;

    write_or_die( $fh, <<'HTML');
            </tbody>
        </table>
    </body>
</html>
HTML

    return;
}

=func dump_html

  Usage       : dump_html( $fh, $definition, $row );
  Purpose     : Dump the data in an HTML table
  Returns     : undef
  Parameters  : Filehandle
                Arrayref (the defintion)
                Arrayref (the row data)
  Throws      : No exceptions
  Comments    : None

=cut

sub dump_html {
    my ( $fh, $definition, $row ) = @_;

    write_or_die( $fh, '<tr>', "\n" );

    my $i = 0;    # Index to definition
    foreach my $cell ( @{$row} ) {
        my ( undef, $type, $link ) = @{ $definition->[$i] };
        my ( $data, $data_to_link ) = @{$cell};

        # Turn into a list of data, even if just one
        if ( ref $data ne 'ARRAY' ) {
            $data         = [$data];
            $data_to_link = [$data_to_link];
        }

        write_or_die( $fh, '<td>' );

        my @data;
        my $j = 0;    # Index to each item when multiple items in one table cell
        foreach my $datum ( @{$data} ) {
            my $datum_to_link = $data_to_link->[$j];

            $datum = defined $datum ? $datum : q{};

            # Make a link if there's a link and all data for the link is defined
            if ( $link && $datum_to_link && all { defined $_ }
                @{$datum_to_link} )
            {
                $datum = sprintf $link, @{$datum_to_link};
            }

            push @data, $datum;

            $j++;
        }

        write_or_die( $fh, join '<br />', @data );

        write_or_die( $fh, '</td>', "\n" );

        $i++;
    }

    write_or_die( $fh, '</tr>', "\n" );

    return;
}

=func parse_table

  Usage       : $regions = DETCT::Misc::Output::parse_table( {
                    analysis     => $analysis,
                    table_file   => 'all.tsv',
                    table_format => 'tsv',
                } );
  Purpose     : Parse table of regions into final internal format
  Returns     : Arrayref [
                    Arrayref [
                        String (region sequence name),
                        Int (region start),
                        Int (region end),
                        undef,
                        undef,
                        String (3' end sequence name) or undef,
                        Int (3' end position) or undef,
                        Int (3' end strand) or undef,
                        Int (3' end read count) or undef,
                        Arrayref [
                            Int (count)
                            ...
                        ],
                        Arrayref [
                            Int (normalised count)
                            ...
                        ],
                        Int (p value) or undef,
                        Int (adjusted p value) or undef,
                        Arrayref [
                            Int (condition fold change) or undef,
                            Int (log2 condition fold change) or undef,
                        ],
                        Arrayref [
                            Arrayref [
                                Int (group fold change) or undef,
                                Int (log2 group fold change) or undef,
                            ],
                            ... (groups)
                        ],
                        Hashref {
                            String (genebuild version) => Arrayref [
                                Arrayref [
                                    String (gene stable id),
                                    String (gene name) or undef,
                                    String (gene description) or undef,
                                    String (gene biotype),
                                    Int (distance to 3' end),
                                    Arrayref [
                                        Arrayref [
                                            String (transcript stable id),
                                            String (transcript biotype),
                                        ],
                                        ... (transcripts)
                                    ],
                                ],
                                ... (genes)
                            ],
                        }
                    ],
                    ... (regions)
                ]
  Parameters  : Hashref {
                    analysis     => DETCT::Analysis,
                    table_file   => String (the table file),
                    table_format => String (the table format),
                }
  Throws      : If analysis is missing
                If table file is missing
                If table format is missing or invalid or HTML
  Comments    : None

=cut

sub parse_table {    ## no critic (ProhibitExcessComplexity)
    my ($arg_ref) = @_;

    confess 'No analysis specified'     if !defined $arg_ref->{analysis};
    confess 'No table file specified'   if !defined $arg_ref->{table_file};
    confess 'No table format specified' if !defined $arg_ref->{table_format};
    confess 'Cannot handle HTML format' if $arg_ref->{table_format} eq 'html';
    if ( all { $_ ne $arg_ref->{table_format} } @FORMATS ) {
        confess sprintf 'Invalid table format (%s) specified',
          $arg_ref->{table_format};
    }

    # Get conditions and groups
    my @samples = @{ $arg_ref->{analysis}->get_all_samples() };
    my @conditions = uniq( nsort( map { $_->condition } @samples ) );
    my @groups =
      uniq( nsort( grep { defined $_ } map { $_->group } @samples ) );

    # Get table
    my @rows = read_file( $arg_ref->{table_file} );

    # Get headings
    my $line = shift @rows;
    my @headings = parse_line( $line, $arg_ref->{table_format} );

    # Get genebuild version from Ensembl Gene ID heading
    ## no critic (ProhibitMagicNumbers)
    my ($genebuild_version) = $headings[9] =~ m/\A (e\d+) \s/xms;
    ## use critic

    # Get ordinal for each count column
    my ( %count_ordinal_for, %normalised_count_ordinal_for );
    ## no critic (ProhibitMagicNumbers)
    my $ordinal = 15;    # First count column
    ## use critic
    while ( defined $headings[$ordinal] ) {
        if ( $headings[$ordinal] =~ m/\A ([\w.-]+) \s count \z/xms ) {
            $count_ordinal_for{$1} = $ordinal;
        }
        elsif ( $headings[$ordinal] =~
            m/\A ([\w.-]+) \s normalised \s count \z/xms )
        {
            $normalised_count_ordinal_for{$1} = $ordinal;
        }
        $ordinal++;
    }

    # Parse rows into regions
    my @regions;
    foreach my $row (@rows) {
        my @row = parse_line( $row, $arg_ref->{table_format} );
        my @region;
        push @region, $row[0];    # Region sequence name
        push @region, $row[1];    # Region start
        push @region, $row[2];    # Region end
        push @region, undef;      # Region maximum read count
        push @region, undef;      # Region log probability sum
        ## no critic (ProhibitMagicNumbers)
        push @region, defined $row[3] ? $row[0] : undef;  # 3' end sequence name
        push @region, $row[3];                            # 3' end position
        push @region, $row[4];                            # 3' end strand
        push @region, $row[5];                            # 3' end read count
        ## use critic

        # Counts and normalised counts, plus store for fold change calculation
        my @counts;
        my @normalised_counts;
        my %counts_for_condition;
        my %counts_for_group_condition;
        foreach my $sample (@samples) {
            push @counts, $row[ $count_ordinal_for{ $sample->name } ];
            my $normalised_count =
              $row[ $normalised_count_ordinal_for{ $sample->name } ];
            push @normalised_counts, $normalised_count;
            if ( scalar @conditions == 2 ) {
                push @{ $counts_for_condition{ $sample->condition } },
                  $normalised_count;
            }
            if ( scalar @conditions == 2 && scalar @groups > 1 ) {
                push @{ $counts_for_group_condition{ $sample->group }
                      { $sample->condition } }, $normalised_count;
            }
        }
        push @region, \@counts;
        push @region, \@normalised_counts;

        ## no critic (ProhibitMagicNumbers)
        push @region, $row[6];    # p value
        push @region, $row[7];    # Adjusted p value
        ## use critic

        # Condition fold change
        push @region,
          DETCT::Misc::R::calc_condition_fold_change( \@conditions,
            \%counts_for_condition );

        # Group fold changes
        push @region,
          DETCT::Misc::R::calc_group_fold_changes( \@conditions, \@groups,
            \%counts_for_group_condition );

        # Skip to next region if no gene annotation
        if ( !defined $row[9] ) {    ## no critic (ProhibitMagicNumbers)
            push @region, {};
            push @regions, \@region;
            next;
        }

        # Gene annotation
        # Not guaranteed to be identical to original data (e.g. no way of
        # knowing which transcripts are associated with which genes), but data
        # will be output correctly as a table (i.e. can roundtrip)
        ## no critic (ProhibitMagicNumbers)
        my @distance             = split /,/xms, $row[8];     # 3' end distance
        my @gene_stable_id       = split /,/xms, $row[9];     # Gene ID
        my @gene_biotype         = split /,/xms, $row[10];    # Gene type
        my @transcript_stable_id = split /,/xms, $row[11];    # Transcript ID
        my @transcript_biotype   = split /,/xms, $row[12];    # Transcript type

        # Gene name (can be missing)
        my @name = defined $row[13] ? split /,\S/xms, $row[13] : ();

        # Gene description (can be missing)
        my @description = defined $row[14] ? split /,\S/xms, $row[14] : ();
        ## use critic

        my @genes;
        while (@gene_stable_id) {
            if ( scalar @gene_stable_id == 1 ) {

                # Last gene, so merge remaining data if necessary
                @name        = ( join q{,}, @name );
                @description = ( join q{,}, @description );
            }
            my @gene;
            push @gene, shift @gene_stable_id;
            push @gene, shift @name;
            push @gene, shift @description;
            push @gene, shift @gene_biotype;
            push @gene, shift @distance;
            if (@gene_stable_id) {

                # Not last gene, so just add one transcript
                push @gene,
                  [ [ shift @transcript_stable_id, shift @transcript_biotype ]
                  ];
            }
            else {
                # Last gene, so add remaining transcripts
                while (@transcript_stable_id) {
                    push @gene,
                      [
                        [
                            shift @transcript_stable_id,
                            shift @transcript_biotype
                        ]
                      ];
                }
            }

            push @genes, \@gene;
        }
        push @region, { $genebuild_version => \@genes };

        push @regions, \@region;
    }

    return \@regions;
}

=func parse_line

  Usage       : my @row = parse_line( $line, 'tsv' );
  Purpose     : Parse line of TSV or CSV into fields
  Returns     : Array
  Parameters  : String (the line from a table)
                String (the table format)
  Throws      : If CSV can't be parsed
  Comments    : None

=cut

sub parse_line {
    my ( $line, $format ) = @_;

    $line =~ s/[\r\n]+ \z//xms;

    my @row;

    if ( $format eq 'tsv' ) {
        @row = split /\t/xms, $line;
        @row = map { $_ eq q{-} || length $_ == 0 ? undef : $_ } @row;
    }
    elsif ( $format eq 'csv' ) {
        confess "Can't parse CSV: $line" if !_csv()->parse($line);
        @row = _csv()->fields();
    }

    return @row;
}

=func convert_table_for_deseq

  Usage       : $regions = DETCT::Misc::Output::convert_table_for_deseq( {
                    regions  => $regions_ary_ref,
                } );
  Purpose     : Convert final internal format to format for DESeq
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
                    regions  => Arrayref (of regions),
                }
  Throws      : If regions are missing
  Comments    : None

=cut

sub convert_table_for_deseq {
    my ($arg_ref) = @_;

    confess 'No regions specified' if !defined $arg_ref->{regions};

    my %regions;

    foreach my $region ( @{ $arg_ref->{regions} } ) {
        my $seq_name = shift @{$region};
        @{$region} = @{$region}[ 0 .. 8 ];   ## no critic (ProhibitMagicNumbers)
        push @{ $regions{$seq_name} }, $region;
    }

    return \%regions;
}

# Cache Text::CSV object
{
    my $text_csv_object = Text::CSV->new( { empty_is_undef => 1 } );

    sub _csv {
        return $text_csv_object;
    }
}

=func dump_duplication_metrics

  Usage       : DETCT::Misc::Output::dump_duplication_metrics( {
                    output_file => 'markdup.tsv',
                    metrics     => \@metrics,
                } );
  Purpose     : Dump duplication metrics in tsv format
  Returns     : undef
  Parameters  : Hashref {
                    output_file => String (the output file),
                    metrics     => Arrayref (of metrics Hashrefs),
                }
  Throws      : If output file is missing
  Comments    : None

=cut

sub dump_duplication_metrics {
    my ($arg_ref) = @_;

    confess 'No output file specified' if !defined $arg_ref->{output_file};

    # Header
    my $metrics_output = q{#}
      . (
        join "\t",
        'Sample',
        'Mapped reads without mapped mate',
        'Mapped read pairs',
        'Mapped reads',
        'Unmapped reads',
        'Duplicate mapped reads without mapped mate',
        'Duplicate mapped read pairs',
        'Optical duplicate mapped read pairs',
        'Duplicate reads',
        'Duplication rate',
        'Estimated library size',
      ) . "\n";

    # Get metrics
    foreach my $output ( @{ $arg_ref->{metrics} } ) {
        $metrics_output .= (
            join "\t",
            $output->{sample_name},
            $output->{mapped_reads_without_mapped_mate},
            $output->{mapped_read_pairs},
            $output->{mapped_reads},
            $output->{unmapped_reads},
            $output->{duplicate_mapped_reads_without_mapped_mate},
            $output->{duplicate_mapped_read_pairs},
            $output->{optical_duplicate_mapped_read_pairs},
            $output->{duplicate_reads},
            $output->{duplication_rate},
            $output->{estimated_library_size},
        ) . "\n";
    }

    write_file( $arg_ref->{output_file}, $metrics_output );

    return;
}

1;
