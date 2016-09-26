#!/usr/bin/env perl

# PODNAME: compare_de_output.pl
# ABSTRACT: Compare output files from differential expression pipeline

## Author         : is1
## Maintainer     : is1
## Created        : 2015-02-26
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
use Statistics::Descriptive;
use DETCT::Misc qw( print_or_die printf_or_die );
use DETCT::Misc::Output;

=head1 DESCRIPTION


=head1 EXAMPLES


=cut

# Default options
my $file1;
my $file2;
my $sig_level = 0.05;    ## no critic (ProhibitMagicNumbers)
my ( $help, $man );

# Get and check command line options
get_and_check_options();

my $data1_for = get_data_by_region( $file1, $sig_level );
my $data2_for = get_data_by_region( $file2, $sig_level );

my $identical_regions =
  get_identical_regions( [ keys %{$data1_for} ], [ keys %{$data2_for} ] );

printf_or_die( "Total region count for %s:\t%d\n",
    $file1, scalar keys %{$data1_for} );
printf_or_die( "Total region count for %s:\t%d\n",
    $file2, scalar keys %{$data2_for} );
printf_or_die( "Identical regions:\t%d\n", scalar @{$identical_regions} );

my ( $stat1, $stat2, $count_for );

# 3' end position

$stat1 = summarise_field_for_all_data( $data1_for, 0 );
$stat2 = summarise_field_for_all_data( $data2_for, 0 );

print_or_die("\n===\n\n3' end position:\n");
printf_or_die( "Regions without 3' end position for %s:\t%d\n",
    $file1, $stat1->{'missing_count'} );
printf_or_die( "Regions without 3' end position for %s:\t%d\n",
    $file2, $stat2->{'missing_count'} );
printf_or_die( "Regions with 3' end position for %s:\t%d\n",
    $file1, $stat1->{'count'} );
printf_or_die( "Regions with 3' end position for %s:\t%d\n",
    $file2, $stat2->{'count'} );

$count_for = compare_numeric_field_for_data_subset( $data1_for, $data2_for,
    $identical_regions, 0 );

print_or_die("\n3' end position for identical regions:\n");
printf_or_die( "Regions without 3' end position for both:\t%d\n",
    $count_for->{'missing_value_both'} );
printf_or_die( "Regions without 3' end position just for %s:\t%d\n",
    $file1, $count_for->{'missing_value1'} );
printf_or_die( "Regions without 3' end position just for %s:\t%d\n",
    $file2, $count_for->{'missing_value2'} );
printf_or_die( "Regions with identical 3' end positions:\t%d\n",
    $count_for->{'identical_value'} );

# 3' end distance

$stat1 = summarise_field_for_all_data( $data1_for, 1 );
$stat2 = summarise_field_for_all_data( $data2_for, 1 );

print_or_die("\n===\n\n3' end distance:\n");
printf_or_die( "Regions without 3' end distance for %s:\t%d\n",
    $file1, $stat1->{'missing_count'} );
printf_or_die( "Regions without 3' end distance for %s:\t%d\n",
    $file2, $stat2->{'missing_count'} );
printf_or_die( "Regions with 3' end distance for %s:\t%d\n",
    $file1, $stat1->{'count'} );
printf_or_die( "Regions with 3' end distance for %s:\t%d\n",
    $file2, $stat2->{'count'} );
printf_or_die( "Mean 3' end distance for %s:\t%.1f\n",
    $file1, $stat1->{'mean'} );
printf_or_die( "Mean 3' end distance for %s:\t%.1f\n",
    $file2, $stat2->{'mean'} );
printf_or_die( "Standard deviation of 3' end distance for %s:\t%.2f\n",
    $file1, $stat1->{'sd'} );
printf_or_die( "Standard deviation of 3' end distance for %s:\t%.2f\n",
    $file2, $stat2->{'sd'} );
printf_or_die( "Median 3' end distance for %s:\t%.1f\n",
    $file1, $stat1->{'median'} );
printf_or_die( "Median 3' end distance for %s:\t%.1f\n",
    $file2, $stat2->{'median'} );
printf_or_die( "Modal 3' end distance for %s:\t%d\n",
    $file1, $stat1->{'mode'} );
printf_or_die( "Modal 3' end distance for %s:\t%d\n",
    $file2, $stat2->{'mode'} );

$count_for = compare_numeric_field_for_data_subset( $data1_for, $data2_for,
    $identical_regions, 1 );

print_or_die("\n3' end position for identical regions:\n");
printf_or_die( "Regions without 3' end distance for both:\t%d\n",
    $count_for->{'missing_value_both'} );
printf_or_die( "Regions without 3' end distance just for %s:\t%d\n",
    $file1, $count_for->{'missing_value1'} );
printf_or_die( "Regions without 3' end distance just for %s:\t%d\n",
    $file2, $count_for->{'missing_value2'} );
printf_or_die( "Regions with identical 3' end distance:\t%d\n",
    $count_for->{'identical_value'} );

# p-value

$stat1 = summarise_field_for_all_data( $data1_for, 2, 'NA' );
$stat2 = summarise_field_for_all_data( $data2_for, 2, 'NA' );

print_or_die("\n===\n\np-value:\n");
printf_or_die( "Regions without p-value for %s:\t%d\n",
    $file1, $stat1->{'missing_count'} );
printf_or_die( "Regions without p-value for %s:\t%d\n",
    $file2, $stat2->{'missing_count'} );
printf_or_die( "Regions with p-value for %s:\t%d\n",
    $file1, $stat1->{'count'} );
printf_or_die( "Regions with p-value for %s:\t%d\n",
    $file2, $stat2->{'count'} );

$count_for = compare_numeric_field_for_data_subset( $data1_for, $data2_for,
    $identical_regions, 2, 'NA' );

print_or_die("\np-value for identical regions:\n");
printf_or_die(
    "Regions without p-value for both:\t%d\n",
    $count_for->{'missing_value_both'}
);
printf_or_die( "Regions without p-value just for %s:\t%d\n",
    $file1, $count_for->{'missing_value1'} );
printf_or_die( "Regions without p-value just for %s:\t%d\n",
    $file2, $count_for->{'missing_value2'} );
printf_or_die( "Regions with identical p-values:\t%d\n",
    $count_for->{'identical_value'} );
printf_or_die( "Regions with lower p-values for %s than %s:\t%d\n",
    $file1, $file2, $count_for->{'lower_value'} );
printf_or_die( "Regions with higher p-values for %s than %s:\t%d\n",
    $file1, $file2, $count_for->{'higher_value'} );

# Adjusted p-value

## no critic (ProhibitMagicNumbers)
$stat1 = summarise_field_for_all_data( $data1_for, 3, 'NA' );
$stat2 = summarise_field_for_all_data( $data2_for, 3, 'NA' );
## use critic

print_or_die("\n===\n\nAdjusted p-value:\n");
printf_or_die( "Regions without adjusted p-value for %s:\t%d\n",
    $file1, $stat1->{'missing_count'} );
printf_or_die( "Regions without adjusted p-value for %s:\t%d\n",
    $file2, $stat2->{'missing_count'} );
printf_or_die( "Regions with adjusted p-value for %s:\t%d\n",
    $file1, $stat1->{'count'} );
printf_or_die( "Regions with adjusted p-value for %s:\t%d\n",
    $file2, $stat2->{'count'} );

$count_for = compare_numeric_field_for_data_subset( $data1_for, $data2_for,
    $identical_regions, 3, 'NA' );    ## no critic (ProhibitMagicNumbers)

print_or_die("\nAdjusted p-value for identical regions:\n");
printf_or_die( "Regions without adjusted p-value for both:\t%d\n",
    $count_for->{'missing_value_both'} );
printf_or_die( "Regions without adjusted p-value just for %s:\t%d\n",
    $file1, $count_for->{'missing_value1'} );
printf_or_die( "Regions without adjusted p-value just for %s:\t%d\n",
    $file2, $count_for->{'missing_value2'} );
printf_or_die( "Regions with identical adjusted p-values:\t%d\n",
    $count_for->{'identical_value'} );
printf_or_die( "Regions with lower adjusted p-values for %s than %s:\t%d\n",
    $file1, $file2, $count_for->{'lower_value'} );
printf_or_die( "Regions with higher adjusted p-values for %s than %s:\t%d\n",
    $file1, $file2, $count_for->{'higher_value'} );

# Significant adjusted p-value

## no critic (ProhibitMagicNumbers)
$stat1 = summarise_field_for_all_data( $data1_for, 4 );
$stat2 = summarise_field_for_all_data( $data2_for, 4 );
## use critic

print_or_die("\n===\n\nSignificant adjusted p-value:\n");
printf_or_die( "Regions without significant adjusted p-value for %s:\t%d\n",
    $file1, $stat1->{'missing_count'} );
printf_or_die( "Regions without significant adjusted p-value for %s:\t%d\n",
    $file2, $stat2->{'missing_count'} );
printf_or_die( "Regions with significant adjusted p-value for %s:\t%d\n",
    $file1, $stat1->{'count'} );
printf_or_die( "Regions with significant adjusted p-value for %s:\t%d\n",
    $file2, $stat2->{'count'} );

$count_for = compare_numeric_field_for_data_subset( $data1_for, $data2_for,
    $identical_regions, 4 );    ## no critic (ProhibitMagicNumbers)

print_or_die("\nSignificant adjusted p-value for identical regions:\n");
printf_or_die( "Regions without significant adjusted p-value for both:\t%d\n",
    $count_for->{'missing_value_both'} );
printf_or_die(
    "Regions without significant adjusted p-value just for %s:\t%d\n",
    $file1, $count_for->{'missing_value1'} );
printf_or_die(
    "Regions without significant adjusted p-value just for %s:\t%d\n",
    $file2, $count_for->{'missing_value2'} );
printf_or_die( "Regions with identical significant adjusted p-values:\t%d\n",
    $count_for->{'identical_value'} );
printf_or_die(
    "Regions with lower significant adjusted p-values for %s than %s:\t%d\n",
    $file1, $file2, $count_for->{'lower_value'} );
printf_or_die(
    "Regions with higher significant adjusted p-values for %s than %s:\t%d\n",
    $file1, $file2, $count_for->{'higher_value'} );

# Get data (currently just 3' end position, distance and p values) by region
sub get_data_by_region {
    my ( $file, $adjusted_p_value_sig_level ) = @_;

    my ($extension) = $file =~ m/[.] ([[:lower:]]{3}) \z/xms;
    if ( !$extension || ( $extension ne 'csv' && $extension ne 'tsv' ) ) {
        confess sprintf '%s is not .csv or .tsv file', $file;
    }

    my %data_for;

    open my $fh, '<', $file;    ## no critic (RequireBriefOpen)
    my $header = <$fh>;
    while ( my $line = <$fh> ) {
        my @fields = DETCT::Misc::Output::parse_line( $line, $extension );

        # Get region ID by joining chr, start, end and strand
        ## no critic (ProhibitMagicNumbers)
        my $region                   = join q{:}, @fields[ 0, 1, 2, 4 ];
        my $three_prime_end_pos      = $fields[3];
        my $three_prime_end_distance = $fields[8];
        my $p_value                  = $fields[6];
        my $adjusted_p_value         = $fields[7];
        ## use critic
        # If multiple 3' ends then take nearest one or most frequest one
        my $min_distance_index;
        if ( $three_prime_end_distance && $three_prime_end_distance =~ m/,/xms )
        {
            my $min_distance;
            my $index = -1;    ## no critic (ProhibitMagicNumbers)
            foreach my $distance ( split /,/xms, $three_prime_end_distance ) {
                $index++;
                next if $distance eq q{-} || !length $distance;
                if (  !defined $min_distance
                    || abs $distance < abs $min_distance )
                {
                    $min_distance       = $distance;
                    $min_distance_index = $index;
                }
            }
            $three_prime_end_distance = $min_distance;
        }
        if ( $three_prime_end_pos && $three_prime_end_pos =~ m/,/xms ) {
            my @pos = split /,/xms, $three_prime_end_pos;
            if ($min_distance_index) {
                $three_prime_end_pos =
                  $pos[ $min_distance_index % scalar @pos ];
            }
            else {
                $three_prime_end_pos = $pos[0];
            }
        }
        my $sig_adjusted_p_value =
          (      $adjusted_p_value ne 'NA'
              && $adjusted_p_value < $adjusted_p_value_sig_level )
          ? $adjusted_p_value
          : undef;
        $data_for{$region} = [
            $three_prime_end_pos, $three_prime_end_distance,
            $p_value,             $adjusted_p_value,
            $sig_adjusted_p_value,
        ];
    }
    close $fh;

    return \%data_for;
}

# Get list of regions identical between two lists
sub get_identical_regions {
    my ( $regions1, $regions2 ) = @_;

    my %count;
    foreach my $region ( @{$regions1}, @{$regions2} ) {
        $count{$region}++;
    }

    my @regions;
    foreach my $region ( sort keys %count ) {
        if ( $count{$region} == 2 ) {
            push @regions, $region;
        }
    }

    return \@regions;
}

# Statistically summarise all data for a specific field
sub summarise_field_for_all_data {
    my ( $data_for, $field, $missing_value ) = @_;

    my $stat          = Statistics::Descriptive::Full->new();
    my $missing_count = 0;

    foreach my $data ( values %{$data_for} ) {
        my $value = $data->[$field];
        if ( !defined $value
            || ( defined $missing_value && $value eq $missing_value ) )
        {
            $missing_count++;
        }
        else {
            $stat->add_data($value);
        }
    }

    return {
        missing_count => $missing_count,
        count         => $stat->count,
        mean          => $stat->mean,
        sd            => $stat->standard_deviation,
        median        => $stat->median,
        mode          => $stat->mode,
    };
}

# Compare specific numeric field for subset of data
sub compare_numeric_field_for_data_subset {
    my ( $data1, $data2, $regions, $field, $missing_value ) = @_;

    my %count_of = (
        missing_value_both => 0,
        missing_value1     => 0,
        missing_value2     => 0,
        identical_value    => 0,
        lower_value        => 0,
        higher_value       => 0,
    );
    foreach my $region ( @{$regions} ) {
        my $value1       = $data1->{$region}->[$field];
        my $value2       = $data2->{$region}->[$field];
        my $value1_undef = !defined $value1
          || ( defined $missing_value && $value1 eq $missing_value );
        my $value2_undef = !defined $value2
          || ( defined $missing_value && $value2 eq $missing_value );

        my $key =
            $value1_undef && $value2_undef ? 'missing_value_both'
          : $value1_undef      ? 'missing_value1'
          : $value2_undef      ? 'missing_value2'
          : $value1 == $value2 ? 'identical_value'
          : $value1 < $value2  ? 'lower_value'
          : $value1 > $value2  ? 'higher_value'
          :                      undef;
        if ( !defined $key ) {
            carp sprintf "Unaccounted for combination of %s and %s\n", $value1,
              $value2;
        }
        $count_of{$key}++;
    }

    return \%count_of;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'file1=s'     => \$file1,
        'file2=s'     => \$file2,
        'sig_level=i' => \$sig_level,
        'help'        => \$help,
        'man'         => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    # Check options
    if ( !$file1 ) {
        pod2usage("--file1 must be specified\n");
    }
    if ( !$file2 ) {
        pod2usage("--file2 must be specified\n");
    }

    return;
}

=head1 USAGE

    compare_de_output.pl
        [--file1 file]
        [--file2 file]
        [--sig_level float]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--file1 FILE>

First differential expression pipeline output file (e.g. 1/all.tsv).

=item B<--file2 FILE>

Second differential expression pipeline output file (e.g. 2/all.tsv).

=item B<--sig_level FLOAT>

Adjusted p-value significance level (defaults to 0.05).

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
