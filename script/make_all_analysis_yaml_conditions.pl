#!/usr/bin/env perl

# PODNAME: make_all_analysis_yaml_conditions.pl
# ABSTRACT: Create analysis YAML files with appropriate pairs of conditions

## Author         : is1
## Maintainer     : is1
## Created        : 2018-09-27
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
use YAML::Tiny;
use Path::Tiny;
use File::Spec;
use File::Path qw( make_path );
use DETCT::Misc qw( write_or_die );
use Algorithm::Combinatorics qw( partitions subsets );

=head1 DESCRIPTION



=head1 EXAMPLES



=cut

# Default options
my $base_yaml = 'analysis.yaml';
my $output_dir;
my $table_file;
my @comparisons;
my @ignore;
my ( $help, $man );

# Get and check command line options
get_and_check_options();

# Make output directory
if ( !-d $output_dir ) {
    make_path($output_dir);
}

my ( $yaml, $yaml_tiny, $is_condition, @all_conditions ) =
  get_conditions( $base_yaml, @ignore );

# Assume some comparisons
if ( !@comparisons ) {
    @comparisons = default_comparisons(@all_conditions);
}

# Make all possible comparisons
if ( !@comparisons ) {
    @comparisons = all_comparisons(@all_conditions);
}

# Write YAML files
foreach my $comparison (@comparisons) {
    my ( $all_exp, $all_con ) = split /:/xms, $comparison;
    confess "Experimental condition missing from $comparison for $base_yaml"
      if !$all_exp;
    confess "Control condition missing from $comparison for $base_yaml"
      if !$all_con;
    my ( $exp, $exp_name ) = split /=/xms, $all_exp;
    my ( $con, $con_name ) = split /=/xms, $all_con;
    if ( !$exp_name ) {
        $exp_name = $exp;
        $exp_name =~ s/,/_/xmsg;
    }
    if ( !$con_name ) {
        $con_name = $con;
        $con_name =~ s/,/_/xmsg;
    }
    my @exp = split /,/xms, $exp;
    my @con = split /,/xms, $con;

    my %rename;
    foreach my $condition (@exp) {
        confess "Unknown condition ($condition) in $comparison for $output_dir"
          if !$is_condition->{$condition};
        $rename{$condition} = $exp_name;
    }
    foreach my $condition (@con) {
        confess "Unknown condition ($condition) in $comparison for $output_dir"
          if !$is_condition->{$condition};
        $rename{$condition} = $con_name;
    }
    my @output;
    foreach my $line ( split /\n/xms, $yaml ) {
        if ( $line =~ m/\A \s+ condition: \s+ (\S+) \z/xms ) {
            my $condition = $1;
            if ( exists $rename{$condition} ) {
                $line =~ s/$condition \z/$rename{$condition}/xms;
            }
        }
        push @output, $line;
    }

    my $output_file = File::Spec->catfile( $output_dir,
        $exp_name . '_vs_' . $con_name . '.yaml' );
    next if -e $output_file;
    if ($table_file) {
        push @output, sprintf 'table_file: %s', $table_file;
    }
    push @output, sprintf 'experimental_condition: %s', $exp_name;
    push @output, sprintf 'control_condition: %s',      $con_name;
    open my $fh, '>', $output_file;
    write_or_die( $fh, ( join "\n", @output ), "\n" );
    close $fh;
}

sub get_conditions {
    my ( $base_yaml, @ignore ) = @_;    ## no critic (ProhibitReusedNames)

    # Read analysis config
    confess "YAML file ($base_yaml) does not exist or cannot be read"
      if !-r $base_yaml;
    ## no critic (ProhibitReusedNames)
    my $yaml      = path($base_yaml)->slurp;
    my $yaml_tiny = YAML::Tiny->read_string($yaml);
    ## use critic
    if ( !$yaml_tiny ) {
        confess sprintf 'YAML file (%s) is invalid: %s', $base_yaml,
          YAML::Tiny->errstr;
    }

    # Get all conditions
    my %is_condition;
    foreach my $sample ( @{ $yaml_tiny->[0]->{'samples'} } ) {
        $is_condition{ $sample->{'condition'} } = 1;
    }
    ## no critic (ProhibitReusedNames)
    my @all_conditions = sort keys %is_condition;
    ## use critic

    # Removed ignored conditions
    foreach my $condition (@ignore) {
        @all_conditions = grep { $_ ne $condition } @all_conditions;
    }

    # Sanity check
    if ( scalar @all_conditions == 1 ) {
        confess "Only one condition (@all_conditions) for $base_yaml";
    }
    if ( scalar @all_conditions == scalar @{ $yaml_tiny->[0]->{'samples'} } ) {
        confess "One condition per sample for $base_yaml";
    }

    return $yaml, $yaml_tiny, \%is_condition, @all_conditions;
}

sub default_comparisons {
    my @conds = @_;

    my @comps;

    my ($wt)  = grep { $_ eq 'wt' } @conds;
    my ($het) = grep { $_ eq 'het' } @conds;
    my ($hom) = grep { $_ eq 'hom' } @conds;
    my ($sib) = grep { $_ eq 'sib' } @conds;
    my ($mut) = grep { $_ eq 'mut' } @conds;

    if ( $wt && $het && $hom ) {
        push @comps, "$het:$wt";
        push @comps, "$hom:$wt";
        push @comps, "$hom:$het";
        push @comps, "$hom:$het,$wt";
        push @comps, "$hom,$het:$wt";
    }
    if ( !@comps ) {
        @comps =
            $wt  && $het ? ("$het:$wt")
          : $wt  && $hom ? ("$hom:$wt")
          : $het && $hom ? ("$hom:$het")
          : $sib && $mut ? ("$mut:$sib")
          :                ();
    }

    return @comps;
}

sub all_comparisons {
    my @conds = @_;

    my @comps;

    if ( scalar @conds < 5 ) {    ## no critic (ProhibitMagicNumbers)
                                  # All possible comparisons
        my $subsets = subsets( \@conds );
        while ( my $subset = $subsets->next ) {
            next if scalar @{$subset} < 2;
            my $partitions = partitions( $subset, 2 );
            while ( my $partition = $partitions->next ) {
                push @comps,
                  ( join q{,}, @{ $partition->[0] } ) . q{:}
                  . ( join q{,}, @{ $partition->[1] } );
                push @comps,
                  ( join q{,}, @{ $partition->[1] } ) . q{:}
                  . ( join q{,}, @{ $partition->[0] } );
            }
        }
    }
    else {
        # Just pairwise comparisons
        my $cond1 = shift @conds;
        foreach my $cond2 (@conds) {
            push @comps, $cond1 . q{:} . $cond2;
        }
    }

    return @comps;
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'base_yaml=s'       => \$base_yaml,
        'output_dir=s'      => \$output_dir,
        'table_file=s'      => \$table_file,
        'comparisons=s@{,}' => \@comparisons,
        'ignore=s@{,}'      => \@ignore,
        'help'              => \$help,
        'man'               => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    if ( !$output_dir ) {
        pod2usage("--output_dir must be specified\n");
    }

    return;
}

=head1 USAGE

    make_all_analysis_yaml_conditions.pl
        [--base_yaml file]
        [--output_dir dir]
        [--table_file file]
        [--comparisons comparison...]
        [--ignore condition...]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--base_yaml FILE>

YAML file to use as base for new YAML files (and to extract conditions from).

=item B<--output_dir DIR>

Directory in which to create YAML files.

=item B<--table_file FILE>

Differential expression pipeline CSV or TSV output file.

=item B<--comparisons COMPARISONS>

Condition comparisons. Each comparison is a pair of experimental and controls
conditions (in that order) separated by a colon (e.g. hom:wt). If multiple
conditions are to be combined then separate them with a comma (e.g. hom:wt,het).
To rename a condition, append an equals sign (e.g. hom=mut:het,wt=sib).

=item B<--ignore CONDITIONS>

Conditions to be ignored.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
