## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Pipeline::Downsample;
## use critic

# ABSTRACT: Object representing a downsampling pipeline

## Author         : is1
## Maintainer     : is1
## Created        : 2014-02-17
## Last commit by : $Author$
## Last modified  : $Date$
## Revision       : $Revision$
## Repository URL : $HeadURL$

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use parent qw(DETCT::Pipeline);

use Class::InsideOut qw( private register id );
use List::Util qw( min );
use YAML qw( DumpFile LoadFile );
use File::Spec;
use DETCT::Misc::BAM qw(
  stats
  downsample
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=method all_parameters_for_stats

  Usage       : all_parameters_for_stats();
  Purpose     : Get all parameters for stats stage
  Returns     : Array of arrayrefs
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_stats {
    my ($self) = @_;

    my @all_parameters;

    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        push @all_parameters, [ $bam_file, @tags ];
    }

    return @all_parameters;
}

=method run_stats

  Usage       : run_stats();
  Purpose     : Run function for stats stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_stats {
    my ( $self, $job ) = @_;

    my ( $bam_file, @tags ) = @{ $job->parameters };

    # Get stats for BAM file
    my $stats = stats(
        {
            bam_file => $bam_file,
            tags     => \@tags,
        }
    );

    my $output_file = $job->base_filename . '.out';

    DumpFile( $output_file, $stats );

    return;
}

=method all_parameters_for_downsample

  Usage       : all_parameters_for_downsample();
  Purpose     : Get all parameters for downsample stage
  Returns     : Array of arrayrefs
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_downsample {
    my ($self) = @_;

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $stats_output_file =
          $self->get_and_check_output_file( 'stats', $component );
        my $stats = LoadFile($stats_output_file);
        my @tags  = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            my $target_read_count =
              $stats->{$tag}->{ $self->analysis->read_count_type };
            push @all_parameters, [ $bam_file, $tag, $target_read_count ];
        }
    }

    # Calculate read count to downsample to
    my $min = min( map { $_->[2] } @all_parameters );
    my $target_read_count = $min - ( $min % $self->analysis->round_down_to );
    ## no critic (ProhibitMagicNumbers)
    map { $_->[3] = $target_read_count } @all_parameters;
    ## use critic

    return @all_parameters;
}

=method run_downsample

  Usage       : run_downsample();
  Purpose     : Run function for downsample stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_downsample {
    my ( $self, $job ) = @_;

    my ( $source_bam_file, $tag, $source_read_count, $target_read_count ) =
      @{ $job->parameters };

    # Downsample BAM file
    my $read_count = downsample(
        {
            source_bam_file   => $source_bam_file,
            source_read_count => $source_read_count,
            tag               => $tag,
            target_bam_file   => $job->base_filename . '.bam',
            target_read_count => $target_read_count,
            read_count_type   => $self->analysis->read_count_type,
        }
    );

    # Only mark as finished if reached target read count
    # Will happen 50% of the time
    if ( $read_count == $target_read_count ) {
        my $output_file = $job->base_filename . '.out';
        DumpFile( $output_file, 1 );
    }

    return;
}

=method input_overview

  Usage       : $pipeline->say_if_verbose($pipeline->input_overview);
  Purpose     : Return textual overview of pipeline's input
  Returns     : Array of Strings
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub input_overview {
    my ($self) = @_;

    my @output;

    push @output, 'Command line:', $self->cmd_line;
    if ( defined $DETCT::VERSION ) {
        push @output, 'DETCT version: ' . $DETCT::VERSION;
    }
    push @output, 'Working directory: ' . $self->analysis_dir;

    push @output, 'BAM files: ' . join q{ },
      $self->analysis->list_all_bam_files();

    push @output, sprintf 'Number of samples: %d',
      scalar @{ $self->analysis->get_all_samples };

    return @output;
}

1;
