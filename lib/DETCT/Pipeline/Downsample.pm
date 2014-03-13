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
  stats_by_tag
  downsample_by_tag
);
use DETCT::Misc::Picard qw(
  mark_duplicates
  merge
);
use DETCT::Misc::SAMtools qw(
  make_index
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=method all_parameters_for_stats_by_tag

  Usage       : all_parameters_for_stats_by_tag();
  Purpose     : Get all parameters for stats_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_stats_by_tag {
    my ($self) = @_;

    my @all_parameters;

    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        push @all_parameters, [ $bam_file, @tags ];
    }

    return @all_parameters;
}

=method run_stats_by_tag

  Usage       : run_stats_by_tag();
  Purpose     : Run function for stats_by_tag stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_stats_by_tag {
    my ( $self, $job ) = @_;

    my ( $bam_file, @tags ) = @{ $job->parameters };

    # Get stats for BAM file
    my $stats = stats_by_tag(
        {
            bam_file => $bam_file,
            tags     => \@tags,
        }
    );

    my $output_file = $job->base_filename . '.out';

    DumpFile( $output_file, $stats );

    return;
}

=method all_parameters_for_downsample_by_tag

  Usage       : all_parameters_for_downsample_by_tag();
  Purpose     : Get all parameters for downsample_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_downsample_by_tag {
    my ($self) = @_;

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $stats_output_file =
          $self->get_and_check_output_file( 'stats_by_tag', $component );
        my $stats = LoadFile($stats_output_file);
        my @tags  = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            my $target_read_count =
              $stats->{$tag}->{ $self->analysis->read_count_type };
            push @all_parameters, [ $bam_file, $tag, $target_read_count ];
        }
    }

    # Calculate maximum possible target read count
    my $max = min( map { $_->[2] } @all_parameters );

    # If no target read count specified or larger than maximum then use maximum
    if (  !$self->analysis->target_read_count
        || $self->analysis->target_read_count > $max )
    {
        $self->analysis->set_target_read_count($max);
    }

    # Round down to specified whole number if required
    if ( $self->analysis->round_down_to ) {
        my $count = $self->analysis->target_read_count;
        $count = $count - ( $count % $self->analysis->round_down_to );
        $self->analysis->set_target_read_count($count);
    }

    return @all_parameters;
}

=method run_downsample_by_tag

  Usage       : run_downsample_by_tag();
  Purpose     : Run function for downsample_by_tag stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_downsample_by_tag {
    my ( $self, $job ) = @_;

    my ( $source_bam_file, $tag, $source_read_count ) = @{ $job->parameters };

    # Downsample BAM file
    my $read_count = downsample_by_tag(
        {
            source_bam_file   => $source_bam_file,
            source_read_count => $source_read_count,
            tag               => $tag,
            target_bam_file   => $job->base_filename . '.bam',
            target_read_count => $self->analysis->target_read_count,
            read_count_type   => $self->analysis->read_count_type,
        }
    );

    # Only mark as finished if reached target read count
    # Will happen 50% of the time
    if ( $read_count == $self->analysis->target_read_count ) {
        my $output_file = $job->base_filename . '.out';
        DumpFile( $output_file, 1 );
    }

    return;
}

=method all_parameters_for_mark_duplicates

  Usage       : all_parameters_for_mark_duplicates();
  Purpose     : Get all parameters for mark_duplicates stage
  Returns     : Array of arrayrefs
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_mark_duplicates {
    my ($self) = @_;

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            $component++;
            my $input_bam_file =
              File::Spec->catfile( $self->analysis_dir, 'downsample_by_tag',
                $component . '.bam' );
            push @all_parameters, [$input_bam_file];
        }
    }

    return @all_parameters;
}

=method run_mark_duplicates

  Usage       : run_mark_duplicates();
  Purpose     : Run function for mark_duplicates stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_mark_duplicates {
    my ( $self, $job ) = @_;

    my ($input_bam_file) = @{ $job->parameters };

    # Mark duplicates
    mark_duplicates(
        {
            dir                 => $job->base_filename,
            input_bam_file      => $input_bam_file,
            output_bam_file     => $job->base_filename . '.bam',
            metrics_file        => $job->base_filename . '.metrics',
            java_binary         => $self->analysis->java_binary,
            mark_duplicates_jar => $self->analysis->mark_duplicates_jar,
            memory              => $job->memory,
            consider_tags       => 1,
        }
    );

    my $output_file = $job->base_filename . '.out';

    DumpFile( $output_file, 1 );

    return;
}

=method all_parameters_for_merge

  Usage       : all_parameters_for_merge();
  Purpose     : Get all parameters for merge stage
  Returns     : Array of arrayrefs
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_merge {
    my ($self) = @_;

    my @all_parameters;

    my @input_bam_files;
    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            $component++;
            my $input_bam_file =
              File::Spec->catfile( $self->analysis_dir, 'mark_duplicates',
                $component . '.bam' );
            push @input_bam_files, $input_bam_file;
        }
    }
    push @all_parameters, \@input_bam_files;

    return @all_parameters;
}

=method run_merge

  Usage       : run_merge();
  Purpose     : Run function for merge stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_merge {
    my ( $self, $job ) = @_;

    my $input_bam_files = $job->parameters;

    my $output_bam_file = File::Spec->catfile( $self->analysis_dir,
        $self->analysis->name . '.bam' );

    # Merge BAM files
    merge(
        {
            dir                 => $job->base_filename,
            input_bam_files     => $input_bam_files,
            output_bam_file     => $output_bam_file,
            java_binary         => $self->analysis->java_binary,
            merge_sam_files_jar => $self->analysis->merge_sam_files_jar,
            memory              => $job->memory,
        }
    );

    my $output_file = $job->base_filename . '.out';

    DumpFile( $output_file, 1 );

    return;
}

=method all_parameters_for_make_index

  Usage       : all_parameters_for_make_index();
  Purpose     : Get all parameters for make_index stage
  Returns     : Array of arrayrefs
  Parameters  : None
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_make_index {
    my ($self) = @_;

    return [];
}

=method run_make_index

  Usage       : run_make_index();
  Purpose     : Run function for make_index stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_make_index {
    my ( $self, $job ) = @_;

    my $bam_file = File::Spec->catfile( $self->analysis_dir,
        $self->analysis->name . '.bam' );

    # Index BAM file
    make_index(
        {
            dir             => $job->base_filename,
            bam_file        => $bam_file,
            samtools_binary => $self->analysis->samtools_binary,
        }
    );

    my $output_file = $job->base_filename . '.out';

    DumpFile( $output_file, 1 );

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
