## no critic (RequireUseStrict, RequireUseWarnings, RequireTidyCode)
package DETCT::Pipeline::Downsample;
## VERSION
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
use List::Util qw( min sum );
use File::Spec;
use File::Slurp;
use DETCT::Misc::BAM qw(
  stats_by_tag
  stats_all_reads
  downsample_by_tag
  downsample_all_reads
);
use DETCT::Misc::Picard qw(
  extract_mark_duplicates_metrics
  merge
);
use DETCT::Misc::SAMtools qw(
  make_index
  flagstats
);
use DETCT::Misc::Output qw(
  dump_duplication_metrics
);

=head1 SYNOPSIS

    # Brief code examples

=cut

=method all_parameters_for_stats_by_tag

  Usage       : all_parameters_for_stats_by_tag();
  Purpose     : Get all parameters for stats_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_stats_by_tag {
    my ( $self, $stage ) = @_;

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
            bam_file       => $bam_file,
            tags           => \@tags,
            skip_sequences => $self->analysis->get_all_skip_sequences(),
        }
    );

    my $output_file = $job->base_filename . '.out';

    $self->dump_serialised( $output_file, $stats );

    return;
}

=method all_parameters_for_stats_all_reads

  Usage       : all_parameters_for_stats_all_reads();
  Purpose     : Get all parameters for stats_all_reads stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_stats_all_reads {
    my ( $self, $stage ) = @_;

    my @all_parameters;

    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        push @all_parameters, [$bam_file];
    }

    return @all_parameters;
}

=method run_stats_all_reads

  Usage       : run_stats_all_reads();
  Purpose     : Run function for stats_all_reads stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_stats_all_reads {
    my ( $self, $job ) = @_;

    my ($bam_file) = @{ $job->parameters };

    # Get stats for BAM file
    my $stats = stats_all_reads(
        {
            bam_file       => $bam_file,
            skip_sequences => $self->analysis->get_all_skip_sequences(),
        }
    );

    my $output_file = $job->base_filename . '.out';

    $self->dump_serialised( $output_file, $stats );

    return;
}

=method all_parameters_for_downsample_by_tag

  Usage       : all_parameters_for_downsample_by_tag();
  Purpose     : Get all parameters for downsample_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_downsample_by_tag {
    my ( $self, $stage ) = @_;

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $stats_output_file =
          $self->get_and_check_output_file( 'stats_by_tag', $component );
        my $stats = $self->load_serialised($stats_output_file);
        my @tags  = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            my $source_read_count =
              $stats->{$tag}->{ $self->analysis->read_count_type };
            push @all_parameters, [ $bam_file, $tag, $source_read_count ];
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
            skip_sequences    => $self->analysis->get_all_skip_sequences(),
        }
    );

    # Only mark as finished if reached target read count
    # Will happen 50% of the time
    if ( $read_count == $self->analysis->target_read_count ) {
        my $output_file = $job->base_filename . '.out';
        my %output      = (
            source_bam_file   => $source_bam_file,
            tag               => $tag,
            read_count_type   => $self->analysis->read_count_type,
            source_read_count => $source_read_count,
            target_read_count => $self->analysis->target_read_count,
        );
        $self->dump_serialised( $output_file, \%output );
    }
    else {
        exit 1;
    }

    return;
}

=method all_parameters_for_downsample_all_reads

  Usage       : all_parameters_for_downsample_all_reads();
  Purpose     : Get all parameters for downsample_all_reads stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_downsample_all_reads {
    my ( $self, $stage ) = @_;

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $stats_output_file =
          $self->get_and_check_output_file( 'stats_all_reads', $component );
        my $stats             = $self->load_serialised($stats_output_file);
        my $source_read_count = $stats->{ $self->analysis->read_count_type };
        push @all_parameters, [ $bam_file, $source_read_count ];
    }

    # Calculate total read count
    my $total_read_count = sum( map { $_->[1] } @all_parameters );

    # If no target read count specified or larger than total then use total
    if (  !$self->analysis->target_read_count
        || $self->analysis->target_read_count > $total_read_count )
    {
        $self->analysis->set_target_read_count($total_read_count);
    }

    # Round down to specified whole number if required
    if ( $self->analysis->round_down_to ) {
        my $count = $self->analysis->target_read_count;
        $count = $count - ( $count % $self->analysis->round_down_to );
        $self->analysis->set_target_read_count($count);
    }

    # Calculate proportion of target read count for each BAM file
    my $count_remaining = $self->analysis->target_read_count;
    foreach my $parameters (@all_parameters) {
        my ( $bam_file, $source_read_count ) = @{$parameters};
        my $target_read_count =
          int( $source_read_count /
              $total_read_count *
              $self->analysis->target_read_count );
        push @{$parameters}, $target_read_count;
        $count_remaining -= $target_read_count;
    }

    # Add any remaining counts to final parameters
    $all_parameters[-1]->[2] += $count_remaining;

    return @all_parameters;
}

=method run_downsample_all_reads

  Usage       : run_downsample_all_reads();
  Purpose     : Run function for downsample_all_reads stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_downsample_all_reads {
    my ( $self, $job ) = @_;

    my ( $source_bam_file, $source_read_count, $target_read_count ) =
      @{ $job->parameters };

    # Downsample BAM file
    my $read_count = downsample_all_reads(
        {
            source_bam_file   => $source_bam_file,
            source_read_count => $source_read_count,
            target_bam_file   => $job->base_filename . '.bam',
            target_read_count => $target_read_count,
            read_count_type   => $self->analysis->read_count_type,
            skip_sequences    => $self->analysis->get_all_skip_sequences(),
        }
    );

    # Only mark as finished if reached target read count
    # Will happen 50% of the time
    if ( $read_count == $target_read_count ) {
        my $output_file = $job->base_filename . '.out';
        my %output      = (
            source_bam_file       => $source_bam_file,
            read_count_type       => $self->analysis->read_count_type,
            source_read_count     => $source_read_count,
            sub_target_read_count => $target_read_count,
            target_read_count     => $self->analysis->target_read_count,
        );
        $self->dump_serialised( $output_file, \%output );
    }
    else {
        exit 1;
    }

    return;
}

=method all_parameters_for_sort_by_queryname_by_tag

  Usage       : all_parameters_for_sort_by_queryname_by_tag();
  Purpose     : Get all parameters for sort_by_queryname_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_sort_by_queryname_by_tag {
    my ( $self, $stage ) = @_;

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

=method run_sort_by_queryname_by_tag

  Usage       : run_sort_by_queryname_by_tag();
  Purpose     : Run function for sort_by_queryname_by_tag stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sort_by_queryname_by_tag {
    my ( $self, $job ) = @_;

    return $self->run_sort_by_queryname($job);
}

=method all_parameters_for_sort_by_queryname_all_reads

  Usage       : all_parameters_for_sort_by_queryname_all_reads();
  Purpose     : Get all parameters for sort_by_queryname_all_reads stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_sort_by_queryname_all_reads {
    my ( $self, $stage ) = @_;

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $input_bam_file =
          File::Spec->catfile( $self->analysis_dir, 'downsample_all_reads',
            $component . '.bam' );
        push @all_parameters, [$input_bam_file];
    }

    return @all_parameters;
}

=method run_sort_by_queryname_all_reads

  Usage       : run_sort_by_queryname_all_reads();
  Purpose     : Run function for sort_by_queryname_all_reads stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sort_by_queryname_all_reads {
    my ( $self, $job ) = @_;

    return $self->run_sort_by_queryname($job);
}

=method run_sort_by_queryname

  Usage       : run_sort_by_queryname();
  Purpose     : Run function for sort_by_queryname stages
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sort_by_queryname {
    my ( $self, $job ) = @_;

    my ($input_bam_file) = @{ $job->parameters };

    if ( $self->analysis->mark_duplicates_method eq 'native' ) {

        # Sort BAM file by read name
        DETCT::Misc::Picard::sort_bam(
            {
                dir             => $job->base_filename,
                input_bam_file  => $input_bam_file,
                output_bam_file => $job->base_filename . '.bam',
                java_binary     => $self->analysis->java_binary,
                sort_bam_jar    => $self->analysis->sort_bam_jar,
                sort_order      => 'queryname',
            }
        );
    }

    my $output_file = $job->base_filename . '.out';

    $self->dump_serialised( $output_file, \1 );

    return;
}

=method all_parameters_for_mark_duplicates_by_tag

  Usage       : all_parameters_for_mark_duplicates_by_tag();
  Purpose     : Get all parameters for mark_duplicates_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_mark_duplicates_by_tag {
    my ( $self, $stage ) = @_;

    # Previous stage will be skipped if using native mark duplicates method
    my $prev_stage_name =
      $self->analysis->mark_duplicates_method eq 'native'
      ? 'sort_by_queryname_by_tag'
      : 'downsample_by_tag';

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            $component++;
            my $input_bam_file =
              File::Spec->catfile( $self->analysis_dir, $prev_stage_name,
                $component . '.bam' );
            push @all_parameters, [ $input_bam_file, $tag ];
        }
    }

    return @all_parameters;
}

=method run_mark_duplicates_by_tag

  Usage       : run_mark_duplicates_by_tag();
  Purpose     : Run function for mark_duplicates_by_tag stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_mark_duplicates_by_tag {
    my ( $self, $job ) = @_;

    return $self->run_mark_duplicates($job);
}

=method all_parameters_for_mark_duplicates_all_reads

  Usage       : all_parameters_for_mark_duplicates_all_reads();
  Purpose     : Get all parameters for mark_duplicates_all_reads stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_mark_duplicates_all_reads {
    my ( $self, $stage ) = @_;

    # Previous stage will be skipped if using native mark duplicates method
    my $prev_stage_name =
      $self->analysis->mark_duplicates_method eq 'native'
      ? 'sort_by_queryname_all_reads'
      : 'downsample_all_reads';

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        my $input_bam_file =
          File::Spec->catfile( $self->analysis_dir, $prev_stage_name,
            $component . '.bam' );
        push @all_parameters, [ $input_bam_file, @tags ];
    }

    return @all_parameters;
}

=method run_mark_duplicates_all_reads

  Usage       : run_mark_duplicates_all_reads();
  Purpose     : Run function for mark_duplicates_all_reads stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_mark_duplicates_all_reads {
    my ( $self, $job ) = @_;

    return $self->run_mark_duplicates($job);
}

=method run_mark_duplicates

  Usage       : run_mark_duplicates();
  Purpose     : Run function for mark_duplicates stages
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_mark_duplicates {
    my ( $self, $job ) = @_;

    my ( $input_bam_file, @tags ) = @{ $job->parameters };

    # No need for tag-specific metrics if just one tag
    if ( scalar @tags == 1 ) {
        @tags = ();
    }

    my $output = \1;

    if ( $self->analysis->mark_duplicates_method eq 'native' ) {

        # Mark Duplicates with native Perl-based method
        $output = DETCT::Misc::BAM::mark_duplicates(
            {
                input_bam_file  => $input_bam_file,
                output_bam_file => $job->base_filename . '.bam',
                consider_tags   => 1,
                tags            => \@tags,
            }
        );
    }
    else {
        # Mark duplicates with Picard MarkDuplicates
        DETCT::Misc::Picard::mark_duplicates(
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
    }

    my $output_file = $job->base_filename . '.out';

    $self->dump_serialised( $output_file, $output );

    return;
}

=method all_parameters_for_mark_duplicates_metrics_by_tag

  Usage       : all_parameters_for_mark_duplicates_metrics_by_tag();
  Purpose     : Get all parameters for mark_duplicates_metrics_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_mark_duplicates_metrics_by_tag {
    my ( $self, $stage ) = @_;

    my @all_parameters;

    my $metrics_output_file =
      File::Spec->catfile( $self->analysis_dir,
        'mark_duplicates_metrics_by_tag',
        'markduplicates.tsv' );

    my @parameters = ($metrics_output_file);
    my $component  = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            $component++;
            my $metrics_file;
            if ( $self->analysis->mark_duplicates_method eq 'native' ) {
                $metrics_file =
                  $self->get_and_check_output_file( 'mark_duplicates_by_tag',
                    $component );
            }
            else {
                $metrics_file = File::Spec->catfile( $self->analysis_dir,
                    'mark_duplicates_by_tag', $component . '.metrics' );
            }
            my $sample_name =
              $self->analysis->get_sample_name_by_bam_file_and_tag( $bam_file,
                $tag );
            push @parameters, [ $sample_name, $metrics_file ];
        }
    }
    push @all_parameters, \@parameters;

    return @all_parameters;
}

=method run_mark_duplicates_metrics_by_tag

  Usage       : run_mark_duplicates_metrics_by_tag();
  Purpose     : Run function for mark_duplicates_metrics_by_tag stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_mark_duplicates_metrics_by_tag {
    my ( $self, $job ) = @_;

    return $self->run_mark_duplicates_metrics($job);
}

=method all_parameters_for_mark_duplicates_metrics_all_reads

  Usage       : all_parameters_for_mark_duplicates_metrics_all_reads();
  Purpose     : Get all parameters for mark_duplicates_metrics_all_reads stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_mark_duplicates_metrics_all_reads {
    my ( $self, $stage ) = @_;

    my @all_parameters;

    my $metrics_output_file =
      File::Spec->catfile( $self->analysis_dir,
        'mark_duplicates_metrics_all_reads',
        'markduplicates.tsv' );

    my @parameters = ($metrics_output_file);
    my $component  = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $sample_names = join q{ },
          $self->analysis->get_sample_names_by_bam_file($bam_file);
        my $metrics_file;
        if ( $self->analysis->mark_duplicates_method eq 'native' ) {
            $metrics_file =
              $self->get_and_check_output_file( 'mark_duplicates_all_reads',
                $component );
        }
        else {
            $metrics_file = File::Spec->catfile( $self->analysis_dir,
                'mark_duplicates_all_reads', $component . '.metrics' );
        }
        push @parameters, [ $sample_names, $metrics_file ];
    }
    push @all_parameters, \@parameters;

    return @all_parameters;
}

=method run_mark_duplicates_metrics_all_reads

  Usage       : run_mark_duplicates_metrics_all_reads();
  Purpose     : Run function for mark_duplicates_metrics_all_reads stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_mark_duplicates_metrics_all_reads {
    my ( $self, $job ) = @_;

    return $self->run_mark_duplicates_metrics($job);
}

=method run_mark_duplicates_metrics

  Usage       : run_mark_duplicates_metrics();
  Purpose     : Run function for mark_duplicates_metrics stages
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_mark_duplicates_metrics {
    my ( $self, $job ) = @_;

    my (@parameters) = @{ $job->parameters };

    my $metrics_output_file = shift @parameters;

    # Get metrics
    my @metrics;
    foreach my $parameter (@parameters) {
        my ( $sample_name, $metrics_file ) = @{$parameter};

        if ( $self->analysis->mark_duplicates_method eq 'native' ) {

            # Get metrics output from native method
            my $output = $self->load_serialised($metrics_file);
            my @tags = sort grep { m/\A_/xms } keys %{$output};
            if ( !@tags ) {

                # Just one set of metrics
                $output->{_all}->{sample_name} = $sample_name;
                push @metrics, $output->{_all};
            }
            else {
                # Metrics for each tag
                $output->{_all}->{sample_name} = $sample_name . ' All';
                push @metrics, $output->{_all};
                foreach my $tag (@tags) {
                    $output->{$tag}->{sample_name} = $sample_name . q{ } . $tag;
                    push @metrics, $output->{$tag};
                }
                $output->{_other}->{sample_name} = $sample_name . ' Other';
                push @metrics, $output->{_other};
            }
        }
        else {
            # Get metrics from Picard MarkDuplicates output file
            my $output = extract_mark_duplicates_metrics(
                { metrics_file => $metrics_file, } );
            $output->{sample_name} = $sample_name;
            push @metrics, $output;
        }
    }

    # Output metrics
    dump_duplication_metrics(
        {
            output_file => $metrics_output_file,
            metrics     => \@metrics,
        }
    );

    my $output_file = $job->base_filename . '.out';

    $self->dump_serialised( $output_file, \1 );

    return;
}

=method all_parameters_for_sort_by_coordinate_by_tag

  Usage       : all_parameters_for_sort_by_coordinate_by_tag();
  Purpose     : Get all parameters for sort_by_coordinate_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_sort_by_coordinate_by_tag {
    my ( $self, $stage ) = @_;

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            $component++;
            my $input_bam_file =
              File::Spec->catfile( $self->analysis_dir,
                'mark_duplicates_by_tag', $component . '.bam' );
            push @all_parameters, [$input_bam_file];
        }
    }

    return @all_parameters;
}

=method run_sort_by_coordinate_by_tag

  Usage       : run_sort_by_coordinate_by_tag();
  Purpose     : Run function for sort_by_coordinate_by_tag stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sort_by_coordinate_by_tag {
    my ( $self, $job ) = @_;

    return $self->run_sort_by_coordinate($job);
}

=method all_parameters_for_sort_by_coordinate_all_reads

  Usage       : all_parameters_for_sort_by_coordinate_all_reads();
  Purpose     : Get all parameters for sort_by_coordinate_all_reads stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_sort_by_coordinate_all_reads {
    my ( $self, $stage ) = @_;

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $input_bam_file =
          File::Spec->catfile( $self->analysis_dir, 'mark_duplicates_all_reads',
            $component . '.bam' );
        push @all_parameters, [$input_bam_file];
    }

    return @all_parameters;
}

=method run_sort_by_coordinate_all_reads

  Usage       : run_sort_by_coordinate_all_reads();
  Purpose     : Run function for sort_by_coordinate_all_reads stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sort_by_coordinate_all_reads {
    my ( $self, $job ) = @_;

    return $self->run_sort_by_coordinate($job);
}

=method run_sort_by_coordinate

  Usage       : run_sort_by_coordinate();
  Purpose     : Run function for sort_by_coordinate stages
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sort_by_coordinate {
    my ( $self, $job ) = @_;

    my ($input_bam_file) = @{ $job->parameters };

    if ( $self->analysis->mark_duplicates_method eq 'native' ) {

        # Sort BAM file by coordinate
        DETCT::Misc::SAMtools::sort_bam(
            {
                dir             => $job->base_filename,
                input_bam_file  => $input_bam_file,
                output_bam_file => $job->base_filename . '.bam',
                samtools_binary => $self->analysis->samtools_binary,
                sort_order      => 'coordinate',
            }
        );
    }

    my $output_file = $job->base_filename . '.out';

    $self->dump_serialised( $output_file, \1 );

    return;
}

=method all_parameters_for_sample_flagstats_by_tag

  Usage       : all_parameters_for_sample_flagstats_by_tag();
  Purpose     : Get all parameters for sample_flagstats_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_sample_flagstats_by_tag {
    my ( $self, $stage ) = @_;

    # Previous stage will be skipped if using native mark duplicates method
    my $prev_stage_name =
      $self->analysis->mark_duplicates_method eq 'native'
      ? 'sort_by_coordinate_by_tag'
      : 'mark_duplicates_by_tag';

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            $component++;
            my $input_bam_file =
              File::Spec->catfile( $self->analysis_dir, $prev_stage_name,
                $component . '.bam' );
            my $prefix =
              $self->analysis->get_sample_name_by_bam_file_and_tag( $bam_file,
                $tag );
            push @all_parameters, [ $input_bam_file, $prefix ];
        }
    }

    return @all_parameters;
}

=method run_sample_flagstats_by_tag

  Usage       : run_sample_flagstats_by_tag();
  Purpose     : Run function for sample_flagstats_by_tag stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sample_flagstats_by_tag {
    my ( $self, $job ) = @_;

    return $self->run_sample_flagstats($job);
}

=method all_parameters_for_sample_flagstats_all_reads

  Usage       : all_parameters_for_sample_flagstats_all_reads();
  Purpose     : Get all parameters for sample_flagstats_all_reads stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_sample_flagstats_all_reads {
    my ( $self, $stage ) = @_;

    # Previous stage will be skipped if using native mark duplicates method
    my $prev_stage_name =
      $self->analysis->mark_duplicates_method eq 'native'
      ? 'sort_by_coordinate_all_reads'
      : 'mark_duplicates_all_reads';

    my @all_parameters;

    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $input_bam_file =
          File::Spec->catfile( $self->analysis_dir, $prev_stage_name,
            $component . '.bam' );
        my $prefix = join q{_},
          $self->analysis->get_sample_names_by_bam_file($bam_file);
        push @all_parameters, [ $input_bam_file, $prefix ];
    }

    return @all_parameters;
}

=method run_sample_flagstats_all_reads

  Usage       : run_sample_flagstats_all_reads();
  Purpose     : Run function for sample_flagstats_all_reads stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sample_flagstats_all_reads {
    my ( $self, $job ) = @_;

    return $self->run_sample_flagstats($job);
}

=method run_sample_flagstats

  Usage       : run_sample_flagstats();
  Purpose     : Run function for sample_flagstats stages
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_sample_flagstats {
    my ( $self, $job ) = @_;

    my ( $bam_file, $prefix ) = @{ $job->parameters };

    my $flagstat_output_file =
      File::Spec->catfile( $job->base_filename, $prefix . '.flagstat.txt' );

    # Get stats
    flagstats(
        {
            dir             => $job->base_filename,
            bam_file        => $bam_file,
            samtools_binary => $self->analysis->samtools_binary,
            output_file     => $flagstat_output_file,
        }
    );

    my $output_file = $job->base_filename . '.out';

    $self->dump_serialised( $output_file, \1 );

    return;
}

=method all_parameters_for_merge_by_tag

  Usage       : all_parameters_for_merge_by_tag();
  Purpose     : Get all parameters for merge_by_tag stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_merge_by_tag {
    my ( $self, $stage ) = @_;

    # Previous stage will be skipped if using native mark duplicates method
    my $prev_stage_name =
      $self->analysis->mark_duplicates_method eq 'native'
      ? 'sort_by_coordinate_by_tag'
      : 'mark_duplicates_by_tag';

    my @all_parameters;

    my @input_bam_files;
    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        my @tags = $self->analysis->list_all_tags_by_bam_file($bam_file);
        foreach my $tag (@tags) {
            $component++;
            my $input_bam_file = File::Spec->catfile( $self->analysis_dir,
                $prev_stage_name, $component . '.bam' );
            push @input_bam_files, $input_bam_file;
        }
    }
    push @all_parameters, \@input_bam_files;

    return @all_parameters;
}

=method run_merge_by_tag

  Usage       : run_merge_by_tag();
  Purpose     : Run function for merge_by_tag stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_merge_by_tag {
    my ( $self, $job ) = @_;

    return $self->run_merge($job);
}

=method all_parameters_for_merge_all_reads

  Usage       : all_parameters_for_merge_all_reads();
  Purpose     : Get all parameters for merge_all_reads stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_merge_all_reads {
    my ( $self, $stage ) = @_;

    # Previous stage will be skipped if using native mark duplicates method
    my $prev_stage_name =
      $self->analysis->mark_duplicates_method eq 'native'
      ? 'sort_by_coordinate_all_reads'
      : 'mark_duplicates_all_reads';

    my @all_parameters;

    my @input_bam_files;
    my $component = 0;
    foreach my $bam_file ( $self->analysis->list_all_bam_files() ) {
        $component++;
        my $input_bam_file = File::Spec->catfile( $self->analysis_dir,
            $prev_stage_name, $component . '.bam' );
        push @input_bam_files, $input_bam_file;
    }
    push @all_parameters, \@input_bam_files;

    return @all_parameters;
}

=method run_merge_all_reads

  Usage       : run_merge_all_reads();
  Purpose     : Run function for merge_all_reads stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_merge_all_reads {
    my ( $self, $job ) = @_;

    return $self->run_merge($job);
}

=method run_merge

  Usage       : run_merge();
  Purpose     : Run function for merge stages
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

    $self->dump_serialised( $output_file, \1 );

    return;
}

=method all_parameters_for_make_index

  Usage       : all_parameters_for_make_index();
  Purpose     : Get all parameters for make_index stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_make_index {
    my ( $self, $stage ) = @_;

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

    $self->dump_serialised( $output_file, \1 );

    return;
}

=method all_parameters_for_merged_flagstats

  Usage       : all_parameters_for_merged_flagstats();
  Purpose     : Get all parameters for merged_flagstats stage
  Returns     : Array of arrayrefs
  Parameters  : DETCT::Pipeline::Stage
  Throws      : No exceptions
  Comments    : None

=cut

sub all_parameters_for_merged_flagstats {
    my ( $self, $stage ) = @_;

    return [];
}

=method run_merged_flagstats

  Usage       : run_merged_flagstats();
  Purpose     : Run function for merged_flagstats stage
  Returns     : undef
  Parameters  : DETCT::Pipeline::Job
  Throws      : No exceptions
  Comments    : None

=cut

sub run_merged_flagstats {
    my ( $self, $job ) = @_;

    my $bam_file = File::Spec->catfile( $self->analysis_dir,
        $self->analysis->name . '.bam' );

    my $flagstat_output_file = File::Spec->catfile( $job->base_filename,
        $self->analysis->name . '.flagstat.txt' );

    # Get stats
    flagstats(
        {
            dir             => $job->base_filename,
            bam_file        => $bam_file,
            samtools_binary => $self->analysis->samtools_binary,
            output_file     => $flagstat_output_file,
        }
    );

    my $output_file = $job->base_filename . '.out';

    $self->dump_serialised( $output_file, \1 );

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
