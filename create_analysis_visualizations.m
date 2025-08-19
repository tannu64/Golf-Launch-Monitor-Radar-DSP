function create_analysis_visualizations(stft_results, signal_quality, trackman_validation, config)
%% Create Analysis Visualizations
% Generates comprehensive plots for Milestone 1 analysis
%
% Inputs:
%   stft_results - Cell array of STFT processing results
%   signal_quality - Signal quality analysis structure
%   trackman_validation - TrackMan validation results
%   config - Configuration structure

fprintf('Generating analysis visualizations...\n');

% Create output directory for figures
output_dir = 'milestone1_figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Figure 1: Representative Spectrograms
create_spectrogram_figure(stft_results, output_dir);

%% Figure 2: Signal Quality Summary
create_signal_quality_figure(signal_quality, output_dir);

%% Figure 3: TrackMan Validation
create_trackman_validation_figure(trackman_validation, output_dir);

%% Figure 4: Velocity Detection Overview
create_velocity_detection_figure(stft_results, output_dir);

%% Figure 5: Data Quality Assessment
create_data_quality_figure(stft_results, signal_quality, output_dir);

fprintf('Visualizations saved to %s directory\n', output_dir);

end

function create_spectrogram_figure(stft_results, output_dir)
%% Create spectrogram visualization for representative shots

fig1 = figure('Position', [100, 100, 1200, 800]);

% Select 4 representative shots for display
num_plots = min(4, length(stft_results));
shot_indices = round(linspace(1, length(stft_results), num_plots));

for i = 1:num_plots
    subplot(2, 2, i);
    
    shot_idx = shot_indices(i);
    stft_data = stft_results{shot_idx}.stft_data;
    shot_info = stft_results{shot_idx}.shot_info;
    
    % Plot spectrogram
    imagesc(stft_data.processing_stats.time_resolution * (0:size(stft_data.magnitude_db_ch1, 2)-1), ...
            stft_data.velocity_axis, stft_data.magnitude_db_ch1);
    
    colormap('jet');
    colorbar;
    axis xy;
    
    title(sprintf('Shot %d: %s\nBall: %.1f mph, Club: %.1f mph', ...
        shot_idx, shot_info.shot_type, shot_info.ball_speed_trackman, shot_info.club_speed_trackman));
    xlabel('Time (s)');
    ylabel('Velocity (mph)');
    ylim([-50, 200]);
    
    % Overlay detected peaks
    hold on;
    if ~isempty(stft_data.peaks.times)
        scatter(stft_data.peaks.times, stft_data.peaks.velocities, 20, 'w', 'filled');
    end
    hold off;
end

sgtitle('Representative Golf Shot Spectrograms');
saveas(fig1, fullfile(output_dir, 'spectrograms.png'));
close(fig1);

end

function create_signal_quality_figure(signal_quality, output_dir)
%% Create signal quality summary figure

fig2 = figure('Position', [200, 200, 1000, 600]);

% SNR distribution
subplot(2, 3, 1);
histogram(signal_quality.detailed.snr_values(~isnan(signal_quality.detailed.snr_values)), 10);
title('SNR Distribution');
xlabel('SNR (dB)');
ylabel('Count');
grid on;

% Peak detection counts
subplot(2, 3, 2);
histogram(signal_quality.detailed.peak_counts, 0:max(signal_quality.detailed.peak_counts));
title('Peak Detection Counts');
xlabel('Number of Peaks');
ylabel('Count');
grid on;

% Saturation rates
subplot(2, 3, 3);
saturation_data = signal_quality.detailed.saturation_rates * 100; % Convert to percentage
bar([mean(saturation_data(:, 1)), mean(saturation_data(:, 2))]);
title('Average Saturation Rates');
xlabel('Channel');
ylabel('Saturation (%)');
xticklabels({'Ch1', 'Ch2'});
grid on;

% Signal bandwidth distribution
subplot(2, 3, 4);
valid_bw = signal_quality.detailed.signal_bandwidths(~isnan(signal_quality.detailed.signal_bandwidths));
if ~isempty(valid_bw)
    histogram(valid_bw, 10);
end
title('Signal Bandwidth Distribution');
xlabel('Bandwidth (Hz)');
ylabel('Count');
grid on;

% Quality flags summary
subplot(2, 3, 5);
flags = signal_quality.quality_flags;
flag_names = fieldnames(flags);
flag_values = cellfun(@(x) flags.(x), flag_names);
bar(flag_values);
title('Quality Assessment Flags');
ylabel('Pass (1) / Fail (0)');
xticklabels(flag_names);
xtickangle(45);
grid on;

% Overall quality score
subplot(2, 3, 6);
pie([signal_quality.overall_quality_score, 1-signal_quality.overall_quality_score], ...
    {'Good Quality', 'Issues Detected'});
title(sprintf('Overall Quality Score: %.2f', signal_quality.overall_quality_score));

sgtitle('Signal Quality Analysis Summary');
saveas(fig2, fullfile(output_dir, 'signal_quality.png'));
close(fig2);

end

function create_trackman_validation_figure(trackman_validation, output_dir)
%% Create TrackMan validation figure

fig3 = figure('Position', [300, 300, 1200, 800]);

% Ball speed comparison
subplot(2, 3, 1);
ball_est = trackman_validation.detailed.ball_speed_estimated;
ball_ref = trackman_validation.detailed.ball_speed_trackman;
valid_mask = ~isnan(ball_est) & ~isnan(ball_ref);

if sum(valid_mask) > 0
    scatter(ball_ref(valid_mask), ball_est(valid_mask), 50, 'b', 'filled');
    hold on;
    plot([0, max(ball_ref(valid_mask))], [0, max(ball_ref(valid_mask))], 'r--', 'LineWidth', 2);
    hold off;
end

title('Ball Speed Validation');
xlabel('TrackMan Ball Speed (mph)');
ylabel('Estimated Ball Speed (mph)');
grid on;
axis equal;

% Club speed comparison
subplot(2, 3, 2);
club_est = trackman_validation.detailed.club_speed_estimated;
club_ref = trackman_validation.detailed.club_speed_trackman;
valid_mask = ~isnan(club_est) & ~isnan(club_ref);

if sum(valid_mask) > 0
    scatter(club_ref(valid_mask), club_est(valid_mask), 50, 'g', 'filled');
    hold on;
    plot([0, max(club_ref(valid_mask))], [0, max(club_ref(valid_mask))], 'r--', 'LineWidth', 2);
    hold off;
end

title('Club Speed Validation');
xlabel('TrackMan Club Speed (mph)');
ylabel('Estimated Club Speed (mph)');
grid on;
axis equal;

% Error distribution
subplot(2, 3, 3);
ball_errors = abs(ball_est(valid_mask) - ball_ref(valid_mask));
club_errors = abs(club_est(valid_mask) - club_ref(valid_mask));

if ~isempty(ball_errors) && ~isempty(club_errors)
    histogram(ball_errors, 10, 'FaceAlpha', 0.5, 'DisplayName', 'Ball Speed Errors');
    hold on;
    histogram(club_errors, 10, 'FaceAlpha', 0.5, 'DisplayName', 'Club Speed Errors');
    hold off;
    legend;
end

title('Speed Estimation Errors');
xlabel('Absolute Error (mph)');
ylabel('Count');
grid on;

% Validation metrics summary
subplot(2, 3, 4);
metrics = [trackman_validation.timing_alignment, ...
           trackman_validation.ball_accuracy_3mph, ...
           trackman_validation.club_accuracy_3mph] * 100;
bar(metrics);
title('Validation Metrics');
ylabel('Percentage (%)');
xticklabels({'Timing Align', 'Ball ±3mph', 'Club ±3mph'});
grid on;
ylim([0, 100]);

% Category analysis
subplot(2, 3, 5);
if isfield(trackman_validation.detailed, 'category_analysis')
    cat_analysis = trackman_validation.detailed.category_analysis;
    categories = fieldnames(cat_analysis);
    
    ball_rms = zeros(length(categories), 1);
    club_rms = zeros(length(categories), 1);
    
    for i = 1:length(categories)
        if isfield(cat_analysis.(categories{i}), 'ball_rms_error')
            ball_rms(i) = cat_analysis.(categories{i}).ball_rms_error;
            club_rms(i) = cat_analysis.(categories{i}).club_rms_error;
        end
    end
    
    x = 1:length(categories);
    bar(x-0.2, ball_rms, 0.4, 'DisplayName', 'Ball RMS Error');
    hold on;
    bar(x+0.2, club_rms, 0.4, 'DisplayName', 'Club RMS Error');
    hold off;
    
    title('RMS Error by Shot Type');
    xlabel('Shot Type');
    ylabel('RMS Error (mph)');
    xticklabels(categories);
    legend;
    grid on;
end

% Correlation scores
subplot(2, 3, 6);
corr_scores = trackman_validation.detailed.correlation_scores;
valid_corr = corr_scores(~isnan(corr_scores));

if ~isempty(valid_corr)
    histogram(valid_corr, 10);
end

title('Correlation Score Distribution');
xlabel('Correlation Score');
ylabel('Count');
grid on;

sgtitle('TrackMan Validation Results');
saveas(fig3, fullfile(output_dir, 'trackman_validation.png'));
close(fig3);

end

function create_velocity_detection_figure(stft_results, output_dir)
%% Create velocity detection overview figure

fig4 = figure('Position', [400, 400, 1000, 600]);

% Compile velocity data from all shots
all_velocities = [];
all_magnitudes = [];
shot_types = {};

for i = 1:length(stft_results)
    if ~isempty(stft_results{i}.stft_data.peaks.velocities)
        velocities = abs(stft_results{i}.stft_data.peaks.velocities);
        magnitudes = stft_results{i}.stft_data.peaks.magnitudes;
        shot_type = stft_results{i}.shot_info.shot_type;
        
        all_velocities = [all_velocities; velocities(:)];
        all_magnitudes = [all_magnitudes; magnitudes(:)];
        shot_types = [shot_types; repmat({shot_type}, length(velocities), 1)];
    end
end

% Velocity distribution by shot type
subplot(2, 2, 1);
if ~isempty(all_velocities)
    wedge_vel = all_velocities(strcmp(shot_types, 'wedge'));
    iron_vel = all_velocities(strcmp(shot_types, 'iron'));
    driver_vel = all_velocities(strcmp(shot_types, 'driver'));
    
    histogram(wedge_vel, 20, 'FaceAlpha', 0.5, 'DisplayName', 'Wedge');
    hold on;
    histogram(iron_vel, 20, 'FaceAlpha', 0.5, 'DisplayName', 'Iron');
    histogram(driver_vel, 20, 'FaceAlpha', 0.5, 'DisplayName', 'Driver');
    hold off;
    legend;
end

title('Detected Velocity Distribution');
xlabel('Velocity (mph)');
ylabel('Count');
grid on;

% Velocity vs Magnitude scatter
subplot(2, 2, 2);
if ~isempty(all_velocities)
    scatter(all_velocities, all_magnitudes, 20, 'filled');
end
title('Velocity vs Signal Magnitude');
xlabel('Velocity (mph)');
ylabel('Magnitude (dB)');
grid on;

% Track statistics
subplot(2, 2, 3);
track_counts = zeros(length(stft_results), 1);
track_durations = [];

for i = 1:length(stft_results)
    tracks = stft_results{i}.stft_data.tracked_objects;
    track_counts(i) = length(tracks);
    
    if ~isempty(tracks)
        track_durations = [track_durations, [tracks.duration]];
    end
end

histogram(track_counts, 0:max(track_counts));
title('Number of Tracks per Shot');
xlabel('Track Count');
ylabel('Number of Shots');
grid on;

% Track duration distribution
subplot(2, 2, 4);
if ~isempty(track_durations)
    histogram(track_durations, 20);
end
title('Track Duration Distribution');
xlabel('Duration (s)');
ylabel('Count');
grid on;

sgtitle('Velocity Detection Analysis');
saveas(fig4, fullfile(output_dir, 'velocity_detection.png'));
close(fig4);

end

function create_data_quality_figure(stft_results, signal_quality, output_dir)
%% Create data quality assessment figure

fig5 = figure('Position', [500, 500, 1000, 600]);

% DC offset analysis
subplot(2, 3, 1);
dc_data = signal_quality.detailed.dc_offsets;
boxplot(dc_data, 'Labels', {'Ch1_I', 'Ch1_Q', 'Ch2_I', 'Ch2_Q'});
title('DC Offset Distribution');
ylabel('DC Offset Value');
grid on;

% File size vs shot type
subplot(2, 3, 2);
file_sizes = zeros(length(stft_results), 1);
shot_types_list = cell(length(stft_results), 1);

for i = 1:length(stft_results)
    file_info = dir(stft_results{i}.shot_info.bin_file);
    if ~isempty(file_info)
        file_sizes(i) = file_info.bytes / 1024 / 1024; % MB
    end
    shot_types_list{i} = stft_results{i}.shot_info.shot_type;
end

unique_types = unique(shot_types_list);
colors = lines(length(unique_types));

for i = 1:length(unique_types)
    mask = strcmp(shot_types_list, unique_types{i});
    scatter(find(mask), file_sizes(mask), 50, colors(i, :), 'filled', 'DisplayName', unique_types{i});
    hold on;
end
hold off;

title('File Size by Shot Type');
xlabel('Shot Index');
ylabel('File Size (MB)');
legend;
grid on;

% Signal duration analysis
subplot(2, 3, 3);
durations = zeros(length(stft_results), 1);
for i = 1:length(stft_results)
    durations(i) = stft_results{i}.iq_data.duration;
end

histogram(durations, 20);
title('Signal Duration Distribution');
xlabel('Duration (s)');
ylabel('Count');
grid on;

% Power levels by channel
subplot(2, 3, 4);
power_ch1 = zeros(length(stft_results), 1);
power_ch2 = zeros(length(stft_results), 1);

for i = 1:length(stft_results)
    power_ch1(i) = 10*log10(stft_results{i}.iq_data.power_ch1);
    power_ch2(i) = 10*log10(stft_results{i}.iq_data.power_ch2);
end

plot(power_ch1, 'b.-', 'DisplayName', 'Channel 1');
hold on;
plot(power_ch2, 'r.-', 'DisplayName', 'Channel 2');
hold off;

title('Signal Power Levels');
xlabel('Shot Index');
ylabel('Power (dB)');
legend;
grid on;

% Temporal analysis
subplot(2, 3, 5);
timestamps = [];
for i = 1:length(stft_results)
    if ~isempty(stft_results{i}.shot_info.timestamp)
        timestamps(end+1) = datenum(stft_results{i}.shot_info.timestamp);
    end
end

if ~isempty(timestamps)
    plot(timestamps - min(timestamps), 1:length(timestamps), 'o-');
    title('Shot Sequence Timeline');
    xlabel('Time from First Shot (days)');
    ylabel('Shot Number');
    grid on;
end

% Summary statistics
subplot(2, 3, 6);
summary_stats = [
    signal_quality.overall_quality_score;
    signal_quality.peak_detection_rate;
    1 - signal_quality.avg_saturation; % Convert to "good" metric
    signal_quality.avg_snr / 30; % Normalize to ~1
];

bar(summary_stats);
title('Overall Data Quality Metrics');
ylabel('Score (0-1)');
xticklabels({'Overall', 'Peak Detect', 'Low Saturation', 'Good SNR'});
xtickangle(45);
grid on;
ylim([0, 1]);

sgtitle('Data Quality Assessment');
saveas(fig5, fullfile(output_dir, 'data_quality.png'));
close(fig5);

end
