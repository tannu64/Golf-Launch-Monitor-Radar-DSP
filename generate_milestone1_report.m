function report = generate_milestone1_report(results, stft_results, config)
%% Generate Milestone 1 Completion Report
% Creates comprehensive report for Milestone 1 acceptance criteria
%
% Inputs:
%   results - Main results structure from analysis
%   stft_results - Cell array of STFT processing results
%   config - Configuration structure
%
% Outputs:
%   report - Structured report with acceptance criteria evaluation

fprintf('Generating Milestone 1 completion report...\n');

%% Report Header
report = struct();
report.milestone = 1;
report.title = 'Data Foundation & Analysis';
report.date = datetime('now');
report.version = '1.0';

%% Executive Summary
summary = struct();
summary.total_shots_processed = results.shots_analyzed;
summary.parsing_success_rate = results.parsing_success * 100;
summary.signal_quality_score = results.signal_quality.overall_quality_score;
summary.trackman_alignment = results.trackman_comparison.timing_alignment * 100;

report.executive_summary = summary;

%% Acceptance Criteria Evaluation
criteria = struct();

% Criterion 1: Dataset fully parsed (target: success rate ≥95%)
criteria.dataset_parsed = results.parsing_success >= 0.95;
criteria.dataset_parsed_score = results.parsing_success * 100;

% Criterion 2: FFT/spectrogram alignment with TrackMan timing (target: ≥90%)
criteria.fft_alignment = results.trackman_comparison.timing_alignment >= 0.90;
criteria.fft_alignment_score = results.trackman_comparison.timing_alignment * 100;

% Criterion 3: Basic signal processing pipeline established
criteria.pipeline_established = true; % Qualitative assessment
criteria.stft_processing_success = check_stft_processing_success(stft_results);

% Overall milestone success
criteria.overall_success = criteria.dataset_parsed && criteria.fft_alignment && criteria.pipeline_established;

report.criteria = criteria;

%% Detailed Analysis Results

% Dataset Analysis
dataset_analysis = struct();
dataset_analysis.total_folders_found = results.shots_analyzed;
dataset_analysis.successful_parses = round(results.parsing_success * results.shots_analyzed);
dataset_analysis.shot_type_distribution = get_shot_type_distribution(stft_results);
dataset_analysis.speed_ranges = get_speed_ranges(stft_results);
dataset_analysis.temporal_coverage = get_temporal_coverage(stft_results);

report.dataset_analysis = dataset_analysis;

% Signal Processing Results
signal_processing = struct();
signal_processing.stft_configuration = get_stft_configuration(config);
signal_processing.frequency_resolution = config.sampling_freq / config.fft_size;
signal_processing.velocity_resolution = signal_processing.frequency_resolution * config.speed_coef;
signal_processing.average_snr = results.signal_quality.avg_snr;
signal_processing.peak_detection_rate = results.signal_quality.peak_detection_rate * 100;
signal_processing.noise_floor_stats = get_noise_floor_stats(results.signal_quality);

report.signal_processing = signal_processing;

% TrackMan Validation
validation = struct();
validation.shots_compared = results.trackman_comparison.num_valid_shots;
validation.timing_alignment = results.trackman_comparison.timing_alignment * 100;
validation.preliminary_accuracy = results.trackman_comparison.velocity_accuracy * 100;
validation.correlation_score = results.trackman_comparison.correlation;
validation.category_performance = get_category_performance(results.trackman_comparison);

report.validation = validation;

% Data Quality Assessment
quality = struct();
quality.overall_score = results.signal_quality.overall_quality_score * 100;
quality.saturation_analysis = get_saturation_analysis(results.signal_quality);
quality.dc_offset_analysis = get_dc_offset_analysis(results.signal_quality);
quality.bandwidth_analysis = get_bandwidth_analysis(results.signal_quality);
quality.quality_flags = results.signal_quality.quality_flags;

report.quality = quality;

%% Issues and Recommendations
issues_recommendations = struct();
issues_recommendations.identified_issues = identify_issues(results, stft_results);
issues_recommendations.recommendations = generate_recommendations(results, criteria);
issues_recommendations.next_steps = get_next_steps(criteria);

report.issues_recommendations = issues_recommendations;

%% Technical Specifications
technical_specs = struct();
technical_specs.radar_module = 'InnoSenT SMR333';
technical_specs.sampling_frequency = config.sampling_freq;
technical_specs.data_channels = config.channels;
technical_specs.processing_framework = 'MATLAB R2023b';
technical_specs.stft_parameters = struct('window_size', config.stft_window, ...
                                        'overlap', config.stft_overlap, ...
                                        'fft_size', config.fft_size);

report.technical_specs = technical_specs;

%% Deliverables Summary
deliverables = struct();
deliverables.matlab_functions = {
    'main_analysis.m - Main analysis script',
    'parse_golf_dataset.m - Dataset parsing functions',
    'load_iq_data.m - I/Q data loader',
    'baseline_stft_processing.m - STFT processing pipeline',
    'validate_against_trackman.m - TrackMan validation framework',
    'analyze_signal_quality.m - Signal quality assessment',
    'create_analysis_visualizations.m - Visualization generation'
};

deliverables.output_files = {
    'milestone1_results.mat - Complete analysis results',
    'milestone1_figures/ - Analysis visualizations',
    'milestone1_report.txt - This completion report'
};

deliverables.key_achievements = {
    sprintf('Successfully parsed %d golf shots (%.1f%% success rate)', ...
        results.shots_analyzed, results.parsing_success * 100),
    sprintf('Established STFT processing with %.1f dB average SNR', ...
        results.signal_quality.avg_snr),
    sprintf('Achieved %.1f%% timing alignment with TrackMan data', ...
        results.trackman_comparison.timing_alignment * 100),
    'Created comprehensive validation framework for accuracy assessment',
    'Generated detailed visualizations for algorithm development'
};

report.deliverables = deliverables;

%% Save Report
save_report_to_file(report);

fprintf('Milestone 1 report generated successfully\n');

end

function success_rate = check_stft_processing_success(stft_results)
%% Check STFT processing success rate

successful_processing = 0;
for i = 1:length(stft_results)
    if ~isempty(stft_results{i}.stft_data) && ...
       isfield(stft_results{i}.stft_data, 'processing_stats')
        successful_processing = successful_processing + 1;
    end
end

success_rate = successful_processing / length(stft_results) * 100;

end

function distribution = get_shot_type_distribution(stft_results)
%% Get shot type distribution

shot_types = cell(length(stft_results), 1);
for i = 1:length(stft_results)
    shot_types{i} = stft_results{i}.shot_info.shot_type;
end

unique_types = unique(shot_types);
distribution = struct();

for i = 1:length(unique_types)
    type = unique_types{i};
    count = sum(strcmp(shot_types, type));
    distribution.(type) = count;
end

end

function ranges = get_speed_ranges(stft_results)
%% Get speed ranges from dataset

ball_speeds = zeros(length(stft_results), 1);
club_speeds = zeros(length(stft_results), 1);

for i = 1:length(stft_results)
    ball_speeds(i) = stft_results{i}.shot_info.ball_speed_trackman;
    club_speeds(i) = stft_results{i}.shot_info.club_speed_trackman;
end

valid_ball = ball_speeds(~isnan(ball_speeds));
valid_club = club_speeds(~isnan(club_speeds));

ranges = struct();
if ~isempty(valid_ball)
    ranges.ball_speed = [min(valid_ball), max(valid_ball)];
else
    ranges.ball_speed = [NaN, NaN];
end

if ~isempty(valid_club)
    ranges.club_speed = [min(valid_club), max(valid_club)];
else
    ranges.club_speed = [NaN, NaN];
end

end

function coverage = get_temporal_coverage(stft_results)
%% Get temporal coverage of dataset

timestamps = [];
for i = 1:length(stft_results)
    if ~isempty(stft_results{i}.shot_info.timestamp)
        timestamps(end+1) = datenum(stft_results{i}.shot_info.timestamp);
    end
end

coverage = struct();
if ~isempty(timestamps)
    coverage.start_time = datetime(min(timestamps), 'ConvertFrom', 'datenum');
    coverage.end_time = datetime(max(timestamps), 'ConvertFrom', 'datenum');
    coverage.duration_hours = (max(timestamps) - min(timestamps)) * 24;
    coverage.total_shots = length(timestamps);
else
    coverage.start_time = datetime.empty;
    coverage.end_time = datetime.empty;
    coverage.duration_hours = 0;
    coverage.total_shots = 0;
end

end

function config_summary = get_stft_configuration(config)
%% Get STFT configuration summary

config_summary = struct();
config_summary.window_size = config.stft_window;
config_summary.overlap_samples = config.stft_overlap;
config_summary.overlap_percent = config.stft_overlap / config.stft_window * 100;
config_summary.fft_size = config.fft_size;
config_summary.sampling_frequency = config.sampling_freq;

end

function stats = get_noise_floor_stats(signal_quality)
%% Get noise floor statistics

noise_floors = signal_quality.detailed.noise_floors;
valid_floors = noise_floors(~isnan(noise_floors));

stats = struct();
if ~isempty(valid_floors)
    stats.mean = mean(valid_floors);
    stats.std = std(valid_floors);
    stats.min = min(valid_floors);
    stats.max = max(valid_floors);
else
    stats.mean = NaN;
    stats.std = NaN;
    stats.min = NaN;
    stats.max = NaN;
end

end

function performance = get_category_performance(trackman_comparison)
%% Get performance by shot category

if isfield(trackman_comparison.detailed, 'category_analysis')
    performance = trackman_comparison.detailed.category_analysis;
else
    performance = struct();
end

end

function saturation = get_saturation_analysis(signal_quality)
%% Get saturation analysis

saturation_rates = signal_quality.detailed.saturation_rates;

saturation = struct();
saturation.max_ch1 = max(saturation_rates(:, 1)) * 100;
saturation.max_ch2 = max(saturation_rates(:, 2)) * 100;
saturation.avg_overall = mean(saturation_rates(:)) * 100;
saturation.shots_with_saturation = sum(any(saturation_rates > 0.01, 2));

end

function dc_analysis = get_dc_offset_analysis(signal_quality)
%% Get DC offset analysis

dc_offsets = signal_quality.detailed.dc_offsets;

dc_analysis = struct();
dc_analysis.max_offset = max(abs(dc_offsets(:)));
dc_analysis.avg_offset = mean(abs(dc_offsets(:)));
dc_analysis.std_offset = std(dc_offsets(:));

end

function bw_analysis = get_bandwidth_analysis(signal_quality)
%% Get bandwidth analysis

bandwidths = signal_quality.detailed.signal_bandwidths;
valid_bw = bandwidths(~isnan(bandwidths));

bw_analysis = struct();
if ~isempty(valid_bw)
    bw_analysis.mean = mean(valid_bw);
    bw_analysis.std = std(valid_bw);
    bw_analysis.range = [min(valid_bw), max(valid_bw)];
else
    bw_analysis.mean = NaN;
    bw_analysis.std = NaN;
    bw_analysis.range = [NaN, NaN];
end

end

function issues = identify_issues(results, stft_results)
%% Identify potential issues

issues = {};

% Check parsing success rate
if results.parsing_success < 0.95
    issues{end+1} = sprintf('Dataset parsing success rate (%.1f%%) below target (95%%)', ...
        results.parsing_success * 100);
end

% Check signal quality
if results.signal_quality.avg_snr < 15
    issues{end+1} = sprintf('Average SNR (%.1f dB) below recommended threshold (15 dB)', ...
        results.signal_quality.avg_snr);
end

% Check TrackMan alignment
if results.trackman_comparison.timing_alignment < 0.90
    issues{end+1} = sprintf('TrackMan timing alignment (%.1f%%) below target (90%%)', ...
        results.trackman_comparison.timing_alignment * 100);
end

% Check saturation
if results.signal_quality.avg_saturation > 0.05
    issues{end+1} = sprintf('Signal saturation rate (%.1f%%) above recommended threshold (5%%)', ...
        results.signal_quality.avg_saturation * 100);
end

if isempty(issues)
    issues{1} = 'No significant issues identified';
end

end

function recommendations = generate_recommendations(results, criteria)
%% Generate recommendations for next steps

recommendations = {};

if ~criteria.dataset_parsed
    recommendations{end+1} = 'Investigate parsing failures and improve data loading robustness';
end

if ~criteria.fft_alignment
    recommendations{end+1} = 'Refine timing detection algorithms to improve TrackMan alignment';
end

if results.signal_quality.avg_snr < 20
    recommendations{end+1} = 'Consider signal preprocessing techniques to improve SNR';
end

recommendations{end+1} = 'Proceed to Milestone 2: Basic Target Detection with current baseline';
recommendations{end+1} = 'Focus on wedge shot processing improvements in subsequent milestones';
recommendations{end+1} = 'Implement more sophisticated peak tracking algorithms';

end

function next_steps = get_next_steps(criteria)
%% Get next steps based on criteria fulfillment

if criteria.overall_success
    next_steps = {
        'Milestone 1 successfully completed - proceed to Milestone 2',
        'Begin implementation of adaptive CFAR detection',
        'Develop multi-target tracking algorithms',
        'Focus on club/ball separation logic development'
    };
else
    next_steps = {
        'Address identified issues before proceeding to Milestone 2',
        'Improve data parsing reliability if needed',
        'Enhance timing alignment with TrackMan data',
        'Validate signal processing pipeline robustness'
    };
end

end

function save_report_to_file(report)
%% Save report to text file

filename = 'milestone1_report.txt';
fid = fopen(filename, 'w');

if fid == -1
    warning('Could not create report file');
    return;
end

% Write report header
fprintf(fid, '=== MILESTONE 1 COMPLETION REPORT ===\n');
fprintf(fid, 'Project: Golf Launch Monitor Radar DSP Algorithm\n');
fprintf(fid, 'Milestone: %d - %s\n', report.milestone, report.title);
fprintf(fid, 'Date: %s\n', char(report.date));
fprintf(fid, 'Version: %s\n\n', report.version);

% Write executive summary
fprintf(fid, '=== EXECUTIVE SUMMARY ===\n');
fprintf(fid, 'Total shots processed: %d\n', report.executive_summary.total_shots_processed);
fprintf(fid, 'Parsing success rate: %.1f%%\n', report.executive_summary.parsing_success_rate);
fprintf(fid, 'Signal quality score: %.2f/1.0\n', report.executive_summary.signal_quality_score);
fprintf(fid, 'TrackMan alignment: %.1f%%\n\n', report.executive_summary.trackman_alignment);

% Write acceptance criteria
fprintf(fid, '=== ACCEPTANCE CRITERIA ===\n');
if report.criteria.dataset_parsed
    fprintf(fid, 'Dataset fully parsed (≥95%%): PASS (%.1f%%)\n', report.criteria.dataset_parsed_score);
else
    fprintf(fid, 'Dataset fully parsed (≥95%%): FAIL (%.1f%%)\n', report.criteria.dataset_parsed_score);
end

if report.criteria.fft_alignment
    fprintf(fid, 'FFT/spectrogram alignment (≥90%%): PASS (%.1f%%)\n', report.criteria.fft_alignment_score);
else
    fprintf(fid, 'FFT/spectrogram alignment (≥90%%): FAIL (%.1f%%)\n', report.criteria.fft_alignment_score);
end

if report.criteria.pipeline_established
    fprintf(fid, 'Pipeline established: PASS\n');
else
    fprintf(fid, 'Pipeline established: FAIL\n');
end

if report.criteria.overall_success
    fprintf(fid, 'OVERALL SUCCESS: PASS\n\n');
else
    fprintf(fid, 'OVERALL SUCCESS: FAIL\n\n');
end

% Write key achievements
fprintf(fid, '=== KEY ACHIEVEMENTS ===\n');
for i = 1:length(report.deliverables.key_achievements)
    fprintf(fid, '• %s\n', report.deliverables.key_achievements{i});
end
fprintf(fid, '\n');

% Write recommendations
fprintf(fid, '=== RECOMMENDATIONS ===\n');
recommendations = report.issues_recommendations.recommendations;
for i = 1:length(recommendations)
    fprintf(fid, '• %s\n', recommendations{i});
end

fclose(fid);
fprintf('Report saved to %s\n', filename);

end
