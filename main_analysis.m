%% Golf Launch Monitor - Radar DSP Analysis
% Milestone 1: Data Foundation & Analysis
% Author: Tanveer Hussain
% Project: Club/Ball Separation Algorithm for InnoSenT SMR333
%
% This script implements the complete data analysis pipeline for golf 
% launch monitor radar data processing and TrackMan validation.

clear; close all; clc;

%% Configuration Parameters
config = struct();
config.data_path = 'D:\MATLAB OR PYTHON  PROJECT\28-01-2024 DATA';
config.sampling_freq = 22700; % Hz (from JSON metadata)
config.speed_coef = 0.1388888888888889; % Frequency to velocity conversion
config.fft_size = 1024; % Points
config.data_step = 64; % Samples
config.channels = 4; % IF1 Upper, IF1 Lower, IF2 Upper, IF2 Lower

% Analysis parameters
config.stft_window = 256; % STFT window size
config.stft_overlap = 128; % STFT overlap
config.velocity_range = [-200, 200]; % mph
config.trackman_tolerance = 3; % mph for initial validation

%% Initialize Analysis Framework
fprintf('=== Golf Launch Monitor Radar DSP Analysis ===\n');
fprintf('Milestone 1: Data Foundation & Analysis\n\n');

% Create results structure
results = struct();
results.config = config;
results.shots_analyzed = 0;
results.parsing_success = 0;
results.trackman_comparison = [];
results.signal_quality = struct();

%% Step 1: Parse Dataset and Extract Shot Information
fprintf('Step 1: Parsing dataset and extracting shot information...\n');

[shot_data, parsing_stats] = parse_golf_dataset(config.data_path);
results.shots_analyzed = length(shot_data);
results.parsing_success = parsing_stats.success_rate;

fprintf('Dataset Analysis Complete:\n');
fprintf('- Total shots found: %d\n', results.shots_analyzed);
fprintf('- Parsing success rate: %.1f%%\n', results.parsing_success * 100);
fprintf('- Speed ranges: Ball %.1f-%.1f mph, Club %.1f-%.1f mph\n\n', ...
    parsing_stats.ball_speed_range, parsing_stats.club_speed_range);

%% Step 2: Baseline STFT/FFT Processing Pipeline
fprintf('Step 2: Establishing baseline STFT/FFT processing pipeline...\n');

% Select representative shots for analysis
test_shots = select_representative_shots(shot_data, 10);

% Process each test shot
stft_results = cell(length(test_shots), 1);
for i = 1:length(test_shots)
    fprintf('Processing shot %d/%d: %s\n', i, length(test_shots), test_shots(i).folder_name);
    
    % Load I/Q data
    iq_data = load_iq_data(test_shots(i).bin_file, config);
    
    % Apply baseline processing
    [stft_result, freq_axis, time_axis] = baseline_stft_processing(iq_data, config);
    
    % Store results
    stft_results{i} = struct();
    stft_results{i}.shot_info = test_shots(i);
    stft_results{i}.stft_data = stft_result;
    stft_results{i}.freq_axis = freq_axis;
    stft_results{i}.time_axis = time_axis;
    stft_results{i}.iq_data = iq_data;
end

fprintf('STFT processing complete for %d test shots\n\n', length(test_shots));

%% Step 3: Signal Quality Assessment and SNR Analysis
fprintf('Step 3: Performing signal quality assessment...\n');

signal_quality = analyze_signal_quality(stft_results, config);
results.signal_quality = signal_quality;

fprintf('Signal Quality Analysis:\n');
fprintf('- Average SNR: %.1f dB\n', signal_quality.avg_snr);
fprintf('- Noise floor: %.1f dB\n', signal_quality.noise_floor);
fprintf('- Peak detection success: %.1f%%\n', signal_quality.peak_detection_rate * 100);
fprintf('- Signal bandwidth: %.1f Hz\n\n', signal_quality.signal_bandwidth);

%% Step 4: TrackMan Validation Framework
fprintf('Step 4: Creating TrackMan validation framework...\n');

trackman_validation = validate_against_trackman(stft_results, config);
results.trackman_comparison = trackman_validation;

fprintf('TrackMan Validation Results:\n');
fprintf('- Shots with timing alignment: %.1f%%\n', trackman_validation.timing_alignment * 100);
fprintf('- Initial velocity accuracy: ±%.1f mph\n', trackman_validation.velocity_accuracy);
fprintf('- Correlation with TrackMan data: %.3f\n\n', trackman_validation.correlation);

%% Step 5: Generate Comprehensive Visualizations
fprintf('Step 5: Generating analysis visualizations...\n');

% Create visualization figures
create_analysis_visualizations(stft_results, signal_quality, trackman_validation, config);

fprintf('Visualizations saved to output directory\n\n');

%% Step 6: Generate Milestone 1 Report
fprintf('Step 6: Generating Milestone 1 completion report...\n');

report = generate_milestone1_report(results, stft_results, config);

% Save results
save('milestone1_results.mat', 'results', 'stft_results', 'config', 'report');

fprintf('=== Milestone 1 Complete ===\n');
fprintf('Acceptance Criteria Status:\n');
if report.criteria.dataset_parsed
    fprintf('✓ Dataset fully parsed: PASS\n');
else
    fprintf('✓ Dataset fully parsed: FAIL\n');
end

if report.criteria.fft_alignment
    fprintf('✓ FFT/spectrogram alignment: PASS\n');
else
    fprintf('✓ FFT/spectrogram alignment: FAIL\n');
end

if report.criteria.fft_alignment
    fprintf('✓ TrackMan timing correlation ≥90%%: PASS\n');
else
    fprintf('✓ TrackMan timing correlation ≥90%%: FAIL\n');
end

if all(struct2array(report.criteria))
    fprintf('\n🎉 MILESTONE 1 SUCCESSFULLY COMPLETED! 🎉\n');
else
    fprintf('\n⚠️  Some acceptance criteria not met. Review required.\n');
end

fprintf('\nNext Steps: Proceed to Milestone 2 - Basic Target Detection\n');
