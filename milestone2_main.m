%% Milestone 2: Enhanced Golf Launch Monitor DSP Pipeline
% Builds on Milestone 1 foundation with advanced CFAR detection,
% improved tracking, and velocity smoothing
%
% Author: Tanveer Hussain
% Project: Golf Launch Monitor Club/Ball Separation & Accuracy Improvement

clear; close all; clc;
warning off;

fprintf('=== MILESTONE 2: ENHANCED DSP PROCESSING ===\n');
fprintf('Building on Milestone 1 foundation...\n\n');

%% Step 1: Load Milestone 1 Results
fprintf('Step 1: Loading Milestone 1 baseline results...\n');

% Check if M1 results exist
if exist('milestone1_results.mat', 'file')
    load('milestone1_results.mat');
    fprintf('✓ Loaded M1 results: %d shots analyzed\n', length(stft_results));
else
    % Run M1 analysis if results don't exist
    fprintf('M1 results not found. Running baseline analysis...\n');
    run('main_analysis.m');
    fprintf('✓ M1 baseline complete\n');
end

%% Step 2: Enhanced Processing Pipeline
fprintf('\nStep 2: Applying Milestone 2 enhancements...\n');

% M2-specific configuration parameters
config.cfar_pfa = 1e-4; % Probability of false alarm
config.cfar_guard_cells = 3;
config.cfar_reference_cells = 16;
config.m2_velocity_tolerance = 5; % mph (tighter than M1)
config.m2_min_track_length = 3; % detections
config.m2_max_time_gap = 0.05; % seconds
config.m2_quality_threshold = 0.4; % track quality threshold
config.m2_process_noise = 25; % Kalman filter Q
config.m2_measurement_noise = 9; % Kalman filter R
config.m2_initial_uncertainty = 100; % Kalman filter P0

% Process each shot through M2 pipeline
num_shots = length(stft_results);
enhanced_results = cell(num_shots, 1);

fprintf('Processing %d shots through M2 pipeline:\n', num_shots);

for i = 1:num_shots
    fprintf('[%d/%d] Processing: %s\n', i, num_shots, stft_results{i}.shot_info.folder_name);
    
    try
        % Step 2a: Apply CFAR detection
        [cfar_result, cfar_stats] = milestone2_cfar_detection(stft_results{i}.stft_data, config);
        
        % Step 2b: Advanced tracking
        [advanced_tracks, tracking_stats] = milestone2_advanced_tracking(cfar_result.peaks_cfar, config);
        
        % Step 2c: Velocity smoothing (if tracks exist)
        if ~isempty(advanced_tracks)
            [smoothed_tracks, smoothing_stats] = milestone2_velocity_smoothing(advanced_tracks, config);
        else
            smoothed_tracks = [];
            smoothing_stats = struct('tracks_processed', 0, 'avg_improvement', 0);
        end
        
        % Store enhanced results
        enhanced_results{i} = struct();
        enhanced_results{i}.shot_info = stft_results{i}.shot_info;
        enhanced_results{i}.m1_result = stft_results{i}.stft_data; % Keep M1 for comparison
        enhanced_results{i}.cfar_result = cfar_result;
        enhanced_results{i}.advanced_tracks = smoothed_tracks;
        enhanced_results{i}.stats = struct();
        enhanced_results{i}.stats.cfar = cfar_stats;
        enhanced_results{i}.stats.tracking = tracking_stats;
        enhanced_results{i}.stats.smoothing = smoothing_stats;
        
    catch ME
        fprintf('   ⚠️ Error processing shot %d: %s\n', i, ME.message);
        
        % Create minimal result structure for failed shots
        enhanced_results{i} = struct();
        enhanced_results{i}.shot_info = stft_results{i}.shot_info;
        enhanced_results{i}.m1_result = stft_results{i}.stft_data;
        enhanced_results{i}.cfar_result = [];
        enhanced_results{i}.advanced_tracks = [];
        enhanced_results{i}.stats = struct();
        enhanced_results{i}.error = ME.message;
    end
end

fprintf('✓ M2 processing pipeline complete\n\n');

%% Step 3: Milestone 2 Validation
fprintf('Step 3: Performing Milestone 2 validation against TrackMan...\n');

m2_validation = milestone2_validation(enhanced_results, config);

%% Step 4: Generate M2 Performance Report
fprintf('\nStep 4: Generating Milestone 2 performance report...\n');

% Create comprehensive M2 results structure
m2_results = struct();
m2_results.config = config;
m2_results.enhanced_results = enhanced_results;
m2_results.validation = m2_validation;
m2_results.processing_date = datetime('now');
m2_results.shots_processed = num_shots;

% Calculate summary statistics
valid_shots = m2_validation.num_valid_shots;
detection_rate = m2_validation.detection_rate;
velocity_accuracy = m2_validation.velocity_accuracy;
ball_rms = m2_validation.ball_rms_error;
club_rms = m2_validation.club_rms_error;

%% Step 5: Display Results Summary
fprintf('\n=== MILESTONE 2 RESULTS SUMMARY ===\n');
fprintf('Processing Summary:\n');
fprintf('- Total shots processed: %d\n', num_shots);
fprintf('- Valid shots for validation: %d\n', valid_shots);
fprintf('- Processing success rate: %.1f%%\n', (num_shots - sum(cellfun(@(x) isfield(x, 'error'), enhanced_results))) / num_shots * 100);

fprintf('\nPerformance Metrics:\n');
fprintf('- Detection rate: %.1f%% (Target: ≥90%%)\n', detection_rate * 100);
fprintf('- Velocity accuracy (±3 mph): %.1f%% (Target: ≥90%%)\n', velocity_accuracy * 100);
fprintf('- Ball speed RMS error: %.2f mph (Target: <3 mph)\n', ball_rms);
fprintf('- Club speed RMS error: %.2f mph (Target: <3 mph)\n', club_rms);

fprintf('\nAccuracy Breakdown:\n');
fprintf('- Ball speed accuracy (±3 mph): %.1f%%\n', m2_validation.ball_accuracy_3mph * 100);
fprintf('- Club speed accuracy (±3 mph): %.1f%%\n', m2_validation.club_accuracy_3mph * 100);
fprintf('- Ball speed accuracy (±1 mph): %.1f%%\n', m2_validation.ball_accuracy_1mph * 100);
fprintf('- Club speed accuracy (±1 mph): %.1f%%\n', m2_validation.club_accuracy_1mph * 100);

%% Step 6: Acceptance Criteria Assessment
fprintf('\n=== MILESTONE 2 ACCEPTANCE CRITERIA ===\n');

criteria = m2_validation.criteria_met;

if criteria.detection_rate
    fprintf('✓ Detection Rate (≥90%%): PASS (%.1f%%)\n', detection_rate * 100);
else
    fprintf('❌ Detection Rate (≥90%%): FAIL (%.1f%%)\n', detection_rate * 100);
end

if criteria.velocity_accuracy
    fprintf('✓ Velocity Accuracy (≥90%% within ±3 mph): PASS (%.1f%%)\n', velocity_accuracy * 100);
else
    fprintf('❌ Velocity Accuracy (≥90%% within ±3 mph): FAIL (%.1f%%)\n', velocity_accuracy * 100);
end

if m2_validation.m2_success
    fprintf('\n🎯 MILESTONE 2 STATUS: ✅ PASSED\n');
    fprintf('All acceptance criteria have been met!\n');
else
    fprintf('\n⚠️  MILESTONE 2 STATUS: ❌ NOT PASSED\n');
    fprintf('Some acceptance criteria require improvement.\n');
end

%% Step 7: Category-Specific Analysis
if isfield(m2_validation.detailed, 'category_analysis')
    fprintf('\n=== PERFORMANCE BY SHOT TYPE ===\n');
    cat_analysis = m2_validation.detailed.category_analysis;
    categories = fieldnames(cat_analysis);
    
    for i = 1:length(categories)
        cat = categories{i};
        if cat_analysis.(cat).count > 0
            fprintf('%s shots (%d total):\n', upper(cat), cat_analysis.(cat).count);
            fprintf('  Ball RMS: %.2f mph, Club RMS: %.2f mph\n', ...
                cat_analysis.(cat).ball_rms_error, cat_analysis.(cat).club_rms_error);
            fprintf('  Accuracy (±3 mph): Ball %.1f%%, Club %.1f%%\n', ...
                cat_analysis.(cat).ball_accuracy_3mph * 100, cat_analysis.(cat).club_accuracy_3mph * 100);
        end
    end
end

%% Step 8: Save Results
fprintf('\nStep 8: Saving Milestone 2 results...\n');

% Save comprehensive results
save('milestone2_results.mat', 'm2_results', 'enhanced_results', 'm2_validation', 'config');

% Save lightweight summary for quick reference
m2_summary = struct();
m2_summary.detection_rate = detection_rate;
m2_summary.velocity_accuracy = velocity_accuracy;
m2_summary.ball_rms_error = ball_rms;
m2_summary.club_rms_error = club_rms;
m2_summary.m2_success = m2_validation.m2_success;
m2_summary.processing_date = datetime('now');

save('milestone2_summary.mat', 'm2_summary');

fprintf('✓ Results saved to milestone2_results.mat\n');

%% Step 9: Debug track structures
fprintf('Step 9: Debugging track structures...\n');
debug_track_structures(enhanced_results);

%% Step 10: Generate Visualizations (Optional)
fprintf('\nStep 10: Creating Milestone 2 visualizations...\n');

try
    % Create M2-specific plots
    create_milestone2_plots(enhanced_results, m2_validation, config);
    fprintf('✓ Visualizations saved to milestone2_figures/\n');
catch ME
    fprintf('⚠️ Visualization generation failed: %s\n', ME.message);
end

fprintf('\n=== MILESTONE 2 ANALYSIS COMPLETE ===\n');
fprintf('Next steps:\n');
if m2_validation.m2_success
    fprintf('- Proceed to Milestone 3: Club/Ball Separation Logic\n');
    fprintf('- Focus on achieving 80%% wedge accuracy, 90%% other clubs\n');
else
    fprintf('- Debug and improve velocity accuracy\n');
    fprintf('- Adjust CFAR and tracking parameters\n');
    fprintf('- Re-run analysis after improvements\n');
end

%% Debugging Functions
function debug_track_structures(enhanced_results)
%% Debug function to analyze track structures and identify data format issues
% This function helps diagnose why velocity extraction is failing

fprintf('\n=== TRACK STRUCTURE DEBUGGING ===\n');

for shot_idx = 1:length(enhanced_results)
    fprintf('\nShot %d: %s\n', shot_idx, enhanced_results{shot_idx}.shot_name);
    
    if isfield(enhanced_results{shot_idx}, 'smoothed_tracks')
        tracks = enhanced_results{shot_idx}.smoothed_tracks;
        fprintf('  Found %d smoothed tracks\n', length(tracks));
        
        if length(tracks) > 0
            % Analyze first track structure
            track1 = tracks(1);
            fprintf('  Track 1 fields: %s\n', strjoin(fieldnames(track1), ', '));
            
            % Check velocity-related fields
            velocity_fields = {'mean_velocity', 'smoothed_velocity', 'velocity_data', 'detections'};
            for field = velocity_fields
                if isfield(track1, field{1})
                    field_data = track1.(field{1});
                    fprintf('  %s: type=%s, size=%s, empty=%s\n', ...
                        field{1}, class(field_data), mat2str(size(field_data)), mat2str(isempty(field_data)));
                    
                    % If it's a struct, show its fields
                    if isstruct(field_data) && ~isempty(field_data)
                        fprintf('    Sub-fields: %s\n', strjoin(fieldnames(field_data), ', '));
                    end
                else
                    fprintf('  %s: NOT FOUND\n', field{1});
                end
            end
            
            % Check quality-related fields
            quality_fields = {'quality_score', 'num_detections', 'duration', 'max_magnitude'};
            for field = quality_fields
                if isfield(track1, field{1})
                    field_data = track1.(field{1});
                    fprintf('  %s: type=%s, value=%s\n', ...
                        field{1}, class(field_data), mat2str(field_data));
                else
                    fprintf('  %s: NOT FOUND\n', field{1});
                end
            end
        end
    else
        fprintf('  No smoothed_tracks field found\n');
    end
end

fprintf('\n=== END TRACK STRUCTURE DEBUGGING ===\n');
end

function create_milestone2_plots(enhanced_results, validation, config)
%% Create Milestone 2 specific visualization plots

% Create output directory
output_dir = 'milestone2_figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Plot 1: M2 Performance Summary
fig1 = figure('Position', [100, 100, 1200, 800]);

subplot(2, 3, 1);
bar([validation.detection_rate, validation.velocity_accuracy] * 100);
title('M2 Key Metrics');
ylabel('Percentage (%)');
xticklabels({'Detection Rate', 'Velocity Accuracy'});
ylim([0, 100]);
yline(90, 'r--', '90% Target');

subplot(2, 3, 2);
bar([validation.ball_rms_error, validation.club_rms_error]);
title('RMS Errors');
ylabel('Error (mph)');
xticklabels({'Ball Speed', 'Club Speed'});
yline(3, 'r--', '3 mph Target');

subplot(2, 3, 3);
accuracy_data = [validation.ball_accuracy_3mph, validation.club_accuracy_3mph, ...
                validation.ball_accuracy_1mph, validation.club_accuracy_1mph] * 100;
bar(accuracy_data);
title('Accuracy Breakdown');
ylabel('Accuracy (%)');
xticklabels({'Ball ±3mph', 'Club ±3mph', 'Ball ±1mph', 'Club ±1mph'});
ylim([0, 100]);

% Add more subplots as needed...

sgtitle('Milestone 2: Enhanced DSP Performance Analysis');
saveas(fig1, fullfile(output_dir, 'm2_performance_summary.png'));

end
