function validation_results = milestone2_validation(enhanced_results, config)
%% Milestone 2: Enhanced Validation Framework
% Validates enhanced algorithms against TrackMan with M2 criteria
%
% Inputs:
%   enhanced_results - Results from M2 processing pipeline
%   config - Configuration structure
%
% Outputs:
%   validation_results - M2-specific validation metrics

fprintf('Performing Milestone 2 validation analysis...\n');

% M2 acceptance criteria
m2_criteria = struct();
m2_criteria.detection_rate_target = 0.90; % ≥90% detection rate
m2_criteria.velocity_accuracy_target = 3; % ±3 mph
m2_criteria.cfar_fa_rate_target = 0.05; % <5% false alarms

num_shots = length(enhanced_results);

% Initialize result arrays
ball_speed_estimated = zeros(num_shots, 1);
club_speed_estimated = zeros(num_shots, 1);
ball_speed_trackman = zeros(num_shots, 1);
club_speed_trackman = zeros(num_shots, 1);
shot_types = cell(num_shots, 1);
detection_success = zeros(num_shots, 1);

% Process each shot
for i = 1:num_shots
    try
        shot_info = enhanced_results{i}.shot_info;
        advanced_tracks = enhanced_results{i}.advanced_tracks;
        
        % Get TrackMan reference values
        ball_speed_ref = shot_info.ball_speed_trackman;
        club_speed_ref = shot_info.club_speed_trackman;
        
        % Store reference values
        ball_speed_trackman(i) = ball_speed_ref;
        club_speed_trackman(i) = club_speed_ref;
        shot_types{i} = shot_info.shot_type;
        
        % Check if we have detections
        if ~isempty(advanced_tracks)
            detection_success(i) = 1;
            
            % Estimate speeds from enhanced tracks
            [ball_est, club_est] = estimate_speeds_from_enhanced_tracks(advanced_tracks, config);
            
            ball_speed_estimated(i) = ball_est;
            club_speed_estimated(i) = club_est;
        else
            detection_success(i) = 0;
            ball_speed_estimated(i) = NaN;
            club_speed_estimated(i) = NaN;
        end
        
    catch ME
        fprintf('Warning: Error processing shot %d: %s\n', i, ME.message);
        fprintf('Error occurred at line: %d in function: %s\n', ME.stack(1).line, ME.stack(1).name);
        if length(ME.stack) > 1
            fprintf('Called from line: %d in function: %s\n', ME.stack(2).line, ME.stack(2).name);
        end
        detection_success(i) = 0;
        ball_speed_estimated(i) = NaN;
        club_speed_estimated(i) = NaN;
    end
end

% Calculate validation metrics
valid_shots = ~isnan(ball_speed_trackman) & ~isnan(club_speed_trackman) & ...
              ~isnan(ball_speed_estimated) & ~isnan(club_speed_estimated);

% Detection rate
detection_rate = mean(detection_success);

if sum(valid_shots) > 0
    % Velocity accuracy metrics
    ball_errors = abs(ball_speed_estimated(valid_shots) - ball_speed_trackman(valid_shots));
    club_errors = abs(club_speed_estimated(valid_shots) - club_speed_trackman(valid_shots));
    
    % Accuracy within ±3 mph
    ball_accuracy_3mph = mean(ball_errors <= 3);
    club_accuracy_3mph = mean(club_errors <= 3);
    overall_accuracy_3mph = mean([ball_errors; club_errors] <= 3);
    
    % RMS errors
    ball_rms_error = sqrt(mean(ball_errors.^2));
    club_rms_error = sqrt(mean(club_errors.^2));
    
    % Accuracy within ±1 mph (ultimate goal)
    ball_accuracy_1mph = mean(ball_errors <= 1);
    club_accuracy_1mph = mean(club_errors <= 1);
    overall_accuracy_1mph = mean([ball_errors; club_errors] <= 1);
    
    % Category-specific analysis
    category_analysis = analyze_by_shot_type_m2(ball_speed_estimated, club_speed_estimated, ...
                                               ball_speed_trackman, club_speed_trackman, shot_types);
    
else
    % No valid shots for comparison
    ball_accuracy_3mph = 0;
    club_accuracy_3mph = 0;
    overall_accuracy_3mph = 0;
    ball_rms_error = inf;
    club_rms_error = inf;
    ball_accuracy_1mph = 0;
    club_accuracy_1mph = 0;
    overall_accuracy_1mph = 0;
    category_analysis = struct();
end

% Check M2 acceptance criteria
criteria_met = struct();
criteria_met.detection_rate = detection_rate >= m2_criteria.detection_rate_target;
criteria_met.velocity_accuracy = overall_accuracy_3mph >= 0.90; % 90% within ±3 mph
criteria_met.ball_accuracy = ball_accuracy_3mph >= 0.90;
criteria_met.club_accuracy = club_accuracy_3mph >= 0.90;

% Overall M2 success
m2_success = criteria_met.detection_rate && criteria_met.velocity_accuracy;

% Compile validation results
validation_results = struct();
validation_results.detection_rate = detection_rate;
validation_results.velocity_accuracy = overall_accuracy_3mph;
validation_results.ball_accuracy_3mph = ball_accuracy_3mph;
validation_results.club_accuracy_3mph = club_accuracy_3mph;
validation_results.ball_rms_error = ball_rms_error;
validation_results.club_rms_error = club_rms_error;

% Additional metrics for M2
validation_results.ball_accuracy_1mph = ball_accuracy_1mph;
validation_results.club_accuracy_1mph = club_accuracy_1mph;
validation_results.overall_accuracy_1mph = overall_accuracy_1mph;

validation_results.num_valid_shots = sum(valid_shots);
validation_results.num_total_shots = num_shots;
validation_results.criteria_met = criteria_met;
validation_results.m2_success = m2_success;

% Store detailed results
validation_results.detailed = struct();
validation_results.detailed.ball_speed_estimated = ball_speed_estimated;
validation_results.detailed.club_speed_estimated = club_speed_estimated;
validation_results.detailed.ball_speed_trackman = ball_speed_trackman;
validation_results.detailed.club_speed_trackman = club_speed_trackman;
validation_results.detailed.shot_types = shot_types;
validation_results.detailed.detection_success = detection_success;
validation_results.detailed.category_analysis = category_analysis;

% Performance comparison with M1 (if available)
if isfield(enhanced_results{1}, 'm1_result')
    validation_results.m1_comparison = compare_with_milestone1(enhanced_results, validation_results);
end

fprintf('Milestone 2 validation complete:\n');
fprintf('- Detection rate: %.1f%% (Target: ≥90%%)\n', detection_rate * 100);
fprintf('- Valid shots analyzed: %d/%d\n', sum(valid_shots), num_shots);
fprintf('- Overall velocity accuracy (±3 mph): %.1f%% (Target: ≥90%%)\n', overall_accuracy_3mph * 100);
fprintf('- Ball speed RMS error: %.2f mph\n', ball_rms_error);
fprintf('- Club speed RMS error: %.2f mph\n', club_rms_error);

if m2_success
    fprintf('🎯 MILESTONE 2 ACCEPTANCE CRITERIA: PASSED\n');
else
    fprintf('⚠️  MILESTONE 2 ACCEPTANCE CRITERIA: NOT MET\n');
    if ~criteria_met.detection_rate
        fprintf('   - Detection rate below 90%%\n');
    end
    if ~criteria_met.velocity_accuracy
        fprintf('   - Velocity accuracy below 90%% at ±3 mph\n');
    end
end

end

function [ball_speed, club_speed] = estimate_speeds_from_enhanced_tracks(tracks, config)
%% Estimate Ball and Club Speeds from Enhanced Tracks - FIXED VERSION
% Fixed version with proper struct handling and error checking
%
% Inputs:
%   tracks - Enhanced tracked objects from M2 pipeline
%   config - Configuration structure
%
% Outputs:
%   ball_speed - Estimated ball speed (mph)
%   club_speed - Estimated club head speed (mph)

% Input validation
if isempty(tracks) || ~isstruct(tracks)
    ball_speed = NaN;
    club_speed = NaN;
    return;
end

% Debug: Print track structure to understand format
fprintf('   Debug: Processing %d tracks\n', length(tracks));
if length(tracks) > 0
    fprintf('   Debug: Track fields: %s\n', strjoin(fieldnames(tracks(1)), ', '));
end

% Filter tracks by quality and golf-relevant velocities
valid_tracks = [];
track_info = []; % Store velocity and quality info

for i = 1:length(tracks)
    track = tracks(i);
    
    % Extract mean velocity with multiple fallback methods
    mean_vel = extract_track_velocity(track);
    
    % Extract quality score with fallback
    quality = extract_track_quality(track);
    
    % Debug output for first few tracks
    if i <= 3
        fprintf('   Debug: Track %d - Velocity: %.1f, Quality: %.2f\n', i, mean_vel, quality);
    end
    
    % Golf-relevant velocity range and quality threshold
    if ~isnan(mean_vel) && ~isnan(quality) && ...
       abs(mean_vel) >= 15 && abs(mean_vel) <= 200 && quality >= 0.3
        
        track_info(end+1, :) = [abs(mean_vel), quality, i]; % [velocity, quality, index]
        valid_tracks(end+1) = track;
    end
end

fprintf('   Debug: Found %d valid tracks from %d total\n', length(valid_tracks), length(tracks));

% Speed estimation logic
if size(track_info, 1) >= 2
    % Sort by combined score (velocity + quality)
    velocities = track_info(:, 1);
    qualities = track_info(:, 2);
    
    % Normalize and combine scores
    max_velocity = max(velocities);
    if max_velocity > 0
        velocity_scores = velocities / max_velocity;
    else
        velocity_scores = ones(size(velocities));
    end
    combined_scores = 0.6 * velocity_scores + 0.4 * qualities;
    
    [~, sort_idx] = sort(combined_scores, 'descend');
    sorted_velocities = velocities(sort_idx);
    
    % Select top 2 candidates
    candidate1 = sorted_velocities(1);
    candidate2 = sorted_velocities(2);
    
    % Apply golf physics: ball speed ≥ club speed
    ball_speed = max(candidate1, candidate2);
    club_speed = min(candidate1, candidate2);
    
    % Sanity check: ensure reasonable separation
    if (ball_speed - club_speed) < 5 % Less than 5 mph difference
        % For very close velocities, apply small offset based on golf physics
        avg_speed = (ball_speed + club_speed) / 2;
        ball_speed = avg_speed + 2.5; % Ball slightly faster
        club_speed = avg_speed - 2.5; % Club slightly slower
    end
    
elseif size(track_info, 1) == 1
    % Single track - assume it's the ball (most common detection)
    ball_speed = track_info(1, 1);
    club_speed = NaN;
else
    % No valid tracks
    ball_speed = NaN;
    club_speed = NaN;
end

fprintf('   Debug: Final speeds - Ball: %.1f mph, Club: %.1f mph\n', ball_speed, club_speed);

end

function velocity = extract_track_velocity(track)
%% Extract velocity from track structure with multiple fallback methods
% Handles different possible data structures

velocity = NaN;

try
    % Method 1: Check for mean_velocity field
    if isfield(track, 'mean_velocity')
        if isnumeric(track.mean_velocity) && ~isempty(track.mean_velocity)
            velocity = track.mean_velocity;
            return;
        end
    end
    
    % Method 2: Calculate from detections array
    if isfield(track, 'detections') && ~isempty(track.detections)
        if isstruct(track.detections) && isfield(track.detections, 'velocity')
            detection_velocities = [track.detections.velocity];
            if ~isempty(detection_velocities) && isnumeric(detection_velocities)
                velocity = mean(detection_velocities);
                return;
            end
        end
    end
    
    % Method 3: Check if velocity data is in a nested structure
    if isfield(track, 'velocity_data')
        if isstruct(track.velocity_data) && isfield(track.velocity_data, 'mean')
            velocity = track.velocity_data.mean;
            return;
        end
    end
    
    % Method 4: Check for smoothed velocity
    if isfield(track, 'smoothed_velocity')
        if isnumeric(track.smoothed_velocity) && ~isempty(track.smoothed_velocity)
            velocity = track.smoothed_velocity;
            return;
        end
    end
    
catch ME
    % If any method fails, continue to next method
    fprintf('   Warning: Velocity extraction method failed: %s\n', ME.message);
end

% If all methods fail, velocity remains NaN
if isnan(velocity)
    fprintf('   Warning: Could not extract velocity from track structure\n');
end

end

function quality = extract_track_quality(track)
%% Extract quality score from track structure with fallback methods

quality = 0.5; % Default quality if not found

try
    % Method 1: Direct quality_score field
    if isfield(track, 'quality_score')
        if isnumeric(track.quality_score) && ~isempty(track.quality_score)
            quality = track.quality_score;
            return;
        end
    end
    
    % Method 2: Calculate basic quality from track statistics
    if isfield(track, 'num_detections') && isfield(track, 'duration')
        if isnumeric(track.num_detections) && isnumeric(track.duration)
            % Simple quality based on track length and duration
            length_score = min(track.num_detections / 5, 1);
            duration_score = min(track.duration / 0.2, 1);
            quality = (length_score + duration_score) / 2;
            return;
        end
    end
    
    % Method 3: Use magnitude-based quality
    if isfield(track, 'max_magnitude')
        if isnumeric(track.max_magnitude) && ~isempty(track.max_magnitude)
            % Normalize magnitude to quality score (rough approximation)
            quality = max(0, min(1, (track.max_magnitude - 80) / 40));
            return;
        end
    end
    
catch ME
    % If extraction fails, use default quality
    fprintf('   Warning: Quality extraction failed: %s\n', ME.message);
end

end

function category_analysis = analyze_by_shot_type_m2(ball_est, club_est, ball_ref, club_ref, shot_types)
%% Enhanced category analysis for M2

categories = {'wedge', 'iron', 'driver'};
category_analysis = struct();

for i = 1:length(categories)
    cat = categories{i};
    cat_mask = strcmp(shot_types, cat) & ~isnan(ball_ref) & ~isnan(club_ref) & ...
               ~isnan(ball_est) & ~isnan(club_est);
    
    if sum(cat_mask) > 0
        ball_errors = abs(ball_est(cat_mask) - ball_ref(cat_mask));
        club_errors = abs(club_est(cat_mask) - club_ref(cat_mask));
        
        category_analysis.(cat) = struct();
        category_analysis.(cat).count = sum(cat_mask);
        category_analysis.(cat).ball_rms_error = sqrt(mean(ball_errors.^2));
        category_analysis.(cat).club_rms_error = sqrt(mean(club_errors.^2));
        category_analysis.(cat).ball_accuracy_3mph = mean(ball_errors <= 3);
        category_analysis.(cat).club_accuracy_3mph = mean(club_errors <= 3);
        category_analysis.(cat).ball_accuracy_1mph = mean(ball_errors <= 1);
        category_analysis.(cat).club_accuracy_1mph = mean(club_errors <= 1);
    else
        category_analysis.(cat) = struct();
        category_analysis.(cat).count = 0;
        category_analysis.(cat).ball_rms_error = NaN;
        category_analysis.(cat).club_rms_error = NaN;
        category_analysis.(cat).ball_accuracy_3mph = NaN;
        category_analysis.(cat).club_accuracy_3mph = NaN;
        category_analysis.(cat).ball_accuracy_1mph = NaN;
        category_analysis.(cat).club_accuracy_1mph = NaN;
    end
end

end

function m1_comparison = compare_with_milestone1(enhanced_results, m2_validation)
%% Compare M2 performance with M1 baseline

% This would compare the M2 results with the original M1 results
% For now, return basic structure - can be enhanced later

m1_comparison = struct();
m1_comparison.improvement_available = true;
m1_comparison.note = 'M1 comparison requires baseline results';

end
