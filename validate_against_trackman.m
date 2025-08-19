function validation_results = validate_against_trackman(stft_results, config)
%% Validate Against TrackMan Data
% Compares algorithm outputs with TrackMan reference measurements
%
% Inputs:
%   stft_results - Cell array of STFT processing results
%   config - Configuration structure
%
% Outputs:
%   validation_results - Structure containing validation metrics

fprintf('Performing TrackMan validation analysis...\n');

num_shots = length(stft_results);
timing_alignments = zeros(num_shots, 1);
velocity_errors = cell(num_shots, 1);
correlation_scores = zeros(num_shots, 1);

% Initialize result arrays
ball_speed_estimated = zeros(num_shots, 1);
club_speed_estimated = zeros(num_shots, 1);
ball_speed_trackman = zeros(num_shots, 1);
club_speed_trackman = zeros(num_shots, 1);
shot_types = cell(num_shots, 1);

for i = 1:num_shots
    try
        shot_info = stft_results{i}.shot_info;
        stft_data = stft_results{i}.stft_data;
        
        % Get TrackMan reference values
        ball_speed_ref = shot_info.ball_speed_trackman;
        club_speed_ref = shot_info.club_speed_trackman;
        
        % Store reference values
        ball_speed_trackman(i) = ball_speed_ref;
        club_speed_trackman(i) = club_speed_ref;
        shot_types{i} = shot_info.shot_type;
        
        % Estimate speeds from radar data
        [ball_est, club_est, timing_quality] = estimate_speeds_from_stft(stft_data, config);
        
        ball_speed_estimated(i) = ball_est;
        club_speed_estimated(i) = club_est;
        timing_alignments(i) = timing_quality;
        
        % Calculate velocity errors
        if ~isnan(ball_speed_ref) && ~isnan(club_speed_ref)
            ball_error = abs(ball_est - ball_speed_ref);
            club_error = abs(club_est - club_speed_ref);
            velocity_errors{i} = [ball_error, club_error];
            
            % Calculate correlation (simplified metric)
            correlation_scores(i) = calculate_correlation_metric(ball_est, club_est, ball_speed_ref, club_speed_ref);
        else
            velocity_errors{i} = [NaN, NaN];
            correlation_scores(i) = NaN;
        end
        
    catch ME
        fprintf('Warning: Error processing shot %d: %s\n', i, ME.message);
        velocity_errors{i} = [NaN, NaN];
        correlation_scores(i) = NaN;
        timing_alignments(i) = 0;
    end
end

% Calculate overall validation metrics
valid_shots = ~isnan(ball_speed_trackman) & ~isnan(club_speed_trackman) & ...
              ~isnan(ball_speed_estimated) & ~isnan(club_speed_estimated);

if sum(valid_shots) > 0
    % Timing alignment percentage (≥90% target)
    timing_alignment_rate = mean(timing_alignments >= 0.7); % 70% threshold for timing quality
    
    % Velocity accuracy (±3 mph target for Milestone 1)
    ball_errors = abs(ball_speed_estimated(valid_shots) - ball_speed_trackman(valid_shots));
    club_errors = abs(club_speed_estimated(valid_shots) - club_speed_trackman(valid_shots));
    
    ball_accuracy_3mph = mean(ball_errors <= 3);
    club_accuracy_3mph = mean(club_errors <= 3);
    overall_accuracy_3mph = mean([ball_errors; club_errors] <= 3);
    
    % Calculate RMS errors
    ball_rms_error = sqrt(mean(ball_errors.^2));
    club_rms_error = sqrt(mean(club_errors.^2));
    
    % Correlation analysis
    valid_correlations = correlation_scores(~isnan(correlation_scores));
    mean_correlation = mean(valid_correlations);
    
    % Speed range analysis
    ball_speed_range_est = [min(ball_speed_estimated(valid_shots)), max(ball_speed_estimated(valid_shots))];
    ball_speed_range_ref = [min(ball_speed_trackman(valid_shots)), max(ball_speed_trackman(valid_shots))];
    
    % Category-specific analysis
    category_analysis = analyze_by_shot_type(ball_speed_estimated, club_speed_estimated, ...
                                            ball_speed_trackman, club_speed_trackman, shot_types);
    
else
    % No valid shots for comparison
    timing_alignment_rate = 0;
    ball_accuracy_3mph = 0;
    club_accuracy_3mph = 0;
    overall_accuracy_3mph = 0;
    ball_rms_error = inf;
    club_rms_error = inf;
    mean_correlation = 0;
    ball_speed_range_est = [NaN, NaN];
    ball_speed_range_ref = [NaN, NaN];
    category_analysis = struct();
end

% Compile validation results
validation_results = struct();
validation_results.timing_alignment = timing_alignment_rate;
validation_results.velocity_accuracy = overall_accuracy_3mph;
validation_results.ball_accuracy_3mph = ball_accuracy_3mph;
validation_results.club_accuracy_3mph = club_accuracy_3mph;
validation_results.ball_rms_error = ball_rms_error;
validation_results.club_rms_error = club_rms_error;
validation_results.correlation = mean_correlation;
validation_results.num_valid_shots = sum(valid_shots);
validation_results.num_total_shots = num_shots;

% Store detailed results for further analysis
validation_results.detailed = struct();
validation_results.detailed.ball_speed_estimated = ball_speed_estimated;
validation_results.detailed.club_speed_estimated = club_speed_estimated;
validation_results.detailed.ball_speed_trackman = ball_speed_trackman;
validation_results.detailed.club_speed_trackman = club_speed_trackman;
validation_results.detailed.velocity_errors = velocity_errors;
validation_results.detailed.timing_alignments = timing_alignments;
validation_results.detailed.correlation_scores = correlation_scores;
validation_results.detailed.shot_types = shot_types;
validation_results.detailed.ball_speed_range_est = ball_speed_range_est;
validation_results.detailed.ball_speed_range_ref = ball_speed_range_ref;
validation_results.detailed.category_analysis = category_analysis;

fprintf('TrackMan validation complete:\n');
fprintf('- Valid shots analyzed: %d/%d\n', sum(valid_shots), num_shots);
fprintf('- Timing alignment rate: %.1f%%\n', timing_alignment_rate * 100);
fprintf('- Overall velocity accuracy (±3 mph): %.1f%%\n', overall_accuracy_3mph * 100);
fprintf('- Ball speed RMS error: %.2f mph\n', ball_rms_error);
fprintf('- Club speed RMS error: %.2f mph\n', club_rms_error);

end

function [ball_speed, club_speed, timing_quality] = estimate_speeds_from_stft(stft_data, config)
%% Estimate Ball and Club Speeds from STFT Data
% Simplified speed estimation for Milestone 1 validation

% Extract tracked objects
tracked_objects = stft_data.tracked_objects;

if isempty(tracked_objects)
    ball_speed = NaN;
    club_speed = NaN;
    timing_quality = 0;
    return;
end

% Sort tracks by magnitude (strongest signals first)
track_magnitudes = [tracked_objects.max_magnitude];
[~, sort_idx] = sort(track_magnitudes, 'descend');
sorted_tracks = tracked_objects(sort_idx);

% Simple heuristic: highest velocity is likely the ball
track_velocities = abs([sorted_tracks.mean_velocity]);
[max_velocity, max_idx] = max(track_velocities);

if length(sorted_tracks) >= 2
    % Two strongest tracks: assume higher velocity is ball
    if max_idx == 1
        ball_speed = abs(sorted_tracks(1).mean_velocity);
        club_speed = abs(sorted_tracks(2).mean_velocity);
    else
        ball_speed = abs(sorted_tracks(max_idx).mean_velocity);
        club_speed = abs(sorted_tracks(1).mean_velocity);
    end
else
    % Only one track detected
    ball_speed = max_velocity;
    club_speed = NaN;
end

% Assess timing quality based on track characteristics
timing_quality = assess_timing_quality(sorted_tracks);

end

function correlation = calculate_correlation_metric(ball_est, club_est, ball_ref, club_ref)
%% Calculate Simple Correlation Metric

% Normalized distance metric
ball_error_norm = abs(ball_est - ball_ref) / ball_ref;
club_error_norm = abs(club_est - club_ref) / club_ref;

% Combined error (lower is better)
combined_error = sqrt(ball_error_norm^2 + club_error_norm^2);

% Convert to correlation-like metric (higher is better)
correlation = exp(-combined_error);

end

function timing_quality = assess_timing_quality(tracks)
%% Assess Timing Quality of Track Detection

if isempty(tracks)
    timing_quality = 0;
    return;
end

% Factors that indicate good timing:
% 1. Track duration (longer tracks are more reliable)
% 2. Number of detections per track
% 3. Velocity consistency within tracks

quality_scores = zeros(length(tracks), 1);

for i = 1:length(tracks)
    track = tracks(i);
    
    % Duration score (0-1)
    duration_score = min(track.duration / 0.5, 1); % Normalize to 0.5 seconds
    
    % Detection density score (0-1)
    detection_score = min(track.num_detections / 10, 1); % Normalize to 10 detections
    
    % Consistency score (0-1) - lower std is better
    if track.velocity_std > 0
        consistency_score = exp(-track.velocity_std / 10); % Normalize to 10 mph std
    else
        consistency_score = 1;
    end
    
    % Combined quality score
    quality_scores(i) = (duration_score + detection_score + consistency_score) / 3;
end

% Overall timing quality is the mean of individual track qualities
timing_quality = mean(quality_scores);

end

function category_analysis = analyze_by_shot_type(ball_est, club_est, ball_ref, club_ref, shot_types)
%% Analyze Performance by Shot Category

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
    else
        category_analysis.(cat) = struct();
        category_analysis.(cat).count = 0;
        category_analysis.(cat).ball_rms_error = NaN;
        category_analysis.(cat).club_rms_error = NaN;
        category_analysis.(cat).ball_accuracy_3mph = NaN;
        category_analysis.(cat).club_accuracy_3mph = NaN;
    end
end

end
