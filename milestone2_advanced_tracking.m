function [advanced_tracks, tracking_stats] = milestone2_advanced_tracking(cfar_peaks, config)
%% Milestone 2: Advanced Multi-Target Tracking
% Improved tracking with tighter tolerances and quality scoring
%
% Inputs:
%   cfar_peaks - Peak detections from CFAR detector
%   config - Configuration with tracking parameters
%
% Outputs:
%   advanced_tracks - Enhanced tracked objects with quality scores
%   tracking_stats - Tracking performance metrics

fprintf('   M2: Applying advanced multi-target tracking...\n');

% Enhanced tracking parameters (tighter than M1)
if ~isfield(config, 'm2_velocity_tolerance')
    config.m2_velocity_tolerance = 5; % mph (reduced from 10)
end
if ~isfield(config, 'm2_min_track_length')
    config.m2_min_track_length = 3; % detections (increased from 2)
end
if ~isfield(config, 'm2_max_time_gap')
    config.m2_max_time_gap = 0.05; % seconds (reduced from 0.2)
end
if ~isfield(config, 'm2_quality_threshold')
    config.m2_quality_threshold = 0.4; % minimum track quality (0-1)
end

% Extract peak data
if isfield(cfar_peaks, 'velocities') && ~isempty(cfar_peaks.velocities)
    velocities = cfar_peaks.velocities;
    times = cfar_peaks.times;
    magnitudes = cfar_peaks.magnitudes;
else
    % No detections
    advanced_tracks = [];
    tracking_stats = struct('total_tracks', 0, 'valid_tracks', 0, 'avg_quality', 0);
    fprintf('   M2 Tracking: No CFAR detections to track\n');
    return;
end

% Tracking parameters
velocity_tolerance = config.m2_velocity_tolerance;
min_track_length = config.m2_min_track_length;
max_time_gap = config.m2_max_time_gap;
quality_threshold = config.m2_quality_threshold;

% Limit detections for computational efficiency
max_detections = 1000;
if length(times) > max_detections
    [~, strong_idx] = sort(magnitudes, 'descend');
    keep_idx = strong_idx(1:max_detections);
    times = times(keep_idx);
    velocities = velocities(keep_idx);
    magnitudes = magnitudes(keep_idx);
end

% Sort detections by time
[sorted_times, sort_idx] = sort(times);
sorted_velocities = velocities(sort_idx);
sorted_magnitudes = magnitudes(sort_idx);

% Initialize tracking
tracks = struct([]);
track_id = 0;

% Advanced nearest-neighbor tracking with gating
for i = 1:length(sorted_times)
    current_time = sorted_times(i);
    current_velocity = sorted_velocities(i);
    current_magnitude = sorted_magnitudes(i);
    
    % Find compatible existing tracks
    compatible_tracks = [];
    for j = 1:length(tracks)
        if ~isempty(tracks(j).detections)
            last_detection = tracks(j).detections(end);
            time_gap = current_time - last_detection.time;
            velocity_diff = abs(current_velocity - last_detection.velocity);
            
            % Velocity-based gating with time consideration
            if time_gap <= max_time_gap && velocity_diff <= velocity_tolerance
                % Additional quality-based gating
                predicted_velocity = predict_velocity(tracks(j), current_time);
                prediction_error = abs(current_velocity - predicted_velocity);
                
                if prediction_error <= velocity_tolerance * 1.5 % Allow 50% more tolerance for prediction
                    compatible_tracks(end+1) = j;
                end
            end
        end
    end
    
    if isempty(compatible_tracks)
        % Create new track
        track_id = track_id + 1;
        new_track = struct();
        new_track.id = track_id;
        new_track.detections = struct('time', current_time, ...
                                     'velocity', current_velocity, ...
                                     'magnitude', current_magnitude);
        
        if isempty(tracks)
            tracks = new_track;
        else
            tracks(end+1) = new_track;
        end
    else
        % Find best matching track using combined distance metric
        best_track_idx = find_best_track_match(compatible_tracks, tracks, ...
                                               current_time, current_velocity, current_magnitude);
        
        % Add detection to best track
        new_detection = struct('time', current_time, ...
                              'velocity', current_velocity, ...
                              'magnitude', current_magnitude);
        tracks(best_track_idx).detections(end+1) = new_detection;
    end
end

% Filter tracks by length and calculate quality scores
valid_tracks = [];
quality_scores = [];

for i = 1:length(tracks)
    if length(tracks(i).detections) >= min_track_length
        track = tracks(i);
        
        % Calculate comprehensive track statistics
        track_times = [track.detections.time];
        track_velocities = [track.detections.velocity];
        track_magnitudes = [track.detections.magnitude];
        
        % Basic statistics
        track.duration = max(track_times) - min(track_times);
        track.mean_velocity = mean(track_velocities);
        track.velocity_std = std(track_velocities);
        track.max_magnitude = max(track_magnitudes);
        track.mean_magnitude = mean(track_magnitudes);
        track.num_detections = length(track.detections);
        
        % Calculate track quality score (0-1)
        track.quality_score = calculate_track_quality(track, config);
        quality_scores(end+1) = track.quality_score;
        
        % Apply quality threshold
        if track.quality_score >= quality_threshold
            if isempty(valid_tracks)
                valid_tracks = track;
            else
                valid_tracks(end+1) = track;
            end
        end
    end
end

advanced_tracks = valid_tracks;

% Calculate tracking statistics
tracking_stats = struct();
tracking_stats.total_tracks_created = track_id;
tracking_stats.tracks_meeting_length = sum(arrayfun(@(t) length(t.detections) >= min_track_length, tracks));
tracking_stats.tracks_meeting_quality = length(valid_tracks);
tracking_stats.tracking_efficiency = length(valid_tracks) / max(track_id, 1);

if ~isempty(quality_scores)
    tracking_stats.avg_quality_score = mean(quality_scores);
    tracking_stats.quality_std = std(quality_scores);
else
    tracking_stats.avg_quality_score = 0;
    tracking_stats.quality_std = 0;
end

if ~isempty(valid_tracks)
    tracking_stats.avg_track_length = mean([valid_tracks.num_detections]);
    tracking_stats.avg_track_duration = mean([valid_tracks.duration]);
    
    % Safely extract velocity range
    velocities = [];
    for i = 1:length(valid_tracks)
        if isfield(valid_tracks(i), 'mean_velocity') && isnumeric(valid_tracks(i).mean_velocity)
            velocities(end+1) = valid_tracks(i).mean_velocity;
        end
    end
    
    if ~isempty(velocities)
        tracking_stats.velocity_range = [min(velocities), max(velocities)];
    else
        tracking_stats.velocity_range = [NaN, NaN];
    end
else
    tracking_stats.avg_track_length = 0;
    tracking_stats.avg_track_duration = 0;
    tracking_stats.velocity_range = [NaN, NaN];
end

fprintf('   M2 Tracking: %d valid tracks (quality≥%.2f) from %d total\n', ...
    length(valid_tracks), quality_threshold, track_id);

end

function predicted_velocity = predict_velocity(track, target_time)
%% Predict velocity at target time using linear extrapolation
% Simple constant velocity model

detections = track.detections;
times = [detections.time];
velocities = [detections.velocity];

if length(times) >= 2
    % Linear extrapolation using last two points
    dt = times(end) - times(end-1);
    dv = velocities(end) - velocities(end-1);
    acceleration = dv / dt;
    
    time_extrapolation = target_time - times(end);
    predicted_velocity = velocities(end) + acceleration * time_extrapolation;
else
    % Single detection - assume constant velocity
    predicted_velocity = velocities(end);
end

end

function best_idx = find_best_track_match(compatible_tracks, tracks, current_time, current_velocity, current_magnitude)
%% Find best matching track using multi-criteria scoring

best_score = inf;
best_idx = compatible_tracks(1);

for i = 1:length(compatible_tracks)
    track_idx = compatible_tracks(i);
    track = tracks(track_idx);
    last_detection = track.detections(end);
    
    % Distance metrics
    time_diff = abs(current_time - last_detection.time);
    velocity_diff = abs(current_velocity - last_detection.velocity);
    magnitude_diff = abs(current_magnitude - last_detection.magnitude);
    
    % Weighted combined score (lower is better)
    score = 0.4 * velocity_diff + 0.3 * time_diff * 1000 + 0.3 * magnitude_diff / 10;
    
    if score < best_score
        best_score = score;
        best_idx = track_idx;
    end
end

end

function quality_score = calculate_track_quality(track, config)
%% Calculate comprehensive track quality score (0-1, higher is better)

% Extract track characteristics
duration = track.duration;
velocity_std = track.velocity_std;
mean_magnitude = track.mean_magnitude;
num_detections = track.num_detections;
mean_velocity = abs(track.mean_velocity);

% Quality components

% 1. Duration score (target: >0.2 seconds for good tracks)
duration_score = min(duration / 0.2, 1);

% 2. Velocity consistency score (penalize erratic motion)
velocity_consistency_score = exp(-velocity_std / 8); % 8 mph std = ~0.37 score

% 3. Signal strength score (normalize to typical golf radar levels)
magnitude_score = max(0, min(1, (mean_magnitude - 80) / 40)); % 80-120 dB range

% 4. Detection density score (more detections = better)
density_score = min(num_detections / 10, 1); % 10 detections = full score

% 5. Golf realism score (reasonable velocities for golf)
if mean_velocity >= 20 && mean_velocity <= 200
    realism_score = 1;
elseif mean_velocity >= 10 && mean_velocity <= 250
    realism_score = 0.7;
else
    realism_score = 0.2;
end

% Combined quality score (weighted average)
quality_score = 0.25 * duration_score + ...
                0.25 * velocity_consistency_score + ...
                0.20 * magnitude_score + ...
                0.15 * density_score + ...
                0.15 * realism_score;

% Ensure score is in [0,1] range
quality_score = max(0, min(1, quality_score));

end
