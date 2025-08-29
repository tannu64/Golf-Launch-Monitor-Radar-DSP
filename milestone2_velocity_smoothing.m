function [smoothed_tracks, smoothing_stats] = milestone2_velocity_smoothing(raw_tracks, config)
%% Milestone 2: Kalman Filter Velocity Smoothing
% Applies temporal smoothing to improve velocity accuracy
%
% Inputs:
%   raw_tracks - Raw tracked objects from advanced tracking
%   config - Configuration with smoothing parameters
%
% Outputs:
%   smoothed_tracks - Velocity-smoothed tracks
%   smoothing_stats - Smoothing performance metrics

fprintf('   M2: Applying Kalman velocity smoothing...\n');

if isempty(raw_tracks)
    smoothed_tracks = [];
    smoothing_stats = struct('tracks_processed', 0, 'avg_improvement', 0);
    fprintf('   M2 Smoothing: No tracks to process\n');
    return;
end

% Kalman filter parameters
if ~isfield(config, 'm2_process_noise')
    config.m2_process_noise = 25; % Q - acceleration variance (mph²/s²)
end
if ~isfield(config, 'm2_measurement_noise')
    config.m2_measurement_noise = 9; % R - velocity measurement variance (mph²)
end
if ~isfield(config, 'm2_initial_uncertainty')
    config.m2_initial_uncertainty = 100; % P0 - initial state uncertainty
end

% Initialize statistics
total_improvement = 0;
tracks_processed = 0;

% Process each track individually
smoothed_tracks = raw_tracks; % Copy structure

for i = 1:length(raw_tracks)
    track = raw_tracks(i);
    
    if length(track.detections) >= 3 % Need minimum points for smoothing
        
        % Extract time series data
        times = [track.detections.time];
        velocities = [track.detections.velocity];
        magnitudes = [track.detections.magnitude];
        
        % Apply Kalman smoothing
        [smoothed_velocities, smoothing_improvement] = kalman_velocity_filter(times, velocities, config);
        
        % Update track with smoothed velocities
        for j = 1:length(track.detections)
            smoothed_tracks(i).detections(j).velocity = smoothed_velocities(j);
            smoothed_tracks(i).detections(j).velocity_raw = velocities(j); % Keep original
        end
        
        % Update track statistics with smoothed data
        smoothed_tracks(i).mean_velocity = mean(smoothed_velocities);
        smoothed_tracks(i).velocity_std = std(smoothed_velocities);
        smoothed_tracks(i).velocity_std_raw = std(velocities); % Keep original std
        
        % Track the improvement
        total_improvement = total_improvement + smoothing_improvement;
        tracks_processed = tracks_processed + 1;
        
    else
        % Too few points for smoothing - keep original
        smoothed_tracks(i).velocity_std_raw = track.velocity_std;
    end
end

% Calculate smoothing statistics
smoothing_stats = struct();
smoothing_stats.tracks_processed = tracks_processed;
smoothing_stats.tracks_total = length(raw_tracks);

if tracks_processed > 0
    smoothing_stats.avg_improvement = total_improvement / tracks_processed;
    
    % Calculate overall velocity standard deviation improvement
    % Safely extract velocity standard deviations
    raw_stds = [];
    smoothed_stds = [];
    
    for i = 1:length(raw_tracks)
        if isfield(raw_tracks(i), 'velocity_std') && isnumeric(raw_tracks(i).velocity_std)
            raw_stds(end+1) = raw_tracks(i).velocity_std;
        else
            raw_stds(end+1) = NaN;
        end
    end
    
    for i = 1:length(smoothed_tracks)
        if isfield(smoothed_tracks(i), 'velocity_std') && isnumeric(smoothed_tracks(i).velocity_std)
            smoothed_stds(end+1) = smoothed_tracks(i).velocity_std;
        else
            smoothed_stds(end+1) = NaN;
        end
    end
    
    valid_mask = ~isnan(raw_stds) & ~isnan(smoothed_stds);
    if sum(valid_mask) > 0
        smoothing_stats.std_improvement = mean(raw_stds(valid_mask)) - mean(smoothed_stds(valid_mask));
    else
        smoothing_stats.std_improvement = 0;
    end
else
    smoothing_stats.avg_improvement = 0;
    smoothing_stats.std_improvement = 0;
end

fprintf('   M2 Smoothing: %d tracks processed, %.2f mph avg std reduction\n', ...
    tracks_processed, smoothing_stats.std_improvement);

end

function [smoothed_velocities, improvement] = kalman_velocity_filter(times, velocities, config)
%% 1D Kalman filter for velocity smoothing
% State: [velocity, acceleration]
% Measurement: velocity

% Kalman filter parameters
Q = config.m2_process_noise; % Process noise (acceleration variance)
R = config.m2_measurement_noise; % Measurement noise (velocity measurement variance)
P0 = config.m2_initial_uncertainty; % Initial uncertainty

n = length(times);
if n < 2
    smoothed_velocities = velocities;
    improvement = 0;
    return;
end

% State: x = [velocity; acceleration]
% State transition model: velocity(k+1) = velocity(k) + acceleration(k) * dt
%                        acceleration(k+1) = acceleration(k) + noise

% Measurement model: z = [1 0] * x (we measure velocity directly)
H = [1, 0]; % Measurement matrix

% Initialize state
x = [velocities(1); 0]; % Initial velocity, zero acceleration
P = [P0, 0; 0, P0]; % Initial covariance

smoothed_velocities = zeros(size(velocities));
smoothed_velocities(1) = velocities(1);

for k = 2:n
    % Time step
    dt = times(k) - times(k-1);
    
    % State transition matrix
    F = [1, dt; 0, 1];
    
    % Process noise covariance
    Q_k = [dt^3/3, dt^2/2; dt^2/2, dt] * Q;
    
    % Prediction step
    x_pred = F * x;
    P_pred = F * P * F' + Q_k;
    
    % Measurement update
    y = velocities(k) - H * x_pred; % Innovation
    S = H * P_pred * H' + R; % Innovation covariance
    K = P_pred * H' / S; % Kalman gain
    
    % Update state
    x = x_pred + K * y;
    P = (eye(2) - K * H) * P_pred;
    
    % Store smoothed velocity
    smoothed_velocities(k) = x(1);
end

% Calculate improvement metric (reduction in velocity standard deviation)
original_std = std(velocities);
smoothed_std = std(smoothed_velocities);
improvement = original_std - smoothed_std;

end
