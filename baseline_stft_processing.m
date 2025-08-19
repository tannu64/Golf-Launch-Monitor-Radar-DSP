function [stft_result, freq_axis, time_axis] = baseline_stft_processing(iq_data, config)
%% Baseline STFT Processing Pipeline
% Performs Short-Time Fourier Transform on I/Q radar data for Doppler analysis
%
% Inputs:
%   iq_data - Structure containing I/Q data from load_iq_data
%   config - Configuration structure with processing parameters
%
% Outputs:
%   stft_result - Structure containing STFT analysis results
%   freq_axis - Frequency axis in Hz
%   time_axis - Time axis in seconds

fprintf('Performing STFT processing...\n');

% Extract parameters
window_size = config.stft_window;
overlap = config.stft_overlap;
sampling_freq = config.sampling_freq;
speed_coef = config.speed_coef;

% Select primary channel for analysis (Channel 1)
signal = iq_data.channel1;

% Apply window function (Hamming window for spectral leakage reduction)
window = hamming(window_size);

% Calculate STFT for both channels
[stft_ch1, freq_stft, time_stft] = stft(signal, sampling_freq, ...
    'Window', window, 'OverlapLength', overlap, 'FFTLength', config.fft_size);

% Also process channel 2 for comparison
[stft_ch2, ~, ~] = stft(iq_data.channel2, sampling_freq, ...
    'Window', window, 'OverlapLength', overlap, 'FFTLength', config.fft_size);

% Convert frequency to velocity using Doppler equation
% For CW radar: v = (λ * f_d) / 2, where λ = c / f_carrier
% speed_coef already incorporates radar-specific conversion factors
velocity_axis = freq_stft * speed_coef;

% Calculate magnitude spectrograms
magnitude_ch1 = abs(stft_ch1);
magnitude_ch2 = abs(stft_ch2);

% Convert to dB scale
magnitude_db_ch1 = 20 * log10(magnitude_ch1 + eps);
magnitude_db_ch2 = 20 * log10(magnitude_ch2 + eps);

% Estimate noise floor
noise_floor_db = estimate_noise_floor(magnitude_db_ch1);

% Apply basic thresholding for peak detection
threshold_db = noise_floor_db + 10; % 10 dB above noise floor
peaks_binary = magnitude_db_ch1 > threshold_db;

% Find peak locations in time-frequency domain
[peak_freq_idx, peak_time_idx] = find(peaks_binary);
peak_frequencies = freq_stft(peak_freq_idx);
peak_velocities = peak_frequencies * speed_coef;
peak_times = time_stft(peak_time_idx);
peak_magnitudes = magnitude_db_ch1(sub2ind(size(magnitude_db_ch1), peak_freq_idx, peak_time_idx));

% Basic peak clustering and tracking
[tracked_objects, tracking_stats] = basic_peak_tracking(peak_velocities, peak_times, peak_magnitudes, config);

% Calculate processing statistics
processing_stats = struct();
processing_stats.window_size = window_size;
processing_stats.overlap = overlap;
processing_stats.fft_size = config.fft_size;
processing_stats.frequency_resolution = sampling_freq / config.fft_size;
processing_stats.time_resolution = (window_size - overlap) / sampling_freq;
processing_stats.velocity_resolution = processing_stats.frequency_resolution * speed_coef;
processing_stats.noise_floor_db = noise_floor_db;
processing_stats.threshold_db = threshold_db;
processing_stats.num_peaks_detected = length(peak_frequencies);

% Store results
stft_result = struct();
stft_result.stft_ch1 = stft_ch1;
stft_result.stft_ch2 = stft_ch2;
stft_result.magnitude_ch1 = magnitude_ch1;
stft_result.magnitude_ch2 = magnitude_ch2;
stft_result.magnitude_db_ch1 = magnitude_db_ch1;
stft_result.magnitude_db_ch2 = magnitude_db_ch2;
stft_result.peaks = struct();
stft_result.peaks.frequencies = peak_frequencies;
stft_result.peaks.velocities = peak_velocities;
stft_result.peaks.times = peak_times;
stft_result.peaks.magnitudes = peak_magnitudes;
stft_result.tracked_objects = tracked_objects;
stft_result.tracking_stats = tracking_stats;
stft_result.processing_stats = processing_stats;

% Output axes
freq_axis = freq_stft;
time_axis = time_stft;

% Add velocity axis for convenience
stft_result.velocity_axis = velocity_axis;

fprintf('STFT processing complete:\n');
fprintf('- Frequency resolution: %.2f Hz (%.2f mph)\n', ...
    processing_stats.frequency_resolution, processing_stats.velocity_resolution);
fprintf('- Time resolution: %.3f s\n', processing_stats.time_resolution);
fprintf('- Noise floor: %.1f dB\n', noise_floor_db);
fprintf('- Peaks detected: %d\n', length(peak_frequencies));
fprintf('- Tracked objects: %d\n', length(tracked_objects));

end

function noise_floor_db = estimate_noise_floor(magnitude_db)
%% Estimate noise floor from magnitude spectrogram
% Uses lower percentile to estimate noise level

% Flatten spectrogram and remove very low values (log artifacts)
magnitude_vector = magnitude_db(:);
magnitude_vector = magnitude_vector(isfinite(magnitude_vector));

% Use 25th percentile as noise floor estimate
noise_floor_db = prctile(magnitude_vector, 25);

end

function [tracked_objects, stats] = basic_peak_tracking(velocities, times, magnitudes, config)
%% Basic Peak Tracking Algorithm
% Groups peaks into velocity tracks over time

% Parameters for tracking
velocity_tolerance = 5; % mph - peaks within this range considered same object
min_track_length = 3; % minimum number of detections for valid track
max_time_gap = 0.1; % seconds - maximum time gap between detections

% Initialize tracking
tracks = struct([]);  % Initialize as empty struct array
track_id = 0;

% Sort detections by time
[sorted_times, sort_idx] = sort(times);
sorted_velocities = velocities(sort_idx);
sorted_magnitudes = magnitudes(sort_idx);

% Simple nearest-neighbor tracking
for i = 1:length(sorted_times)
    current_time = sorted_times(i);
    current_velocity = sorted_velocities(i);
    current_magnitude = sorted_magnitudes(i);
    
    % Find existing tracks that could match this detection
    matching_tracks = [];
    for j = 1:length(tracks)
        last_detection = tracks(j).detections(end);
        time_gap = current_time - last_detection.time;
        velocity_diff = abs(current_velocity - last_detection.velocity);
        
        if time_gap <= max_time_gap && velocity_diff <= velocity_tolerance
            matching_tracks(end+1) = j;
        end
    end
    
    if isempty(matching_tracks)
        % Create new track
        track_id = track_id + 1;
        new_track = struct();
        new_track.id = track_id;
        new_track.detections = struct('time', current_time, 'velocity', current_velocity, 'magnitude', current_magnitude);
        
        if isempty(tracks)
            tracks = new_track;
        else
            tracks(end+1) = new_track;
        end
    else
        % Add to closest existing track (by velocity)
        best_track_idx = matching_tracks(1);
        best_velocity_diff = inf;
        
        for j = matching_tracks
            last_detection = tracks(j).detections(end);
            velocity_diff = abs(current_velocity - last_detection.velocity);
            if velocity_diff < best_velocity_diff
                best_velocity_diff = velocity_diff;
                best_track_idx = j;
            end
        end
        
        % Add detection to best matching track
        new_detection = struct('time', current_time, 'velocity', current_velocity, 'magnitude', current_magnitude);
        tracks(best_track_idx).detections(end+1) = new_detection;
    end
end

% Filter tracks by minimum length and calculate statistics
valid_tracks = [];
for i = 1:length(tracks)
    if length(tracks(i).detections) >= min_track_length
        track = tracks(i);
        
        % Calculate track statistics
        track_times = [track.detections.time];
        track_velocities = [track.detections.velocity];
        track_magnitudes = [track.detections.magnitude];
        
        track.duration = max(track_times) - min(track_times);
        track.mean_velocity = mean(track_velocities);
        track.velocity_std = std(track_velocities);
        track.max_magnitude = max(track_magnitudes);
        track.mean_magnitude = mean(track_magnitudes);
        track.num_detections = length(track.detections);
        
        if isempty(valid_tracks)
            valid_tracks = track;
        else
            valid_tracks(end+1) = track;
        end
    end
end

tracked_objects = valid_tracks;

% Calculate tracking statistics
stats = struct();
stats.total_tracks_created = track_id;
stats.valid_tracks = length(valid_tracks);

if track_id > 0
    stats.tracking_efficiency = length(valid_tracks) / track_id;
else
    stats.tracking_efficiency = 0;
end

if ~isempty(valid_tracks)
    stats.avg_track_length = mean([valid_tracks.num_detections]);
    stats.velocity_range = [min([valid_tracks.mean_velocity]), max([valid_tracks.mean_velocity])];
else
    stats.avg_track_length = 0;
    stats.velocity_range = [NaN, NaN];
end

end
