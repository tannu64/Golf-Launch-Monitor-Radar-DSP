function [enhanced_stft_result, cfar_stats] = milestone2_cfar_detection(stft_result, config)
%% Milestone 2: Adaptive CFAR Detection
% Enhances the baseline STFT result with adaptive thresholding
%
% Inputs:
%   stft_result - Output from baseline_stft_processing()
%   config - Configuration structure with CFAR parameters
%
% Outputs:
%   enhanced_stft_result - Enhanced result with CFAR peaks
%   cfar_stats - CFAR performance statistics

fprintf('   M2: Applying adaptive CFAR detection...\n');

% Add CFAR-specific config parameters if not present
if ~isfield(config, 'cfar_pfa')
    config.cfar_pfa = 1e-4; % Probability of false alarm (0.01%)
end
if ~isfield(config, 'cfar_guard_cells')
    config.cfar_guard_cells = 3; % Guard cells around test cell
end
if ~isfield(config, 'cfar_reference_cells')
    config.cfar_reference_cells = 16; % Reference cells for noise estimation
end

% Extract data from baseline STFT result
magnitude_db = stft_result.magnitude_db_ch1;
[nFreq, nTime] = size(magnitude_db);

% Create frequency mask for golf-relevant range (50-1500 Hz)
freq_resolution = config.sampling_freq / config.fft_size;
freq_axis = (0:nFreq-1) * freq_resolution;
golf_freq_mask = freq_axis >= 50 & freq_axis <= 1500;

% Initialize CFAR threshold matrix
cfar_threshold = zeros(size(magnitude_db));
cfar_detections = false(size(magnitude_db));

% CFAR parameters
guard_cells = config.cfar_guard_cells;
ref_cells = config.cfar_reference_cells;
pfa = config.cfar_pfa;

% Calculate CFAR scaling factor
alpha = ref_cells * (pfa^(-1/ref_cells) - 1);
cfar_factor_db = 10 * log10(alpha);

% Apply CFAR detection
detection_count = 0;
for t = 1:nTime
    for f = 1:nFreq
        if golf_freq_mask(f) % Only process golf-relevant frequencies
            
            % Define reference window boundaries
            f_start = max(1, f - guard_cells - ref_cells/2);
            f_end = min(nFreq, f + guard_cells + ref_cells/2);
            
            % Exclude guard cells and test cell
            f_guard_start = f - guard_cells;
            f_guard_end = f + guard_cells;
            
            % Extract reference cells (excluding guard region)
            ref_indices = [f_start:(f_guard_start-1), (f_guard_end+1):f_end];
            ref_indices = ref_indices(ref_indices >= 1 & ref_indices <= nFreq);
            
            if length(ref_indices) >= ref_cells/2 % Minimum reference cells
                % Calculate noise estimate from reference cells
                reference_values = magnitude_db(ref_indices, t);
                noise_estimate = mean(reference_values);
                
                % Calculate adaptive threshold
                cfar_threshold(f, t) = noise_estimate + cfar_factor_db;
                
                % Apply threshold test
                if magnitude_db(f, t) > cfar_threshold(f, t)
                    cfar_detections(f, t) = true;
                    detection_count = detection_count + 1;
                end
            else
                % Fallback to fixed threshold if insufficient reference cells
                cfar_threshold(f, t) = median(magnitude_db(:, t)) + 6; % 6 dB above median
                if magnitude_db(f, t) > cfar_threshold(f, t)
                    cfar_detections(f, t) = true;
                    detection_count = detection_count + 1;
                end
            end
        end
    end
end

% Extract CFAR peak information
[peak_freq_idx, peak_time_idx] = find(cfar_detections);

if ~isempty(peak_freq_idx)
    % Convert to physical units
    peak_frequencies = freq_axis(peak_freq_idx)';
    peak_velocities = peak_frequencies * config.speed_coef;
    
    % Get time axis from STFT result
    if isfield(stft_result, 'processing_stats')
        time_resolution = stft_result.processing_stats.time_resolution;
        peak_times = (peak_time_idx - 1) * time_resolution;
    else
        % Fallback time calculation
        window_step = config.stft_window - config.stft_overlap;
        time_resolution = window_step / config.sampling_freq;
        peak_times = (peak_time_idx - 1) * time_resolution;
    end
    
    % Get peak magnitudes
    peak_magnitudes = magnitude_db(sub2ind(size(magnitude_db), peak_freq_idx, peak_time_idx));
    
    % Limit number of peaks for processing efficiency
    max_peaks = 2000;
    if length(peak_frequencies) > max_peaks
        [~, sort_idx] = sort(peak_magnitudes, 'descend');
        keep_idx = sort_idx(1:max_peaks);
        peak_frequencies = peak_frequencies(keep_idx);
        peak_velocities = peak_velocities(keep_idx);
        peak_times = peak_times(keep_idx);
        peak_magnitudes = peak_magnitudes(keep_idx);
    end
else
    peak_frequencies = [];
    peak_velocities = [];
    peak_times = [];
    peak_magnitudes = [];
end

% Calculate CFAR statistics
cfar_stats = struct();
cfar_stats.total_detections = detection_count;
cfar_stats.detection_rate = detection_count / sum(golf_freq_mask(:)) / nTime;
cfar_stats.pfa_setting = pfa;
cfar_stats.cfar_factor_db = cfar_factor_db;
cfar_stats.avg_threshold = mean(cfar_threshold(golf_freq_mask));
cfar_stats.peaks_after_limit = length(peak_frequencies);

% Create enhanced STFT result (copy from baseline and add CFAR data)
enhanced_stft_result = stft_result;

% Update peaks with CFAR detections
enhanced_stft_result.peaks_cfar = struct();
enhanced_stft_result.peaks_cfar.frequencies = peak_frequencies;
enhanced_stft_result.peaks_cfar.velocities = peak_velocities;
enhanced_stft_result.peaks_cfar.times = peak_times;
enhanced_stft_result.peaks_cfar.magnitudes = peak_magnitudes;

% Add CFAR-specific data
enhanced_stft_result.cfar_data = struct();
enhanced_stft_result.cfar_data.threshold_matrix = cfar_threshold;
enhanced_stft_result.cfar_data.detection_matrix = cfar_detections;
enhanced_stft_result.cfar_data.config = config;

% Keep original peaks for comparison
enhanced_stft_result.peaks_baseline = stft_result.peaks;

fprintf('   M2 CFAR: %d detections (PFA=%.0e, threshold=%.1f dB)\n', ...
    detection_count, pfa, cfar_stats.avg_threshold);

end
