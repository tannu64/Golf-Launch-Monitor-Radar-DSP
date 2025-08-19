function signal_quality = analyze_signal_quality(stft_results, config)
%% Analyze Signal Quality
% Performs comprehensive signal quality assessment on processed radar data
%
% Inputs:
%   stft_results - Cell array of STFT processing results
%   config - Configuration structure
%
% Outputs:
%   signal_quality - Structure containing signal quality metrics

fprintf('Analyzing signal quality across %d shots...\n', length(stft_results));

% Initialize metric arrays
snr_values = zeros(length(stft_results), 1);
noise_floors = zeros(length(stft_results), 1);
peak_counts = zeros(length(stft_results), 1);
signal_bandwidths = zeros(length(stft_results), 1);
saturation_rates = zeros(length(stft_results), 2); % [Ch1, Ch2]
dc_offsets = zeros(length(stft_results), 4); % [Ch1_I, Ch1_Q, Ch2_I, Ch2_Q]

for i = 1:length(stft_results)
    try
        stft_data = stft_results{i}.stft_data;
        iq_data = stft_results{i}.iq_data;
        
        % Extract signal quality metrics
        snr_values(i) = calculate_snr(stft_data);
        noise_floors(i) = stft_data.processing_stats.noise_floor_db;
        peak_counts(i) = stft_data.processing_stats.num_peaks_detected;
        signal_bandwidths(i) = estimate_signal_bandwidth(stft_data);
        
        % I/Q data quality metrics
        saturation_rates(i, :) = [iq_data.saturated_samples_ch1, iq_data.saturated_samples_ch2];
        dc_offsets(i, :) = [iq_data.dc_offset_ch1, iq_data.dc_offset_ch2];
        
    catch ME
        fprintf('Warning: Error analyzing shot %d: %s\n', i, ME.message);
        snr_values(i) = NaN;
        noise_floors(i) = NaN;
        peak_counts(i) = 0;
        signal_bandwidths(i) = NaN;
        saturation_rates(i, :) = [NaN, NaN];
        dc_offsets(i, :) = [NaN, NaN, NaN, NaN];
    end
end

% Calculate overall quality metrics
valid_shots = ~isnan(snr_values);
num_valid = sum(valid_shots);

if num_valid > 0
    avg_snr = mean(snr_values(valid_shots));
    snr_std = std(snr_values(valid_shots));
    avg_noise_floor = mean(noise_floors(valid_shots));
    avg_peak_count = mean(peak_counts(valid_shots));
    avg_bandwidth = mean(signal_bandwidths(valid_shots));
    
    % Peak detection success rate
    peak_detection_rate = mean(peak_counts > 0);
    
    % Saturation analysis
    max_saturation_ch1 = max(saturation_rates(:, 1));
    max_saturation_ch2 = max(saturation_rates(:, 2));
    avg_saturation = mean(saturation_rates(:));
    
    % DC offset analysis
    max_dc_offset = max(abs(dc_offsets(:)));
    avg_dc_offset = mean(abs(dc_offsets(:)));
    
else
    avg_snr = NaN;
    snr_std = NaN;
    avg_noise_floor = NaN;
    avg_peak_count = 0;
    avg_bandwidth = NaN;
    peak_detection_rate = 0;
    max_saturation_ch1 = NaN;
    max_saturation_ch2 = NaN;
    avg_saturation = NaN;
    max_dc_offset = NaN;
    avg_dc_offset = NaN;
end

% Signal quality assessment
quality_flags = struct();
quality_flags.good_snr = avg_snr > 15; % dB
quality_flags.low_saturation = avg_saturation < 0.01; % <1%
quality_flags.reasonable_dc = max_dc_offset < 1000; % Arbitrary threshold
quality_flags.good_peak_detection = peak_detection_rate > 0.8; % >80%

% Overall quality score (0-1)
quality_score = sum(struct2array(quality_flags)) / length(fieldnames(quality_flags));

% Compile results
signal_quality = struct();
signal_quality.avg_snr = avg_snr;
signal_quality.snr_std = snr_std;
signal_quality.noise_floor = avg_noise_floor;
signal_quality.signal_bandwidth = avg_bandwidth;
signal_quality.peak_detection_rate = peak_detection_rate;
signal_quality.avg_peak_count = avg_peak_count;
signal_quality.max_saturation_ch1 = max_saturation_ch1;
signal_quality.max_saturation_ch2 = max_saturation_ch2;
signal_quality.avg_saturation = avg_saturation;
signal_quality.max_dc_offset = max_dc_offset;
signal_quality.avg_dc_offset = avg_dc_offset;
signal_quality.quality_flags = quality_flags;
signal_quality.overall_quality_score = quality_score;
signal_quality.num_valid_shots = num_valid;

% Detailed arrays for further analysis
signal_quality.detailed = struct();
signal_quality.detailed.snr_values = snr_values;
signal_quality.detailed.noise_floors = noise_floors;
signal_quality.detailed.peak_counts = peak_counts;
signal_quality.detailed.signal_bandwidths = signal_bandwidths;
signal_quality.detailed.saturation_rates = saturation_rates;
signal_quality.detailed.dc_offsets = dc_offsets;

fprintf('Signal quality analysis complete:\n');
fprintf('- Valid shots: %d/%d\n', num_valid, length(stft_results));
fprintf('- Average SNR: %.1f ± %.1f dB\n', avg_snr, snr_std);
fprintf('- Peak detection success: %.1f%%\n', peak_detection_rate * 100);
fprintf('- Overall quality score: %.2f/1.0\n', quality_score);

end

function snr_db = calculate_snr(stft_data)
%% Calculate Signal-to-Noise Ratio

magnitude_db = stft_data.magnitude_db_ch1;
noise_floor = stft_data.processing_stats.noise_floor_db;

% Find signal regions (above threshold)
threshold = noise_floor + 6; % 6 dB above noise floor
signal_mask = magnitude_db > threshold;

if sum(signal_mask(:)) > 0
    % Calculate mean signal power in signal regions
    signal_power = mean(magnitude_db(signal_mask));
    snr_db = signal_power - noise_floor;
else
    snr_db = 0; % No signal detected
end

end

function bandwidth = estimate_signal_bandwidth(stft_data)
%% Estimate Signal Bandwidth

magnitude_db = stft_data.magnitude_db_ch1;
freq_axis = stft_data.processing_stats.frequency_resolution * (0:size(magnitude_db, 1)-1);

% Average over time to get frequency profile
avg_spectrum = mean(magnitude_db, 2);

% Find -3dB bandwidth around peak
[max_power, max_idx] = max(avg_spectrum);
threshold_3db = max_power - 3;

% Find frequency range above threshold
above_threshold = avg_spectrum > threshold_3db;
freq_indices = find(above_threshold);

if length(freq_indices) > 1
    bandwidth = freq_axis(freq_indices(end)) - freq_axis(freq_indices(1));
else
    bandwidth = 0;
end

end
