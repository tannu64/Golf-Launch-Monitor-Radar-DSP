function iq_data = load_iq_data(bin_file, config)
%% Load I/Q Data from Binary File
% Reads 4-channel I/Q radar data from binary file
%
% Inputs:
%   bin_file - Path to binary data file
%   config - Configuration structure with radar parameters
%
% Outputs:
%   iq_data - Structure containing processed I/Q data

% Validate inputs
if ~exist(bin_file, 'file')
    error('Binary file not found: %s', bin_file);
end

% Get file information
file_info = dir(bin_file);
file_size_bytes = file_info.bytes;

% Calculate expected data format
% Assuming 4 channels (IF1 Upper, IF1 Lower, IF2 Upper, IF2 Lower)
% Data format: 16-bit signed integers (2 bytes per sample)
bytes_per_sample = 2;
channels = 4;
bytes_per_frame = channels * bytes_per_sample;
total_samples = file_size_bytes / bytes_per_frame;

fprintf('Loading I/Q data from: %s\n', bin_file);
fprintf('File size: %.2f MB\n', file_size_bytes / 1024 / 1024);
fprintf('Expected samples: %d\n', total_samples);

try
    % Open binary file
    fid = fopen(bin_file, 'r');
    if fid == -1
        error('Cannot open binary file: %s', bin_file);
    end
    
    % Read all data as 16-bit signed integers
    raw_data = fread(fid, [channels, total_samples], 'int16');
    fclose(fid);
    
    % Separate channels
    if1_upper = raw_data(1, :);
    if1_lower = raw_data(2, :);
    if2_upper = raw_data(3, :);
    if2_lower = raw_data(4, :);
    
    % Form complex I/Q pairs
    % Channel 1: IF1 Upper (I) + j*IF1 Lower (Q)
    % Channel 2: IF2 Upper (I) + j*IF2 Lower (Q)
    iq_ch1 = complex(if1_upper, if1_lower);
    iq_ch2 = complex(if2_upper, if2_lower);
    
    % Create time axis
    time_axis = (0:length(iq_ch1)-1) / config.sampling_freq;
    
    % Calculate signal statistics
    power_ch1 = mean(abs(iq_ch1).^2);
    power_ch2 = mean(abs(iq_ch2).^2);
    
    % Store in output structure
    iq_data = struct();
    iq_data.channel1 = iq_ch1;
    iq_data.channel2 = iq_ch2;
    iq_data.raw_data = raw_data;
    iq_data.time_axis = time_axis;
    iq_data.sampling_freq = config.sampling_freq;
    iq_data.duration = length(iq_ch1) / config.sampling_freq;
    iq_data.num_samples = length(iq_ch1);
    iq_data.power_ch1 = power_ch1;
    iq_data.power_ch2 = power_ch2;
    iq_data.file_path = bin_file;
    
    % Data quality checks
    iq_data.saturated_samples_ch1 = check_saturation(if1_upper, if1_lower);
    iq_data.saturated_samples_ch2 = check_saturation(if2_upper, if2_lower);
    iq_data.dc_offset_ch1 = [mean(if1_upper), mean(if1_lower)];
    iq_data.dc_offset_ch2 = [mean(if2_upper), mean(if2_lower)];
    
    fprintf('I/Q data loaded successfully:\n');
    fprintf('- Duration: %.3f seconds\n', iq_data.duration);
    fprintf('- Samples: %d\n', iq_data.num_samples);
    fprintf('- Channel 1 power: %.2e\n', power_ch1);
    fprintf('- Channel 2 power: %.2e\n', power_ch2);
    fprintf('- Saturation Ch1: %.2f%%, Ch2: %.2f%%\n', ...
        iq_data.saturated_samples_ch1 * 100, iq_data.saturated_samples_ch2 * 100);
    
catch ME
    if exist('fid', 'var') && fid ~= -1
        fclose(fid);
    end
    error('Error reading I/Q data: %s', ME.message);
end

end

function saturation_rate = check_saturation(i_data, q_data)
%% Check for signal saturation
% Looks for samples at maximum/minimum values indicating ADC saturation

max_val = 2^15 - 1;  % Maximum value for 16-bit signed integer
min_val = -2^15;     % Minimum value for 16-bit signed integer

saturated_i = sum(i_data == max_val | i_data == min_val);
saturated_q = sum(q_data == max_val | q_data == min_val);
total_samples = length(i_data) + length(q_data);

saturation_rate = (saturated_i + saturated_q) / total_samples;

end
