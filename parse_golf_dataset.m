function [shot_data, stats] = parse_golf_dataset(data_path)
%% Parse Golf Dataset
% Extracts shot information from folder structure and JSON metadata
%
% Inputs:
%   data_path - Path to dataset directory
%
% Outputs:
%   shot_data - Array of structures containing shot information
%   stats - Parsing statistics and data summary

fprintf('Parsing golf dataset from: %s\n', data_path);

% Get all log folders
folder_pattern = fullfile(data_path, 'LOG*');
log_folders = dir(folder_pattern);
log_folders = log_folders([log_folders.isdir]);

fprintf('Found %d log folders\n', length(log_folders));

% Initialize output
shot_data = [];
parsing_errors = 0;
ball_speeds = [];
club_speeds = [];

for i = 1:length(log_folders)
    try
        folder_name = log_folders(i).name;
        folder_path = fullfile(data_path, folder_name);
        
        % Parse folder name to extract TrackMan reference speeds
        [ball_speed, club_speed] = parse_folder_name(folder_name);
        
        % Find JSON and BIN files
        json_files = dir(fullfile(folder_path, '*.json'));
        bin_files = dir(fullfile(folder_path, '*.bin'));
        
        if isempty(json_files) || isempty(bin_files)
            fprintf('Warning: Missing files in folder %s\n', folder_name);
            parsing_errors = parsing_errors + 1;
            continue;
        end
        
        % Parse JSON metadata
        json_file = fullfile(folder_path, json_files(1).name);
        metadata = parse_json_metadata(json_file);
        
        % Create shot data structure
        shot_info = struct();
        shot_info.folder_name = folder_name;
        shot_info.folder_path = folder_path;
        shot_info.json_file = json_file;
        shot_info.bin_file = fullfile(folder_path, bin_files(1).name);
        shot_info.ball_speed_trackman = ball_speed;
        shot_info.club_speed_trackman = club_speed;
        shot_info.metadata = metadata;
        shot_info.timestamp = extract_timestamp(folder_name);
        
        % Determine shot type based on ball speed
        shot_info.shot_type = classify_shot_type(ball_speed);
        
        % Add to dataset
        if isempty(shot_data)
            shot_data = shot_info;
        else
            shot_data(end+1) = shot_info;
        end
        
        % Collect speed statistics (only for valid shots with reasonable values)
        if ~isnan(ball_speed) && ~isnan(club_speed) && ...
           ball_speed > 0 && ball_speed < 300 && ...  % Reasonable ball speed range
           club_speed > 0 && club_speed < 200        % Reasonable club speed range
            ball_speeds(end+1) = ball_speed;
            club_speeds(end+1) = club_speed;
        end
        
    catch ME
        fprintf('Error parsing folder %s: %s\n', folder_name, ME.message);
        parsing_errors = parsing_errors + 1;
    end
end

% Calculate statistics
total_folders = length(log_folders);
successful_parses = length(shot_data);
success_rate = successful_parses / total_folders;

stats = struct();
stats.total_folders = total_folders;
stats.successful_parses = successful_parses;
stats.parsing_errors = parsing_errors;
stats.success_rate = success_rate;

if ~isempty(ball_speeds)
    stats.ball_speed_range = [min(ball_speeds), max(ball_speeds)];
    stats.club_speed_range = [min(club_speeds), max(club_speeds)];
    stats.ball_speed_mean = mean(ball_speeds);
    stats.club_speed_mean = mean(club_speeds);
    
    % Shot type distribution
    shot_types = {shot_data.shot_type};
    stats.wedge_count = sum(strcmp(shot_types, 'wedge'));
    stats.iron_count = sum(strcmp(shot_types, 'iron'));
    stats.driver_count = sum(strcmp(shot_types, 'driver'));
else
    stats.ball_speed_range = [NaN, NaN];
    stats.club_speed_range = [NaN, NaN];
    stats.ball_speed_mean = NaN;
    stats.club_speed_mean = NaN;
    stats.wedge_count = 0;
    stats.iron_count = 0;
    stats.driver_count = 0;
end

fprintf('Dataset parsing complete:\n');
fprintf('- Successfully parsed: %d/%d folders (%.1f%%)\n', ...
    successful_parses, total_folders, success_rate * 100);
fprintf('- Wedge shots: %d, Iron shots: %d, Driver shots: %d\n', ...
    stats.wedge_count, stats.iron_count, stats.driver_count);

end

function [ball_speed, club_speed] = parse_folder_name(folder_name)
%% Parse TrackMan speeds from folder name
% Expected format: LOG[timestamp]_B[ball_speed]-H[club_speed]

ball_speed = NaN;
club_speed = NaN;

% Handle special case for missed shots
if contains(folder_name, 'BMISS') || contains(folder_name, 'HMISS')
    return;
end

% Extract ball speed
ball_match = regexp(folder_name, 'B([\d.]+)', 'tokens');
if ~isempty(ball_match)
    ball_speed = str2double(ball_match{1}{1});
end

% Extract club head speed
club_match = regexp(folder_name, 'H([\d.]+)', 'tokens');
if ~isempty(club_match)
    club_speed = str2double(club_match{1}{1});
end

end

function metadata = parse_json_metadata(json_file)
%% Parse JSON metadata file

try
    % Read JSON file and handle BOM
    fid = fopen(json_file, 'r', 'n', 'UTF-8');
    if fid == -1
        error('Cannot open JSON file: %s', json_file);
    end
    json_text = fread(fid, '*char')';
    fclose(fid);
    
    % Remove BOM if present
    if length(json_text) >= 3 && strcmp(json_text(1:3), char([239, 187, 191]))
        json_text = json_text(4:end);
    end
    
    metadata = jsondecode(json_text);
    
    % Ensure all expected fields are present
    required_fields = {'samplingFreq', 'speedCoef', 'fftNum', 'dataStep', ...
                      'binFile_DataArray', 'dataMode'};
    
    for i = 1:length(required_fields)
        if ~isfield(metadata, required_fields{i})
            warning('Missing field %s in JSON metadata', required_fields{i});
        end
    end
    
catch ME
    warning('Error parsing JSON file %s: %s', json_file, ME.message);
    metadata = struct();
end

end

function timestamp = extract_timestamp(folder_name)
%% Extract timestamp from folder name

timestamp_match = regexp(folder_name, 'LOG(\d{14})', 'tokens');
if ~isempty(timestamp_match)
    timestamp_str = timestamp_match{1}{1};
    % Convert to MATLAB datetime
    timestamp = datetime(timestamp_str, 'InputFormat', 'yyyyMMddHHmmss');
else
    timestamp = datetime.empty;
end

end

function shot_type = classify_shot_type(ball_speed)
%% Classify shot type based on ball speed

if isnan(ball_speed)
    shot_type = 'unknown';
elseif ball_speed < 100
    shot_type = 'wedge';
elseif ball_speed < 140
    shot_type = 'iron';
else
    shot_type = 'driver';
end

end
