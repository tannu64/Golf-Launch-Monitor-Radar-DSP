function selected_shots = select_representative_shots(shot_data, num_shots)
%% Select Representative Shots for Analysis
% Selects a diverse set of shots covering different club types and speeds
%
% Inputs:
%   shot_data - Array of shot structures from parse_golf_dataset
%   num_shots - Number of shots to select
%
% Outputs:
%   selected_shots - Array of selected shot structures

if isempty(shot_data)
    selected_shots = [];
    return;
end

% Filter out invalid shots (NaN speeds)
valid_shots = shot_data(~isnan([shot_data.ball_speed_trackman]));

if length(valid_shots) <= num_shots
    selected_shots = valid_shots;
    return;
end

% Categorize shots by type
wedge_shots = valid_shots(strcmp({valid_shots.shot_type}, 'wedge'));
iron_shots = valid_shots(strcmp({valid_shots.shot_type}, 'iron'));
driver_shots = valid_shots(strcmp({valid_shots.shot_type}, 'driver'));

% Calculate target distribution
wedge_target = ceil(num_shots * 0.4);  % 40% wedges (focus area)
iron_target = ceil(num_shots * 0.4);   % 40% irons
driver_target = num_shots - wedge_target - iron_target; % Remaining for drivers

% Select shots from each category
selected_wedges = select_from_category(wedge_shots, wedge_target);
selected_irons = select_from_category(iron_shots, iron_target);
selected_drivers = select_from_category(driver_shots, driver_target);

% Combine selections
selected_shots = [selected_wedges, selected_irons, selected_drivers];

% If we don't have enough shots in some categories, fill from others
if length(selected_shots) < num_shots
    remaining_shots = setdiff(valid_shots, selected_shots);
    additional_needed = num_shots - length(selected_shots);
    
    if length(remaining_shots) >= additional_needed
        % Select additional shots spanning the speed range
        additional_shots = select_by_speed_distribution(remaining_shots, additional_needed);
        selected_shots = [selected_shots, additional_shots];
    else
        selected_shots = [selected_shots, remaining_shots];
    end
end

% Ensure we don't exceed the requested number
if length(selected_shots) > num_shots
    selected_shots = selected_shots(1:num_shots);
end

% Sort by timestamp for logical processing order
if ~isempty(selected_shots) && ~isempty([selected_shots.timestamp])
    [~, sort_idx] = sort([selected_shots.timestamp]);
    selected_shots = selected_shots(sort_idx);
end

fprintf('Selected %d representative shots:\n', length(selected_shots));
if ~isempty(selected_shots)
    fprintf('- Wedges: %d\n', sum(strcmp({selected_shots.shot_type}, 'wedge')));
    fprintf('- Irons: %d\n', sum(strcmp({selected_shots.shot_type}, 'iron')));
    fprintf('- Drivers: %d\n', sum(strcmp({selected_shots.shot_type}, 'driver')));
    
    ball_speeds = [selected_shots.ball_speed_trackman];
    fprintf('- Ball speed range: %.1f - %.1f mph\n', min(ball_speeds), max(ball_speeds));
end

end

function selected = select_from_category(category_shots, target_count)
%% Select shots from a specific category

if isempty(category_shots) || target_count <= 0
    selected = [];
    return;
end

if length(category_shots) <= target_count
    selected = category_shots;
    return;
end

% Select shots to cover the speed range within the category
ball_speeds = [category_shots.ball_speed_trackman];
[~, sort_idx] = sort(ball_speeds);
sorted_shots = category_shots(sort_idx);

% Select evenly distributed across the speed range
indices = round(linspace(1, length(sorted_shots), target_count));
selected = sorted_shots(indices);

end

function selected = select_by_speed_distribution(shots, target_count)
%% Select shots to achieve good speed distribution

if isempty(shots) || target_count <= 0
    selected = [];
    return;
end

if length(shots) <= target_count
    selected = shots;
    return;
end

% Sort by ball speed
ball_speeds = [shots.ball_speed_trackman];
[~, sort_idx] = sort(ball_speeds);
sorted_shots = shots(sort_idx);

% Select evenly distributed across the entire speed range
indices = round(linspace(1, length(sorted_shots), target_count));
selected = sorted_shots(indices);

end
