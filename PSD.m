%% Power Spectral Density (PSD) Analysis Script
% This script performs spectral analysis of turbulent flow data using different methods:
% - Welch's method: Overlapping segments with windowing
% - Periodogram method: Single FFT of entire signal
% - Bartlett method: Non-overlapping segments with rectangular windows
%
% Required Input:
% - savepoints : Contains probe location coordinates (X,Y,Z)
% - Flow data files in format: Flow0_X_Y_Z_dt_*.dat
% - IntegralScales_*.txt files for integral length scales
%
% Author: Hossein Seyedzadeh
% Email: hossein.seyyedzadeh@stonybrook.edu
% Last modified: April 25, 2025

%% Initialization and Setup
clc;
clear;
close all;
tic     % Start timing execution

%% Input Parameters and Configuration
dt = 0.0005;                  % Time step (seconds)
savepoints = load('savepoints'); % Probe locations file
lp = length(savepoints);      % Number of probe locations

% Analysis configuration
case_name = 'HS';            % Options: 'HR', 'HS', 'MR', 'MS'
component = 'w';             % Options: 'u', 'v', 'w', 'vmag', 'tke'
psd_method = 'bartlett';     % Options: 'welch', 'periodogram', 'bartlett'
Lf = 0.06;                   % Fixed fish characteristic length (m)

% PSD calculation parameters
window_size = 4096;          % Window size in samples
overlap_ratio = 0.5;         % Overlap ratio (0-1, used in Welch's method)
nfft = 4096;                 % Number of FFT points (frequency resolution)
cutoff_frequency = 30;      % Maximum frequency to analyze (Hz)

% Frequency band analysis parameters
band_width_factor = 0.2;     % Bandwidth around characteristic frequencies (±20%)

% Set input/output paths based on selected case
path_input = fullfile("Input", case_name, "");
path_output = fullfile("Output", case_name, "");

% Create output directory if it doesn't exist
if ~exist(path_output, 'dir')
    mkdir(path_output);
end

% PSD output subfolder
output_dir = fullfile(path_output, 'PSD_Results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Setup log file
log_filename = fullfile(output_dir, ['PSD_Analysis_Log_', case_name, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.txt']);
log_file = fopen(log_filename, 'w');

% Write analysis parameters to log
fprintf(log_file, '=================================================================\n');
fprintf(log_file, 'PSD ANALYSIS LOG - %s\n', datestr(now));
fprintf(log_file, '=================================================================\n');
fprintf(log_file, 'Case: %s\n', case_name);
fprintf(log_file, 'Input path: %s\n', path_input);
fprintf(log_file, 'Output path: %s\n', output_dir);
fprintf(log_file, 'Component analyzed: %s\n', upper(component));
fprintf(log_file, 'PSD method: %s\n', psd_method);
fprintf(log_file, 'Window size: %d samples\n', window_size);
fprintf(log_file, 'Overlap ratio: %.2f\n', overlap_ratio);
fprintf(log_file, 'NFFT: %d\n', nfft);
fprintf(log_file, 'Probe count: %d\n', lp);
fprintf(log_file, 'Fish length (Lf): %.2f cm\n', Lf*100);
fprintf(log_file, '=================================================================\n\n');

%% Data Loading and Preprocessing
% Create probe mapping table
probe_info = cell(lp,1);

% Define coordinates-to-label mapping
for i = 1:lp
    coords = savepoints(i,:);
    
    % Assign correct label based on coordinates
    if abs(coords(1) - 2.08) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 2.5) < 0.01
        label = 'P1 (2.08,0.5,2.5)';
    elseif abs(coords(1) - 1.25) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 2.5) < 0.01
        label = 'P2 (1.25,0.5,2.5)';
    elseif abs(coords(1) - 0.417) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 2.5) < 0.01
        label = 'P3 (0.42,0.5,2.5)';    
    elseif abs(coords(1) - 2.08) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 5.0) < 0.01
        label = 'P4 (2.08,0.5,5.0)';
    elseif abs(coords(1) - 1.25) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 5.0) < 0.01
        label = 'P5 (1.25,0.5,5.0)';
    elseif abs(coords(1) - 0.417) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 5.0) < 0.01
        label = 'P6 (0.42,0.5,5.0)';
    elseif abs(coords(1) - 2.08) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 7.5) < 0.01
        label = 'P7 (2.08,0.5,7.5)';
    elseif abs(coords(1) - 1.25) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 7.5) < 0.01
        label = 'P8 (1.25,0.5,7.5)';
    elseif abs(coords(1) - 0.417) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 7.5) < 0.01
        label = 'P9 (0.42,0.5,7.5)';
    elseif abs(coords(1) - 2.08) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 10.0) < 0.01
        label = 'P10 (2.08,0.5,10.0)';
    elseif abs(coords(1) - 1.25) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 10.0) < 0.01
        label = 'P11 (1.25,0.5,10.0)';
    elseif abs(coords(1) - 0.417) < 0.01 && abs(coords(2) - 0.5) < 0.01 && abs(coords(3) - 10.0) < 0.01
        label = 'P12 (0.42,0.5,10.0)';
    else
        label = sprintf('Unknown (%.2f,%.2f,%.2f)', coords(1), coords(2), coords(3));
    end
    
    probe_info{i} = label;
end

% Initialize arrays
[all_data, time] = deal(cell(lp,1));

% Log data loading process
fprintf(log_file, 'LOADING DATA\n');
fprintf(log_file, 'Loading and preprocessing data from %d probe locations...\n\n', lp);

% Load and preprocess data for each probe location
for i = 1:lp
    % Load file with probe coordinates
    filename = sprintf('Flow0_%1.2e_%1.2e_%1.2e_dt_%g.dat',...
        savepoints(i,1), savepoints(i,2), savepoints(i,3), dt);
    f = load(fullfile(path_input, filename));

    % Data cleanup: remove unused columns and duplicates
    f(:,2:5) = [];  % Remove unnecessary velocity components
    f(:,5:7) = [];  % Remove additional unused data
    
    % Remove duplicate timestamps for clean time series
    [~, idx] = unique(f(:,1), 'stable');
    f = f(idx,:);
    
    % Store processed data
    all_data{i} = f;
    time{i} = f(:,1);
    
    fprintf(log_file, 'Loaded probe %d: %s\n', i, filename);
end

% Find common time length for all probes
min_length = min(cellfun(@(x) size(x,1), all_data));
fprintf(log_file, '\nCommon time series length: %d samples\n\n', min_length);

%% Velocity Fluctuations and Mean Velocity Calculation
fprintf(log_file, 'VELOCITY ANALYSIS\n');
fprintf(log_file, 'Calculating %s fluctuations and mean velocities...\n\n', upper(component));

% Initialize arrays
fluctuations = zeros(min_length, lp);
mean_velocities = zeros(lp, 4); % [u_mean, v_mean, w_mean, vmag_mean]
integral_length_scales = zeros(lp, 1);

% Process each probe
for i = 1:lp
    data = all_data{i}(1:min_length,:);
    
    % Calculate mean of each velocity component
    u_mean = mean(data(:,2));
    v_mean = mean(data(:,3));
    w_mean = mean(data(:,4));
    
    % Calculate mean velocity magnitude
    vmag_mean = mean(sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2));
    
    % Store mean velocities
    mean_velocities(i,:) = [u_mean, v_mean, w_mean, vmag_mean];
    
    % Select and process velocity component for fluctuations
    switch lower(component)
        case 'u'
            vel = data(:,2);
        case 'v'
            vel = data(:,3);
        case 'w'
            vel = data(:,4);
        case 'vmag'
            % Calculate velocity magnitude
            vel = sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2);
        case 'tke'
            % Calculate turbulent kinetic energy
            u_fluc = data(:,2) - u_mean;
            v_fluc = data(:,3) - v_mean;
            w_fluc = data(:,4) - w_mean;
            vel = 0.5 * (u_fluc.^2 + v_fluc.^2 + w_fluc.^2);
        otherwise
            error('Invalid component. Use ''u'', ''v'', ''w'', ''vmag'', or ''tke''');
    end
    
    % Remove mean and store fluctuations
    fluctuations(:,i) = vel - mean(vel);
end

% Create probe order and sort index for consistent ordering
probe_order = zeros(lp,1);

% Extract probe numbers for sorting
for i = 1:lp
    if ischar(probe_info{i}) || isstring(probe_info{i})
        if contains(probe_info{i}, 'P')
            probe_num = sscanf(probe_info{i}(2:end), '%d');
            probe_order(i) = probe_num;
        else
            probe_order(i) = inf;
        end
    else
        probe_info{i} = sprintf('Unknown %d', i);
        probe_order(i) = inf;
    end
end

% Create sorted index
[~, sort_idx] = sort(probe_order);

% Write mean velocities to log file in sorted order
fprintf(log_file, '---------- Mean Velocity Summary ----------\n');
fprintf(log_file, 'Probe #\t|\tLocation\t|\tu_mean\t|\tv_mean\t|\tw_mean\t|\tMagnitude\n');
fprintf(log_file, '--------------------------------------------------------------------------------\n');

for i = 1:lp
    idx = sort_idx(i);
    loc_info = probe_info{idx};
    coords_str = regexprep(loc_info, '.*\((.*)\)', '$1');
    
    fprintf(log_file, 'P%d\t|\t%s\t|\t%.3f\t|\t%.3f\t|\t%.3f\t|\t%.3f\n', ...
        i, coords_str, mean_velocities(idx,1), mean_velocities(idx,2), mean_velocities(idx,3), mean_velocities(idx,4));
end
fprintf(log_file, '--------------------------------------------------------------------------------\n\n');

%% Check Taylor's Frozen Turbulence Hypothesis Validity
fprintf(log_file, 'TAYLOR HYPOTHESIS VERIFICATION\n');
fprintf(log_file, 'Checking validity of Taylor''s frozen turbulence hypothesis (u''_rms/u ≤ 0.3)...\n\n');

% Calculate ratio for each probe
taylor_ratio = zeros(lp, 1);
taylor_validity = false(lp, 1);

fprintf(log_file, '---------- Taylor Hypothesis Validity Check ----------\n');
fprintf(log_file, 'Probe #\t|\tLocation\t|\tMean Vel.\t|\tRMS Fluct.\t|\tRatio u''_rms/u\t|\tValid?\n');
fprintf(log_file, '--------------------------------------------------------------------------------------------------\n');

for i = 1:lp
    idx = sort_idx(i);
    
    % Get mean velocity for the component being analyzed
    switch lower(component)
        case 'u'
            u_mean = mean_velocities(idx,1);
        case 'v'
            u_mean = mean_velocities(idx,2);
        case 'w'
            u_mean = mean_velocities(idx,3);
        otherwise
            u_mean = mean_velocities(idx,4);  % Use magnitude for other components
    end
    
    % Calculate RMS of fluctuations
    u_rms = sqrt(mean(fluctuations(:,idx).^2));
    
    % Calculate ratio
    ratio = u_rms / abs(u_mean);
    taylor_ratio(idx) = ratio;
    
    % Check if ratio is below threshold
    is_valid = ratio <= 0.3;
    taylor_validity(idx) = is_valid;
    
    % Extract location info for display
    loc_info = probe_info{idx};
    coords_str = regexprep(loc_info, '.*\((.*)\)', '$1');
    
    % Display results
    % Determine validity string
    if is_valid
        validity_str = 'Yes';
    else
        validity_str = 'No';
    end
    
    fprintf(log_file, 'P%d\t|\t%s\t|\t%.3f\t|\t%.3f\t|\t%.3f\t\t|\t%s\n', ...
        i, coords_str, u_mean, u_rms, ratio, validity_str);
end
fprintf(log_file, '--------------------------------------------------------------------------------------------------\n\n');

%% Alternative Analysis for Non-Taylor Points Using Cross-Correlation
fprintf(log_file, '\nALTERNATIVE ANALYSIS FOR NON-TAYLOR POINTS\n');
fprintf(log_file, 'Applying cross-correlation method for probes where Taylor''s hypothesis is invalid...\n\n');

% Identify probes that don't satisfy Taylor's hypothesis
non_taylor_probes = find(~taylor_validity);

if isempty(non_taylor_probes)
    fprintf(log_file, 'All probes satisfy Taylor''s hypothesis. No alternative analysis needed.\n\n');
else
    fprintf(log_file, 'Probes requiring alternative analysis: %s\n\n', mat2str(non_taylor_probes));
    
    % Create a spatial map of probe locations
    probe_locations = zeros(lp, 3);
    for i = 1:lp
        probe_locations(i,:) = savepoints(i,:);
    end
    
    % For each non-Taylor probe, find nearest neighbors
    for i = 1:length(non_taylor_probes)
        probe_idx = non_taylor_probes(i);
        current_loc = probe_locations(probe_idx,:);
        
        % Calculate distance to all other probes
        distances = sqrt(sum((probe_locations - repmat(current_loc, lp, 1)).^2, 2));
        distances(probe_idx) = inf;  % Exclude self
        
        % Find nearest probes (up to 3)
        [sorted_dist, sorted_idx] = sort(distances);
        nearest_probes = sorted_idx(1:min(3, lp-1));
        
        fprintf(log_file, 'Probe %d: Nearest probes are %s (distances: %s meters)\n', ...
            probe_idx, mat2str(nearest_probes), mat2str(sorted_dist(1:length(nearest_probes)), 3));
        
        % Get fluctuation data for current probe
        u_current = fluctuations(:, probe_idx);
        
        % Initialize results storage
        cross_corr_results = zeros(length(nearest_probes), 3);  % [max_corr, lag_time, spatial_scale]
        
        % Process each nearest neighbor
        for j = 1:length(nearest_probes)
            neighbor_idx = nearest_probes(j);
            u_neighbor = fluctuations(:, neighbor_idx);
            
            % Calculate cross-correlation
            [xcorr_result, lags] = xcorr(u_current, u_neighbor, 'coeff');
            
            % Find peak correlation and corresponding lag
            [max_corr, max_idx] = max(xcorr_result);
            lag_samples = lags(max_idx);
            lag_time = lag_samples * dt;
            
            % Calculate spatial distance between probes
            spatial_dist = sorted_dist(j);
            
            % If lag_time is non-zero, calculate effective advection velocity and spatial scale
            if lag_time ~= 0
                effective_velocity = spatial_dist / abs(lag_time);
                
                % Calculate integral time scale using autocorrelation
                [acf, acf_lags] = xcorr(u_current, 'coeff');
                acf = acf(acf_lags >= 0);  % Only positive lags
                acf_lags = acf_lags(acf_lags >= 0);
                
                % Find first zero crossing or use decay to 1/e
                zero_idx = find(acf <= 1/exp(1), 1, 'first');
                if isempty(zero_idx)
                    zero_idx = find(acf <= 0, 1, 'first');
                    if isempty(zero_idx)
                        zero_idx = length(acf);
                    end
                end
                
                % Integral time scale
                integral_time_scale = trapz(acf_lags(1:zero_idx) * dt, acf(1:zero_idx));
                
                % Spatial scale using effective velocity
                spatial_scale = effective_velocity * integral_time_scale;
            else
                effective_velocity = NaN;
                spatial_scale = NaN;
            end
            
            % Store results
            cross_corr_results(j, :) = [max_corr, lag_time, spatial_scale];
            
            fprintf(log_file, '  - With Probe %d: Corr=%.3f, Lag=%.4fs, Spatial Scale=%.4fm\n', ...
                neighbor_idx, max_corr, lag_time, spatial_scale);
        end
        
        % Determine best estimate based on highest correlation
        [best_corr, best_idx] = max(cross_corr_results(:, 1));
        best_scale = cross_corr_results(best_idx, 3);
        
        if ~isnan(best_scale) && best_corr > 0.3  % Only use if correlation is reasonable
            fprintf(log_file, '  => Estimated length scale for Probe %d: %.4f m (from cross-correlation)\n', ...
                probe_idx, best_scale);
            
            % Update integral length scale for this probe
            integral_length_scales(probe_idx) = best_scale;
            fprintf(log_file, '     (Previous value from file: %.4f m)\n', integral_length_scales(probe_idx));
        else
            fprintf(log_file, '  => Cannot estimate reliable length scale from cross-correlation\n');
            fprintf(log_file, '     (Using original value: %.4f m)\n', integral_length_scales(probe_idx));
        end
    end
    
    fprintf(log_file, '\n---------- Updated Integral Length Scales ----------\n');
    fprintf(log_file, 'Probe #\t|\tLocation\t\t|\tLength Scale (m)\t|\tMethod\n');
    fprintf(log_file, '----------------------------------------------------------------------------\n');
    
    % Display in ordered form with method indication
    for i = 1:length(sort_idx)
        idx = sort_idx(i);
        loc_info = probe_info{idx};
        coords_str = regexprep(loc_info, '.*\((.*)\)', '$1'); 
        if taylor_validity(idx)
            method = 'Taylor';
        else
            method = 'Cross-correlation';
        end
        
        fprintf(log_file, 'P%d\t|\t%s\t|\t%.4f\t\t|\t%s\n', ...
            i, coords_str, integral_length_scales(idx), method);
    end
    fprintf(log_file, '----------------------------------------------------------------------------\n\n');
end

%% Read Integral Length Scales from File
fprintf(log_file, 'INTEGRAL LENGTH SCALES\n');
fprintf(log_file, 'Reading integral length scales from file...\n');

% Look for integral scales files
integral_scale_files = dir(fullfile(path_output, 'IntegralScales_*.txt'));

if isempty(integral_scale_files)
    fprintf(log_file, 'Error: No integral scales files found in %s\n', path_input);
    error('No integral scales files found in %s', path_input);
end

% Initialize array for integral length scales
integral_length_scales = zeros(lp, 1);
found_scales = false(lp, 1);

% Process each integral scale file
for file_idx = 1:length(integral_scale_files)
    file_path = fullfile(integral_scale_files(file_idx).folder, integral_scale_files(file_idx).name);
    fprintf(log_file, 'Processing file: %s\n', integral_scale_files(file_idx).name);
    
    try
        fileData = readtable(file_path, 'ReadVariableNames', true);
        
        if ~any(strcmp('L', fileData.Properties.VariableNames))
            fprintf(log_file, 'Warning: File %s does not contain column "L", skipping...\n', ...
                integral_scale_files(file_idx).name);
            continue;
        end
        
        % Process each entry in the file
        for i = 1:height(fileData)
            if any(strcmp('Point', fileData.Properties.VariableNames))
                if iscell(fileData.Point)
                    pointName = fileData.Point{i};
                else
                    pointName = char(fileData.Point(i));
                end
                
                % Extract coordinates from point name
                parts = strsplit(pointName, '_');
                if length(parts) >= 4
                    x = str2double(parts{2});
                    y = str2double(parts{3});
                    z = str2double(parts{4});
                    
                    % Find matching probe
                    for probe_idx = 1:lp
                        coords = savepoints(probe_idx,:);
                        if abs(coords(1) - x) < 0.01 && abs(coords(2) - y) < 0.01 && abs(coords(3) - z) < 0.01
                            integral_length_scales(probe_idx) = fileData.L(i);
                            found_scales(probe_idx) = true;
                            break;
                        end
                    end
                end
            end
        end
    catch ME
        fprintf(log_file, 'Warning: Error reading file %s: %s\n', ...
            integral_scale_files(file_idx).name, ME.message);
    end
end

% Check if we found scales for all probes
if ~all(found_scales)
    missing_probes = find(~found_scales);
    fprintf(log_file, 'Warning: Missing integral scales for probes: %s\n', mat2str(missing_probes));
    
    % Fill missing values with mean
    if any(found_scales)
        mean_scale = mean(integral_length_scales(found_scales));
        integral_length_scales(~found_scales) = mean_scale;
        fprintf(log_file, 'Using mean value (%.4f m) for missing probes\n', mean_scale);
    else
        integral_length_scales(:) = 0.1;  % Default value
        fprintf(log_file, 'No scales found. Using default value (0.1 m) for all probes\n');
    end
end

% Display organized length scales
fprintf(log_file, '\n---------- Integral Length Scale Verification ----------\n');
fprintf(log_file, 'Probe #\t|\tLocation\t\t|\tLength Scale (m)\n');
fprintf(log_file, '------------------------------------------------------\n');

% Display in ordered form (P1 to P12)
for i = 1:length(sort_idx)
    idx = sort_idx(i);
    loc_info = probe_info{idx};
    coords_str = regexprep(loc_info, '.*\((.*)\)', '$1'); % Extract just the coordinates
    fprintf(log_file, 'P%d\t|\t%s\t|\t%.4f\n', i, coords_str, integral_length_scales(idx));
end
fprintf(log_file, '------------------------------------------------------\n\n');

% %% Integral Length Scale Calculation using MATLAB funstions
% fprintf('Calculating integral length scale for each probe...\n')
% integral_length_scale = zeros(lp,1);
% 
% for i = 1:lp
%     % Use the velocity fluctuations for the selected component
%     u_fluc = fluctuations(:,i);
% 
%     % Compute autocorrelation (biased)
%     [acf, lags] = xcorr(u_fluc, 'biased');
%     mid = ceil(length(acf)/2);
%     acf = acf(mid:end); % Only positive lags
%     lags = lags(mid:end);
% 
%     % Normalize autocorrelation
%     acf = acf / acf(1);
% 
%     % Find first zero crossing
%     zero_idx = find(acf <= 0, 1, 'first');
%     if isempty(zero_idx)
%         zero_idx = length(acf);
%     end
% 
%     % Integrate autocorrelation up to first zero crossing
%     integral_length_scale(i) = trapz(lags(1:zero_idx)*dt, acf(1:zero_idx));
% end
% 
% % Display results
% for i = 1:lp
%     fprintf('Probe %d: Integral time scale = %.4f s\n', i, integral_length_scale(i));
% end
% 
% % Add after displaying the results
% fprintf('\nConverting to integral length scales...\n')
% for i = 1:lp
%     data = all_data{i}(1:min_length,:);
%     % mean_vel = mean(sqrt(data(:,2).^2 + data(:,3).^2 + data(:,4).^2)); % Mean velocity magnitude
%     mean_vel=mean(data(:,4));
%     length_scale = mean_vel * integral_length_scale(i);
%     fprintf('Probe %d: Integral length scale = %.4f m\n', i, length_scale);
% end

%% Calculate Fish and Integral Scale Frequencies
fprintf(log_file, 'CHARACTERISTIC FREQUENCY ANALYSIS\n');

% Fish frequency calculation (Umag/Lf)
fish_frequencies = zeros(lp, 1);
fprintf(log_file, '\n---------- Fish Frequency Analysis ----------\n');
fprintf(log_file, 'Probe #\t|\tLocation\t|\tMean Vel. (m/s)\t|\tFrequency (Hz)\n');
fprintf(log_file, '------------------------------------------------------------------\n');

% Calculate fish frequency for each probe
for i = 1:lp
    idx = sort_idx(i);
    velocity_mag = mean_velocities(idx,4); % Using velocity magnitude
    freq = velocity_mag / Lf;
    fish_frequencies(idx) = freq;
    
    loc_info = probe_info{idx};
    coords_str = regexprep(loc_info, '.*\((.*)\)', '$1');
    
    fprintf(log_file, 'P%d\t|\t%s\t|\t%.3f\t\t|\t%.2f\n', ...
        i, coords_str, velocity_mag, freq);
end
fprintf(log_file, '------------------------------------------------------------------\n\n');

% Integral scale frequency calculation (Umag/L)
integral_scale_frequencies = zeros(lp, 1);

fprintf(log_file, '---------- Integral Scale Frequency Analysis ----------\n');
fprintf(log_file, 'Probe #\t|\tLocation\t|\tUmag (m/s)\t|\tL (m)\t|\tFrequency (Hz)\n');
fprintf(log_file, '--------------------------------------------------------------------------------\n');

% Calculate integral scale frequency for each probe
for i = 1:lp
    idx = sort_idx(i);
    velocity_mag = mean_velocities(idx,4); % Using velocity magnitude
    length_scale = integral_length_scales(idx);
    
    if length_scale > 0
        freq = velocity_mag / length_scale;
    else
        freq = NaN; % Avoid division by zero
    end
    
    integral_scale_frequencies(idx) = freq;
    
    loc_info = probe_info{idx};
    coords_str = regexprep(loc_info, '.*\((.*)\)', '$1');
    
    fprintf(log_file, 'P%d\t|\t%s\t|\t%.3f\t|\t%.4f\t|\t%.2f\n', ...
        i, coords_str, velocity_mag, length_scale, freq);
end
fprintf(log_file, '--------------------------------------------------------------------------------\n\n');

%% PSD Calculation
fprintf(log_file, 'POWER SPECTRAL DENSITY CALCULATION\n');
fprintf(log_file, 'Calculating Power Spectral Density using %s method...\n\n', upper(psd_method));

fs = 1/dt;  % Sampling frequency
[f_psd, psd] = deal(zeros(nfft/2+1, lp));

% Configure PSD calculation based on selected method
switch lower(psd_method)
    case 'welch'
        window = hann(window_size);
        noverlap = round(window_size*overlap_ratio);
        
        % Calculate PSD for each probe location
        for i = 1:lp
            [psd(:,i), f_psd(:,i)] = pwelch(fluctuations(:,i), window, noverlap, nfft, fs);
        end
        
    case 'periodogram'
        window = hann(min_length);
        
        % Calculate PSD for each probe location
        for i = 1:lp
            [psd(:,i), f_psd(:,i)] = periodogram(fluctuations(:,i), window, nfft, fs, 'onesided');
        end
        
    case 'bartlett'
        window = rectwin(window_size);
        noverlap = 0;
        
        % Calculate PSD for each probe location
        for i = 1:lp
            [psd(:,i), f_psd(:,i)] = pwelch(fluctuations(:,i), window, noverlap, nfft, fs);
        end
        
    otherwise
        fprintf(log_file, 'Error: Invalid PSD method. Use ''welch'', ''periodogram'', or ''bartlett''\n');
        error('Invalid PSD method. Use ''welch'', ''periodogram'', or ''bartlett''');
end

%% Exclude data after cutoff frequency
indices = f_psd(:,1) <= cutoff_frequency;
f_psd = f_psd(indices, :);
psd = psd(indices, :);

fprintf(log_file, 'PSD calculation complete. Using frequencies up to %d Hz.\n\n', cutoff_frequency);

%% Energy Analysis for Fish and Integral Scale Frequencies
fprintf(log_file, 'ENERGY ANALYSIS\n');
fprintf(log_file, 'Analyzing energy content at characteristic frequencies...\n\n');

% Initialize arrays
fish_scale_energy = zeros(lp,1);
fish_scale_percentage = zeros(lp,1);
integral_scale_energy = zeros(lp,1);
integral_scale_percentage = zeros(lp,1);
total_energy = zeros(lp,1);

% Fish frequency energy analysis
fprintf(log_file, '---------- Fish Scale Energy Analysis ----------\n');
fprintf(log_file, 'Probe #\t|\tFish Freq.\t|\tBand Range\t\t|\tEnergy\t\t|\tPercentage\n');
fprintf(log_file, '----------------------------------------------------------------------------------\n');

for i = 1:lp
    idx = sort_idx(i);
    fish_freq = fish_frequencies(idx);
    
    % Define frequency band around fish frequency
    freq_min = fish_freq * (1 - band_width_factor);
    freq_max = fish_freq * (1 + band_width_factor);
    
    % Find indices of frequency band in PSD data
    band_indices = f_psd(:,1) >= freq_min & f_psd(:,1) <= freq_max;
    
    % Check if any frequencies are in the band
    if any(band_indices)
        band_energy = trapz(f_psd(band_indices,1), psd(band_indices,idx));
    else
        band_energy = 0;
        fprintf(log_file, 'Warning: No frequencies found in fish band %.2f-%.2f Hz for probe %d\n', ...
            freq_min, freq_max, i);
    end
    
    % Calculate total energy across all frequencies
    total = trapz(f_psd(:,1), psd(:,idx));
    total_energy(idx) = total;
    
    % Calculate percentage
    if total > 0
        percentage = (band_energy / total) * 100;
    else
        percentage = 0;
    end
    
    fish_scale_energy(idx) = band_energy;
    fish_scale_percentage(idx) = percentage;
    
    fprintf(log_file, 'P%d\t|\t%.2f Hz\t|\t%.2f - %.2f Hz\t|\t%.3e\t|\t%.2f%%\n', ...
        i, fish_freq, freq_min, freq_max, band_energy, percentage);
end
fprintf(log_file, '----------------------------------------------------------------------------------\n\n');

% Integral scale energy analysis
fprintf(log_file, '---------- Integral Scale Energy Analysis ----------\n');
fprintf(log_file, 'Probe #\t|\tInt. Freq.\t|\tBand Range\t\t|\tEnergy\t\t|\tPercentage\n');
fprintf(log_file, '----------------------------------------------------------------------------------\n');

for i = 1:lp
    idx = sort_idx(i);
    int_freq = integral_scale_frequencies(idx);
    
    % Skip if NaN
    if isnan(int_freq)
        integral_scale_energy(idx) = NaN;
        integral_scale_percentage(idx) = NaN;
        fprintf(log_file, 'P%d\t|\tN/A\t\t|\tN/A\t\t\t|\tN/A\t\t|\tN/A\n', i);
        continue;
    end
    
    % Define frequency band around integral scale frequency
    freq_min = int_freq * (1 - band_width_factor);
    freq_max = int_freq * (1 + band_width_factor);
    
    % Find indices of frequency band in PSD data
    band_indices = f_psd(:,1) >= freq_min & f_psd(:,1) <= freq_max;
    
    % Check if any frequencies are in the band
    if any(band_indices)
        band_energy = trapz(f_psd(band_indices,1), psd(band_indices,idx));
    else
        band_energy = 0;
        fprintf(log_file, 'Warning: No frequencies found in integral band %.2f-%.2f Hz for probe %d\n', ...
            freq_min, freq_max, i);
    end
    
    % Calculate percentage
    if total_energy(idx) > 0
        percentage = (band_energy / total_energy(idx)) * 100;
    else
        percentage = 0;
    end
    
    integral_scale_energy(idx) = band_energy;
    integral_scale_percentage(idx) = percentage;
    
    fprintf(log_file, 'P%d\t|\t%.2f Hz\t|\t%.2f - %.2f Hz\t|\t%.3e\t|\t%.2f%%\n', ...
        i, int_freq, freq_min, freq_max, band_energy, percentage);
end
fprintf(log_file, '----------------------------------------------------------------------------------\n\n');

%% Final Frequency and Energy Comparison
fprintf(log_file, '---------- Frequency and Energy Comparison ----------\n');
fprintf(log_file, 'Probe #\t|\tFish Freq.\t|\tInt. Freq.\t|\tRatio\t|\tFish Energy %%\t|\tInt. Energy %%\t|\tRatio\n');
fprintf(log_file, '-----------------------------------------------------------------------------------------------------------\n');

for i = 1:lp
    idx = sort_idx(i);
    fish_freq = fish_frequencies(idx);
    int_freq = integral_scale_frequencies(idx);
    
    if ~isnan(int_freq) && int_freq > 0
        freq_ratio = fish_freq / int_freq;
    else
        freq_ratio = NaN;
    end
    
    fish_energy = fish_scale_percentage(idx);
    int_energy = integral_scale_percentage(idx);
    
    if ~isnan(int_energy) && int_energy > 0
        energy_ratio = fish_energy / int_energy;
    else
        energy_ratio = NaN;
    end
    
    fprintf(log_file, 'P%d\t|\t%.2f Hz\t|\t%.2f Hz\t|\t%.2f\t|\t%.2f%%\t\t|\t%.2f%%\t\t|\t%.2f\n', ...
        i, fish_freq, int_freq, freq_ratio, fish_energy, int_energy, energy_ratio);
end
fprintf(log_file, '-----------------------------------------------------------------------------------------------------------\n\n');

%% Set Default Font
set(0, 'defaultAxesFontName', 'Times New Roman');
set(0, 'defaultTextFontName', 'Times New Roman');
set(0, 'defaultLegendFontName', 'Times New Roman');
set(0, 'defaultColorbarFontName', 'Times New Roman');

%% Visualization
% Main PSD plot
fprintf(log_file, 'VISUALIZATIONS\n');
fprintf(log_file, 'Generating PSD and analysis plots...\n');

% Create new figure with adequate size
fig = figure('Position', [100 100 950 700]);
hold on;

% Define distinct visual styles for each probe
line_styles = {'-', '--', ':', '-.', '-', '--', ':', '-.', '-', '--', ':', '-.'};
marker_styles = {'none', 'none', 'none', 'none', 'o', 'o', 'o', 'o', 's', 's', 's', 's'};
marker_spacing = [0, 0, 0, 0, 50, 50, 50, 50, 50, 50, 50, 50];

% Define distinct color palette with high contrast
colors = [
    0.0000, 0.4470, 0.7410;  % blue
    0.8500, 0.3250, 0.0980;  % red
    0.4660, 0.6740, 0.1880;  % green
    0.9290, 0.6940, 0.1250;  % yellow
    0.3010, 0.7450, 0.9330;  % light blue
    0.6350, 0.0780, 0.1840;  % burgundy
    0.4940, 0.1840, 0.5560;  % purple
    0.0000, 0.0000, 0.0000;  % black
    0.8000, 0.8000, 0.0000;  % olive
    0.0000, 0.5000, 0.0000;  % dark green
    1.0000, 0.0000, 1.0000;  % magenta
    0.5000, 0.5000, 0.5000;  % gray
];

% Plot each probe with distinct style
h_lines = gobjects(lp,1);

% Use consistent probe order
for i = 1:lp
    idx = sort_idx(i);
    style_idx = mod(i-1, length(line_styles)) + 1;
    color_idx = mod(i-1, size(colors,1)) + 1;
    marker = marker_styles{style_idx};
    
    if strcmp(marker, 'none')
        % No markers
        h_lines(i) = loglog(f_psd(:,1), psd(:,idx), line_styles{style_idx}, ...
            'LineWidth', 2, 'Color', colors(color_idx,:));
    else
        % Add markers with specified spacing
        spacing = marker_spacing(style_idx);
        h_lines(i) = loglog(f_psd(:,1), psd(:,idx), line_styles{style_idx}, ...
            'LineWidth', 2, 'Color', colors(color_idx,:));
        hold on;
        loglog(f_psd(1:spacing:end,1), psd(1:spacing:end,idx), marker, ...
            'MarkerSize', 6, 'Color', colors(color_idx,:), 'MarkerFaceColor', 'white');
    end
end

% Force logarithmic scale on both axes
set(gca, 'XScale', 'log', 'YScale', 'log');

xlim([0.5, cutoff_frequency]);

% Enhance grid and ticks for better readability
grid on;
set(gca, 'GridLineStyle', ':');
set(gca, 'FontSize', 12);

% Set appropriate y-axis label based on component
ylabel_text = 'PSD (m^2/s^2/Hz)';  % Default value
switch lower(component)
    case 'tke'
        ylabel_text = 'PSD (m^2/s^2/Hz)';
    case 'vmag'
        ylabel_text = 'PSD (m^2/s^2/Hz)';
    case 'u'
        ylabel_text = 'PSD (m^2/s^2/Hz)';
    case 'v'
        ylabel_text = 'PSD (m^2/s^2/Hz)';
    case 'w'
        ylabel_text = 'PSD (m^2/s^2/Hz)';
end

ylabel(ylabel_text, 'FontSize', 12, 'FontWeight', 'bold')
% title(['PSD Analysis - ', upper(component), ' Component'], 'FontSize', 14)
grid on

% Sort labels by probe number
sorted_labels = probe_info(sort_idx);

% Add Kolmogorov -5/3 reference line
hold on
f_ref = f_psd(f_psd(:,1) > 1 & f_psd(:,1) < 10, 1);  % Select inertial subrange
ref_line = 0.1*(f_ref).^(-5/3);                       % Reference slope
kolm_line = loglog(f_ref, ref_line, 'k-', 'LineWidth', 1.5);

% Calculate average fish frequency and add vertical line
mean_fish_freq = mean(fish_frequencies(~isnan(fish_frequencies)));
min_psd = min(min(psd));
max_psd = max(max(psd));
fish_line = loglog([mean_fish_freq mean_fish_freq], [min_psd max_psd], 'r--', 'LineWidth', 2);

% Add band indicators for fish frequency range
lower_freq = mean_fish_freq * (1 - band_width_factor);
upper_freq = mean_fish_freq * (1 + band_width_factor);
light_red = [1 0.8 0.8];
fish_band = fill([lower_freq upper_freq upper_freq lower_freq], [min_psd min_psd max_psd max_psd], light_red, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
uistack(fish_band, 'bottom');  % Put the band behind other elements

simplified_labels = cell(size(sorted_labels));
for i = 1:length(sorted_labels)
    % Extract just the probe number (P1, P2, etc.) from the full label
    if contains(sorted_labels{i}, 'P')
        % Extract text before the space or opening parenthesis
        parts = strsplit(sorted_labels{i}, ' ');
        simplified_labels{i} = parts{1};  % Take just "P1", "P2", etc.
    else
        simplified_labels{i} = sorted_labels{i};  % Keep original if not a probe
    end
end

% Add legend with all reference elements
%legend([h_lines; kolm_line; fish_line], [simplified_labels; {'-5/3 slope'; 'Fish Scale Frequency'}], 'Location', 'southwest', 'Interpreter', 'none')
hold off

% Create comparison bar chart of percentage energy
fig_compare = figure('Position', [150 150 900 500]);

% Combine data for comparison
comparison_data = [fish_scale_percentage(sort_idx), integral_scale_percentage(sort_idx)];

% Create grouped bar chart
h = bar(1:lp, comparison_data);
set(h(1), 'FaceColor', [0.2 0.6 1])
set(h(2), 'FaceColor', [1 0.6 0.2])

xlabel('Probe Number', 'FontSize', 12);
ylabel('Energy Percentage (%)', 'FontSize', 12);
% title('Comparison of Energy Content: Fish vs. Integral Scales', 'FontSize', 14);
legend('Fish Scale (L_f = 6 cm)', 'Integral Length Scale', 'Location', 'best');
xticks(1:lp);
xticklabels(arrayfun(@(i) ['P' num2str(i)], 1:lp, 'UniformOutput', false));
grid on;

%% Save Results
saveas(fig, fullfile(output_dir, ['PSD_', upper(component), '_', upper(psd_method), '.png']));
saveas(fig_compare, fullfile(output_dir, 'Fish_vs_Integral_Energy_Comparison.png'));

% Save analysis results and parameters
save(fullfile(output_dir, ['PSD_data_', lower(psd_method), '.mat']), 'f_psd', 'psd', 'component',...
    'psd_method', 'window_size', 'overlap_ratio', 'nfft', 'dt', 'integral_length_scales',...
    'mean_velocities', 'fish_frequencies', 'integral_scale_frequencies', 'Lf',...
    'fish_scale_energy', 'fish_scale_percentage', 'total_energy',...
    'integral_scale_energy', 'integral_scale_percentage', 'taylor_ratio', 'taylor_validity');

% Record runtime
runtime = toc;
fprintf(log_file, '\nSAVED RESULTS\n');
fprintf(log_file, 'Analysis saved to directory: %s\n', output_dir);
fprintf(log_file, 'PSD plot saved as: %s\n', ['PSD_', upper(component), '_', upper(psd_method), '.png']);
fprintf(log_file, 'Data file saved as: %s\n', ['PSD_data_', lower(psd_method), '.mat']);
fprintf(log_file, 'Log file saved as: %s\n', log_filename);
fprintf(log_file, 'Script runtime: %.2f seconds\n', runtime);

% Close the log file
fclose(log_file);

% Print only minimal info to console
fprintf('PSD Analysis complete for case %s, component %s\n', case_name, upper(component));
fprintf('Results and detailed log saved to %s\n', output_dir);
fprintf('Runtime: %.2f seconds\n', runtime);