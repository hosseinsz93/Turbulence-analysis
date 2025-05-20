%% Turbulence Analysis Script
% This script analyzes turbulent flow data to calculate:
% - Reynolds stresses
% - Turbulent Kinetic Energy (TKE)
% - Autocorrelation Functions (ACF)
% - Integral Time and Length Scales
%
% Author: Hossein Seyedzadeh
% Email: hossein.seyyedzadeh@stonybrook.edu
% Last modified: March 30, 2025

%% Initialization
clc;
clear;
close all;
tic;  % Start measuring runtime

%% Parameters and Data Loading
% Load measurement points coordinates from file
savepoints = load('savepoints');
dt = 0.0005;  % Time step [s]
lp = length(savepoints);  % Number of measurement points
% Select case to analyze
case_name = 'MR';  % Options: 'HR', 'HS', 'MR', 'MS'
% Set input/output paths based on selected case
path_input = fullfile("Input", case_name, "");  % Input directory with case subfolder
path_output = fullfile("Output", case_name, "");  % Output directory with case subfolder

% Create output directory if it doesn't exist
if ~exist(path_output, 'dir')
    mkdir(path_output);
    fprintf('Created output directory: %s\n', path_output);
end

fprintf('Processing case: %s\n', case_name);
fprintf('Input path: %s\n', path_input);
fprintf('Output path: %s\n', path_output);

%% Main Processing Loop
for i = 1:lp
    %% File Loading and Preprocessing
    % Generate filename based on coordinates and load data
    filename = sprintf('Flow0_%1.2e_%1.2e_%1.2e_dt_%g.dat', ...
        savepoints(i,1), savepoints(i,2), savepoints(i,3), dt);
    full_filename = fullfile(path_input, filename);
    f = load(full_filename);
    
    % Initialize arrays on first iteration
    if (i == 1)
        Tl = zeros(1,lp);  % Integral time scale
        L = zeros(1,lp);   % Integral length scale
        fileLabels = cell(1,lp);  % Store filenames for labeling
    end
    
    % Create label for current point
    fileLabels{i} = sprintf('Flow0_%1.2e_%1.2e_%1.2e', ...
        savepoints(i,1), savepoints(i,2), savepoints(i,3));
    
    % Data cleanup: remove unused columns and duplicates
    f(:,2:5) = [];  % Remove unnecessary velocity components
    f(:,5:7) = [];  % Remove additional unused data
    
    % Remove time-duplicated entries
    j = 2;
    while (j ~= length(f)+1)
        if (f(j,1) <= f(j-1,1))
            f(j,:) = [];
        else
            j = j+1;
        end
    end
    
    %% Initialize Statistical Arrays
    if (i == 1)
        % Mean velocities
        u_mean = zeros(1,lp);
        v_mean = zeros(1,lp);
        w_mean = zeros(1,lp);
        
        % Fluctuating components
        uf = zeros(length(f),lp);  % u' fluctuation
        vf = zeros(length(f),lp);  % v' fluctuation
        wf = zeros(length(f),lp);  % w' fluctuation
        
        % Reynolds stresses
        uuf = zeros(length(f),lp);  % u'u'
        vvf = zeros(length(f),lp);  % v'v'
        wwf = zeros(length(f),lp);  % w'w'
        
        % Running means
        uuf_mean = zeros(length(f),lp);
        vvf_mean = zeros(length(f),lp);
        wwf_mean = zeros(length(f),lp);
        
        % Turbulent Kinetic Energy
        TKE = zeros(length(f),lp);
    end
    
    %% Statistical Calculations
    % Calculate mean velocities
    u_mean(1,i) = mean(f(:,2));
    v_mean(1,i) = mean(f(:,3));
    w_mean(1,i) = mean(f(:,4));

    vmag(1,i) = sqrt(u_mean(1,i).^2+v_mean(1,i).^2+w_mean(1,i).^2);
    
    % Calculate velocity fluctuations
    uf(:,i) = f(:,2) - u_mean(1,i);  % u' = u - ū
    vf(:,i) = f(:,3) - v_mean(1,i);  % v' = v - v̄
    wf(:,i) = f(:,4) - w_mean(1,i);  % w' = w - w̄
    
    % Calculate Reynolds stresses
    uuf(:,i) = uf(:,i).*uf(:,i);  % u'u'
    vvf(:,i) = vf(:,i).*vf(:,i);  % v'v'
    wwf(:,i) = wf(:,i).*wf(:,i);  % w'w'
    
    % Calculate running means
    for j = 1:length(f)
        uuf_mean(j,i) = mean(uuf(1:j,i));
        vvf_mean(j,i) = mean(vvf(1:j,i));
        wwf_mean(j,i) = mean(wwf(1:j,i));
    end
    
    % Calculate Turbulent Kinetic Energy
    TKE(:,i) = (uuf_mean(:,i) + vvf_mean(:,i) + wwf_mean(:,i))/2;
    
    %% Autocorrelation Function Calculation
    lagMax = 800;  % Maximum time lag
    if (i == 1)
        ACF = zeros(lagMax,lp+1);
    end
    
    % Calculate autocorrelation using streamwise fluctuations
    s0 = mean(wf(:,i).^2);  % Variance of streamwise fluctuations
    for lag = 1:lagMax
        sk = 0;
        for k = 1:(length(wf)-lag)
            sk = sk + wf(k,i)*wf(k+lag,i);
        end
        sk = sk/(length(wf)-lag);
        ACF(lag,i+1) = sk/s0;  % Normalized autocorrelation
        if (i == 1)
            ACF(lag,1) = dt*lag;  % Time lag array
        end
    end
    
    %% Integral Scale Calculations
    % Calculate integral time scale using trapezoidal integration
    flag = 1;
    sum = 0;
    n = 1;
    while flag > 0 && n <= lagMax
        if (n == 1)
            sum = (1 + ACF(n,i+1))/2 * dt;
        else
            sum = sum + (ACF(n,i+1) + ACF(n-1,i+1))/2 * dt;
        end
        flag = ACF(n,i+1);
        n = n + 1;
    end
    Tl(1,i) = sum;  % Integral time scale
    L(1,i) = Tl(1,i) * w_mean(1,i);  % Integral length scale
    % L(1,i) = Tl(1,i) * vmag(1,i);  % Integral length scale
end

%% Post-Processing
% Add unity correlation at zero lag
ACF1 = zeros(1,lp+1);
ACF1(1,2:end) = 1;
ACF = [ACF1; ACF];
Ts = Tl*600;  % Calculate sampling time

%% Visualization
% Setup time array for plotting
time = dt:dt:dt*length(f);

% Figure 1: Reynolds Stress Evolution
fig1 = figure(1);
plot(time, uuf_mean)
legend(fileLabels, 'Location', 'best')
xlabel("Time(s)", 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
ylabel("<u'u'>", 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
xlim([0, length(f)*dt])
set(gcf, 'position', [10 10 900 550])

% Figure 2: Turbulent Kinetic Energy
fig2 = figure(2);
plot(time, TKE)
legend(fileLabels, 'Location', 'best')
xlabel("Time(s)", 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
ylabel("Turbulence Kinetic Energy", 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
xlim([0, length(f)*dt])
set(gcf, 'position', [10 10 900 550])

% Figure 3: Autocorrelation Function
fig3 = figure(3);
plot(ACF(:,1), ACF(:,2:end))
legend(fileLabels, 'Location', 'best')
xlabel("Time(s)", 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
ylabel("Autocorrelation Function", 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman')
xlim([0, (length(ACF)-1)*dt])
set(gcf, 'position', [10 10 900 550])

%% Save Results
% Save figures
saveas(fig1, fullfile(path_output, "ReS.tif"))
saveas(fig2, fullfile(path_output, "TKE.tif"))
saveas(fig3, fullfile(path_output, "ACF.tif"))

% Save numerical results to file
result_filename = fullfile(path_output, sprintf("IntegralScales_%s.txt", case_name));
fw = fopen(result_filename, "wb");
fprintf(fw, 'Point\tTl\tL\tTs\n');  % Write header
for i = 1:lp
    fprintf(fw, '%s\t%f\t%f\t%f\n', fileLabels{i}, Tl(1,i), L(1,i), Ts(1,i));
end
fclose(fw);

%% Display Results
% Display integral scales
fprintf('\nIntegral Scales Summary:\n');
fprintf('%-20s %-12s %-12s %-12s\n', 'Point', 'Tl', 'L', 'Ts');

% Print data row by row
for i = 1:lp
  fprintf('%-20s %-12.6f %-12.6f %-12.6f\n', ...
      fileLabels{i}, Tl(1,i), L(1,i), Ts(1,i));
end

% Display runtime
runtime = toc;
fprintf('\nScript runtime: %.2f seconds\n', runtime);