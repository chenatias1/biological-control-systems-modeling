clc ; clear all ; close all ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Simulation of the Cardiovascular System %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Complete the code in places where the "..." string appears %%

%% Flags:
Plot_flag = 1; % 0 = off , 1 = on
 
%% system parameters:
% Time parameters
HR           = 60+3+7 ; % [BPM] % 60 + sum of last digits from all members
dt           = 5e-4        ; % [sec]
Heart_cycles = 20          ; % total heart cycles
N_per_cycle      = round(60 / (HR * dt))     ; % number of steps per heart cycle

% Heart parameters
V0           = 15  ; % [ml] V0 for Plv calculation
Emax         = 2.0 ; % contractility
N_Systole    = round(N_per_cycle * 1/3); % number of points per ventricle Systole
N_Diastole   = round(N_per_cycle * 2/3); % number of points per ventricle Diastole
n_s = 1:N_Systole;
E_norm_systole(n_s) = 0.5 * (1 + sin((2*pi*n_s)/N_Systole - pi/2));
E_norm_systole(N_Systole+1:N_per_cycle) = 0;
E_diastole = 10/(120-V0);
E = max(E_diastole, Emax * E_norm_systole);% Heart Elasticity: combine Systole (elasticity function) & Diastole (uniform value) 
 
% Vascular constants:
Ra = 0.1;  % arterial resistance 
Rp = 1.0;  % peripheral resistance
Rv = 0.01; % venous filling resistance
Ca = 2.0;  % arterial compliance
Cv = 300.0; % venous compliance 
 
% Initiate variables:
%Volume [ml]
Vlv(1)  = 120;  % left ventricle
Va(1)   = 270;  % arteries
Vv(1)   = 2700; % veins 
%Pressure [mmHg]
Plv(1)  = 0;    % left ventricle
Pa(1)   = 70;   % arterial capacitor
Pv(1)   = 9;    % venous filling 
Pao(1)  = 100;   % aorta
%Flow [ml/sec]
Qlv(1)  = 0;    % left ventricle (outflow)
Qp(1)   = 0;    % peripheral resistance
Qv(1)   = 0;    % ventricle filling (inflow)
 
%% Main Program
% Initialization of continuous variables:
Vlv_full = []; Va_full = []; Vv_full = [];
Plv_full = []; Pa_full = []; Pv_full = []; Pao_full = [];
Qlv_full = []; Qp_full = []; Qv_full = [];

% Main Program:
for CycleIdx = 1:Heart_cycles % Main loop for each heart cycle  
    for StepInCycle = 2:N_per_cycle
        % Volumes [ml]
        Vlv(StepInCycle) = Vlv(StepInCycle-1) + (Qv(StepInCycle-1) - Qlv(StepInCycle-1)) * dt;
        Va(StepInCycle) = Va(StepInCycle-1) + (Qlv(StepInCycle-1) - Qp(StepInCycle-1)) * dt;
        Vv(StepInCycle) = Vv(StepInCycle-1) + (Qp(StepInCycle-1) - Qv(StepInCycle-1)) * dt;

        % Pressures [mmHg]
        Plv(StepInCycle) = E(StepInCycle) * (Vlv(StepInCycle) - V0);
        Pa(StepInCycle) = Va(StepInCycle) / Ca;
        Pv(StepInCycle) = Vv(StepInCycle) / Cv;

        % Update aortic pressure
        if Plv(StepInCycle) > Pa(StepInCycle)
            Pao(StepInCycle) = Plv(StepInCycle);
        else
            Pao(StepInCycle) = Pa(StepInCycle) + Qlv(StepInCycle-1) * Ra;
        end

        % Flows [ml/sec]
        Qlv(StepInCycle) = max(0, (Pao(StepInCycle) - Pa(StepInCycle)) / Ra);  
        Qp(StepInCycle) = max(0, (Pa(StepInCycle) - Pv(StepInCycle)) / Rp); 
        Qv(StepInCycle) = max(0, (Pv(StepInCycle) - Plv(StepInCycle)) / Rv);  
    end

    % Save each variable from the current cycle to a continuous variable:
    Vlv_full = [Vlv_full, Vlv];
    Va_full = [Va_full, Va];
    Vv_full = [Vv_full, Vv];
    Plv_full = [Plv_full, Plv];
    Pa_full = [Pa_full, Pa];
    Pv_full = [Pv_full, Pv];
    Pao_full = [Pao_full, Pao];
    Qlv_full = [Qlv_full, Qlv];
    Qp_full = [Qp_full, Qp];
    Qv_full = [Qv_full, Qv];

    % Update the initial variables before the next cycle:
    % Volumes [ml]
    Vlv(1) = Vlv(end); % Left ventricle
    Va(1) = Va(end); % Arteries
    Vv(1) = Vv(end); % Veins

    % Pressures [mmHg]
    Plv(1) = Plv(end); % Left ventricle
    Pa(1) = Pa(end); % Arterial capacitor
    Pv(1) = Pv(end); % Venous filling
    Pao(1) = Pao(end); % Aorta

    % Flows [ml/sec]
    Qlv(1) = Qlv(end); % Left ventricle (outflow)
    Qp(1) = Qp(end); % Peripheral resistance
    Qv(1) = Qv(end); % Ventricle filling (inflow)
end

 
%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time vector for plotting
total_steps = length(Pao_full);
t = (0:(total_steps - 1)) * dt;

% Plot aortic pressure (Pao) over time
figure;
plot(t, Pao_full, 'b-', 'LineWidth', 1.5); % Plot Pao as a function of time
hold on;

% Initialize arrays to store maxima values and their corresponding times
Pao_max_values = zeros(1, Heart_cycles);
time_max_values = zeros(1, Heart_cycles);

% Find and mark the maxima of Pao for each heart cycle
for cycle = 1:Heart_cycles
    % Indices for the current cycle
    start_idx = (cycle - 1) * N_per_cycle + 1;
    end_idx = cycle * N_per_cycle;

    % Find the maximum Pao value within the current cycle
    [Pao_max_values(cycle), max_idx] = max(Pao_full(start_idx:end_idx));
    time_max_values(cycle) = t(start_idx + max_idx - 1);

    % Mark the maximum point on the plot
    scatter(time_max_values(cycle), Pao_max_values(cycle), 'r', 'filled'); % Red dot for maxima
end

% Add title and labels to the plot
title('Aortic Pressure (Pao) Over Time with Max', 'FontSize', 16);
xlabel('Time [sec]', 'FontSize', 14);
ylabel('Pressure [mmHg]', 'FontSize', 14);

% Add legend
legend('Aortic Pressure', 'Max Points', 'FontSize', 12);

% Add grid for better visualization
grid on;

% Hold off to finish the plot
hold off;

% Selected cycle for analysis (cycle number 20)
selected_cycle = 20;

% Find start and end indices for cycle 20
start_idx_20 = (selected_cycle - 1) * N_per_cycle + 1;
end_idx_20 = selected_cycle * N_per_cycle;

% 1.2.1 Find the maximum value in cycle 20
[max_Pao_20, max_idx_20] = max(Pao_full(start_idx_20:end_idx_20));
fprintf('1.2.1 The maximum value in cycle 20 is: %.2f mmHg\n', max_Pao_20);

%1.2.2
% Establish the acceptable range for stabilization as ±1% around the maximum value in cycle 20
Pao_lower_bound = 0.99 * Pao_max_values(20);
Pao_upper_bound = 1.01 * Pao_max_values(20); 

% Find the first cycle where the peak pressures fall within the stabilization range
stabilization_cycle_idx = find(Pao_max_values >= Pao_lower_bound & Pao_max_values <= Pao_upper_bound, 1, 'first');

% Assign the stabilization cycle
stabilization_cycle = stabilization_cycle_idx;

% Print the result
fprintf('1.2.2 The pressure stabilizes after %d heartbeats.\n', stabilization_cycle);


% 1.2.3 Calculate the time difference between successive maxima, starting from cycle 6
% Extract time maxima from cycle 6 onwards
time_max_values_stabilized = time_max_values(stabilization_cycle:end);

% Calculate differences between successive maxima for stabilized region
time_difference = diff(time_max_values_stabilized);

% Print the result
fprintf('1.2.3 Average time between successive max (from cycle 6 onwards): %.4f seconds\n', mean(time_difference));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define time vector for plotting
t = (0:(length(Pao_full) - 1)) * dt;

% Define the range for the selected single cardiac cycle
selected_cycle = 10; %Selected cycle after stabilization
start_idx = (selected_cycle - 1) * N_per_cycle + 1;
end_idx = selected_cycle * N_per_cycle;
t_cycle = t(start_idx:end_idx); % Time vector for the selected cycle

% Identify valve events
% Aortic valve events
aortic_open_idx = find(Plv_full(start_idx:end_idx) >= Pao_full(start_idx:end_idx), 1, 'first');
% Aortic valve closing: happens when Pao > Plv after opening
aortic_close_idx = find(Pao_full(start_idx:end_idx) > Plv_full(start_idx:end_idx) & ...
                        (1:length(Pao_full(start_idx:end_idx))') > aortic_open_idx, 1, 'first');

% AV valve events
av_close_idx = find(Plv_full(start_idx:end_idx) > Pv_full(start_idx:end_idx), 1, 'first');
av_open_idx = find(Pv_full(start_idx:end_idx) > Plv_full(start_idx:end_idx) & ...
                   (1:length(Pv_full(start_idx:end_idx))') > aortic_close_idx, 1, 'first');

% Plot all parameters together in one figure
figure;

% Plot Plv, Pao, Pv, and Pa on the same graph
subplot(2, 1, 1);
hold on;
plot(t_cycle, Plv_full(start_idx:end_idx), 'r', 'LineWidth', 1.5, 'DisplayName', 'Left Ventricular Pressure (Plv)');
plot(t_cycle, Pao_full(start_idx:end_idx), 'b', 'LineWidth', 1.5, 'DisplayName', 'Aortic Pressure (Pao)');
plot(t_cycle, Pv_full(start_idx:end_idx), 'g', 'LineWidth', 1.5, 'DisplayName', 'Venous Pressure (Pv)');
plot(t_cycle, Pa_full(start_idx:end_idx), 'm', 'LineWidth', 1.5, 'DisplayName', 'Arterial Pressure (Pa)');

% Mark valve events
scatter(t_cycle(av_close_idx), Plv_full(start_idx + av_close_idx - 1), 100, 'filled', 'MarkerFaceColor', 'y', 'DisplayName', 'AV Valve Closes');
scatter(t_cycle(aortic_open_idx), Plv_full(start_idx + aortic_open_idx - 1), 100, 'filled', 'MarkerFaceColor', 'c', 'DisplayName', 'Aortic Valve Opens');
scatter(t_cycle(aortic_close_idx), Pa_full(start_idx + aortic_close_idx - 1), 100, 'filled', 'MarkerFaceColor', 'k', 'DisplayName', 'Aortic Valve Closes');
scatter(t_cycle(av_open_idx), Pv_full(start_idx + av_open_idx - 1), 100, 'filled', 'MarkerFaceColor', 'm', 'DisplayName', 'AV Valve Opens');

% Add systole and diastole text
[~, systole_peak_idx] = max(Pao_full(start_idx:end_idx)); % Systole peak index
text(t_cycle(systole_peak_idx), max(Pao_full(start_idx:end_idx)) * 1.05, 'Systole', 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
text(t_cycle(round(length(t_cycle) * 0.7)), min(Pa_full(start_idx:end_idx)) * 1.1, ...
    'Diastole', 'FontSize', 12, 'Color', 'b', 'HorizontalAlignment', 'center');
% Add title, legend, and labels
title('Pressure Changes During a Single Cardiac Cycle', 'FontSize', 14);
xlabel('Time [sec]', 'FontSize', 12);
ylabel('Pressure [mmHg]', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid on;
hold off;

% Plot Left Ventricular Volume (Vlv) in the second subplot
subplot(2, 1, 2);
plot(t_cycle, Vlv_full(start_idx:end_idx), 'k', 'LineWidth', 1.5);
hold on;

% Highlight EDV and ESV
[~, edv_idx] = max(Vlv_full(start_idx:end_idx)); % End-diastolic volume (EDV)
EDV = Vlv_full(start_idx + edv_idx - 1);
scatter(t_cycle(edv_idx), EDV, 100, 'r', 'filled');

[~, esv_idx] = min(Vlv_full(start_idx:end_idx)); % End-systolic volume (ESV)
ESV = Vlv_full(start_idx + esv_idx - 1);
scatter(t_cycle(esv_idx), ESV, 100, 'b', 'filled');

% Calculate and display Stroke Volume (SV)
SV = EDV - ESV;
legend('Left Ventricular Volume', sprintf('EDV: %.2f ml', EDV), sprintf('ESV: %.2f ml', ESV), 'FontSize', 10, 'Location', 'best');

% Add title, labels, and legend
title('Left Ventricular Volume (Vlv)', 'FontSize', 14);
xlabel('Time [sec]', 'FontSize', 12);
ylabel('Volume [ml]', 'FontSize', 12);
grid on;
hold off;

% Print values
fprintf('2 ');
fprintf('End-Diastolic Volume (EDV): %.2f ml\n', EDV);
fprintf('End-Systolic Volume (ESV): %.2f ml\n', ESV);
fprintf('Stroke Volume (SV): %.2f ml\n', SV);

% Print valve event times
if ~isempty(aortic_open_idx)
    fprintf('Aortic Valve Opens at Time: %.4f sec\n', t_cycle(aortic_open_idx));
end
if ~isempty(aortic_close_idx)
    fprintf('Aortic Valve Closes at Time: %.4f sec\n', t_cycle(aortic_close_idx));
end
if ~isempty(av_open_idx)
    fprintf('AV Valve Opens at Time: %.4f sec\n', t_cycle(av_open_idx));
end
if ~isempty(av_close_idx)
    fprintf('AV Valve Closes at Time: %.4f sec\n', t_cycle(av_close_idx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3.1 

% Define key indices based on cardiac events
mvc_idx = av_close_idx; % MVC: Mitral Valve Closes
avo_idx = aortic_open_idx; % AVO: Aortic Valve Opens
avc_idx = aortic_close_idx; % AVC: Aortic Valve Closes
mvo_idx = av_open_idx; % MVO: Mitral Valve Opens

% Calculate time for each point relative to t2
t2 = 0; % t2 is the starting point (MVC)
t3 = t_cycle(avo_idx) - t_cycle(mvc_idx); % Time for t3 relative to t2
t4 = t_cycle(avc_idx) - t_cycle(mvc_idx); % Time for t4 relative to t2
t1 = t_cycle(mvo_idx) - t_cycle(mvc_idx); % Time for t1 relative to t2

% Assign volumes and pressures
V2 = Vlv_full(start_idx + mvc_idx - 1); % Volume at t2 (MVC)
P2 = Plv_full(start_idx + mvc_idx - 1); % Pressure at t2 (MVC)

V3 = Vlv_full(start_idx + avo_idx - 1); % Volume at t3 (AVO)
P3 = Plv_full(start_idx + avo_idx - 1); % Pressure at t3 (AVO)

V4 = Vlv_full(start_idx + avc_idx - 1); % Volume at t4 (AVC)
P4 = Plv_full(start_idx + avc_idx - 1); % Pressure at t4 (AVC)

V1 = Vlv_full(start_idx + mvo_idx - 1); % Volume at t1 (MVO)
P1 = Plv_full(start_idx + mvo_idx - 1); % Pressure at t1 (MVO)

% Display results
fprintf('3.1 ');
fprintf('Point t2 (MVC): t = %.2f s, V = %.2f ml, P = %.2f mmHg\n', t2, V2, P2);
fprintf('Point t3 (AVO): t = %.2f s, V = %.2f ml, P = %.2f mmHg\n', t3, V3, P3);
fprintf('Point t4 (AVC): t = %.2f s, V = %.2f ml, P = %.2f mmHg\n', t4, V4, P4);
fprintf('Point t1 (MVO): t = %.2f s, V = %.2f ml, P = %.2f mmHg\n', t1, V1, P1);

%% 3.2.1 - Linear Elastance Calculation

%{
% Define the time points (in s)
time_points = [0.00, 0.03, 0.16, 0.24]; % Time in s
%}

time_points = [t2, t3, t4, t1]; % Time in s

% Calculate elastance values 
E1 = P1 / (V1 - V0);
E2 = P2 / (V2 - V0);
E3 = P3 / (V3 - V0);
E4 = P4 / (V4 - V0);

% Create a list of elastance values
E_values = [E2, E3, E4, E1]; % Start from t2 to match the sequence t2 -> t3 -> t4 -> t1

% Generate a dense time vector for interpolation
time_dense = linspace(time_points(1), time_points(end), N_per_cycle); 

% Interpolate elastance values linearly across the given time points
E_dense = interp1(time_points, E_values, time_dense, 'linear');

% Calculate t2_new
t1 = time_points(4); % Time of t1
t2 = time_points(1); % Time of t2
cycle_duration = N_per_cycle * dt; % Full heart cycle duration
t2_new =  t2 + cycle_duration; % New t2 (end of diastole)

% Define the linear connection between t1 and t2_new
connection_time = linspace(t1, t2_new, 100);
connection_elastance = linspace(E_values(end), E_values(1), 100);

figure;
hold on;

% Plot the linear elastance curve
plot(time_dense, E_dense, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Elastance');

% Mark specific elastance points (t2, t3, t4, t1) in blue
scatter(time_points, E_values, 100, 'b', 'filled', 'DisplayName', 'Key Points (t2, t3, t4, t1)')

% Mark t2_new in blue
scatter(t2_new, connection_elastance(end), 100, 'm', 'filled', 'DisplayName', 't2_{new}');

% Plot the linear connection to t2_new in red
plot(connection_time, connection_elastance, 'r-', 'LineWidth', 2, 'DisplayName', 'Connection to t2_{new}');

% Add labels, legend, and grid
title('Changing Elasticity Over Time', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12); % Time in seconds
ylabel('Elastance (mmHg/ml)', 'FontSize', 12);
legend('Location', 'Best');
grid on;
hold off;

%% 3.2.2 - Midpoint Elastance Calculation 

% Define t_d as the exact midpoint for ventricular ejection (between t3 and t4)
t3 = time_points(2); % Time of t3 (AVO) in s
t4 = time_points(3); % Time of t4 (AVC) in s
t_d = t3 + 0.5 * (t4 - t3); % Exact midpoint in s

% Calculate the elastance value at t_d (using the interpolated linear curve)
E_td = interp1(time_dense, E_dense, t_d, 'linear'); % Use interpolated dense curve

% Calculate t2_new
t1 = time_points(4); % Time of t1
t2 = time_points(1); % Time of t2
cycle_duration = N_per_cycle * dt; % Full heart cycle duration
t2_new =  t2 + cycle_duration; % New t2 (end of diastole)

% Define the linear connection between t1 and t2_new
connection_time = linspace(t1, t2_new, 100);
connection_elastance = linspace(E_values(end), E_values(1), 100);

figure;
hold on;

% Plot the linear elastance curve
plot(time_dense, E_dense, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Elastance');

% Mark specific elastance points (t2, t3, t4, t1) in blue
scatter(time_points, E_values, 100, 'b', 'filled', 'DisplayName', 'Key Points (t2, t3, t4, t1)');

% Mark the midpoint t_d in green
scatter(t_d, E_td, 150, 'g', 'filled', 'DisplayName', 't_d');

% Mark t2_new in purple
scatter(t2_new, connection_elastance(end), 150, 'm', 'filled', 'DisplayName', 't2_{new}');

% Plot the linear connection to t2_new in red
plot(connection_time, connection_elastance, 'r-', 'LineWidth', 2, 'DisplayName', 'Connection to t2_{new}');

% Add labels, legend, and grid
title('Changing Elasticity Over Time (Including td)', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12); % Time in milliseconds
ylabel('Elastance (mmHg/ml)', 'FontSize', 12);
legend('Location', 'Best');
grid on;
hold off;

% Print Results for t_d
fprintf('Midpoint t_d and Elastance:\n');
fprintf('t_d = %.4f s\n', t_d);
fprintf('E(t_d) = %.4f mmHg/ml\n', E_td);
fprintf('t2_new = %.4f s\n', t2_new);

%% 3.2.3 

% Define the time vector for E (sinusoidal timeline in s)
t_model = (0:(length(E) - 1)) * dt; % Time in s

time_offset = 0.02; 
time_dense_shifted = time_dense + time_offset;

t1_original = time_points(4); % Original t1 in ms
t2 = time_points(1); % Time of t2 in ms
t1_shifted = t1_original + time_offset; % Shifted t1

% Define t2_new (end of the diastole)
% Verify linear model duration matches cycle duration
cycle_duration = N_per_cycle * dt; % Full heart cycle duration
t2_new =  time_points(1) + cycle_duration; %  t2_new based on full cycle duration


% Define the time and elastance vectors for the linear connection between t1_shifted and t2_new
connection_time = linspace(t1_shifted, t2_new, 100); % Linear time vector
connection_elastance = linspace(E_dense(end), E_dense(1), 100); % Linear elastance decrease

% Plot the graph
figure;
hold on;

% Plot the original shifted elastance curve (red, for all the original points)
plot(time_dense_shifted, E_dense, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Elastance');

% Plot the linear connection between t1_shifted and t2_new (red)
plot(connection_time, connection_elastance, 'r-', 'LineWidth', 2, 'DisplayName', 'Connection to t2_{new}');

% Plot the sinusoidal elastance curve (black)
plot(t_model, E, 'k-', 'LineWidth', 2, 'DisplayName', 'Sinusoidal Elastance');

% Add the key points in blue
scatter(time_points + time_offset, E_values, 100, 'b', 'filled', 'DisplayName', 'Key Points (t2, t3, t4, t1)');

% Add the t_d point (adjusted for the linear curve shift) in green
scatter(t_d + time_offset, E_td, 150, 'g', 'filled', 'DisplayName', 't_d');

% Add t2_new point in magenta
scatter(t2_new, connection_elastance(end), 150, 'm', 'filled', 'DisplayName', 't2_{new}');

% Add title, labels, and legend
title('Changing Elasticity Over Time (With Sinusodial Elasticity)', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12); % Time in seconds
ylabel('Elastance (mmHg/ml)', 'FontSize', 12);
legend('Location', 'Best');
grid on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for 4.1 and 4.2
BV_percentages = 50:25:300; % Blood volume percentages
Pao_means = zeros(size(BV_percentages)); % Store mean Pao for each BV percentage
stable_cycle = stabilization_cycle; % Start of stable cycles
stable_start_idx = (stable_cycle - 1) * N_per_cycle + 1; % Index for stable phase

%% Main Program
for i = 1:length(BV_percentages)
    % Adjust initial volumes for each BV percentage
    BV_ratio = BV_percentages(i) / 100;
    Vlv(1) = 120 * BV_ratio;
    Va(1) = 270 * BV_ratio;
    Vv(1) = 2700 * BV_ratio;
    % Initialization of continuous variables:
    Vlv_full = []; Va_full = []; Vv_full = [];
    Plv_full = []; Pa_full = []; Pv_full = []; Pao_full = [];
    Qlv_full = []; Qp_full = []; Qv_full = [];
    for CycleIdx = 1:Heart_cycles % Main loop for each heart cycle  
        for StepInCycle = 2:N_per_cycle
            % Volumes [ml]
            Vlv(StepInCycle) = Vlv(StepInCycle-1) + (Qv(StepInCycle-1) - Qlv(StepInCycle-1)) * dt;
            Va(StepInCycle) = Va(StepInCycle-1) + (Qlv(StepInCycle-1) - Qp(StepInCycle-1)) * dt;
            Vv(StepInCycle) = Vv(StepInCycle-1) + (Qp(StepInCycle-1) - Qv(StepInCycle-1)) * dt;
            % Pressures [mmHg]
            Plv(StepInCycle) = E(StepInCycle) * (Vlv(StepInCycle) - V0);
            Pa(StepInCycle) = Va(StepInCycle) / Ca;
            Pv(StepInCycle) = Vv(StepInCycle) / Cv;
            % Update aortic pressure
            if Plv(StepInCycle) > Pa(StepInCycle)
                Pao(StepInCycle) = Plv(StepInCycle);
            else
                Pao(StepInCycle) = Pa(StepInCycle) + Qlv(StepInCycle-1) * Ra;
            end
            % Flows [ml/sec]
            Qlv(StepInCycle) = max(0, (Pao(StepInCycle) - Pa(StepInCycle)) / Ra);  
            Qp(StepInCycle) = max(0, (Pa(StepInCycle) - Pv(StepInCycle)) / Rp); 
            Qv(StepInCycle) = max(0, (Pv(StepInCycle) - Plv(StepInCycle)) / Rv);  
        end
        % Save each variable from the current cycle to a continuous variable:
        Vlv_full = [Vlv_full, Vlv];
        Va_full = [Va_full, Va];
        Vv_full = [Vv_full, Vv];
        Plv_full = [Plv_full, Plv];
        Pa_full = [Pa_full, Pa];
        Pv_full = [Pv_full, Pv];
        Pao_full = [Pao_full, Pao];
        Qlv_full = [Qlv_full, Qlv];
        Qp_full = [Qp_full, Qp];
        Qv_full = [Qv_full, Qv];
        % Update the initial variables before the next cycle:
        % Volumes [ml]
        Vlv(1) = Vlv(end); % Left ventricle
        Va(1) = Va(end); % Arteries
        Vv(1) = Vv(end); % Veins
        % Pressures [mmHg]
        Plv(1) = Plv(end); % Left ventricle
        Pa(1) = Pa(end); % Arterial capacitor
        Pv(1) = Pv(end); % Venous filling
        Pao(1) = Pao(end); % Aorta
        % Flows [ml/sec]
        Qlv(1) = Qlv(end); % Left ventricle (outflow)
        Qp(1) = Qp(end); % Peripheral resistance
        Qv(1) = Qv(end); % Ventricle filling (inflow)
    end
    % Calculate mean Pao for the stable cycles
    Pao_means(i) = mean(Pao_full(stable_start_idx:end));
end

%% 4.1: Print results
fprintf('4.1: Mean Pao values for different BV percentages:\n');
for i = 1:length(BV_percentages)
    fprintf('BV = %d%%: Mean Pao = %.2f mmHg\n', BV_percentages(i), Pao_means(i));
end

%% 4.2: Plot results
figure;
plot(BV_percentages, Pao_means, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
idx_100 = find(BV_percentages == 100); % Find index for 100% BV
scatter(BV_percentages(idx_100), Pao_means(idx_100), 100, 'r', 'filled', 'DisplayName', '100% BV');
title('Mean Pao vs. BV Percentage', 'FontSize', 16);
xlabel('BV [% of original]', 'FontSize', 14); % Updated label
ylabel('Mean Pao [mmHg]', 'FontSize', 14);
legend('Mean Pao', '100% BV', 'Location', 'Best');
grid on;
hold off;