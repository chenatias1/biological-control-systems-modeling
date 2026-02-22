clc ; clear all ; close all ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Simulation of the Cardiovascular System %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% New parameters
Heart_cycles = 200;
disturbance = 0.15; % disturbance simulating bleeding, meaning a sudden 15% decrease in BV
cycle_num = 1:Heart_cycles;

%PID controlers 
Kp = 0.00005;
Ki = 0.005;
Kd = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot
Pao_wo_control_system = Pao_part2_proj(Heart_cycles, disturbance, 0, 0, 0);
Pao_mean = Pao_part2_proj(Heart_cycles, disturbance, Kp, Ki, Kd);

figure
plot(cycle_num, Pao_mean, 'b', 'LineWidth', 1.5); 
hold on
plot(cycle_num, Pao_wo_control_system, 'r', 'LineWidth', 1.5); 
ax = gca;
ax.FontSize = 12; 
title('Average Pao As Function Of Cycle Number', 'FontSize', 14);
xlabel('Heart Cycle Number', 'FontSize', 12); 
ylabel('Average Pao [mmHg]', 'FontSize', 12); 
legend('with PID', 'without PID', 'FontSize', 12); 
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pao_in = 77.828 ; 
disturbance_beat = 100;
% Define stable cycles before and after disturbance
stable_cycles_before_start = 50;  % First stable cycle before disturbance
stable_cycles_before_end = 90;    % Last stable cycle before disturbance
stable_cycles_after_start = 180;  % First stable cycle after disturbance
stable_cycles_after_end = 200;    % Last stable cycle after disturbance

%% Error calculation

% System with PID control
% Calculate the error as the absolute difference between the input value and the mean output
error_with_control_before = abs(Pao_in - mean(Pao_mean(stable_cycles_before_start:stable_cycles_before_end)));
error_with_control_after = abs(Pao_in - mean(Pao_mean(stable_cycles_after_start:stable_cycles_after_end)));
% System without PID control
% Calculate the error as the absolute difference between the input value and the mean output
error_without_control_before = abs(Pao_in - mean(Pao_wo_control_system(stable_cycles_before_start:stable_cycles_before_end)));
error_without_control_after = abs(Pao_in - mean(Pao_wo_control_system(stable_cycles_after_start:stable_cycles_after_end)));
% Display results in scientific notation if values are very small
fprintf('Error with PID control before disturbance: %.5e mmHg\n', error_with_control_before);
fprintf('Error with PID control after disturbance: %.5e mmHg\n', error_with_control_after);
fprintf('Error without PID control before disturbance: %.5e mmHg\n', error_without_control_before);
fprintf('Error without PID control after disturbance: %.5e mmHg\n', error_without_control_after);

%% O.S (undershoot in our case)

% Before disturbance (all cycles before cycle 100)
min_value_before = min(Pao_mean(1:disturbance_beat-1)); % Minimum value before disturbance
us_before = ((Pao_in - min_value_before)/Pao_in) * 100; % Undershoot percentage before disturbance
% After disturbance (all cycles after cycle 100)
min_value_after = min(Pao_mean(disturbance_beat:end)); % Minimum value after disturbance
us_after = ((Pao_in - min_value_after)/Pao_in) * 100; % Undershoot percentage after disturbance
% Display results
fprintf('Minimum value before disturbance: %.5f mmHg\n', min_value_before);
fprintf('Undershoot before disturbance: %.5e%%\n', us_before);
fprintf('Minimum value after disturbance: %.5f mmHg\n', min_value_after);
fprintf('Undershoot after disturbance: %.5e%%\n', us_after);

%% Settling Time - in Heartbeats and Seconds

start_cycle_before = 5; % Start considering values from cycle 5
% Time parameters
HR = 80; % Heart rate in BPM
dt = 5e-4; % Time step in seconds
N_per_cycle = round(60 / (HR * dt)); % Number of steps per heart cycle
time_per_cycle = N_per_cycle * dt; % Duration of one cycle in seconds
% Define stable range
stable_min = Pao_in * 0.99; % Lower bound of 1% range
stable_max = Pao_in * 1.01; % Upper bound of 1% range
% Initialize variables
settling_heartbeats_before = NaN; % Default if not found
settling_heartbeats_after = NaN; % Default if not found
% Before disturbance (cycle 5 to 99)
for cycle_idx = start_cycle_before:disturbance_beat-1
    % Check if the next 5 cycles are within the stable range
    if all(Pao_mean(cycle_idx:cycle_idx+4) >= stable_min & Pao_mean(cycle_idx:cycle_idx+4) <= stable_max)
        settling_heartbeats_before = cycle_idx - start_cycle_before; % Heartbeats from the start
        break;
    end
end
% After disturbance (cycle 101 to end)
for cycle_idx = disturbance_beat+1:Heart_cycles-5
    % Check if the next 5 cycles are within the stable range
    if all(Pao_mean(cycle_idx:cycle_idx+4) >= stable_min & Pao_mean(cycle_idx:cycle_idx+4) <= stable_max)
        settling_heartbeats_after = cycle_idx - disturbance_beat; % Heartbeats from disturbance
        break;
    end
end
% Convert heartbeats to seconds
settling_time_seconds_before = settling_heartbeats_before * time_per_cycle;
settling_time_seconds_after = settling_heartbeats_after * time_per_cycle;
% Display results
if isnan(settling_heartbeats_before)
    fprintf('Settling heartbeats before disturbance: Not found (system did not stabilize)\n');
else
    fprintf('Settling heartbeats before disturbance: %d heartbeats\n', settling_heartbeats_before);
    fprintf('Settling time before disturbance: %.5f seconds\n', settling_time_seconds_before);
end

if isnan(settling_heartbeats_after)
    fprintf('Settling heartbeats after disturbance: Not found (system did not stabilize)\n');
else
    fprintf('Settling heartbeats after disturbance: %d heartbeats\n', settling_heartbeats_after);
    fprintf('Settling time after disturbance: %.5f seconds\n', settling_time_seconds_after);
end

%% Rise Time 

% Define bounds for 10% and 90%
lower_bound = 0.1 * Pao_in; % 10% of Pao_in
upper_bound = 0.9 * Pao_in; % 90% of Pao_in
% Initialize variables
cycle_10_percent = NaN; % Cycle index where 10% is reached
cycle_90_percent = NaN; % Cycle index where 90% is reached
% After disturbance (cycles 101 to end)
for cycle_idx = disturbance_beat+1:Heart_cycles
    % Check for 10% boundary
    if isnan(cycle_10_percent) && Pao_mean(cycle_idx) >= lower_bound
        cycle_10_percent = cycle_idx; % First cycle reaching 10%
    end
    % Check for 90% boundary
    if isnan(cycle_90_percent) && Pao_mean(cycle_idx) >= upper_bound
        cycle_90_percent = cycle_idx; % First cycle reaching 90%
        break; % Exit loop once 90% is reached
    end
end

% Calculate response time in seconds
if ~isnan(cycle_10_percent) && ~isnan(cycle_90_percent)
    response_time_seconds = (cycle_90_percent - cycle_10_percent) * time_per_cycle;
else
    response_time_seconds = NaN; % Default if bounds are not reached
end

% Display results
if isnan(cycle_10_percent)
    fprintf('Cycle for 10%% of desired value: Not found\n');
else
    fprintf('Cycle for 10%% of desired value: %d\n', cycle_10_percent);
end

if isnan(cycle_90_percent)
    fprintf('Cycle for 90%% of desired value: Not found\n');
else
    fprintf('Cycle for 90%% of desired value: %d\n', cycle_90_percent);
end

if isnan(response_time_seconds)
    fprintf('rise time: Not found (bounds not reached)\n');
else
    fprintf('rise time: %.5f seconds\n', response_time_seconds);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define PID controller parameters
Kp_original = 0.00005;
Ki_original = 0.005;
Kd_original = 0.001;

% Varying Kp
Kp_values =  [0.5, 10, 50, 100, 500]* Kp_original;
colors = lines(length(Kp_values));
figure;
hold on;

Pao_Average_original = Pao_part2_proj(Heart_cycles, disturbance, Kp_original, Ki_original, Kd_original);

% Plot the original Kp
plot(cycle_num, Pao_Average_original, '--k', 'LineWidth', 3, 'DisplayName', 'Original Kp');

for i = 1:length(Kp_values)
    Kp = Kp_values(i);
    Pao_Average = Pao_part2_proj(Heart_cycles, disturbance, Kp, Ki_original, Kd_original);
    relative_value = Kp / Kp_original; % Calculate relative value
    plot(cycle_num, Pao_Average, 'Color', colors(i, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%.2f Kp_{original}', relative_value));
end
title('Effect of Kp on Mean Pao', 'FontSize', 14);
xlabel('Heart Cycle Number', 'FontSize', 12);
ylabel('Average Pao [mmHg]', 'FontSize', 12);
legend show;
grid on;

% Varying Ki
Ki_values = [0.1, 0.5, 5, 10] * Ki_original;
figure;
hold on;

% Plot the original Ki
plot(cycle_num, Pao_Average_original, '--k', 'LineWidth', 1.5, 'DisplayName', 'Original Ki');

for Ki = Ki_values
    Pao_Average = Pao_part2_proj(Heart_cycles, disturbance, Kp_original, Ki, Kd_original);
    relative_value = Ki / Ki_original; % Calculate relative value
    plot(cycle_num, Pao_Average, 'DisplayName', sprintf('%.2f Ki_{original}', relative_value));
end
title('Effect of Ki on Mean Pao', 'FontSize', 14);
xlabel('Heart Cycle Number', 'FontSize', 12);
ylabel('Average Pao [mmHg]', 'FontSize', 12);
legend show;
grid on;

% Varying Kd
Kd_values = [0.5, 10, 50,100] * Kd_original;
figure;
hold on;

% Plot the original Kd
plot(cycle_num, Pao_Average_original, '--k', 'LineWidth', 1.5, 'DisplayName', 'Original Kd');

for Kd = Kd_values
    Pao_Average = Pao_part2_proj(Heart_cycles, disturbance, Kp_original, Ki_original, Kd);
    relative_value = Kd / Kd_original; % Calculate relative value
    plot(cycle_num, Pao_Average, 'DisplayName', sprintf('%.2f Kd_{original}', relative_value));
end
title('Effect of Kd on Mean Pao', 'FontSize', 14);
xlabel('Heart Cycle Number', 'FontSize', 12);
ylabel('Average Pao [mmHg]', 'FontSize', 12);
legend show;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calulating the mean Pao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Pao_mean = Pao_part2_proj(Heart_cycles, disturbance, Kp, Ki, Kd)

    %% system parameters:
    % Time parameters
    HR           = 80 ; % [BPM] 
    dt           = 5e-4        ; % [sec]
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
     
    % Initialization of continuous variables:
    Vlv_full = []; Va_full = []; Vv_full = [];
    Plv_full = []; Pa_full = []; Pv_full = []; Pao_full = []; Pao_mean = zeros([1, Heart_cycles]);
    Qlv_full = []; Qp_full = []; Qv_full = [];

    % Additional constants
    Rp_original = 1.0; % peripheral resistance
    Cv_original = 300.0; % venous compliance
    disturbance_beat = 100;
    Pao_in = 77.828;

    % Error
    error = zeros([1, Heart_cycles+1]);
    error_sum = 0  ;

%% Main Program:
    for CycleIdx = 1:Heart_cycles % Main loop for each heart cycle
        % Generating a disturbance after a typical settling beat
        if CycleIdx == disturbance_beat
            reducing_percentage = disturbance; % Try to start with smaller or no reduction
            % Volumes [ml]
            Vlv(1)  = Vlv(1) * (1 - reducing_percentage);  % Left ventricle
            Va(1)   = Va(1) * (1 - reducing_percentage);   % Arteries
            Vv(1)   = Vv(1) * (1 - reducing_percentage);   % Veins
        end
    
        for StepInCycle = 2:N_per_cycle
            % Calculating all variables for each cycle at N points:
            % Volumes [ml]
            Vlv(StepInCycle) = Vlv(StepInCycle-1) + (Qv(StepInCycle-1) - Qlv(StepInCycle-1)) * dt;
            Va(StepInCycle)  = Va(StepInCycle-1) + (Qlv(StepInCycle-1) - Qp(StepInCycle-1)) * dt;
            Vv(StepInCycle)  = Vv(StepInCycle-1) + (Qp(StepInCycle-1) - Qv(StepInCycle-1)) * dt;
    
            % Pressures [mmHg]
            Plv(StepInCycle) = E(StepInCycle) * (Vlv(StepInCycle) - V0);
            Pa(StepInCycle)  = Va(StepInCycle) / Ca;
            Pv(StepInCycle)  = Vv(StepInCycle) / Cv;
            Pao(StepInCycle) = max(Plv(StepInCycle), Pa(StepInCycle));
    
            % Flows [ml/sec]
            Qlv(StepInCycle) = max(0, (Pao(StepInCycle) - Pa(StepInCycle)) / Ra);
            Qp(StepInCycle)  = max(0, (Pa(StepInCycle) - Pv(StepInCycle)) / Rp);
            Qv(StepInCycle)  = max(0, (Pv(StepInCycle) - Plv(StepInCycle)) / Rv);
        end
    
        % Mean Aortic Pressure
        Pao_mean(CycleIdx) = mean(Pao);
    
        % Error Calculation
        error(CycleIdx) = Pao_in - Pao_mean(CycleIdx);
        error_sum = error_sum + error(CycleIdx);
    
        % PID Controller
        P = error(CycleIdx) * Kp; % Proportional term
        I = error_sum * Ki; % Integral term
        if CycleIdx == 1
            D = 0; % No derivative term in the first cycle
        else
            D = (error(CycleIdx) - error(CycleIdx-1)) * Kd; % Derivative term
        end
        PID = P + I + D; % Total PID control output
        
        % Parameter Update
        Emax = max(0.1, Emax + PID); % Update Emax based on PID, ensuring a minimum value
        Cv = max(0.001, Cv_original + PID); % Update Cv based on PID, ensuring a minimum value
        Rp = max(0.001, Rp_original + PID); % Update Rp based on PID, ensuring a minimum value

        % Save variables from the current cycle to continuous variables:
        % Volumes [ml]
        Vlv_full = [Vlv_full, Vlv];
        Va_full = [Va_full, Va];
        Vv_full = [Vv_full, Vv];
    
        % Pressures [mmHg]
        Plv_full = [Plv_full, Plv];
        Pa_full = [Pa_full, Pa];
        Pv_full = [Pv_full, Pv];
        Pao_full = [Pao_full, Pao];
    
        % Flows [ml/sec]
        Qlv_full = [Qlv_full, Qlv];
        Qp_full = [Qp_full, Qp];
        Qv_full = [Qv_full, Qv];
    
        % Update the initial variables before the next cycle:
        % Volumes [ml]
        Vlv(1) = Vlv(end);
        Va(1)  = Va(end);
        Vv(1)  = Vv(end);
    
        % Pressures [mmHg]
        Plv(1) = Plv(end);
        Pa(1)  = Pa(end);
        Pv(1)  = Pv(end);
        Pao(1) = Pao(end);
    
        % Flows [ml/sec]
        Qlv(1) = Qlv(end);
        Qp(1)  = Qp(end);
        Qv(1)  = Qv(end);
    end
end
