# Cardiovascular System Simulation + PID Control (MATLAB)

MATLAB project that simulates a simplified cardiovascular system and analyzes aortic pressure dynamics.
Includes:
1) Baseline hemodynamic simulation and event-based analysis (valve timing, EDV/ESV/SV, elastance modeling).
2) Closed-loop control using a PID controller to stabilize mean aortic pressure under a blood-volume disturbance.

## Features
- Time-domain simulation of LV / arterial / venous compartments (pressures, volumes, flows)
- Aortic pressure (Pao) tracking across multiple cardiac cycles
- Automatic detection of valve events (opening/closing)
- Stroke volume calculation (EDV, ESV, SV)
- Elastance modeling (linear reconstruction + comparison to sinusoidal elastance)
- Disturbance simulation (sudden blood-volume drop)
- PID control to reduce steady-state error and improve recovery dynamics

## Repository structure
- `proj_cardio_part1.m` - Baseline cardiovascular simulation + analysis (Q1-Q4 style outputs)
- `final_proj_part2.m` - Disturbance + PID control experiments (Q5-Q6 style outputs)

## Requirements
- MATLAB R2020b+ (should also work on nearby versions)
- No external toolboxes required (unless your MATLAB setup needs one for specific plotting)

## How to run

### Part 1 - Baseline simulation and analysis
1. Open MATLAB and set the working directory to the repository folder.
2. Run:
   - `proj_cardio_part1`
3. Outputs:
   - Plots for Pao over time (with cycle maxima)
   - Single-cycle pressure plots with valve events
   - LV volume plot with EDV/ESV and stroke volume
   - Elastance plots (linear model + comparison)

### Part 2 - Disturbance + PID control
1. Run:
   - `final_proj_part2`
2. Outputs:
   - Mean Pao vs. heart cycle number (with PID vs without PID)
   - Metrics: steady-state error, undershoot, settling time, rise time
   - Sensitivity plots showing the effect of Kp / Ki / Kd on mean Pao

## Notes
- The model is a simplified lumped-parameter representation (not a full physiological model).
- PID control is demonstrated as a feedback approach to compensate for a simulated bleeding disturbance.

## Technologies
- MATLAB
- Numerical simulation (ODE-like discrete time stepping)
- Signal/metrics extraction (peaks, mean values, thresholds)
- Control systems concepts (PID: proportional, integral, derivative)

## License
Add a license if needed (MIT is common for course projects).
