# Cardiovascular System Simulation + PID Control (MATLAB)

MATLAB project that simulates a lumped-parameter cardiovascular system (left ventricle, arterial and venous compartments) and analyzes the aortic pressure response.  
Part 2 extends the simulation with a disturbance (blood volume drop) and a PID controller to restore the mean aortic pressure.

## Technologies
- MATLAB (scripts + local function)
- Control concepts: PID (Proportional-Integral-Derivative)
- Signal/physiology modeling: lumped-parameter cardiovascular model

## Requirements
- MATLAB R2016b+ recommended (supports local functions inside scripts reliably)
- No external toolboxes required (basic MATLAB is enough)

## Part 1 - Cardiovascular System Simulation
**Goal:** Simulate pressures, volumes, and flows over multiple heart cycles and analyze Pao (aortic pressure).

### What it does
- Builds elastance curve for the left ventricle (systole/diastole).
- Simulates:
  - Volumes: `Vlv`, `Va`, `Vv`
  - Pressures: `Plv`, `Pa`, `Pv`, `Pao`
  - Flows: `Qlv`, `Qp`, `Qv`
- Produces plots:
  - Pao over time with peak markers
  - Pressures in one cardiac cycle + valve events
  - Vlv over one cardiac cycle (EDV/ESV/SV)
  - Mean Pao vs blood volume percentage (BV%)

### Run
1. Open MATLAB.
2. Go to `part1/`.
3. Run:
   - `proj_cardio_part1.m`

Outputs: figures + printed metrics (max Pao, stabilization cycle, time between peaks, EDV/ESV/SV, valve event times, BV sweep).

  ## Part 2 - Disturbance + PID Control
  **Goal:** Apply a disturbance (sudden blood volume reduction) and compare mean Pao with and without PID.
  
  ### What it does
  - Simulates 200 heart cycles.
  - Introduces a disturbance at cycle 100: `disturbance` (e.g., 0.15 = 15% BV decrease).
  - Computes mean Pao per cycle:
    - Without control: `Kp=0, Ki=0, Kd=0`
    - With PID: `Kp, Ki, Kd` tuned values
  - Produces:
    - Plot of mean Pao across cycles (with vs without PID)
    - Performance metrics:
      - steady-state error (before/after disturbance)
      - undershoot (as defined in the analysis)
      - settling time (heartbeats + seconds)
      - rise time (10% to 90% after disturbance)
    - Sensitivity plots for varying `Kp`, `Ki`, `Kd`
  
  ### Run
  1. Open MATLAB.
  2. Go to `part2/`.
  3. Run:
     - `final_proj_part2.m`
  
  Note: `final_proj_part2.m` includes a local function `Pao_part2_proj(...)` at the bottom. Keep it in the same file.
  
  ## Parameters You Can Edit
  - Heart rate: `HR`
  - Time step: `dt`
  - Number of cycles: `Heart_cycles`
  - Disturbance magnitude: `disturbance`
  - PID gains: `Kp`, `Ki`, `Kd`
  - Stabilization windows used for metrics (cycle ranges)
  
  ## Results 
  - Part 1 shows convergence to a stable cardiac cycle and supports valve timing + PV-related analysis.
  - Part 2 demonstrates how PID control improves recovery of mean aortic pressure after a sudden BV decrease, compared to the uncontrolled system.
  
  ## Notes
  - If you get errors about local functions, ensure you are using a MATLAB version that supports functions at the end of scripts (R2016b+). Otherwise, move `Pao_part2_proj` into its own file `Pao_part2_proj.m`.
