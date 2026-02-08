# Optimizing Fractional-Polynomial Powered Descent Guidance Laws

This repository contains the MATLAB implementation for the **"Optimizing Fractional-Polynomial Powered Descent Guidance Laws"** paper. It implements an optimization framework for the **Fractional Polynomial-Powered Descent Guidance (FP2DG)** law, optimizing three parameters to minimize fuel and throttle usage while satisfying hard constraints (glideslope, pointing, thrust).

This codebase demonstrates that by optimizing just three key variables, two polynomial powers ($\gamma_1, \gamma_2$) and the time-to-go ($t_{go}$), the **FP2DG** law can achieve near-optimal fuel performance comparable to complex numerical methods, while retaining the computational efficiency required for onboard flight computers.

## Theoretical Background

**FP2DG** enhances classical analytical methods (like Apollo Powered Descent Guidance) with the benefits of modern numerical optimizers into a single flexible framework.

*   **Guidance Law:** The commanded thrust acceleration is defined as a function of time:
    $$ \mathbf{a}_T(t) = \mathbf{a}_{T_f}^* + \mathbf{c}_1 t^{\gamma_1} + \mathbf{c}_2 t^{\gamma_2} $$
    Where $\mathbf{c}_1$ and $\mathbf{c}_2$ are coefficient vectors computed analytically to satisfy boundary conditions (Target Position $\mathbf{r}_f$ and Velocity $\mathbf{v}_f$).

*   **Optimization:** The algorithm solves a Nonlinear Programming (NLP) problem to find the optimal set $X = [\gamma_1, \gamma_2, t_{go}]$ that minimizes a combined cost function:
    $$ J = \beta \int ||\mathbf{a}_T|| dt + (1-\beta) \int ||\mathbf{a}_T||^2 dt $$
    *   $\beta \to 1$: Maximizes fuel efficiency (Minimum Control).
    *   $\beta \to 0$: Maximizes throttle smoothness (Minimum Control Effort).

    The optimizer enforces the following constraints on $X$:

    *   $\gamma_1, \gamma_2 > 0$: Prevents singularities in the $t^{\gamma}$ basis functions at $t = 0$ (landing).
    *   $\gamma_2 > \gamma_1$: Ensures the two polynomial basis functions are distinct so that $\delta \neq 0$ and $\gamma_2$ is always greater than $\gamma_1$.
    *   $t_{go} > 0$: Required for a real trajectory.
    *   Additionally, thrust magnitude is constrained at every trajectory node: $T_{min} \leq ||\mathbf{a}_T|| \cdot m \leq T_{max}$. Optional glideslope (position vector angle from vertical) and pointing (thrust vector angle from vertical) constraints can also be enforced.

*   **Simulation:** The optimized $\gamma_1, \gamma_2$ and $t_{go}$ are used in the closed-loop tracking form of FP2DG, which recomputes the commanded acceleration at each integration step using current state feedback:

    $$\mathbf{a}_T(t) = \gamma_1\!\left(\frac{k_r}{2\gamma_1+4}-1\right)\mathbf{a}_{T_f}^* + \left(\frac{\gamma_1 k_r}{2\gamma_1+4}-\gamma_1-1\right)\mathbf{g} + \frac{\gamma_1+1}{t_{go}}\left(1-\frac{k_r}{\gamma_1+2}\right)(\mathbf{v}_f^*-\mathbf{v}) + \frac{k_r}{t_{go}^2}(\mathbf{r}_f^*-\mathbf{r}-\mathbf{v}\,t_{go})$$

    where $k_r = (\gamma_2+2)(\gamma_1+2)$. Unlike the open-loop planning form, this tracking law feeds back the current position $\mathbf{r}$ and velocity $\mathbf{v}$, making it robust to perturbations. Time-to-go decrements in real time as the vehicle flies. If thrust magnitude exceeds engine limits, the command is clamped to $T_{min}$ or $T_{max}$. Near landing ($t_{go} \cdot T_{ref} < 0.2$ s), the acceleration command is frozen at its last computed value to avoid the $t_{go} \to 0$ singularity.

## Functions

### `userTest.m`
**Purpose:** Primary entry point. Defines the mission scenario, vehicle/planet parameters, constraint flags, and optimization settings, then calls `getParams.m` (or `getParamsDIVERT.m` for divert scenarios) to run the full optimization + simulation pipeline.

**User-configured variables:**
*   `beta`: Cost weighting factor (1.0 = fuel optimal, 0.0 = smoothest throttle).
*   `glideSlopeEnabled`, `pointingEnabled`, `reOptimizationEnabled`, `divertEnabled`: Constraint and mode toggles.
*   `PDIState`: Initial state at Powered Descent Initiation — altitude (m), longitude/latitude (deg), inertial velocity (m/s), flight path angle (deg), azimuth (deg).
*   `planetaryParams`: Moon radius (m), surface gravity (m/s^2), Earth gravity (m/s^2).
*   `vehicleParams`: Initial mass (kg), dry mass (kg), Isp (s), max/min thrust (N).
*   `targetState`: Landing site coordinates (deg), terminal position/velocity/acceleration targets (ENU, m or m/s), divert configuration.
*   `optimizationParams`: Initial guess `[gamma1, gamma2, tgo]`, node count ( default of 301, must be odd), glideslope/pointing settings, re-optimization frequency/stop time, gamma epsilon bounds.

**Outputs:**
*   Console: Optimal parameters (gamma1, gamma2, kr, tgo), fuel cost from optimizer and simulation.
*   Figures: Full trajectory visualization via `plotting.m`.

### `getParams.m`
**Purpose:** Central managing function. Transforms the PDI state from spherical to MCMF Cartesian coordinates, non-dimensionalizes all quantities, runs the optimizer, runs the trajectory simulation (static or re-optimizing), and post-processes results back to dimensional units.

**Inputs:**
*   `PDIState`, `planetaryParams`, `targetState`, `vehicleParams`, `optimizationParams`: Mission definition structs (from `userTest.m`).
*   `betaParam`: Cost weighting (scalar).
*   `doPlots`: Boolean, enables trajectory plots via `plotting.m`.
*   `verboseOutput`: Boolean, enables detailed console output.
*   `dispersion`: Boolean, adjusts control flow when called from `dispersionStudy.m`.
*   `runSimulation`: Boolean, whether to run ODE45 simulation after optimization.
*   `monteCarloSeed` (optional, 11th arg): Acceleration scaling matrix for `accelMonteCarlo.m`.

**Outputs:**
*   `gammaOpt`, `gamma2Opt`: Optimal fractional polynomial exponents.
*   `krOpt`: Gain for alternate form of FP2DG, `(gamma2+2)*(gamma1+2)`.
*   `tgoOpt`: Optimal time-to-go (seconds, dimensional).
*   `aTOptim`: Commanded acceleration profile from the optimizer (3xN, non-dimensional).
*   `exitflag`: `fmincon` convergence flag for the initial optimization.
*   `optFuelCost`: Fuel consumed predicted by the optimizer (kg).
*   `simFuelCost`: Fuel consumed in the ODE45 simulation (kg).
*   `aTSim`: Commanded acceleration history from the simulation (3xM).
*   `finalPosSim`: Final position error vector in ENU (m), from the simulation endpoint.
*   `optHistory`: Parameter history during re-optimization — columns: `[t_elapsed, gamma1, gamma2, kr, tgo]` (non-dimensional). Empty if re-optimization is disabled.
*   `ICstates`: Table of initial condition states at each re-optimization segment (7xK: r, v, m). Empty if re-optimization is disabled.
*   `exitFlags`: Vector of `fmincon` exit flags for each re-optimization call.
*   `problemParams`: Struct of dimensional problem constants (positions, velocities, thrust limits, etc.), passed through for downstream use.
*   `nonDimParams`: Struct of non-dimensional quantities (r0ND, v0ND, ispND, etc.).
*   `refVals`: Reference values used for non-dimensionalization (L_ref, T_ref, V_ref, A_ref, M_ref).
*   `optTable`, `simTable`: 2x1 vectors `[landing error (m); fuel cost (kg)]` for optimizer vs simulation comparison.

### `getParamsDIVERT.m`
**Purpose:** Alternate of `getParams.m` for divert scenarios. Forces glideslope and pointing constraints off and re-optimization on, then runs the same pipeline as `getParams.m`.

**Inputs/Outputs:** Same as `getParams.m`, minus `optTable`/`simTable`.

### `optimizationLoop.m`
**Purpose:** Solves the NLP using `fmincon` (SQP algorithm). Finds optimal `[gamma1, gamma2, tgo]` that minimize the cost function subject to thrust limit constraints and optional glideslope/pointing constraints. Contains the objective function and nonlinear constraint function as local functions.

**Inputs:**
*   `paramsX0`: Initial guess `[gamma1, gamma2, tgo]` (non-dimensional tgo).
*   `betaParam`: Cost weighting (scalar).
*   `problemParams`: Dimensional problem constants (used for glideslope/pointing constraint evaluation).
*   `nonDimParams`: Non-dimensional state/vehicle values.
*   `optimizationParams`: Solver settings — node count, constraint flags, glideslope/pointing parameters, gamma epsilon bounds.
*   `refVals`: Reference scaling values.
*   `delta_t`: Beyond-termination targeting time offset (non-dimensional).
*   `verboseOutput`: Boolean for detailed console output.
*   `dispersion`: Boolean (adjusts behavior for dispersion study calls).

**Outputs:**
*   `optParams`: Optimized `[gamma1, gamma2, tgo]` (non-dimensional tgo).
*   `optCost`: Final objective function value.
*   `aTOptim`: Optimal acceleration profile (3x301, non-dimensional).
*   `mOptim`: Mass profile along the optimal trajectory (1x301, non-dimensional).
*   `rdOptim`: Position profile along the optimal trajectory (3x301, non-dimensional).
*   `vdOptim`: Velocity profile along the optimal trajectory (3x301, non-dimensional).
*   `exitflag`: `fmincon` convergence status (1 = converged, <=0 = failed).

**Constraints:** Linear inequality enforces `gamma1 >= eps`, `gamma2 >= gamma1 + eps`, `tgo >= 0.01`. Nonlinear constraints include up to 2x301 thrust limit constraints, plus optional glideslope (301) and pointing (301) constraints per node.

### `calculateCoeffs.m`
**Purpose:** Analytically solves for the guidance law coefficients. Given the current state, target state, gravity, and the polynomial exponents, computes `c1` and `c2` so that the trajectory satisfies the position and velocity boundary conditions.

**Inputs:**
*   `r`: Current position (3x1, non-dimensional).
*   `v`: Current velocity (3x1, non-dimensional).
*   `tgo`: Time-to-go (scalar, non-dimensional).
*   `gamma1`, `gamma2`: Fractional polynomial exponents.
*   `afStar`: Target terminal acceleration (3x1, non-dimensional).
*   `rfStar`: Target terminal position (3x1, non-dimensional).
*   `vfStar`: Target terminal velocity (3x1, non-dimensional).
*   `g`: Constant gravity vector at landing site (3x1, non-dimensional).

**Outputs:**
*   `c1`, `c2`: Guidance coefficient vectors (3x1, non-dimensional).
*   `c1_num`, `c2_num` (optional): Numerator vectors before division by delta, used for diagnostics.

### `closedLoopSim.m`
**Purpose:** Propagates the trajectory using ODE45 with fixed guidance parameters (no re-optimization during flight). At each integration step, the FP2DG guidance law recomputes the commanded acceleration using the current state and decrementing time-to-go. Thrust is clamped to engine limits if exceeded.

**Inputs:**
*   `gamma`, `gamma2`: Fixed exponents from the optimizer.
*   `tgo0`: Initial time-to-go (non-dimensional).
*   `problemParams`: Dimensional problem constants.
*   `nonDimParams`: Non-dimensional state/vehicle values (initial state, thrust limits, Isp, gravity, targets).
*   `refVals`: Reference scaling values.
*   `delta_t`: Beyond-termination targeting time offset (non-dimensional). Not implemented currently.
*   `monteCarloSeed` (optional): Acceleration scaling matrix applied to thrust commands, used by `accelMonteCarlo.m`.

**Outputs:**
*   `tTraj`: Time history (Mx1, non-dimensional).
*   `stateTraj`: State history (Mx7: columns are r_x, r_y, r_z, v_x, v_y, v_z, m, all non-dimensional).
*   `aTList`: Commanded acceleration at each time step (3xM, non-dimensional).
*   `flag_thrustGotLimited`: Boolean, true if thrust was clamped at any point during the simulation.

### `simReOpt.m`
**Purpose:** Propagates the trajectory with periodic re-optimization. Integrates forward in segments of `updateFreq` seconds, re-solves the NLP at the end of each segment using the current state as new initial conditions, then continues. Stops re-optimizing at `updateStop` seconds before landing and freezes the thrust command near touchdown to avoid the tgo→0 singularity. Also handles divert scenarios — if a divert altitude trigger is reached, switches the target position and re-optimizes once.

**Inputs:**
*   `gamma0`, `gamma20`: Initial exponents.
*   `tgo0`: Initial time-to-go (non-dimensional).
*   `problemParams`, `nonDimParams`, `refVals`: Problem definition structs.
*   `delta_t`: Beyond-termination targeting offset (non-dimensional). Not implemented.
*   `optimizationParams`: Re-optimization settings — `updateFreq` (s), `updateStop` (s), constraint flags.
*   `betaParam`: Cost weighting.
*   `verboseOutput`: Boolean for console output.
*   `divertPoint` (optional): New target position in MCMF (3x1, non-dimensional) for divert scenarios.

**Outputs:**
*   `tTraj`: Full time history across all segments (Mx1, non-dimensional).
*   `stateTraj`: Full state history (Mx7, non-dimensional).
*   `aTList`: Commanded acceleration across all segments (3xM, non-dimensional).
*   `flag_thrustGotLimited`: Boolean, true if thrust was clamped at any point.
*   `optHistory`: Re-optimization parameter log (Kx5: columns are `t_elapsed, gamma1, gamma2, kr, tgo`, all non-dimensional).
*   `ICstates`: Initial condition states at each segment boundary (7xK, non-dimensional).
*   `exitFlags`: Vector of `fmincon` exit flags for each re-optimization (Kx1).

### `simpsonComp13Integral.m`
**Purpose:** Composite Simpson's 1/3 rule numerical integration. Requires an even number of intervals (odd number of points, hence the 301 default node count).

**Inputs:**
*   `t`: Time vector (1xN).
*   `y`: Data vector to integrate (1xN).

**Outputs:**
*   `cost`: Computed integral value (scalar).

### `dispersionStudy.m`
**Purpose:** Monte Carlo dispersion analysis. Perturbs the PDI initial conditions (altitude, lat, lon, velocity, FPA, azimuth, mass) using pre-generated seed files from `Seeds/` and 3-sigma ranges, then runs `getParams.m` for each case via `parfor`. Saves results to `Dispersion301/` as .mat and .csv, generates statistical plots via `statsPlotting.m`, and prints summary tables (mean, std, min, max, 3-sigma bounds) for all key parameters.

**User-configured variables:**
*   Nominal state structs (same as `userTest.m`).
*   3-sigma dispersion ranges: `r_disp` (m), `lon_disp`/`lat_disp` (deg), `v_disp` (m/s), `fpa_disp` (deg), `azmth_disp` (deg), `mass_disp` (fraction).

**Outputs:**
*   `Dispersion301/<run_name>/` — .mat results, .csv summary, exit flag table, success table.
*   Console: Summary statistics table and exit flag breakdown.
*   Figures: Statistical distribution plots from `statsPlotting.m`.

### `rerunDispersionCase.m`
**Purpose:** Re-runs a single dispersion case by index with plotting and verbose output enabled. Used for investigating specific cases from the dispersion study.

**Inputs:**
*   `idx`: Case index into the seed arrays.
*   `PDINom`: Nominal PDI state struct.
*   `planetaryParams`, `targetState`, `vehicleNom`, `optimizationParams`: Standard mission structs.
*   `beta`: Cost weighting.
*   `seeds`: Struct of loaded seed arrays from `Seeds/`.

**Outputs:**
*   `outputSingle`: Struct with `gamma`, `gamma2`, `kr`, `tgo`, `optFuel`, `simFuel`, `finalError`.

### `accelMonteCarlo.m`
**Purpose:** Acceleration uncertainty Monte Carlo study. Applies random scaling factors (from `Seeds/accel_seeds.dat`) to the commanded acceleration during simulation to assess guidance robustness against thrust execution errors. Structure mirrors `dispersionStudy.m` but perturbs acceleration rather than initial conditions.

**User-configured variables:**
*   `accel_monte_carlo`: 3-sigma acceleration perturbation as a fraction (e.g., 0.10 = 10%).

**Outputs:**
*   `AccelMonteCarlo301/<run_name>/` — .mat results, .csv summary, exit flag/success tables.

### `gammaSweep.m`
**Purpose:** Parameter sweep over gamma1 and gamma2 at a fixed tgo. Evaluates the cost function over a grid of gamma values to visualize the cost landscape and identify the optimal set. Constraints are disabled.

**User-configured variables:**
*   `fixedTgo`: Fixed time-to-go for the sweep (seconds, dimensional).
*   `gammaRange`: `[min, max]` range for both gamma1 and gamma2.
*   `betaVal`: Cost weighting.
*   `gridResolution`: Number of grid points per axis.

**Outputs:**
*   Figures: Cost surface/contour plots over the gamma1-gamma2 space. Saved to `Gamma Sweep/`.

### `tgoSweep.m`
**Purpose:** Parameter sweep over tgo at fixed gamma1 and gamma2. Evaluates cost and constraint feasibility across a range of time-to-go values to visualize how tgo affects optimality and constraint activity.

**User-configured variables:**
*   `fixedGamma1`, `fixedGamma2`: Fixed exponents.
*   `tgoRange`: `[min, max]` range for tgo (seconds, dimensional).
*   `betaVal`: Cost weighting.

**Outputs:**
*   Figures: Cost vs tgo curves with constraint feasibility markers.

### `plotting.m`
**Purpose:** Generates trajectory visualization figures for a single run. Produces plots of throttle, altitude vs ground range, 3D trajectory, velocity components, acceleration components, mass history, and constraint boundaries. Supports overlay of multiple runs on the same figures with auto-colored legend entries.

**Inputs:**
*   `tTraj`, `stateTraj`: Time and state histories from simulation (non-dimensional).
*   `optParams`: Optimized `[gamma1, gamma2, tgo]`.
*   `optCost`: Objective function value.
*   `aTOptim`, `mOptim`, `rdOptim`, `vdOptim`: Optimizer-predicted profiles (3x301 or 1x301).
*   `aTSim`: Simulation acceleration history (3xM).
*   `refVals`: Reference values for redimensionalization.
*   `problemParams`: Dimensional constants for plotting limit lines.
*   `nonDimParams`: Non-dimensional constants.
*   `optimParams`: Constraint settings (for glideslope/pointing overlay).
*   `flag_thrustGotLimited`: Boolean for labeling saturated runs.
*   `optHistory`, `ICstates`: Re-optimization data (can be empty).
*   `betaParam`: For labeling.

### `statsPlotting.m`
**Purpose:** Generates statistical distribution plots for dispersion and Monte Carlo studies. Filters to converged cases (exitflag 1 or 2), then produces histograms and scatter plots for gamma, gamma2, tgo, fuel cost, and guidance coefficients.

**Inputs:**
*   `Results`: Struct array from `dispersionStudy.m` or `accelMonteCarlo.m`.
*   `beta`: Cost weighting (for labeling).
*   `optimizationParams`: Constraint settings (for labeling).

### `CoordinateFunctions/PDI2MCMF.m`
**Purpose:** Converts the PDI state from spherical coordinates (altitude, lat, lon, velocity magnitude, flight path angle, azimuth) to Moon-Centered Moon-Fixed Cartesian position and velocity vectors.

**Inputs:**
*   `altitude_km`: Altitude above surface (km).
*   `lonInitDeg`, `latInitDeg`: PDI longitude and latitude (deg).
*   `landingLonDeg`, `landingLatDeg`: Landing site coordinates (deg), used for heading calculation.
*   `inertialVel_mps`: Inertial velocity magnitude (m/s).
*   `flightPathAngleDeg`: Flight path angle (deg), angle above local horizontal.
*   `azimuth`: Heading angle (rad), clockwise from North.
*   `radMoon`: Moon radius (m).

**Outputs:**
*   `r_mcmf`: Position in MCMF (3x1, m).
*   `v_mcmf`: Velocity in MCMF (3x1, m/s).

### `CoordinateFunctions/ENU2MCMF.m`
**Purpose:** Rotates a vector from East-North-Up coordinates (anchored at a given lat/lon) into MCMF. For position vectors, also translates from the local ENU origin to the MCMF origin.

**Inputs:**
*   `enu`: Vector in ENU frame (3x1, non-dimensional).
*   `anchorLatDeg`: Anchor point latitude (deg).
*   `anchorLonDeg`: Anchor point longitude (deg).
*   `isPosition`: Boolean — if true, applies both rotation and translation; if false, rotation only (for velocity/acceleration vectors).

**Outputs:**
*   `mcmf`: Vector in MCMF frame (3x1).

### `CoordinateFunctions/MCMF2ENU.m`
**Purpose:** Rotates a vector from MCMF into East-North-Up coordinates anchored at a given lat/lon. For position vectors, subtracts the anchor point location before rotating.

**Inputs:**
*   `X`: Vector in MCMF frame (3x1 or 3xN).
*   `landingLatDeg`: Anchor point latitude (deg).
*   `landingLonDeg`: Anchor point longitude (deg).
*   `isPosition`: Boolean — if true, applies translation then rotation; if false, rotation only.
*   `isDim`: Boolean — if true, uses dimensional Moon radius (m); if false, uses non-dimensional.

**Outputs:**
*   `enu`: Vector in ENU frame (3x1 or 3xN).
*   `alt` (optional): Altitude above the Moon surface (scalar or 1xN), only computed when `isPosition` is true.

### `enuBasis.m`
**Purpose:** Computes the East, North, and Up unit vectors in MCMF coordinates for a given geodetic latitude and longitude. Used internally by all coordinate transform functions.

**Inputs:**
*   `lat`: Geodetic latitude (rad).
*   `lon`: Geodetic longitude (rad).

**Outputs:**
*   `E`: East unit vector (3x1).
*   `N`: North unit vector (3x1).
*   `U`: Up unit vector (3x1).