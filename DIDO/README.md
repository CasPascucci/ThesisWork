# DIDO Optimal Control - Lunar Powered Descent

Solves the lunar powered descent problem using DIDO's pseudospectral method to produce an optimal solution for comparison against FP2DG.

## Optimal Control Problem

- **States (7):** r_x, r_y, r_z, v_x, v_y, v_z, m
- **Controls (3):** aT_x, aT_y, aT_z (commanded acceleration)
- **Cost:** J = integral of [beta*||aT|| + (1-beta)*||aT||^2] dt
- **Dynamics:** dr/dt = v, dv/dt = aT + g(r), dm/dt = -||aT||*m / isp
- **Path constraint:** minThrust <= m*||aT|| <= maxThrust
- **Boundary conditions:** r(0), v(0), m(0) fixed; r(tf), v(tf) fixed; m(tf) free; tf free

All quantities non-dimensionalized with L_ref=10000m, same as `getParams.m`.

## Files

### `FP2DG_problem.m`
**Purpose:** Main entry point. Sets up the DIDO problem with the same parameters as `userTest.m`, runs the solver, post-processes and plots results. Includes V&V via Hamiltonian check and ODE45 propagation.

### `FP2DG_preamble.m`
**Purpose:** Extracts named variables (r, v, m, aT, t) from DIDO's `primal` structure. Called by dynamics, cost, and path functions.

### `FP2DG_dynamics.m`
**Purpose:** Equations of motion. Returns 7xN state derivatives [drdt; dvdt; dmdt] at all collocation nodes.

### `FP2DG_cost.m`
**Purpose:** Objective function. Returns endpoint cost (0) and running cost (beta*||aT|| + (1-beta)*||aT||^2) at each node.

### `FP2DG_events.m`
**Purpose:** Boundary conditions. Returns 13-element vector: initial r, v, m (7) and final r, v (6). All set to equality in bounds.

### `FP2DG_path.m`
**Purpose:** Path constraint. Returns thrust magnitude h = m*||aT|| at each node, bounded by engine thrust limits.
