# 1D Diffusion with Degradation (FEniCS)

This project simulates **1D diffusion of a concentration** using the Finite Element Method (FEniCS).

It compares two cases:
- **Susceptible**  diffusion only (no degradation)
- **Resistant**  diffusion + first-order degradation

The goal is to see how degradation changes the final concentration profile.

## Model

We solve the PDE:

∂u/∂t = D ∂²u/∂x² − k u

Where:
- u(x,t) = concentration
- D = diffusion coefficient
- k = degradation rate
- Domain: x ∈ [0, L]

### Boundary condition
- At x = 0 → fixed concentration (Dirichlet condition)
- Right boundary → natural zero-flux condition

### Initial condition
- u(x,0) = 0 everywhere


## What the Code Does

1. Builds a 1D mesh
2. Defines a CG1 (linear) finite element space
3. Uses implicit Euler for time stepping
4. Runs two simulations:
   - k = 0 (Susceptible)
   - k = 1 (Resistant)
5. Plots the final concentration profile at time T


## How to Run

Make sure you have FEniCS installed.

Then run:

python your_script_name.py

A matplotlib window will open showing both concentration profiles.


## Expected Result

- The Susceptible case diffuses further into the domain.
- The Resistant case has lower concentration due to degradation.
- Increasing k increases suppression of concentration.


## Notes

You can change parameters at the top of the script:
- D (diffusion coefficient)
- k_res (degradation rate)
- dt (time step)
- T (final time)
- N (mesh resolution)

Higher N → smoother curve but slower computation.
Smaller dt → more accurate but slower simulation.
