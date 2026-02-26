from fenics import *
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Parameters (change these to play around)
# -------------------------
L = 1.0          # length of the 1D domain (basically from x=0 to x=L)
N = 200          # how fine the mesh is (more = smoother result but slower)

D = 1e-2         # diffusion coefficient (how fast stuff spreads)
dt = 0.01        # time step size
T = 0.6          # total simulation time

C_left = 1.0     # concentration fixed at the left boundary (x=0)

# Degradation rates
# k = 0 means no degradation
# larger k means stuff disappears faster
k_susc = 0.0     # Susceptible case → no degradation
k_res  = 1.0     # Resistant case → degradation happening

# -------------------------
# Mesh and function space
# -------------------------

# Create a 1D mesh from 0 to L split into N pieces
mesh = IntervalMesh(N, 0.0, L)

# Use linear finite elements (continuous, degree 1)
V = FunctionSpace(mesh, "CG", 1)

# -------------------------
# Boundary condition at x=0
# -------------------------

# This function tells FEniCS which boundary we mean
def left_boundary(x, on_boundary):
    # Only apply on the boundary AND where x is basically 0
    return on_boundary and near(x[0], 0.0)

# Fix concentration at x=0 to C_left
bc = DirichletBC(V, Constant(C_left), left_boundary)

# -------------------------
# Initial condition
# -------------------------

# Start with zero concentration everywhere
# So the source at x=0 will slowly fill the domain
u0 = Function(V)
u0.assign(Constant(0.0))

# Trial and test functions (standard FEM setup)
u = TrialFunction(V)
v = TestFunction(V)

# -------------------------
# Function that runs the simulation
# -------------------------

def run_simulation(k_deg, label="case"):
    """
    Runs diffusion + degradation simulation for a given k_deg.
    Returns x values and final concentration profile.
    """

    # u_n = solution from previous time step
    u_n = Function(V)
    u_n.assign(u0)

    # u_sol = solution at current time step
    u_sol = Function(V)

    # Weak form using implicit Euler:
    #
    # (u - u_n)/dt = D*u_xx - k*u
    #
    # After rearranging and writing in weak form, we get:
    #
    # u*v + dt*(D*grad(u)·grad(v) + k*u*v) = u_n*v
    #
    # Everything on left = matrix
    # Everything on right = known from previous step

    a = (u*v + dt*(D*dot(grad(u), grad(v)) + k_deg*u*v))*dx
    L_form = (u_n*v)*dx

    t = 0.0
    while t < T - 1e-12:
        t += dt

        # Solve the linear system for this time step
        solve(a == L_form, u_sol, bc)

        # Update for next step
        u_n.assign(u_sol)

    # Get coordinates of nodes
    x = V.tabulate_dof_coordinates().reshape(-1)

    # Get concentration values
    c = u_sol.vector().get_local()

    # Sort them so the plot isn't weirdly scrambled
    idx = np.argsort(x)
    return x[idx], c[idx]

# -------------------------
# Run both cases
# -------------------------

x_s, c_s = run_simulation(k_susc, label="Susceptible")
x_r, c_r = run_simulation(k_res,  label="Resistant")

# -------------------------
# Plot results
# -------------------------

plt.figure()
plt.plot(x_s, c_s, label=f"Susceptible (k={k_susc})")
plt.plot(x_r, c_r, label=f"Resistant (k={k_res})")
plt.xlabel("x")
plt.ylabel("Concentration")
plt.title(f"1D diffusion (+ degradation) at t={T:.3f}")
plt.legend()
plt.tight_layout()
plt.show()