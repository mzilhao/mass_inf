import numpy as np

with open('data.dat') as f:
    metadata = [line for line in f if line.startswith('#') and 'Simulation' in line]
    v_slices = [line for line in f if 'v =' in line]

data = np.loadtxt('data.dat')  # Shape: (N_points, 7) → [u, r, phi, sigma, mass, drdv, ricci]

# Extract v values from comment lines and reshape into grid if needed

