
# %%
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# %%
# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# Add the workspace root to sys.path for package imports
workspace_root = os.path.abspath(os.path.join(script_dir, '..'))
if workspace_root not in sys.path:
	sys.path.insert(0, workspace_root)

from analysis.read_data import MassInflationData

# %%
run_dir = os.path.join(script_dir, 'RN_config00')
reader = MassInflationData(run_dir)

# %%
U, V, R = reader['r']

plt.figure()
plt.contour(V, U, R, levels=50, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
plt.colorbar()
plt.title('r(u,v)')
plt.tight_layout()
plt.show()

levels = np.linspace(1.2, 1.4, 10)

plt.figure()
contour = plt.contour(V, U, R, levels=levels, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
cbar = plt.colorbar(contour)
plt.title('r(u,v)')
plt.tight_layout()
plt.show()

# %%
U, V, f = reader['drdv']

f_abs = np.abs(f)

# Overlay: contour where |drdv| ~ 0
epsilon = 1e-5  # Adjust as needed for your data's noise floor
plt.figure()
cf = plt.contourf(V, U, np.log10(f_abs), levels=50, cmap='viridis')
# Overlay the zero contour
contour_zero = plt.contour(V, U, f_abs, levels=[epsilon], colors='red', linewidths=2, linestyles='dashed')
plt.xlabel('v')
plt.ylabel('u')
plt.colorbar(cf)
plt.title(f'drdv(u,v) with |drdv|={epsilon} contour')
plt.tight_layout()
plt.show()

# %%
U, V, f = reader['Guu']

f_abs = np.abs(f)

plt.figure()
cf = plt.contourf(V, U, np.log10(f_abs), levels=50, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
plt.colorbar(cf)
plt.title(f'Guu(u,v)')
plt.tight_layout()
plt.show()

# %%
