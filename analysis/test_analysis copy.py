
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
run_dir = os.path.join(workspace_root, 'runs', 'RN_config_Ana')
reader = MassInflationData(run_dir)



# %%
U, V, R = reader['r']

plt.figure()
plt.contour(V, U, R, levels=100, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
plt.ylim(22, 40)
plt.colorbar()
plt.title('r(u,v)')
plt.tight_layout()
plt.show()

levels = np.linspace(1.45, 1.53, 50)

plt.figure()
contour = plt.contour(V, U, R, levels=levels, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
plt.ylim(22, 40)
cbar = plt.colorbar(contour)
plt.title('r(u,v)')
plt.tight_layout()
plt.show()

levels = np.linspace(0.53, 0.66, 50)
plt.figure()
contour = plt.contour(V, U, R, levels=levels, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
plt.ylim(10, 30)
plt.xlim(5, 30)
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





U, V, M = reader['mass']
M_abs = np.abs(M)
plt.figure()
cf = plt.contour(V, U, np.log10(M_abs), levels=50, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
plt.ylim(22, 40)
plt.xlim(5, 30)
plt.colorbar(cf)
plt.title('log10 |m(u,v)|')
plt.tight_layout()
plt.show()



# %%
U, V, Ricci = reader['Ricci']

Ricci_abs = np.abs(Ricci)
Ricci_safe = np.maximum(Ricci_abs, 1e-16)  # avoid log(0)

plt.figure()
cf = plt.contour(V, U, np.log10(Ricci_safe), levels=80, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
plt.ylim(10, 40)
plt.colorbar(cf)
plt.title('log10 |R(u,v)|')
plt.tight_layout()
plt.show()

# %%

U, V, M = reader['mass']

u_indices = [-100,-20,-50, -75,-140]  # try different slices

plt.figure()

for ui in u_indices:
    v_vals = V[:, ui]
    m_vals = M[:, ui]
    m_safe = np.maximum(np.abs(m_vals), 1e-16)
    label = f'u = {U[0, ui]:.3f}'
    plt.plot(v_vals, np.log(m_safe),label=label)
plt.xlabel('v')
plt.ylabel('log |m|')
plt.title('Mass growth for different u')
plt.legend()
plt.tight_layout()
plt.show()
# %%




U, V, M = reader['mass']

u_indices = [-100,-20,-50, -75,-140]  # try different slices

plt.figure()

for ui in u_indices:
    v_vals = V[:, ui]
    m_vals = M[:, ui]
    m_safe = np.maximum(np.abs(m_vals), 1e-16)
    label = f'u = {U[0, ui]:.3f}'
    plt.plot(v_vals, m_safe,label=label)
plt.xlabel('v')
plt.ylabel('|m|')
plt.title('Mass growth for different u')
plt.legend()
plt.tight_layout()
plt.show()
# %%
