import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# Add the target module directory (e.g., ../analysis) to sys.path
analysis_dir = os.path.join(script_dir, '..', 'analysis')
sys.path.insert(0, os.path.abspath(analysis_dir))

from read_data import MassInflationData

run_dir = os.path.join(script_dir, 'RN_config00')

reader = MassInflationData(run_dir)

U,V,R = reader.get_data('fields.dat', col=1)

_,_,drdu = reader.get_data('derivatives.dat', col=1)
_,_,drdv = reader.get_data('derivatives.dat', col=2)


plt.figure()
plt.contour(V, U, R, levels=50, cmap='viridis')
plt.xlabel('v')
plt.ylabel('u')
plt.colorbar()
plt.title('r(u,v)')
plt.tight_layout()
plt.show()
