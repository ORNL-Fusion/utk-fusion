# -*- coding: utf-8 -*-
"""
Created on Thu May 19 08:55:48 2022

@author: jmateja
"""

import oedge_plots
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter



ncpath = "/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/d3d-187111-inj-39.nc"
   
op = oedge_plots.OedgePlots(ncpath)
    
fig, ax = plt.subplots()
for ring in range(56, 57): # Change this
    s, nz = op.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    nz = savgol_filter(nz, window_length = 15, polyorder = 3) # Might also want to play around with this
    ax.plot(s, nz, label=ring)
ax.set_xlabel("Distance from Outer target (m)") # Which target depends on configuration.
ax.set_ylabel("Impurity Density (m-3)") # Will be m-3 only if ABSFAC was supplied.
ax.set_title("OSP, Drifts off, Dperp = 0.3 m^2/s")
fig.tight_layout()
fig.show()