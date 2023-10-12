"""
Generates a plot with multiple pdfs to more easily show effects of DIVIMP parameters such as
the cross field diffusion coefficient and drifts.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.signal import savgol_filter

#Initiate flux tube parameters
z_fluxtube     = -0.188
r_fluxtube     = 2.3175 - 0.03
show_plot      = True
reverse        = False
smooth         = True
smooth_window  = 15

#Set up nc paths for runs you want to compare. (Is there a more efficient way to do this?)
divimp_nc_path1 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/Million_Particle_Runs/d3d-187111-inj-42-mil.nc"
divimp_nc_path2 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/Million_Particle_Runs/d3d-187111-inj-44-mil.nc"
divimp_nc_path3 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/Million_Particle_Runs/d3d-187111-inj-39-mil.nc"

divimp_nc_paths = [divimp_nc_path1, divimp_nc_path2, divimp_nc_path3]

sys.path.append(os.getcwd().split("lim")[0] + "oedge")
import oedge_plots

#hsv = plt.get_cmap('hsv')
#colors = iter(hsv(np.linspace(0,100,1000)))
# First load the DIVIMP run.
for divimp_nc_path in divimp_nc_paths:
    op = oedge_plots.OedgePlots(divimp_nc_path)
    cell = op.find_ring_knot(r_fluxtube, z_fluxtube, return_cell=True)
    ring, knot = op.find_ring_knot(r_fluxtube, z_fluxtube)

    # Grab the impurity density and cell volume along our selected ring.
    s, imp = op.along_ring(int(ring), "DDLIMS", charge="all", plot_it=False)
    s, vol = op.along_ring(int(ring), "KVOLS", plot_it=False)
    tot_imp = imp * vol
    if smooth:
        if type(smooth_window) == type(None):
            smooth_window = int(len(tot_imp)/15)
            if smooth_window % 2 == 0:
                smooth_window += 1
        tot_imp_smooth = savgol_filter(tot_imp, smooth_window, 2)
    if reverse:

        # Need to offset the S values by the min otherwise the S values will
        # not be the same.
        s   = s.min() + s.max() - s[::-1]
        tot_imp = tot_imp[::-1]
        if smooth:
            tot_imp_smooth = tot_imp_smooth[::-1]

    #fig, ax = plt.subplots()
    #plt.plot(s, tot_imp, color=next(colors))
    if smooth:
       norm = tot_imp_smooth/max(tot_imp_smooth) 
       plt.plot(s, norm)
#plt.legend(["Drift Mult. = 0", "Drift Mult. = 0.5", "Drift Mult. = 1.0"], loc = "upper right")
#plt.title("OSP Drift Multiplier Comparison with Dperp = 1.0 m^2/s")
plt.title("ISP Cross Field Diffusion Comparison with Drift Multiplier = 1.0")
plt.legend(["Dperp = 0.3 m^2/s", "Dperp = 0.6 m^2/s", "Dperp = 1.0 m^2/s"], loc = "upper right")
plt.xlabel("Distance from target")
plt.ylabel("Nomralized number of impurity ions")
#fig.tight_layout()
plt.show()