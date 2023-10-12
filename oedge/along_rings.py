"""
Along ring profiles for multiple rings
"""

import oedge_plots
import matplotlib.pyplot as plt
from scipy.signal import medfilt
import matplotlib as mpl


def multiple_rings(ncpath,  ymax = 1e14):
    #ncpath = "/mnt/c/Users/jmateja/Documents/DIVIMP Related/SAS-VW/190424/424_inj_files/d3d-190424-inj-10_nd_mil.nc"
    
    op = oedge_plots.OedgePlots(ncpath)
    mpl.rc('font', size = 20)
    fig, ax = plt.subplots()
    for ring in range(47, 51): # Change this
        s, nz = op.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
        nz = medfilt(nz, 5) # Might also want to play around with this
        ax.plot(s, nz, label="Ring: "+str(ring))
    ax.set_xlabel("Distance from Outer target (m)") # Which target depends on configuration.
    ax.set_ylabel("Impurity Density (m-3)") # Will be m-3 only if ABSFAC was supplied.
    #ax.set_title(str(launch)+", Drifts off, Dperp = 0.3 m^2/s")
    ax.set_title("Along Ring Impurity Density")
    ax.set_ylim(0, ymax)
    ax.legend()
    fig.tight_layout()
    fig.show()
    #print(s)
    #print(nz)
    return

def launch_comp(ring, ymin = 0, ymax = 5e14):
    ncpath1 = "/mnt/c/Users/jmateja/Documents/DIVIMP Related/SAS-VW/190424/424_inj_files/drifts_launches/d3d-190424-inj-sasv_dp1.nc"
    ncpath2 = "/mnt/c/Users/jmateja/Documents/DIVIMP Related/SAS-VW/190424/424_inj_files/drifts_launches/d3d-190424-inj-sasv_dp2.nc"
    ncpath3 = "/mnt/c/Users/jmateja/Documents/DIVIMP Related/SAS-VW/190424/424_inj_files/drifts_launches/d3d-190424-inj-sasv_dp3.nc"
    ncpath4 = "/mnt/c/Users/jmateja/Documents/DIVIMP Related/SAS-VW/190424/424_inj_files/drifts_launches/d3d-190424-inj-sasv_dp4.nc"
    #ncpath5 = "/mnt/c/Users/jmateja/Documents/DIVIMP Related/SAS-VW/190424/424_inj_files/d3d-190424-inj-13_nd_mil.nc"
    #ncpath6 = "/mnt/c/Users/jmateja/Documents/DIVIMP Related/SAS-VW/190424/424_inj_files/d3d-190424-inj-14_nd_mil.nc"
    
    op1 = oedge_plots.OedgePlots(ncpath1)
    op2 = oedge_plots.OedgePlots(ncpath2)
    op3 = oedge_plots.OedgePlots(ncpath3)
    op4 = oedge_plots.OedgePlots(ncpath4)
    #op5 = oedge_plots.OedgePlots(ncpath5)
    #op6 = oedge_plots.OedgePlots(ncpath6)
    mpl.rc('font', size = 16)
    fig, ax = plt.subplots()
    s1, nz1 = op1.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    s2, nz2 = op2.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    s3, nz3 = op3.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    s4, nz4 = op4.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    #s5, nz5 = op5.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    #s6, nz6 = op6.along_ring(ring, "DDLIMS", charge="all", plot_it=False)
    nz1 = medfilt(nz1, 5) # Might also want to play around with this
    nz2 = medfilt(nz2, 5)
    nz3 = medfilt(nz3, 5)
    nz4 = medfilt(nz4, 5)
    #nz5 = medfilt(nz5, 5)
    #nz6 = medfilt(nz6, 5)
    ax.plot(s1, nz1, label="Progressive Angle", linewidth = 2)
    ax.plot(s2, nz2, label="Inner Target", linewidth = 2)
    ax.plot(s3, nz3, label="Inner Wall", linewidth = 2)
    ax.plot(s4, nz4, label="Leading Edge", linewidth = 2)
    #ax.plot(s5, nz5, label="Outside slot with leading edge")
    #ax.plot(s6, nz6, label="Leading edge only")
    ax.set_xlabel("Distance from Outer target (m)") # Which target depends on configuration.
    ax.set_ylabel("Impurity Density (m-3)") # Will be m-3 only if ABSFAC was supplied.
    ax.set_title("Along Ring Impurity Density For Ring:" + str(ring))
    ax.set_ylim([ymin, ymax])
    ax.legend()
    fig.tight_layout()
    fig.show()
    #print(s)
    #print(nz)
    return

def adding():
    ncpath1 = "/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/d3d-187111-inj-30.nc"
    ncpath2 = "/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/d3d-187111-inj-31.nc"
    ncpath3 = "/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/d3d-187111-inj-32.nc"
   
    op1 = oedge_plots.OedgePlots(ncpath1)
    op2 = oedge_plots.OedgePlots(ncpath2)
    op3 = oedge_plots.OedgePlots(ncpath3)
    
    s1, nz1 = op1.along_ring(56, "DDLIMS", charge="all", plot_it=False)
    norm_nz1 = nz1/max(nz1)
    outer = savgol_filter(norm_nz1, window_length = 15, polyorder = 3)
    
    s2,nz2 = op2.along_ring(56, "DDLIMS", charge="all", plot_it=False)
    norm_nz2 = nz2/max(nz2)
    inner = savgol_filter(norm_nz2, window_length = 15, polyorder = 3)
    
    s3,nz3 = op3.along_ring(56, "DDLIMS", charge="all", plot_it=False)
    norm_nz3 = nz3/max(nz3)
    wall = savgol_filter(norm_nz3, window_length = 15, polyorder = 3)
    
    total = 0*wall + 1.0*inner
    
    plt.plot(s1, total)
    plt.title("Normalized Impurity Distribution vs S")
    plt.xlabel("S (m)")
    plt.ylabel("Normalized Impurity Distribution")
    plt.text(40, 0.9, "Ring "+str(56))

    plt.show()
    return

def normalize():
    ncpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/d3d-187111-inj-34.nc"

    op = oedge_plots.OedgePlots(ncpath)

    s, imp = op.along_ring(56, "DDLIMS", charge ="all", plot_it=False)
    
    norm_imp = imp/max(imp)

    plt.plot(s, norm_imp)
    plt.show()
    
    
    return s, norm_imp

def charges(ring, max_charge = 10):

    ncpath1 = "/mnt/c/Users/jmateja/Documents/DIVIMP Related/SAS-VW/190424/424_inj_files/drifts_launches/d3d-190424-inj-14d.nc"
    
    op1 = oedge_plots.OedgePlots(ncpath1)
    
    mpl.rc('font', size = 16)
    fig, ax = plt.subplots()
    
    for i in range(0, max_charge + 1):
        s1, nz1 = op1.along_ring(ring, "DDLIMS", charge=i, plot_it=False)
        nz1 = medfilt(nz1, 5) # Might also want to play around with this
        ax.plot(s1, nz1, label=i, linewidth = 2)
    
    ax.set_xlabel("Distance from Outer target (m)") # Which target depends on configuration.
    ax.set_ylabel("Impurity Density (m-3)") # Will be m-3 only if ABSFAC was supplied.
    ax.set_title("Charge State Along Ring Impurity Density For Ring:" + str(ring))
    #ax.set_ylim([ymin, ymax])
    ax.legend()
    fig.tight_layout()
    fig.show()
    return