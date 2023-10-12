"""
Overlays centerline plots on same graph to show comparison between different Dperp values
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import savgol_filter
import lim_plots
from LimWallToolkit import LimWallToolkit

def mul_pdfs():

    ncpath1 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/Million_Particle_Runs/d3d-187111-inj-42-mil.nc"
    ncpath2 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/Million_Particle_Runs/d3d-187111-inj-44-mil.nc"
    ncpath3 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/Million_Particle_Runs/d3d-187111-inj-39-mil.nc"

    z_fluxtube     = -0.188
    r_fluxtube     = 2.3175 - 0.03
    show_plot      = False
    reverse        = False
    smooth         = True
    smooth_window  = 15
    lwt = LimWallToolkit()
    ncpath_array = [ncpath1, ncpath2, ncpath3]
    dperp = ['Dperp = 0.3', 'Dperp = 0.6', 'Dperp = 1.0']
    mpl.rc('font', size = 24)
    fig, ax = plt.subplots()
    i = 0
    for ncpath in ncpath_array:
        prob_dict = lwt.divimp_prob_dist(divimp_nc_path=ncpath, r_fluxtube=r_fluxtube, z_fluxtube=z_fluxtube, show_plot=show_plot, 
        reverse=reverse, smooth=smooth, smooth_window=smooth_window)
        s = prob_dict['s']
        tot_imp = prob_dict['tot_imp_smooth']
        norm = tot_imp/max(tot_imp)
        ax.plot(s, norm, label = dperp[i],linewidth = 4.2)
        i = i+1
    ax.set_title("Interface Profile", fontsize = 25, fontweight = "bold")
    ax.set_xlabel("Distance from target", fontsize = 25, fontweight = "bold")
    ax.set_ylabel("Normalized # of impurity ions", fontsize = 25, fontweight = "bold")
    ax.legend(loc = 'upper left', fontsize = 23)
    fig.tight_layout()
    plt.show()

    return


def mul_centerline():
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)

    ncpath1 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/lim/Methane Stuff/lim187111/187111_lim_4400_DCP_inj42.nc"
    ncpath2 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/lim/Methane Stuff/lim187111/187111_lim_4400_DCP_inj44_10mil_updated.nc"
    ncpath3 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/lim/Methane Stuff/lim187111/187111_lim_4400_DCP_inj39.nc"
    fontsize = 25
    ms = 2
    ncpath_array = [ncpath1, ncpath2, ncpath3]
    mpl.rc('font', size = 24)
    fig, ax = plt.subplots()
    dperp = ['Dperp = 0.3', 'Dperp = 0.6', 'Dperp = 1.0']
    i = 0
    for ncpath in ncpath_array:
        lp = lim_plots.LimPlots(ncpath)
        cp_dict = lp.centerline(log = True, show_plot = False)

        itf_x = cp_dict['itf_x']
        otf_x = cp_dict['otf_x']
        itf_y = cp_dict['itf_y_filtered']
        otf_y = cp_dict['otf_y_filtered']

        ax.semilogy(itf_x*100, itf_y, '-', label=dperp[i], ms=ms, color=tableau20[i+2])
        ax.semilogy(otf_x*100, otf_y, '-', ms=ms, color=tableau20[i+2])
        
        i = i+1
    ax.legend(fontsize=fontsize)
    ax.set_xlabel('Distance along probe (cm)', fontsize=fontsize, fontweight = "bold")
    ax.set_ylabel('Deposition (arbitrary units)', fontsize=fontsize, fontweight = "bold")
    ax.set_title("MCP 3DLIM Simulated Deposition", fontsize = 25, fontweight = "bold")   
    fig.tight_layout()
    plt.show()



def combo():

    ncpath1 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/Million_Particle_Runs/d3d-187111-inj-34-mil.nc"
    ncpath2 = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/Million_Particle_Runs/d3d-187111-inj-39-mil.nc"

    z_fluxtube     = -0.188
    r_fluxtube     = 2.3175 - 0.03
    show_plot      = False
    reverse        = False
    smooth         = True
    smooth_window  = 15
    lwt = LimWallToolkit()
    prob_dict_osp = lwt.divimp_prob_dist(divimp_nc_path=ncpath1, r_fluxtube=r_fluxtube, z_fluxtube=z_fluxtube, show_plot=show_plot, 
        reverse=reverse, smooth=smooth, smooth_window=smooth_window)
    s = prob_dict_osp['s']
    tot_imp_osp = prob_dict_osp['tot_imp_smooth']
    prob_dict_isp = lwt.divimp_prob_dist(divimp_nc_path=ncpath2, r_fluxtube=r_fluxtube, z_fluxtube=z_fluxtube, show_plot=show_plot, 
        reverse=reverse, smooth=smooth, smooth_window=smooth_window)
    tot_imp_isp = prob_dict_isp['tot_imp_smooth']
    
    osp_norm = tot_imp_osp / max(tot_imp_osp)
    isp_norm = tot_imp_isp / max(tot_imp_isp)
    
    combo = 0.3 * osp_norm + 0.7 * isp_norm
    
    mpl.rc('font', size = 25)
    fig, ax = plt.subplots()

    ax.plot(s, combo,linewidth = 4.2, color = "tab:red")
    ax.set_title("Interface Profile", fontsize = 27, fontweight = "bold")
    ax.set_xlabel("Distance from target", fontsize = 27, fontweight = "bold")
    ax.set_ylabel("Normalized # of impurity ions", fontsize = 27, fontweight = "bold")
    fig.tight_layout()
    plt.show()
    for i in range(0, len(combo)):
        print("{:<6.2f} {:.3e}".format(s[i], combo[i]))
    return




