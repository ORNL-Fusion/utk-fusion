"""
Trying to plot the far SOL radial impurity density from a 3DLIM run with no probe
"""

import matplotlib.pyplot as plt
import LimPlots
import netCDF4
import matplotlib
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


def get_lim_data(lim_path, lim_absfac):
    lp = LimPlots.LimPlots(lim_path)
    lpdata = lp.plot_par_rad("nz", 21, charge="all")

    # Let's not choose right in the middle (Y = 0) since flows are zero there, and it actually leads to distracting
    # valleys in the radial profile due to no impurities there.
    mididx = np.argmin(np.abs(lpdata["X"][:, 0])) - 5
    lim_nz = lpdata["Z"][mididx].data * lim_absfac

    # This value is chosen in lim_a2_again.py. It is that needed to match the experimental deposition.
    cpexcelpath = "/mnt/c/Users/jmateja/Documents/SASVWcp/SASVW Updated CP Spreadsheets/MCP/xlsx/MCPL01W.xlsx"
    xl = pd.ExcelFile(cpexcelpath)

    df = pd.read_excel(xl)
    tip = df['R tip (cm)'].loc[df.index[0]]
    rad_locs = df['R (cm)'].to_numpy()
    psin = df['Psin Shot: 190422'].to_numpy()

    f_psin = interp1d(rad_locs, psin, fill_value = "extrapolate")
    xlocs = (tip - lpdata["Y"][0]*100)
    new_psins = f_psin(xlocs)
    return new_psins, lim_nz

def plot_radial_profile(ncpath, absfac):
    
    psin, nz = get_lim_data(ncpath, absfac)

    matplotlib.rc('font', size = 24)
    fig, ax = plt.subplots()
    ax.semilogy(psin, nz)
    ax.set_xlabel("Psin", fontsize = 20, fontweight = "bold")
    ax.set_ylabel("W Density (m-3)", fontsize = 20, fontweight = "bold")
    ax.set_title("190422 Radial W Density", fontsize = 20, fontweight = "bold")
    plt.show()
    
    return