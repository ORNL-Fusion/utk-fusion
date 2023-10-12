"""
Author: Jeremy Mateja (With some functions taken from lim_plots authored
by Shawn Zamperini)
Email: jmateja@vols.utk.edu

This is my own version of lim plots. Mostly I did this so I can 
learn about the different stuff in the nc file and add my own kinds 
of plots, like overlaying the 3DLIM result on the LAMS plots.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import netCDF4
import matplotlib
#import pickle
#import sys
#from scipy.optimize import curve_fit
from scipy.signal import medfilt
from scipy.interpolate import interp1d
#from scipy.interpolate import griddata
#from matplotlib import colors


# Some plot properties to make them a bit nicer.
plt.ion()
#plt.rcParams['font.family'] = 'serif'
fontsize = 27
ms = 4
lw = 5
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the tableau20 RGBs to numbers between (0,1) since this is how mpl accepts them.

for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)
class LimPlots:
    """
    Learning about the use of classes. It seems to save having 
    to define a lot of the same variables and what not for
    different functions in the class. Not sure how necessary it is.
    """
    def __init__(self, ncpath): #potentially add combine_repeat_runs
        """
        Initiallized with only the nc file. No need for the lim and dat files

        ncpath: Path to the NetCDF file
        """
        self.ncpath = ncpath
        self.nc = netCDF4.Dataset(ncpath)


    def get_dep_array(self, num_runs=399):
        """
        This is designed to load the deposition arrays for the
        collector probes

        To-do:
        Shawn talks about adding option for repeat runs but idk
        how important that is
        """

        # You only need to load it in once, so check to see if it 
        # has been defined yet.
        try:
            self.dep_arr

        except AttributeError:
            # This creates the array for the deposition file
            dep_arr = np.array(self.nc.variables['NERODS3'][0] * -1)
            print(len(dep_arr))
        
        # Define dep_arr so next time you won't have to choose all of the
        # file locations.
        self.dep_arr = dep_arr

        return self.dep_arr
    
    def centerline(self, filter = True, log=False, fit_exp=False, plotnum=0, show_plot=False):
        """
        Probably the most used function. Plots the OTF and ITF deposition
        along the centerlines on the same plot. 

        log:        Option to make y axis a log scale
        fit_exp:    Do an exponential fit onto the data and get the lambdas

        To-Do:
        Maybe add option so that ITF/OTF is only over a certain range (like
        first 5 cm)
        """
        matplotlib.rc('font', size = 25)
        #The deposition array
        dep_arr = self.get_dep_array()

        # Location of each P bin, and its width
        ps     = np.array(self.nc.variables['PS'][:].data)
        pwids  = np.array(self.nc.variables['PWIDS'][:].data)

        # Array of poloidal locations (i.e. the center of each p bin)
        pol_locs = ps - pwids/2.0

        # Distance cell centers along surfaces (i.e. the radial locations)
        rad_locs = np.array(self.nc.variables['ODOUTS'][:].data)

        # Get the centerline index (or closest to it).
        cline = np.abs(pol_locs).min()

        # Index the deposition array at the centerline for plotting
        itf_x = rad_locs[np.where(rad_locs > 0.0)[0]]
        itf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs > 0.0)[0]]
        otf_x = rad_locs[np.where(rad_locs < 0.0)[0]] * -1
        otf_y = dep_arr[np.where(pol_locs == cline)[0], np.where(rad_locs < 0.0)[0]]
        
        #plotting stuff
        if plotnum == 0:
            if show_plot:
                fig = plt.figure()
                ax = fig.add_subplot(111)
        else:
            ax = self.master_fig.axes[plotnum-1]
        
        #set up filter arrays
        if filter:
            itf_y_filtered = medfilt(itf_y, 11)
            otf_y_filtered = medfilt(otf_y, 11)
        
        if show_plot:
            if log:
                if filter:
                    ax.semilogy(itf_x*100, itf_y_filtered, '-', label='ITF', ms=ms, color = tableau20[6])
                    ax.semilogy(otf_x*100, otf_y_filtered, '-', label='OTF', ms=ms, color = tableau20[8])
                else:
                    ax.semilogy(itf_x*100, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
                    ax.semilogy(otf_x*100, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])
            else:
                if filter:
                    ax.plot(itf_x*100, itf_y_filtered, '-', label='ITF', ms=ms, color=tableau20[6])
                    ax.plot(otf_x*100, otf_y_filtered, '-', label='OTF', ms=ms, color=tableau20[8])
                else:
                    ax.plot(itf_x*100, itf_y, '-', label='ITF', ms=ms, color=tableau20[6])
                    ax.plot(otf_x*100, otf_y, '-', label='OTF', ms=ms, color=tableau20[8])
            
            ax.legend(fontsize=fontsize)
            ax.set_xlabel('Distance along probe (cm)', fontsize=fontsize)
            ax.set_ylabel('Deposition (arbitrary units)', fontsize=fontsize)

        if plotnum == 0:
            if show_plot:
                fig.tight_layout()
                fig.show()
        
        return itf_x, otf_x, itf_y_filtered, otf_y_filtered
        

    def lams_overlay(self, absfac = 1.0, time = 4.0, shift = 0):

        #Getting the 3dlim arrays. I wonder if there is a better way of writing this.
        lim_itfx, lim_otfx, lim_itfy_filtered, lim_otfy_filtered = self.centerline()
        colors = matplotlib.cm.nipy_spectral(np.linspace(0.2,0.9,2))
        #Now we have to load in LAMS profiles
        cpexcelpath = "/mnt/c/Users/jmateja/Documents/SASVWcp/Messer_Sinclair_mastersheets/MCP01W.xlsx"
        xl = pd.ExcelFile(cpexcelpath)
        # Maybe make new excel spreadsheet to help streamline this more
        #probe = '01W'
        df_itf = pd.read_excel(xl, sheet_name = 'L')  ######
        df_otf = pd.read_excel(xl, sheet_name = 'R')  ######
        #df_itf = pd.read_excel(xl, sheet_name = 'MCPL02W_190424_190425')  ######
        #df_otf = pd.read_excel(xl, sheet_name = 'MCPR02W_190424_190425')  ######
        lams_itfx = df_itf['Distance along probe (cm)'].to_numpy()
        lams_otfx = df_otf['Distance along probe (cm)'].to_numpy()
        lams_itfy = df_itf['W areal density (1e15cm-2)'].to_numpy() * 1e19/2
        lams_otfy = df_otf['W areal density (1e15cm-2)'].to_numpy() * 1e19/2
        lams_itfy_error = df_itf['W areal density error (1e15cm-2)'].to_numpy() * 1e19/2
        lams_otfy_error = df_otf['W areal density error (1e15cm-2)'].to_numpy() * 1e19/2
        #Convert 3dlim profiles to r-rsep (Do this at later date)

        

        #Normalize the lams and 3DLIM profiles
        limtot = np.concatenate((lim_itfy_filtered, lim_otfy_filtered))
        lamstot = np.concatenate((lams_itfy, lams_otfy))

        lim_itfnorm = lim_itfy_filtered/max(limtot)
        lim_otfnorm = lim_otfy_filtered/max(limtot)

        lams_itfnorm = lams_itfy/max(lamstot)
        lams_otfnorm = lams_otfy/max(lamstot)
        lams_itfnormerror =  lams_itfy_error/max(lamstot)
        lams_otfnormerror =  lams_otfy_error/max(lamstot)

        #Plot 3dlim and lams results together
        #matplotlib.use('TKAgg')
        matplotlib.rc('font', size = 24)
        fig, ax = plt.subplots(figsize = (4, 3.5))
        plt.rcParams['xtick.labelsize']=18
        plt.rcParams['ytick.labelsize']=18
        x = [0, 2, 4, 6, 8]
        ticks = x

        ax.plot(lim_itfx*100-shift, lim_itfnorm, '-.', color='r', ms = ms, linewidth = 3, label = "3DLIM ITF")
        ax.plot(lim_otfx*100-shift, lim_otfnorm, '-.', color='b', ms = ms, linewidth = 3, label = "3DLIM OTF")
        #ax.plot(lams_itfx, lams_itfnorm, color = 'r', alpha = 0.7)
        ax.fill_between(lams_itfx, lams_itfnorm-lams_itfnormerror, lams_itfnorm+lams_itfnormerror,  
            facecolor = colors[1], alpha = 0.65, label = "LAMS ITF")
        #ax.plot(lams_otfx, lams_otfnorm, color = 'b', alpha = 0.7)
        ax.fill_between(lams_otfx, lams_otfnorm-lams_otfnormerror, lams_otfnorm+lams_otfnormerror,  
            facecolor = colors[0], alpha = 0.8, label = "LAMS OTF")
        ax.set_xlabel("Distance along probe [cm]", fontsize = 24)
        ax.set_ylabel("Normalized W Deposition", fontsize = 24)
        ax.set_ylim(bottom = 0.0, top = 1.1)
        ax.set_xlim(left = 0-shift, right = 10)
        ax.set_xticks(x)
        ax.set_xticklabels(ticks)
        ax.legend(loc = 'upper right', fontsize = 16)
        #ax.set_title("\"Low-Density\" Unfavorable Deposition Profiles", fontsize = 14, fontweight = "bold")

        plt.tight_layout()
        plt.show()

        #Find ABSFAC
        # absfac_text = f'{absfac:.2E}'
        # matplotlib.rc('font', size = 24)
        
        # fig, ax = plt.subplots(figsize = (10, 8))
        # x = [0, 2, 4, 6, 8]
        # ticks = x

        # ax.plot(lim_itfx*100-shift, absfac*lim_itfy_filtered, '-', color=tableau20[6], ms = ms, linewidth = 4, label = "3DLIM ITF")
        # ax.plot(lim_otfx*100-shift, absfac*lim_otfy_filtered, '-', color=tableau20[8], ms = ms, linewidth = 4, label = "3DLIM OTF")
        # ax.fill_between(lams_itfx, (lams_itfy-lams_itfy_error)/ time, (lams_itfy+lams_itfy_error) / time,  
            # facecolor = tableau20[6], alpha = 0.4, label = "LAMS ITF")
        # #ax.plot(lams_otfx, lams_otfnorm, color = 'b', alpha = 0.5)
        # ax.fill_between(lams_otfx, (lams_otfy-lams_otfy_error) / time, (lams_otfy+lams_otfy_error) / time,  
            # facecolor = tableau20[8], alpha = 0.5, label = "LAMS OTF")
        # ax.set_xlabel("Distance along probe [cm]", fontsize = 20, fontweight = "bold")
        # ax.set_ylabel("W Deposition (m-2 s-1)", fontsize = 20, fontweight = "bold")
        # ax.set_xlim(left = 0-shift, right = 10)
        # ax.set_xticks(x)
        # ax.text(6, absfac*2, "ABSFAC="+str(absfac_text))
        # ax.set_xticklabels(ticks)
        # ax.legend(loc = 'upper right', fontsize = 15)
        # ax.set_title("MCP 01W Deposition Profiles", fontsize = 20, fontweight = "bold")

        # plt.tight_layout()
        # plt.show()


        #Radial profiles of W density using ddlim3 after running 3dlim again without probe
        #plot 411 and 407 RCP together. Hopefully they are similar

        return
"""
    def plot_radial_density(self, pol_idx = 21, absfac = 1.0, charge = "all"):
        
        nc = self.nc
        xouts = nc['XOUTS'][:].data
        youts = nc['YOUTS'][:].data
        xwids = nc["XWIDS"][:].data
        mask = xwids != 0
        xouts = xouts[mask]
        xkeep_min = np.nonzero(xouts)[0].min()
        xkeep_max = np.nonzero(xouts)[0].max()+1
        ykeep_min = np.nonzero(youts)[0].min()
        ykeep_max = np.nonzero(youts)[0].max()+1
        xouts = xouts[xkeep_min:xkeep_max+1]

        try:
            yabsorb1a = float(nc["yabsorb1a"][:].data)
            yabsorb2a = float(nc["yabsorb2a"][:].data)
        except:
            yabsorb1a = -99
            yabsorb2a = 99
        if yabsorb1a > yabsorb2a:
            yabsorb_max = yabsorb1a
            yabsorb_min = yabsorb2a
        else:
            yabsorb_max = yabsorb2a
            yabsorb_min = yabsorb1a

        ykeep = np.where(np.logical_and(youts>=yabsorb_min, youts<=yabsorb_max))[0]
        youts = youts[ykeep]
        print(len(ykeep))
        Y, X = np.meshgrid(xouts, youts)

        if charge == "all":
            print("Summing across all charge states...")
            nz = nc["DDLIM3"][:].data[pol_idx, 1:, :, :].sum(axis=0)
        else:
            nz = nc["DDLIM3"][:].data[pol_idx, charge, :, :]
            nz = nz * absfac
        nz = nz[ykeep_min:ykeep_max, xkeep_min:xkeep_max]
        #print(nz)
        nz = nz[ykeep, :]
        mididx = np.argmin(np.abs(X[:, 0])) - 5
        nz = nz[mididx].data*absfac


        #Convert xouts to psin coordinates at some point
        # Read in excel file, use psin and distance along probe columns to convert

        cpexcelpath = "/mnt/c/Users/jmateja/Documents/SASVWcp/SASVW Updated CP Spreadsheets/MCP/xlsx/MCPL01W.xlsx"
        xl = pd.ExcelFile(cpexcelpath)

        df = pd.read_excel(xl)
        tip = df['R tip (cm)'].loc[df.index[0]]
        rad_locs = df['R (cm)'].to_numpy()
        psin = df['Psin Shot: 190422']

        f_psin = interp1d(rad_locs, psin, fill_value = "extrapolate")
        xlocs = (tip - Y[0])*100
        #ew_xs = xouts*-100 + tip
        new_psins = f_psin(xlocs)
        #gfile_path = "/mnt/c/Users/jmateja/Documents/SASVWlim/lim190422/190422_3000.pickle"
        #with open(gfile_path, "rb") as f:
        #    gfile = pickle.load(f)
        #R = gfile["R"]
        #Z = gfile["Z"]
        #Rs, Zs = np.meshgrid(R, Z)
        #psin = gfile["PSIRZ_NORM"]

        #lim_rzs = zip(rad_locs, np.full(len(rad_locs), -0.188))
        #print(list(lim_rzs))
        #lim_psins = griddata((Rs.flatten(), Zs.flatten()), psin.flatten(), list(lim_rzs))        
        #print(lim_psins)
        
        matplotlib.rc('font', size = 24)
        fig, ax = plt.subplots()
        ax.semilogy(new_psins, nz)
        ax.set_xlabel("Psin", fontsize = 20, fontweight = "bold")
        ax.set_ylabel("W Density (m-3)", fontsize = 20, fontweight = "bold")
        ax.set_title("190422 Radial W Density", fontsize = 20, fontweight = "bold")
        plt.show()

        
        
        return
"""

#ncpath = '/mnt/c/research/3dlim/190422/Results/190422_3000_MCP_v65.nc'
#lp = LimPlots(ncpath)
#lp.lams_overlay()