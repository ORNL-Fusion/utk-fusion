"""
Program: updated_cp_plots.py
Author: Jeremy Mateja
Date: Aug 17, 2022

This program creates plots of the collector probe data as a function of both distance along probe and r-rsep
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import scipy.interpolate as scinter
import scipy.integrate as scinteg
import math
# Be mindful of which side of the CP you are trying to plot


#df = pd.read_excel(cpexcelpath)

def comparison_plot():
    
    fig, ax = plt.subplots(nrows = 2, ncols = 2, sharex= 'col', sharey = 'row', figsize=(12, 6))
    
    x = [4.5, 5, 5.5, 6, 6.5]
    ticks = x
    for j in range(0,2):    
        for k in range(0,2):
            if k == 0 and j == 0:
                #cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_sinclair_mastersheets/fwd bt dens right.xlsx"
                cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_sinclair_mastersheets/dcp unfav otf scan.xlsx"
                shots = [190423, 190425, 190426, 190492]
            elif k ==1 and j == 0:
                #cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_sinclair_mastersheets/fwd bt dens left.xlsx"
                cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_sinclair_mastersheets/dcp unfav itf scan.xlsx"
                shots = [190423, 190425, 190426, 190492]
            elif k == 0 and j == 1:
                #cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_sinclair_mastersheets/rev bt dens left.xlsx"
                cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_sinclair_mastersheets/dcp fav otf scan.xlsx"
                shots = [190454, 190456, 190461, 190482]
            else:
                #cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_sinclair_mastersheets/rev bt dens right.xlsx"
                cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_sinclair_mastersheets/dcp fav itf scan.xlsx"
                shots = [190454, 190456, 190461, 190482]
            
            xl = pd.ExcelFile(cpexcelpath)
            i = 0
            for sheet_name in xl.sheet_names:
                df = pd.read_excel(xl, sheet_name = sheet_name)
                shot = shots[i]
                r_rsep = df['R-Rsep OMP (cm) Shot:'+str(shot)].to_numpy()
                #psin = df["Psin"].to_numpy()
                if shot == 190454:
                    w_arealdens = df["W areal density (1e15cm-2)"].to_numpy()/3
                else:
                    w_arealdens = df["W areal density (1e15cm-2)"].to_numpy()/2
                ax[k, j].scatter(r_rsep, w_arealdens, label = "Shot:" + " " +str(shot)+", Line Average Density: "+sheet_name, s = 10, alpha = 0.3)
                ax[k, j].set_xlabel("R-Rsep OMP [cm]", fontsize = 12, fontweight = "bold")
                ax[k, j].set_ylabel("W Areal Density [1e15 cm-2]", fontsize = 12, fontweight = "bold")
                ax[k, j].set_ylim(bottom = 0.0, top = 0.04)
                ax[k, j].set_xlim(left = 4.0, right = 7.0)
                ax[k,j].set_xticks(x)
                ax[k,j].set_xticklabels(ticks)
                average = []
                newrs = []
                window = 5
                for ind in range(window-1, len(w_arealdens) - window, 5):
                    average.append(np.mean(w_arealdens[ind-window:ind+window]))
                    newrs.append(r_rsep[ind])
                ax[k,j].plot(newrs, average, label = "Shot "+str(shot)+ " averaged profile", linewidth = 4.5)
                ax[k,j].legend(loc = 'upper left', fontsize = 7)
                i = i + 1
        #norm = w_arealdens/max(w_arealdens)
        #ax1.plot(newrs, average, label = sheet_name+' Averaged', linewidth = 4)

    ax[0, 0].set_title("Unfavorable Bt L-mode CP Profiles (Pinj = 2.20 MW)", fontsize = 13, fontweight = "bold")
    ax[0, 1].set_title("Favorable Bt L-mode CP Profiles (Pinj = 1.25 MW)", fontsize = 13, fontweight = "bold")
    ax[0,0].text(4.1, 0.022, "(a)", fontsize = 14, fontweight = "bold")
    ax[0,1].text(4.1, 0.022, "(b)", fontsize = 14, fontweight = "bold")
    ax[1,0].text(4.1, 0.022, "(c)", fontsize = 14, fontweight = "bold")
    ax[1,1].text(4.1, 0.022, "(d)", fontsize = 14, fontweight = "bold")
    ax[0,0].text(6.7, 0.037, "OTF", fontsize = 14, fontweight = "bold")
    ax[0,1].text(6.7, 0.037, "OTF", fontsize = 14, fontweight = "bold")
    ax[1,0].text(6.7, 0.037, "ITF", fontsize = 14, fontweight = "bold")
    ax[1,1].text(6.7, 0.037, "ITF", fontsize = 14, fontweight = "bold")
   

    #new_tick_locations = np.array([8, 10, 12, 14, 16])

    #def tick_function(r_rsep, psin, ticks):
    #    f = scinter.interp1d(r_rsep, psin, fill_value = 'extrapolate')
    #    new_ticks = f(ticks)
    #    return ["%.3f" % z for z in new_ticks]

    #ax2.set_xlim(ax.get_xlim())
    #ax2.set_xticks(new_tick_locations)
    #ax2.set_xticklabels(tick_function(r_rsep, psin, new_tick_locations))
    #ax2.set_xlabel("Psin", fontsize = 'large', fontweight = "bold")
    #ax2.xaxis.set_minor_locator(MultipleLocator(5))
    plt.tight_layout()
    plt.show()

    return

def totcontent_plot():
    cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/CP MiMES Fwd Dens Scan.xlsx"
    xl = pd.ExcelFile(cpexcelpath)

    sums = []
    mcprs = []
    mcpls = []
    shots = [190428, 190450]
    for shot in shots:
        for sheet_name in xl.sheet_names:
            if str(shot) in sheet_name:    
                df = pd.read_excel(cpexcelpath, sheet_name = sheet_name)
                
                psins = df['Psin Shot: '+str(shot)].to_numpy()
                ws = df['W areal density (1e15cm-2)'].to_numpy()
                rrsep = df['R-Rsep OMP (cm) Shot:'+str(shot)].to_numpy()
                #psin_crit_lb = 1.2 #r-rsep (8.44 cm - 15.4 cm)
                wstrunc = []
                psis_trunc = []
                rrsep_trunc = []
                for i in range(0, len(ws)-1):
                    if psins[i] > 1.20 and psins[i] < 1.35:
                        wstrunc.append(ws[i])
                        psis_trunc.append(psins[i])
                        rrsep_trunc.append(rrsep[i])
                
                tot_cont_rs = scinteg.simps(wstrunc, rrsep_trunc)
                print(tot_cont_rs)
                plt.scatter(rrsep_trunc, wstrunc)
                #sums.append(tot_cont_rs)
                if 'L' in sheet_name:
                    mcpls.append(tot_cont_rs)
                elif 'R' in sheet_name:
                    mcprs.append(tot_cont_rs)

    plt.show()
    """
    totals = []
    for i in range(0, len(sums)-1, 2): #seperate probe faces
        totals.append(sums[i] + sums[i+1])
    """
    print(mcpls)
    print(mcprs)
    #ticks = ["MCP01W", "MCP02W", "MCP03W", "MCP10W"]
    densticks = ["3.15e19 Shot 190423", "3.51e19 Shot 190425", "3.8e19 Shot 190426", "4.40e19 Shot 190492"]
    
    #fig, ax = plt.subplots(nrows = 1, ncols = 2, sharex=True, sharey = True, figsize=(12, 6))
    fig, ax = plt.subplots()
    
    x = [1,2,3,4]
   
    ax.scatter(x, mcpls, label = "ITF")
    ax.scatter(x, mcprs, label = "OTF")
    ax.set_xlim(left = 0, right = 5)

    #ax2.set_xlim(ax.get_xlim())
    #ax2.set_xticks(x)
    #ax2.set_xticklabels(ticks)
    #ax3.set_xlim(ax.get_xlim())
    #ax3.set_xticks(x)
    #ax3.set_xticklabels(ticks)

    ax.set_ylabel("Collector Probe Total Content (1e15)", fontsize = 'large', fontweight = "bold")
    #ax[0,0].set_title("ITF")
    #ax[0,1].set_title("OTF")
    ax.set_xticks(x)
    ax.set_xticklabels(densticks)
    ax.set_xlabel("Density (m-3)", fontsize = 'large', fontweight = "bold")
    ax.set_xlabel("Density (m-3)", fontsize = 'large', fontweight = "bold")
    
    #fig.suptitle("Unfavorable Bt Density Scan MCP Total W Content", fontweight = "bold")
    ax.set_title("Unfavorable Bt Density Scan MCP Total W Content", fontweight = "bold")
    ax.legend(loc = "upper left")
    plt.tight_layout()
    plt.show()

    return

def mimes_dimes():
    
    fig, ax = plt.subplots(nrows = 1, ncols = 2, sharex= True, sharey = True, figsize=(12, 6))
    
    xdimes = [5.5, 6, 6.5, 7, 7.5]
    dimes_ticks = xdimes
    xmimes = [9, 12, 15, 18, 21]
    mimes_ticks = xmimes
    shots = [190428, 190450]
    cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/h-mode sinclair messer cp data.xlsx"
    for shot in shots:
        xl = pd.ExcelFile(cpexcelpath)
        #i = 0
        for sheet_name in xl.sheet_names:
            if str(shot) in sheet_name:
                df = pd.read_excel(xl, sheet_name = sheet_name)
                #shot = shots[i]
                r_rsep = df['R-Rsep OMP (cm) Shot:'+str(shot)].to_numpy()
                #psin = df["Psin"].to_numpy()
                #if shot == 190454:
                #    w_arealdens = df["W areal density (1e15cm-2)"].to_numpy()/3
                #else:
                #    w_arealdens = df["W areal density (1e15cm-2)"].to_numpy()/2
                w_arealdens = df["W areal density (1e15cm-2)"].to_numpy()
                if "fwd" in sheet_name:    
                    label = "Shot:" + " " +str(shot)+ ", Unfavorable Bt, Line average density = 6.48e19 (m-3)"
                    if 'M' in sheet_name and 'R' in sheet_name:
                        ax[0].scatter(r_rsep, w_arealdens, label = label, s = 10, alpha = 0.3)
                        #j = 0
                        k = 0
                    elif 'M' in sheet_name and 'L' in sheet_name:
                        ax[1].scatter(r_rsep, w_arealdens, label = label, s = 10, alpha = 0.3)
                        #j = 0
                        k = 1
                    """
                    elif 'D' in sheet_name and 'R' in sheet_name:
                        ax[0, 1].scatter(r_rsep, w_arealdens, label = label, s = 10, alpha = 0.3)
                        j = 0
                        k = 1
                    elif 'D' in sheet_name and 'L' in sheet_name:
                        ax[1, 1].scatter(r_rsep, w_arealdens, label = label, s = 10, alpha = 0.3)
                        j = 1
                        k = 1
                    """
                else:
                    label = "Shot:" + " " +str(shot)+ ", Favorable Bt, Line average density = 5.8e19 (m-3)"
                    if 'M' in sheet_name and 'L' in sheet_name:
                        ax[0].scatter(r_rsep, w_arealdens, label = label, s = 10, alpha = 0.3)
                        #j = 0
                        k = 0
                    elif 'M' in sheet_name and 'R' in sheet_name:
                        ax[1].scatter(r_rsep, w_arealdens, label = label, s = 10, alpha = 0.3)
                        #j = 0
                        k = 1
                    """
                    elif 'D' in sheet_name and 'L' in sheet_name:
                        ax[0, 1].scatter(r_rsep, w_arealdens, label = label, s = 10, alpha = 0.3)
                        j = 0
                        k = 1
                    elif 'D' in sheet_name and 'R' in sheet_name:
                        ax[1, 1].scatter(r_rsep, w_arealdens, label = label, s = 10, alpha = 0.3)
                        j = 1
                        k = 1
                    """

                """Need to fix this. R for favorble is ITF, OTF for unfavorable"""
                #if k == 0:
                ax[k].set_xlim(left = 7.0, right = 22.0)
                ax[k].set_xticks(xmimes)
                ax[k].set_xticklabels(mimes_ticks)
                #else:
                    #ax[j,k].set_xlim(left = 4.5, right = 8.0)
                    #ax[j,k].set_xticks(xdimes)
                    #ax[j,k].set_xticklabels(dimes_ticks)
                ax[k].set_xlabel("R-Rsep OMP [cm]", fontsize = 12, fontweight = "bold")
                ax[k].set_ylabel("W Areal Density [1e15 cm-2]", fontsize = 12, fontweight = "bold")
                ax[k].set_ylim(bottom = 0.0, top = 0.022)
                
                average = []
                newrs = []
                window = 5
                if 'M' in sheet_name:
                    for ind in range(window-1, len(w_arealdens) - window, 5):
                        average.append(np.mean(w_arealdens[ind-window:ind+window]))
                        newrs.append(r_rsep[ind])
                    ax[k].plot(newrs, average, label = "Shot "+str(shot)+ " averaged profile", linewidth = 4.5)
                    ax[k].legend(loc = 'upper right', fontsize = 10)

    ax[0].set_title("MCP H-mode OTF CP Profiles", fontsize = 13, fontweight = "bold")
    ax[1].set_title("MCP H-mode ITF CP Profiles", fontsize = 13, fontweight = "bold")
    #ax[0, 1].set_title("DCP H-mode CP Profiles", fontsize = 13, fontweight = "bold")
    #ax[0,0].text(7.2, 0.063, "(a)", fontsize = 15, fontweight = "bold")
    #x[0,1].text(4.6, 0.063, "(c)", fontsize = 15, fontweight = "bold")
    #ax[1,0].text(7.2, 0.063, "(b)", fontsize = 15, fontweight = "bold")
    #ax[1,1].text(4.6, 0.063, "(d)", fontsize = 15, fontweight = "bold")
    #ax[0,0].text(20.8, 0.04, "OTF", fontsize = 14, fontweight = "bold")
    #ax[0,1].text(7.7, 0.04, "OTF", fontsize = 14, fontweight = "bold")
    #ax[1,0].text(20.8, 0.04, "ITF", fontsize = 14, fontweight = "bold")
    #ax[1,1].text(7.7, 0.04, "ITF", fontsize = 14, fontweight = "bold")
    ax[0].text(21, 0.0175, "(a)", fontsize = 15, fontweight = "bold")
    ax[1].text(21, 0.0175, "(b)", fontsize = 15, fontweight = "bold")
   
    plt.tight_layout()
    plt.show()
    return

def profile_plot():
    
    cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_Sinclair_mastersheets/CP MiMES Fwd Dens Scan.xlsx"
    xl = pd.ExcelFile(cpexcelpath)
    mpl.rc('font', size = 24)
    fig, ax = plt.subplots(figsize=(10, 8))
    x = [7, 10, 13, 16, 19]
    #x = [5.5, 6, 6.5, 7, 7.5]
    ticks = x
    
    shots = [190422, 190422]
    labels = ["ITF", "OTF"]
    #labels = [' Mode Reference', 'SASOK Ne Puffing', 'GASC Ne Puffing', 'RF10 Ne Puffing', 'CPMID NE Puffing', 'IPD B Injection'] #Parsons

    i = 0
    for sheet_name in xl.sheet_names:
        if "190422" in sheet_name: #or "11x" in sheet_name:
            #if "R" in sheet_name:
            df = pd.read_excel(xl, sheet_name = sheet_name)
            shot = shots[i]
            if "R" in sheet_name:
                r_rsep = df['R-Rsep OMP (cm) Shot:'+str(shot)].to_numpy()-1.0
            else:
                r_rsep = df['R-Rsep OMP (cm) Shot:'+str(shot)].to_numpy()
            w_arealdens = df['W areal density (1e15cm-2)'].to_numpy()
            ax.scatter(r_rsep, w_arealdens, s = 10, alpha = 0.3)
            ax.set_xlabel("R-Rsep OMP [cm]", fontsize = 20, fontweight = "bold")
            ax.set_ylabel("W Areal Density [1e15 cm-2]", fontsize = 20, fontweight = "bold")
            ax.set_ylim(bottom = 0.0, top = 0.03)
            ax.set_xlim(left = 6, right = 20)
            #ax.set_xlim(left = 5.0, right = 8.0)
            ax.set_xticks(x)
            ax.set_xticklabels(ticks)
            average = []
            newrs = []
            window = 5
            for ind in range(window-1, len(w_arealdens) - window, 5):
                average.append(np.mean(w_arealdens[ind-window:ind+window]))
                newrs.append(r_rsep[ind])
            ax.plot(newrs, average, label = labels[i], linewidth = 4.5)
            ax.legend(loc = 'upper right', fontsize = 15)
            i = i + 1
    #norm = w_arealdens/max(w_arealdens)
    #ax1.plot(newrs, average, label = sheet_name+' Averaged', linewidth = 4)

    ax.set_title("MCP 01W Deposition", fontsize = 20, fontweight = "bold")
   

    #new_tick_locations = np.array([8, 10, 12, 14, 16])

    #def tick_function(r_rsep, psin, ticks):
    #    f = scinter.interp1d(r_rsep, psin, fill_value = 'extrapolate')
    #    new_ticks = f(ticks)
    #    return ["%.3f" % z for z in new_ticks]

    #ax2.set_xlim(ax.get_xlim())
    #ax2.set_xticks(new_tick_locations)
    #ax2.set_xticklabels(tick_function(r_rsep, psin, new_tick_locations))
    #ax2.set_xlabel("Psin", fontsize = 'large', fontweight = "bold")
    #ax2.xaxis.set_minor_locator(MultipleLocator(5))
    plt.tight_layout()
    plt.show()

    return