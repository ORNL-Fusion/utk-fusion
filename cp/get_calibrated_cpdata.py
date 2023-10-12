"""
Program: get_calibrated_cpdata.py
Author: Jeremy Mateja
Email:  jmateja@vols.utk.edu
Date: 8/9/2022

This program pulls W areal densities and distance along probe from a master Excel file. It then
outputs corrected machine coordinate data based on probe angle offset and separatrix location. Returns
probe locations in terms of psin, R probe mapped to midplane, R-Rmidout
"""

import numpy as np
import pandas as pd
import meas_locations
import Calculate_rho_gt1 as crgt
import scipy.interpolate as scinter
#import sys
"""
The Excel sheet should contain the R value of the tip. Hopefully this keeps things more consistent and
less confusing. The tip position will be on the CP consolidated spreadsheet and is a value given by
Dmitry. Note this is only something you have to do when dealing with MiMES. DiMES will always have the
same value.
"""
#Specify your Master Sheet path
cpexcelpath = "/mnt/c/Users/jmateja/Documents/SASVWcp/Cacheris_Abrams_mastersheets/raw_mcps_21_23.xlsx"

# In the future it may be better to just specify the sheet beforehand instead of redoing all of them
time = 3000 # This should maybe never change. Check reviewplus. Can pick anytime on the flat top region

# This is a set up for reading through all of the sheets. Might have to change this when only looking at one sheet
xl = pd.ExcelFile(cpexcelpath)


#Note: The following function should only be used for MiMES probes. This function reads
#the master excel sheet and calls the meas_locations function to correct for the 13 degree
#offset in the MiMES probe
def offset_angle_correction(df, sheet):
    
    r_probe = df['R tip (cm)'].loc[df.index[0]]
    locations = df['Distance along probe (cm)'].to_numpy()

    # Hopefully the naming convention doesn't change
    # r_probe is the r location (cm) of the tip of the probe. This should be given by Dmitry/Shawn
    # locations is "Distance along probe (cm)"
    if 'L' in sheet:
        r_corr = meas_locations.calc_R_measAD(r_probe, locations) #AD is for the left side (Think this is only for USN)
    else:
        r_corr = meas_locations.calc_R_measAU(r_probe, locations) #AU is for the right side
    return r_corr

def get_psin(shot, tree, MiMES = True, time = time): #used to include Rin, Zin

    # crgt is an python file in the EFIT class from the ORNL GitHub. Need to ask Zeke for access I believe
    '''For Seth: crgt.standalone is the function you want to change'''
    if MiMES:
        f_psiN = crgt.standalone(shot, time, tree = tree)
    else:
        f_psiN = crgt.standalone(shot, time, tree = tree, MiMES = False)
    return f_psiN

def get_rmapped(probepsi, zm, tree, shot):

    # This gets r-rsep OMP (cm). First get the psin values at the midplane, then interpolate with probe psi values
    rs = np.linspace(2, 2.5, 100)
    zs = [zm]*len(rs)
    fpsin = get_psin(shot, tree) # psin function at midplane
    
    psin = fpsin(rs, zs)
    #Stuff below is for the one DCP that was in for three shots. Shouldn't have to uncomment this too much

    fromp = scinter.interp1d(psin, rs, fill_value = "extrapolate")
    rmapped = fromp(probepsi) # Gets rprobe mapped to midplane
    rmid = fromp(1.0)
    rfinal = rmapped - rmid # r-rsep OMP 
    return rfinal

def rsep_probe(probesci, r_corr, fpsin = None, MiMES = True):

    #This function gets r-rsep at the probe location
    if MiMES:
        f = scinter.interp1d(probesci, r_corr, fill_value="extrapolate")
        probe_rsep = f(1.0) # Find rsep at probe location
        probe_seprs = abs(r_corr) - abs(probe_rsep)
    else: #For DCP
        zs = np.linspace(-1.4, -0.75, 250)
        rs = [1.4859]*len(zs)
        new_psis = fpsin(rs, zs)
        #new_psis2 = fpsin2(rs, zs)
        psi_trunc = new_psis < 1.001
        new_zs = zs[psi_trunc]
        z_sep = new_zs[0]*100
        probe_seprs = abs(r_corr) - abs(z_sep)

    return probe_seprs

def get_werrors(df):
    #Uses error prop to get W count errors
    W180 = df['W180'].to_numpy()
    W182 = df['W182'].to_numpy()
    W183 = df['W183'].to_numpy()
    W184 = df['W184'].to_numpy()
    W186 = df['W186'].to_numpy()

    W180_error = np.sqrt(W180)
    W182_error = np.sqrt(W182)
    W183_error = np.sqrt(W183)
    W184_error = np.sqrt(W184)
    W186_error = np.sqrt(W186)

    WTOT_errors = np.sqrt(np.power(W180_error, 2)+np.power(W182_error, 2)+np.power(W183_error, 2)+
    np.power(W184_error, 2)+np.power(W186_error, 2))

    return WTOT_errors

def calibration(df, calconstant, slope_error, errors):

    #Gets W areal densities. Be very careful in making sure you are using the right calibration
    #Will have to make changes as new probes and calibrations continue to be analyzed. There's 
    #definitely a better way of doing this but it was nice to keep track of everything this way.
    WTOT = df['WTOT'].to_numpy()
   
    areal_densities = calconstant*WTOT
    areal_dens_error = areal_densities*np.sqrt(np.power(errors/WTOT, 2)+np.power(slope_error/calconstant, 2))
    
    return areal_densities, areal_dens_error


#This function should call all of the functions in the class to make a comprehensive spreadsheet with the new coordinate data
def make_spreadsheet(xl = xl):

    # Cycle through all sheet names and apply the right shot and rmidout based on the probe number
    # Admittedly this current method seems a bit ridiculous, but will be cleaned up for new probes
    for sheet_name in xl.sheet_names:
        df = pd.read_excel(cpexcelpath, sheet_name = sheet_name)
        
        #Listing all the shots like this is probably ridiculous but it was useful for keeping track of stuff

        calconstant = df['calconstant'].loc[df.index[0]]
        tree = str(df['Tree'].loc[df.index[0]])
        slope_error = df['slope_error'].loc[df.index[0]]
        shot = int(df['shot'].loc[df.index[0]])
        zm = df['zm'].loc[df.index[0]]
        tree2 = str(df['Tree2'].loc[df.index[0]])
        shot2 = int(df['shot2'].loc[df.index[0]])
        zm2 = df['zm2'].loc[df.index[0]]
        
        #Check if MiMES or DiMES probe in excel sheet
        if 'M' in sheet_name:
            r_corr = offset_angle_correction(df, sheet_name)
            z_corr = [-0.188]*len(r_corr)
            fpsin = get_psin(shot, tree)
            fpsin2 = get_psin(shot2, tree2)
            psis = fpsin(r_corr/100., z_corr)
            psis2 = fpsin2(r_corr/100., z_corr)
            #psis = (psis1 + psis2)/2
            probe_seprs_shot = rsep_probe(psis, r_corr)
            probe_seprs_shot2 = rsep_probe(psis2, r_corr)
            df['R (cm)'] = r_corr
            df['R-Rsep (cm) Shot: '+ str(shot)] = probe_seprs_shot
            df['R-Rsep (cm) Shot: '+ str(shot2)] = probe_seprs_shot2
        else: # For DCP
            z_probe = -120.5
            locations = df['Distance along probe (cm)'].to_numpy()
            z_corr = z_probe - locations
            r_corr = [1.4859]*len(z_corr)
            
            fpsin = get_psin(shot, tree, MiMES = False)
            fpsin2 = get_psin(shot2, tree2, MiMES = False)
            psis = fpsin(r_corr, z_corr/100.)
            psis2 = fpsin2(r_corr, z_corr/100.)
            probe_seprs_shot = rsep_probe(psis, z_corr, fpsin, MiMES=False)
            probe_seprs_shot2 = rsep_probe(psis2, z_corr, fpsin2, MiMES = False)
            df['Z (cm)'] = z_corr
            df['Z-Zsep (cm) Shot: '+str(shot)] = probe_seprs_shot
            df['Z-Zsep (cm) Shot: '+str(shot2)] = probe_seprs_shot2


        rfinal_shot = get_rmapped(psis, zm, tree, shot)
        rfinal_shot2 = get_rmapped(psis2, zm2, tree2, shot2)  
        w_errors = get_werrors(df)
        areal_densities, areal_dens_error = calibration(df, calconstant, slope_error, w_errors)
        df.pop('WTOT') # Don't need total tungsten numbers in final sheet
        df.pop('W180')
        df.pop('W182')
        df.pop('W183')
        df.pop('W184')
        df.pop('W186')
        df.pop('calconstant')
        df.pop('slope_error')
        df.pop('shot')
        df.pop('zm')
        df.pop('Tree')
        df.pop('shot2')
        df.pop('zm2')
        df.pop('Tree2')
        #df['Psin Shot: '+str(shot3)] = psis3
        df['Psin Shot: '+str(shot)] = psis
        df['Psin Shot: '+str(shot2)] = psis2
        df['R-Rsep OMP (cm) Shot:'+str(shot)] = rfinal_shot*100
        df['R-Rsep OMP (cm) Shot:'+str(shot2)] = rfinal_shot2*100
        #df['R-Rsep OMP Uncertainty (cm)'] = r_errors*100
        df['W areal density (1e15cm-2)'] = areal_densities
        df['W areal density error (1e15cm-2)'] = areal_dens_error

        #Makes output files. Make sure to change the path
        if 'M' in sheet_name:
            df.to_excel("/mnt/c/Users/jmateja/Documents/SASVWcp/SASVW Updated CP Spreadsheets/MCP/xlsx/"+sheet_name+".xlsx")
            df.to_csv("/mnt/c/Users/jmateja/Documents/SASVWcp/SASVW Updated CP Spreadsheets/MCP/csv/"+sheet_name+".csv")
        else:
            df.to_excel("/mnt/c/Users/jmateja/Documents/SASVWcp/SASVW Updated CP Spreadsheets/DCP/xlsx/"+sheet_name+".xlsx")
            df.to_csv("/mnt/c/Users/jmateja/Documents/SASVWcp/SASVW Updated CP Spreadsheets/DCP/csv/"+sheet_name+"new.csv")
    return