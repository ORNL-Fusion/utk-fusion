"""
Program: get_cp_coordinates.py
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
#import matplotlib.pyplot as plt
"""
The Excel sheet should contain the R value of the tip. Hopefully this keeps things more consistent and
less confusing. The tip position will be on the CP consolidated spreadsheet and is a value given by
Dmitry. Note this is only something you have to do when dealing with MiMES. DiMES will always have the
same value.
"""
#Specify your Master Sheet path
cpexcelpath = "/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/Messer_Sinclair_mastersheets/SASVW CP Master Sheet 10W.xlsx"

# In the future it may be better to just specify the sheet beforehand instead of redoing all of them
time = 3500 # This should maybe never change. Check reviewplus. Can pick anytime on the flat top region

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
    """
    if shot3:
        fpsin3 = get_psin(shot3, tree3)
        psin3 = fpsin3(rs, zs)
        psin_av = (psin1 + psin2 + psin3) / 3
        fromp3 = scinter.interp1d(psin3, rs, fill_value = "extrapolate")
        rmid3 = fromp3(1.0)
    else:
        psin_av = (psin1 + psin2) / 2
    """
    fromp = scinter.interp1d(psin, rs, fill_value = "extrapolate")
    """Ignore these commented out lines for now. From a different point when I was
    writing. Could probably delete them"""
    #fromp2 = scinter.interp1d(psin2, rs, fill_value = "extrapolate")
    
    #fromp_av = scinter.interp1d(psin_av, rs, fill_value = "extrapolate")
    
    #rmid1 = fromp1(1.0)
    #rmid2 = fromp1(1.0)
    #rmid_av = fromp_av(1.0)

    #rmapped1 = fromp1(probepsi)
    #rmapped2 = fromp2(probepsi)

    rmapped = fromp(probepsi) # Gets rprobe mapped to midplane
    rmid = fromp(1.0)
    #rmapped1 = fromp_av(probepsi1)
    #rmapped2 = fromp_av(probepsi2)
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
"""
def get_space_errors(psi_av, psis1, psis2, rs, psis3 = None):

    #Need to find errors in psi since I am averaging them. Best way currently would be to take difference between
    # average psi and individual psis. square then square root (std deviation basically)
    #fromp1 = scinter.interp1d(psis1, rs, fill_value = "extrapolate")
    #intercept1 = fromp1(0)
    #print(intercept1)
    #slope1 = (fromp1(1)-intercept1)
    #print(slope1)
    #fromp2 = scinter.interp1d(psis2, rs, fill_value = "extrapolate")
    #intercept2 = fromp2(0)
    #print(intercept2)
    #slope2 = (fromp2(1)-intercept2)
    #intercept_av = function(0)
    #print(intercept_av)
    #slope_av = function(1)-intercept_av
    #print(slope_av)
    if psis3:
        psi_error = np.sqrt(1/3*np.power(psis1-psi_av, 2)+np.power(psis2-psi_av, 2)+np.power(psis3-psi_av, 2))

        fromp3 = scinter.interp1d(psis3, rs, fill_value = "extrapolate")
        intercept3 = fromp2(0)
        slope3 = (fromp1(1)-intercept3)
        slope_error = np.sqrt(1/3*(np.power(slope1 - slope_av, 2)+np.power(slope2 - slope_av, 2)+np.power(slope3 - slope_av, 2)))
        intercept_error = np.sqrt(1/3*(np.power(intercept1 - intercept_av, 2)+np.power(intercept2 - intercept_av, 2)+np.power(intercept3 - intercept_av, 2)))

    else:
        psi_error = np.sqrt(1/2*(np.power(psis1-psi_av, 2)+np.power(psis2-psi_av, 2)))
        print(psi_error)
    
        print(psi_error)
        slope_error = np.sqrt(1/2*(np.power(slope1 - slope_av, 2)+np.power(slope2 - slope_av, 2)))
        print(slope_error)
        intercept_error = np.sqrt(1/2*(np.power(intercept1 - intercept_av, 2)+np.power(intercept2 - intercept_av, 2)))
        print(intercept_error)
    
    #psis_slope_error = rs*np.sqrt(np.power(psi_error/psi_av,2)+np.power(slope_error/slope_av,2))
    #print(psis_slope_error)
    #sys.exit()
    #r_rsep_error = np.sqrt(np.power(rs*np.sqrt(np.power(psi_error/psi_av,2)+np.power(slope_error/slope_av,2)),2)+np.power(intercept_error,2))

    
    #Now find error in r-rsep slope using same premise? 

    return psi_error
"""
def calibration(sheet, df, errors):

    #Gets W areal densities. Be very careful in making sure you are using the right calibration
    #Will have to make changes as new probes and calibrations continue to be analyzed. There's 
    #definitely a better way of doing this but it was nice to keep track of everything this way.
    WTOT = df['WTOT'].to_numpy()
    """
    if "1" in sheet and "L" in sheet:
        calconstant = 6.75e-7 # For LW01 DiMES and MiMES only
        slope_error = 3.23e-8 # For LW01 DiMES and MiMES only
    elif "1" in sheet and "R" in sheet:
        calconstant = 9.25e-7 # For RW01 DiMES and MiMES only
        slope_error = 1.62e-8 # For RW01 DiMES and MiMES only
    elif "R" in sheet and "1" not in sheet:
        calconstant = 8.05e-7 # For RW02-RW05 DiMES and MiMES only
        slope_error = 2.23e-8 # For RW02-RW05 DiMES and MiMES only
    elif "L" in sheet and "1" not in sheet:
        calconstant = 6.01e-7 # For LW02-LW05 DiMES and MiMES only
        slope_error = 2.34e-8 # For LW02-LW05 DiMES and MiMES only
    
    if "R" in sheet:
        calconstant = 9.11e-7 # For RW06-RW09 DiMES and MiMES only
        slope_error = 1.95e-8
    elif "L" in sheet:
        calconstant = 6.94e-7 # For LW06-LW09 DiMES and MiMES only
        slope_error = 2.41e-8
    """
    if "L" in sheet:
        calconstant = 7.79E-7  # For LW10 DiMES and MiMES only
        slope_error = 2.51E-8
    else:
        calconstant = 8.13E-7 # For RW10 DiMES and MiMES only
        slope_error = 2.33E-8
    
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
        """
        if '1W' in sheet_name:
            shot1 = 190422
            shot2 = 190423
            zm1 = -0.083
            zm2 = -0.083
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
        elif '2W' in sheet_name:
            shot1 = 190424
            shot2 = 190425
            zm1 = -0.075
            zm2 = -0.080
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
        elif '3W' in sheet_name:
            shot1 = 190426
            shot2 = 190427
            zm1 = -0.081
            zm2 = -0.080
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
        elif '4W' in sheet_name:
            shot1 = 190428
            shot2 = 190429
            zm1 = -0.090
            zm2 = -0.086
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
        elif '5W' in sheet_name:
            shot1 = 190450
            shot2 = 190451
            zm1 = -0.084
            zm2 = -0.083
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
        
        if '6W' in sheet_name:
            shot1 = 190454
            shot2 = 190455
            shot3 = 190453
            zm1 = -0.075
            zm2 = -0.081
            zm3 = -0.081
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
            tree3 = 'EFIT02'
        elif '7W' in sheet_name:
            shot1 = 190456
            shot2 = 190457
            zm1 = -0.08
            zm2 = -0.082
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
        elif '8W' in sheet_name:
            shot1 = 190459
            shot2 = 190461
            zm1 = -0.077
            zm2 = -0.081
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
        elif '9W' in sheet_name:
            shot1 = 190481
            shot2 = 190482
            zm1 = -0.084
            zm2 = -0.077
            tree1 = 'EFIT02'
            tree2 = 'EFIT02'
        """
        if '10W' in sheet_name:
            shot1 = 190491
            shot2 = 190492
            zm1 = -0.076
            zm2 = -0.076
            tree1 = 'EFIT01'
            tree2 = 'EFIT02'
        
        
        #Check if MiMES or DiMES probe in excel sheet
        if 'M' in sheet_name:
            r_corr = offset_angle_correction(df, sheet_name)
            z_corr = [-0.188]*len(r_corr)
            fpsin1 = get_psin(shot1, tree1)
            fpsin2 = get_psin(shot2, tree2)
            psis1 = fpsin1(r_corr/100., z_corr)
            psis2 = fpsin2(r_corr/100., z_corr)
            #psis = (psis1 + psis2)/2
            probe_seprs_shot1 = rsep_probe(psis1, r_corr)
            probe_seprs_shot2 = rsep_probe(psis2, r_corr)
            df['R (cm)'] = r_corr
            df['R-Rsep (cm) Shot: '+ str(shot1)] = probe_seprs_shot1
            df['R-Rsep (cm) Shot: '+ str(shot2)] = probe_seprs_shot2
        else: # For DCP
            z_probe = -120.5
            locations = df['Distance along probe (cm)'].to_numpy()
            z_corr = z_probe - locations
            r_corr = [1.4859]*len(z_corr)
            
            fpsin1 = get_psin(shot1, tree1, MiMES = False)
            fpsin2 = get_psin(shot2, tree2, MiMES = False)
            #fpsin3 = get_psin(shot3, tree3, MiMES = False)
            psis1 = fpsin1(r_corr, z_corr/100.)
            psis2 = fpsin2(r_corr, z_corr/100.)
            #psis3 = fpsin3(r_corr, z_corr/100.)
            #if shot3:
            #    fpsin3 = get_psin(shot3, tree3, MiMES = False)
            #    psis3 = fpsin3(r_corr, z_corr/100.)
            #    probe_seprs_shot3 = rsep_probe(psis3, z_corr)
            probe_seprs_shot1 = rsep_probe(psis1, z_corr, fpsin1, MiMES=False)
            probe_seprs_shot2 = rsep_probe(psis2, z_corr, fpsin2, MiMES = False)
            #probe_seprs_shot3 = rsep_probe(psis3, z_corr, fpsin3, MiMES = False)
            df['Z (cm)'] = z_corr
            #df['Z-Zsep (cm) Shot: '+str(shot3)] = probe_seprs_shot3
            df['Z-Zsep (cm) Shot: '+str(shot1)] = probe_seprs_shot1
            df['Z-Zsep (cm) Shot: '+str(shot2)] = probe_seprs_shot2
            

            """
            if shot3:
                fpsin3 = get_psin(shot3, tree3, MiMES = False)
                psis = (fpsin1(r_corr, z_corr/100.) + fpsin2(r_corr, z_corr/100.) + fpsin3(r_corr, z_corr/100.)) / 3
                zm = (zm1 + zm2 + zm3) / 3
                probe_seprs = rsep_probe(psis, z_corr, MiMES = False, fpsin1 = fpsin1, fpsin2 = fpsin2, fpsin3 = fpsin3)
            else:
                psis = (fpsin1(r_corr, z_corr/100.) + fpsin2(r_corr, z_corr/100.)) / 2
                zm = (zm1 + zm2) / 2
                probe_seprs = rsep_probe(psis, z_corr, MiMES = False, fpsin1 = fpsin1, fpsin2 = fpsin2)
            """
        
        
        #if shot3:
            #rfinal, r_errors = get_rmapped(psis, zm, tree1, tree2, shot1, shot2, shot3 = shot3, tree3 = tree3)*100
        rfinal_shot1 = get_rmapped(psis1, zm1, tree1, shot1)
        rfinal_shot2 = get_rmapped(psis2, zm2, tree2, shot2)
        #rfinal_shot3 = get_rmapped(psis3, zm3, tree3, shot3)   
        w_errors = get_werrors(df)
        areal_densities, areal_dens_error = calibration(sheet_name, df, w_errors)
        df.pop('WTOT') # Don't need total tungsten numbers in final sheet
        df.pop('W180')
        df.pop('W182')
        df.pop('W183')
        df.pop('W184')
        df.pop('W186')
        #df['Psin Shot: '+str(shot3)] = psis3
        df['Psin Shot: '+str(shot1)] = psis1
        df['Psin Shot: '+str(shot2)] = psis2
        #df['R-Rsep OMP (cm) Shot:'+str(shot3)] = rfinal_shot3*100
        df['R-Rsep OMP (cm) Shot:'+str(shot1)] = rfinal_shot1*100
        df['R-Rsep OMP (cm) Shot:'+str(shot2)] = rfinal_shot2*100
        #df['R-Rsep OMP Uncertainty (cm)'] = r_errors*100
        df['W areal density (1e15cm-2)'] = areal_densities
        df['W areal density error (1e15cm-2)'] = areal_dens_error

        #Makes output files. Make sure to change the path
        if 'M' in sheet_name:
            df.to_excel("/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/SASVW Updated CP Spreadsheets/MCP/xlsx/"+sheet_name+"new.xlsx")
            df.to_csv("/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/SASVW Updated CP Spreadsheets/MCP/csv/"+sheet_name+"new.csv")
        else:
            df.to_excel("/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/SASVW Updated CP Spreadsheets/DCP/xlsx/"+sheet_name+"new.xlsx")
            df.to_csv("/mnt/c/Users/jmateja/Documents/GitHub/utk-fusion/cp/SASVW Updated CP Spreadsheets/DCP/csv/"+sheet_name+"new.csv")
    return