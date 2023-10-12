# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:07:35 2023

@author: Jake Maeker
"""

import oedge_plots
import matplotlib.pyplot as plt
import sys
import scipy
import scipy.io
import numpy as np
import netCDF4 as nc
import matplotlib as mpl
import numpy as np

plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "Times New Roman"

sys.path.append("C:/Users/jmateja/Documents/DIVIMP Related/")

#bg is the path to the background .nc file, b is for baffle .nc file, d is divertor .ncfile
ncpath_bg = "C:/Users/jmateja/Documents/DIVIMP Related/w-56854-ext-improved.nc"
ncpath_b = "C:/Users/jmateja/Documents/DIVIMP Related/w-56854-ext-B-100k-baffle.nc"
ncpath_d = "C:/Users/jmateja/Documents/DIVIMP Related/w-56854-ext-B-100k-lodivou.nc"

op = oedge_plots.OedgePlots(ncpath_bg)

bg_data = nc.Dataset(ncpath_bg)
b_data = nc.Dataset(ncpath_b)
d_data = nc.Dataset(ncpath_d)

idx_cell = bg_data['KVOLS'][:] > 0
index = bg_data['KORPG'][:][idx_cell] - 1

#hydrogen radiated power = HPOWLS (units are W/m3)
#total impurity radiated power = POWLS


h_prad = b_data['HPOWLS'][:,:]
b_prad = b_data['POWLS'][:,:,:] * b_data['ABSFAC'][:]
d_prad = d_data['POWLS'][:,:,:] * d_data['ABSFAC'][:]
h_dens = b_data['KNBS'][:,:]
b_dens = b_data['DDLIMS'][:,:,:] * b_data['ABSFAC'][:]
d_dens = d_data['DDLIMS'][:,:,:] * d_data['ABSFAC'][:]
i_temp = b_data['KTIBS'][:,:]
e_temp = b_data['KTEBS'][:,:]
nrs = b_data['NRS'][:]
nks = b_data['NKS'][:]
rvertp = b_data['RVERTP'][:]
zvertp = b_data['ZVERTP'][:]
korpg = b_data['KORPG'][:]
area = b_data['KAREAS'][:]
rs = b_data['RS'][:]
zs = b_data['ZS'][:]
#print(rs)

mesh = []
num_cells = 0

gridrs = []
gridzs = []
for ir in range(nrs):

    # Scan through the knots.
    for ik in range(nks[ir]):

        # Get the cell index of this knot on this ring.
        index = korpg[ir, ik] - 1

        # Only if the area of this cell is not zero append the corners.
        if area[ir, ik] != 0.0:
            vertices = list(zip(rvertp[index][0:4], zvertp[index][0:4]))
            mesh.append(vertices)
            num_cells = num_cells + 1

            # Print out a warning is the cell center is not within the vertices.
            cell = mpl.path.Path(list(vertices))
            r = rs[ir, ik]
            z = zs[ir, ik]
            gridrs.append(r)
            gridzs.append(z)
            if not cell.contains_point([r, z]):
                print("Error: Cell center not within vertices.")
                print("  (ir, ik)    = ({}, {})".format(ir, ik))
                print("  Vertices    = {}".format(vertices))
                print("  Cell center = ({}, {})".format(r, z))

# Save the results in the class.
print(max(gridzs))
num_cells = num_cells
mesh = mesh

mask = np.copy(rs)
#z_mask = np.copy(zs)
for i in range(0, 190):
    for j in range(0, 260):
        if rs[i][j] in gridrs and zs[i][j] in gridzs:
            mask[i][j] = 1
        else:
            mask[i][j] = 0



sum_prad = np.sum(h_prad, axis=0) + np.sum(b_prad, axis=0) + np.sum(d_prad, axis=0)
sum_pradW = np.sum(b_prad, axis=0) + np.sum(d_prad, axis=0) 
sum_pradWd = np.sum(d_prad, axis=0)
sum_pradWb = np.sum(b_prad, axis=0)
sum_dens = h_dens + np.sum(b_dens, axis=0) + np.sum(d_dens, axis=0)
sum_densW = np.sum(b_dens, axis=0) + np.sum(d_dens, axis=0)
sum_densWd = np.sum(d_dens, axis=0)
sum_densWb = np.sum(b_dens, axis=0)


# print(np.shape(multispecies_prad))
total_prad = sum_prad[idx_cell]
total_pradW = sum_pradW[idx_cell]
total_pradWb = sum_pradWb[idx_cell]
total_pradWd = sum_pradWd[idx_cell]
total_dens = sum_dens[idx_cell]
total_densW = sum_densW[idx_cell]
total_densWb = sum_densWb[idx_cell]
total_densWd = sum_densWd[idx_cell]
total_h_dens =h_dens[idx_cell]
total_i_temp = i_temp[idx_cell]
total_e_temp = e_temp[idx_cell]


op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Total radiated power (W m-3)", cmap='nipy_spectral', own_data=total_prad)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Impurity radiated power (W m-3)", cmap='nipy_spectral', own_data=total_pradW)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Impurity radiated power baffle (W m-3)", cmap='nipy_spectral', own_data=total_pradWb)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Impurity radiated power divertor (W m-3)", cmap='nipy_spectral', own_data=total_pradWd)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Total density (m-3)", cmap='nipy_spectral', own_data=total_dens)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Impurity density (m-3)", cmap='nipy_spectral', own_data=total_densW)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Impurity density baffle (m-3)", cmap='nipy_spectral', own_data=total_densWb)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Impurity density divertor (m-3)", cmap='nipy_spectral', own_data=total_densWd)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Plasma density (m-3)", cmap='nipy_spectral', own_data=total_h_dens)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Electron Temperature (eV)", cmap='nipy_spectral', own_data=total_e_temp)
op.plot_contour_polygon(dataname='none',scaling=op.absfac, normtype='log', cbar_label="Ion Temperature (eV)", cmap='nipy_spectral', own_data=total_i_temp)

scipy.io.savemat('PRAD_56854.mat',{'h_prad': h_prad, 'b_prad': b_prad, 'd_prad': d_prad,'h_dens': h_dens, 'b_dens': b_dens, 'd_dens': d_dens, 'i_temp': i_temp, 'e_temp': e_temp, 'R': rs,'Z': zs, 'mask': mask},)



