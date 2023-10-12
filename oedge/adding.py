import oedge_plots

ncpath1 = "/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/d3d-187111-inj-24.nc"
ncpath2 = "/Users/jmateja/Documents/GitHub/utk-fusion/DIVIMP Related/187111/d3d-187111-inj-25.nc"

op1 = oedge_plots.OedgePlots(ncpath1)
op2 = oedge_plots.OedgePlots(ncpath2)
ddlim1 = op1.read_data_2d("DDLIMS", scaling=op1.absfac, charge="all")
ddlim2 = op2.read_data_2d("DDLIMS", scaling=op2.absfac, charge="all")

comb_ddlim = ddlim1 + ddlim2

op1.plot_contour_polygon(dataname="DDLIMS", own_data=comb_ddlim, normtype="log", cbar_label="C13 Density (m-3)", cmap="nipy_spectral")