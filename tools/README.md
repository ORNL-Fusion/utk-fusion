# tools
This is a collection of tools mostly related to DIII-D data access.

## gadata.py
This is a DIII-D script found on the internal website, i.e. Shawn did not write this. Example usage:
```
$ python3
>>> from gadata import gadata
>>> import MDSplus as mds
>>> conn = mds.Connection('localhost')   # Needed if linking iris to localhost. See utk-fusion README for details.
>>> conn = mds.Connection('atlas.gat.com')  # If working on NoMachine or on the DIII-D network (needs testing).
>>> ga = gadata("DENSV2", 167196, connection=conn) 
>>> time = ga.xdata
>>> dens = ga.zdata
```

## get_lp.py
This allows one to pull target Langmuir probe data into a Python session and optionally save it as a .csv file. Example usage:
```
$ python3
>>> import get_lp

# Get LPs plotted against R-Rsep, bin each LP data into 5 bins equally spaced in indicated 
# time range and report the median values in each bin.
>>> lps = get_lp.plot_lps(167195, 2000, 5000, xtype="rminrsep", bins=5, filter="median")
>>> lps.keys()

# Plot against psin, report average values of each bin.
>>> lps = get_lp.plot_lps(167195, 2000, 5000, xtype="psin", bins=5, filter="average", csv_path="/path/to/save/lp_167196.csv")

```
<img src="lp_ex.png" width="400" height="300">
