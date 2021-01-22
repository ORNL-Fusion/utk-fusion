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
<p align="center">
  <img src="lp_ex.png" width="450" height="200">
</p>

Note the inner and outer target LPs are plotted on the same graph here. 

## elm_frequency.py
A relatively simple script to estimate the ELM frequency. There are more rigourous ways to do this, but this is a simple starting point I suppose. Example usage:

```
$ python3
>>> import elm_frequency as elm
>>> time, fs = elm.get_fs_data(167247, "FS04")  # Load filterscope data.

# Pass in the filterscope data and see if we correctly identify the peaks. distance and height are
# the same as distance, height in the scipy function find_peaks, so best read the documentation for more info.
>>> elm.calc_freq(time, fs, "FS04", time_window=[3000, 4500], distance=200, height=2e16)
```
It will probably require some trial and error playing with distance and height. The red stars indicate a detected ELM. You want to make sure you aren't double counting ELMs, or counting tiny spikes in the data as ELMs (unless you think they might be!).

<p align="center">
  <img src="elm_ex.png" width="300" height="225">
</p>

## ThomsonClass.py
To do. Unsure if I want to document this since generally it may be better to use OMFIT when it comes to getting clean filtered data, but if you want just the "raw" ne, Te data this will work. Also the elementary 2D plot capability, but still maybe best to hunt down Adam Mclean's IDL scripts on iris for that.
