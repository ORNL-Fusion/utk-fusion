# utk-fusion
This repository is a collection of Python scripts related to accessing and analyzing collector probe data, OEDGE/3DLIM/SOLPS plotting GUIs and routines, common DIII-D data analysis such as loading and filtering Thomson scattering and Langmuir probe data. See the README in each folder for descriptions of how to use the scripts.

## cp
Contains files related to accessing collector probe data and collector probe relevant tasks.

## lim
Contains the 3DLIM plotting GUI to help vizualize 3DLIM output.

## oedge
Contains the OEDGE (DIVIMP) plotting GUI to help vizualize OEDGE output.

## solps
Contains the SOLPS-ITER plotting scripts. This is still a work in progress, but the basic bones are there.

## tools
Contains useful tools related to DIII-D data access, such as pulling and filtering Langmuir probe or Thomson scattering data.

## misc
Contains collection of miscellaneous scripts related to structuring LAMS experimental data.


## Instructions on linking localhost to atlas
When not on the DIII-D network, you will need to link localhost (i.e., your computer) to atlas so that you can use MDSplus to pull data from the server. This can easily be done if you are running Linux or Ubuntu on Windows 10 (WSL). Simply open up another terminal and enter:
```
ssh -Y -p 2039 -L 8000:atlas.gat.com:8000 username@cybele.gat.com
```
It will prompt you for your cybele password. Just leave this terminal open and you should be able to access the MDSplus data stored on atlas. This has not been tested with Putty, but try and it update this README if it works. 
