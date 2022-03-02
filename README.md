# MPM_Space_Sciences

## INTRODUCTION ##
These codes were developed during my Thesis work where I performed some simulations with the Planet Simulator Model (PlaSim). The Climate Model is open source and you can find the original version [here](https://github.com/HartmutBorth/PlaSim). However, in order to develop this project, I started from another version which can be found [here](https://github.com/jhardenberg/plasim). Even if the code in this repository was used for the PlaSim model I believe that with slight modifications it can be applied to other models that have their output data in a Netcdf format. In the following part I will explain how to use the code. This is a beta version. Consequently, bugs might occur. Feel free to open issues in that case.

**NOTE 1: The following instructions are thought for ubuntu.**

**NOTE 2: Some of the linked files can be too heavy to be opened on github. Therefore, you have to donwnload them and open them with jupyter.**

## PRE-REQUISITES ##
1. Download cdo `sudo apt-get update && sudo apt install cdo`;
2. Install netcdf `sudo apt-get install libnetcdf-dev libnetcdff-dev`;
3. Install [Anaconda distribution](https://noviello.it/come-installare-anaconda-su-ubuntu-20-04-lts/) with Jupyter Notebook (except for the file in CDO and application) with the following important packages:
   - netcdf:
    `conda install -c anaconda netcdf4`;
    - cdo:
    `pip install cdo`;
    - basemap:
    `conda install -c anaconda basemap`;
    - cartopy:
    `conda install -c conda-forge cartopy` or see [here](https://anaconda.org/conda-forge/cartopy);
    - xarray:
    `conda install -c anaconda xarray`
    
The basemap package is deprecated and in future a complete transition to cartopy is needed. However, basemap holds some features which are still not present in cartopy. As a consequence, I avoided an update regarding cartopy.
