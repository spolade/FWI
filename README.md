# DT CLIMATE - Wildfire Use Case (FWI Index)

![latest_release: v2.2.5](https://github.com/spolade/FWI/releases/tag/v2.2.5)

This repository contains the scripts related to the Wildfire FWI use case of the Climate Adaptation Digital Twin (Climate DT). All the work is being developed in the frame of the [Destination Earth initiative](https://destination-earth.eu/) from the European Commission, where [ECMWF](https://destine.ecmwf.int/) is one of the Entrusted Entities.

LICENSE NOTE: the European Union, represented by the European Commission is the direct and sole owner of the intellectual property rights of these Results. 

## Description

This repository calculates the Canadian Fire Weather Index (FWI), a widely used index for evaluating and assessing fire weather conditions. The code is adapted from the original algorithm provided by Van Wagner and Pickett in 1985: Van Wagner, C.E.; Pickett, T.L. (1985). Equations and FORTRAN program for the Canadian Forest Fire Weather Index System. Canadian Forestry Service, Petawawa National Forestry Institute, Chalk River, Ontario. Forestry Technical Report 33, 18 p.
Several improvements have been made to adapt the code for global use by reducing regional dependencies.
The FWI application consists of two main files:

functions_cal_fwi.py: the submodule containing calculation functions
run_wildfires_fwi.py: the main script that calls the submodule and performs the calculations

Instructions on the required packages and how to use the code are provided below.

## Implemented indicators



- **FWI Indexs**:
    - Input:
        - `tas:   xarray.DataArray ; (time,lat,lon)` -> temperature.
        - `pr:    xarray.DataArray ; (time,lat,lon)` -> precipitation.
        - `uwind: xarray.DataArray ; (time,lat,lon)` -> u-component of wind.
        - `vwind: xarray.DataArray ; (time,lat,lon)` -> v-component of wind.
        - `d2m:   xarray.DataArray ; (time,lat,lon)` -> dew point temperature.

    - Output:
        - `fwi: xarray.DataArray ; (time,lat,lon)` -> FWI index



## Version
Current version can be found at the latest publised tags in the git information.

## How to run

Each submodule includes a description of its aim, inputs, outputs and corresponding referenences. The following is an example of how to run the `FWI applications` in a Python environment:

```
import xarray as xr
from fwi.functions_cal_fwi import FFMC, DMC, DC, ISI, BUI, FWI



## Post-Processing

The different annual statistics are calculated in the same run script using the daily FWI values.

   1. Annual average
   2. Annual maximum
   3. Annual minimum
   4. Annual Percentiles FWI values for 50, 95, and 99 %
   5. Daily Severity Rating (DSR)
   6. Threshold exceedance days

            VLow -- wi < 5.2
            Low -- fwi >= 5.2 & fwi < 11.2
            Mod -- fwi >= 11.2 & fwi < 21.3
            High -- fwi >= 21.3 & fwi < 38.0
            VHigh -- fwi >= 38.0 & fwi < 50
            Ext -- fwi >= 50





## Support

For any feedback, comments and/or issues you can contact me through Gitlab 

-------
To install the necessary dependencies for the package:
```
pip install git+https://github.com/spolade/FWI.git@main
```

To copy the repository to your local directory:
```
git clone https://github.com/spolade/FWI.git
```

To install the package locally from the root directory (where the `setup.py` file is located):
```
pip install .
```

