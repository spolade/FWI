import os
import numpy as np
import subprocess

# Provide the input and output DIR path
os.environ['indir_path'] = '/scratch/project_465000454/poladesu/FWI/'
os.environ['outdir_path'] = '/scratch/project_465000454/poladesu/FWI/'
os.environ['ct_path'] = '/scratch/project_465000454/poladesu/FWI/'

# Provide the data file name for all variables 
temp_name   ="2m_temperature_DMIN_era5Land.nc"
pr_name     ="total_precipitation_DSUM_era5Land.nc"
uwind_name  ="10m_u_component_of_wind_DMIN_era5Land.nc"
vwind_name  ="10m_v_component_of_wind_DMIN_era5Land.nc"
d2m_name    ="2m_dewpoint_temperature_DMIN_era5Land.nc"
out_name    ="FWI_1900-2010_ERA5LandNEW.nc"
ct_name     ="FWI_Const.nc"

# Combine data path and file names
input_names = 'indir_path ' + 'outdir_path ' + ' ' + temp_name + ' ' + pr_name + ' ' + uwind_name + ' ' + vwind_name + ' ' + d2m_name  + ' ' + out_name + ' ' + 'ct_path '+ ct_name 




# Call F95 excutable and provide the dir path and data file names

subprocess.run(["./main", "f"], text=True, input=str(input_names))
