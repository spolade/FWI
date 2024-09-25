import os
import sys  # used for running all at once
import time
import numpy as np
import xarray as xr
import pandas as pd
import argparse

from functions_cal_FWI import *


def cli():
    # First step, create a parser:
    parser = argparse.ArgumentParser(description="Runscript for the wildfires_fwi application")


    # Second step, add positional arguments or
    # https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.add_argument
    parser.add_argument('-year', required=True, help="Input year for the wildfires_fwi app", default=1)
    parser.add_argument('-month', required=True, help="Input month for the wildfires_fwi app", default=2)
    parser.add_argument('-day', required=True, help="Input day for the wildfires_fwi app", default=3)
    parser.add_argument('-hpcrootdir', required=True, help="ROOT directory of the experiment", default=4)
    parser.add_argument('-hpcprojdir', required=True, help="PROJECT directory of the HPC", default=8)
    parser.add_argument('-hpctmpdir', required=True, help="Input expid for the wildfires_fwi app", default=4)

    # Third step, parse arguments.
    # The default args list is taken from sys.args
    args = parser.parse_args()

    year = args.year
    month = args.month
    day = args.day

    hpctmpdir = args.hpctmpdir
    hpcrootdir = args.hpcrootdir
    hpcprojdir = args.hpcprojdir

    # Provide the input and output DIR path
    in_path  = f'{hpctmpdir}/'
    out_path = f'{hpctmpdir}/'
    ct_path  = f'{hpctmpdir}/'

    # Provide the input and output DIR path
    #in_path  = f'{hpcrootdir}/tmp'
    #out_path = f'{hpcrootdir}/tmp'
    #ct_path  = f'{hpcrootdir}/tmp'

    #in_path  = '/scratch/project_465000454/tmp/a15m'
    #out_path = '/scratch/project_465000454/tmp/a15m'
    #ct_path  = '/scratch/project_465000454/tmp/a15m'

    #in_path  = f'{hpcrootdir}/tmp'
    #out_path = f'{hpcrootdir}/tmp'
    #ct_path  = f'{hpcrootdir}/tmp'


    # Provide the data file name for all variables 
    previous_day = int(args.day) - 1

    temp_name   = f'{year}_{month}_{day}_T12_00_2t_raw_data.nc' 
    #pr_name     = f'{year}_{month}_{day}_T13_tp_timestep_60_daily_noon_sum.nc'
    pr_name     = f'{year}_{month}_{str(previous_day)}_T13_tp_timestep_60_daily_noon_sum.nc'
    uwind_name  = f'{year}_{month}_{day}_T12_00_10u_raw_data.nc'
    vwind_name  = f'{year}_{month}_{day}_T12_00_10v_raw_data.nc' 
    d2m_name    = f'{year}_{month}_{day}_T12_00_2d_raw_data.nc'
    out_name    = f'fwi_output_{year}{month}{day}.nc'
    ct_name     = f"FWI_Const1.nc"

    file_t2m   = os.path.join(in_path, temp_name)
    file_pr    = os.path.join(in_path, pr_name)
    file_d2m   = os.path.join(in_path, d2m_name)
    file_10u   = os.path.join(in_path, uwind_name)
    file_10v   = os.path.join(in_path, vwind_name)
    out_file   = os.path.join(in_path, out_name)
    const_file = os.path.join(in_path, ct_name)

    # Read data

    ds_d2m     = xr.open_dataset(file_d2m)
    ds_10u     = xr.open_dataset(file_10u)
    ds_10v     = xr.open_dataset(file_10v)
    ds_t2m     = xr.open_dataset(file_t2m)

    # vwind =ds_vwind.sel(time=ds_vwind['time.hour'] == 12)
    # Calculate wind speed 	

    # Extract Data only for 12 UTC 
    ds_10uN =ds_10u.sel(time=ds_10u['time.hour'] == 12)
    ds_10vN =ds_10v.sel(time=ds_10v['time.hour'] == 12)
    ds_t2mN =ds_t2m.sel(time=ds_t2m['time.hour'] == 12)
    ds_d2mN =ds_d2m.sel(time=ds_d2m['time.hour'] == 12)


    wspd = np.sqrt((ds_10uN['10u'] * ds_10uN['10u']) + (ds_10vN['10v'] * ds_10vN['10v']))
    del  ds_10u, ds_10v, ds_10uN,ds_10vN 
        
        
    # Access temperature and dew point temperature variables
    tas = ds_t2mN['2t']  # Current temperature (in Kelvin)
    dev = ds_d2mN['2d']  # Dew point temperature (in Kelvin)


    if os.path.exists(file_pr):
        # Perform your operation here
        ds_pr      = xr.open_dataset(file_pr)
        ds_pr['tp'] = ds_pr['tp']* 1000
        pr          = ds_pr['tp']
    else:
        pr = xr.full_like(tas, np.nan)


    # Calculate saturation vapor pressure from dew point temperature
    sat_dew = 6.11 * np.exp(53.49 - 6808 / dev - 5.09 * np.log(dev))
        
    # Calculate saturation vapor pressure from current temperature
    sat_tas = 6.11 * np.exp(53.49 - 6808 / tas - 5.09 * np.log(tas))
        
    # Calculate relative humidity
    rhum = (sat_dew / sat_tas) * 100
            
    #Correct the units 
        
    #Temperature K == Deg C
    tas = tas - 273.15
        
    #"Wind speed m/s == km/h
    wspd = wspd * 3.6
        
    rhum = xr.where(rhum > 100, 100, rhum)
        
        
    FWI_all = xr.full_like(tas, 0.0)
        
        
    numb_day = sum([pd.Period(f'{year}-{i}-1').daysinmonth for i in range(1,13)])
    numb_day_aligned = xr.full_like(tas.isel(time=0), numb_day)
        
    MONTH_in = tas.time.dt.month
    MONTH_aligned = MONTH_in.broadcast_like(tas)
        
    LAT_in    = tas.lat
    LAT_aligned = LAT_in.broadcast_like(tas)
        
    # Define conf variables
    # Which variables are used: 'hursmin-tasmax', 'hurs-tasmax'
    type_variables = "hursmin-tasmax"
    # Which type of drying factor: 'original', 'NSH', 'NSHeq'
    adjust_DryingFactor = "NSHeq"
    # Which type of DayLength: 'original', 'bins', 'continuous'
    adjust_DayLength = "continuous"
    # Adjustment with overwintering DC: 'original', 'wDC'
    adjust_overwinterDC = "original"
            
            
    class Config:
        def __init__(self, type_variables, adjust_DryingFactor, adjust_DayLength, adjust_overwinterDC):
            self.type_variables = type_variables
            self.adjust_DryingFactor = adjust_DryingFactor
            self.adjust_DayLength = adjust_DayLength
            self.adjust_overwinterDC = adjust_overwinterDC
                
    # Create an instance of the Config class with the desired values
    cfg = Config("hursmin-tasmax", "NSHeq", "continuous", "original")
        
        
        
    for i in range(len(tas.time)):
        #First initialize FFMC, DMC, and DC
        tas_slice = tas.isel(time=i)
        rhum_slice = rhum.isel(time=i)
        wspd_slice = wspd.isel(time=i)
        pr_slice = pr.isel(time=i)
        lat_slice = LAT_aligned.isel(time=i)
        month_slice = MONTH_aligned.isel(time=i)
        
        if int(tas_slice.time.dt.month) == 1 and int(tas_slice.time.dt.day) ==1:
            FFMCPrev = xr.full_like(tas.isel(time=0), 0.00001)
            DMCPrev = xr.full_like(tas.isel(time=0), 0.00001)
            DCPrev = xr.full_like(tas.isel(time=0), 0.00001)
            dataset = xr.Dataset({
            'FFMCPrev': FFMCPrev,
            'DMCPrev': DMCPrev,
            'DCPrev': DCPrev})
            dataset.to_netcdf(const_file)			
        else:
            ds_con      = xr.open_dataset(const_file)
            FFMCPrev = ds_con.FFMCPrev
            DMCPrev = ds_con.DMCPrev
            DCPrev = ds_con.DCPrev
            os.remove(const_file)
                    
        FFMCPrev_in = FFMCPrev.values
        DMCPrev_in  = DMCPrev.values
        DCPrev_in   = DCPrev.values	  
          

        FFMC_in = FFMC(tas_slice.values, rhum_slice.values, wspd_slice.values, pr_slice.values, FFMCPrev_in)      
        DMC_in =  DMC(tas_slice.values, rhum_slice.values, pr_slice.values, DMCPrev_in, lat_slice.values, numb_day_aligned.values, month_slice.values, cfg) 
        DC_in = DC(tas_slice.values, pr_slice.values, DCPrev_in, lat_slice.values, month_slice.values, cfg) 
        ISI_in = ISI(wspd_slice.values, FFMC_in)
        BUI_in =  BUI(DMC_in, DC_in)
        FWI_in = FWI(ISI_in, BUI_in)
        FWI_all[i,:,:] = FWI_in
            
            
        FFMC_in_xr = xr.DataArray(FFMC_in, coords={'lat': tas.lat, 'lon': tas.lon}, dims=['lat', 'lon'])
        DMC_in_xr  = xr.DataArray(DMC_in, coords={'lat': tas.lat, 'lon': tas.lon}, dims=['lat', 'lon'])
        DC_in_xr   = xr.DataArray(DC_in, coords={'lat': tas.lat, 'lon': tas.lon}, dims=['lat', 'lon'])
        
        dataset    = xr.Dataset({
            'FFMCPrev': FFMC_in_xr,
            'DMCPrev': DMC_in_xr,
            'DCPrev': DC_in_xr})
        dataset.to_netcdf(const_file)
        dataset.close()
                
        del FFMC_in, DMC_in, DC_in, ISI_in, BUI_in, FWI_in  
            
    FWI_all.attrs = {'long_name': 'Fire Weather Index', 'units': 'unit'}
    FWI_all = FWI_all.to_dataset(name='fwi')

    FWI_all .to_netcdf(path=out_file)	

    print("Finished!")


    checkOutputFile = f'{hpcrootdir}/tmp/fwi_output_{year}{month}{day}.nc'



    if os.path.exists(checkOutputFile):
        print('Output for wildfires_fwi application has been generated!')
        print('Inspect output for validity.')
    else:
        print('Output for wildfires_fwi application has NOT been generated!')

    # =====================================================================


if __name__ == "__main__":
    cli()
