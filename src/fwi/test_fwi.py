# Code to unit testing the FWI calculation modules with the data provided by Wang et al. 2015. The code tests each module output with the expected output.
# How to test
# pytest test_fwi.py --maxfail=5 --disable-warnings > test_fwi_output.log 2>&1
# To write the success message
# python add_success_message.py

import pytest
import fwi
import numpy as np
import xarray as xr
import pandas as pd

from fwi.functions_cal_fwi import FFMC, DMC, DC, ISI, BUI, FWI


### Sample inputs and expected outputs
## Which variables are used: 'hursmin-tasmax', 'hurs-tasmax'
## Which type of drying factor: 'original', 'NSH', 'NSHeq'
## Which type of DayLength: 'original', 'bins', 'continuous'
## Adjustment with overwintering DC: 'original', 'wDC'


# Read the sample data file 
df = pd.read_csv('/Users/poladesu/D_Drive/Prog_FMI/dmin_DtEarth/FWI_test_data.csv', delimiter=';', skiprows=0, names=["Lat", "Month", "Day", "Temp", "RH", "Wind", "Rain", "FFMC", "DMC", "DC", "ISI", "BUI", "FWI"])

# Var LIST: "Lat", "Month", "Day", "Temp.", "RH", "Wind", "Rain", "FFMC", "DMC", "DC", "ISI", "BUI", "FWI"])

LAT_1d      = df['Lat'].values
MONTH_1d    = df['Month'].values
numb_day_1d = df['Day'].values
TEMP_1d     = df['Temp'].values
RH_1d       = df['RH'].values
WIND_1d     = df['Wind'].values
RAIN_1d     = df['Rain'].values

LAT      = LAT_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
MONTH    = MONTH_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
numb_day = numb_day_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
TEMP     = TEMP_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
RH       = RH_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
WIND     = WIND_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
RAIN     = RAIN_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)

expected_ffmc_1d  = df['FFMC'].values
expected_dmc_1d   = df['DMC'].values
expected_dc_1d    = df['DC'].values
expected_isi_1d   = df['ISI'].values
expected_bui_1d   = df['BUI'].values
expected_fwi_1d   = df['FWI'].values
   
expected_ffmc      = expected_ffmc_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
expected_dmc      = expected_dmc_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
expected_dc      = expected_dc_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
expected_isi      = expected_isi_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
expected_bui      = expected_bui_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)
expected_fwi      = expected_fwi_1d.reshape(len(TEMP_1d), 1).repeat(2, axis=1)

    
class Config:
    def __init__(self, type_variables, adjust_DryingFactor, adjust_DayLength, adjust_overwinterDC):
        self.type_variables = type_variables
        self.adjust_DryingFactor = adjust_DryingFactor
        self.adjust_DayLength = adjust_DayLength
        self.adjust_overwinterDC = adjust_overwinterDC

# Create an instance of the Config class with the desired values
#cfg = Config("hursmin-tasmax", "NSHeq", "continuous", "original")
cfg = Config("hursmin-tasmax", "NSHeq", "original", "original")


FFMCPrev = np.full_like(TEMP, 85.0)
DMCPrev = np.full_like(TEMP, 6.0)
DCPrev = np.full_like(TEMP, 15.0)


def test_ffmc():
    results = []  # To store calculated FFMC values
    ffmc_prev = FFMCPrev[0]  # Initialize with the first FFMCPrev value

    for i in range(len(TEMP)):
        result = FFMC(TEMP[i,], RH[i,], WIND[i,], RAIN[i,], ffmc_prev)
        results.append(result)
        ffmc_prev = result  # Use the current result as FFMCPrev for the next step
        
        # Check the result with pytest
        assert pytest.approx(result, 0.1) == expected_ffmc[i], f"FFMC test failed at index {i}, got {result}"  
        
    print("All FFMC tests passed.")

    
def test_dmc():
    results = []  # To store calculated DMC values
    dmc_prev = DMCPrev[0]  # Initialize with the first DMCPrev value

    for i in range(len(TEMP)):
        result = DMC(TEMP[i,], RH[i,], RAIN[i,], dmc_prev, LAT[i,], numb_day[i,1], MONTH[i,1], cfg)
        results.append(result)
        dmc_prev = result  # Use the current result as DMCPrev for the next step

        # Check the result with pytest
        assert pytest.approx(result, 0.1) == expected_dmc[i], f"DMC test failed at index {i}, got {result}"

    print("All DMC tests passed.")


def test_dc():
    results = []  # To store calculated DC values
    dc_prev = DCPrev[0]  # Initialize with the first DCPrev value

    for i in range(len(TEMP)):
        result = DC(TEMP[i,], RAIN[i,], dc_prev, LAT[i,], MONTH[i,], cfg)
        results.append(result)
        dc_prev = result  # Use the current result as DCPrev for the next step

        # Check the result with pytest
        assert pytest.approx(result, 0.1) == expected_dc[i], f"DC test failed at index {i}, got {result}"

    print("All DC tests passed.")

    
def test_bui():
    results = []  # To store calculated BUI values
    dc_prev = DCPrev[0]  # Initialize with the first DCPrev value
    dmc_prev = DMCPrev[0]  # Initialize with the first DMCPrev value

    for i in range(len(TEMP)):
        # Calculate DMC and DC for the current time step
        dmc1 = DMC(TEMP[i,], RH[i,], RAIN[i,], dmc_prev, LAT[i,], numb_day[i, 1], MONTH[i, 1], cfg)
        dc1  = DC(TEMP[i,], RAIN[i,], dc_prev, LAT[i,], MONTH[i,], cfg)
        dmc = np.round(dmc1, 1)
        dc = np.round(dc1, 1)
        
        # Calculate BUI using DMC and DC
        result1 = BUI(dmc, dc)
        result = np.round(result1, 1)
        results.append(result)
            
        # Update prev values for the next iteration
        dmc_prev = dmc  # Use the current DMC result as the next dmc_prev
        dc_prev = dc   # Use the current DC result as the next dc_prev
        print(dmc, dc, result)
        
       # Check the result with pytest
        assert pytest.approx(result, 0.1) == expected_bui[i], f"BUI test failed at index {i}, got {result}"

    print("All BUI tests passed.")


def test_isi():
    results = []  # To store calculated BUI values
    ffmc_prev = FFMCPrev[0]  # Initialize with the first FFMCPrev value

    for i in range(len(TEMP)):
        ffmc = FFMC(TEMP[i,], RH[i,], WIND[i,], RAIN[i,], ffmc_prev)
        wind_in = WIND[i,]
        # Calculate ISS using Wind and FFMC
        result = ISI(wind_in, ffmc)
        results.append(result)

        # Update prev values for the next iteration
        ffmc_prev = ffmc
        print(round(result[0],1))
        print(expected_isi[i,1])
        
    	# Check the result with pytest
        assert pytest.approx(round(result[0],1), 0.1) == expected_isi[i,1], f"ISI test failed at index {i}, got {result}"

    print("All ISI tests passed.")


def test_fwi():
    results = []  # To store calculated BUI values
    dc_prev = DCPrev[0]  # Initialize with the first DCPrev value
    dmc_prev = DMCPrev[0]  # Initialize with the first DMCPrev value
    ffmc_prev = FFMCPrev[0]  # Initialize with the first FFMCPrev value

    for i in range(len(TEMP)):
        # Calculate DMC, DC, and FFMC for the current time step
        dmc = DMC(TEMP[i,], RH[i,], RAIN[i,], dmc_prev, LAT[i,], numb_day[i, 1], MONTH[i, 1], cfg)
        dc  = DC(TEMP[i,], RAIN[i,], dc_prev, LAT[i,], MONTH[i,], cfg)
        ffmc = FFMC(TEMP[i,], RH[i,], WIND[i,], RAIN[i,], ffmc_prev)
        
        # Calculate BUI using DMC and DC
        bui = BUI(dmc, dc)
        wind_in = WIND[i,]
        # Calculate ISS using Wind and FFMC
        isi = ISI(wind_in, ffmc)
        result = FWI(isi,bui)
        
        results.append(result)

        # Update prev values for the next iteration
        dmc_prev = dmc  # Use the current DMC result as the next dmc_prev
        dc_prev = dc   # Use the current DC result as the next dc_prev
        ffmc_prev = ffmc # Use the current FFMC result as the next ffmc_prev
        
        # Check the result with pytest
        assert pytest.approx(round(result[0],1), 0.1) == expected_fwi[i,1], f"FWI test failed at index {i}, got {result}"

    print("All FWI tests passed.")

