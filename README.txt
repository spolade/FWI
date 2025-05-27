DT CLIMATE - Wildfire FWI USE CASE
This repository calculates the Canadian Fire Weather Index (FWI), a commonly used index to calculate and assess fire weather. This code is modified from the initial algorithm provided by Van Wagner and Pickett in 1985 (Van Wagner, C.E.; Pickett, T.L. (1985) Equations and FORTRAN program for the Canadian Forest Fire Weather Index SystemCanadian Forestry Service, Petawawa National Forestry Institute, Chalk River, Ontario. Forestry Technical Report 33. 18 p.). Several improvements are applied to use this code globally (reduce the regional dependencies). The FWI code is tested with the ERA-5 Land dataset. However, it can be easily adapted for any model, such as CMIP6 models. These scripts are developed in Fortran 95. Instructions on required packages and how to use the code are provided below.


Computation of the FWI:
Functions used for the calculation of the Canadian Fire Weather Index (FWI)

The FWI System outputs three moisture codes: Fine Fuel Moisture Code (FFMC); Duff Moisture Code (DMC); Drought Code (DC)
 & three fire behavior outputs: Initial Spread Index (ISI); Build-Up Index (BUI); Fire Weather Index (FWI)


Moisture codes:

Function FFMC-  Input: Temp, relative humidity, wind speed, rain, and a previous FFMC value 
                       Output Fine Fuel Moisture Code (FFMC)
                       call FFMCcalc(17,42,25,0,85) = 87.692980092774448

Function DMC-    Input: Temp, relative humidity, rain, previous DMC value, latitude, and month 
                       OutPut: Duff Moisture Code (DMC)
                       call DMCcalc(17,42,0,6,45.98,4) = 8.5450511359999997
Function DC-     Input: Temp, rain, the previous DC value, latititude, month 
                       Output: Drought Code (DC)
                       call DCcalc(17,0,15,45.98,4) = 19.013999999999999


Fire behavior outputs:

Function ISI-    Input: wind speed and current FFMC value 
                       Output: Initial Spread Index (ISI)
                       call ISICalc(25,87.692980092774448) = 10.853661073655068
Function BUI-    Input: the current DMC and DC values 
                      Output: Build-Up Index (BUI)
                       call BUICalc(8.5450511359999997,19.013999999999999) = 8.4904265358371838

Function FWI-    Input: the current ISI and BUI values 
                       Output: Fire Weather Index (FWI)
