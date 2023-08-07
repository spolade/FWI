This repository calculates the Canadian Fire Weather Index (FWI), a commonly used index to calculate and assess fire weather. This code is modified from the initial algorithm provided by Simard (1970). Several improvements are applied to use this code globally (reduce the regional dependencies). The FWI code is tested with the ERA-5 Land dataset. However, it can be easily adapted for any model, such as CMIP6 models. These scripts are developed in Fortran 95. Instructions on required packages and how to use the code are provided below.

Requirements on packages:
•	CDI library: is a C and Fortran Interface to access Climate and NWP model Data.
•	NetCDF-fortran/4.5.4

Computation of the FWI:
Moisture codes:  

Subroutine FFMCcalc-  Input: Temp, relative humidity, wind speed, rain, and a previous v     FFMC value 
                      Output: Fine Fuel Moisture Code (FFMC)


Subroutine DMCcalc-    Input: Temp, relative humidity, rain, previous DMC value, latitude, and month 
                       Output: Duff Moisture Code (DMC)


Subroutine DCcalc-     Input: Temp, rain, the previous DC value, latitude, month 
                       Output: Drought Code (DC)


Fire behaviour outputs:

Subroutine ISICalc-   Input: wind speed and current FFMC value 
                      Output: Initial Spread Index (ISI)

Subroutine BUICalc-   Input: the current DMC and DC values 
                      Output: Build-Up Index (BUI)

Subroutine FWICalc-  Input: the current ISI and BUI values 
                     Output: Fire Weather Index (FWI)

