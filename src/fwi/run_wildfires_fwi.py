import os
import numpy as np
import xarray as xr
import pandas as pd
import argparse
import calendar
import fwi

from fwi.functions_cal_fwi import FFMC, DMC, DC, ISI, BUI, FWI

print(
    fwi.__init__
)  # this is to overcome issue with the pip installation of the package

# First step, create a parser:
parser = argparse.ArgumentParser(
    description="Runscript for the wildfires_fwi application"
)


# Second step, add positional arguments or
# https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser.add_argument
parser.add_argument(
    "--year", required=True, help="Input year for the wildfires_fwi app", default=1
)
parser.add_argument(
    "--month", required=True, help="Input month for the wildfires_fwi app", default=2
)
parser.add_argument(
    "--day", required=True, help="Input day for the wildfires_fwi app", default=3
)
parser.add_argument(
    "--hpcrootdir", required=True, help="ROOT directory of the experiment", default=4
)
parser.add_argument(
    "--hpcprojdir", required=True, help="PROJECT directory of the HPC", default=8
)
parser.add_argument(
    "--hpctmpdir",
    required=True,
    help="Input expid for the wildfires_fwi app",
    default=4,
)

# Third step, parse arguments.
# The default args list is taken from sys.args
args = parser.parse_args()

year  = args.year
month = args.month
day   = args.day

hpctmpdir  = args.hpctmpdir
hpcrootdir = args.hpcrootdir
hpcprojdir = args.hpcprojdir

# Provide the input and output DIR path
in_path   = f"{hpctmpdir}/"
out_path  = f"{hpctmpdir}/"
ct_path   = f"{hpctmpdir}/"

# Provide the data file name for all variables
previous_day   = int(args.day) - 1
previous_month = int(args.month) - 1

temp_name    = f"{year}_{month}_{day}_T12_00_2t_raw_data.nc"
pr_name      = f"{year}_{month}_{str(previous_day).zfill(2)}_T13_avg_tprate_timestep_60_daily_noon_mean.nc"
uwind_name   = f"{year}_{month}_{day}_T12_00_10u_raw_data.nc"
vwind_name   = f"{year}_{month}_{day}_T12_00_10v_raw_data.nc"
d2m_name     = f"{year}_{month}_{day}_T12_00_2d_raw_data.nc"
out_name     = f"{year}_{month}_{day}_T12_00_fwi.nc"
restart_name = f"{year}_{month}_{day}_T12_00_fwi_restart_file.nc"
restart_prv  = f"{year}_{str(previous_month).zfill(2)}_{day}_T12_00_fwi_restart_file.nc"
ct_name      = "FWI_Const1.nc"


file_t2m         = os.path.join(in_path, temp_name)
file_pr          = os.path.join(in_path, pr_name)
file_d2m         = os.path.join(in_path, d2m_name)
file_10u         = os.path.join(in_path, uwind_name)
file_10v         = os.path.join(in_path, vwind_name)
out_file         = os.path.join(in_path, out_name)
restart_file     = os.path.join(in_path, restart_name)
restart_file_prv = os.path.join(in_path, restart_prv)
const_file       = os.path.join(in_path, ct_name)

# Read data

ds_d2m = xr.open_dataset(file_d2m)
ds_10u = xr.open_dataset(file_10u)
ds_10v = xr.open_dataset(file_10v)
ds_t2m = xr.open_dataset(file_t2m)

# Read Variable Name
d2m_var_name = list(ds_d2m.data_vars.keys())[0]
u10_var_name = list(ds_10u.data_vars.keys())[0]
v10_var_name = list(ds_10v.data_vars.keys())[0]
t2m_var_name = list(ds_t2m.data_vars.keys())[0]

# Calculate wind speed
# Extract Data only for 12 UTC
ds_10uN = ds_10u.sel(time=ds_10u["time.hour"] == 12)
ds_10vN = ds_10v.sel(time=ds_10v["time.hour"] == 12)
ds_t2mN = ds_t2m.sel(time=ds_t2m["time.hour"] == 12)
ds_d2mN = ds_d2m.sel(time=ds_d2m["time.hour"] == 12)


wspd = np.sqrt(
    (ds_10uN[u10_var_name] * ds_10uN[u10_var_name])
    + (ds_10vN[v10_var_name] * ds_10vN[v10_var_name])
)
wspd.attrs["units"] = ds_10uN[u10_var_name].attrs.get("units", "").lower()
wspd_units = wspd.attrs.get("units", "").lower()

del ds_10u, ds_10v, ds_10uN, ds_10vN


# Access temperature and dew point temperature variables
tas = ds_t2mN[t2m_var_name]  # Current temperature (in Kelvin)
dev = ds_d2mN[d2m_var_name]  # Dew point temperature (in Kelvin)


# Read the PR data
if os.path.exists(file_pr):
    ds_pr = xr.open_dataset(file_pr)
    pr_var_name = list(ds_pr.data_vars.keys())[0]
    # Get the variable's unit
    pr_units = ds_pr[pr_var_name].attrs.get("units", "").lower()

    # Convert precipitation  to mm/day
    if pr_units in ["kg m-2 s-1", "kg/m2/s", "kg m**-2 s**-1"]:  # kg/m²/s to mm/day
        ds_pr[pr_var_name] = ds_pr[pr_var_name] * 86400
        pr_new_units = "mm/day"
    elif pr_units in ["M", "m"]:  # m to mm/day
        ds_pr[pr_var_name] = ds_pr[pr_var_name] * 1000
        pr_new_units = "mm/day"
    elif pr_units in ["mm/day", "mm d-1"]:  # Already in mm/day
        pr_new_units = "mm/day"
    else:
        pr_new_units = "Unknown"  # Handle unknown units
    pr = ds_pr[pr_var_name]
    pr.attrs["units"] = pr_new_units
else:
    pr = xr.full_like(tas, np.nan)


# Calculate saturation vapor pressure from dew point temperature
sat_dew = 6.11 * np.exp(53.49 - 6808 / dev - 5.09 * np.log(dev))

# Calculate saturation vapor pressure from current temperature
sat_tas = 6.11 * np.exp(53.49 - 6808 / tas - 5.09 * np.log(tas))

# Calculate relative humidity
rhum = (sat_dew / sat_tas) * 100



# Correct the units

# Temperature
tas_units = tas.attrs.get("units", "").lower()
if tas_units in ["K", "k", "K deg", "k deg"]:  # #Temperature K == Deg C
    tas = tas - 273.15
    tas_new_units = "C"
elif tas_units in ["C", "C deg"]:  # Already in c
    tas_new_units = "C"
else:
    tas_new_units = "Unknown"  # Handle unknown units

# Wind Speed
if wspd_units in ["m/s", "m s**-1"]:  # "Wind speed m/s == km/h
    wspd = wspd * 3.6
    wspd_new_units = "km/h"
elif wspd_units in ["km/h", "km**-h"]:  # Already in km/h
    wspd_new_units = "km/h"
else:
    wspd_new_units = "Unknown"  # Handle unknown units


rhum = xr.where(rhum > 100, 100, rhum)


FWI_all = xr.full_like(tas, 0.0)


numb_day = sum([pd.Period(f"{year}-{i}-1").daysinmonth for i in range(1, 13)])
numb_day_aligned = xr.full_like(tas.isel(time=0), numb_day)

MONTH_in      = tas.time.dt.month
MONTH_aligned = MONTH_in.broadcast_like(tas)

LAT_in      = tas.lat
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
    def __init__(
        self, type_variables, adjust_DryingFactor, adjust_DayLength, adjust_overwinterDC
    ):
        self.type_variables = type_variables
        self.adjust_DryingFactor = adjust_DryingFactor
        self.adjust_DayLength = adjust_DayLength
        self.adjust_overwinterDC = adjust_overwinterDC


# Create an instance of the Config class with the desired values
cfg = Config("hursmin-tasmax", "NSHeq", "continuous", "original")


for i in range(len(tas.time)):
    # First initialize FFMC, DMC, and DC
    tas_slice = tas.isel(time=i)
    rhum_slice = rhum.isel(time=i)
    wspd_slice = wspd.isel(time=i)
    pr_slice = pr.isel(time=i)
    lat_slice = LAT_aligned.isel(time=i)
    month_slice = MONTH_aligned.isel(time=i)

    if os.path.exists(restart_file):
	    with xr.open_dataset(restart_file) as ds:
		    FFMCPrev = ds["FFMCPrev"].load()
		    DMCPrev  = ds["DMCPrev"].load()
		    DCPrev   = ds["DCPrev"].load()
    elif os.path.exists(const_file):
	    with xr.open_dataset(const_file) as ds:
		    FFMCPrev = ds["FFMCPrev"].load()
		    DMCPrev  = ds["DMCPrev"].load()
		    DCPrev   = ds["DCPrev"].load()
    else:
	    template = tas_slice  # match current timestep shape
	    FFMCPrev = xr.full_like(template, 0.0)
	    DMCPrev  = xr.full_like(template, 0.0)
	    DCPrev   = xr.full_like(template, 0.0)
	    xr.Dataset(
		    {"FFMCPrev": FFMCPrev, "DMCPrev": DMCPrev, "DCPrev": DCPrev}
	    ).to_netcdf(const_file)
    FFMCPrev_in = FFMCPrev.values
    DMCPrev_in = DMCPrev.values
    DCPrev_in = DCPrev.values

    FFMC_in = FFMC(
        tas_slice.values,
        rhum_slice.values,
        wspd_slice.values,
        pr_slice.values,
        FFMCPrev_in,
    )
    DMC_in = DMC(
        tas_slice.values,
        rhum_slice.values,
        pr_slice.values,
        DMCPrev_in,
        lat_slice.values,
        numb_day_aligned.values,
        month_slice.values,
        cfg,
    )
    DC_in = DC(
        tas_slice.values,
        pr_slice.values,
        DCPrev_in,
        lat_slice.values,
        month_slice.values,
        cfg,
    )
    ISI_in = ISI(wspd_slice.values, FFMC_in)
    BUI_in = BUI(DMC_in, DC_in)
    FWI_in = FWI(ISI_in, BUI_in)
    FWI_all[i, :, :] = FWI_in

    FFMC_in_xr = xr.DataArray(
        FFMC_in, coords={"lat": tas.lat, "lon": tas.lon}, dims=["lat", "lon"]
    )
    DMC_in_xr = xr.DataArray(
        DMC_in, coords={"lat": tas.lat, "lon": tas.lon}, dims=["lat", "lon"]
    )
    DC_in_xr = xr.DataArray(
        DC_in, coords={"lat": tas.lat, "lon": tas.lon}, dims=["lat", "lon"]
    )

    dataset = xr.Dataset(
        {"FFMCPrev": FFMC_in_xr, "DMCPrev": DMC_in_xr, "DCPrev": DC_in_xr}
    )
    dataset.to_netcdf(const_file)
    # Write the restart file for the 1st of each month
    if int(tas_slice.time.dt.day) == 1:
        dataset.to_netcdf(restart_file)
        if os.path.exists(restart_file_prv):
            print("File deleted successfully")
            os.remove(restart_file_prv)
    dataset.close()

    del FFMC_in, DMC_in, DC_in, ISI_in, BUI_in, FWI_in

FWI_all.attrs = {"long_name": "Fire Weather Index", "units": " "}
FWI_all = FWI_all.to_dataset(name="fwi")


# Add global attributes
# Copy only the "dataset" attribute if it exists
if "activity" in ds_t2m.attrs:
    FWI_all.attrs["activity"] = ds_t2m.attrs["activity"]
if "dataset" in ds_t2m.attrs:
    FWI_all.attrs["dataset"] = ds_t2m.attrs["dataset"]
if "experiment" in ds_t2m.attrs:
    FWI_all.attrs["experiment"] = ds_t2m.attrs["experiment"]
if "generation" in ds_t2m.attrs:
    FWI_all.attrs["generation"] = ds_t2m.attrs["generation"]
if "type" in ds_t2m.attrs:
    FWI_all.attrs["type"] = ds_t2m.attrs["type"]
if "levtype" in ds_t2m.attrs:
    FWI_all.attrs["levtype"] = ds_t2m.attrs["levtype"]
if "model" in ds_t2m.attrs:
    FWI_all.attrs["model"] = ds_t2m.attrs["model"]
if "class" in ds_t2m.attrs:
    FWI_all.attrs["class"] = ds_t2m.attrs["class"]
if "realization" in ds_t2m.attrs:
    FWI_all.attrs["realization"] = ds_t2m.attrs["realization"]
if "stream" in ds_t2m.attrs:
    FWI_all.attrs["stream"] = ds_t2m.attrs["stream"]
if "resolution" in ds_t2m.attrs:
    FWI_all.attrs["resolution"] = ds_t2m.attrs["resolution"]
if "expver" in ds_t2m.attrs:
    FWI_all.attrs["expver"] = ds_t2m.attrs["expver"]

FWI_all.attrs["history"] = (
    "Created on " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
)

FWI_all.to_netcdf(path=out_file)




checkOutputFile = os.path.join(in_path, f"{year_in}_*_T12_00_fwi.nc")

if os.path.exists(checkOutputFile):
    print("Output for wildfires_fwi application has been generated!")
    print("Inspect output for validity.")
else:
    print("Output for wildfires_fwi application has NOT been generated!")

print("Finished!")




# Create The yearly Stats
# Parse date inputs
year_in  = int(args.year)
month_in = int(args.month)
day_in   = int(args.day)


OUTPUT_DIR    = in_path
FWI_THRESHOLD = 15
WINDOW        = 14




if (month_in == 12 and day_in == 31) or (month_in == 12 and day_in == 30 and not calendar.isleap(year_in)):
    print(f"Triggering yearly statistics for year {year_in}...")

    file_pattern = os.path.join(in_path, f"{year_in}_*_T12_00_fwi.nc")
    file_list = sorted(glob.glob(file_pattern))
    expected_days = 366 if calendar.isleap(year_in) else 365

    if len(file_list) == expected_days:
        ds_combined = xr.open_mfdataset(file_list, combine="by_coords", chunks={"time": 30}).chunk({"time": -1})

        # Mean, Max, Min
        ds_combined.fwi.mean(dim="time").to_dataset(name="fwi").to_netcdf(os.path.join(OUTPUT_DIR, f"{year_in}_fwi_mean.nc"), compute=False).compute()
        ds_combined.fwi.max(dim="time").to_dataset(name="fwi").to_netcdf(os.path.join(OUTPUT_DIR, f"{year_in}_fwi_max.nc"), compute=False).compute()
        ds_combined.fwi.min(dim="time").to_dataset(name="fwi").to_netcdf(os.path.join(OUTPUT_DIR, f"{year_in}_fwi_min.nc"), compute=False).compute()

        # Percentiles
        for perc in [50, 95, 99]:
            percentile = ds_combined.fwi.quantile(perc / 100, dim="time")
            percentile.to_dataset(name="fwi").to_netcdf(os.path.join(OUTPUT_DIR, f"{year_in}_fwi_{perc}per.nc"), compute=False).compute()

        # DSR Monthly
        dsr = 0.0272 * np.power(ds_combined.fwi, 1.77)
        dsr.resample(time='1M').mean().to_dataset(name="fwi").to_netcdf(os.path.join(OUTPUT_DIR, f"{year_in}_fwi_DSR_monthly.nc"), compute=False).compute()

        # Threshold exceedance days
        thresholds = {
            "VLow": (ds_combined.fwi < 5.2),
            "Low": ((ds_combined.fwi >= 5.2) & (ds_combined.fwi < 11.2)),
            "Mod": ((ds_combined.fwi >= 11.2) & (ds_combined.fwi < 21.3)),
            "High": ((ds_combined.fwi >= 21.3) & (ds_combined.fwi < 38.0)),
            "VHigh": ((ds_combined.fwi >= 38.0) & (ds_combined.fwi < 50)),
            "Ext": (ds_combined.fwi >= 50)
        }

        for label, mask in thresholds.items():
            mask.sum(dim="time").to_dataset(name="fwi").to_netcdf(os.path.join(OUTPUT_DIR, f"{year_in}_fwi_days_{label}.nc"), compute=False).compute()


        ds_combined.close()
    else:
        print(f"Incomplete data: {len(file_list)} of {expected_days} files found.")
else:
    print(f"Not the last day of the year ({year_in}-{month_in:02d}-{day_in:02d}). Skipping yearly stats.")


# =====================================================================
