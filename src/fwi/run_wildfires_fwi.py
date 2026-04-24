import os
import glob
import argparse
import calendar

import dask
import fwi
import numpy as np
import pandas as pd
import xarray as xr

from fwi.functions_cal_fwi import FFMC, DMC, DC, ISI, BUI, FWI


# ==========================================================
# Controls
# ==========================================================
DO_ANNUAL = True

print(fwi.__init__)  # Helps work around package import issues


# ==========================================================
# Argument parsing
# ==========================================================
parser = argparse.ArgumentParser(
    description="Runscript for the wildfires_fwi application"
)

parser.add_argument("--year", required=True, help="Input year")
parser.add_argument("--month", required=True, help="Input month")
parser.add_argument("--day", required=True, help="Input day")
parser.add_argument("--hpcrootdir", required=True, help="ROOT directory of the experiment")
parser.add_argument("--hpcprojdir", required=True, help="PROJECT directory of the HPC")
parser.add_argument("--hpctmpdir", required=True, help="Directory containing raw meteo input files")

args = parser.parse_args()

year = str(args.year)
month = str(args.month).zfill(2)
day = str(args.day).zfill(2)

year_in = int(args.year)
month_in = int(args.month)
day_in = int(args.day)

hpcrootdir = args.hpcrootdir
hpcprojdir = args.hpcprojdir
hpctmpdir = args.hpctmpdir


# ==========================================================
# Helpers
# ==========================================================
def copy_selected_global_attrs(src_ds, dst_ds):
    keys = [
        "activity",
        "dataset",
        "experiment",
        "generation",
        "type",
        "levtype",
        "model",
        "class",
        "realization",
        "stream",
        "resolution",
        "expver",
    ]
    for key in keys:
        if key in src_ds.attrs:
            dst_ds.attrs[key] = src_ds.attrs[key]
    return dst_ds


def add_history_attr(ds):
    ds.attrs["history"] = "Created on " + pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    return ds


def get_first_var_name(ds):
    return list(ds.data_vars.keys())[0]


def normalize_temperature_to_celsius(da):
    units = da.attrs.get("units", "").strip().lower()
    if units in ["k", "k deg", "kelvin"]:
        out = da - 273.15
        out.attrs["units"] = "C"
        return out
    if units in ["c", "c deg", "degc", "celsius"]:
        da.attrs["units"] = "C"
        return da
    da.attrs["units"] = "Unknown"
    return da


def normalize_wind_to_kmh(da):
    units = da.attrs.get("units", "").strip().lower()
    if units in ["m/s", "m s-1", "m s**-1"]:
        out = da * 3.6
        out.attrs["units"] = "km/h"
        return out
    if units in ["km/h", "km h-1", "km**-h"]:
        da.attrs["units"] = "km/h"
        return da
    da.attrs["units"] = "Unknown"
    return da


def normalize_precip_to_mm_day(da):
    units = da.attrs.get("units", "").strip().lower()
    if units in ["kg m-2 s-1", "kg/m2/s", "kg m**-2 s**-1"]:
        out = da * 86400
        out.attrs["units"] = "mm/day"
        return out
    if units in ["m"]:
        out = da * 1000
        out.attrs["units"] = "mm/day"
        return out
    if units in ["mm/day", "mm d-1"]:
        da.attrs["units"] = "mm/day"
        return da
    da.attrs["units"] = "Unknown"
    return da


def build_output_dir_from_input(raw_input_dir):
    if "/opa/" not in raw_input_dir:
        raise ValueError(
            f"Expected '/opa/' in hpctmpdir, but got: {raw_input_dir}"
        )
    return raw_input_dir.replace("/opa/", "/output/", 1)


class Config:
    def __init__(self, type_variables, adjust_DryingFactor, adjust_DayLength, adjust_overwinterDC):
        self.type_variables = type_variables
        self.adjust_DryingFactor = adjust_DryingFactor
        self.adjust_DayLength = adjust_DayLength
        self.adjust_overwinterDC = adjust_overwinterDC


# ==========================================================
# Paths
# ==========================================================
raw_input_dir = hpctmpdir.rstrip("/")
output_dir = build_output_dir_from_input(raw_input_dir)

os.makedirs(output_dir, exist_ok=True)

in_path = raw_input_dir
out_path = output_dir
ct_path = output_dir

previous_day = day_in - 1
previous_month = month_in - 1

temp_name = f"{year}_{month}_{day}_T12_00_2t_raw_data.nc"
pr_name = f"{year}_{month}_{str(previous_day).zfill(2)}_T13_avg_tprate_timestep_60_daily_noon_mean.nc"
uwind_name = f"{year}_{month}_{day}_T12_00_10u_raw_data.nc"
vwind_name = f"{year}_{month}_{day}_T12_00_10v_raw_data.nc"
d2m_name = f"{year}_{month}_{day}_T12_00_2d_raw_data.nc"
out_name = f"{year}_{month}_{day}_T12_00_fwi.nc"
restart_name = f"{year}_{month}_{day}_T12_00_fwi_restart_file.nc"
restart_prv = f"{year}_{str(previous_month).zfill(2)}_{day}_T12_00_fwi_restart_file.nc"
ct_name = "FWI_Const1.nc"

file_t2m = os.path.join(in_path, temp_name)
file_pr = os.path.join(in_path, pr_name)
file_d2m = os.path.join(in_path, d2m_name)
file_10u = os.path.join(in_path, uwind_name)
file_10v = os.path.join(in_path, vwind_name)

out_file = os.path.join(out_path, out_name)
restart_file = os.path.join(ct_path, restart_name)
restart_file_prv = os.path.join(ct_path, restart_prv)
const_file = os.path.join(ct_path, ct_name)


# ==========================================================
# Read input data
# ==========================================================
ds_d2m = xr.open_dataset(file_d2m)
ds_10u = xr.open_dataset(file_10u)
ds_10v = xr.open_dataset(file_10v)
ds_t2m = xr.open_dataset(file_t2m)

d2m_var_name = get_first_var_name(ds_d2m)
u10_var_name = get_first_var_name(ds_10u)
v10_var_name = get_first_var_name(ds_10v)
t2m_var_name = get_first_var_name(ds_t2m)

ds_10uN = ds_10u.sel(time=ds_10u["time.hour"] == 12)
ds_10vN = ds_10v.sel(time=ds_10v["time.hour"] == 12)
ds_t2mN = ds_t2m.sel(time=ds_t2m["time.hour"] == 12)
ds_d2mN = ds_d2m.sel(time=ds_d2m["time.hour"] == 12)

wspd = np.sqrt(ds_10uN[u10_var_name] ** 2 + ds_10vN[v10_var_name] ** 2)
wspd.attrs["units"] = ds_10uN[u10_var_name].attrs.get("units", "")

tas = ds_t2mN[t2m_var_name]
dev = ds_d2mN[d2m_var_name]

tas = normalize_temperature_to_celsius(tas)
wspd = normalize_wind_to_kmh(wspd)

if os.path.exists(file_pr):
    ds_pr = xr.open_dataset(file_pr)
    pr_var_name = get_first_var_name(ds_pr)
    pr = normalize_precip_to_mm_day(ds_pr[pr_var_name])
else:
    pr = xr.full_like(tas, np.nan)
    pr.attrs["units"] = "mm/day"


# ==========================================================
# Derived variables
# ==========================================================
sat_dew = 6.11 * np.exp(53.49 - 6808 / dev - 5.09 * np.log(dev))
sat_tas = 6.11 * np.exp(53.49 - 6808 / tas - 5.09 * np.log(tas))
rhum = (sat_dew / sat_tas) * 100
rhum = xr.where(rhum > 100, 100, rhum)

FWI_all = xr.full_like(tas, 0.0)

numb_day = sum(pd.Period(f"{year_in}-{i}-1").daysinmonth for i in range(1, 13))
numb_day_aligned = xr.full_like(tas.isel(time=0), numb_day)

MONTH_aligned = tas.time.dt.month.broadcast_like(tas)
LAT_aligned = tas.lat.broadcast_like(tas)

cfg = Config("hursmin-tasmax", "NSHeq", "continuous", "original")


# ==========================================================
# Main daily FWI loop
# ==========================================================
for i in range(len(tas.time)):
    tas_slice = tas.isel(time=i)
    rhum_slice = rhum.isel(time=i)
    wspd_slice = wspd.isel(time=i)
    pr_slice = pr.isel(time=i)
    lat_slice = LAT_aligned.isel(time=i)
    month_slice = MONTH_aligned.isel(time=i)

    if os.path.exists(restart_file):
        with xr.open_dataset(restart_file) as ds_restart:
            FFMCPrev = ds_restart["FFMCPrev"].load()
            DMCPrev = ds_restart["DMCPrev"].load()
            DCPrev = ds_restart["DCPrev"].load()
    elif os.path.exists(const_file):
        with xr.open_dataset(const_file) as ds_const:
            FFMCPrev = ds_const["FFMCPrev"].load()
            DMCPrev = ds_const["DMCPrev"].load()
            DCPrev = ds_const["DCPrev"].load()
    else:
        template = tas_slice
        FFMCPrev = xr.full_like(template, 0.0)
        DMCPrev = xr.full_like(template, 0.0)
        DCPrev = xr.full_like(template, 0.0)
        xr.Dataset(
            {"FFMCPrev": FFMCPrev, "DMCPrev": DMCPrev, "DCPrev": DCPrev}
        ).to_netcdf(const_file)

    FFMC_in = FFMC(
        tas_slice.values,
        rhum_slice.values,
        wspd_slice.values,
        pr_slice.values,
        FFMCPrev.values,
    )
    DMC_in = DMC(
        tas_slice.values,
        rhum_slice.values,
        pr_slice.values,
        DMCPrev.values,
        lat_slice.values,
        numb_day_aligned.values,
        month_slice.values,
        cfg,
    )
    DC_in = DC(
        tas_slice.values,
        pr_slice.values,
        DCPrev.values,
        lat_slice.values,
        month_slice.values,
        cfg,
    )
    ISI_in = ISI(wspd_slice.values, FFMC_in)
    BUI_in = BUI(DMC_in, DC_in)
    FWI_in = FWI(ISI_in, BUI_in)

    FWI_all[i, :, :] = FWI_in

    ffmc_da = xr.DataArray(FFMC_in, coords={"lat": tas.lat, "lon": tas.lon}, dims=["lat", "lon"])
    dmc_da = xr.DataArray(DMC_in, coords={"lat": tas.lat, "lon": tas.lon}, dims=["lat", "lon"])
    dc_da = xr.DataArray(DC_in, coords={"lat": tas.lat, "lon": tas.lon}, dims=["lat", "lon"])

    restart_ds = xr.Dataset({"FFMCPrev": ffmc_da, "DMCPrev": dmc_da, "DCPrev": dc_da})
    restart_ds.to_netcdf(const_file)

    if int(tas_slice.time.dt.day) == 1:
        restart_ds.to_netcdf(restart_file)
        if os.path.exists(restart_file_prv):
            os.remove(restart_file_prv)

    del FFMC_in, DMC_in, DC_in, ISI_in, BUI_in, FWI_in


# ==========================================================
# Write daily FWI file
# ==========================================================
FWI_all.attrs = {"long_name": "Fire Weather Index", "units": " "}
FWI_all = FWI_all.astype("float32")
FWI_all = FWI_all.to_dataset(name="fwi")

FWI_all = copy_selected_global_attrs(ds_t2m, FWI_all)
FWI_all = add_history_attr(FWI_all)

encoding = {
    "fwi": {
        "dtype": "float32",
        "zlib": True,
        "complevel": 5,
        "shuffle": True
    }
}

FWI_all.to_netcdf(
    path=out_file,
    format="NETCDF4",
    engine="netcdf4",
    encoding=encoding
)


# ==========================================================
# Annual statistics
# ==========================================================
OUTPUT_DIR = out_path

dask.config.set({"array.slicing.split_large_chunks": True})

if DO_ANNUAL and month_in == 12 and day_in == 31:
    print(f"Triggering yearly statistics for year {year_in}...")

    file_pattern = os.path.join(out_path, f"{year_in}_*_T12_00_fwi.nc")
    file_list = sorted(glob.glob(file_pattern))
    expected_days = 366 if calendar.isleap(year_in) else 365

    if len(file_list) != expected_days:
        print(f"Incomplete data: {len(file_list)} of {expected_days} files found.")
    else:
        ds_year = xr.open_mfdataset(
            file_list,
            combine="by_coords",
            parallel=False,
            data_vars="minimal",
            coords="minimal",
            compat="override",
            chunks={"time": 15},
            engine="netcdf4",
        )

        fwi_year = ds_year["fwi"]
        spatial_dims = [d for d in fwi_year.dims if d != "time"]

        chunk_map = {"time": 15}
        for dim in spatial_dims:
            chunk_map[dim] = 200

        fwi_year = fwi_year.chunk(chunk_map).astype("float32")

        def enc(da):
            return {
                "zlib": True,
                "complevel": 4,
                "dtype": "float32",
                "chunksizes": tuple(da.chunksizes[d][0] for d in da.dims),
            }

        def enc_int(_da):
            return {
                "zlib": True,
                "complevel": 4,
                "dtype": "int16",
            }

        annual_basic = xr.Dataset(
            {
                "fwi_mean": fwi_year.mean("time"),
                "fwi_max": fwi_year.max("time"),
                "fwi_min": fwi_year.min("time"),
                "days_VLow": (fwi_year < 5.2).sum("time").astype("int16"),
                "days_Low": ((fwi_year >= 5.2) & (fwi_year < 11.2)).sum("time").astype("int16"),
                "days_Mod": ((fwi_year >= 11.2) & (fwi_year < 21.3)).sum("time").astype("int16"),
                "days_High": ((fwi_year >= 21.3) & (fwi_year < 38.0)).sum("time").astype("int16"),
                "days_VHigh": ((fwi_year >= 38.0) & (fwi_year < 50.0)).sum("time").astype("int16"),
                "days_Ext": (fwi_year >= 50.0).sum("time").astype("int16"),
            }
        )

        annual_basic.attrs.update(FWI_all.attrs)
        annual_basic = add_history_attr(annual_basic)

        annual_basic_path = os.path.join(OUTPUT_DIR, f"{year_in}_fwi_basic_and_days.nc")
        annual_basic.to_netcdf(
            annual_basic_path,
            format="NETCDF4",
            engine="netcdf4",
            encoding={
                "fwi_mean": enc(annual_basic["fwi_mean"]),
                "fwi_max": enc(annual_basic["fwi_max"]),
                "fwi_min": enc(annual_basic["fwi_min"]),
                "days_VLow": enc_int(annual_basic["days_VLow"]),
                "days_Low": enc_int(annual_basic["days_Low"]),
                "days_Mod": enc_int(annual_basic["days_Mod"]),
                "days_High": enc_int(annual_basic["days_High"]),
                "days_VHigh": enc_int(annual_basic["days_VHigh"]),
                "days_Ext": enc_int(annual_basic["days_Ext"]),
            },
            compute=True,
        )

        quantile_path = os.path.join(OUTPUT_DIR, f"{year_in}_fwi_quantiles.nc")
        print(f"Computing quantiles -> {quantile_path}", flush=True)

        chunk_map_q = {"time": -1}
        for dim in [d for d in fwi_year.dims if d != "time"]:
            chunk_map_q[dim] = 100

        fwi_q = fwi_year.chunk(chunk_map_q).astype("float32")
        q = fwi_q.quantile([0.50, 0.95, 0.99], dim="time").astype("float32")

        q_ds = q.to_dataset(name="fwi")
        q_ds.attrs.update(FWI_all.attrs)
        q_ds = add_history_attr(q_ds)

        q_ds.to_netcdf(
            quantile_path,
            format="NETCDF4",
            engine="netcdf4",
            encoding={
                "fwi": {
                    "zlib": True,
                    "complevel": 4,
                    "dtype": "float32",
                    "chunksizes": tuple(q_ds["fwi"].chunksizes[d][0] for d in q_ds["fwi"].dims),
                }
            },
            compute=True,
        )

        print("Quantiles written OK.", flush=True)

        ds_year.close()
else:
    if not DO_ANNUAL:
        print("Annual statistics disabled in script.")
    else:
        print(
            f"Not the last day of the year ({year_in}-{month_in:02d}-{day_in:02d}). "
            "Skipping yearly stats."
        )
