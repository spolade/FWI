"""
FWI run script for ClimateDT 

"""
import argparse
import calendar
import glob
import os
import re
import sys
import tempfile
from time import perf_counter
import fcntl
import hashlib
import json
from collections import OrderedDict
from contextlib import contextmanager
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# Allow execution directly from this repository without installing it.
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from fwi.functions_cal_fwi import FFMC, DMC, DC, ISI, BUI, FWI, wDC

# Regular FWI and components are independent of this additional product.
DO_OVERWINTERING = True  # False: regular products only; True: also *_fwi_OverWintering.nc.
OVERWINTER_CARRY_OVER = 0.75
OVERWINTER_WETTING_EFFICIENCY = 0.75
OVERWINTER_START_TEMP = 12.0  # Strictly above, in degrees Celsius.
OVERWINTER_END_TEMP = 5.0  # Strictly below, in degrees Celsius.
OVERWINTER_CONSECUTIVE_DAYS = 3  # Completed local noon-to-noon maxima; no future days.

EXCLUDE_FIRST_DAY = True  # Automatically exclude the first date of a fresh integration once.
MONTHLY_CHUNK = False
DO_ANNUAL = True  # FWI annual statistics, independent of component annual statistics.
DO_ANNUAL_OVERWINTERING = True  # Independent annual statistics for the additional seasonal FWI.
DO_ANNUAL_COMPONENTS = True  # Component annual statistics, independent of DO_ANNUAL.
WRITE_COMPONENTS = True  # True: also write YYYY_MM_DD_fwi_comp.nc.
LSM_FILE = None  # Initial mask only; saved state takes priority on subsequent days.
LAND_THRESHOLD = 0.5  # Land fraction >= threshold; 0 means any positive land fraction.
FINAL_LOCAL_DATE = None  # Optional last complete local date, e.g. "2050-12-30".
# An explicit terminal date permits labeled partial-year statistics; never infer it from missing data.
PRECIP_TIME_LABEL = "end"  # Confirmed input convention: 10:00 means 09:00–10:00 UTC.
# Used only when CF time bounds are absent; None requires bounds or a CLI override.
NOON_CONVENTION = "local_mean_solar_noon_longitude_v1"
VARIABLES = {"tas": "2t", "dew": "2d", "u": "10u", "v": "10v", "pr": "avg_tprate"}

# Serial I/O only: bounded rectangular reads in ORIGINAL file longitude order.
# Budget covers one decoded read block, not total application memory; decoding
# copies, prepared global fields and FWI arrays also require RAM.
IO_TILE_COLUMNS = 128
IO_TILE_ROWS = 64  # Contiguous (time,lat,lon) rain/Tmax files: read full-width latitude strips.
IO_BUFFER_MB = 128.0
HANDOVER_CONTEXT = None  # Set only while preparing the transition weather.
COMPATIBILITY_KEYS = ("model", "realization", "resolution")

REPORT_TIMINGS = True  # Console timings only; no output metadata changes.


def report_timing(stage, started, detail=''):
    if REPORT_TIMINGS:
        print(f'[timing] {stage}: {perf_counter() - started:.2f} s{detail}', flush=True)

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
    raise ValueError(f"Unsupported units: {units!r}")


def normalize_wind_to_kmh(da):
    units = da.attrs.get("units", "").strip().lower()
    if units in ["m/s", "m s-1", "m s**-1"]:
        out = da * 3.6
        out.attrs["units"] = "km/h"
        return out
    if units in ["km/h", "km h-1", "km**-h"]:
        da.attrs["units"] = "km/h"
        return da
    raise ValueError(f"Unsupported units: {units!r}")


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
    raise ValueError(f"Unsupported units: {units!r}")


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



def temperature_kelvin(da):
    units = da.attrs.get("units", "").strip().lower()
    if units in ("k", "k deg", "kelvin"):
        return da
    if units in ("c", "c deg", "degc", "celsius"):
        return da + 273.15
    raise ValueError(f"Unsupported temperature units: {units!r}")


def relative_humidity(tas, dew):
    t, td = temperature_kelvin(tas), temperature_kelvin(dew)
    if (t.values <= 0).any() or (td.values <= 0).any():
        raise ValueError("Temperature and dewpoint must be above absolute zero")
    sat = lambda x: 6.11 * np.exp(53.49 - 6808 / x - 5.09 * np.log(x))
    return (100 * sat(td) / sat(t)).clip(min=0, max=100)


def chunk_dates(year, month, day, monthly):
    if monthly:
        start = 1 if day is None else day
        date(year, month, start)  # Validate an optional resume day.
        return [date(year, month, d) for d in range(start, calendar.monthrange(year, month)[1] + 1)]
    if day is None:
        raise ValueError("--day is required in daily mode")
    return [date(year, month, day)]



def metadata_path(path):
    """Legacy companion path, used only when reading old products."""
    return Path(str(path) + '.metadata.json')


def file_digest(path):
    """Validate a legacy sidecar against the exact old NetCDF contents."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _validated_metadata(payload, context):
    """Require a finite JSON dictionary; reject malformed private state."""
    if not isinstance(payload, dict) or not all(isinstance(key, str) for key in payload):
        raise ValueError(f'Invalid FWI runtime metadata dictionary: {context}')
    try:
        encoded = json.dumps(payload, default=_json_value, allow_nan=False, sort_keys=True)
        return json.loads(encoded)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'Invalid FWI runtime metadata: {context}: {exc}') from exc


def read_product_metadata(path, digest=None):
    """Read metadata stored atomically inside NetCDF, or a validated old sidecar."""
    from netCDF4 import Dataset
    with Dataset(path, 'r') as root:
        if '_fwi_runtime' in root.groups:
            group = root.groups['_fwi_runtime']
            if getattr(group, 'format', None) != 'fwi_runtime_metadata_v1':
                raise ValueError(f'Unknown FWI runtime metadata format: {path}')
            try:
                payload = json.loads(group.getncattr('metadata_json'))
            except (AttributeError, TypeError, ValueError) as exc:
                raise ValueError(f'Malformed FWI runtime metadata: {path}') from exc
            return _validated_metadata(payload, path)
    sidecar = metadata_path(path)
    if not sidecar.exists():
        return {}  # Legacy products may keep operational fields in root attributes.
    try:
        with sidecar.open() as stream:
            payload = json.load(stream)
    except (OSError, ValueError) as exc:
        raise ValueError(f'Malformed legacy FWI metadata: {sidecar}') from exc
    if not isinstance(payload, dict) or payload.get('format') != 'fwi_product_metadata_v1':
        raise ValueError(f'Unknown FWI metadata format: {sidecar}')
    versions = payload.get('versions')
    digest = file_digest(path) if digest is None else digest
    metadata = versions.get(digest) if isinstance(versions, dict) else None
    if not isinstance(metadata, dict):
        raise ValueError(f'FWI metadata does not match the NetCDF file: {sidecar}')
    return _validated_metadata(metadata, sidecar)


@contextmanager
def open_product_dataset(path, **kwargs):
    """Restore runtime attributes and optional seasonal state only in memory."""
    from netCDF4 import Dataset
    metadata = read_product_metadata(path)
    with Dataset(path, 'r') as root:
        has_overwinter_state = '_fwi_overwintering' in root.groups
    with xr.open_dataset(path, **kwargs) as ds:
        if not metadata and re.fullmatch(r'\d{4}_\d{2}_\d{2}_fwi\.nc', Path(path).name):
            metadata = daily_metadata_from_constant(path, ds)
        if any(key in ds.attrs and ds.attrs[key] != value for key, value in metadata.items()):
            raise ValueError(f'Conflicting public and private FWI metadata: {path}')
        ds.attrs.update(metadata)
        if has_overwinter_state:
            with xr.open_dataset(path, engine='netcdf4', group='_fwi_overwintering') as state:
                ds.encoding['_fwi_overwinter_state'] = state.load()
        yield ds


def daily_metadata_from_constant(path, ds):
    """Recover annual-validation settings without adding groups to daily FWI.

    The integration constant fixes sampling and mask threshold for the run.
    Verify the daily file's actual grid and UTC noon coordinates rather than
    assuming a filename implies a particular sampling convention.
    """
    constant = Path(path).parent / 'FWI_Const1.nc'
    if not constant.is_file():
        raise FileNotFoundError(f'Annual/resume validation of {path} requires {constant}; '
                                'the daily FWI file intentionally contains no runtime group')
    settings = read_product_metadata(constant)
    stamp = annual_file_date(path)
    with xr.open_dataset(constant) as saved:
        settings = dict(saved.attrs, **settings)
        if settings.get('noon_convention') != NOON_CONVENTION or settings.get('sampling') not in ('linear', 'nearest'):
            raise ValueError(f'Invalid local-noon runtime settings in {constant}')
        if date.fromisoformat(settings['state_date']) < stamp:
            raise ValueError(f'Constant precedes daily FWI product: {path}')
        for name in ('lat', 'lon'):
            if name not in ds or not np.array_equal(ds[name].values, saved[name].values):
                raise ValueError(f'Daily FWI grid differs from integration constant: {path}')
        if ('lsm' in ds) != ('lsm' in saved):
            raise ValueError(f'Daily FWI mask presence differs from integration constant: {path}')
    if 'local_noon_utc' not in ds or ds.local_noon_utc.dims != ('lon',):
        raise ValueError(f'Missing local-noon UTC coordinates in daily FWI: {path}')
    hours = local_noon_hours(stamp, ds.lon.values, settings['sampling'])
    expected = np.rint(hours * 3.6e12).astype(np.int64).astype('datetime64[ns]')
    if not np.array_equal(ds.local_noon_utc.values, expected):
        raise ValueError(f'Daily FWI noon coordinates differ from integration sampling: {path}')
    return {key: settings[key] for key in ('noon_convention', 'sampling', 'land_threshold') if key in settings}


def _json_value(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f'Unsupported metadata value: {type(value).__name__}')


def atomic_write(ds, path, overwinter_state=None, **kwargs):
    """Atomically publish root products and private runtime groups in one NetCDF.

    Public global attributes and variables retain their original contract.
    Constant/restart metadata and optional overwintering arrays live in private
    NetCDF groups. Daily *_fwi.nc omits the runtime group and annual validation
    obtains its integration settings from FWI_Const1.nc in the same directory.
    Old JSON sidecars are accepted for migration but are never created here.
    """
    started = perf_counter()
    from netCDF4 import Dataset
    path = Path(path)
    if kwargs.get('engine', 'netcdf4') != 'netcdf4':
        raise ValueError('Private FWI state requires the netcdf4 engine')
    if kwargs.get('format', 'NETCDF4') != 'NETCDF4':
        raise ValueError('Private FWI state requires NETCDF4 format')
    if kwargs.get('compute', True) is not True:
        raise ValueError('Atomic FWI writes require compute=True')
    if kwargs.get('mode', 'w') != 'w' or kwargs.get('group') is not None:
        raise ValueError('Atomic FWI products must replace the NetCDF root')
    kwargs.update(engine='netcdf4', format='NETCDF4')
    public = ds.copy(deep=False)
    public.encoding = {key: value for key, value in public.encoding.items()
                       if not key.startswith('_fwi_')}
    public.attrs = {}
    copy_selected_global_attrs(ds, public)
    if 'lsm' in public:
        encoding = {name: dict(options) for name, options in kwargs.get('encoding', {}).items()}
        encoding['lsm'] = {'zlib': True, 'complevel': 4, 'dtype': 'float64'}
        kwargs['encoding'] = encoding
    operational = _validated_metadata(
        {key: value for key, value in ds.attrs.items()
         if key not in public.attrs and key != 'history'}, path)
    # History remains the only generated public global attribute.
    public.attrs['history'] = 'Created on ' + pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S.%f')
    if overwinter_state is not None and not isinstance(overwinter_state, xr.Dataset):
        raise TypeError('overwinter_state must be an xarray Dataset')
    fd, name = tempfile.mkstemp(prefix=path.name + '.', suffix='.tmp', dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        public.to_netcdf(temp, **kwargs)
        if overwinter_state is not None:
            state = overwinter_state.copy(deep=False)
            state.encoding = {key: value for key, value in state.encoding.items()
                              if not key.startswith('_fwi_')}
            state_encoding = {name: {'zlib': True, 'complevel': 4} for name in state.data_vars}
            state.to_netcdf(temp, engine='netcdf4', format='NETCDF4', mode='a',
                            group='_fwi_overwintering', encoding=state_encoding)
        # Daily FWI is a clean public product. Integration settings remain in
        # the constant/restarts and are recovered for annual validation there.
        if not re.fullmatch(r'\d{4}_\d{2}_\d{2}_fwi\.nc', path.name):
            with Dataset(temp, 'a') as root:
                group = root.createGroup('_fwi_runtime')
                group.setncattr('format', 'fwi_runtime_metadata_v1')
                group.setncattr('metadata_json', json.dumps(operational, allow_nan=False, sort_keys=True))
        os.replace(temp, path)
        report_timing(f'write {path.name}', started)
    finally:
        if temp.exists():
            temp.unlink()


# --- Hourly local-mean-solar-noon processing ---

def local_noon_hours(stamp, longitude, sampling="linear"):
    """UTC hours since epoch for 12:00 local mean solar time.

    Canonical longitude is [-180, 180); +180 follows the -180 date convention.
    This is longitude time, with no political time zones or equation of time.
    """
    lon = (np.asarray(longitude, dtype=float) + 180.0) % 360.0 - 180.0
    midnight = pd.Timestamp(stamp).value / 3.6e12
    hours = midnight + 12.0 - lon / 15.0
    if sampling == "nearest":
        hours = np.floor(hours + 0.5)  # Half-hour ties go to the later UTC hour.
    return hours


def check_grid(da, path):
    if set(da.dims) != {"time", "lat", "lon"}:
        raise ValueError(f"Expected only time, lat, lon dimensions: {path}")
    for coord in ("lat", "lon"):
        if coord not in da.coords or da[coord].dims != (coord,):
            raise ValueError(f"Explicit 1D {coord} coordinates required: {path}")
        values = da[coord].values
        if len(values) == 0 or not np.isfinite(values).all() or len(np.unique(values)) != len(values):
            raise ValueError(f"Invalid or duplicate {coord}: {path}")
    if (np.abs(da.lat.values) > 90).any():
        raise ValueError(f"Latitude outside [-90,90]: {path}")
    if ((da.lon.values < -180) | (da.lon.values > 360)).any():
        raise ValueError(f"Longitude outside supported [-180,360] convention: {path}")
    normalized = (da.lon.values + 180) % 360 - 180
    if len(np.unique(normalized)) != len(normalized):
        raise ValueError(f"Duplicate cyclic longitude (e.g. both 0 and 360): {path}")


def utc_hours(values, path):
    arr = np.asarray(values)
    if not np.issubdtype(arr.dtype, np.datetime64) or np.isnat(arr).any():
        raise ValueError(f"Decoded Gregorian UTC datetime timestamps required: {path}")
    return arr.astype('datetime64[ns]').astype(np.int64) / 3.6e12


BOUNDARY_WAIT_HOUR = None


class BoundaryInputNotReady(FileNotFoundError):
    """Only the next UTC day's boundary input is unavailable."""


def missing_weather_hour(hour):
    message = f'Missing required UTC weather hour {pd.to_datetime(hour, unit="h")}'
    if BOUNDARY_WAIT_HOUR is not None and hour >= BOUNDARY_WAIT_HOUR:
        return BoundaryInputNotReady(message)
    return FileNotFoundError(message)


def missing_rain_window(source, start, end, message):
    if BOUNDARY_WAIT_HOUR is not None:
        cursor = start
        gaps = []
        for left, right, _ in source.intervals:
            if right <= cursor or left >= end:
                continue
            if left > cursor + 1e-7:
                gaps.append((cursor, min(left, end)))
            cursor = max(cursor, min(right, end))
        if cursor < end - 1e-7:
            gaps.append((cursor, end))
        # End-labelled next-day 00:00 supplies the preceding 23:00–00:00 hour.
        if gaps and all(left >= BOUNDARY_WAIT_HOUR - 1 - 1e-7 for left, _ in gaps):
            return BoundaryInputNotReady(message)
    return FileNotFoundError(message)


def open_hourly_input(path, variable):
    """Retain xarray's reader selection, but identify files that cannot open."""
    try:
        return xr.open_dataset(path)
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        if BOUNDARY_WAIT_HOUR is not None:
            match = re.match(r'(\d{4})[_-]?(\d{2})[_-]?(\d{2})', Path(path).name)
            if match:
                try:
                    file_date = date(*map(int, match.groups()))
                except ValueError:
                    file_date = None
                boundary_date = pd.to_datetime(BOUNDARY_WAIT_HOUR, unit='h').date()
                if file_date is not None and file_date >= boundary_date:
                    raise BoundaryInputNotReady(f'Next-day input is not readable yet: {path}: {exc}') from exc
        try:
            size = f'{Path(path).stat().st_size} bytes'
        except OSError:
            size = 'size unavailable'
        raise RuntimeError(
            f'Cannot open hourly {variable} input: {path} ({size}). '
            'Check that the file is a complete, readable NetCDF file and that '
            'its download/write has finished. No input is skipped or replaced '
            f'with zero. Reader error: {exc}'
        ) from exc


def hourly_file_overlaps_day(path, local_date):
    """Reject unrelated dated files before touching their NetCDF contents.

    A global local-noon day can require UTC dates D-1, D and D+1. Filename
    dates are a coarse discovery filter only; decoded times remain authoritative
    for sampling, intervals, missing-hour and duplicate checks. Date ranges
    (including multi-day/monthly files) are retained when they overlap. Undated
    names remain supported and must be opened because their coverage is unknown.
    """
    if local_date is None:
        return True
    matches = re.findall(r'(?<![0-9])(\d{4})[_-]?(\d{2})[_-]?(\d{2})(?![0-9])', Path(path).name)
    dates = []
    for year, month, day in matches:
        try:
            dates.append(date(int(year), int(month), int(day)))
        except ValueError:
            # An ambiguous name must never silently hide required weather.
            return True
    if not dates:
        return True
    return (min(dates) <= local_date + timedelta(days=1)
            and max(dates) >= local_date - timedelta(days=1))


class HourlySource:
    """Index timestamps first; load only the requested time/longitude slices.

    Filename dates exclude unrelated files before opening; decoded UTC timestamps
    identify the data. Multi-day files must name their full date range. Sources are
    rescanned each calculation day to see newly downloaded data.
    """
    def __init__(self, directories, token, precipitation=False, time_label=PRECIP_TIME_LABEL, local_date=None):
        self.records = {}
        self.intervals = []
        self.lat = self.lon = None
        self.attrs = {}
        self.used_paths = set()
        self._datasets = None
        self.inferred_bounds = False
        self.contiguous_row_storage = True
        pattern = re.compile(rf"(?<![A-Za-z0-9]){re.escape(token)}(?![A-Za-z0-9])")
        paths = sorted({p.resolve() for directory in directories for p in Path(directory).glob('*.nc')
                        if pattern.search(p.name) and hourly_file_overlaps_day(p, local_date) and p.is_file()})
        if not paths:
            raise FileNotFoundError(f"No hourly {token} files in {directories}")
        for path in paths:
            with open_hourly_input(path, token) as ds:
                names = [n for n, da in ds.data_vars.items() if 'time' in da.dims and 'lat' in da.dims and 'lon' in da.dims]
                if len(names) != 1:
                    raise ValueError(f"Expected one hourly weather variable in {path}: {names}")
                da = ds[names[0]]
                self.contiguous_row_storage = (self.contiguous_row_storage
                    and da.dims == ('time', 'lat', 'lon') and da.encoding.get('chunksizes') is None)
                check_grid(da, path)
                if self.lat is None:
                    self.lat, self.lon = da.lat.values.copy(), da.lon.values.copy()
                elif not np.array_equal(self.lat, da.lat.values) or not np.array_equal(self.lon, da.lon.values):
                    raise ValueError(f"Grid changes between {token} files: {path}")
                metadata = copy_selected_global_attrs(da, xr.Dataset())
                copy_selected_global_attrs(ds, metadata)  # Dataset attributes take precedence.
                for key, value in metadata.attrs.items():
                    if key in COMPATIBILITY_KEYS and key in self.attrs and str(self.attrs[key]).lower() != str(value).lower():
                        raise ValueError(f'Conflicting source metadata {key!r}: {path}')
                    self.attrs[key] = value
                times = utc_hours(da.time.values, path)
                if HANDOVER_CONTEXT is not None:
                    check_boundary_source(path, times, metadata.attrs)

                if len(times) == 0:
                    raise ValueError(f"Empty hourly file: {path}")
                if len(times) > 1 and (np.diff(times) <= 0).any():
                    raise ValueError(f"Duplicate/unsorted timestamps: {path}")
                if not precipitation and not np.allclose(times, np.round(times), atol=1e-7, rtol=0):
                    raise ValueError(f"Meteorological timestamps must be on UTC whole hours: {path}")
                bounds = None
                if precipitation:
                    bounds_name = da.time.attrs.get('bounds')
                    if bounds_name:
                        if bounds_name not in ds:
                            raise ValueError(f"Missing time-bounds variable {bounds_name}: {path}")
                        bounds = utc_hours(ds[bounds_name].values, path)
                        if bounds.shape != (len(times), 2):
                            raise ValueError(f"Expected (time,2) precipitation bounds: {path}")
                    elif time_label:
                        self.inferred_bounds = True
                        bounds = np.column_stack((times - 1, times)) if time_label == 'end' else np.column_stack((times, times + 1))
                    else:
                        raise ValueError(f"No precipitation time bounds in {path}. Supply --precip-time-label end or start only after checking hourly interval semantics.")
                    if not np.allclose(bounds[:, 1] - bounds[:, 0], 1.0, atol=1e-7, rtol=0):
                        raise ValueError(f"Precipitation must be independent ONE-hour means/amounts, not daily or forecast-cumulative values: {path}")
                    if not np.isfinite(bounds).all() or ((times < bounds[:, 0]) | (times > bounds[:, 1])).any():
                        raise ValueError(f"Invalid precipitation bounds/timestamps: {path}")
                for i, hour in enumerate(times):
                    record = (path, names[0], i)
                    if precipitation:
                        self.intervals.append((bounds[i, 0], bounds[i, 1], record))
                    else:
                        h = int(round(hour))
                        if h in self.records:
                            previous_path, previous_variable, previous_index = self.records[h]
                            raise ValueError(
                                f"Duplicate {token} UTC hour {pd.to_datetime(h, unit='h')}. "
                                f"First record: {previous_path} variable={previous_variable} time_index={previous_index}; "
                                f"second record: {path} variable={names[0]} time_index={i}. "
                                "These are decoded NetCDF timestamps, not filename dates. "
                                "Check both files and the automatically discovered adjacent input directories; "
                                "no duplicate is silently selected or discarded.")
                        self.records[h] = record
        if precipitation:
            self.intervals.sort(key=lambda x: (x[0], x[1]))
            for prev, curr in zip(self.intervals, self.intervals[1:]):
                if curr[0] < prev[1] - 1e-7:
                    raise ValueError(f"Overlapping/duplicate precipitation intervals: {prev[2][0]} and {curr[2][0]}")

    def read(self, record, columns, kind, rows=None):
        path, name, index = record
        def read_slice(ds):
            if rows is None:
                return ds[name].isel(time=index, lon=columns, drop=True).transpose('lat', 'lon').load().astype('float64')
            # Shared dimension gives paired land points, not a Cartesian product.
            return ds[name].isel(time=index, lat=xr.DataArray(rows, dims='point'),
                                 lon=xr.DataArray(columns, dims='point'), drop=True).transpose('point').load().astype('float64')
        if self._datasets is None:
            with open_hourly_input(path, name) as ds:
                da = read_slice(ds)
        else:
            # Reuse daily/monthly files during this band. Bound handles even if
            # input files each hold just one hour; no dataset crosses processes.
            if path not in self._datasets:
                if len(self._datasets) >= 4:
                    _, oldest = self._datasets.popitem(last=False)
                    oldest.close()
                self._datasets[path] = open_hourly_input(path, name)
            self._datasets.move_to_end(path)
            da = read_slice(self._datasets[path])
        if not np.isfinite(da.values).all():
            raise ValueError(f"Non-finite {kind} input: {path}")
        units = da.attrs.get('units', '').strip().lower()
        if kind == 'temperature':
            da = temperature_kelvin(da)
            if (da.values <= 0).any():
                raise ValueError(f"Invalid Kelvin temperature: {path}")
        elif kind == 'wind':
            if units in ('m/s', 'm s-1', 'm s**-1'):
                pass
            elif units in ('km/h', 'km h-1', 'km h**-1', 'km**-h'):
                da = da / 3.6
            else:
                raise ValueError(f"Unsupported wind units {units!r}: {path}")
        else:  # mm per hour: one-hour mean rates or independent one-hour amounts
            if units in ('kg m-2 s-1', 'kg/m2/s', 'kg m**-2 s**-1', 'mm/s', 'mm s-1'):
                da = da * 3600
            elif units in ('m/s', 'm s-1', 'm s**-1'):
                da = da * 3600 * 1000
            elif units in ('m',):
                da = da * 1000
            elif units in ('mm', 'kg m-2', 'kg m**-2', 'kg/m2', 'mm/h', 'mm h-1', 'mm hr-1'):
                pass
            else:
                raise ValueError(f"Unsupported hourly precipitation units {units!r}: {path}")
            if (da.values < 0).any():
                raise ValueError(f"Negative hourly precipitation: {path}")
        self.used_paths.add(str(path))
        return da.values

    @contextmanager
    def cached_reads(self):
        self._datasets = OrderedDict()
        try:
            yield
        finally:
            datasets, self._datasets = self._datasets, None
            for ds in datasets.values():
                ds.close()

    def sample(self, targets, kind, columns=None):
        # Linear interpolation, never extrapolation; exact-hour targets need one record.
        columns = np.arange(len(self.lon)) if columns is None else np.asarray(columns)
        if np.shape(targets) != (len(columns),):
            raise ValueError('One target hour is required per selected longitude')
        lo = np.floor(targets + 1e-9).astype(np.int64)
        fraction = targets - lo
        fraction[np.abs(fraction) < 1e-7] = 0.0
        result = np.zeros((len(self.lat), len(columns)), dtype=np.float64)
        for hours, weights in ((lo, 1 - fraction), (lo + 1, fraction)):
            for hour in np.unique(hours[weights > 0]):
                positions = np.where((hours == hour) & (weights > 0))[0]
                if int(hour) not in self.records:
                    raise FileNotFoundError(f"Missing required UTC weather hour {pd.to_datetime(hour, unit='h')}")
                result[:, positions] += self.read(self.records[int(hour)], columns[positions], kind) * weights[positions]
        return result

    def accumulate(self, targets, columns=None):
        columns = np.arange(len(self.lon)) if columns is None else np.asarray(columns)
        if np.shape(targets) != (len(columns),):
            raise ValueError('One target hour is required per selected longitude')
        starts = targets - 24.0
        coverage = np.zeros(len(columns))
        result = np.zeros((len(self.lat), len(columns)), dtype=np.float64)
        for begin, end, record in self.intervals:
            overlap = np.maximum(0.0, np.minimum(end, targets) - np.maximum(begin, starts))
            positions = np.where(overlap > 1e-8)[0]
            if not len(positions):
                continue
            result[:, positions] += self.read(record, columns[positions], 'rain') * overlap[positions]
            coverage += overlap
        if not np.allclose(coverage, 24, rtol=0, atol=1e-6):
            bad = np.where(np.abs(coverage - 24) > 1e-6)[0]
            raise FileNotFoundError(f"Incomplete 24-hour rain window at {len(bad)} longitudes; first lon={self.lon[columns[bad[0]]]}, coverage={coverage[bad[0]]:.6f} hours")
        return result


    def sample_points(self, targets, rows, columns, kind):
        lo = np.floor(targets + 1e-9).astype(np.int64)
        fraction = targets - lo
        fraction[np.abs(fraction) < 1e-7] = 0.0
        result = np.zeros(len(rows), dtype=np.float64)
        for hours, weights in ((lo, 1 - fraction), (lo + 1, fraction)):
            for hour in np.unique(hours[weights > 0]):
                positions = np.flatnonzero((hours == hour) & (weights > 0))
                if int(hour) not in self.records:
                    raise missing_weather_hour(hour)
                result[positions] += self.read(self.records[int(hour)], columns[positions], kind,
                                               rows=rows[positions]) * weights[positions]
        return result

    def accumulate_points(self, targets, rows, columns):
        starts = targets - 24.0
        coverage = np.zeros(len(rows))
        result = np.zeros(len(rows), dtype=np.float64)
        if not len(rows):
            return result
        earliest, latest = starts.min(), targets.max()
        for begin, end, record in self.intervals:
            if end <= earliest or begin >= latest:
                continue
            overlap = np.maximum(0.0, np.minimum(end, targets) - np.maximum(begin, starts))
            positions = np.flatnonzero(overlap > 1e-8)
            if not len(positions):
                continue
            result[positions] += self.read(record, columns[positions], 'rain', rows=rows[positions]) * overlap[positions]
            coverage += overlap
        if not np.allclose(coverage, 24, rtol=0, atol=1e-6):
            bad = np.flatnonzero(np.abs(coverage - 24) > 1e-6)[0]
            raise FileNotFoundError(f'Incomplete 24-hour rain window at land point '
                                    f'lat={self.lat[rows[bad]]}, lon={self.lon[columns[bad]]}; '
                                    f'coverage={coverage[bad]:.6f} hours')
        return result



def hourly_directories(input_dir, stamp, extra_dirs=(), monthly=False):
    """Use only explicitly supplied integration directories.

    Autosubmit's YYYYMMDD directory identifies the integration start, not the
    current chunk. It stays fixed across months and years. Never infer sibling
    directories from the calculation date: they may belong to other runs.
    stamp/monthly remain accepted for compatibility with existing callers.
    Filename dates and decoded timestamps select the required weather inside
    these directories. Historical handover directories are supplied explicitly.
    """
    candidates = [Path(input_dir)] + [Path(directory) for directory in extra_dirs]
    return list(dict.fromkeys(path.resolve() for path in candidates if path.is_dir()))


def active_land(mask, threshold=LAND_THRESHOLD):
    if not np.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError('Land threshold must be between 0 and 1')
    values = np.asarray(mask)
    return values > 0 if threshold == 0 else values >= threshold


def normalize_land_mask(mask, lat, lon):
    """Align a static fractional mask by coordinate identity; never interpolate."""
    mask = mask.copy(deep=False)
    units = str(mask.attrs.get('units', '')).strip().lower()
    if units not in ('', '1', '1.0', 'fraction', 'dimensionless', '(0 - 1)', '0-1', '0 - 1', '~'):
        raise ValueError(f'lsm must use fractional units, not {units!r}; convert the mask explicitly first')
    for dim in list(mask.dims):
        if dim not in ('lat', 'lon'):
            if mask.sizes[dim] != 1:
                raise ValueError(f'lsm must be static: extra dimension {dim!r} is not singleton')
            mask = mask.isel({dim: 0}, drop=True)
    if set(mask.dims) != {'lat', 'lon'}:
        raise ValueError('lsm must have lat and lon dimensions')
    for name in ('lat', 'lon'):
        if name not in mask.coords or mask[name].dims != (name,):
            raise ValueError(f'Explicit one-dimensional {name} coordinates required for lsm')
    # Reuse the weather coordinate checks, including duplicate cyclic longitudes.
    check_grid(mask.expand_dims(time=[pd.Timestamp('2000-01-01')]), 'lsm')
    target_lat, target_lon = np.asarray(lat), np.asarray(lon)
    normalized_lon = (mask.lon.values.astype(float) + 180) % 360 - 180
    target_normalized = (target_lon.astype(float) + 180) % 360 - 180
    tolerance = 1e-10  # Only roundoff in equivalent longitude conventions; no spatial interpolation.
    normalized_lon[np.isclose(np.abs(normalized_lon), 180, rtol=0, atol=tolerance)] = -180
    target_normalized[np.isclose(np.abs(target_normalized), 180, rtol=0, atol=tolerance)] = -180

    def coordinate_indices(values, target):
        if len(values) != len(target):
            raise ValueError('lsm grid does not match the weather grid; regrid the mask explicitly first')
        order = np.argsort(values)
        sorted_values = np.asarray(values)[order]
        if len(values) > 1 and (np.diff(sorted_values) <= 2 * tolerance).any():
            raise ValueError('Ambiguous lsm coordinates within rounding tolerance')
        indices = pd.Index(sorted_values).get_indexer(target, method='nearest', tolerance=tolerance)
        if (indices < 0).any() or len(np.unique(indices)) != len(target):
            raise ValueError('lsm grid does not match the weather grid; regrid the mask explicitly first')
        return order[indices]

    lat_index = coordinate_indices(mask.lat.values, target_lat)
    lon_index = coordinate_indices(normalized_lon, target_normalized)
    aligned = mask.isel(lat=lat_index, lon=lon_index).transpose('lat', 'lon').load()
    values = aligned.values.astype(np.float64)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError('lsm must contain finite land fractions between 0 and 1')
    # Keep fractions at float64 precision through saved state and annual reuse.
    return xr.DataArray(values, dims=('lat', 'lon'), coords={'lat': target_lat, 'lon': target_lon},
                        name='lsm', attrs={'long_name': 'Land-sea mask (land fraction)', 'units': '1'})


def discover_land_mask(directories, lat, lon, lsm_file=None, fallback=None, local_date=None):
    """Use the integration's saved static mask; discover inputs only initially."""
    if fallback is not None:
        saved = normalize_land_mask(fallback, lat, lon)
        print('Reusing land-sea mask from saved FWI state; input LSM discovery skipped.', flush=True)
        return saved, []
    if lsm_file is not None:
        paths = [Path(lsm_file).expanduser().resolve()]
        if not paths[0].is_file():
            raise FileNotFoundError(f'Land-sea mask file missing: {paths[0]}')
    else:
        def generated(path):
            return (path.name.startswith('FWI_Const') or path.name.endswith(
                ('_LNOON_inputs.nc', '_fwi.nc', '_fwi_comp.nc', '_fwi_restart_file.nc',
                 '_fwi_basic_and_days.nc', '_fwi_quantiles.nc', '_fwi_comp_annual.nc', '_fwi_OverWintering.nc')))
        weather_pattern = re.compile(r'(?<![A-Za-z0-9])(?:' + '|'.join(map(re.escape, VARIABLES.values())) + r')(?![A-Za-z0-9])')
        paths = sorted({path.resolve() for directory in directories for path in Path(directory).glob('*.nc')
                        if not generated(path)
                        and (not weather_pattern.search(path.name) or hourly_file_overlaps_day(path, local_date))
                        and path.is_file()})
    selected = None
    mask_paths = []
    for path in paths:
        with xr.open_dataset(path) as ds:
            if 'lsm' not in ds:
                if lsm_file is not None:
                    raise ValueError(f'No variable named lsm in {path}')
                continue
            candidate = normalize_land_mask(ds.lsm, lat, lon)
        if selected is not None and not np.array_equal(selected.values, candidate.values):
            raise ValueError(f'Conflicting lsm masks in input directories: {path}')
        selected = candidate
        mask_paths.append(str(path))
    if selected is None:
        print('WARNING: Land-sea mask (lsm) not found and no cached mask is available. '
              'Land-sea masking is NOT applied; FWI will be calculated on all grid cells, '
              'including ocean cells.', file=sys.stderr, flush=True)
    return selected, mask_paths


def temperature_maximum_24h(source, targets, columns, rows=None):
    """Maximum of piecewise-linear hourly temperature over [noon-24h, noon].

    Include UTC-hour knots and interpolated endpoints. Sampling only at 24
    equally offset times could miss a peak between those times.
    """
    targets = np.asarray(targets, dtype=float)
    columns = np.asarray(columns, dtype=int)
    if not len(targets):
        return np.empty((len(source.lat), 0)) if rows is None else np.empty(0)
    starts = targets - 24.0
    sample = (lambda hours: source.sample(hours, 'temperature', columns)) if rows is None else (
        lambda hours: source.sample_points(hours, rows, columns, 'temperature'))
    values = np.maximum(sample(starts), sample(targets))
    for hour in range(int(np.floor(starts.min())) + 1, int(np.ceil(targets.max()))):
        positions = np.flatnonzero((hour > starts + 1e-7) & (hour < targets - 1e-7))
        if not len(positions):
            continue
        if hour not in source.records:
            raise FileNotFoundError(f'Missing required UTC temperature hour for overwintering maximum: {pd.to_datetime(hour, unit="h")}')
        if rows is None:
            sample_values = source.read(source.records[hour], columns[positions], 'temperature')
            values[:, positions] = np.maximum(values[:, positions], sample_values)
        else:
            sample_values = source.read(source.records[hour], columns[positions], 'temperature', rows=rows[positions])
            values[positions] = np.maximum(values[positions], sample_values)
    return values


def weather_output_keys(sources):
    return list(sources) + (['tasmax'] if getattr(sources['tas'], 'with_maximum_temperature', False) else [])


def prepare_temperature_maxima(source, targets):
    values = np.full((len(source.lat), len(source.lon)), np.nan)
    hours = np.floor(targets + 1e-9).astype(np.int64)
    for hour in np.unique(hours):
        columns = np.flatnonzero(hours == hour)
        with source.cached_reads():
            values[:, columns] = temperature_maximum_24h(source, targets[columns], columns)
    return values


def prepare_masked_band(sources, columns, targets, land):
    """Read and calculate paired land points, then restore this band's shape."""
    rows, local_columns = np.nonzero(land[:, columns])
    values = {key: np.full((land.shape[0], len(columns)), np.nan) for key in weather_output_keys(sources)}
    if not len(rows):
        return values
    point_columns = columns[local_columns]
    point_targets = targets[local_columns]
    for key, source in sources.items():
        with source.cached_reads():
            points = (source.accumulate_points(point_targets, rows, point_columns) if key == 'pr' else
                      source.sample_points(point_targets, rows, point_columns,
                                           'temperature' if key in ('tas', 'dew') else 'wind'))
        values[key][rows, local_columns] = points
        if key == 'tas' and getattr(source, 'with_maximum_temperature', False):
            with source.cached_reads():
                maxima = temperature_maximum_24h(source, point_targets, point_columns, rows)
            values['tasmax'][rows, local_columns] = maxima
    return values


def prepare_masked_values(sources, targets, land):
    """Serial band preparation keeps temporary point arrays bounded in size."""
    values = {key: np.full(land.shape, np.nan) for key in weather_output_keys(sources)}
    active_columns = np.flatnonzero(land.any(axis=0))
    hours = np.floor(targets + 1e-9).astype(np.int64)
    for hour in np.unique(hours[active_columns]):
        columns = active_columns[hours[active_columns] == hour]
        band = prepare_masked_band(sources, columns, targets[columns], land)
        for key in weather_output_keys(sources):
            values[key][:, columns] = band[key]
    return values


def add_derived_weather(fields, land=None):
    if land is None:
        fields['relative_humidity'] = relative_humidity(fields.tas, fields.dewpoint)
        fields['wind_speed'] = np.hypot(fields.u10, fields.v10) * 3.6
    else:
        rh = np.full(land.shape, np.nan)
        wind = np.full(land.shape, np.nan)
        if land.any():
            tas = xr.DataArray(fields.tas.values[land], dims='point', attrs={'units': 'K'})
            dew = xr.DataArray(fields.dewpoint.values[land], dims='point', attrs={'units': 'K'})
            rh[land] = relative_humidity(tas, dew).values
            wind[land] = np.hypot(fields.u10.values[land], fields.v10.values[land]) * 3.6
        fields['relative_humidity'] = fields.tas.copy(data=rh)
        fields['wind_speed'] = fields.tas.copy(data=wind)
    fields['relative_humidity'].attrs = {'units': '%'}
    fields['wind_speed'].attrs = {'units': 'km h-1'}



def validate_io_settings(tile_columns, buffer_mb):
    if isinstance(tile_columns, bool) or not isinstance(tile_columns, (int, np.integer)) or tile_columns < 1:
        raise ValueError('io_tile_columns must be a positive integer')
    if not np.isfinite(buffer_mb) or buffer_mb <= 0:
        raise ValueError('io_buffer_mb must be finite and positive')
    return int(buffer_mb * 1024 * 1024)


def read_tile_records(source, records, first_column, stop_column, buffer_bytes,
                      first_row=0, stop_row=None):
    """Yield decoded hours from bounded rectangular reads, with handles reused.

    Input records are already in scientific processing order. Only adjacent
    indices within the SAME file are combined. Never load an entire month or
    use paired/fancy indexing against the NetCDF backend. The caller releases
    each yielded view before requesting the next block.
    """
    stop_row = len(source.lat) if stop_row is None else stop_row
    position = 0
    while position < len(records):
        path, name, index = records[position]
        if path not in source._datasets:
            if len(source._datasets) >= 4:
                _, oldest = source._datasets.popitem(last=False)
                oldest.close()
            source._datasets[path] = open_hourly_input(path, name)
        source._datasets.move_to_end(path)
        variable = source._datasets[path][name]
        bytes_per_hour = ((stop_row - first_row) * (stop_column - first_column)
                          * max(8, variable.dtype.itemsize))
        count_limit = buffer_bytes // bytes_per_hour
        if count_limit < 1:
            raise ValueError('I/O buffer cannot hold one hourly tile; reduce io_tile_columns or increase io_buffer_mb')
        end = position + 1
        while end < len(records) and end - position < count_limit:
            next_path, next_name, next_index = records[end]
            if next_path != path or next_name != name or next_index != index + end - position:
                break
            end += 1
        block = variable.isel(time=slice(index, index + end - position),
                              lat=slice(first_row, stop_row),
                              lon=slice(first_column, stop_column)).transpose('time', 'lat', 'lon').load()
        source.block_reads += 1
        attrs = dict(block.attrs)
        for offset in range(end - position):
            yield records[position + offset], block.values[offset], attrs
        del block
        position = end


def convert_tile_values(source, record, raw, rows, columns, attrs, kind):
    """Match HourlySource.read arithmetic, checking only contributing cells."""
    values = xr.DataArray(raw[rows, columns], dims='point', attrs=attrs).astype('float64')
    path = record[0]
    if not np.isfinite(values.values).all():
        raise ValueError(f'Non-finite {kind} input: {path}')
    units = values.attrs.get('units', '').strip().lower()
    if kind == 'temperature':
        values = temperature_kelvin(values)
        if (values.values <= 0).any():
            raise ValueError(f'Invalid Kelvin temperature: {path}')
    elif kind == 'wind':
        if units in ('m/s', 'm s-1', 'm s**-1'):
            pass
        elif units in ('km/h', 'km h-1', 'km h**-1', 'km**-h'):
            values = values / 3.6
        else:
            raise ValueError(f'Unsupported wind units {units!r}: {path}')
    else:
        if units in ('kg m-2 s-1', 'kg/m2/s', 'kg m**-2 s**-1', 'mm/s', 'mm s-1'):
            values = values * 3600
        elif units in ('m/s', 'm s-1', 'm s**-1'):
            values = values * 3600 * 1000
        elif units in ('m',):
            values = values * 1000
        elif units in ('mm', 'kg m-2', 'kg m**-2', 'kg/m2', 'mm/h', 'mm h-1', 'mm hr-1'):
            pass
        else:
            raise ValueError(f'Unsupported hourly precipitation units {units!r}: {path}')
        if (values.values < 0).any():
            raise ValueError(f'Negative hourly precipitation: {path}')
    source.used_paths.add(str(path))
    return values.values


def interpolation_hours(targets):
    lower = np.floor(targets + 1e-9).astype(np.int64)
    fraction = targets - lower
    fraction[np.abs(fraction) < 1e-7] = 0.0
    return lower, fraction


def prepare_weather_tile(source, targets, first, stop, land, kind, buffer_bytes, maximum=False,
                         first_row=0, stop_row=None):
    lower, fraction = interpolation_hours(targets)
    active_columns = land.any(axis=0)
    needed = set(lower[active_columns]) | set((lower + 1)[active_columns & (fraction > 0)])
    starts = targets - 24.0
    if maximum:
        start_lower, start_fraction = interpolation_hours(starts)
        needed.update(start_lower[active_columns])
        needed.update((start_lower + 1)[active_columns & (start_fraction > 0)])
        for hour in range(int(np.floor(starts[active_columns].min())) + 1,
                          int(np.ceil(targets[active_columns].max()))):
            if np.any(active_columns & (hour > starts + 1e-7) & (hour < targets - 1e-7)):
                needed.add(hour)
    hours = sorted(needed)
    for hour in hours:
        if hour not in source.records:
            raise missing_weather_hour(hour)
    records = [source.records[hour] for hour in hours]
    result = np.where(land, 0.0, np.nan)
    if maximum:
        start_values = np.where(land, 0.0, np.nan)
        interior_maximum = np.full(land.shape, -np.inf)
    position = 0
    for record, raw, attrs in read_tile_records(source, records, first, stop, buffer_bytes, first_row, stop_row):
        hour = hours[position]
        weights = np.where(lower == hour, 1 - fraction, 0.0) + np.where(lower + 1 == hour, fraction, 0.0)
        selected_columns = weights > 0
        if maximum:
            start_weights = (np.where(start_lower == hour, 1 - start_fraction, 0.0)
                             + np.where(start_lower + 1 == hour, start_fraction, 0.0))
            interior = (hour > starts + 1e-7) & (hour < targets - 1e-7)
            selected_columns = selected_columns | (start_weights > 0) | interior
        rows, columns = np.nonzero(land & selected_columns[None, :])
        values = convert_tile_values(source, record, raw, rows, columns, attrs, kind)
        used = weights[columns] > 0
        result[rows[used], columns[used]] += values[used] * weights[columns[used]]
        if maximum:
            used = start_weights[columns] > 0
            start_values[rows[used], columns[used]] += values[used] * start_weights[columns[used]]
            used = interior[columns]
            r, c = rows[used], columns[used]
            interior_maximum[r, c] = np.maximum(interior_maximum[r, c], values[used])
        del raw
        position += 1
    maxima = np.maximum(np.maximum(start_values, result), interior_maximum) if maximum else None
    return result, maxima


def prepare_rain_tile(source, targets, first, stop, land, buffer_bytes, first_row=0, stop_row=None):
    starts = targets - 24.0
    active_columns = land.any(axis=0)
    plan = []
    coverage = np.zeros(len(targets))
    for begin, end, record in source.intervals:
        overlap = np.maximum(0.0, np.minimum(end, targets) - np.maximum(begin, starts))
        if np.any((overlap > 1e-8) & active_columns):
            plan.append((record, overlap))
            coverage += overlap
    if not np.allclose(coverage[active_columns], 24, rtol=0, atol=1e-6):
        bad = np.flatnonzero(active_columns & (np.abs(coverage - 24) > 1e-6))[0]
        raise missing_rain_window(source, starts[bad], targets[bad],
                                  f'Incomplete 24-hour rain window at lon={source.lon[first + bad]}; '
                                  f'coverage={coverage[bad]:.6f} hours')
    result = np.where(land, 0.0, np.nan)
    records = [record for record, _ in plan]
    position = 0
    for record, raw, attrs in read_tile_records(source, records, first, stop, buffer_bytes, first_row, stop_row):
        weights = plan[position][1]
        rows, columns = np.nonzero(land & (weights > 1e-8)[None, :])
        values = convert_tile_values(source, record, raw, rows, columns, attrs, 'rain')
        # Preserve the original chronological additions: no cumsum subtraction
        # or reordered reductions, including at fractional-hour endpoints.
        result[rows, columns] += values * weights[columns]
        del raw
        position += 1
    return result


def prepare_tiled_values(sources, targets, land, tile_columns, buffer_mb, tile_rows):
    buffer_bytes = validate_io_settings(tile_columns, buffer_mb)
    validate_io_settings(tile_rows, buffer_mb)
    source = sources['tas']
    shape = (len(source.lat), len(source.lon))
    tile_columns = min(tile_columns, buffer_bytes // (shape[0] * 8))
    if tile_columns < 1:
        raise ValueError('I/O buffer cannot hold one full-latitude longitude column')
    values = {}
    for key, source in sources.items():
        started = perf_counter()
        maximum = key == 'tas' and source.with_maximum_temperature
        values[key] = np.full(shape, np.nan)
        if maximum:
            values['tasmax'] = np.full(shape, np.nan)
        source.block_reads = 0
        # A contiguous file stores a longitude row as one byte span. Narrow
        # full-latitude bands cause thousands of strided reads per hour. For
        # 24-hour reductions read full-width latitude strips instead. This can
        # read extra non-contributing hours, but every used sample is unchanged.
        # Chunked/compressed inputs retain longitude tiles to limit decoding of
        # unrelated chunks. Weather needing only noon also retains narrow tiles.
        row_reads = ((key == 'pr' or maximum) and source.contiguous_row_storage
                     and buffer_bytes >= shape[1] * 8)
        rows_per_tile = min(tile_rows, buffer_bytes // (shape[1] * 8)) if row_reads else shape[0]
        columns_per_tile = shape[1] if row_reads else tile_columns
        tiles = ((row, min(row + rows_per_tile, shape[0]), first, min(first + columns_per_tile, shape[1]))
                 for row in range(0, shape[0], rows_per_tile)
                 for first in range(0, shape[1], columns_per_tile))
        # Handles survive across every tile of this source, then close on both
        # success and exceptions. A new day still rebuilds the timestamp index.
        with source.cached_reads():
            for first_row, stop_row, first, stop in tiles:
                active = (np.ones((stop_row - first_row, stop - first), dtype=bool) if land is None
                          else land[first_row:stop_row, first:stop])
                if not active.any():
                    continue
                if key == 'pr':
                    values[key][first_row:stop_row, first:stop] = prepare_rain_tile(
                        source, targets[first:stop], first, stop, active, buffer_bytes, first_row, stop_row)
                else:
                    result, maxima = prepare_weather_tile(source, targets[first:stop], first, stop, active,
                        'temperature' if key in ('tas', 'dew') else 'wind', buffer_bytes, maximum, first_row, stop_row)
                    values[key][first_row:stop_row, first:stop] = result
                    if maximum:
                        values['tasmax'][first_row:stop_row, first:stop] = maxima
        report_timing(f'hourly {key}' + (' including 24h temperature maximum' if maximum else ''),
                      started, f' ({source.block_reads} block reads, {"latitude strips" if row_reads else "longitude tiles"})')
    return values


def prepare_local_inputs(stamp, directories, sampling='linear', precip_time_label=PRECIP_TIME_LABEL,
                         lsm_file=None, land_threshold=LAND_THRESHOLD, land_mask=None, overwintering=DO_OVERWINTERING,
                         io_tile_columns=IO_TILE_COLUMNS, io_buffer_mb=IO_BUFFER_MB, io_tile_rows=IO_TILE_ROWS):
    started = perf_counter()
    validate_io_settings(io_tile_columns, io_buffer_mb)
    validate_io_settings(io_tile_rows, io_buffer_mb)
    sources = {key: HourlySource(directories, token, key == 'pr', precip_time_label, local_date=stamp)
               for key, token in VARIABLES.items()}
    first = sources['tas']
    first.with_maximum_temperature = bool(overwintering)
    for source in sources.values():
        if not np.array_equal(source.lat, first.lat) or not np.array_equal(source.lon, first.lon):
            raise ValueError('All hourly variables must have identical lat/lon coordinates and ordering')
    active_land(np.zeros(1), land_threshold)  # Validate the configured cutoff even without a mask.
    mask, mask_paths = discover_land_mask(directories, first.lat, first.lon, lsm_file, land_mask, local_date=stamp)
    report_timing('input indexing and mask discovery', started)
    land = None if mask is None else active_land(mask, land_threshold)
    targets = local_noon_hours(stamp, first.lon, sampling)
    coords = {'lat': first.lat, 'lon': first.lon}
    fields = xr.Dataset()
    output_names = {'tas':'tas', 'dew':'dewpoint', 'u':'u10', 'v':'v10', 'pr':'precip_24h'}
    prepared_values = prepare_tiled_values(sources, targets, land, io_tile_columns, io_buffer_mb, io_tile_rows)
    for key, source in sources.items():
        values = prepared_values.pop(key)
        da = xr.DataArray(values, dims=('lat','lon'), coords=coords)
        da.attrs['units'] = 'mm' if key == 'pr' else ('K' if key in ('tas','dew') else 'm s-1')
        fields[output_names[key]] = da
    if overwintering:
        maxima = prepared_values.pop('tasmax')
        fields['tasmax_24h'] = xr.DataArray(maxima, dims=('lat', 'lon'), coords=coords,
            attrs={'units': 'K', 'long_name': 'Maximum hourly temperature over previous 24 hours ending at local noon'})
    add_derived_weather(fields, land)
    for name, hours in (('local_noon_utc', targets), ('precip_window_start_utc', targets - 24)):
        values = np.rint(hours * 3.6e12).astype(np.int64).astype('datetime64[ns]')
        fields[name] = xr.DataArray(values, dims=('lon',), coords={'lon': first.lon})
        fields[name].attrs['long_name'] = name.replace('_', ' ')
    # Original ds_t2m attribute list. OPA can omit Dataset metadata on its
    # aggregated temperature file: use variable metadata or matching companion
    # sources to fill omissions, never hardcode an experiment/model identity.
    for source in sources.values():
        for key, value in source.attrs.items():
            if key in COMPATIBILITY_KEYS and key in fields.attrs and str(fields.attrs[key]).lower() != str(value).lower():
                raise ValueError(f'Hourly variables have conflicting source metadata {key!r}')
            fields.attrs[key] = value
    if HANDOVER_CONTEXT is not None:
        branch = 'historical' if stamp < HANDOVER_CONTEXT['start'] else 'projection'
        fields.attrs.update(HANDOVER_CONTEXT['attrs'][branch])
    fields.attrs.update(noon_convention=NOON_CONVENTION, sampling=sampling,
                        local_date=stamp.isoformat(),
                        precipitation_intervals='inferred one hour, timestamp '+str(precip_time_label) if sources['pr'].inferred_bounds else 'validated one-hour time bounds',
                        precipitation_partial_hours='uniform rate within each hourly interval')
    fields.attrs['input_files'] = '; '.join(sorted({p for source in sources.values() for p in source.used_paths}))
    if mask is not None:
        fields['lsm'] = mask
        fields.attrs.update(land_threshold=float(land_threshold), land_mask_files='; '.join(mask_paths))
        print(f'Land-sea mask: {int(land.sum())} of {land.size} grid cells selected.', flush=True)
    return fields


def existing_daily_product(output, stamp):
    """Read earlier local-noon filenames while writing only the requested plain names."""
    paths = [Path(output) / f'{stamp:%Y_%m_%d}_fwi.nc',
             Path(output) / f'{stamp:%Y_%m_%d}_LNOON_fwi.nc']
    found = [path for path in paths if path.exists()]
    if len(found) > 1:
        raise ValueError(f'Ambiguous daily FWI products for {stamp}: {found}')
    return found[0] if found else paths[0]


def write_prepared_inputs(stamp, fields, input_dir):
    """Prepared daily weather belongs alongside the hourly OPA inputs."""
    prepared = fields.expand_dims(time=[pd.Timestamp(stamp) + pd.Timedelta(hours=12)])
    if 'lsm' in fields:
        prepared['lsm'] = fields.lsm  # Static mask, without a repeated time dimension.
    prepared.time.attrs['comment'] = 'Nominal local-date label, not a common UTC valid time; see local_noon_utc(lon).'
    atomic_write(prepared, Path(input_dir) / f'{stamp:%Y_%m_%d}_LNOON_inputs.nc', engine='netcdf4',
                 encoding={k:{'zlib':True,'complevel':4} for k,v in prepared.data_vars.items() if np.issubdtype(v.dtype,np.floating)})


def validate_overwintering_settings():
    """Validate configurable DC-only overwintering parameters before a run."""
    parameters = {
        'overwintering_carry_over': OVERWINTER_CARRY_OVER,
        'overwintering_wetting_efficiency': OVERWINTER_WETTING_EFFICIENCY,
        'overwintering_start_temp': OVERWINTER_START_TEMP,
        'overwintering_end_temp': OVERWINTER_END_TEMP,
    }
    for name, value in parameters.items():
        if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
            raise ValueError(f'{name} must be a finite number')
        try:
            finite = np.isfinite(float(value))
        except (TypeError, ValueError):
            finite = False
        if not finite:
            raise ValueError(f'{name} must be a finite number')
        parameters[name] = float(value)
    if not 0 < parameters['overwintering_carry_over'] <= 1:
        raise ValueError('OVERWINTER_CARRY_OVER must be greater than 0 and at most 1')
    if not 0 <= parameters['overwintering_wetting_efficiency'] <= 1:
        raise ValueError('OVERWINTER_WETTING_EFFICIENCY must be between 0 and 1')
    if parameters['overwintering_start_temp'] <= parameters['overwintering_end_temp']:
        raise ValueError('OVERWINTER_START_TEMP must exceed OVERWINTER_END_TEMP')
    consecutive = OVERWINTER_CONSECUTIVE_DAYS
    if isinstance(consecutive, (bool, np.bool_)) or not isinstance(consecutive, (int, np.integer)) or consecutive < 1:
        raise ValueError('OVERWINTER_CONSECUTIVE_DAYS must be a positive integer')
    parameters['overwintering_consecutive_days'] = int(consecutive)
    parameters['overwintering_method'] = 'trailing_noon_to_noon_tasmax_v1'
    return parameters


def advance_overwintering(stamp, fields, previous_ow_state, initial_dc_previous,
                         regular_dc, dmc, isi, cfg, land=None):
    """Advance an independent DC stream; regular FFMC/DMC/ISI remain unchanged.

    Season transitions use the last N completed local-noon-to-noon temperature
    maxima, including today, with strict > start and < end thresholds. Each
    cell is initialized from ordinary DC until its first season is detected.
    The first shutdown then establishes a complete winter history. Shutdown
    precipitation belongs to winter; spring startup precipitation belongs to
    the ordinary startup-day DC update. No interval is omitted or counted twice.
    Winter DCPrev is held finite for restart, while winter FWI output is NaN.
    """
    settings = validate_overwintering_settings()
    stamp = pd.Timestamp(stamp).date()
    if 'lat' not in fields.coords or 'lon' not in fields.coords:
        raise ValueError('Overwintering fields require latitude and longitude coordinates')
    latitude = np.asarray(fields.lat.values)
    longitude = np.asarray(fields.lon.values)
    if latitude.ndim != 1 or longitude.ndim != 1:
        raise ValueError('Overwintering requires one-dimensional latitude and longitude coordinates')
    shape = (latitude.size, longitude.size)
    active_land = np.ones(shape, dtype=bool) if land is None else np.asarray(land, dtype=bool)
    if active_land.shape != shape:
        raise ValueError('Overwintering land mask grid differs from the weather grid')

    def packed(value, name, *, finite=True, nonnegative=False):
        if isinstance(value, xr.DataArray):
            if set(value.dims) != {'lat', 'lon'}:
                raise ValueError(f'{name} must have only latitude and longitude dimensions')
            if not np.array_equal(value.lat.values, latitude) or not np.array_equal(value.lon.values, longitude):
                raise ValueError(f'{name} grid differs from the weather grid')
            value = value.transpose('lat', 'lon').values
        values = np.asarray(value, dtype=np.float64)
        if values.shape != shape:
            raise ValueError(f'{name} shape differs from the weather grid')
        values = values[active_land].copy()
        if finite and not np.all(np.isfinite(values)):
            raise ValueError(f'{name} contains non-finite values on land')
        if nonnegative and np.any(values < 0):
            raise ValueError(f'{name} contains negative values on land')
        return values

    for name in ('tas', 'tasmax_24h', 'precip_24h'):
        if name not in fields:
            raise ValueError(f'Missing {name} needed for overwintering')
    temperature = packed(fields.tas, 'tas')
    maximum = packed(fields.tasmax_24h, 'tasmax_24h')
    if np.any(temperature <= 0) or np.any(maximum <= 0):
        raise ValueError('Overwintering temperatures must be positive Kelvin values')
    temperature -= 273.15
    maximum -= 273.15
    rain = packed(fields.precip_24h, 'precip_24h', nonnegative=True)
    baseline_dc = packed(regular_dc, 'regular DC', nonnegative=True)
    dmc_values = packed(dmc, 'DMC', nonnegative=True)
    isi_values = packed(isi, 'ISI', nonnegative=True)
    n = baseline_dc.size
    if previous_ow_state is None:
        own_dc_previous = packed(initial_dc_previous, 'initial previous DC', nonnegative=True)
        season_active = np.zeros(n, dtype=bool)
        season_seen = np.zeros(n, dtype=bool)
        fall_dc = np.full(n, np.nan, dtype=np.float64)
        winter_rain = np.zeros(n, dtype=np.float64)
        warm_days = np.zeros(n, dtype=np.int32)
        cold_days = np.zeros(n, dtype=np.int32)
    else:
        expected = (stamp - timedelta(days=1)).isoformat()
        if previous_ow_state.attrs.get('state_date') != expected:
            raise ValueError(f'Overwintering state must represent {expected}, immediately before {stamp}')
        for name, expected_value in settings.items():
            if previous_ow_state.attrs.get(name) != expected_value:
                raise ValueError(f'Overwintering state setting {name} differs; cannot change method or parameters mid-run')
        names = ('DCPrev', 'SeasonActive', 'DCf', 'rw', 'season_seen', 'warm_days', 'cold_days')
        for name in names:
            if name not in previous_ow_state:
                raise ValueError(f'Overwintering state lacks {name}')
        own_dc_previous = packed(previous_ow_state.DCPrev, 'overwintering DCPrev', nonnegative=True)
        fall_dc = packed(previous_ow_state.DCf, 'overwintering DCf', finite=False)
        winter_rain = packed(previous_ow_state.rw, 'overwintering rw', nonnegative=True)
        flags = {}
        for name in ('SeasonActive', 'season_seen'):
            values = packed(previous_ow_state[name], name)
            if np.any((values != 0) & (values != 1)):
                raise ValueError(f'Overwintering state {name} must contain only 0 or 1 on land')
            flags[name] = values.astype(bool)
        season_active = flags['SeasonActive']
        season_seen = flags['season_seen']
        counters = {}
        for name in ('warm_days', 'cold_days'):
            values = packed(previous_ow_state[name], name, nonnegative=True)
            if np.any(values != np.floor(values)) or np.any(values > OVERWINTER_CONSECUTIVE_DAYS):
                raise ValueError(f'Overwintering state {name} contains invalid consecutive-day counters')
            counters[name] = values.astype(np.int32)
        warm_days = counters['warm_days']
        cold_days = counters['cold_days']
        winter = season_seen & ~season_active
        if np.any(season_active & ~season_seen):
            raise ValueError('Overwintering active-season flag is inconsistent with seasonal history')
        if not np.all(np.isfinite(fall_dc[winter])) or np.any(fall_dc[winter] < 0):
            raise ValueError('Overwintering winter state requires a finite, nonnegative autumn DC')
        if np.any(~np.isnan(fall_dc[~winter])) or np.any(winter_rain[~winter] != 0):
            raise ValueError('Overwintering autumn DC and winter rain are inconsistent with seasonal state')

    warm_days = np.where(maximum > OVERWINTER_START_TEMP,
                         np.minimum(warm_days + 1, OVERWINTER_CONSECUTIVE_DAYS), 0)
    cold_days = np.where(maximum < OVERWINTER_END_TEMP,
                         np.minimum(cold_days + 1, OVERWINTER_CONSECUTIVE_DAYS), 0)
    starting = ~season_active & (warm_days >= OVERWINTER_CONSECUTIVE_DAYS)
    ending = season_active & (cold_days >= OVERWINTER_CONSECUTIVE_DAYS)
    restarting = starting & season_seen
    continuing = season_active & ~ending
    # Prior to the first detected active season, retain the ordinary DC exactly.
    own_dc = own_dc_previous.copy()
    own_dc[~season_seen] = baseline_dc[~season_seen]
    lat_grid = np.broadcast_to(latitude[:, None], shape)[active_land].astype(np.float64)
    if np.any(continuing):
        own_dc[continuing] = DC(temperature[continuing], rain[continuing],
                               own_dc_previous[continuing].copy(), lat_grid[continuing], stamp.month, cfg)
    if np.any(restarting):
        startup = wDC(fall_dc[restarting].copy(), winter_rain[restarting].copy(),
                      a=OVERWINTER_CARRY_OVER, b=OVERWINTER_WETTING_EFFICIENCY)
        own_dc[restarting] = DC(temperature[restarting], rain[restarting],
                               startup.copy(), lat_grid[restarting], stamp.month, cfg)
    # Yesterday's DC ends at yesterday's noon: today's full interval is winter.
    fall_dc[ending] = own_dc_previous[ending]
    winter_rain[ending] = rain[ending]
    inactive_continuation = season_seen & ~season_active & ~starting
    winter_rain[inactive_continuation] += rain[inactive_continuation]
    season_active = (season_active | starting) & ~ending
    season_seen = season_seen | starting
    fall_dc[starting] = np.nan
    winter_rain[starting] = 0.0
    if not np.all(np.isfinite(own_dc)) or np.any(own_dc < 0):
        raise ValueError('Overwintering DC calculation produced invalid values on land')
    if not np.all(np.isfinite(winter_rain)):
        raise ValueError('Overwintering winter precipitation accumulation is not finite')
    valid_fwi = season_active | ~season_seen
    result = np.full(n, np.nan, dtype=np.float64)
    if np.any(valid_fwi):
        bui = BUI(dmc_values[valid_fwi], own_dc[valid_fwi])
        result[valid_fwi] = FWI(isi_values[valid_fwi], bui)
        if not np.all(np.isfinite(result[valid_fwi])):
            raise ValueError('Overwintering FWI calculation produced non-finite values on active land')

    def full(values, flag=False):
        array = np.zeros(shape, dtype=np.int32) if flag else np.full(shape, np.nan, dtype=np.float64)
        array[active_land] = values
        return (('lat', 'lon'), array)

    state = xr.Dataset({
        'DCPrev': full(own_dc),
        'SeasonActive': full(season_active, flag=True),
        'DCf': full(fall_dc),
        'rw': full(winter_rain),
        'season_seen': full(season_seen, flag=True),
        'warm_days': full(warm_days, flag=True),
        'cold_days': full(cold_days, flag=True),
    }, coords={'lat': fields.lat, 'lon': fields.lon},
        attrs=dict(settings, state_date=stamp.isoformat()))
    return state, full(result)[1]


def process_local_day(stamp, fields, output, state, input_dir, write_components=WRITE_COMPONENTS,
                      overwintering=DO_OVERWINTERING):
    started = perf_counter()
    validate_overwintering_mode(state, overwintering)
    previous_ow = None if state is None else state.encoding.get('_fwi_overwinter_state')
    tas = fields.tas - 273.15
    mask = normalize_land_mask(fields.lsm, tas.lat.values, tas.lon.values) if 'lsm' in fields else None
    threshold = float(fields.attrs.get('land_threshold', LAND_THRESHOLD))
    land = None if mask is None else active_land(mask, threshold)
    select = (lambda values: values) if land is None else (lambda values: values[land])
    if state is None:
        print(f'Initializing zero moisture codes for {stamp}; spin-up required.', flush=True)
        state = xr.Dataset({k: xr.zeros_like(tas) for k in ('FFMCPrev', 'DMCPrev', 'DCPrev')})
    else:
        check_identity(state.attrs, fields.attrs, 'restart and weather')
        if state.attrs.get('noon_convention') != NOON_CONVENTION or state.attrs.get('sampling') != fields.attrs['sampling']:
            raise ValueError('Restart convention differs: do not reuse 12 UTC state or change sampling during a run')
        if state.attrs.get('state_date') != (stamp - timedelta(days=1)).isoformat():
            raise ValueError(f'Restart must represent the preceding local date for {stamp}')
        if 'lsm' in state:
            saved_mask = normalize_land_mask(state.lsm, tas.lat.values, tas.lon.values)
            if (mask is None or not np.array_equal(saved_mask.values, mask.values)
                    or state.attrs.get('land_threshold') != threshold):
                raise ValueError('Land-sea mask or threshold differs from saved constant')
    for name in ('FFMCPrev', 'DMCPrev', 'DCPrev'):
        xr.align(tas, state[name], join='exact')
        if state[name].dims != tas.dims or not np.isfinite(select(state[name].values)).all():
            raise ValueError(f'Invalid restart grid/values on active cells: {name}')
    for name in ('tas', 'relative_humidity', 'wind_speed', 'precip_24h'):
        xr.align(tas, fields[name], join='exact')
        if fields[name].dims != tas.dims or not np.isfinite(select(fields[name].values)).all():
            raise ValueError(f'Non-finite prepared weather on active cells: {name}')
    cfg = Config('hursmin-tasmax', 'NSHeq', 'continuous', 'original')
    temperatures = select(tas.values)
    lat = select(np.broadcast_to(tas.lat.values[:, None], tas.shape)).astype(float)
    day_of_year = stamp.timetuple().tm_yday
    doy = np.full(temperatures.shape, day_of_year, dtype=int)
    month = np.full(temperatures.shape, stamp.month, dtype=int)
    rain, rh, wind = (select(fields[name].values) for name in
                      ('precip_24h', 'relative_humidity', 'wind_speed'))
    if temperatures.size:
        ffmc = FFMC(temperatures, rh, wind, rain, select(state.FFMCPrev.values).copy())
        dmc = DMC(temperatures, rh, rain, select(state.DMCPrev.values).copy(), lat, doy, month, cfg)
        dc = DC(temperatures, rain, select(state.DCPrev.values).copy(), lat, month, cfg)
        isi, bui = ISI(wind, ffmc), BUI(dmc, dc)
        values = FWI(isi, bui)
    else:
        ffmc, dmc, dc, isi, bui, values = (np.empty(0, dtype=float) for _ in range(6))
    if not all(np.isfinite(v).all() for v in (ffmc, dmc, dc, isi, bui, values)):
        raise ValueError(f'Non-finite FWI/component output on active cells for {stamp}')

    def restore(values):
        if land is None:
            return values
        grid = np.full(tas.shape, np.nan, dtype=np.float64)
        grid[land] = values
        return grid

    next_ow = ow_values = None
    if overwintering:
        if previous_ow is None:
            print(f'Starting overwintering history on {stamp}; regular DC is used until a complete seasonal cycle is available.', flush=True)
        next_ow, ow_values = advance_overwintering(
            stamp, fields, previous_ow, state.DCPrev.values.copy(), restore(dc),
            restore(dmc), restore(isi), cfg, land=land)
    state = xr.Dataset({k: tas.copy(data=restore(v)) for k, v in
                        zip(('FFMCPrev', 'DMCPrev', 'DCPrev'), (ffmc, dmc, dc))})
    for name in state.data_vars:
        state[name].attrs = {'units': '1'}
    copy_selected_global_attrs(fields, state)
    state.attrs.update(state_date=stamp.isoformat(), noon_convention=NOON_CONVENTION,
                       sampling=fields.attrs['sampling'], overwintering_enabled=bool(overwintering))
    if next_ow is not None:
        state.encoding['_fwi_overwinter_state'] = next_ow
    if mask is not None:
        state['lsm'] = mask
        state.attrs['land_threshold'] = threshold
    nominal = pd.Timestamp(stamp) + pd.Timedelta(hours=12)
    out = tas.copy(data=restore(values)).astype('float32').expand_dims(time=[nominal]).to_dataset(name='fwi')
    out.fwi.attrs = {'units': '1', 'long_name': 'Fire Weather Index at local mean solar noon'}
    out['local_noon_utc'] = fields.local_noon_utc
    if mask is not None:
        out['lsm'] = mask
    out.attrs.update(fields.attrs)
    out = add_history_attr(out)
    out.time.attrs['comment'] = 'Nominal local-date label, not a common UTC valid time; see local_noon_utc(lon).'
    components = None
    if write_components:
        codes = {'FFMC': ffmc, 'DMC': dmc, 'DC': dc, 'ISI': isi, 'BUI': bui}
        long_names = {'FFMC': 'Fine Fuel Moisture Code', 'DMC': 'Duff Moisture Code',
                      'DC': 'Drought Code', 'ISI': 'Initial Spread Index', 'BUI': 'Buildup Index'}
        components = xr.Dataset({name: tas.copy(data=restore(code)).astype('float32').expand_dims(time=[nominal])
                                 for name, code in codes.items()})
        for name in codes:
            components[name].attrs = {'units': '1', 'long_name': long_names[name]}
        components['local_noon_utc'] = fields.local_noon_utc
        if mask is not None:
            components['lsm'] = mask
        components.attrs.update(out.attrs)
        components.time.attrs.update(out.time.attrs)
    report_timing('FWI calculations and output assembly', started)
    write_prepared_inputs(stamp, fields, input_dir)
    atomic_write(out, output / f'{stamp:%Y_%m_%d}_fwi.nc', engine='netcdf4',
                 encoding={'fwi': {'zlib': True, 'complevel': 5}})
    if components is not None:
        # All requested products must be written before advancing the constant.
        atomic_write(components, output / f'{stamp:%Y_%m_%d}_fwi_comp.nc', engine='netcdf4',
                     encoding={name: {'zlib': True, 'complevel': 5} for name in codes})
    if overwintering:
        overwintered = out.copy(deep=False)
        overwintered['fwi'] = out.fwi.copy(data=ow_values.astype('float32')[None, ...])
        overwintered.attrs = dict(out.attrs)
        overwintered.attrs.update(overwintering_enabled=True, **next_ow.attrs)
        atomic_write(overwintered, output / f'{stamp:%Y_%m_%d}_fwi_OverWintering.nc',
                     engine='netcdf4', encoding={'fwi': {'zlib': True, 'complevel': 5}})
    if stamp.day == 1:
        atomic_write(state, output / f'{stamp:%Y_%m_%d}_fwi_restart_file.nc', overwinter_state=next_ow)
    atomic_write(state, output / 'FWI_Const1.nc', overwinter_state=next_ow)
    rain_range = f'{rain.min():.4g}..{rain.max():.4g} mm' if rain.size else 'no selected land cells'
    print(f'Completed local date {stamp}; DOY={day_of_year}; rain range {rain_range}', flush=True)
    return state, out


def validate_overwintering_mode(state, enabled):
    if not isinstance(enabled, (bool, np.bool_)):
        raise ValueError('overwintering must be True or False')
    if enabled:
        validate_overwintering_settings()
    if state is None:
        return
    previous = state.attrs.get('overwintering_enabled')
    if previous is not None and previous != bool(enabled):
        raise ValueError('Overwintering on/off setting differs from saved constant; use a separate output integration for a different mode')
    stored = state.encoding.get('_fwi_overwinter_state')
    if previous is True and stored is None:
        raise ValueError('Saved constant is missing its overwintering state group')
    if not enabled and stored is not None:
        raise ValueError('Overwintering state exists but overwintering is disabled')



def check_identity(left, right, context):
    """Only model, member (realization), resolution are identity constraints."""
    for key in COMPATIBILITY_KEYS:
        if key not in left or key not in right:
            raise ValueError(f'Missing {key} metadata: {context}')
        if str(left[key]).strip().lower() != str(right[key]).strip().lower():
            raise ValueError(f'{key} mismatch: {context}: {left[key]!r} vs {right[key]!r}')


def check_boundary_source(path, hours, attrs):
    context = HANDOVER_CONTEXT
    branch = 'historical' if path.parent.resolve() in context['historical_dirs'] else 'projection'
    check_identity(context['identity'], attrs, str(path))
    cut = (np.datetime64(context['start'], 'h') - np.datetime64('1970-01-01', 'h')).astype(int)
    # Timestamp ownership, including end-labelled precipitation: January 1 00
    # is supplied by the projection, even though its interval starts in history.
    valid = hours < cut if branch == 'historical' else hours >= cut
    if not np.all(valid):
        raise ValueError(f'{branch} UTC timestamps overlap the other experiment at {context["start"]}: {path}')
    context['attrs'].setdefault(branch, dict(attrs))


@contextmanager
def boundary_weather(args, stamp, start, identity):
    global HANDOVER_CONTEXT
    historical = set()
    for directory in args.historical_input_dir:
        historical.update(hourly_directories(directory, stamp, monthly=args.monthly_chunk))
    projection = hourly_directories(args.hpctmpdir, stamp, args.extra_input_dir, args.monthly_chunk)
    if historical.intersection(projection):
        raise ValueError('Historical and projection input directories must be separate')
    previous = HANDOVER_CONTEXT
    HANDOVER_CONTEXT = dict(start=start, identity=identity, historical_dirs=historical, attrs={})
    try:
        yield projection + sorted(historical)
    finally:
        HANDOVER_CONTEXT = previous


def prepare_handover_fields(args, stamp, start, state):
    with boundary_weather(args, stamp, start, state.attrs) as directories:
        return prepare_local_inputs(stamp, directories, args.sampling, args.precip_time_label,
            lsm_file=args.lsm_file, land_threshold=args.land_threshold, land_mask=state.get('lsm'),
            overwintering=args.overwintering, io_tile_columns=args.io_tile_columns,
            io_buffer_mb=args.io_buffer_mb, io_tile_rows=args.io_tile_rows)


def link_historical_year(source, destination, last_day, identity):
    """Read-only references avoid copying a full year of global fields per branch."""
    for suffix in ('fwi', 'fwi_comp', 'fwi_OverWintering'):
        for path in sorted(source.glob(f'{last_day.year}_??_??_{suffix}.nc')):
            stamp = date.fromisoformat(path.name[:10].replace('_', '-'))
            if stamp >= last_day:
                continue
            with xr.open_dataset(path) as ds:
                check_identity(identity, ds.attrs, str(path))
            target = destination / path.name
            if target.is_symlink():
                if target.resolve() != path.resolve():
                    raise ValueError(f'Conflicting historical link: {target}')
            elif target.exists():
                raise ValueError(f'Refusing to replace existing historical product: {target}')
            else:
                target.symlink_to(path.resolve())


def scenario_output_paths(args, dates, base):
    """Route by the projection's own input metadata, not historical companions."""
    token = re.compile(r'(?<![A-Za-z0-9])2t(?![A-Za-z0-9])')
    paths = sorted(p for directory in hourly_directories(args.hpctmpdir, dates[0], monthly=args.monthly_chunk)
                   for p in Path(directory).glob('*.nc')
                   if token.search(p.name) and hourly_file_overlaps_day(p, dates[0]))
    if not paths:
        raise FileNotFoundError('Cannot determine experiment: no temperature input for requested date')
    with open_hourly_input(paths[0], '2t') as ds:
        attrs = {}
        for da in ds.data_vars.values():
            if 'time' in da.dims and 'lat' in da.dims and 'lon' in da.dims:
                attrs.update(da.attrs)
        attrs.update(ds.attrs)
    experiment = str(attrs.get('experiment', '')).strip()
    if not experiment or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', experiment):
        raise ValueError('Missing or unsafe experiment label in temperature metadata')
    if experiment.lower() == 'hist':
        if getattr(args, 'historical_output_dir', None):
            raise ValueError('Handover requires projection inputs, not historical inputs')
        return base, None
    output = base / 'projections' / experiment
    boundary = base / 'historical_boundary' / experiment
    return output, boundary


def execute_with_handover(args, dates, directories, output):
    base = Path(output)
    output, boundary = scenario_output_paths(args, dates, base)
    args = argparse.Namespace(**vars(args))
    args.scenario_boundary_dir = boundary
    legacy = base / 'FWI_Const1.nc'
    if output != base and legacy.exists():
        with xr.open_dataset(legacy) as ds:
            legacy_experiment = str(ds.attrs.get('experiment', '')).lower()
        if legacy_experiment != 'hist':
            raise ValueError('Shared output root already contains a non-historical constant; '
                             'migrate that existing projection explicitly before using scenario folders')
        if not (output / 'FWI_Const1.nc').exists() and not getattr(args, 'historical_output_dir', None):
            raise ValueError('Historical state exists in the shared root: supply handover arguments '
                             'to start this projection without reinitializing')
    output.mkdir(parents=True, exist_ok=True)
    print(f'Scenario output directory: {output}', flush=True)
    # Main holds the shared root lock; this also protects direct programmatic calls.
    if output == base:
        return execute_scenario_handover(args, dates, directories, output)
    with integration_lock(output):
        return execute_scenario_handover(args, dates, directories, output)


def execute_scenario_handover(args, dates, directories, output):
    """Branch completion is durable and idempotent; history is never modified."""
    saved_start = output / 'FWI_projection_start_date.txt'
    requested = getattr(args, 'historical_output_dir', None)
    if not requested:
        if saved_start.exists() and dates[0] <= date.fromisoformat(saved_start.read_text().strip()):
            raise ValueError('Boundary retry requires --historical-output-dir and --historical-input-dir')
        return execute_chunk(args, dates, directories, output)
    if args.initialize or args.restart_file or args.prepare_only:
        raise ValueError('Handover supplies the restart; do not combine with initialize/restart-file/prepare-only')
    if not args.historical_input_dir:
        raise ValueError('Handover requires --historical-input-dir')
    for directory in args.historical_input_dir:
        if not directory.is_dir():
            raise FileNotFoundError(directory)
    historical = requested.resolve()
    if historical == output.resolve():
        raise ValueError('Projection output must be separate from the historical output directory')
    if saved_start.exists():
        start = date.fromisoformat(saved_start.read_text().strip())
        if args.projection_start_date and args.projection_start_date != start:
            raise ValueError('Projection start date differs from saved handover')
    else:
        start = args.projection_start_date or dates[0]
        if dates[0] != start:
            raise ValueError('Initial handover must run the first projection date')
        if (output / 'FWI_Const1.nc').exists() or list(output.glob('*_fwi.nc')):
            raise ValueError('Initial handover needs a fresh projection output directory')
        saved_start.write_text(start.isoformat()+'\n')
    if dates[0] < start:
        raise ValueError('Requested projection date precedes its start')
    if args.final_local_date is not None and args.final_local_date < start:
        raise ValueError('Projection final local date precedes projection start')
    last_day = start - timedelta(days=1)
    boundary = getattr(args, 'scenario_boundary_dir', None) or output / 'historical_boundary'
    boundary.mkdir(parents=True, exist_ok=True)
    constant = output / 'FWI_Const1.nc'
    branch_constant = boundary / 'FWI_Const1.nc'
    if dates[0] > start and constant.exists():
        # Ordinary continuation must not reread the large historical state.
        return execute_chunk(args, dates, directories, output)
    with integration_lock(boundary):
        if not branch_constant.exists():
            seed = historical / 'FWI_Const1.nc'
            with open_product_dataset(seed) as ds:
                state = ds.load()
            if state.attrs.get('state_date') != (last_day-timedelta(days=1)).isoformat():
                raise ValueError(f'Historical constant must end on {last_day-timedelta(days=1)}; final day {last_day} is completed per branch')
            validate_overwintering_mode(state, args.overwintering)
            fields = prepare_handover_fields(args, last_day, start, state)
            # No shared OPA prepared file is overwritten by another scenario:
            # prepared boundary weather is written inside this projection OPA dir.
            state, _ = process_local_day(last_day, fields, boundary, state, args.hpctmpdir,
                                        args.write_components, args.overwintering)
        with open_product_dataset(branch_constant) as ds:
            boundary_state = ds.load()
        if boundary_state.attrs.get('state_date') != last_day.isoformat():
            raise ValueError('Invalid saved historical boundary date')
        validate_overwintering_mode(boundary_state, args.overwintering)
        # Rebuild statistics on initial/retried boundary only, never each later day.
        if not constant.exists() or dates[0] == start:
            link_historical_year(historical, boundary, last_day, boundary_state.attrs)
            if DO_ANNUAL and not args.no_annual:
                annual_statistics(str(boundary), last_day.year, boundary_state, coverage_end=last_day)
            if args.annual_components:
                annual_component_statistics(str(boundary), last_day.year, boundary_state, coverage_end=last_day)
            if args.overwintering and getattr(args, 'annual_overwintering', DO_ANNUAL_OVERWINTERING):
                annual_overwintering_statistics(str(boundary), last_day.year, boundary_state, coverage_end=last_day)
        if not constant.exists():
            atomic_write(boundary_state, constant,
                         overwinter_state=boundary_state.encoding.get('_fwi_overwinter_state'))
    # The first projection day needs historical rain and, for OW, temperatures.
    # It is processed explicitly; later days retain the ordinary chunk pipeline.
    if start in dates:
        with open_product_dataset(constant) as ds:
            state = ds.load()
        saved = date.fromisoformat(state.attrs['state_date'])
        if saved == last_day:
            fields = prepare_handover_fields(args, start, start, state)
            process_local_day(start, fields, output, state, args.hpctmpdir,
                              args.write_components, args.overwintering)
        elif saved < start:
            raise ValueError('Projection state is not continuous with historical boundary')
    return execute_chunk(args, dates, directories, output)


def main():
    started = perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    ow_switch = parser.add_mutually_exclusive_group()
    ow_switch.add_argument('--overwintering', dest='overwintering', action='store_true',
                           help='Also write *_fwi_OverWintering.nc; regular FWI always continues')
    ow_switch.add_argument('--no-overwintering', dest='overwintering', action='store_false',
                           help='Write regular FWI only, without reading the additional temperature history')
    parser.set_defaults(overwintering=DO_OVERWINTERING)
    parser.add_argument('--year', type=int, required=True)
    parser.add_argument('--month', type=int, required=True)
    parser.add_argument('--day', type=int)
    parser.add_argument('--hpctmpdir', '--input-dir', dest='hpctmpdir', required=True)
    parser.add_argument('--hpcrootdir', default='')
    parser.add_argument('--hpcprojdir', default='')
    parser.add_argument('--output-dir')
    parser.add_argument('--extra-input-dir', action='append', default=[], help='Additional adjacent-day hourly directory; repeat if needed')
    parser.add_argument('--monthly-chunk','--monthly_chunk', action='store_true', default=MONTHLY_CHUNK)
    parser.add_argument('--sampling', choices=('linear','nearest'), default='linear')
    parser.add_argument('--precip-time-label', choices=('start','end'), default=PRECIP_TIME_LABEL,
                        help='Timestamp convention when precipitation lacks time bounds: end = preceding hour, start = following hour (default: %(default)s); CF bounds take precedence')
    parser.add_argument('--restart-file', type=Path)
    parser.add_argument('--initialize', action='store_true', help='New integration only; start on given day if previous 24h inputs exist')
    parser.add_argument('--no-annual', action='store_true',
                        help='Disable annual FWI statistics only; component annual statistics have a separate switch')
    annual_component_switch = parser.add_mutually_exclusive_group()
    annual_component_switch.add_argument('--annual-components', dest='annual_components', action='store_true',
                                         help='Enable the single annual component file (script default: DO_ANNUAL_COMPONENTS)')
    annual_component_switch.add_argument('--no-annual-components', dest='annual_components', action='store_false',
                                         help='Disable annual component statistics independently of FWI')
    parser.set_defaults(annual_components=DO_ANNUAL_COMPONENTS)
    ow_annual = parser.add_mutually_exclusive_group()
    ow_annual.add_argument('--annual-overwintering', dest='annual_overwintering', action='store_true',
                           help='Write annual statistics of the daily overwintered FWI stream')
    ow_annual.add_argument('--no-annual-overwintering', dest='annual_overwintering', action='store_false')
    parser.set_defaults(annual_overwintering=DO_ANNUAL_OVERWINTERING)
    parser.add_argument('--lsm-file', type=Path, default=LSM_FILE,
                        help='Initial lsm file; otherwise discover in input directories. Once saved, reuse the state mask without checking input masks.')
    parser.add_argument('--land-threshold', type=float, default=LAND_THRESHOLD,
                        help='Land fraction cutoff (default: %(default)s); 0 retains any positive land fraction')
    component_switch = parser.add_mutually_exclusive_group()
    component_switch.add_argument('--write-components', dest='write_components', action='store_true',
                              help='Write all five daily components in YYYY_MM_DD_fwi_comp.nc')
    component_switch.add_argument('--no-write-components', dest='write_components', action='store_false',
                              help='Disable new daily component files; annual component control is separate')
    parser.set_defaults(write_components=WRITE_COMPONENTS)
    parser.add_argument('--final-local-date', type=date.fromisoformat, default=FINAL_LOCAL_DATE,
                        help='Last complete local date (YYYY-MM-DD) of the entire integration, not each chunk; allows labeled partial final-year statistics and is saved for retries')
    parser.add_argument('--prepare-only', action='store_true', help='Write local-noon weather fields in hpctmpdir (OPA) only; do not read/write FWI state')
    parser.add_argument('--io-tile-columns', type=int, default=IO_TILE_COLUMNS,
                        help='Maximum contiguous longitude columns in one serial read tile (default: %(default)s)')
    parser.add_argument('--io-buffer-mb', type=float, default=IO_BUFFER_MB,
                        help='Decoded hourly read-block budget in MiB, excluding global fields and decoding copies (default: %(default)s)')
    parser.add_argument('--io-tile-rows', type=int, default=IO_TILE_ROWS,
                        help='Maximum latitude rows per strip for contiguous rain/Tmax input files (default: %(default)s)')
    parser.add_argument('--historical-output-dir', type=Path,
                        help='Read-only historical daily products and December 30 (or boundary-minus-one) constant')
    parser.add_argument('--historical-input-dir', type=Path, action='append', default=[],
                        help='Historical hourly boundary directory; repeat if required')
    parser.add_argument('--projection-start-date', type=date.fromisoformat,
                        help='First projection UTC date; defaults to first requested date on initial handover, then saved')
    args = parser.parse_args()
    validate_io_settings(args.io_tile_columns, args.io_buffer_mb)
    validate_io_settings(args.io_tile_rows, args.io_buffer_mb)
    active_land(np.zeros(1), args.land_threshold)
    validate_overwintering_mode(None, args.overwintering)
    # Validate an explicit input mask during initial discovery only. Once a
    # saved mask exists, continuation must not depend on that input file.
    dates = chunk_dates(args.year,args.month,args.day,args.monthly_chunk)
    directories = [args.hpctmpdir] + args.extra_input_dir
    for directory in directories:
        if not Path(directory).is_dir():
            raise FileNotFoundError(f'Input directory missing: {directory}')
    output = Path(args.output_dir or build_output_dir_from_input(args.hpctmpdir.rstrip('/')))
    output.mkdir(parents=True, exist_ok=True)
    with integration_lock(output):
        execute_with_handover(args, dates, directories, output)
    report_timing('complete application', started)


@contextmanager
def integration_lock(output):
    # Shared constant file must have only one writer. flock is released on exit/crash.
    with (output / '.fwi_local_noon.lock').open('a') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Another local-noon job is using this output directory') from exc
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def write_pending_day(path, stamp):
    """Publish the small pending-date marker atomically, without JSON sidecars."""
    fd, temporary = tempfile.mkstemp(prefix=path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(stamp.isoformat() + '\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def execute_chunk(args, dates, directories, output):
    """Monthly workflow: finish pending predecessor, attempt this month end.

    Defer only when required next-day boundary input is unavailable.
    The state is authoritative; the marker records work still owed, not success
    for that local day. No meteorological exception is suppressed.
    """
    global BOUNDARY_WAIT_HOUR
    if args.prepare_only:
        return _execute_chunk(args, dates, directories, output)
    dates = list(dates)
    pending_path = output / 'FWI_pending_local_date.txt'
    pending = date.fromisoformat(pending_path.read_text().strip()) if pending_path.exists() else None
    requested_first = dates[0]
    if pending is not None and args.final_local_date == pending - timedelta(days=1):
        # Explicitly close a terminal integration without inventing its unavailable day.
        _execute_chunk(args, dates, directories, output)
        pending_path.unlink()
        return
    if pending is not None and requested_first > pending:
        next_day = pending + timedelta(days=1)
        if (requested_first.year, requested_first.month) != (next_day.year, next_day.month):
            raise ValueError(f'Pending local date {pending} must be completed by the immediately following month')
        constant = output / 'FWI_Const1.nc'
        with open_product_dataset(constant) as ds:
            saved = date.fromisoformat(ds.attrs['state_date'])
        if saved not in (pending - timedelta(days=1), pending):
            raise ValueError(f'Cannot complete pending {pending}: state is {saved}; earlier days are incomplete')
        if args.initialize or args.restart_file:
            raise ValueError('Pending continuation must use its saved constant, not initialize/restart-file')
        finish = argparse.Namespace(**vars(args))
        finish.initialize = False
        finish.restart_file = None
        # Keep monthly directory discovery, but bypass this wrapper for exactly
        # the pending date. It also retries annual output if state already exists.
        print(f'Completing pending local date {pending} before {requested_first:%Y-%m}.', flush=True)
        _execute_chunk(finish, [pending], directories, output)
        pending_path.unlink()
        pending = None
    elif pending is not None and (requested_first.year, requested_first.month) != (pending.year, pending.month):
        raise ValueError(f'Requested chunk conflicts with pending local date {pending}')

    # An explicit final date retains the original finalization semantics, e.g.
    # final-local-date=December 30 writes disclosed partial final-year statistics.
    final_marker = output / 'FWI_final_local_date.txt'
    declared_end = args.final_local_date
    if declared_end is None and final_marker.exists():
        declared_end = date.fromisoformat(final_marker.read_text().strip())
    last = dates[-1]
    defer = (args.monthly_chunk and (declared_end is None or last < declared_end)
             and last.day == calendar.monthrange(last.year, last.month)[1])
    if defer:
        if pending is not None and pending != last:
            raise ValueError('Different month-end date already pending')
        # A successful retry must not turn an already completed day into pending.
        constant = output / 'FWI_Const1.nc'
        saved = None
        if constant.exists():
            with open_product_dataset(constant) as ds:
                saved = date.fromisoformat(ds.attrs['state_date'])
        if saved is not None and saved >= last:
            result = _execute_chunk(args, dates, directories, output)
            if pending_path.exists():
                pending_path.unlink()
            return result
        write_pending_day(pending_path, last)
        work_dates = dates[:-1]
        if work_dates:
            _execute_chunk(args, work_dates, directories, output)
        with open_product_dataset(constant) as ds:
            if ds.attrs.get('state_date') != (last - timedelta(days=1)).isoformat():
                raise ValueError(f'Cannot defer {last}: preceding local day is not complete')
        finish = argparse.Namespace(**vars(args))
        finish.initialize = False
        finish.restart_file = None
        previous_wait = BOUNDARY_WAIT_HOUR
        BOUNDARY_WAIT_HOUR = int((np.datetime64(last + timedelta(days=1), 'h')
                                 - np.datetime64('1970-01-01', 'h')).astype(int))
        try:
            _execute_chunk(finish, [last], directories, output)
        except BoundaryInputNotReady as exc:
            print(f'Boundary input pending: {exc}', flush=True)
        else:
            pending_path.unlink()
            print(f'Completed month-end {last}: required boundary data are available.', flush=True)
            return
        finally:
            BOUNDARY_WAIT_HOUR = previous_wait
        print(f'Deferred local date {last}; saved state through {last-timedelta(days=1)}. '
              'The next monthly job will complete it before its own dates. '
              'Annual statistics for this year remain pending if this is December 31.', flush=True)
        return
    if pending is not None:
        raise ValueError('Pending month-end retry requires the original monthly chunk; '
                         'do not change final-date settings until pending work is resolved')
    return _execute_chunk(args, dates, directories, output)


def _execute_chunk(args, dates, directories, output):
    if args.initialize and args.restart_file:
        raise ValueError('Choose --initialize or --restart-file')
    marker = output / 'FWI_start_date_local_noon.txt'
    final_marker = output / 'FWI_final_local_date.txt'
    constant = output / 'FWI_Const1.nc'
    legacy_constant = output / 'FWI_Const_local_noon.nc'
    default_restart = constant if constant.exists() else legacy_constant
    products = list(output.glob('*_fwi.nc')) + list(output.glob('*_fwi_comp.nc')) + list(output.glob('*_fwi_OverWintering.nc'))
    checkpoints = list(output.glob('*_fwi_restart_file.nc'))
    if args.initialize and (constant.exists() or legacy_constant.exists() or products or checkpoints):
        raise ValueError('--initialize cannot reset an existing local-noon integration')
    saved_end = date.fromisoformat(final_marker.read_text().strip()) if final_marker.exists() else None
    if saved_end is not None and args.final_local_date is not None and saved_end != args.final_local_date:
        raise ValueError(f'Final local date conflicts with saved declaration {saved_end}')
    declared_end = args.final_local_date or saved_end
    if declared_end is not None and dates[0].year > declared_end.year:
        raise ValueError(f'Requested year is beyond the declared final local date {declared_end}')
    finalization_only = declared_end is not None and dates[0] > declared_end
    if declared_end is not None:
        excluded = [stamp for stamp in dates if stamp > declared_end]
        if excluded:
            print(f'Excluding {len(excluded)} requested date(s) after declared final local date {declared_end}.', flush=True)
        dates = [stamp for stamp in dates if stamp <= declared_end]
        if finalization_only:
            if args.prepare_only:
                return
            if args.initialize:
                raise ValueError('Cannot initialize after the declared final local date')
            # A final daily job (e.g. Dec31) can finalize an already completed Dec30.
            # It must never manufacture the missing final date or advance the state.
            dates = [declared_end]
    state = None
    if not args.prepare_only and not args.initialize:
        restart = args.restart_file or default_restart
        if restart.exists():
            with open_product_dataset(restart) as ds:
                state = ds.load()
    if state is not None:
        validate_overwintering_mode(state, args.overwintering)
        if 'lsm' in state and state.attrs.get('land_threshold') != args.land_threshold:
            raise ValueError('Land threshold differs from saved constant')
        if state.attrs.get('noon_convention') != NOON_CONVENTION or state.attrs.get('sampling') != args.sampling:
            raise ValueError('Restart noon convention/sampling mismatch')
        saved_date = date.fromisoformat(state.attrs['state_date'])
        if not args.restart_file and constant.exists() and legacy_constant.exists():
            with open_product_dataset(legacy_constant) as older:
                if older.attrs.get('noon_convention') != NOON_CONVENTION or older.attrs.get('sampling') != args.sampling:
                    raise ValueError('Legacy and current constant conventions disagree')
                if date.fromisoformat(older.attrs['state_date']) > saved_date:
                    raise ValueError('Legacy constant is newer than FWI_Const1.nc; refusing to resume older state')
        if declared_end is not None and saved_date > declared_end:
            raise ValueError(f'Constant is later than declared final local date {declared_end}; refusing to truncate completed work')
    if finalization_only and (state is None or state.attrs['state_date'] != declared_end.isoformat()):
        raise FileNotFoundError(f'Final local date {declared_end} has not been completed; run it before finalizing statistics')
    if declared_end is not None and not args.prepare_only and not final_marker.exists():
        with final_marker.open('x') as stream:
            stream.write(declared_end.isoformat()+'\n')
    start = date.fromisoformat(marker.read_text().strip()) if marker.exists() else None
    if not args.prepare_only and not args.restart_file and not args.initialize and not finalization_only and EXCLUDE_FIRST_DAY and start is None and not constant.exists() and not legacy_constant.exists() and not products and not checkpoints:
        start = dates[0]
        with marker.open('x') as stream:
            stream.write(start.isoformat()+'\n')
    if start is not None and not args.prepare_only:
        if dates[0] < start:
            raise ValueError(
                f'Requested date {dates[0].isoformat()} precedes saved simulation start '
                f'{start.isoformat()} in {marker}. Check the requested chunk dates and '
                'output directory. A fresh integration needs a separate output directory; '
                'do not remove an existing constant or change the start marker to bypass this check.')
        if start in dates:
            print(f'Skipping only initial local date {start}: no preceding 24h rainfall.', flush=True)
            dates.remove(start)
        if not dates:
            return
    if not args.prepare_only and not args.initialize:
        if state is None and (args.restart_file or products or checkpoints or (start is not None and dates[0] != start+timedelta(days=1))):
            raise FileNotFoundError('Local-noon constant missing for continuation; do not substitute 12 UTC state')
    last_output = None
    final_date = dates[-1]
    if state is not None:
        saved_date = date.fromisoformat(state.attrs['state_date'])
        if saved_date > final_date:
            raise ValueError('Constant is later than the requested chunk; refusing to rewind it')
        completed = [stamp for stamp in dates if stamp <= saved_date]
        available_completed = []
        for stamp in completed:
            product = existing_daily_product(output, stamp)
            # The validated constant determines recurrence continuity. Removed daily
            # outputs are disclosed as missing by annual aggregation, not recomputed.
            if product.exists():
                available_completed.append(product)
        if completed:
            print(f'Resuming after completed local date {saved_date}', flush=True)
            dates = [stamp for stamp in dates if stamp > saved_date]
            last_output = xr.Dataset(attrs=dict(state.attrs))
            if available_completed:
                with open_product_dataset(available_completed[-1]) as ds:
                    last_output = xr.Dataset(attrs=dict(ds.attrs))
    for stamp in dates:
        # Discover again each day because Autosubmit produces chunks incrementally.
        day_directories = hourly_directories(args.hpctmpdir, stamp, args.extra_input_dir, args.monthly_chunk)
        print(f"Hourly inputs for local date {stamp}: " + "; ".join(map(str, day_directories)), flush=True)
        fields = prepare_local_inputs(stamp,day_directories,args.sampling,args.precip_time_label,
                                      lsm_file=args.lsm_file, land_threshold=args.land_threshold,
                                      land_mask=state.get('lsm') if state is not None else None,
                                      overwintering=args.overwintering,
                                      io_tile_columns=args.io_tile_columns, io_buffer_mb=args.io_buffer_mb,
                                      io_tile_rows=args.io_tile_rows)
        if args.prepare_only:
            write_prepared_inputs(stamp, fields, args.hpctmpdir)
        else:
            state,last_output = process_local_day(stamp,fields,output,state,input_dir=args.hpctmpdir,
                                                  write_components=args.write_components, overwintering=args.overwintering)
    annual_boundary = (final_date.month == 12 and final_date.day == 31) or final_date == declared_end
    if not args.prepare_only and annual_boundary and DO_ANNUAL and not args.no_annual:
        if final_date == declared_end:
            annual_statistics(str(output),final_date.year,last_output,coverage_end=declared_end)
        else:
            annual_statistics(str(output),final_date.year,last_output)
    if not args.prepare_only and annual_boundary and args.annual_components:
        if final_date == declared_end:
            annual_component_statistics(str(output), final_date.year, last_output, coverage_end=declared_end)
        else:
            annual_component_statistics(str(output), final_date.year, last_output)

    if (not args.prepare_only and annual_boundary and args.overwintering
            and getattr(args, 'annual_overwintering', DO_ANNUAL_OVERWINTERING)):
        annual_overwintering_statistics(str(output), final_date.year, last_output,
                                       coverage_end=declared_end if final_date == declared_end else None)


def annual_overwintering_file_date(path):
    match = re.fullmatch(r'(\d{4})_(\d{2})_(\d{2})_fwi_OverWintering\.nc', Path(path).name)
    if not match:
        raise ValueError(f'Malformed daily overwintered FWI filename: {path}')
    return date(*map(int, match.groups()))


def validate_annual_overwintering_files(out_path, year_in, expected_times):
    allowed = {stamp.date() for stamp in expected_times}
    products = {}
    for path in sorted(Path(out_path).glob(f'{year_in}*_fwi_OverWintering.nc')):
        if not path.is_file():
            continue
        stamp = annual_overwintering_file_date(path)
        if stamp not in allowed or stamp in products:
            raise ValueError(f'Duplicate or outside-coverage overwintered FWI date: {path}')
        products[stamp] = str(path)
    if not products:
        raise ValueError(f'No available daily overwintered FWI products for {year_in} in {out_path}')
    return [products[stamp] for stamp in sorted(products)]


def validate_annual_overwintering_dataset(ds):
    validate_annual_daily_dataset(ds, annual_overwintering_file_date)
    if not ds.attrs.get('overwintering_enabled', False):
        raise ValueError('Daily overwintered FWI is missing its seasonal-stream metadata')
    return ds


@contextmanager
def open_annual_overwintering_dataset(file_list, chunks=None, expected_sampling=None):
    with _open_annual_product_dataset(file_list, ('fwi',), annual_overwintering_file_date,
                                     validate_annual_overwintering_dataset, chunks, expected_sampling,
                                     allow_seasonal_nan=True) as ds:
        yield ds


def reduce_overwintering_statistics(fwi_year, metadata, quantile_input=None, land_mask=None):
    # Winter NaNs are intentional: they must never become zero-FWI observations.
    basic, quantiles = reduce_annual_statistics(fwi_year, metadata, quantile_input, land_mask)
    counts = fwi_year.notnull().sum('time').astype('int16')
    if land_mask is not None:
        mask = normalize_land_mask(land_mask, fwi_year.lat.values, fwi_year.lon.values)
        active = xr.DataArray(active_land(mask, float(metadata.get('land_threshold', LAND_THRESHOLD))),
                              dims=('lat', 'lon'), coords={'lat': mask.lat, 'lon': mask.lon})
        counts = counts.where(active)
    basic['days_valid'] = counts
    basic.days_valid.attrs = {'units': '1', 'long_name': 'Days with a finite overwintered FWI value',
                              'comment': 'Seasonal winter NaNs are excluded; file-date coverage is recorded separately.'}
    return basic, quantiles


def annual_overwintering_statistics(out_path, year_in, FWI_all, coverage_end=None):
    """Summarize available local dates and disclose every absent calendar date."""
    allowed_times, _ = annual_coverage(year_in, coverage_end)
    file_list = validate_annual_overwintering_files(out_path, year_in, allowed_times)
    expected_times, coverage = annual_coverage(
        year_in, coverage_end, available_dates=[annual_overwintering_file_date(path) for path in file_list])
    import dask

    label = 'partial-year' if coverage['is_partial_year'] else 'full-year'
    print(f"Computing {label} overwintered FWI statistics for {year_in}: {coverage['days_used']} available days, "
          f"{coverage['missing_days']} missing calendar days", flush=True)
    # Both the source handles and temporary Dask settings are released on errors.
    with dask.config.set({'array.slicing.split_large_chunks': True}):
        with open_annual_overwintering_dataset(file_list, chunks={'time': 15},
                                 expected_sampling=FWI_all.attrs.get('sampling')) as ds_year:
            fwi_year = ds_year['fwi']
            validate_annual_timestamps(fwi_year, expected_times)
            spatial_dims = [dim for dim in fwi_year.dims if dim != 'time']
            fwi_year = fwi_year.chunk({'time': 15, **{dim: 200 for dim in spatial_dims}}).astype('float32')
            fwi_quantile = fwi_year.chunk({'time': -1, **{dim: 100 for dim in spatial_dims}})
            metadata = dict(ds_year.attrs)
            metadata.update(FWI_all.attrs)
            metadata.update(coverage)
            metadata.update(noon_convention=NOON_CONVENTION, sampling=ds_year.attrs['sampling'])
            if 'lsm' in ds_year:
                metadata['land_threshold'] = ds_year.attrs['land_threshold']
            basic, quantiles = reduce_overwintering_statistics(fwi_year, metadata, fwi_quantile,
                                                       land_mask=ds_year.get('lsm'))
            write_annual_products(out_path, year_in, basic, quantiles, suffix='_OverWintering')
    print(f'{label.capitalize()} overwintered FWI statistics written OK.', flush=True)


def annual_statistics(out_path, year_in, FWI_all, coverage_end=None):
    """Summarize available local dates and disclose every absent calendar date."""
    if not DO_ANNUAL:
        print('Annual statistics disabled in script.', flush=True)
        return
    allowed_times, _ = annual_coverage(year_in, coverage_end)
    file_list = validate_annual_files(out_path, year_in, allowed_times)
    expected_times, coverage = annual_coverage(
        year_in, coverage_end, available_dates=[annual_file_date(path) for path in file_list])
    import dask

    label = 'partial-year' if coverage['is_partial_year'] else 'full-year'
    print(f"Computing {label} statistics for {year_in}: {coverage['days_used']} available days, "
          f"{coverage['missing_days']} missing calendar days", flush=True)
    # Both the source handles and temporary Dask settings are released on errors.
    with dask.config.set({'array.slicing.split_large_chunks': True}):
        with open_annual_dataset(file_list, chunks={'time': 15},
                                 expected_sampling=FWI_all.attrs.get('sampling')) as ds_year:
            fwi_year = ds_year['fwi']
            validate_annual_timestamps(fwi_year, expected_times)
            spatial_dims = [dim for dim in fwi_year.dims if dim != 'time']
            fwi_year = fwi_year.chunk({'time': 15, **{dim: 200 for dim in spatial_dims}}).astype('float32')
            fwi_quantile = fwi_year.chunk({'time': -1, **{dim: 100 for dim in spatial_dims}})
            metadata = dict(ds_year.attrs)
            metadata.update(FWI_all.attrs)
            metadata.update(coverage)
            metadata.update(noon_convention=NOON_CONVENTION, sampling=ds_year.attrs['sampling'])
            if 'lsm' in ds_year:
                metadata['land_threshold'] = ds_year.attrs['land_threshold']
            basic, quantiles = reduce_annual_statistics(fwi_year, metadata, fwi_quantile,
                                                       land_mask=ds_year.get('lsm'))
            write_annual_products(out_path, year_in, basic, quantiles)
    print(f'{label.capitalize()} statistics written OK.', flush=True)


def annual_coverage(year_in, coverage_end=None, available_dates=None):
    """Return observed timestamps and missing dates over the FULL calendar year."""
    first, last = date(year_in, 1, 1), date(year_in, 12, 31)
    end = last if coverage_end is None else coverage_end
    if not isinstance(end, date) or end < first or end > last:
        raise ValueError(f'Annual coverage end must be a date within {year_in}')
    if available_dates is None:
        available_dates = [stamp.date() for stamp in pd.date_range(first, end)]
    observed = list(available_dates)
    if not observed:
        raise ValueError(f'No available daily FWI products for {year_in}')
    if any(type(stamp) is not date or stamp < first or stamp > end for stamp in observed):
        raise ValueError(f'Available dates must be within declared coverage for {year_in}')
    observed_set = set(observed)
    if len(observed) != len(observed_set):
        raise ValueError('Duplicate dates in annual coverage')
    observed.sort()
    expected = pd.DatetimeIndex(observed) + pd.Timedelta(hours=12)
    excluded = pd.date_range(end + timedelta(days=1), last)
    missing = [stamp.date() for stamp in pd.date_range(first, last) if stamp.date() not in observed_set]
    days_in_year = 366 if calendar.isleap(year_in) else 365
    return expected, {
        'coverage_start': observed[0].isoformat(),
        'coverage_end': observed[-1].isoformat(),
        'declared_coverage_end': end.isoformat(),
        'days_used': len(expected),
        'days_in_year': days_in_year,
        'coverage_fraction': len(expected) / days_in_year,
        'is_partial_year': int(bool(missing)),
        'missing_days': len(missing),
        'missing_dates': ','.join(stamp.isoformat() for stamp in missing),
        'excluded_dates': ','.join(stamp.date().isoformat() for stamp in excluded),
    }


def validate_annual_files(out_path, year_in, expected_times):
    """Discover available dates, rejecting ambiguous names and excluded products."""
    root = Path(out_path)
    allowed = {stamp.date() for stamp in expected_times}
    products = {}
    for path in sorted(root.glob(f'{year_in}*_fwi.nc')):
        if not path.is_file():
            continue
        stamp = annual_file_date(path)
        if stamp not in allowed:
            raise ValueError(f'Daily FWI product outside declared annual coverage: {path}')
        if stamp in products:
            raise ValueError(f'Ambiguous duplicate daily FWI date {stamp}: {products[stamp]} and {path}')
        products[stamp] = path
    if not products:
        raise ValueError(f'No available daily FWI products for {year_in} in {root}')
    return [str(products[stamp]) for stamp in sorted(products)]


def annual_file_date(path):
    match = re.fullmatch(r'(\d{4})_(\d{2})_(\d{2})(?:_LNOON)?_fwi\.nc', Path(path).name)
    if not match:
        raise ValueError(f'Malformed daily FWI filename: {path}')
    try:
        return date(*(int(part) for part in match.groups()))
    except ValueError as exc:
        raise ValueError(f'Invalid date in daily FWI filename: {path}') from exc


def validate_annual_daily_dataset(ds, date_parser=annual_file_date):
    """Check a file's single nominal local-date timestamp against its filename."""
    path = Path(ds.encoding.get('source', ''))
    stamp = date_parser(path)
    if 'fwi' not in ds or 'time' not in ds.fwi.dims:
        raise ValueError(f'Invalid daily FWI product: {path}')
    if ds.attrs.get('noon_convention') != NOON_CONVENTION:
        raise ValueError(f'Daily FWI noon convention differs; cannot mix 12 UTC and local-noon products: {path}')
    if ds.attrs.get('sampling') not in ('linear', 'nearest'):
        raise ValueError(f'Missing or invalid daily FWI sampling metadata: {path}')
    check_grid(ds.fwi, path)
    expected = pd.DatetimeIndex([pd.Timestamp(stamp) + pd.Timedelta(hours=12)])
    validate_annual_timestamps(ds.fwi, expected)
    return ds


def validate_annual_timestamps(fwi, expected_times):
    if 'time' not in fwi.coords or not np.array_equal(fwi.time.values, expected_times.values):
        raise ValueError('Annual daily timestamps must match every expected nominal local date exactly once at 12:00')


@contextmanager
def _open_annual_product_dataset(file_list, names, date_parser, validator,
                                 chunks=None, expected_sampling=None, allow_seasonal_nan=False):
    """Shared daily validation and mask continuity for FWI and its components."""
    from contextlib import ExitStack
    with ExitStack() as sources:
        datasets, paths = [], []
        land_mask = None
        threshold = None
        for path in sorted(file_list, key=date_parser):
            ds = sources.enter_context(open_product_dataset(path, engine='netcdf4', chunks=chunks))
            validator(ds)
            if expected_sampling is None:
                expected_sampling = ds.attrs['sampling']
            if ds.attrs['sampling'] != expected_sampling:
                raise ValueError(f'Annual FWI products use inconsistent sampling: {path}')
            if datasets and any(not np.array_equal(ds[coord].values, datasets[0][coord].values)
                                for coord in ('lat', 'lon')):
                raise ValueError(f'Annual FWI products use conflicting grids: {path}')
            if 'lsm' in ds:
                current_mask = normalize_land_mask(ds.lsm, ds.lat.values, ds.lon.values)
                current_threshold = float(ds.attrs.get('land_threshold', LAND_THRESHOLD))
                if not np.isfinite(current_threshold) or not 0 <= current_threshold <= 1:
                    raise ValueError(f'Invalid annual land-mask threshold: {path}')
                if land_mask is None:
                    land_mask = current_mask
                    threshold = current_threshold
                elif threshold != current_threshold or not np.array_equal(land_mask.values, current_mask.values):
                    raise ValueError(f'Annual FWI products use conflicting land masks or thresholds: {path}')
            elif land_mask is not None:
                raise ValueError(f'Annual daily product has lost the previously applied land mask: {path}')
            datasets.append(ds)
            paths.append(path)
        if not datasets:
            raise ValueError('No available daily FWI products to combine')
        active = None
        if land_mask is not None:
            active = xr.DataArray(active_land(land_mask, threshold), dims=('lat', 'lon'),
                                  coords={'lat': land_mask.lat, 'lon': land_mask.lon})
        checks = []
        for ds in datasets:
            field_checks = []
            for name in names:
                finite = (np.isfinite(ds[name]) | ds[name].isnull()) if allow_seasonal_nan else np.isfinite(ds[name])
                if active is not None:
                    finite = finite.where(active, True)
                field_checks.append(finite.all())
            checks.append(xr.concat(field_checks, dim='_annual_field'))
        # Dask executes reductions to only nfiles x nfields boolean values.
        valid = np.asarray(xr.concat(checks, dim='_annual_source').compute().values, dtype=bool)
        if not valid.all():
            row, column = np.argwhere(~valid)[0]
            raise ValueError(f'Non-finite annual {names[column].upper()} on an active grid cell: {paths[row]}')
        if land_mask is not None:
            prepared = []
            for ds in datasets:
                item = ds.copy(deep=False)
                item['lsm'] = land_mask
                for name in names:
                    item[name] = item[name].where(active)
                prepared.append(item)
            datasets = prepared
        with xr.combine_by_coords(datasets, data_vars='minimal', coords='minimal',
                                  compat='override', join='exact', combine_attrs='override') as combined:
            if land_mask is not None:
                combined['lsm'] = land_mask
                combined.attrs['land_threshold'] = threshold
            yield combined


@contextmanager
def open_annual_dataset(file_list, chunks=None, expected_sampling=None):
    """Open FWI days using the shared annual mask and validity checks."""
    with _open_annual_product_dataset(file_list, ('fwi',), annual_file_date,
                                     validate_annual_daily_dataset, chunks, expected_sampling) as ds:
        yield ds


def annual_component_file_date(path):
    match = re.fullmatch(r'(\d{4})_(\d{2})_(\d{2})_fwi_comp\.nc', Path(path).name)
    if not match:
        raise ValueError(f'Malformed daily FWI component filename: {path}')
    try:
        return date(*(int(part) for part in match.groups()))
    except ValueError as exc:
        raise ValueError(f'Invalid date in daily FWI component filename: {path}') from exc


def validate_annual_component_files(out_path, year_in, expected_times):
    """Find available component dates independently of the daily FWI files."""
    root = Path(out_path)
    allowed = {stamp.date() for stamp in expected_times}
    products = {}
    for path in sorted(root.glob(f'{year_in}*_fwi_comp.nc')):
        if not path.is_file():
            continue
        stamp = annual_component_file_date(path)
        if stamp not in allowed:
            raise ValueError(f'Daily FWI component product outside declared annual coverage: {path}')
        if stamp in products:
            raise ValueError(f'Ambiguous duplicate daily FWI component date {stamp}: {products[stamp]} and {path}')
        products[stamp] = path
    if not products:
        raise ValueError(f'No available daily FWI component products for {year_in} in {root}')
    return [str(products[stamp]) for stamp in sorted(products)]


def validate_annual_component_dataset(ds):
    path = Path(ds.encoding.get('source', ''))
    stamp = annual_component_file_date(path)
    names = ('FFMC', 'DMC', 'DC', 'ISI', 'BUI')
    missing = [name for name in names if name not in ds]
    if missing:
        raise ValueError(f'Missing daily FWI components {missing}: {path}')
    if ds.attrs.get('noon_convention') != NOON_CONVENTION:
        raise ValueError(f'Daily FWI component noon convention differs: {path}')
    if ds.attrs.get('sampling') not in ('linear', 'nearest'):
        raise ValueError(f'Missing or invalid daily FWI component sampling metadata: {path}')
    expected = pd.DatetimeIndex([pd.Timestamp(stamp) + pd.Timedelta(hours=12)])
    for name in names:
        if set(ds[name].dims) != {'time', 'lat', 'lon'}:
            raise ValueError(f'Invalid daily {name} dimensions: {path}')
        check_grid(ds[name], path)
        validate_annual_timestamps(ds[name], expected)
    return ds


@contextmanager
def open_annual_component_dataset(file_list, chunks=None, expected_sampling=None):
    names = ('FFMC', 'DMC', 'DC', 'ISI', 'BUI')
    with _open_annual_product_dataset(file_list, names, annual_component_file_date,
                                     validate_annual_component_dataset, chunks, expected_sampling) as ds:
        yield ds


def reduce_component_statistics(components_year, metadata, quantile_input=None, land_mask=None):
    """Mean/min/max and 50th/95th/99th percentiles for all five FWI components.

    Components have different scales, so the FWI danger-category thresholds do
    not apply to them. All statistics use only the available daily products.
    """
    names = ('FFMC', 'DMC', 'DC', 'ISI', 'BUI')
    values = components_year[list(names)].astype('float32')
    source = values if quantile_input is None else quantile_input[list(names)].astype('float32')
    active = None
    if land_mask is not None:
        land_mask = normalize_land_mask(land_mask, values.lat.values, values.lon.values)
        threshold = float(metadata.get('land_threshold', LAND_THRESHOLD))
        if not np.isfinite(threshold) or not 0 <= threshold <= 1:
            raise ValueError('Invalid annual component land-mask threshold')
        active = xr.DataArray(active_land(land_mask, threshold), dims=('lat', 'lon'),
                              coords={'lat': land_mask.lat, 'lon': land_mask.lon})
        # Avoid all-NaN quantile warnings; mask every final result below.
        values = values.where(active, 0.0)
        source = source.where(active, 0.0)
    product = xr.Dataset()
    for name in names:
        for statistic in ('mean', 'min', 'max'):
            key = f'{name}_{statistic}'
            product[key] = getattr(values[name], statistic)('time').astype('float32')
            product[key].attrs = {'units': '1', 'long_name': f'Annual {statistic} of {name} on available days'}
        key = f'{name}_quantile'
        product[key] = source[name].quantile([0.50, 0.95, 0.99], dim='time').astype('float32')
        product[key].attrs = {'units': '1', 'long_name': f'Annual quantiles of {name} on available days'}
    if active is not None:
        product = product.where(active)
        product['lsm'] = land_mask.copy(deep=False)
    product.attrs.update({key: value for key, value in metadata.items()
                          if key not in {'local_date', 'state_date', 'input_files', 'history'}})
    add_history_attr(product)
    return product


def write_annual_component_product(out_path, year_in, product):
    """Write all component statistics together, preserving the original globals."""
    path = Path(out_path) / f'{year_in}_fwi_comp_annual.nc'
    encoding = {}
    for name, da in product.data_vars.items():
        entry = {'zlib': True, 'complevel': 4, 'dtype': 'float64' if name == 'lsm' else 'float32'}
        if da.chunks is not None:
            entry['chunksizes'] = tuple(da.chunksizes[dim][0] for dim in da.dims)
        encoding[name] = entry
    atomic_write(product, path, format='NETCDF4', engine='netcdf4', encoding=encoding, compute=True)
    return path


def annual_component_statistics(out_path, year_in, FWI_all, coverage_end=None):
    """Write one annual file for available component days, including partial years.

    The caller controls the independent annual-component switch. Direct calls
    calculate statistics regardless of the script's default switch value.
    """
    allowed_times, _ = annual_coverage(year_in, coverage_end)
    file_list = validate_annual_component_files(out_path, year_in, allowed_times)
    expected_times, coverage = annual_coverage(
        year_in, coverage_end, available_dates=[annual_component_file_date(path) for path in file_list])
    import dask

    label = 'partial-year' if coverage['is_partial_year'] else 'full-year'
    print(f'Computing {label} component statistics for {year_in}: '
          f'{coverage["days_used"]} available days, {coverage["missing_days"]} missing calendar days', flush=True)
    names = ('FFMC', 'DMC', 'DC', 'ISI', 'BUI')
    with dask.config.set({'array.slicing.split_large_chunks': True}):
        with open_annual_component_dataset(file_list, chunks={'time': 15},
                                           expected_sampling=FWI_all.attrs.get('sampling')) as ds_year:
            values = ds_year[list(names)]
            validate_annual_timestamps(values.FFMC, expected_times)
            spatial_dims = [dim for dim in values.FFMC.dims if dim != 'time']
            values = values.chunk({'time': 15, **{dim: 200 for dim in spatial_dims}}).astype('float32')
            quantile_values = values.chunk({'time': -1, **{dim: 100 for dim in spatial_dims}})
            metadata = dict(ds_year.attrs)
            metadata.update(FWI_all.attrs)
            metadata.update(coverage)
            metadata.update(noon_convention=NOON_CONVENTION, sampling=ds_year.attrs['sampling'])
            if 'lsm' in ds_year:
                metadata['land_threshold'] = ds_year.attrs['land_threshold']
            product = reduce_component_statistics(values, metadata, quantile_values, ds_year.get('lsm'))
            write_annual_component_product(out_path, year_in, product)
    print(f'{label.capitalize()} component statistics written OK.', flush=True)



def reduce_annual_statistics(fwi_year, metadata, quantile_input=None, land_mask=None):
    """Observed-day reductions; masked ocean cells remain missing in all products."""
    fwi_year = fwi_year.astype('float32')
    source = fwi_year if quantile_input is None else quantile_input.astype('float32')
    active = None
    if land_mask is not None:
        land_mask = normalize_land_mask(land_mask, fwi_year.lat.values, fwi_year.lon.values)
        threshold = float(metadata.get('land_threshold', LAND_THRESHOLD))
        if not np.isfinite(threshold) or not 0 <= threshold <= 1:
            raise ValueError('Invalid annual land-mask threshold')
        active = xr.DataArray(active_land(land_mask, threshold), dims=('lat', 'lon'),
                              coords={'lat': land_mask.lat, 'lon': land_mask.lon})
        # Ocean values are excluded from every published reduction. Temporary
        # zeros avoid all-NaN quantile warnings; the final ocean results are NaN.
        fwi_year = fwi_year.where(active, 0.0)
        source = source.where(active, 0.0)
    annual_basic = xr.Dataset({
        'fwi_mean': fwi_year.mean('time'),
        'fwi_max': fwi_year.max('time'),
        'fwi_min': fwi_year.min('time'),
        'days_VLow': (fwi_year < 5.2).sum('time').astype('int16'),
        'days_Low': ((fwi_year >= 5.2) & (fwi_year < 11.2)).sum('time').astype('int16'),
        'days_Mod': ((fwi_year >= 11.2) & (fwi_year < 21.3)).sum('time').astype('int16'),
        'days_High': ((fwi_year >= 21.3) & (fwi_year < 38.0)).sum('time').astype('int16'),
        'days_VHigh': ((fwi_year >= 38.0) & (fwi_year < 50.0)).sum('time').astype('int16'),
        'days_Ext': (fwi_year >= 50.0).sum('time').astype('int16'),
    })
    quantiles = source.quantile([0.50, 0.95, 0.99], dim='time').astype('float32').to_dataset(name='fwi')
    if active is not None:
        annual_basic = annual_basic.where(active)
        quantiles = quantiles.where(active)
        for product in (annual_basic, quantiles):
            product['lsm'] = land_mask.copy(deep=False)
    attrs = {key: value for key, value in metadata.items()
             if key not in {'local_date', 'state_date', 'input_files', 'history'}}
    for product in (annual_basic, quantiles):
        product.attrs.update(attrs)
        add_history_attr(product)
    return annual_basic, quantiles


def write_annual_products(out_path, year_in, basic, quantiles, suffix=""):
    """Use the requested annual names and preserve missing ocean day counts."""
    root = Path(out_path)
    basic_path = root / f'{year_in}_fwi{suffix}_basic_and_days.nc'
    quantile_path = root / f'{year_in}_fwi{suffix}_quantiles.nc'

    def encoding(name, da):
        is_count = name.startswith('days_')
        dtype = 'int16' if is_count else ('float64' if name == 'lsm' else 'float32')
        enc = {'zlib': True, 'complevel': 4, 'dtype': dtype}
        if is_count:
            enc['_FillValue'] = np.int16(-32767)
        if da.chunks is not None:
            enc['chunksizes'] = tuple(da.chunksizes[dim][0] for dim in da.dims)
        return enc

    for path, product in ((basic_path, basic), (quantile_path, quantiles)):
        atomic_write(product, path, format='NETCDF4', engine='netcdf4',
                     encoding={name: encoding(name, da) for name, da in product.data_vars.items()}, compute=True)
    return basic_path, quantile_path



if __name__ == "__main__":
    main()
