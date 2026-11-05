"""
FWI run script for ClimateDT 

"""
import argparse
from datetime import date
from pathlib import Path

from fwi.run_fwi import run_fwi
from fwi.run_fwi import (
    DO_ANNUAL_COMPONENTS,
    DO_ANNUAL_OVERWINTERING,
    DO_OVERWINTERING,
    FINAL_LOCAL_DATE,
    IO_BUFFER_MB,
    IO_TILE_COLUMNS,
    IO_TILE_ROWS,
    LAND_THRESHOLD,
    LSM_FILE,
    MONTHLY_CHUNK,
    PRECIP_TIME_LABEL,
    WRITE_COMPONENTS,
)


def _get_parser():
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
    return parser


def main(argv=None):
    args = _get_parser().parse_args(argv)
    return run_fwi(args)


if __name__ == "__main__":
    main()
