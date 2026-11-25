#!/usr/bin/env python3
# -----------------------------------------------------------------------------------------------------------------------------
#   Copyright 2026 Didier M. Roche (a.k.a. dmr) | iLOVECLIM / FRATRES coding group
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
# -----------------------------------------------------------------------------------------------------------------------------
"""Make the initial land isotope map of the LBM (land_init_mode = 'file', wisoatm.nml; ROADMAP D69, wiso:3a:0007).

Nothing fractionates on land (D64): a well-mixed snow layer at steady state has the composition of the snowfall that fed
it, weighted by the snowfall amount. This script takes that composition from the yearly ECBilt output of a run (WISO = 1)
and writes it on the model grid, so that a cold start begins with the snow (and the soil water) near equilibrium instead
of at uniform namelist values (snow over the ice sheets: decades to centuries of spin-up otherwise, O57).

Input : one or more ECBilt yearly output files (e.g. outputdata/atmos/atmym*.nc), with snow, d18s, d17s, dds (snowfall)
        and pp, d18p, d17p, ddp (precipitation), in permil; time records = years.
Output: netCDF with lat, lon (as in the input) and d18o_snow, dexc_snow, o17exc_snow, d18o_soil, dexc_soil, o17exc_soil
        (permil, per meg), fill -1e30 where the amount is below --min-amount (the model then uses the namelist values).

Averaging: deltas and d-excess are linear in the ratios, so their amount-weighted means over the years are exact; the
17O-excess is not, so it is recomputed from the averaged d17O and d18O: ln(1 + d17) - 0.528 ln(1 + d18).

Usage: make_wisolbm_init.py OUT.nc IN.nc [IN2.nc ...] [--last N] [--min-amount A]
"""
import argparse
import sys

import numpy as np
import netCDF4 as nc

FILL_IN = 1.0e19        # output fill values of ECBilt are 1e20
FILL_OUT = -1.0e30      # missing for the model reader (|value| > 1e10)
THETA = 0.528


def read(files, names, last):
    """Concatenate the variables over the time records of the files; keep the last N records."""
    out = {n: [] for n in names}
    lat = lon = None
    for f in files:
        with nc.Dataset(f) as d:
            if lat is None:
                lat, lon = d['lat'][:].data.copy(), d['lon'][:].data.copy()
            for n in names:
                if n not in d.variables:
                    sys.exit(f'{f}: variable {n} not found (yearly output of a WISO = 1 run needed)')
                x = d[n][:]
                x = np.ma.masked_invalid(np.ma.masked_where(np.abs(np.ma.filled(x, 0.0)) > FILL_IN, x))
                out[n].append(x if x.ndim == 3 else x[None])
    out = {n: np.ma.concatenate(v, axis=0) for n, v in out.items()}
    if last:
        out = {n: v[-last:] for n, v in out.items()}
    return lat, lon, out


def mean_composition(amount, d18, d17, dd, min_amount):
    """Amount-weighted means over the years (permil); 17O-excess from the averaged deltas (per meg)."""
    w = np.ma.masked_where(np.ma.getmaskarray(d18) | np.ma.getmaskarray(dd) | np.ma.getmaskarray(d17), amount)
    w = np.ma.filled(np.ma.maximum(w, 0.0), 0.0)
    wsum = w.sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        m = {k: (w*np.ma.filled(v, 0.0)).sum(axis=0)/wsum for k, v in (('d18', d18), ('d17', d17), ('dd', dd))}
    ok = wsum/amount.shape[0] >= min_amount
    dexc = m['dd'] - 8.0*m['d18']
    o17 = (np.log1p(m['d17']*1e-3) - THETA*np.log1p(m['d18']*1e-3))*1e6
    return [np.where(ok, x, FILL_OUT) for x in (m['d18'], dexc, o17)], ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('inputs', nargs='+')
    ap.add_argument('--last', type=int, default=0, help='use the last N years only (default: all)')
    ap.add_argument('--min-amount', type=float, default=0.1,
                    help='mean amount (cm/yr, as the output) below which a cell is missing (default 0.1)')
    a = ap.parse_args()

    lat, lon, v = read(a.inputs, ['snow', 'd18s', 'd17s', 'dds', 'pp', 'd18p', 'd17p', 'ddp'], a.last)
    snow, oks = mean_composition(v['snow'], v['d18s'], v['d17s'], v['dds'], a.min_amount)
    soil, okp = mean_composition(v['pp'], v['d18p'], v['d17p'], v['ddp'], a.min_amount)

    with nc.Dataset(a.out, 'w') as d:
        d.title = 'Initial land water isotopes of the LBM (land_init_mode = file, wiso:3a:0007)'
        d.source = 'make_wisolbm_init.py: amount-weighted snowfall (snow) and precipitation (soil) composition of ' + \
                   ', '.join(a.inputs) + (f', last {a.last} years' if a.last else f', {v["snow"].shape[0]} years')
        d.min_amount_cm_per_yr = a.min_amount
        d.createDimension('lat', len(lat))
        d.createDimension('lon', len(lon))
        x = d.createVariable('lat', 'f8', ('lat',)); x[:] = lat; x.units = 'degrees_north'
        x = d.createVariable('lon', 'f8', ('lon',)); x[:] = lon; x.units = 'degrees_east'
        for (name, units, long), data in zip(
                [('d18o_snow', 'permil', 'snow delta 18O'), ('dexc_snow', 'permil', 'snow d-excess'),
                 ('o17exc_snow', 'per meg', 'snow 17O-excess'), ('d18o_soil', 'permil', 'soil water delta 18O'),
                 ('dexc_soil', 'permil', 'soil water d-excess'), ('o17exc_soil', 'per meg', 'soil water 17O-excess')],
                snow + soil):
            x = d.createVariable(name, 'f8', ('lat', 'lon'), fill_value=FILL_OUT)
            x[:] = data; x.units = units; x.long_name = long + ' (amount-weighted mean)'

    print(f'{a.out}: snow defined in {int(oks.sum())} cells, soil in {int(okp.sum())} of {oks.size}; '
          f'snow d18O {np.min(snow[0][oks]):.1f} to {np.max(snow[0][oks]):.1f} permil')


if __name__ == '__main__':
    main()
