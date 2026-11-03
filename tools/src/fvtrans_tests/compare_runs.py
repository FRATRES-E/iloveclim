#!/usr/bin/env python3
"""Compare the moisture climate of ECBilt runs (ROADMAP step 1.3, wiso:1e:0004).

Usage: python3 compare_runs.py [--report-figure] [--mask <fractoc.dat>] <outdir> <label>=<rundir> [<label>=<rundir> ...]
  e.g. python3 compare_runs.py --report-figure cmp spectral=run_spec fv=run_fv nohdiff=run_spec_tdifq

The first run is the reference. In each <rundir> (searched recursively) the script reads
  atmmmyl*.nc    ECBilt monthly means (one or several files, taken in name order, concatenated in time);
  wbudget*.txt   the water budget of atm_wbudget_mod (optional, wiso:1e:0003).
Only complete years are used, and the runs are cut to the shortest one.

Writes <outdir>/compare_report.md and <outdir>/figures/*.png:
  - water budget per run (global means, mm/yr): E, P (large-scale, convective), clipping C1/C2/C3, transport and
    diffusion integrals, residual;
  - per variable (pp, lsp, cp, snow, evap, eminp, q when present): global mean and interannual spread, difference to
    the reference, area with a significant difference (Welch t-test on annual means, 5 %), zonal means, difference maps;
  - monthly precipitation distribution (area-weighted percentiles) and regional means of pp, evap, eminp, q with a
    significance mark (polar caps, dry regions, tropics, Sahel box SAHEL);
  - clipping maps (user_t1 levels 1-3 = C1, C2, C3 in cm/yr, wiso:1e:0003) as zonal means.
With --report-figure, also <outdir>/figures/report_differences.png: one row per variable (snow, pp, q), one column per
run other than the reference, differences to the reference only, cells not significant at 5 % greyed, one colour scale
per row (symmetric, at the 99th percentile of the absolute differences, colour bar extended beyond).
Column water q is converted from m (model output) to mm.
Coastlines: every map shows the model's own coastline, the 0.5 contour of the ocean fraction of ECBilt's fractoc.dat
(T21, as the model sees it). The file is taken from --mask, otherwise found under the reference run directory
(e.g. inputdata/fractoc.dat); without it the maps are drawn without coastlines.
Note: ECBilt monthly means are built from one 4 h step per day, so monthly precipitation extremes are noisy;
multi-year means are the robust quantity. The budget file integrates every step.
Needs numpy, matplotlib, netCDF4 (scipy optional, for the exact t quantile).
"""

import glob
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4  # noqa: E402
import numpy as np  # noqa: E402

import nicopal as ncp

try:
    from scipy import stats as _stats
except ImportError:  # pragma: no cover
    _stats = None

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlesize": 10, "axes.titlecolor": INK, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "legend.frameon": False})

VARS = [("pp", "total precipitation", "cm/yr"), ("lsp", "large-scale precipitation", "cm/yr"),
        ("cp", "convective precipitation", "cm/yr"), ("snow", "snowfall", "cm/yr"),
        ("evap", "evaporation", "cm/yr"), ("eminp", "evaporation minus precipitation", "cm/yr"),
        ("q", "column water", "mm")]
SCALE = {"q": 1000.0}  # model output units -> report units (q: m -> mm)
REPORT_VARS = ["snow", "pp", "q"]  # rows of --report-figure
CLIP = ["C1 (clipping in ec_moisture)", "C2 (removed by ec_moisbalance)", "C3 (clipping in ec_convec_two)"]
BUDGET_COLS = ["year", "month", "nstep", "W_end", "dW", "E", "Pls", "Pcv", "T", "D", "C1", "C2", "C3", "R", "absR"]
DRY_CMYR = 25.0      # dry region: reference climatological precipitation below this (cm/yr)
POLAR_LAT = 66.5     # polar caps
SAHEL = (-20.0, 40.0, 5.0, 20.0)  # West Africa / Sahel box: lon west, lon east (degrees east), lat south, lat north
PCTS = [1, 5, 25, 50, 75, 95, 99]

CMAP_SEQ = ncp.pal("Uranium_r")
CMAP_DIV = ncp.pal("Sulfur")
CMAP_DIV_R = ncp.pal("Sulfur_r")
CMAP_DIV_ALT = ncp.pal("Vanadium")
CMAP_DIV_ALT_R = ncp.pal("Vanadium_r")

def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def find(rundir, pattern):
    return sorted(glob.glob(os.path.join(rundir, "**", pattern), recursive=True))


def read_run(rundir):
    """Monthly fields (time, lat, lon) per variable, user_t1 (time, 3, lat, lon), coordinates, budget table."""
    files = find(rundir, "atmmmyl*.nc")
    if not files:
        sys.exit(f"compare_runs: no atmmmyl*.nc under {rundir}")
    fields, lat, lon = {}, None, None
    for f in files:
        with netCDF4.Dataset(f) as d:
            if lat is None:
                lat, lon = np.array(d["lat"][:], float), np.array(d["lon"][:], float)
            for name in [v[0] for v in VARS] + ["user_t1"]:
                if name in d.variables:
                    a = d[name][:]
                    a = np.ma.filled(a.astype(float), np.nan)*SCALE.get(name, 1.0)
                    fields.setdefault(name, []).append(a)
    fields = {k: np.concatenate(v, axis=0) for k, v in fields.items()}
    budget = None
    bfiles = find(rundir, "wbudget*.txt")
    if bfiles:
        rows = []
        for f in bfiles:
            rows += [list(map(float, ln.split())) for ln in open(f) if ln.strip() and not ln.startswith("#")]
        if rows:
            budget = dict(zip(BUDGET_COLS, np.array(rows).T))
    return dict(files=files, bfiles=bfiles, fields=fields, lat=lat, lon=lon, budget=budget)


def lat_weights(lat):
    """Gaussian weights (sum 1) matched to the model latitudes; cos(lat) if they are not Gaussian."""
    x, w = np.polynomial.legendre.leggauss(len(lat))
    mu = np.sin(np.deg2rad(lat))
    if np.allclose(np.sort(mu), x, atol=1e-4):
        w = w[np.argsort(np.argsort(mu))]
    else:
        w = np.cos(np.deg2rad(lat))
    return w/w.sum()


def gmean(f, w):
    """Area mean over the last two axes (lat, lon)."""
    return np.nansum(np.nanmean(f, axis=-1)*w, axis=-1)


def annual(f, nyr):
    return f[:nyr*12].reshape((nyr, 12) + f.shape[1:]).mean(axis=1)


def tcrit(dof):
    if _stats is not None:
        return _stats.t.ppf(0.975, dof)
    return 2.0


def significant(a, b):
    """Welch t-test on annual means a, b (years, ...): mask where the difference is significant at 5 %."""
    na, nb = a.shape[0], b.shape[0]
    va, vb = a.var(axis=0, ddof=1)/na, b.var(axis=0, ddof=1)/nb
    se = np.sqrt(va + vb)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (b.mean(0) - a.mean(0))/se
        dof = (va + vb)**2/(va**2/(na - 1) + vb**2/(nb - 1))
    dof = np.where(np.isfinite(dof), dof, 1.0)
    return np.abs(t) > tcrit(np.maximum(dof, 1.0))


def edges(c, period=None):
    e = np.concatenate([[c[0] - (c[1] - c[0])/2], (c[1:] + c[:-1])/2, [c[-1] + (c[-1] - c[-2])/2]])
    if period is None:
        e[0], e[-1] = max(e[0], -90.0), min(e[-1], 90.0)
    return e


def wpercentiles(x, w, pcts):
    """Weighted percentiles of x (flattened) with weights w."""
    x, w = x.ravel(), w.ravel()
    ok = np.isfinite(x)
    x, w = x[ok], w[ok]
    i = np.argsort(x)
    cw = np.cumsum(w[i])
    cw /= cw[-1]
    return [x[i][np.searchsorted(cw, p/100.0)] for p in pcts]


COAST = None  # (lon, lat, ocean fraction) with a periodic column, set by read_coast


def pm180(lon):
    """Index shift and longitudes for maps centred on 0 (-180..180): roll by -k, then lon_plot."""
    k = int(np.argmax(lon >= 180.0)) if np.any(lon >= 180.0) else 0
    return k, np.roll(np.where(lon >= 180.0, lon - 360.0, lon), -k)


def mapmesh(ax, lat, lon, f, **kw):
    """pcolormesh of a (lat, lon) field on a map centred on longitude 0 (-180..180), with axis ticks."""
    k, lp = pm180(lon)
    m = ax.pcolormesh(edges(lp, 360.0), edges(lat), np.roll(f, -k, axis=-1), shading="flat", **kw)
    ax.set_xlim(-180.0, 180.0)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    return m


def read_coast(path, lat, lon):
    """Ocean fraction from ECBilt's fractoc.dat: one header line, nlat*nlon values, longitude fastest, row 1 south."""
    global COAST
    try:
        vals = np.array(open(path).read().split("\n", 1)[1].split(), float)
        frac = vals[:len(lat)*len(lon)].reshape(len(lat), len(lon))
    except (OSError, IndexError, ValueError):
        print("compare_runs: could not read the land-sea mask", path, "- maps without coastlines")
        return
    k, lp = pm180(lon)
    frac = np.roll(frac, -k, axis=1)
    COAST = (np.append(lp, lp[0] + 360.0), lat, np.concatenate([frac, frac[:, :1]], axis=1))
    print("compare_runs: coastlines from", path)


def add_coast(ax):
    """Model coastline (ocean fraction 0.5) on a longitude-latitude map."""
    if COAST is not None:
        ax.contour(COAST[0], COAST[1], COAST[2], levels=[0.5], colors=INK, linewidths=0.6)


def report_figure(figdir, runs, labels, ref, nyr, lat, lon):
    """Differences to the reference for REPORT_VARS: rows = variables, columns = runs other than the reference."""
    others = [lab for lab in labels if lab != ref]
    long = {v[0]: (v[1], v[2]) for v in VARS}
    rows = [v for v in REPORT_VARS if all(v in runs[lab]["fields"] for lab in labels)]
    if not others or not rows:
        print("compare_runs: --report-figure needs a second run and at least one of", ", ".join(REPORT_VARS))
        return
    fig, axs = plt.subplots(len(rows), len(others), figsize=(5.0*len(others) + 1.2, 2.9*len(rows)), squeeze=False,
                            constrained_layout=True)
    for r, name in enumerate(rows):
        ann = {lab: annual(runs[lab]["fields"][name], nyr) for lab in labels}
        diffs = [ann[lab].mean(0) - ann[ref].mean(0) for lab in others]
        dmax = np.nanpercentile(np.abs(np.concatenate([d.ravel() for d in diffs])), 99) or 1.0
        for c, (lab, d) in enumerate(zip(others, diffs)):
            a = axs[r, c]
            sig = significant(ann[ref], ann[lab]) if nyr > 1 else np.ones_like(d, bool)
            m = mapmesh(a, lat, lon, np.where(sig, d, np.nan), cmap=CMAP_DIV if name == "q" else CMAP_DIV_R,
                             vmin=-dmax, vmax=dmax)
            mapmesh(a, lat, lon, np.where(sig, np.nan, 0.0), cmap="Greys", vmin=-1, vmax=4)
            add_coast(a)
            a.set_title(f"{long[name][0]}: {lab} - {ref}", loc="left")
            a.set_ylabel("latitude" if c == 0 else "")
            a.set_xlabel("longitude" if r == len(rows) - 1 else "")
        fig.colorbar(m, ax=axs[r, :], shrink=0.9, extend="both", label=long[name][1])
    fig.savefig(os.path.join(figdir, "report_differences.png"), dpi=200)
    plt.close(fig)


def main():
    argv = sys.argv[1:]
    want_report_fig = "--report-figure" in argv
    argv = [a for a in argv if a != "--report-figure"]
    mask = None
    if "--mask" in argv:
        k = argv.index("--mask")
        if k + 1 >= len(argv):
            sys.exit(__doc__)
        mask = argv[k + 1]
        del argv[k:k + 2]
    args = argv
    if len(args) < 2 or any("=" not in a for a in args[1:]):
        sys.exit(__doc__)
    outdir = args[0]
    figdir = os.path.join(outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    labels = [a.split("=", 1)[0] for a in args[1:]]
    runs = {lab: read_run(a.split("=", 1)[1]) for lab, a in zip(labels, args[1:])}
    ref = labels[0]
    lat, lon = runs[ref]["lat"], runs[ref]["lon"]
    if mask is None:
        found = find(args[1].split("=", 1)[1], "fractoc.dat")
        mask = found[0] if found else None
    if mask is not None:
        read_coast(mask, lat, lon)
    w = lat_weights(lat)
    w2 = np.repeat(w[:, None], len(lon), axis=1)/len(lon)
    color = dict(zip(labels, SERIES))

    nmon = min(min(f.shape[0] for f in r["fields"].values()) for r in runs.values())
    nyr = nmon//12
    if nyr < 1:
        sys.exit("compare_runs: less than one complete year of monthly means")
    lines = ["# ECBilt moisture: run comparison", "",
             f"Reference: **{ref}**. {nyr} complete years per run (shortest run). Monthly means are one 4 h step per "
             "day; the budget integrates every step.", "", "| run | netCDF files | budget files |", "|---|---|---|"]
    for lab in labels:
        r = runs[lab]
        lines.append(f"| {lab} | {len(r['files'])} | {len(r['bfiles'])} |")

    # 1. water budget
    if any(runs[lab]["budget"] is not None for lab in labels):
        lines += ["", "## Water budget (global means, mm/yr; W in mm)", "",
                  "| term | " + " | ".join(labels) + " |", "|---|" + "---|"*len(labels)]
        terms = [("W_end", "column water W (mean)"), ("E", "evaporation E"), ("Pls", "large-scale precipitation"),
                 ("Pcv", "convective precipitation"), ("T", "transport integral T"), ("D", "diffusion integral D"),
                 ("C1", "C1 clipping in ec_moisture"), ("C2", "C2 removed by ec_moisbalance"),
                 ("C1-C2", "C1 - C2 (net water created)"), ("C3", "C3 clipping in convection"),
                 ("R", "residual R"), ("absR", "sum of per-step abs residuals")]
        for key, name in terms:
            row = []
            for lab in labels:
                b = runs[lab]["budget"]
                if b is None:
                    row.append("-")
                    continue
                n = min(len(b["E"]), nyr*12)
                if key == "W_end":
                    v = b["W_end"][:n].mean()
                elif key == "C1-C2":
                    v = (b["C1"][:n] - b["C2"][:n]).sum()/(n/12)
                else:
                    v = b[key][:n].sum()/(n/12)
                row.append(f"{v:.4g}")
            lines.append(f"| {name} | " + " | ".join(row) + " |")

    # 2. variables: global means, significance, zonal means, maps
    common = [v for v in VARS if all(v[0] in runs[lab]["fields"] for lab in labels)]
    if common:
        lines += ["", "## Climatology and differences to the reference", "",
                  "| variable | run | global mean | interannual sd | diff to ref | pattern rms diff | "
                  "area with significant diff |", "|---|---|---|---|---|---|---|"]
    for name, long, unit in common:
        ann = {lab: annual(runs[lab]["fields"][name], nyr) for lab in labels}
        clim = {lab: ann[lab].mean(0) for lab in labels}
        gm = {lab: gmean(ann[lab], w) for lab in labels}
        for lab in labels:
            if lab == ref:
                sd = gm[lab].std(ddof=1) if nyr > 1 else 0.0
                lines.append(f"| {name} ({unit}) | {lab} | {gm[lab].mean():.4g} | {sd:.2g} | - | - | - |")
                continue
            d = clim[lab] - clim[ref]
            rms = np.sqrt(np.nansum(w2*d**2))
            sig = significant(ann[ref], ann[lab]) if nyr > 1 else np.zeros_like(d, bool)
            lines.append(f"| {name} ({unit}) | {lab} | {gm[lab].mean():.4g} | "
                         f"{gm[lab].std(ddof=1) if nyr > 1 else 0:.2g} | {gm[lab].mean() - gm[ref].mean():+.3g} | "
                         f"{rms:.3g} | {100*np.sum(w2*sig):.0f} % |")

        # zonal means and their differences (two panels, one quantity each)
        fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), sharex=True)
        for lab in labels:
            zm = np.nanmean(clim[lab], axis=1)
            ax[0].plot(lat, zm, "-", color=color[lab], lw=2, label=lab)
            if lab != ref:
                ax[1].plot(lat, zm - np.nanmean(clim[ref], axis=1), "-", color=color[lab], lw=2, label=f"{lab} - {ref}")
        ax[1].axhline(0.0, color=INK2, lw=1)
        ax[0].set_title(f"{long}, zonal mean ({unit})", loc="left")
        ax[1].set_title(f"difference to {ref} ({unit})", loc="left")
        for a in ax:
            a.set_xlabel("latitude (degrees)")
            a.legend()
            style(a)
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, f"zonal_{name}.png"), dpi=150)
        plt.close(fig)

        # difference maps, non-significant cells greyed
        others = [lab for lab in labels if lab != ref]
        if others:
            fig, axs = plt.subplots(1, len(others) + 1, figsize=(5.2*(len(others) + 1), 3.2), squeeze=False)
            axs = axs[0]
            m = mapmesh(axs[0], lat, lon, clim[ref], cmap=CMAP_SEQ if name != "eminp" else CMAP_DIV_ALT)
            fig.colorbar(m, ax=axs[0], shrink=0.85, label=unit)
            axs[0].set_title(f"{name}: {ref}", loc="left")
            dmax = max(np.nanmax(np.abs(clim[lab] - clim[ref])) for lab in others) or 1.0
            for a, lab in zip(axs[1:], others):
                d = clim[lab] - clim[ref]
                sig = significant(ann[ref], ann[lab]) if nyr > 1 else np.ones_like(d, bool)
                m = mapmesh(a, lat, lon, np.where(sig, d, np.nan), cmap=CMAP_DIV_R, # if name != "eminp" else CMAP_DIV,
                                 vmin=-dmax, vmax=dmax)
                mapmesh(a, lat, lon, np.where(sig, np.nan, 0.0), cmap="Greys", vmin=-1, vmax=4)
                fig.colorbar(m, ax=a, shrink=0.85, label=unit)
                a.set_title(f"{lab} - {ref} (grey: not significant)", loc="left")
            for a in axs:
                add_coast(a)
                a.set_xlabel("longitude")
                a.set_ylabel("latitude")
            fig.tight_layout()
            fig.savefig(os.path.join(figdir, f"map_{name}.png"), dpi=150)
            plt.close(fig)

    # 3. monthly precipitation distribution, regions
    if all("pp" in runs[lab]["fields"] for lab in labels):
        lines += ["", "## Monthly precipitation distribution (area-weighted percentiles, cm/yr)", "",
                  "| run | " + " | ".join(f"p{p}" for p in PCTS) + " |", "|---|" + "---|"*len(PCTS)]
        fig, ax = plt.subplots(figsize=(6, 3.6))
        for lab in labels:
            pp = runs[lab]["fields"]["pp"][:nyr*12]
            wt = np.broadcast_to(w2, pp.shape)
            pv = wpercentiles(pp, wt, PCTS)
            lines.append(f"| {lab} | " + " | ".join(f"{v:.3g}" for v in pv) + " |")
            qs = np.linspace(0.5, 99.5, 199)
            ax.plot(qs, wpercentiles(pp, wt, qs), "-", color=color[lab], lw=2, label=lab)
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_xlabel("percentile of area-weighted monthly values")
        ax.set_ylabel("precipitation (cm/yr)")
        ax.set_title("Monthly precipitation, quantile function", loc="left")
        ax.legend()
        style(ax)
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, "pp_quantiles.png"), dpi=150)
        plt.close(fig)

        clim_ref = annual(runs[ref]["fields"]["pp"], nyr).mean(0)
        latg = np.repeat(lat[:, None], len(lon), axis=1)
        long_ = np.repeat(lon[None, :], len(lat), axis=0)
        lon180 = np.where(long_ > 180.0, long_ - 360.0, long_)
        regions = [("south polar cap (< -66.5)", latg < -POLAR_LAT), ("north polar cap (> 66.5)", latg > POLAR_LAT),
                   (f"dry ({ref} pp < {DRY_CMYR:g} cm/yr)", clim_ref < DRY_CMYR),
                   ("tropics (|lat| < 30)", np.abs(latg) < 30.0),
                   (f"Sahel ({SAHEL[0]:g} to {SAHEL[1]:g} E, {SAHEL[2]:g} to {SAHEL[3]:g} N)",
                    (lon180 >= SAHEL[0]) & (lon180 <= SAHEL[1]) & (latg >= SAHEL[2]) & (latg <= SAHEL[3]))]
        lines += ["", "## Regional means (time mean; * = difference to the reference significant at 5 %)", "",
                  "| region | area | variable | " + " | ".join(labels) + " |", "|---|---|---|" + "---|"*len(labels)]
        for rname, mask in regions:
            area = np.sum(w2*mask)
            if area == 0:
                continue
            for name in ("pp", "evap", "eminp", "q"):
                if not all(name in runs[lab]["fields"] for lab in labels):
                    continue
                # annual regional means (years), then time mean and Welch test against the reference
                series = {lab: np.nansum(w2*mask*annual(runs[lab]["fields"][name], nyr), axis=(1, 2))/area
                          for lab in labels}
                cells = []
                for lab in labels:
                    mark = ""
                    if lab != ref and nyr > 1 and bool(significant(series[ref], series[lab])):
                        mark = "*"
                    cells.append(f"{series[lab].mean():.4g}{mark}")
                lines.append(f"| {rname} | {100*area:.1f} % | {name} | " + " | ".join(cells) + " |")

    # 4. clipping maps (user_t1)
    have = [lab for lab in labels if "user_t1" in runs[lab]["fields"]]
    if have:
        lines += ["", "## Clipping (user_t1, cm/yr, time mean)", "",
                  "| run | field | global mean | max | latitude of the zonal maximum |", "|---|---|---|---|---|"]
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.4), sharex=True)
        for lab in have:
            u = np.nanmean(runs[lab]["fields"]["user_t1"][:nyr*12], axis=0)
            for k in range(3):
                zm = np.nanmean(u[k], axis=1)
                lines.append(f"| {lab} | {CLIP[k]} | {np.nansum(w2*u[k]):.4g} | {np.nanmax(u[k]):.4g} | "
                             f"{lat[np.nanargmax(zm)]:.1f} |")
                axs[k].plot(lat, zm, "-", color=color[lab], lw=2, label=lab)
        for k, a in enumerate(axs):
            a.set_title(CLIP[k] + ", zonal mean", loc="left")
            a.set_xlabel("latitude (degrees)")
            a.set_ylabel("cm/yr")
            a.legend()
            style(a)
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, "clipping_zonal.png"), dpi=150)
        plt.close(fig)

        # maps of C1 and C2 for the runs that clip at all
        for lab in have:
            u = np.nanmean(runs[lab]["fields"]["user_t1"][:nyr*12], axis=0)
            if not np.nanmax(u[0]) > 0.0:
                continue
            fig, axs = plt.subplots(1, 2, figsize=(11, 3.4))
            vmax = max(np.nanmax(u[0]), np.nanmax(u[1]))
            for k, a in enumerate(axs):
                m = mapmesh(a, lat, lon, u[k], cmap=CMAP_SEQ, vmin=0.0, vmax=vmax)
                add_coast(a)
                fig.colorbar(m, ax=a, shrink=0.85, label="cm/yr")
                a.set_title(f"{lab}: {CLIP[k]}", loc="left")
                a.set_xlabel("longitude")
                a.set_ylabel("latitude")
            fig.tight_layout()
            fig.savefig(os.path.join(figdir, f"clipping_map_{lab}.png"), dpi=150)
            plt.close(fig)

    if want_report_fig:
        report_figure(figdir, runs, labels, ref, nyr, lat, lon)

    open(os.path.join(outdir, "compare_report.md"), "w").write("\n".join(lines) + "\n")
    print("compare_runs: report and figures in", outdir)


if __name__ == "__main__":
    main()
