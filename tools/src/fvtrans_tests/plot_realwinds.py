#!/usr/bin/env python3
"""Figures and summary of fvt_realwinds (ROADMAP step 1.1, wiso:1d:0003).

Usage: python3 plot_realwinds.py <outdir>
Reads <outdir>/realwinds_steps.csv and realwinds_rows.csv; writes <outdir>/figures/realwinds_*.png and
<outdir>/realwinds_report.md. Needs numpy and matplotlib.
"""

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlesize": 10, "axes.titlecolor": INK, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "legend.frameon": False})


def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def load(path):
    rows = list(csv.DictReader(open(path)))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: plot_realwinds.py <outdir>")
    outdir = sys.argv[1]
    os.makedirs(os.path.join(outdir, "figures"), exist_ok=True)
    st = load(os.path.join(outdir, "realwinds_steps.csv"))
    rw = load(os.path.join(outdir, "realwinds_rows.csv"))
    nlat = int(rw["row"].max())
    lat = rw["lat"][:nlat]
    months = sorted(set(rw["month"].astype(int)))
    cxmax = np.array([rw["cx_max"][rw["month"] == m] for m in months]).max(0)
    cxmean = np.array([rw["cx_mean"][rw["month"] == m] for m in months]).mean(0)
    cymax = np.array([rw["cy_max"][rw["month"] == m] for m in months]).max(0)
    cymean = np.array([rw["cy_mean"][rw["month"] == m] for m in months]).mean(0)

    # Courant numbers by latitude
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
    for a, (mx, mn, name) in zip(ax, ((cxmax, cxmean, "zonal"), (cymax, cymean, "meridional"))):
        a.plot(lat, mx, "-o", color=SERIES[0], lw=2, ms=4, mec=SURF, label="maximum over the year")
        a.plot(lat, mn, "-o", color=SERIES[1], lw=2, ms=4, mec=SURF, label="mean over the year")
        a.axhline(1.0, color=INK2, lw=1, ls="--")
        a.set_title(f"{name} Courant number per 4 h step (umoisr x 800 hPa winds)", loc="left")
        a.set_xlabel("latitude (degrees)")
        style(a)
    ax[0].set_ylabel("Courant number")
    ax[0].legend(loc="upper center")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "figures", "realwinds_courant.png"), dpi=150)
    plt.close(fig)

    # substeps and clipping through the year
    day = (st["istep"] - 1)/6.0
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    for k, (key, name) in enumerate((("nsub_1", "out_max = 1"), ("nsub_075", "out_max = 0.75"))):
        v = st[key].astype(int)
        vals = np.arange(1, v.max() + 1)
        ax[0].bar(vals + (k - 0.5)*0.38, [(v == x).mean()*100 for x in vals], width=0.36, color=SERIES[k], label=name)
    ax[0].set_xlabel("FV substeps per 4 h step")
    ax[0].set_ylabel("share of steps (%)")
    ax[0].set_title("Substeps needed over one model year", loc="left")
    ax[0].legend()
    style(ax[0])
    ax[1].plot(day, st["clipped_frac"]*100, color=SERIES[0], lw=1, label="current spectral chain")
    ax[1].set_xlabel("day of the year")
    ax[1].set_ylabel("mass clipped into cormois (% of total)")
    ax[1].set_title("Clipping per step, current spectral chain", loc="left")
    style(ax[1])
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "figures", "realwinds_substeps_clipping.png"), dpi=150)
    plt.close(fig)

    def q(v, p):
        return np.percentile(v, p)

    lines = ["# FV transport with the real ECBilt 800 hPa winds (one model year)", "",
             f"{len(st['istep'])} steps of 4 h; fluxes from psi(:,3) and chi(:,3) times umoisr.", "",
             "| quantity | median | 99th pct | max |", "|---|---|---|---|"]
    for key, name in (("err_u", "u from psi, chi vs model u800+udivg (relative)"),
                      ("err_v", "v from psi, chi vs model v800+vdivg (relative)"),
                      ("err_div", "divergence from chi vs model divg (relative)"),
                      ("cx_max", "max zonal Courant number"), ("cy_max", "max meridional Courant number"),
                      ("out_max_nonshift", "max outflow fraction (cells with cx < 1)"),
                      ("clipped_frac", "spectral chain: mass clipped into cormois / total"),
                      ("clipped_frac_nodif", "same without ec_hdiff"),
                      ("rms_incr_fv", "rms one-step increment, FV pos.-def."),
                      ("rms_incr_spec", "rms one-step increment, spectral chain"),
                      ("rms_diff_nodif", "rms difference FV - spectral (no diffusion)")):
        v = st[key]
        lines.append(f"| {name} | {np.median(v):.3g} | {q(v, 99):.3g} | {v.max():.3g} |")
    lines += ["", f"FV minimum after one step over the year: {st['qmin_fv'].min():.3g}; "
              f"max relative mass change FV: {np.abs(st['mass_change_fv']).max():.2g}, "
              f"spectral: {np.abs(st['mass_change_spec']).max():.2g}.", ""]
    for key, name in (("nsub_1", "out_max = 1"), ("nsub_075", "out_max = 0.75")):
        v = st[key].astype(int)
        lines.append(f"Substeps with {name}: mean {v.mean():.2f}, distribution " +
                     ", ".join(f"{x}: {(v == x).mean()*100:.1f} %" for x in range(1, v.max() + 1)) + ".")
    lines += ["", "| row | latitude | cx max | cx mean | cy max | cy mean |", "|---|---|---|---|---|---|"]
    for i in range(nlat):
        lines.append(f"| {i + 1} | {lat[i]:.1f} | {cxmax[i]:.2f} | {cxmean[i]:.3f} | {cymax[i]:.2f} | {cymean[i]:.3f} |")
    open(os.path.join(outdir, "realwinds_report.md"), "w").write("\n".join(lines) + "\n")
    print("plot_realwinds: figures and realwinds_report.md in", outdir)


if __name__ == "__main__":
    main()
