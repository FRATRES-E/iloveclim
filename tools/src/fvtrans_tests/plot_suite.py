#!/usr/bin/env python3
"""Figures and summary tables of the FV transport evaluation suite (ROADMAP step 1.1, wiso:1c:0001).

Usage: python3 plot_suite.py <outdir>
Reads <outdir>/summary.csv and the netCDF files written by fvt_suite; writes <outdir>/figures/*.png and
<outdir>/report.md. Needs numpy, matplotlib, netCDF4.
"""

import csv
import os
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from netCDF4 import Dataset  # noqa: E402

# reference palette (dataviz skill): categorical slots in fixed order, one-hue sequential ramp, text inks
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SCHEMES = ["ppm_none", "ppm_mono", "ppm_pdef", "vanleer", "pdef_mono"]
LABEL = {"ppm_none": "PPM unlimited", "ppm_mono": "PPM monotone", "ppm_pdef": "PPM pos.-def.",
         "vanleer": "van Leer", "pdef_mono": "pos.-def. carrier + monotone ratios"}
COLOR = dict(zip(SCHEMES, SERIES))
CONV_RUNS = [(32, 72), (64, 144), (128, 288), (256, 576)]   # constant Courant numbers
DLON = {32: 5.625, 64: 2.8125, 128: 1.40625, 256: 0.703125}

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlesize": 10, "axes.titlecolor": INK, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "legend.frameon": False})
SEQ_CMAP = LinearSegmentedColormap.from_list("seq_blue", SEQ)


def load(outdir):
    rows = list(csv.DictReader(open(os.path.join(outdir, "summary.csv"))))
    table = {}
    for r in rows:
        key = (int(r["nlat"]), int(r["nstep"]), r["timing"], r["case"], r["mode"], r["scheme"], r["tracer"], r["metric"])
        table[key] = float(r["value"])
    return table


def get(t, nlat, nstep, case, mode, scheme, tracer, metric, timing="mid"):
    return t.get((nlat, nstep, timing, case, mode, scheme, tracer, metric), np.nan)


def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_convergence(t, outdir):
    panels = [("deform", "density", "hills", "Gaussian hills, non-divergent flow"),
              ("sbr_a90", "density", "wg_bell", "Cosine bell, rotation over the pole"),
              ("divergent", "ratio", "bells", "Cosine bells, divergent flow (C')")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, (case, mode, tracer, title) in zip(axes, panels):
        avail = [(n, s) for n, s in CONV_RUNS if not np.isnan(get(t, n, s, case, mode, "ppm_none", tracer, "l2"))]
        x = np.array([DLON[n] for n, _ in avail])
        for sch in SCHEMES[:4]:
            y = np.array([get(t, n, s, case, mode, sch, tracer, "l2") for n, s in avail])
            ax.loglog(x, y, "-o", color=COLOR[sch], lw=2, ms=5, mec=SURF, mew=1.5, label=LABEL[sch])
        if len(x) > 1:
            ref = get(t, avail[0][0], avail[0][1], case, mode, "ppm_none", tracer, "l2")
            ax.loglog(x, ref*(x/x[0])**2, "--", color=INK2, lw=1, label="2nd-order slope")
        ax.set_xlabel("grid spacing (degrees)")
        ax.set_title(title, loc="left")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{v:.2g}" for v in x])
        ax.minorticks_off()
        ax.invert_xaxis()
        style(ax)
    axes[0].set_ylabel("normalised l2 error at t = T")
    axes[0].legend(loc="lower left", fontsize=8)
    fig.suptitle("Convergence at constant Courant number (T21/72, T42/144, T85/288, T170/576 steps per period)",
                 x=0.01, ha="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "figures", "convergence.png"), dpi=150)
    plt.close(fig)


def read_nc(path, var):
    with Dataset(path) as nc:
        lonb = nc["lon_bnds"][:]
        latb = nc["lat_bnds"][:]
        f = nc[var][:]
    lon_edges = np.append(lonb[:, 0], lonb[-1, 1])
    lat_edges = np.append(latb[:, 0], latb[-1, 1])
    return lon_edges, lat_edges, np.asarray(f)


def fig_maps(outdir, nlat, nstep, case, mode, var, schemes, title, fname):
    nrow = len(schemes)
    fig, axes = plt.subplots(nrow, 3, figsize=(11, 2.15*nrow + 0.7), sharex=True, sharey=True, squeeze=False,
                             constrained_layout=True)
    times = ["t = 0", "t = T/2", "t = T"]
    under_cmap = LinearSegmentedColormap.from_list("under", [SERIES[1], SERIES[1]])
    pm = None
    for row, sch in enumerate(schemes):
        path = os.path.join(outdir, f"{case}_{mode}_{sch}_n{nlat}_mid.nc")
        if not os.path.exists(path):
            continue
        lone, late, f = read_nc(path, var)
        for col in range(3):
            ax = axes[row, col]
            pm = ax.pcolormesh(lone, late, f[col], cmap=SEQ_CMAP, vmin=0.0, vmax=1.0, shading="flat")
            # cells below the initial minimum (0.1) shaded orange
            und = np.ma.masked_where(f[col] >= 0.1 - 1.0e-6, np.ones_like(f[col]))
            ax.pcolormesh(lone, late, und, cmap=under_cmap, vmin=0, vmax=1, shading="flat")
            ax.set_title(f"{LABEL[sch]}, {times[col]}", loc="left", fontsize=8, color=INK2)
            ax.set_xlim(0, 360)
            ax.set_ylim(-90, 90)
            ax.set_aspect("equal")
            ax.set_yticks([-60, 0, 60])
            ax.set_xticks([0, 90, 180, 270, 360])
    if pm is not None:
        cb = fig.colorbar(pm, ax=axes, shrink=0.8, pad=0.01)
        cb.set_label("mixing ratio")
    fig.suptitle(title + "; orange cells: below the initial minimum 0.1", x=0.0, ha="left", color=INK, fontsize=10)
    fig.savefig(os.path.join(outdir, "figures", fname), dpi=150)
    plt.close(fig)


def fig_mixing(outdir, nlat, t):
    chi = np.linspace(0.1, 1.0, 200)
    configs = [("deform", "density", "ppm_mono", "non-divergent, PPM monotone"),
               ("divergent", "density", "ppm_mono", "divergent, option C (independent densities)"),
               ("divergent", "ratio", "ppm_mono", "divergent, option C' (ratios to the carrier)")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharex=True, sharey=True)
    for ax, (case, mode, sch, title) in zip(axes, configs):
        path = os.path.join(outdir, f"{case}_{mode}_{sch}_n{nlat}_mid.nc")
        if not os.path.exists(path):
            continue
        with Dataset(path) as nc:
            c = np.asarray(nc["bells"][1]).ravel()
            x = np.asarray(nc["correlated"][1]).ravel()
        ax.plot(chi, -0.8*chi**2 + 0.9, color=INK2, lw=1)
        ax.plot([0.1, 1.0], [-0.8*0.01 + 0.9, 0.1], color=INK2, lw=1, ls="--")
        ax.plot([0.1, 1.0, 1.0, 0.1, 0.1], [0.1, 0.1, 0.892, 0.892, 0.1], color=GRID, lw=1)
        ax.scatter(c, x, s=8, color=SERIES[0], edgecolors=SURF, linewidths=0.3)
        lr, lu, lo = (get(t, nlat, CONV_STEPS[nlat], case, mode, sch, "pair", m) for m in ("lr", "lu", "lo"))
        ax.set_title(title, loc="left", fontsize=9)
        ax.text(0.03, 0.04, f"l_r {lr:.2e}\nl_u {lu:.2e}\nl_o {lo:.2e}", transform=ax.transAxes, fontsize=8,
                color=INK, va="bottom")
        ax.set_xlabel("cosine bells")
        style(ax)
    axes[0].set_ylabel("correlated cosine bells")
    fig.suptitle(f"Mixing at t = T/2 (T{ {32: 21, 64: 42, 128: 85, 256: 170}[nlat] }): curve = initial relation, "
                 "dashed = chord (real-mixing region in between)", x=0.01, ha="left", color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "figures", f"mixing_n{nlat}.png"), dpi=150)
    plt.close(fig)


CONV_STEPS = dict(CONV_RUNS)


def report(t, outdir):
    lines = ["# FV transport suite: summary", "",
             "All errors normalised as in Williamson et al. (1992) / Lauritzen et al. (2012); exact solution at t = T is "
             "the initial state. `under`/`over`: largest excursion below/above the initial range at any step, "
             "relative to that range.", ""]
    for nlat, nstep in CONV_RUNS:
        if np.isnan(get(t, nlat, nstep, "deform", "density", "ppm_none", "hills", "l2")):
            continue
        lines += [f"## T{ {32: 21, 64: 42, 128: 85, 256: 170}[nlat] } ({nlat} x {2*nlat}), {nstep} steps per period", ""]
        lines += ["| scheme | SBR a=0 l2 | SBR a=90 l2 | hills l2 | bells l2 | cyl. l2 | cyl. under | cyl. over | "
                  "lf(0.5) | l_r | l_u | l_o |", "|" + "---|"*12]
        for sch in SCHEMES[:4]:
            g = lambda c, tr, m, mode="density": get(t, nlat, nstep, c, mode, sch, tr, m)  # noqa: E731
            lines.append(f"| {sch} | {g('sbr_a00','wg_bell','l2'):.3g} | {g('sbr_a90','wg_bell','l2'):.3g} | "
                         f"{g('deform','hills','l2'):.3g} | {g('deform','bells','l2'):.3g} | "
                         f"{g('deform','cylinders','l2'):.3g} | {g('deform','cylinders','trans_under'):.1e} | "
                         f"{g('deform','cylinders','trans_over'):.1e} | {g('deform','bells','lf_0.50'):.1f} | "
                         f"{g('deform','pair','lr'):.2e} | {g('deform','pair','lu'):.2e} | {g('deform','pair','lo'):.2e} |")
        lines += ["", "Divergent flow, C (density) vs C' (ratio):", "",
                  "| mode | scheme | bells l2 | cyl. under | cyl. over | corr. under | corr. over | l_r | l_u | l_o |",
                  "|" + "---|"*10]
        for mode in ("density", "ratio"):
            for sch in SCHEMES:
                g = lambda tr, m: get(t, nlat, nstep, "divergent", mode, sch, tr, m)  # noqa: E731
                if np.isnan(g("bells", "l2")):
                    continue
                lines.append(f"| {mode} | {sch} | {g('bells','l2'):.3g} | {g('cylinders','trans_under'):.1e} | "
                             f"{g('cylinders','trans_over'):.1e} | {g('correlated','trans_under'):.1e} | "
                             f"{g('correlated','trans_over'):.1e} | {g('pair','lr'):.2e} | {g('pair','lu'):.2e} | "
                             f"{g('pair','lo'):.2e} |")
        lines.append("")
    # sensitivity at T21: time step and wind timing
    lines += ["## T21 sensitivity: step and wind timing (PPM pos.-def. / monotone)", "",
              "| run | SBR a=90 l2 (pdef) | hills l2 (pdef) | bells l2 (pdef) | cyl. under (mono) | "
              "div. C' cyl. over (mono) |", "|" + "---|"*6]
    for nstep, timing in ((72, "mid"), (72, "start"), (144, "mid")):
        g = lambda c, mo, s, tr, m: get(t, 32, nstep, c, mo, s, tr, m, timing)  # noqa: E731
        if np.isnan(g("deform", "density", "ppm_pdef", "hills", "l2")):
            continue
        lines.append(f"| {nstep} steps, winds at {timing} | {g('sbr_a90','density','ppm_pdef','wg_bell','l2'):.3g} | "
                     f"{g('deform','density','ppm_pdef','hills','l2'):.3g} | "
                     f"{g('deform','density','ppm_pdef','bells','l2'):.3g} | "
                     f"{g('deform','density','ppm_mono','cylinders','trans_under'):.1e} | "
                     f"{g('divergent','ratio','ppm_mono','cylinders','trans_over'):.1e} |")
    open(os.path.join(outdir, "report.md"), "w").write("\n".join(lines) + "\n")


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: plot_suite.py <outdir>")
    outdir = sys.argv[1]
    os.makedirs(os.path.join(outdir, "figures"), exist_ok=True)
    t = load(outdir)
    fig_convergence(t, outdir)
    for nlat in (32, 64):
        fig_maps(outdir, nlat, CONV_STEPS[nlat], "deform", "density", "cylinders", ["ppm_mono", "ppm_pdef"],
                 f"Slotted cylinders, non-divergent flow, T{ {32: 21, 64: 42}[nlat] }", f"cylinders_deform_n{nlat}.png")
        fig_maps(outdir, nlat, CONV_STEPS[nlat], "divergent", "density", "cylinders", ["ppm_mono"],
                 f"Slotted cylinders, divergent flow, option C, T{ {32: 21, 64: 42}[nlat] }", f"cylinders_div_C_n{nlat}.png")
        fig_maps(outdir, nlat, CONV_STEPS[nlat], "divergent", "ratio", "cylinders", ["ppm_mono"],
                 f"Slotted cylinders, divergent flow, option C', T{ {32: 21, 64: 42}[nlat] }",
                 f"cylinders_div_Cp_n{nlat}.png")
        fig_mixing(outdir, nlat, t)
    report(t, outdir)
    print("plot_suite: figures and report.md in", outdir)


if __name__ == "__main__":
    main()
