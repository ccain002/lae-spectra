"""Reionization histories and CMB optical depth of the three models -> docs/reionization_histories.png

    python docs/make_reionization_plot.py                       # plot from docs/reionization_histories.csv
    python docs/make_reionization_plot.py --otf-dir DIR         # first rebuild the CSV from FlexRT output

DIR holds full_mono_<model>/output_files/otf_<model>.txt (FlexRT's on-the-fly global averages),
e.g. /expanse/lustre/projects/ucr138/ccain002/AGN_reionization/JWST_paper_models.

tau(<z) = c sigma_T int_0^z n_e(z') / ((1+z') H(z')) dz', with n_e from the mass-weighted ionized
H fraction, He II ionized alongside H, He III fully ionized below z = 3.5, and gas fully ionized
below the end of the simulations (z = 4.8). Cosmology as in FlexRT: Om = 0.305147, Ob = 0.048,
h = 0.68, Y = 0.246.
"""
import argparse
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "reionization_histories.csv")
PNG = os.path.join(HERE, "reionization_histories.png")

# model label in the data files -> name in Cohon et al. 2025
MODELS = [("late_end_late_start_4", "late start"),
          ("late_end_early_start_13", "early start"),
          ("late_end_very_early_start_3", "very early start")]

OM, OB, HH, Y = 0.305147, 0.048, 0.68, 0.246
SIGMA_T, C, G, M_H, MPC = 6.6524587e-25, 2.99792458e10, 6.6743e-8, 1.6735575e-24, 3.0856776e24
Z_HE3 = 3.5
PLANCK_TAU, PLANCK_ERR = 0.058, 0.006  # value quoted in Cohon et al. 2025


def read_otf(path):
    txt = open(path).read()
    z = np.array([float(v) for v in re.findall(r"Redshift\s*:\s*(\S+)", txt)])
    f = np.array([[float(a), float(b)] for a, b in re.findall(r"\nfH2\s*:\s*(\S+)\s+(\S+)", txt)])
    assert len(z) == len(f) and np.all(np.diff(z) <= 0)
    return z[::-1], f[::-1, 0], f[::-1, 1]  # ascending z; volume- and mass-weighted x_HII


def tau_cumulative(z_sim, x_mass, zmax=30.0):
    """tau(<z) on z_sim, for the history x_mass(z_sim) (ascending z)."""
    H0 = HH * 100e5 / MPC
    nH0 = (1 - Y) * OB * 3 * H0 ** 2 / (8 * np.pi * G) / M_H
    y = Y / (4 * (1 - Y))
    z = np.unique(np.concatenate([np.linspace(0, zmax, 30001), z_sim]))
    x = np.interp(z, z_sim, x_mass, left=1.0, right=0.0)
    x[z < z_sim[0]] = 1.0
    ne = nH0 * (1 + z) ** 3 * (x * (1 + y) + y * (z < Z_HE3))
    integrand = C * SIGMA_T * ne / ((1 + z) * H0 * np.sqrt(OM * (1 + z) ** 3 + 1 - OM))
    cum = np.concatenate([[0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(z))])
    return np.interp(z_sim, z, cum), cum[-1]


def build_csv(otf_dir):
    cols, header = None, ["z"]
    for model, _ in MODELS:
        z, xv, xm = read_otf(os.path.join(otf_dir, "full_mono_" + model, "output_files", "otf_" + model + ".txt"))
        tau, _ = tau_cumulative(z, xm)
        if cols is None:
            cols, zref = [z], z
        assert np.allclose(z, zref), "models are on different redshift grids"
        cols += [xv, xm, tau]
        header += [model + "_xHII_vol", model + "_xHII_mass", model + "_tau_lt_z"]
    keep = np.unique(np.concatenate([np.arange(0, len(zref), 4), [len(zref) - 1]]))
    np.savetxt(CSV, np.column_stack(cols)[keep], delimiter=",", header=",".join(header), comments="",
               fmt="%.6g")
    print("wrote", CSV)


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = np.genfromtxt(CSV, delimiter=",", names=True)
    z = d["z"]
    # default categorical slots 1-3 (blue, orange, aqua), plus line style as a second cue
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]
    styles = ["-", "--", "-."]
    ink, ink2, grid, surface = "#1f1f1e", "#5c5b55", "#e4e3df", "#fcfcfb"

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "text.color": ink})
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9), facecolor=surface)
    for a in ax:
        a.set_facecolor(surface)
        a.grid(color=grid, lw=0.8)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.set_xlim(4.8, 16)
        a.set_xlabel("redshift z")
        a.axvspan(5.0, 13.04, color="#efeeea", lw=0, zorder=0)

    for (model, name), c, ls in zip(MODELS, colors, styles):
        x = d[model + "_xHII_vol"]
        tau = d[model + "_tau_lt_z"]
        tau_tot = tau_cumulative(z, d[model + "_xHII_mass"])[1]
        ax[0].plot(z, x, color=c, ls=ls, lw=2, label=name)
        ax[1].plot(z, tau, color=c, ls=ls, lw=2, label="{}: τ = {:.3f}".format(name, tau_tot))
        zl, off, ha = (7.4, (-6, -4), "right") if name == "late start" else (10.5, (6, 9), "left")
        i = np.argmin(abs(z - zl))
        ax[0].annotate(name, (z[i], x[i]), xytext=off, textcoords="offset points", fontsize=9, color=ink,
                       ha=ha, va="top" if ha == "right" else "baseline")

    ax[0].set(ylabel="volume-weighted ionized fraction", ylim=(0, 1.02),
              title="Reionization histories")
    ax[0].text(9.0, 0.96, "range covered by the spectra (z = 5–13)", fontsize=8, color=ink2, ha="center")

    ax[1].axhspan(PLANCK_TAU - PLANCK_ERR, PLANCK_TAU + PLANCK_ERR, color="#d9d8d2", lw=0, zorder=0.5)
    ax[1].axhline(PLANCK_TAU, color=ink2, lw=1, zorder=0.6)
    ax[1].text(15.9, PLANCK_TAU + PLANCK_ERR + 0.0015, "Planck τ = 0.058 ± 0.006", fontsize=8, color=ink2,
               ha="right", va="bottom")
    ax[1].set(ylabel="CMB optical depth τ(<z)", ylim=(0.03, 0.09), title="Thomson optical depth")
    ax[1].legend(loc="lower right", frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(PNG, dpi=150, facecolor=surface)
    print("wrote", PNG)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--otf-dir", help="rebuild the CSV from FlexRT otf files in this directory first")
    a = p.parse_args()
    if a.otf_dir:
        build_csv(a.otf_dir)
    plot()
