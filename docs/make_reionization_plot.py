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
PLANCK_TAU, PLANCK_ERR = 0.058, 0.006  # Planck+18 band as plotted in Cohon et al. 2025
DEBELSUNCE_TAU, DEBELSUNCE_LO, DEBELSUNCE_HI = 0.0627, 0.0058, 0.0050  # de Belsunce et al. 2021


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
    """Same style as Fig. 1 of Cohon et al. 2025."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = np.genfromtxt(CSV, delimiter=",", names=True)
    z = d["z"]
    styles = [dict(color="k", ls="-"), dict(color="r", ls="--"), dict(color="b", ls=":")]

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for a in ax:
        a.minorticks_on()
        a.tick_params(which="both", direction="in", top=True, right=True, labelsize=14)
        a.tick_params(which="major", length=8, width=1)
        a.tick_params(which="minor", length=4, width=1)
        a.set_xlim(4.3, 18.3)
        a.set_xlabel("Redshift", fontsize=18, family="serif")

    for (model, name), st in zip(MODELS, styles):
        ax[0].plot(z, d[model + "_xHII_vol"], lw=2.5, label=name.upper(), **st)
        ax[1].plot(z, d[model + "_tau_lt_z"], lw=2.5, **st)

    ax[0].set_ylabel("Ionized Fraction", fontsize=18, family="serif")
    ax[0].set_ylim(-0.03, 1.05)
    ax[0].legend(loc="upper right", frameon=False, prop={"family": "serif", "size": 11}, handlelength=3)

    ax[1].axhspan(PLANCK_TAU - PLANCK_ERR, PLANCK_TAU + PLANCK_ERR, color="g", alpha=0.3, lw=0, label="Planck+18")
    ax[1].axhline(PLANCK_TAU, color="darkgreen", lw=2.5)
    ax[1].axhspan(DEBELSUNCE_TAU - DEBELSUNCE_LO, DEBELSUNCE_TAU + DEBELSUNCE_HI, color="m", alpha=0.3, lw=0,
                  label="de Belsunce+21")
    ax[1].axhline(DEBELSUNCE_TAU, color="m", ls="-.", lw=2.5)
    ax[1].set_ylabel(r"$\tau_{\rm es}$", fontsize=18)
    ax[1].set_ylim(0.041, 0.09)
    ax[1].legend(loc="upper left", fontsize=14)

    fig.tight_layout()
    fig.savefig(PNG, dpi=150)
    print("wrote", PNG)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--otf-dir", help="rebuild the CSV from FlexRT otf files in this directory first")
    a = p.parse_args()
    if a.otf_dir:
        build_csv(a.otf_dir)
    plot()
