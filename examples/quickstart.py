"""Quick tour of laespec.

    python examples/quickstart.py                              # bundled example (z = 6.97)
    python examples/quickstart.py /path/to/LAE_Catalogs.tar    # full data (or a directory of it)

With no argument and no $LAESPEC_DATA, this uses the spectrum file shipped
with the package. Writes quickstart.png in the current directory.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import laespec as ls  # noqa: E402

if len(sys.argv) > 1 or os.environ.get("LAESPEC_DATA"):
    data = ls.LAEData(*sys.argv[1:])
else:
    data = ls.example()
print(data.summary())

model = "late_end_early_start_13"
z = 7.0

typical = data.spectra(model, z, sampling="uniform")  # every k-th halo of the catalog
print("loaded {} sightlines at z = {} (skip = {})".format(len(typical), typical.z, typical.skip))

fig, ax = plt.subplots(1, 3, figsize=(15, 4))

# 1. Mean transmission and the 16-84% scatter between sightlines
ax[0].plot(typical.wavelength, typical.mean(), "k-", label="uniform (skip={})".format(typical.skip))
lo, hi = ls.transmission_percentiles(typical.flux, (16, 84))
ax[0].fill_between(typical.wavelength, lo, hi, color="k", alpha=0.15, label="16-84%")
if len(data.redshifts(model, sampling="mass")):  # full data only: the most massive halos
    bright = data.spectra(model, z, sampling="mass")
    ax[0].plot(bright.wavelength, bright.mean(), "r-", label="most massive (skip=1)")
ax[0].set(xlabel=r"rest wavelength [$\AA$]", ylabel="IGM transmission",
          title="{}, z = {}".format(model, typical.z))
ax[0].legend()

# 2. Mean transmission in bins of host-halo mass, vs velocity offset from line centre
logm = np.log10(typical.mass)
bins = [8.5, 9.5, 10.5, 11.5]  # the uniform sample has ~1 halo above 10^11.5
stacked, counts = ls.stack(typical.flux, logm, bins)
v = typical.velocity
for i, n in enumerate(counts):
    if n:
        ax[1].plot(v, stacked[i], label=r"$\log M$ = {}-{} (N={})".format(bins[i], bins[i + 1], n))
ax[1].set(xlabel="velocity offset [km/s]", ylabel="mean transmission", xlim=(-200, 3000))
ax[1].legend(fontsize=8)

# 3. Transmitted fraction of a Gaussian Lya line (offset 200 km/s, sigma 100 km/s)
T = ls.lya_transmission(typical.wavelength, typical.flux, dv=200, sigma_kms=100)
ax[2].plot(logm, T, ",", color="0.6", alpha=0.5)
med = [np.median(T[(logm >= a) & (logm < b)]) if np.any((logm >= a) & (logm < b)) else np.nan
       for a, b in zip(bins[:-1], bins[1:])]
ax[2].plot(0.5 * (np.array(bins[:-1]) + bins[1:]), med, "ko-", label="median")
ax[2].set(xlabel=r"$\log_{10} M_{\rm halo}$ [$M_\odot/h$]", ylabel=r"$T_{\rm Ly\alpha}$", ylim=(0, 1))
ax[2].legend()

fig.tight_layout()
fig.savefig("quickstart.png", dpi=120)
print("wrote quickstart.png")
