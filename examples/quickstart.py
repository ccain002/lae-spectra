"""Quick tour of laespec.

    python examples/quickstart.py /path/to/LAE_Catalogs.tar     # or a fetched data directory

Writes quickstart.png in the current directory.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import laespec as ls  # noqa: E402

data = ls.LAEData(*sys.argv[1:])  # no argument: uses $LAESPEC_DATA
print(data.summary())

model = "late_end_early_start_13"
z = 8.4194

# The two samplings at the same redshift (this reproduces the original example.py):
bright = data.spectra(model, z, sampling="mass")      # 10000 most massive halos
typical = data.spectra(model, z, sampling="uniform")  # every k-th halo of the catalog

fig, ax = plt.subplots(1, 3, figsize=(15, 4))

ax[0].plot(bright.wavelength, bright.mean(), "r-", label="most massive (skip=1)")
ax[0].plot(typical.wavelength, typical.mean(), "k-", label="uniform (skip={})".format(typical.skip))
lo, hi = ls.transmission_percentiles(typical.flux, (16, 84))
ax[0].fill_between(typical.wavelength, lo, hi, color="k", alpha=0.15, label="16-84%")
ax[0].set(xlabel=r"rest wavelength [$\AA$]", ylabel="IGM transmission", title="z = {}".format(z))
ax[0].legend()

# Mean transmission in bins of halo mass, vs. velocity offset from line centre
logm = np.log10(typical.mass)
bins = [8.5, 9.5, 10.5, 11.5]
stacked, counts = ls.stack(typical.flux, logm, bins)
v = typical.velocity
for i, n in enumerate(counts):
    if n:
        ax[1].plot(v, stacked[i], label=r"$\log M$ = {}-{} (N={})".format(bins[i], bins[i + 1], n))
ax[1].set(xlabel="velocity offset [km/s]", ylabel="mean transmission", xlim=(-200, 3000))
ax[1].legend(fontsize=8)

# Transmitted fraction of a Gaussian Lya line (offset 200 km/s, sigma 100 km/s) vs redshift
for m in data.models:
    zs = data.redshifts(m, sampling="mass")
    if len(zs) < 3:
        continue
    zs = zs[::8]  # each file is 40 MB; a subset of redshifts keeps this quick
    tl = [np.median(ls.lya_transmission(s.wavelength, s.flux, dv=200, sigma_kms=100))
          for s in (data.spectra(m, zz, sampling="mass") for zz in zs)]
    ax[2].plot(zs, tl, label=m)
ax[2].set(xlabel="z", ylabel=r"median $T_{\rm Ly\alpha}$ (10$^4$ most massive halos)")
ax[2].legend(fontsize=8)

fig.tight_layout()
fig.savefig("quickstart.png", dpi=120)
print("wrote quickstart.png")
