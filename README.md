# laespec

Helpers for reading and manipulating the simulated **Lyman-alpha emitter (LAE)
damping-wing spectra**: IGM transmission along sightlines toward halos at
5 < z < 13, for several reionization histories, plus the matching halo
catalogs.

- Reads the data straight out of the original `LAE_Catalogs.tar`, the
  downloadable per-model archives, or extracted files. No unpacking needed:
  files are memory-mapped in place.
- Matches each sightline to its host halo (mass, position, velocity).
- Provides vectorised helpers for velocity offsets, resampling, rebinning,
  instrumental smoothing, stacking in halo-mass bins, bootstrap errors, and the
  transmitted fraction of an intrinsic Lyα line (T_Lyα).
- Ships with one real spectrum file (z ≈ 7), so examples and tests run out of the box.

Only `numpy` is required (Python ≥ 3.7). `matplotlib` is needed only for the example.
The install includes a ~41 MB example spectrum file, so you can try everything right away.

## Install

### Quick install (not on Ubuntu)

If your Python lets you `pip install` directly (macOS, conda, most HPC
modules), this is all you need:

```bash
pip install git+https://github.com/ccain002/lae-spectra.git
# or, for development
git clone https://github.com/ccain002/lae-spectra.git && cd lae-spectra && pip install -e .
```

**This won't work with the system Python on Ubuntu 23.04 or later** (or
Debian 12+). There, pip refuses with `error: externally-managed-environment`.
Use the venv steps below. They also work everywhere else, and keep `laespec`
separate from your other Python packages.

### Ubuntu (20.04, 22.04, 24.04) and anywhere else: virtual environment (venv)

```bash
# 1. One-time system packages (venv is not installed by default on Ubuntu)
sudo apt update
sudo apt install -y python3-venv python3-pip git

# 2. Get the code
git clone https://github.com/ccain002/lae-spectra.git
cd lae-spectra

# 3. Create and activate a venv inside the repo
python3 -m venv .venv
source .venv/bin/activate

# 4. Install laespec plus matplotlib (for the example) and pytest
python -m pip install --upgrade pip
python -m pip install -e ".[plot,test]"

# 5. Check it works
pytest tests
python examples/quickstart.py      # writes quickstart.png
```

Run `source .venv/bin/activate` again in every new terminal before using
`laespec`, and `deactivate` to leave the venv. `.venv/` is already in
`.gitignore`. If step 3 says `ensurepip is not available`, install the
versioned package instead, e.g. `sudo apt install python3.12-venv`, where the
version matches `python3 --version`. Then delete the half-made `.venv` and
redo step 3.

### Other systems (macOS, HPC clusters, Windows/WSL)

The same steps work without the `apt` line. On a cluster, load a Python
module first (e.g. `module load python`). On native Windows, activate with
`.venv\Scripts\activate`.

### Without a clone

To install the package into an active venv without keeping a checkout:

```bash
python -m pip install --upgrade pip   # older pip rejects the "name @ url" form
python -m pip install "laespec[plot] @ git+https://github.com/ccain002/lae-spectra.git"
```

The bundled example data is included either way, but the tests and
`examples/` are only in a clone.

### Using it in Jupyter

```bash
python -m pip install ipykernel
python -m ipykernel install --user --name laespec --display-name "Python (laespec)"
```

Then pick the "Python (laespec)" kernel in Jupyter.

## Try it without downloading anything

The package ships with one real spectrum file (about 41 MB): model
`late_end_early_start_13`, z = 6.9746, `skip=539`, so 10,000 sightlines that
uniformly sample the halo catalog. It includes each sightline's host-halo
properties.

```python
import laespec as ls
spec = ls.example().spectra("late_end_early_start_13", z=7.0)
spec.flux.shape, spec.mass[:3]
```

```bash
python examples/quickstart.py     # plots from the bundled file -> quickstart.png
pytest tests                      # includes checks against the bundled file
```

## Get the full data

The full dataset is 45.5 GB: 3 reionization models, 85 redshifts, and the
halo catalogs. It is shared through **Globus**, as a single file,
`LAE_Catalogs.tar`:

**Globus collection:**
<https://app.globus.org/file-manager?origin_id=e06caa06-533b-47fb-b874-7c120d8cf896&origin_path=%2F>

1. Open the link and sign in to Globus. Your university login usually works.
   If it says you don't have access, ask the maintainer.
2. In the other pane, pick where the file should go: your cluster's Globus
   collection, or [Globus Connect Personal](https://www.globus.org/globus-connect-personal)
   on a laptop or workstation.
3. Select `LAE_Catalogs.tar` and click **Start**. Globus verifies the
   transfer's checksum and emails you when it's done.

There's no need to extract anything:

```python
data = ls.LAEData("/path/to/LAE_Catalogs.tar")   # or the directory you transferred into
```

```bash
export LAESPEC_DATA=/path/to/LAE_Catalogs.tar    # optional: lets LAEData() find it by default
laespec info /path/to/LAE_Catalogs.tar           # summary of models / redshifts / files
```

If you have an Expanse account, skip the transfer and read it in place:
`LAEData("/expanse/lustre/projects/ucr138/ccain002/LAE_Catalogs/LAE_Catalogs.tar")`.

## Usage

```python
import numpy as np
import laespec as ls

data = ls.LAEData("LAE_Catalogs.tar")    # or a directory; LAEData() uses $LAESPEC_DATA; ls.example() for the bundled file
data.models                               # ['late_end_early_start_13', 'late_end_late_start_4', 'late_end_very_early_start_3']
data.redshifts("late_end_early_start_13")

# One file = 10000 sightlines at one redshift
spec = data.spectra("late_end_early_start_13", z=7.0, sampling="uniform")  # nearest snapshot within dz = 0.05 (here 6.9746)
spec.wavelength      # (1000,) rest-frame Angstrom, ascending
spec.flux            # (10000, 1000) transmission, lazy read-only memmap
spec.velocity        # km/s from Lya line centre
spec.mass            # (10000,) host halo mass [Msun/h]
spec.halos["x"]      # also y, z [cMpc/h], vx, vy, vz [km/s]

bright = spec.select(spec.mass > 1e11)            # subset of sightlines
mean = bright.mean()
lo, med, hi = ls.transmission_percentiles(spec.flux, (16, 50, 84))

# Transmitted fraction of a Gaussian Lya line, 200 km/s redward, sigma 100 km/s
T_lya = ls.lya_transmission(spec.wavelength, spec.flux, dv=200, sigma_kms=100)

# Mean spectrum in halo-mass bins
stacked, counts = ls.stack(spec.flux, np.log10(spec.mass), bins=[9, 10, 11, 12])

# Instrument-like processing
sm = ls.smooth(spec.wavelength, spec.flux, resolution=1000)       # or fwhm_kms=...
wl2, f2 = ls.rebin(spec.wavelength, sm, 4)
T_500 = ls.transmission_at(spec.wavelength, spec.flux, dv=500)     # per sightline
obs_wl = spec.observed_wavelength                                   # (1+z) * rest
```

`sampling="mass"` gives the file with `skip=1` (the 10⁴ most massive halos).
`sampling="uniform"` gives `skip=k>1`, every k-th halo of the mass-sorted
catalog. That is a uniform sample of the halo population, which is dominated
by faint, low-mass halos. Pass `skip=` to choose a file exactly.

See [`examples/quickstart.py`](examples/quickstart.py) for a complete script
that makes plots. It runs on the bundled file by default.

### Function reference

| | |
|---|---|
| `example()` | `LAEData` for the bundled z = 6.97 file |
| `LAEData(*paths)` | index the data; `.models`, `.redshifts()`, `.table()`, `.summary()`, `.find()` |
| `LAEData.spectra(model, z, sampling, skip, ascending, mmap)` | load one file → `Spectra` |
| `LAEData.halo_catalog(z)` | full halo catalog (structured array) |
| `LAEData.matched_halos(file)` | halo rows for each sightline |
| `velocity_offset`, `wavelength_from_velocity`, `to_observed`, `to_rest` | coordinates |
| `mean_transmission`, `transmission_percentiles`, `optical_depth`, `bootstrap`, `stack` | statistics |
| `interpolate`, `transmission_at`, `rebin`, `smooth`, `ensure_ascending` | resampling |
| `line_profile`, `lya_transmission` | Lyα line transmission |

All helper functions in `laespec.spectra` work on plain arrays with the
spectral axis last, so you can use them on any spectra.

## Data format

All files are raw little-endian binary with no header. You can read them
without this package: see `example.py` in the original tarball and the
docstring of [`laespec/catalog.py`](laespec/catalog.py).

| file | dtype | shape | contents |
|---|---|---|---|
| `wav_out_LAE` | float32 | (1000,) | rest wavelength [Å], **descending** 1315.67 → 1210.78, Δλ ≈ 0.105 Å |
| `spec_out_mono_<model>.z=<z>_skip<k>_LAE` | float32 | (N, 1000) | transmission e^-τ on that grid; N = 10000 |
| `halo_xvm_catalog.z=<z>` | float64 | (N_halo, 7) | x, y, z [cMpc/h], vx, vy, vz [km/s], M [Msun/h], sorted by descending M |
| `halos_for_<spectrum file>` | float64 | (N, 7) | the halo rows matching each sightline (in the downloadable archives) |

Sightline `i` of a `skip<k>` file belongs to halo row `i*k` of the halo
catalog at the same redshift. With `skip` = N_halo // 10000, the sample
spans the whole catalog. The box is 200 cMpc/h on a side (positions span 0–200).

**Models** (reionization histories; `laespec` drops the `mono_` prefix from
file names). All three come from the simulations of
[Cohon et al. 2025](https://arxiv.org/abs/2508.05739):

- **Ionization fields:** the FlexRT radiative-transfer code
  ([Cain & D'Aloisio 2024](https://arxiv.org/abs/2409.04521)), run in a
  200 h⁻¹ cMpc box.
- **Gas density and velocity:** a 2048³ RadHydro simulation of the same volume
  (Δx ≈ 97 h⁻¹ kpc).
- **Shared calibration:** all three histories end reionization late, at
  z ≈ 5–5.5, so they match the Lyα forest mean transmission at z ≤ 6. The
  ionizing emissivity is shared among halos in proportion to their UV
  luminosity.
- **What differs:** how early reionization starts.

| model | name in Cohon et al. 2025 | ionized fraction at z = 13 | τ_CMB |
|---|---|---|---|
| `late_end_late_start_4` | late start | < 1% | 0.049, slightly more than 1σ below Planck's 0.058 ± 0.006 |
| `late_end_early_start_13` | early start | ≈ 5% | 0.063, within 1σ |
| `late_end_very_early_start_3` | very early start | ≈ 25% | 0.084, in tension with Planck |

![Reionization histories and CMB optical depth of the three models](docs/reionization_histories.png)

The histories come from FlexRT's global averages. τ is computed from them with
He II reionization at z = 3.5. The histories and τ(<z) are in
[`docs/reionization_histories.csv`](docs/reionization_histories.csv), and
[`docs/make_reionization_plot.py`](docs/make_reionization_plot.py) remakes the
figure.

Each model has 85 snapshots, z = 5.00–13.03, with a `skip=1` file and a
uniform file at each redshift. These spectra use the same simulations as that
paper, but select their sightlines differently: one per halo, ranked by mass
or sampled uniformly.

## Tests

With the venv active (see [Install](#install)):

```bash
python -m pip install -e ".[test]"   # already done if you followed the install steps
pytest tests
```

Most tests build a tiny synthetic dataset. `test_example_data` checks the
bundled real file against reference values.

## For the data maintainer

[`docs/SHARING.md`](docs/SHARING.md) explains how to share the full data over
Globus. It also covers the optional `laespec pack` / `fetch` route through
Zenodo.
