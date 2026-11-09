"""Locate and load LAE damping-wing spectra and halo catalogs.

The data can live in any mix of

* the original ``LAE_Catalogs.tar`` archive,
* the per-model archives produced by ``laespec pack``, or
* a directory of extracted files,

and :class:`LAEData` reads all of them the same way. Files inside an
*uncompressed* tar are memory-mapped in place, so nothing has to be
extracted first.

File formats (all raw binary, native little-endian, no header)
---------------------------------------------------------------
``wav_out_LAE``
    float32, 1000 values: rest-frame wavelength in Angstrom, stored in
    *descending* order (1315.67 -> 1210.78 A).
``spec_out_<model>.z=<zz.zzzz>_skip<k>_LAE``
    float32, shape (n_sightlines, 1000): IGM transmission exp(-tau) on the
    ``wav_out_LAE`` grid. Row ``i`` is halo number ``i * k`` in the halo
    catalog at the same redshift (``k = 1``: the 10000 most massive halos;
    ``k > 1``: uniform sampling of the whole catalog).
``halo_xvm_catalog.z=<zz.zzzz>``
    float64, shape (n_halos, 7): x, y, z [cMpc/h], vx, vy, vz [km/s],
    M [Msun/h], sorted by descending mass.
``halos_for_<spectrum file name>``  (written by ``laespec pack``)
    Same layout as a halo catalog, holding only the rows matched to the
    sightlines of one spectrum file.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tarfile
import warnings
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

N_PIX = 1000
N_HALO_DEFAULT = 10000
SPEC_DTYPE = np.dtype("<f4")
HALO_DTYPE = np.dtype("<f8")
HALO_FIELDS = ("x", "y", "z", "vx", "vy", "vz", "mass")
HALO_RECORD = np.dtype([(f, HALO_DTYPE) for f in HALO_FIELDS])

WAVELENGTH_FILE = "wav_out_LAE"
SPEC_RE = re.compile(r"^spec_out_(?P<model>.+)\.z=(?P<zstr>\d+\.\d+)_skip(?P<skip>\d+)_LAE$")
HALO_RE = re.compile(r"^halo_xvm_catalog\.z=(?P<zstr>\d+\.\d+)$")
MATCHED_PREFIX = "halos_for_"

_CACHE_DIR = os.environ.get(
    "LAESPEC_CACHE", os.path.join(os.path.expanduser("~"), ".cache", "laespec")
)


# --------------------------------------------------------------------------
# file index
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Entry:
    """Where the bytes of one data file live (a plain file or a tar member)."""

    path: str      # file on disk (the tar itself for archive members)
    offset: int    # byte offset of the data within ``path``
    size: int      # size in bytes
    member: Optional[str] = None  # member name if inside a tar

    def memmap(self, dtype, shape=None, mmap=True):
        dtype = np.dtype(dtype)
        count = self.size // dtype.itemsize
        if shape is None:
            shape = (count,)
        if count == 0:  # np.memmap cannot map empty files
            return np.zeros(shape, dtype=dtype)
        if mmap:
            return np.memmap(self.path, dtype=dtype, mode="r", offset=self.offset, shape=shape)
        with open(self.path, "rb") as fh:
            fh.seek(self.offset)
            return np.fromfile(fh, dtype=dtype, count=count).reshape(shape)

    def read_bytes(self):
        with open(self.path, "rb") as fh:
            fh.seek(self.offset)
            return fh.read(self.size)


def _index_tar(path: str, use_cache: bool = True) -> List[Entry]:
    """List the regular-file members of an uncompressed tar, with data offsets.

    Scanning the 48 GB original archive takes ~10 s, so the result is cached
    in ``~/.cache/laespec`` (override with ``$LAESPEC_CACHE``), keyed on the
    archive's path, size and mtime.
    """
    st = os.stat(path)
    key = hashlib.sha1(
        "{}|{}|{}".format(os.path.abspath(path), st.st_size, int(st.st_mtime)).encode()
    ).hexdigest()
    cache = os.path.join(_CACHE_DIR, "tarindex_{}.json".format(key))
    if use_cache and os.path.exists(cache):
        try:
            with open(cache) as fh:
                return [Entry(path, o, s, m) for m, o, s in json.load(fh)]
        except (OSError, ValueError):
            pass

    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic in (b"\x1f\x8b", b"BZ", b"\xfd7"):
        raise ValueError(
            "{} is compressed; laespec can only memory-map uncompressed .tar "
            "files. Decompress it first (e.g. `gunzip`).".format(path)
        )

    rows = []
    with tarfile.open(path, "r:") as tf:
        for m in tf:
            if m.isfile():
                rows.append((m.name, m.offset_data, m.size))
    if use_cache:
        try:
            os.makedirs(_CACHE_DIR, exist_ok=True)
            with open(cache, "w") as fh:
                json.dump(rows, fh)
        except OSError:
            pass
    return [Entry(path, o, s, m) for m, o, s in rows]


def build_index(paths: Sequence[str], use_cache: bool = True) -> Dict[str, Entry]:
    """Map data-file basename -> :class:`Entry` for every file reachable
    from ``paths`` (files, directories, and uncompressed .tar archives,
    searched recursively)."""
    index: Dict[str, Entry] = {}

    def add(name, entry):
        base = os.path.basename(name)
        if base in index and index[base].size != entry.size:
            warnings.warn(
                "duplicate file {!r} with different sizes; using the one in {}".format(
                    base, entry.path
                )
            )
        index[base] = entry

    def visit(p):
        if os.path.isdir(p):
            for root, dirs, files in os.walk(p):
                dirs.sort()
                for f in sorted(files):
                    visit(os.path.join(root, f))
        elif p.endswith(".tar"):
            for e in _index_tar(p, use_cache=use_cache):
                add(e.member, e)
        elif os.path.isfile(p):
            add(p, Entry(p, 0, os.path.getsize(p)))

    for p in paths:
        if not os.path.exists(p):
            raise FileNotFoundError(p)
        visit(p)
    return index


# --------------------------------------------------------------------------
# records
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class SpectrumFile:
    """Metadata for one spectrum file, parsed from its name and size."""

    name: str
    model: str
    zstr: str
    skip: int
    n_sightlines: int

    @property
    def z(self) -> float:
        return float(self.zstr)

    @property
    def sampling(self) -> str:
        """``'mass'`` (skip=1: the most massive halos) or ``'uniform'``."""
        return "mass" if self.skip == 1 else "uniform"

    @classmethod
    def from_name(cls, name: str, size: int) -> Optional["SpectrumFile"]:
        m = SPEC_RE.match(os.path.basename(name))
        if m is None:
            return None
        row = N_PIX * SPEC_DTYPE.itemsize
        if size % row:
            warnings.warn("{}: size {} is not a multiple of {} bytes".format(name, size, row))
        return cls(os.path.basename(name), m["model"], m["zstr"], int(m["skip"]), size // row)


@dataclass
class Spectra:
    """A set of sightline spectra from one file.

    Attributes
    ----------
    wavelength : (n_pix,) rest-frame wavelength [A]
    flux : (n_sightlines, n_pix) IGM transmission (a read-only memmap unless
        loaded with ``mmap=False``; call ``np.array(spec.flux)`` for a copy)
    file : the :class:`SpectrumFile` it came from
    halo_index : (n_sightlines,) row of each sightline in the full halo
        catalog at this redshift, or None if the mapping is unknown
    """

    wavelength: np.ndarray
    flux: np.ndarray
    file: SpectrumFile
    halo_index: Optional[np.ndarray] = None
    _halo_loader: Optional[Callable[[], np.ndarray]] = field(default=None, repr=False)
    _halos: Optional[np.ndarray] = field(default=None, repr=False)

    # convenience pass-throughs
    @property
    def z(self) -> float:
        return self.file.z

    @property
    def model(self) -> str:
        return self.file.model

    @property
    def skip(self) -> int:
        return self.file.skip

    def __len__(self):
        return self.flux.shape[0]

    @property
    def halos(self) -> np.ndarray:
        """Structured array (fields x, y, z, vx, vy, vz, mass) of the halo
        hosting each sightline."""
        if self._halos is None:
            if self._halo_loader is None:
                raise LookupError("no halo information available for {}".format(self.file.name))
            self._halos = self._halo_loader()
        return self._halos

    @property
    def mass(self) -> np.ndarray:
        return self.halos["mass"]

    @property
    def velocity(self) -> np.ndarray:
        """Velocity offset from Lyman-alpha line centre [km/s]."""
        from .spectra import velocity_offset

        return velocity_offset(self.wavelength)

    @property
    def observed_wavelength(self) -> np.ndarray:
        """Observed-frame wavelength [A] at the source redshift."""
        return self.wavelength * (1.0 + self.z)

    def mean(self) -> np.ndarray:
        """Mean transmission over sightlines."""
        return np.asarray(self.flux, dtype=np.float64).mean(axis=0)

    def select(self, mask) -> "Spectra":
        """Return the subset of sightlines picked by a boolean mask or index array.

        Example: ``spec.select(spec.mass > 1e11)``.
        """
        mask = np.asarray(mask)
        halos = self._halos[mask] if self._halos is not None else None
        loader = None
        if halos is None and self._halo_loader is not None:
            parent = self._halo_loader
            loader = lambda: parent()[mask]  # noqa: E731
        return Spectra(
            wavelength=self.wavelength,
            flux=self.flux[mask],
            file=self.file,
            halo_index=None if self.halo_index is None else self.halo_index[mask],
            _halo_loader=loader,
            _halos=halos,
        )


# --------------------------------------------------------------------------
# main entry point
# --------------------------------------------------------------------------

class LAEData:
    """Access the LAE spectra and halo catalogs.

    Parameters
    ----------
    *paths : str
        One or more of: ``LAE_Catalogs.tar``, a directory of extracted
        files, or a directory containing the archives from ``laespec
        fetch``. Defaults to ``$LAESPEC_DATA`` if no path is given.

    Examples
    --------
    >>> data = LAEData("LAE_Catalogs.tar")
    >>> data.models
    >>> spec = data.spectra("late_end_early_start_13", z=7.0, sampling="uniform")
    >>> spec.flux.shape, spec.mass[:5]
    """

    def __init__(self, *paths: str, use_cache: bool = True):
        if not paths:
            env = os.environ.get("LAESPEC_DATA")
            if not env:
                raise ValueError("give a data path or set $LAESPEC_DATA")
            paths = tuple(env.split(os.pathsep))
        paths = tuple(os.path.expandvars(os.path.expanduser(str(p))) for p in paths)
        self.paths = paths
        self._index = build_index(paths, use_cache=use_cache)

        self.files: List[SpectrumFile] = []
        self._halo_cats: Dict[str, Entry] = {}
        for base, e in self._index.items():
            sf = SpectrumFile.from_name(base, e.size)
            if sf is not None:
                self.files.append(sf)
                continue
            m = HALO_RE.match(base)
            if m:
                self._halo_cats[m["zstr"]] = e
        self.files.sort(key=lambda f: (f.model, f.z, f.skip))
        if not self.files:
            warnings.warn("no spectrum files found under {}".format(paths))

    # ---- discovery -------------------------------------------------------

    @staticmethod
    def _norm_model(model: str) -> str:
        return model if model.startswith("mono_") else "mono_" + model

    @property
    def models(self) -> List[str]:
        """Reionization-model labels (without the ``mono_`` prefix)."""
        return sorted({f.model[5:] if f.model.startswith("mono_") else f.model for f in self.files})

    def redshifts(self, model: Optional[str] = None, sampling: Optional[str] = None) -> np.ndarray:
        """Sorted redshifts with spectra (optionally for one model / sampling)."""
        fs = self.files
        if model is not None:
            mm = self._norm_model(model)
            fs = [f for f in fs if f.model == mm]
        if sampling is not None:
            fs = [f for f in fs if f.sampling == sampling]
        return np.array(sorted({f.z for f in fs}))

    @property
    def halo_redshifts(self) -> np.ndarray:
        return np.array(sorted(float(z) for z in self._halo_cats))

    def table(self):
        """List of dicts describing every spectrum file (pandas-friendly)."""
        return [
            dict(model=f.model, z=f.z, skip=f.skip, sampling=f.sampling,
                 n_sightlines=f.n_sightlines, name=f.name)
            for f in self.files
        ]

    def summary(self) -> str:
        lines = ["LAE data from: " + ", ".join(self.paths)]
        lines.append("  {} spectrum files, {} halo catalogs, wavelength grid: {}".format(
            len(self.files), len(self._halo_cats),
            "yes" if WAVELENGTH_FILE in self._index else "MISSING"))
        n_matched = sum(1 for b in self._index if b.startswith(MATCHED_PREFIX))
        if n_matched:
            lines.append("  {} matched-halo files".format(n_matched))
        for model in self.models:
            mm = self._norm_model(model)
            fs = [f for f in self.files if f.model == mm]
            zs = sorted({f.z for f in fs})
            nsl = sorted({f.n_sightlines for f in fs})
            lines.append("  {:<32s} {:3d} files, {:3d} redshifts ({:.2f}-{:.2f}), "
                         "sightlines/file: {}".format(model, len(fs), len(zs), zs[0], zs[-1],
                                                      ",".join(map(str, nsl))))
        return "\n".join(lines)

    def __repr__(self):
        return "<LAEData: {} spectrum files, {} models>".format(len(self.files), len(self.models))

    def find(self, model: str, z: float, sampling: str = "uniform",
             skip: Optional[int] = None, tol: float = 0.02) -> SpectrumFile:
        """Return the spectrum file for ``model`` closest to redshift ``z``.

        ``sampling`` is ``'uniform'`` (every k-th halo, k > 1) or ``'mass'``
        (the most massive halos, k = 1); give ``skip`` to pick an exact file.
        Raises ``LookupError`` if nothing lies within ``tol`` of ``z``.
        """
        mm = self._norm_model(model)
        fs = [f for f in self.files if f.model == mm]
        if not fs:
            raise LookupError("unknown model {!r}; available: {}".format(model, self.models))
        if skip is not None:
            fs = [f for f in fs if f.skip == skip]
        elif sampling is not None:
            if sampling not in ("uniform", "mass"):
                raise ValueError("sampling must be 'uniform' or 'mass'")
            fs = [f for f in fs if f.sampling == sampling]
        if not fs:
            raise LookupError("no {} files for {} (skip={})".format(sampling, model, skip))
        best = min(fs, key=lambda f: abs(f.z - z))
        if abs(best.z - z) > tol:
            avail = ", ".join("{:.4f}".format(x) for x in sorted({f.z for f in fs}))
            raise LookupError("no file for {} within dz={} of z={}; available: {}".format(
                model, tol, z, avail))
        return best

    # ---- loading ---------------------------------------------------------

    def _entry(self, base: str) -> Entry:
        try:
            return self._index[base]
        except KeyError:
            raise FileNotFoundError("{} not found under {}".format(base, self.paths)) from None

    def wavelength(self, ascending: bool = True) -> np.ndarray:
        """Rest-frame wavelength grid [A] (ascending by default; the raw
        file is descending)."""
        w = np.array(self._entry(WAVELENGTH_FILE).memmap(SPEC_DTYPE), dtype=np.float64)
        return w[::-1].copy() if ascending and w[0] > w[-1] else w

    def load_flux(self, f: SpectrumFile, mmap: bool = True) -> np.ndarray:
        """Raw transmission array (n_sightlines, 1000) in file (descending-λ) order."""
        return self._entry(f.name).memmap(SPEC_DTYPE, (f.n_sightlines, N_PIX), mmap=mmap)

    def spectra(self, model: str, z: float, sampling: str = "uniform",
                skip: Optional[int] = None, ascending: bool = True,
                mmap: bool = True, tol: float = 0.02) -> Spectra:
        """Load one spectrum file as a :class:`Spectra` object.

        With ``mmap=True`` (default) the flux is a lazy read-only memmap: cheap
        to open, read from disk on access. With ``ascending=True`` wavelength
        and flux are flipped to increasing wavelength (a view, no copy).
        """
        f = self.find(model, z, sampling=sampling, skip=skip, tol=tol)
        wl = self.wavelength(ascending=False)
        flux = self.load_flux(f, mmap=mmap)
        if ascending and wl[0] > wl[-1]:
            wl, flux = wl[::-1].copy(), flux[:, ::-1]

        idx = self._halo_indices(f)
        loader = None
        if (MATCHED_PREFIX + f.name) in self._index or (idx is not None and f.zstr in self._halo_cats):
            loader = lambda: self.matched_halos(f)  # noqa: E731
        return Spectra(wl, flux, f, idx, loader)

    def _halo_indices(self, f: SpectrumFile) -> Optional[np.ndarray]:
        idx = np.arange(f.n_sightlines, dtype=np.int64) * f.skip
        n_cat = self.n_halos(f.zstr) if f.zstr in self._halo_cats else None
        if f.n_sightlines != N_HALO_DEFAULT or (n_cat is not None and idx[-1] >= n_cat):
            return None
        return idx

    def n_halos(self, z) -> int:
        """Number of halos in the catalog at redshift ``z``."""
        return self._halo_entry(z).size // (len(HALO_FIELDS) * HALO_DTYPE.itemsize)

    def _halo_entry(self, z) -> Entry:
        if isinstance(z, str) and z in self._halo_cats:
            return self._halo_cats[z]
        zf = float(z)
        if not self._halo_cats:
            raise LookupError("no halo catalogs found under {}".format(self.paths))
        key = min(self._halo_cats, key=lambda s: abs(float(s) - zf))
        if abs(float(key) - zf) > 1e-3:
            raise LookupError("no halo catalog at z={} (nearest: {})".format(z, key))
        return self._halo_cats[key]

    def halo_catalog(self, z, mmap: bool = True) -> np.ndarray:
        """Full halo catalog at redshift ``z`` as a structured array with
        fields x, y, z [cMpc/h], vx, vy, vz [km/s], mass [Msun/h], sorted by
        descending mass."""
        e = self._halo_entry(z)
        n = e.size // HALO_RECORD.itemsize
        return e.memmap(HALO_RECORD, (n,), mmap=mmap)

    def matched_halos(self, f: SpectrumFile) -> np.ndarray:
        """Halo properties for each sightline of spectrum file ``f``.

        Uses a pre-extracted ``halos_for_*`` file when present (as shipped by
        ``laespec pack``/``fetch``), otherwise slices the full catalog.
        """
        base = MATCHED_PREFIX + f.name
        if base in self._index:
            return np.array(self._index[base].memmap(HALO_RECORD))
        idx = self._halo_indices(f)
        if idx is None:
            raise LookupError(
                "{}: sightline->halo mapping is unknown for files with {} sightlines "
                "(row i = halo i*skip only holds for {}-sightline files)".format(
                    f.name, f.n_sightlines, N_HALO_DEFAULT))
        return np.array(self.halo_catalog(f.zstr)[idx])
