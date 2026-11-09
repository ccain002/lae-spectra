"""Array helpers for manipulating damping-wing transmission spectra.

Every function takes plain numpy arrays: ``wavelength`` is a 1-D rest-frame
grid in Angstrom (ascending or descending), and ``flux`` has the spectral
axis last, i.e. shape ``(n_pix,)`` or ``(n_sightlines, n_pix)``. Use them on
``Spectra.wavelength`` / ``Spectra.flux`` or on any other spectra.
"""

from __future__ import annotations

from typing import Callable, Optional, Sequence, Tuple

import numpy as np

LYA_REST = 1215.67  # Lyman-alpha rest wavelength [A]
C_KMS = 299792.458  # speed of light [km/s]

_trapz = getattr(np, "trapezoid", None) or np.trapz


# ---- coordinates ---------------------------------------------------------

def velocity_offset(wavelength, lambda0: float = LYA_REST):
    """Velocity offset [km/s] of rest-frame ``wavelength`` from ``lambda0``
    (positive = redward)."""
    return C_KMS * (np.asarray(wavelength, dtype=np.float64) / lambda0 - 1.0)


def wavelength_from_velocity(dv, lambda0: float = LYA_REST):
    """Inverse of :func:`velocity_offset`."""
    return lambda0 * (1.0 + np.asarray(dv, dtype=np.float64) / C_KMS)


def to_observed(wavelength, z: float):
    """Rest-frame -> observed-frame wavelength."""
    return np.asarray(wavelength) * (1.0 + z)


def to_rest(wavelength_obs, z: float):
    """Observed-frame -> rest-frame wavelength."""
    return np.asarray(wavelength_obs) / (1.0 + z)


def ensure_ascending(wavelength, flux):
    """Return ``(wavelength, flux)`` ordered by increasing wavelength (views)."""
    wavelength = np.asarray(wavelength)
    if wavelength[0] > wavelength[-1]:
        return wavelength[::-1], np.asarray(flux)[..., ::-1]
    return wavelength, np.asarray(flux)


# ---- statistics ----------------------------------------------------------

def mean_transmission(flux, weights=None, axis: int = 0):
    """(Weighted) mean transmission over sightlines (in float64)."""
    flux = np.asarray(flux, dtype=np.float64)
    if weights is None:
        return flux.mean(axis=axis)
    return np.average(flux, axis=axis, weights=weights)


def transmission_percentiles(flux, q: Sequence[float] = (16, 50, 84), axis: int = 0):
    """Percentiles of the transmission over sightlines; returns shape (len(q), n_pix)."""
    return np.percentile(np.asarray(flux, dtype=np.float64), q, axis=axis)


def optical_depth(flux, floor: float = 1e-30):
    """tau = -ln(T), with T clipped at ``floor`` to keep it finite."""
    return -np.log(np.clip(np.asarray(flux, dtype=np.float64), floor, None))


def bootstrap(flux, n_boot: int = 200, statistic: Callable = mean_transmission,
              seed: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Bootstrap over sightlines. Returns ``(statistic(flux), std over resamples)``."""
    flux = np.asarray(flux)
    rng = np.random.default_rng(seed)
    n = flux.shape[0]
    reps = np.array([statistic(flux[rng.integers(0, n, n)]) for _ in range(n_boot)])
    return statistic(flux), reps.std(axis=0)


def stack(flux, values, bins, statistic: Callable = mean_transmission):
    """Stack sightlines in bins of a per-sightline quantity (e.g. halo mass).

    Returns ``(stacked, counts)`` where ``stacked`` has shape
    ``(len(bins) - 1, n_pix)`` (NaN for empty bins).

    Example: ``stack(spec.flux, np.log10(spec.mass), [9, 10, 11, 12])``.
    """
    values = np.asarray(values)
    flux = np.asarray(flux)
    which = np.digitize(values, bins) - 1
    nb = len(bins) - 1
    out = np.full((nb, flux.shape[-1]), np.nan)
    counts = np.zeros(nb, dtype=int)
    for i in range(nb):
        sel = which == i
        counts[i] = sel.sum()
        if counts[i]:
            out[i] = statistic(flux[sel])
    return out, counts


# ---- resampling ----------------------------------------------------------

def interpolate(wavelength, flux, new_wavelength, fill_value=np.nan):
    """Linearly interpolate ``flux`` onto ``new_wavelength`` (vectorised over
    sightlines). Points outside the grid get ``fill_value``."""
    wl, fl = ensure_ascending(wavelength, flux)
    new = np.asarray(new_wavelength, dtype=np.float64)
    i = np.clip(np.searchsorted(wl, new, side="right") - 1, 0, len(wl) - 2)
    t = (new - wl[i]) / (wl[i + 1] - wl[i])
    out = fl[..., i] * (1.0 - t) + fl[..., i + 1] * t
    outside = (new < wl[0]) | (new > wl[-1])
    if np.any(outside):
        out = np.where(outside, fill_value, out)
    return out


def transmission_at(wavelength, flux, at_wavelength=None, dv: Optional[float] = None):
    """Transmission of each sightline at one rest wavelength [A] or at
    velocity offset ``dv`` [km/s] from Lyman-alpha."""
    if (at_wavelength is None) == (dv is None):
        raise ValueError("give exactly one of at_wavelength or dv")
    if dv is not None:
        at_wavelength = wavelength_from_velocity(dv)
    return interpolate(wavelength, flux, np.atleast_1d(at_wavelength))[..., 0]


def rebin(wavelength, flux, factor: int):
    """Average ``factor`` adjacent pixels together (trailing pixels dropped)."""
    wavelength = np.asarray(wavelength)
    flux = np.asarray(flux)
    n = (len(wavelength) // factor) * factor
    wl = wavelength[:n].reshape(-1, factor).mean(axis=1)
    fl = flux[..., :n].reshape(flux.shape[:-1] + (-1, factor)).mean(axis=-1)
    return wl, fl


def smooth(wavelength, flux, fwhm_kms: Optional[float] = None,
           resolution: Optional[float] = None):
    """Convolve with a Gaussian line-spread function.

    Give either ``fwhm_kms`` or a spectral resolving power ``resolution``
    (R = lambda / dlambda, i.e. FWHM = c / R). The kernel width is set in
    velocity at the grid's mean pixel spacing (the grid is uniform in
    wavelength, so this is accurate to a few per cent across 1210-1316 A).
    Edges are padded with the end values.
    """
    if (fwhm_kms is None) == (resolution is None):
        raise ValueError("give exactly one of fwhm_kms or resolution")
    if resolution is not None:
        fwhm_kms = C_KMS / resolution
    wavelength = np.asarray(wavelength, dtype=np.float64)
    flux = np.asarray(flux, dtype=np.float64)
    dv_pix = np.mean(np.abs(np.diff(velocity_offset(wavelength))))
    sigma = fwhm_kms / (2.0 * np.sqrt(2.0 * np.log(2.0))) / dv_pix
    if sigma < 0.1:
        return flux.copy()
    half = int(np.ceil(4 * sigma))
    k = np.exp(-0.5 * (np.arange(-half, half + 1) / sigma) ** 2)
    k /= k.sum()
    pad = [(0, 0)] * (flux.ndim - 1) + [(half, half)]
    fp = np.pad(flux, pad, mode="edge")
    n = flux.shape[-1]
    out = np.zeros_like(flux)
    for j, w in enumerate(k):
        out += w * fp[..., j:j + n]
    return out


# ---- Lyman-alpha line transmission ---------------------------------------

def line_profile(wavelength, dv: float = 0.0, sigma_kms: float = 100.0, fwhm_kms=None):
    """Gaussian intrinsic Lyman-alpha emission line on ``wavelength``,
    centred ``dv`` km/s from line centre and normalised so that its
    integral over wavelength is 1."""
    if fwhm_kms is not None:
        sigma_kms = fwhm_kms / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    v = velocity_offset(wavelength)
    prof = np.exp(-0.5 * ((v - dv) / sigma_kms) ** 2)
    norm = abs(_trapz(prof, wavelength))
    if norm == 0:
        raise ValueError("line profile falls entirely outside the wavelength grid")
    return prof / norm


def lya_transmission(wavelength, flux, profile=None, dv: float = 0.0,
                     sigma_kms: float = 100.0):
    """Fraction of an intrinsic Lyman-alpha line transmitted by the IGM,
    T_Lya = int T(lambda) L(lambda) dlambda / int L(lambda) dlambda,
    for each sightline.

    ``profile`` is the intrinsic line on the same grid; if omitted a
    Gaussian with offset ``dv`` and width ``sigma_kms`` is used (see
    :func:`line_profile`). Note the grid starts at ~1210.8 A (-1200 km/s),
    so very blue lines are truncated.
    """
    wavelength = np.asarray(wavelength, dtype=np.float64)
    if profile is None:
        profile = line_profile(wavelength, dv=dv, sigma_kms=sigma_kms)
    profile = np.asarray(profile, dtype=np.float64)
    num = _trapz(np.asarray(flux, dtype=np.float64) * profile, wavelength, axis=-1)
    return num / _trapz(profile, wavelength)
