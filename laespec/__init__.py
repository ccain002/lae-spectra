"""laespec: read and manipulate LAE Lyman-alpha damping-wing spectra.

Quick start::

    import laespec as ls
    spec = ls.example().spectra("late_end_early_start_13", z=7.0)  # bundled sample, no download

    from laespec import LAEData
    data = LAEData("/path/to/LAE_Catalogs.tar")   # or a directory of archives/files
    spec = data.spectra("late_end_early_start_13", z=7.0, sampling="uniform")
    spec.wavelength, spec.flux, spec.mass
"""

from .catalog import (EXAMPLE_DATA, EXAMPLE_MODEL, EXAMPLE_Z, HALO_FIELDS,  # noqa: F401
                      LAEData, Spectra, SpectrumFile, example)
from .spectra import (C_KMS, LYA_REST, bootstrap, ensure_ascending,  # noqa: F401
                      interpolate, line_profile, lya_transmission, mean_transmission,
                      optical_depth, rebin, smooth, stack, to_observed, to_rest,
                      transmission_at, transmission_percentiles, velocity_offset,
                      wavelength_from_velocity)
from .distribute import fetch  # noqa: F401

__version__ = "0.1.0"
