"""laespec: read and manipulate LAE Lyman-alpha damping-wing spectra.

Quick start::

    from laespec import LAEData
    data = LAEData("/path/to/LAE_Catalogs.tar")   # or a directory of archives/files
    spec = data.spectra("late_end_early_start_13", z=7.0, sampling="uniform")
    spec.wavelength, spec.flux, spec.mass
"""

from .catalog import (HALO_FIELDS, LAEData, Spectra, SpectrumFile)  # noqa: F401
from .spectra import (C_KMS, LYA_REST, bootstrap, ensure_ascending,  # noqa: F401
                      interpolate, line_profile, lya_transmission, mean_transmission,
                      optical_depth, rebin, smooth, stack, to_observed, to_rest,
                      transmission_at, transmission_percentiles, velocity_offset,
                      wavelength_from_velocity)
from .distribute import fetch  # noqa: F401

__version__ = "0.1.0"
