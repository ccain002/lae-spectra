import os
import tarfile

import numpy as np
import pytest

import laespec
from laespec import LAEData
from laespec import catalog, distribute
from laespec import spectra as sp

N = catalog.N_HALO_DEFAULT


@pytest.fixture(autouse=True)
def _cache(tmp_path, monkeypatch):
    monkeypatch.setattr(catalog, "_CACHE_DIR", str(tmp_path / "cache"))


def _make_dataset(root, n_cat=3 * N):
    """Tiny fake dataset with the real file layout."""
    spec_dir = root / "LAE_Catalogs" / "LAE_Spectra"
    halo_dir = root / "LAE_Catalogs" / "halo_catalogs"
    spec_dir.mkdir(parents=True)
    halo_dir.mkdir()
    wl = np.linspace(1315.67, 1210.775, catalog.N_PIX).astype("<f4")  # descending, like the real file
    wl.tofile(str(spec_dir / "wav_out_LAE"))
    halos = np.zeros((n_cat, 7))
    halos[:, 6] = np.logspace(12, 8, n_cat)  # descending mass
    halos[:, 0] = np.arange(n_cat)  # x = row number, to check the mapping
    halos.tofile(str(halo_dir / "halo_xvm_catalog.z=07.0445"))
    np.zeros(0).tofile(str(halo_dir / "halo_xvm_catalog.z=28.1959"))  # empty, as in the real data
    rng = np.random.default_rng(0)
    for model in ("late_end_early_start_13", "late_end_late_start_4"):
        for skip in (1, 3):
            f = rng.uniform(0, 1, (N, catalog.N_PIX)).astype("<f4")
            f[:, 0] = np.arange(N)  # tag rows (first pixel = reddest in file order)
            f.tofile(str(spec_dir / "spec_out_mono_{}.z=07.0445_skip{}_LAE".format(model, skip)))
    # a file whose halo mapping is unknown (20x the sightlines)
    np.ones((2 * N, catalog.N_PIX), "<f4").tofile(
        str(spec_dir / "spec_out_mono_early_end_early_start_3.z=07.0445_skip3_LAE"))
    return wl


@pytest.fixture
def dataset(tmp_path):
    root = tmp_path / "data"
    root.mkdir()
    wl = _make_dataset(root)
    return root, wl


@pytest.fixture
def tarred(dataset, tmp_path):
    root, wl = dataset
    tpath = tmp_path / "LAE_Catalogs.tar"
    with tarfile.open(str(tpath), "w") as tf:
        tf.add(str(root / "LAE_Catalogs"), arcname="LAE_Catalogs")
    return tpath, wl


@pytest.mark.parametrize("kind", ["dir", "tar"])
def test_discovery_and_loading(kind, dataset, tarred):
    path = dataset[0] if kind == "dir" else tarred[0]
    wl_raw = dataset[1]
    d = LAEData(str(path))
    assert d.models == ["early_end_early_start_3", "late_end_early_start_13", "late_end_late_start_4"]
    assert np.allclose(d.redshifts("late_end_late_start_4"), [7.0445])
    assert len(d.files) == 5
    assert d.n_halos(7.0445) == 3 * N
    assert len(d.halo_catalog(28.1959)) == 0

    w = d.wavelength()
    assert w[0] < w[-1] and np.allclose(w, wl_raw[::-1])

    s = d.spectra("late_end_early_start_13", 7.04, sampling="uniform")
    assert s.skip == 3 and s.flux.shape == (N, catalog.N_PIX)
    assert np.allclose(s.flux[:, -1], np.arange(N))  # flipped: file pixel 0 is now last
    assert np.array_equal(s.halo_index, np.arange(N) * 3)
    assert np.array_equal(s.halos["x"], np.arange(N) * 3)
    assert s.mass[0] == pytest.approx(1e12)

    sub = s.select(s.mass > 1e11)
    assert len(sub) == np.sum(s.mass > 1e11)
    assert np.all(sub.mass > 1e11)

    raw = d.spectra("mono_late_end_early_start_13", 7.04, sampling="mass", ascending=False, mmap=False)
    assert raw.skip == 1 and np.allclose(raw.flux[:, 0], np.arange(N))

    big = d.spectra("early_end_early_start_3", 7.04)
    assert big.flux.shape[0] == 2 * N and big.halo_index is None
    with pytest.raises(LookupError):
        big.halos

    with pytest.raises(LookupError):
        d.find("late_end_early_start_13", 8.0)
    with pytest.raises(LookupError):
        d.find("nope", 7.0)


def test_pack_fetch_roundtrip(tarred, tmp_path):
    tpath, _ = tarred
    out = tmp_path / "packed"
    man = distribute.pack([str(tpath)], str(out), full_halos=True, split_gb=1e-9,
                          base_url="file://" + str(out), log=lambda *a: None)
    kinds = sorted(e["kind"] for e in man["files"])
    assert kinds.count("spectra") == 3 and "matched_halos" in kinds
    assert kinds.count("halo_catalogs") == 2  # split into parts

    dest = tmp_path / "fetched"
    distribute.fetch(str(dest), manifest=str(out / "manifest.json"),
                     models=["late_end_early_start_13"], log=lambda *a: None)
    got = sorted(os.listdir(str(dest)))
    assert got == ["LAE_matched_halos.tar", "LAE_spectra_late_end_early_start_13.tar"]

    # collaborator view: spectra + matched halos only, no full catalogs
    d = LAEData(str(dest))
    s = d.spectra("late_end_early_start_13", 7.0445)
    assert np.array_equal(s.halos["x"], np.arange(N) * 3)
    full = LAEData(str(tpath))
    assert np.array_equal(full.spectra("late_end_early_start_13", 7.0445).flux, s.flux)


def test_spectra_helpers():
    wl = np.linspace(1210.775, 1315.67, 1000)
    assert sp.velocity_offset(sp.LYA_REST) == 0
    assert sp.wavelength_from_velocity(sp.velocity_offset(1250.0)) == pytest.approx(1250.0)
    flux = np.vstack([np.full(1000, 0.5), (wl > 1216).astype(float)])

    assert np.allclose(sp.mean_transmission(flux)[wl > 1216], 0.75)
    assert np.allclose(sp.interpolate(wl[::-1], flux[:, ::-1], [1250.0]), [[0.5], [1.0]])
    assert np.isnan(sp.interpolate(wl, flux, [1000.0])).all()
    assert np.allclose(sp.transmission_at(wl, flux, dv=1000.0), [0.5, 1.0])

    sm = sp.smooth(wl, flux, fwhm_kms=300.0)
    assert sm.shape == flux.shape and np.allclose(sm[0], 0.5)
    assert np.allclose(sm[1, [0, -1]], [0, 1])

    w2, f2 = sp.rebin(wl, flux, 3)
    assert f2.shape == (2, 333) and len(w2) == 333

    prof = sp.line_profile(wl, dv=300.0, sigma_kms=50.0)
    assert np.trapz(prof, wl) == pytest.approx(1.0)
    T = sp.lya_transmission(wl, flux, dv=300.0, sigma_kms=50.0)
    assert np.allclose(T, [0.5, 1.0], atol=1e-6)

    st, counts = sp.stack(flux, [9.5, 10.5], [9, 10, 11, 12])
    assert list(counts) == [1, 1, 0] and np.isnan(st[2]).all()

    m, err = sp.bootstrap(flux, n_boot=20, seed=1)
    assert m.shape == err.shape == (1000,)
    assert laespec.optical_depth(np.array([1.0]))[0] == 0
