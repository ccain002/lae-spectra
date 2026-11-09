"""Package the data for sharing (``pack``) and download it again (``fetch``).

``pack`` splits the monolithic ``LAE_Catalogs.tar`` into

* ``LAE_spectra_<model>.tar``      - one per reionization model, with the
  wavelength grid included so each archive is usable on its own;
* ``LAE_matched_halos.tar``        - for every spectrum file, just the halo
  rows matched to its sightlines (~0.3 GB in total instead of ~26 GB of
  full halo catalogs);
* ``LAE_halo_catalogs[_partN].tar`` - optional, the full halo catalogs;

and writes ``manifest.json`` with the size, SHA-256 and download URL of
each. Archives are uncompressed tars, which :class:`laespec.LAEData` reads
in place, so collaborators never need to extract them.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import sys
import tarfile
import time
import urllib.request
from typing import Iterable, List, Optional

from .catalog import (HALO_RE, MATCHED_PREFIX, N_HALO_DEFAULT, WAVELENGTH_FILE,
                      LAEData)

MANIFEST_NAME = "manifest.json"
_DEFAULT_MANIFEST = os.path.join(os.path.dirname(__file__), MANIFEST_NAME)
_CHUNK = 8 << 20


# ---------------------------------------------------------------------------
# pack
# ---------------------------------------------------------------------------

class _HashingWriter:
    def __init__(self, fh):
        self.fh, self.sha, self.n = fh, hashlib.sha256(), 0

    def write(self, b):
        self.sha.update(b)
        self.n += len(b)
        return self.fh.write(b)


def _write_tar(path, members, log):
    """Write an uncompressed tar from (arcname, entry_or_bytes) pairs, hashing
    as we go. Returns (size, sha256)."""
    tmp = path + ".part"
    with open(tmp, "wb") as raw:
        hw = _HashingWriter(raw)
        with tarfile.open(fileobj=hw, mode="w|", format=tarfile.GNU_FORMAT) as tf:
            for arcname, src in members:
                ti = tarfile.TarInfo(arcname)
                ti.mtime = int(time.time())
                ti.mode = 0o644
                if isinstance(src, bytes):
                    ti.size = len(src)
                    tf.addfile(ti, io.BytesIO(src))
                else:
                    ti.size = src.size
                    with open(src.path, "rb") as fh:
                        fh.seek(src.offset)
                        tf.addfile(ti, fh)
    os.replace(tmp, path)
    log("  wrote {} ({:.2f} GB)".format(os.path.basename(path), hw.n / 1e9))
    return hw.n, hw.sha.hexdigest()


def _url(base_url, name):
    if not base_url:
        return None
    if "{name}" in base_url:
        return base_url.format(name=name)
    return base_url.rstrip("/") + "/" + name


def pack(sources: Iterable[str], outdir: str, models: Optional[List[str]] = None,
         matched_halos: bool = True, full_halos: bool = False, split_gb: float = 20.0,
         base_url: Optional[str] = None, log=print) -> dict:
    """Split the data into shareable archives + ``manifest.json`` (see module doc).

    ``base_url`` is where the archives will be downloadable from, e.g.
    ``https://zenodo.org/records/1234567/files/{name}?download=1``; it can also
    be filled in later with :func:`set_urls`.
    """
    data = LAEData(*sources)
    os.makedirs(outdir, exist_ok=True)
    wl = data._entry(WAVELENGTH_FILE)
    entries = []

    for model in models or data.models:
        fs = [f for f in data.files if f.model == data._norm_model(model)]
        name = "LAE_spectra_{}.tar".format(model)
        log("packing {} ({} files)".format(name, len(fs)))
        members = [("LAE_Spectra/" + WAVELENGTH_FILE, wl)]
        members += [("LAE_Spectra/" + f.name, data._entry(f.name)) for f in fs]
        size, sha = _write_tar(os.path.join(outdir, name), members, log)
        entries.append(dict(name=name, kind="spectra", model=model, size=size, sha256=sha))

    if matched_halos:
        members, skipped = [], []
        for f in data.files:
            if models and f.model[5:] not in models and f.model not in models:
                continue
            try:
                members.append(("matched_halos/" + MATCHED_PREFIX + f.name,
                                data.matched_halos(f).tobytes()))
            except LookupError:
                skipped.append(f.name)
        if members:
            name = "LAE_matched_halos.tar"
            log("packing {} ({} files)".format(name, len(members)))
            size, sha = _write_tar(os.path.join(outdir, name), members, log)
            entries.append(dict(name=name, kind="matched_halos", size=size, sha256=sha))
        for s in skipped:
            log("  no halo match for {} (not {} sightlines or no catalog)".format(s, N_HALO_DEFAULT))

    if full_halos:
        cats = sorted(((z, e) for z, e in data._halo_cats.items()), key=lambda t: -float(t[0]))
        parts, cur, cur_size = [], [], 0
        for z, e in cats:
            if cur and cur_size + e.size > split_gb * 1e9:
                parts.append(cur)
                cur, cur_size = [], 0
            cur.append(("halo_catalogs/halo_xvm_catalog.z=" + z, e))
            cur_size += e.size
        if cur:
            parts.append(cur)
        for i, members in enumerate(parts):
            name = ("LAE_halo_catalogs.tar" if len(parts) == 1
                    else "LAE_halo_catalogs_part{}.tar".format(i + 1))
            log("packing {} ({} catalogs)".format(name, len(members)))
            size, sha = _write_tar(os.path.join(outdir, name), members, log)
            zs = [HALO_RE.match(os.path.basename(a))["zstr"] for a, _ in members]
            entries.append(dict(name=name, kind="halo_catalogs", size=size, sha256=sha,
                                z_range=[float(min(zs)), float(max(zs))]))

    for e in entries:
        e["url"] = _url(base_url, e["name"])
    manifest = dict(dataset="LAE Lyman-alpha damping-wing spectra",
                    created=time.strftime("%Y-%m-%d"), files=entries)
    with open(os.path.join(outdir, MANIFEST_NAME), "w") as fh:
        json.dump(manifest, fh, indent=2)
    log("wrote " + os.path.join(outdir, MANIFEST_NAME))
    return manifest


def set_urls(manifest_path: str, base_url: str) -> dict:
    """Fill in download URLs in an existing manifest (after uploading)."""
    with open(manifest_path) as fh:
        man = json.load(fh)
    for e in man["files"]:
        e["url"] = _url(base_url, e["name"])
    with open(manifest_path, "w") as fh:
        json.dump(man, fh, indent=2)
    return man


# ---------------------------------------------------------------------------
# fetch
# ---------------------------------------------------------------------------

def load_manifest(src: Optional[str] = None) -> dict:
    """Read a manifest from a path or URL (default: the one bundled with laespec)."""
    src = src or os.environ.get("LAESPEC_MANIFEST") or _DEFAULT_MANIFEST
    if src.startswith(("http://", "https://", "file://")):
        with urllib.request.urlopen(src) as r:
            return json.loads(r.read().decode())
    with open(src) as fh:
        return json.load(fh)


def select(manifest: dict, models: Optional[List[str]] = None, matched_halos: bool = True,
           halo_catalogs: bool = False) -> List[dict]:
    out = []
    for e in manifest["files"]:
        k = e["kind"]
        if k == "spectra" and (not models or e["model"] in models):
            out.append(e)
        elif k == "matched_halos" and matched_halos:
            out.append(e)
        elif k == "halo_catalogs" and halo_catalogs:
            out.append(e)
    return out


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(_CHUNK), b""):
            h.update(b)
    return h.hexdigest()


def _download(url, dest, size, log):
    part = dest + ".part"
    have = os.path.getsize(part) if os.path.exists(part) else 0
    req = urllib.request.Request(url, headers={"User-Agent": "laespec"})
    if have:
        req.add_header("Range", "bytes={}-".format(have))
    with urllib.request.urlopen(req) as r:
        if have and getattr(r, "status", 200) != 206:
            have = 0  # server ignored the Range header; start over
        mode = "ab" if have else "wb"
        done, t0, last = have, time.time(), 0.0
        with open(part, mode) as out:
            for b in iter(lambda: r.read(_CHUNK), b""):
                out.write(b)
                done += len(b)
                now = time.time()
                if now - last > 2 and sys.stderr.isatty():
                    rate = (done - have) / max(now - t0, 1e-9) / 1e6
                    pct = 100.0 * done / size if size else 0
                    sys.stderr.write("\r  {:5.1f}%  {:8.1f} MB/s".format(pct, rate))
                    last = now
    if sys.stderr.isatty():
        sys.stderr.write("\n")
    return part


def fetch(dest: str, manifest: Optional[str] = None, models: Optional[List[str]] = None,
          matched_halos: bool = True, halo_catalogs: bool = False, verify: bool = True,
          log=print) -> List[str]:
    """Download the selected archives into ``dest`` (resuming partial
    downloads, skipping ones already present) and check their SHA-256.

    Returns the list of local archive paths; pass ``dest`` to
    :class:`laespec.LAEData` to use them.
    """
    dest = os.path.expanduser(dest)
    man = load_manifest(manifest)
    todo = select(man, models, matched_halos, halo_catalogs)
    if not todo:
        raise ValueError("nothing to download: the manifest lists no matching files "
                         "(is it still the empty placeholder?)")
    missing = [e["name"] for e in todo if not e.get("url")]
    if missing:
        raise ValueError("manifest has no URL for: " + ", ".join(missing))
    os.makedirs(dest, exist_ok=True)
    total = sum(e["size"] for e in todo)
    log("{} archive(s), {:.1f} GB total -> {}".format(len(todo), total / 1e9, dest))
    paths = []
    for e in todo:
        path = os.path.join(dest, e["name"])
        if os.path.exists(path) and os.path.getsize(path) == e["size"]:
            log("  {} already present".format(e["name"]))
            paths.append(path)
            continue
        log("  downloading {} ({:.2f} GB)".format(e["name"], e["size"] / 1e9))
        part = _download(e["url"], path, e["size"], log)
        if os.path.getsize(part) != e["size"]:
            raise IOError("{}: got {} bytes, expected {}; rerun to resume".format(
                e["name"], os.path.getsize(part), e["size"]))
        if verify and _sha256(part) != e["sha256"]:
            os.remove(part)
            raise IOError("{}: checksum mismatch (download corrupted); rerun".format(e["name"]))
        os.replace(part, path)
        paths.append(path)
    log("done. Use:  LAEData({!r})".format(dest))
    return paths
