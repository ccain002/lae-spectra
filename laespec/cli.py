"""Command-line interface: ``laespec {info,fetch,pack,set-urls}``."""

from __future__ import annotations

import argparse
import sys


def main(argv=None):
    p = argparse.ArgumentParser(prog="laespec", description="LAE damping-wing spectra tools")
    sub = p.add_subparsers(dest="cmd")

    s = sub.add_parser("info", help="summarise the data under one or more paths")
    s.add_argument("paths", nargs="+", help="LAE_Catalogs.tar, data directories or archives")
    s.add_argument("--files", action="store_true", help="list every spectrum file")

    s = sub.add_parser("fetch", help="download the shared archives")
    s.add_argument("--dest", required=True, help="download directory")
    s.add_argument("--manifest", help="manifest path or URL (default: bundled manifest)")
    s.add_argument("--model", action="append", help="only this model (repeatable; default all)")
    s.add_argument("--no-matched-halos", action="store_true")
    s.add_argument("--halo-catalogs", action="store_true", help="also get full halo catalogs (large)")
    s.add_argument("--list", action="store_true", help="show what would be downloaded and exit")
    s.add_argument("--no-verify", action="store_true", help="skip SHA-256 check")

    s = sub.add_parser("pack", help="(maintainer) split the data into shareable archives")
    s.add_argument("sources", nargs="+")
    s.add_argument("--out", required=True)
    s.add_argument("--model", action="append")
    s.add_argument("--full-halos", action="store_true", help="also package full halo catalogs")
    s.add_argument("--split-gb", type=float, default=20.0, help="max size of each halo-catalog part")
    s.add_argument("--base-url", help="download URL prefix or template containing {name}")

    s = sub.add_parser("set-urls", help="(maintainer) fill in download URLs in a manifest")
    s.add_argument("manifest")
    s.add_argument("--base-url", required=True)

    a = p.parse_args(argv)
    if a.cmd == "info":
        from .catalog import LAEData

        d = LAEData(*a.paths)
        print(d.summary())
        if a.files:
            for r in d.table():
                print("  {model:<36s} z={z:7.4f} skip={skip:<5d} n={n_sightlines:<6d} {name}".format(**r))
    elif a.cmd == "fetch":
        from . import distribute

        if a.list:
            man = distribute.load_manifest(a.manifest)
            todo = distribute.select(man, a.model, not a.no_matched_halos, a.halo_catalogs)
            for e in todo:
                print("{:<45s} {:7.2f} GB".format(e["name"], e["size"] / 1e9))
            print("total {:.2f} GB".format(sum(e["size"] for e in todo) / 1e9))
            return
        distribute.fetch(a.dest, a.manifest, a.model, not a.no_matched_halos,
                         a.halo_catalogs, verify=not a.no_verify)
    elif a.cmd == "pack":
        from . import distribute

        distribute.pack(a.sources, a.out, a.model, full_halos=a.full_halos,
                        split_gb=a.split_gb, base_url=a.base_url)
    elif a.cmd == "set-urls":
        from . import distribute

        distribute.set_urls(a.manifest, a.base_url)
        print("updated " + a.manifest)
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
