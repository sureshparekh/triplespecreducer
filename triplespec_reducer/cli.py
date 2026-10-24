"""
The command line.

    triplespec-reducer                         the desktop interface
    triplespec-reducer gui [RAW]               the same, with a raw folder filled in
    triplespec-reducer organize RAW OUT        read a night and write OUT/night.json
    triplespec-reducer reduce OUT [SOURCE ...] run every step from OUT/night.json
    triplespec-reducer products OUT SOURCE     list what exists for SOURCE
"""

import argparse
import json
import os
import sys
import time


def _progress(t0):
    return lambda step, msg="": print(f"[{time.time() - t0:7.1f} s] {step:9s} {msg}", flush=True)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] == "gui":
        from triplespec_reducer.gui.launch import main as gui
        return gui([sys.argv[0]] + argv[1:2])
    p = argparse.ArgumentParser(prog="triplespec-reducer",
                                description="SOAR TripleSpec 4.1 reduction (Spextool TS4_SOAR, in Python)")
    sub = p.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("organize", help="read a night's headers and propose the plan")
    o.add_argument("raw")
    o.add_argument("out")
    o.add_argument("--no-simbad", action="store_true", help="do not look the standards up on SIMBAD")
    r = sub.add_parser("reduce", help="run the steps on a plan")
    r.add_argument("out", help="the output folder holding night.json")
    r.add_argument("sources", nargs="*", help="science targets (default: all)")
    r.add_argument("--steps", default="flat,wavecal,extract,combine,telluric,merge")
    r.add_argument("--force", action="store_true", help="redo steps whose products exist")
    r.add_argument("--idl-defaults", action="store_true",
                   help="every setting exactly as the IDL panels start (kernel continuum degree 1)")
    r.add_argument("--xs-background", action="store_true",
                   help="subtract the extended-source background (the IDL tool never does)")
    r.add_argument("--combine-images", action="store_true",
                   help="average all pairs into one image before extracting")
    r.add_argument("--options", help="a JSON file of Options fields to set")
    q = sub.add_parser("products", help="every file made for a source")
    q.add_argument("out")
    q.add_argument("source")
    a = p.parse_args(argv)

    from triplespec_reducer.night import organize, Night
    if a.cmd == "organize":
        t0 = time.time()
        n = organize(a.raw, a.out, lookup=not a.no_simbad, progress=lambda m: print(m))
        path = n.save()
        print(f"{len(n.frames)} frames, {len(n.sequences)} targets  ({time.time() - t0:.1f} s)")
        for s in n.sequences:
            print(f"  {s.name:34s} {s.role:9s} {s.source:9s} {s.pattern:10s} pairs {len(s.pairs):3d}  "
                  f"airmass {s.airmass:.2f}  wavecal {s.wavecal_from}  standard {s.standard or '-'}")
            for note in s.notes:
                print(f"      {note}")
        for w in n.warnings:
            print("  warning:", w)
        print("plan:", path)
        return 0
    n = Night.load(os.path.join(a.out, "night.json"))
    if a.cmd == "products":
        for step, path in n.products(a.source).items():
            print(f"{step:12s} {path}")
        return 0
    from triplespec_reducer.pipeline import Reducer, Options
    opts = Options.idl_defaults() if a.idl_defaults else Options()
    if a.xs_background:
        opts.xs_background = True
    if a.combine_images:
        opts.combine_images = True
    if a.options:
        with open(a.options) as fh:
            for k, v in json.load(fh).items():
                setattr(opts, k, v)
    t0 = time.time()
    red = Reducer(n, opts, progress=_progress(t0), force=a.force)
    prods = red.run(sources=a.sources or None, steps=tuple(a.steps.split(",")))
    for src, files in prods.items():
        print(src)
        for step, path in files.items():
            print(f"   {step:12s} {path}")
    return 1 if red.failures else 0


if __name__ == "__main__":
    sys.exit(main())
