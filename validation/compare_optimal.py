"""
Point-source extraction as xspextool's Do All Steps makes it for TS4_SOAR
(optimal extraction, the panel's other find modes and options), compared
with IDL (validation/idl/ref_optimal.pro, and ref_1d.pro for the 1D extraction made
without a wavecal file). Every case goes through
extract.reduce_pair with the panel's inputs; IDL's flat and wavecal are
used, so only extraction is being tested.

    python validation/compare_optimal.py [case ...]
"""

import os
import sys
import tempfile

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import flat as F, wavecal as W, extract as E   # noqa: E402
from triplespec_reducer.instrument import read_atran                    # noqa: E402

RAW = os.path.join(os.path.dirname(HERE), "data", "raw", "2025-11-02")
REF = os.path.join(HERE, "idl_out")


def raw(n):
    return os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits")


# case: (frames, mode, reduce_pair settings)
CASES = {
    "op": ([116, 117], "A-B", dict()),
    "om": ([116, 117], "A-B", dict(medprofile=True)),
    "nb": ([116, 117], "A-B", dict(psbgsub=False)),
    "gs": ([116, 117], "A-B", dict(findmode="guess", guess=[9.5, 28.3])),
    "fx": ([116, 117], "A-B", dict(findmode="fixed", guess=[9.5, 28.3])),
    "it": ([116, 117], "A-B", dict(trace_in=os.path.join(REF, "op_trace.dat"))),
    "os": ([116, 117], "A-B", dict(ps_orders=[3, 5, 6])),
    "am": ([116], "A", dict(findmode="fixed", guess=[9.8, 20.0])),
    # no wavecal file: MC_EXTPSSPEC / MC_EXTXSSPEC (validation/idl/ref_1d.pro)
    "p1o": ([116, 117], "A-B", dict(one_d=True)),
    "p1s": ([116, 117], "A-B", dict(one_d=True, psbgsub=False)),
    "x1f": ([120, 121], "A-B", dict(one_d=True, source="extended", xs_positions=[20.1],
                                    xs_radii=[1.0], xsbg=((2., 6.), (32., 36.)), xsbgdeg=0)),
    "x1n": ([120, 121], "A-B", dict(one_d=True, source="extended", xs_positions=[20.1],
                                    xs_radii=[1.0], psbgsub=False, xsbgsub=False)),
}


def run(case, fl, wc, atran, outdir):
    frames, mode, kw = CASES[case]
    kw = dict(kw)
    if kw.pop("one_d", False):
        wc = W.no_wavecal(fl)                 # the Full Wavecal Name field left blank
    op = [os.path.join(outdir, f"{case}_{k}.fits") for k in range(2)]
    return E.reduce_pair([raw(n) for n in frames], fl, wc, op, atran=atran, mode=mode, **kw)


def trace_text(case, fl, wc, atran, outdir):
    """The trace file the case stores, against IDL's text."""
    path = os.path.join(outdir, f"{case}_trace.dat")
    frames, mode, k = CASES[case]
    r = E.reduce_pair([raw(n) for n in frames], fl, wc, [os.path.join(outdir, "t.fits")] * 2,
                      atran=atran, mode=mode, trace_out=path, **k)
    return open(path).read() == open(os.path.join(REF, f"{case}_trace.dat")).read(), r


def compare(case, r):
    """(values that differ, values compared, max |difference|) against IDL."""
    nd = nt = 0
    dmax = 0.0
    k = 0
    for o in r["orders"]:
        for a in range(1, len(r["signs"]) + 1):
            ref = fits.getdata(os.path.join(REF, f"{case}_spec_{k}.fits"))
            p = r["spectra"][(o, a)]
            if p.shape != ref.shape:
                raise AssertionError(f"{case} order {o} ap {a}: shape {p.shape} vs IDL {ref.shape}")
            ok = (p == ref) | (np.isnan(p) & np.isnan(ref))
            nd += int((~ok).sum())
            nt += p.size
            if (~ok).any():
                dmax = max(dmax, float(np.nanmax(np.abs(p - ref)[~ok])))
            k += 1
    return nd, nt, dmax


if __name__ == "__main__":
    fl = F.read_flat(os.path.join(REF, "flat.fits"))
    wc = W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"), rotation=fl["rotation"])
    atran = read_atran(fl["rp"])
    with tempfile.TemporaryDirectory() as tmp:
        for case in sys.argv[1:] or CASES:
            r = run(case, fl, wc, atran, tmp)
            nd, nt, dmax = compare(case, r)
            print(f"  {case}: {nd} of {nt} values differ from IDL (max |d| {dmax:.3g}); "
                  f"signs {[int(v) for v in r['signs']]}")
            if case == "it":
                same = np.array_equal(r["positions"], fits.getdata(os.path.join(REF, "it_appos.fits")),
                                      equal_nan=True)
                print(f"      Import Trace positions identical to IDL: {same}")
        for case in ("op", "fx", "os"):
            same, _ = trace_text(case, fl, wc, atran, tmp)
            print(f"  trace file {case}: text identical to IDL's: {same}")
