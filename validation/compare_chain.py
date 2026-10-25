"""
The reduction as a chain: each Python stage reads only the files the Python
stage before it wrote (flat -> wavelength calibration -> extraction ->
combination -> telluric correction -> merged orders -> cleaned spectra), and
every file is compared with IDL's file from IDL's own chain
(validation/idl/ref_flat.pro, ref_wavecal.pro, ref_chain.pro): the data bit
for bit, the header card for card and its history line for line.

Two cases (validation/idl_out/chain_<case>/, made by ref_chain.pro):

    gui   the point sources extracted as the TS4_SOAR panel does by default
          (optimal extraction, PSF radius 2.2, background) -- the Python
          defaults
    sum   the aperture sum with the background (optimal=False), what the
          earlier references (ref_dg.pro) used

The telluric kernel is the one number IDL makes with its FFT, which is not
reproducible to the last bit (docs/validation.md): ours is compared with a
tolerance, and the chain then goes on from IDL's kernels (and the scale and
deviations made with them), so that every file after it can be exact.

    python validation/compare_chain.py [gui] [sum]

The stage-by-stage comparisons feed each stage IDL's file, so they cannot see
what a file's header carries to the next stage: Spextool writes numbers as
text (FXADDPAR) and reads that text back (FXPAR), rounded and typed by it.
This script follows the numbers through that round trip.
"""

import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import (flat as F, wavecal as W, extract as E,   # noqa: E402
                                combine as C, telluric as T, merge as M, clean as CL)
from triplespec_reducer.instrument import read_atran                      # noqa: E402

RAW = os.path.join(os.path.dirname(HERE), "data", "raw", "2025-11-02")
REF = os.path.join(HERE, "idl_out")
OUT = os.path.join(HERE, "py_chain")
os.makedirs(OUT, exist_ok=True)

# Cards that cannot agree: the software version, the date a file was written,
# and the FITS structure cards astropy writes itself.
SKIP = {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "EXTEND", "COMMENT",
        "HISTORY", "", "DATE", "VERSION",
        # MC_ADDHISTLABEL's banner starts in column 2, not a valid FITS record;
        # written here as a blank-keyword card with the same text
        "#######"}
# Ours, not IDL's: the wavecal's record of the arc type, offset and 1D rms.
OURS = {"TSRTYPE", "TSROFFS", "TSRRMS1D"}
failures = []


def rd(path, ext=0):
    with fits.open(path) as h:
        d = h[ext].data
        return d.astype(d.dtype.newbyteorder("="))


def same_data(name, a, b):
    ok = a.shape == b.shape and a.dtype == b.dtype and \
        bool(np.all((a == b) | (np.isnan(a.astype(float)) & np.isnan(b.astype(float)))))
    print(f"  {name:40s} data {'identical' if ok else 'DIFFERENT'}"
          + ("" if ok else f"  {a.shape}{a.dtype} vs IDL {b.shape}{b.dtype}"))
    if not ok:
        failures.append(name + " data")


def history(path):
    """
    The HISTORY lines and MC_ADDHISTLABEL banners, in order, from the raw
    80-character records (IDL writes a banner from column 2, which astropy
    does not read as a card; ours is a blank-keyword card): banners by text.
    """
    with open(path, "rb") as fh:
        b = fh.read(2880 * 40)
    out, prev = [], False
    for i in range(0, len(b), 80):
        r = b[i:i + 80].decode("ascii", "replace")
        if r.startswith("END     "):
            break
        if r.startswith("HISTORY "):
            out.append(("H", r[8:].rstrip()))
        elif "####" in r:
            if out and out[-1][0] == "B" and prev:         # a banner wrapped onto a 2nd record
                out[-1] = ("B", out[-1][1] + r.strip())
            else:
                out.append(("B", r.strip()))
            prev = True
            continue
        prev = False
    return out


# Ours, not IDL's: wavecal.py records the frames it used as HISTORY lines.
OURS_HISTORY = {"wavecal_NGC134.fits"}


def same_header(name, ours, idl):
    a, b = fits.getheader(ours), fits.getheader(idl)
    ha, hb = history(ours), history(idl)
    if ha != hb and name not in OURS_HISTORY:
        print(f"  {name:40s} history DIFFERENT ({len(ha)} vs IDL {len(hb)} lines)")
        failures.append(name + " history")
    diff = [k for k in b if k not in SKIP and (k not in a or a.cards[k].image != b.cards[k].image)]
    extra = [k for k in a if k not in b and k not in SKIP | OURS]
    print(f"  {name:40s} header {'identical' if not diff and not extra else 'DIFFERENT'}")
    for k in diff[:8]:
        print(f"      ours {a.cards[k].image.rstrip() if k in a else '(missing)'}\n"
              f"      IDL  {b.cards[k].image.rstrip()}")
    if extra:
        print("      only ours:", extra)
    if diff or extra:
        failures.append(name + " header")


SR = [(2.00, 2.06), (1.43, 1.48), (1.11, 1.15), (0.93, 0.97), (0.94, 1.00)]
CASES = {"gui": dict(), "sum": dict(optimal=False)}


def chain(case, fl, wc, outroot=OUT):
    """
    One case from extraction to the cleaned spectra, every stage reading the
    files of the stage before; each file checked against chain_<case>.
    Returns the merged spectrum; differences go to `failures`.
    """
    atran = read_atran(fl["rp"])
    ref = os.path.join(REF, f"chain_{case}")
    out = os.path.join(outroot, case)
    os.makedirs(out, exist_ok=True)
    log = {}
    for line in open(os.path.join(ref, "ref_chain.log")):
        t = line.split()
        if t and t[0].isupper() and t[0].isalpha():
            log.setdefault(t[0], t[1:])

    def check(name):
        same_data(name, rd(os.path.join(out, name)), rd(os.path.join(ref, name)))
        same_header(name, os.path.join(out, name), os.path.join(ref, name))

    print(f"\n=== case '{case}' (validation/idl_out/chain_{case}) ===")
    print("extraction, from our flat and wavecal files:")
    for (a, b), kind in (((116, 117), "point"), ((118, 119), "point"),
                         ((120, 121), "extended"), ((126, 127), "extended")):
        files = [os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits") for n in (a, b)]
        op = [os.path.join(out, f"spectra{a:04d}.fits"), os.path.join(out, f"spectra{b:04d}.fits")]
        r = E.reduce_pair(files, fl, wc, op, source=kind, atran=atran, xs_positions=[20.1],
                          flatname="flat.fits", wavecalname="wavecal_NGC134.fits", **CASES[case])
        for p in r["written"]:
            check(os.path.basename(p))

    print("combination, from our spectra files:")
    for names, oname in ((["spectra0116", "spectra0117", "spectra0118", "spectra0119"], "comb_HD225200"),
                         (["spectra0120", "spectra0126"], "comb_NGC134")):
        L = C.load([os.path.join(out, f + ".fits") for f in names])
        C.scale(L, 5)
        C.combine(L)
        p = os.path.join(out, oname + ".fits")
        C.write(p, L)
        check(oname + ".fits")

    print("telluric correction, from our combined files:")
    st = T.load(os.path.join(out, "comb_HD225200.fits"), os.path.join(out, "comb_NGC134.fits"), 6.33, 6.31)
    T.conkernel(st, order=6, normreg=[(0.990, 0.997), (1.013, 1.020)], normdegree=3,
                linereg=(0.997, 1.013))
    f32 = np.float32
    ok = st["vshift"] == f32(float(log["VSHIFT"][0]))
    print(f"  {'Vega velocity shift':40s} {'identical' if ok else 'DIFFERENT'}")
    if not ok:
        failures.append("vshift")
    kd = max(float(np.max(np.abs(st["kernels"][o] - rd(os.path.join(ref, f"tc_kern_{o}.fits"))))
                   / np.max(np.abs(rd(os.path.join(ref, f"tc_kern_{o}.fits"))))) for o in st["stdorders"])
    print(f"  {'kernel (IDL FFT)':40s} max |d|/max {kd:.1e}"
          f"  scale {st['scale']:.8f} vs IDL {float(log['KERNEL'][1]):.8f}")
    if kd > 1e-5:
        failures.append("kernel")
    st["kernels"] = {o: rd(os.path.join(ref, f"tc_kern_{o}.fits")).astype(f32) for o in st["stdorders"]}
    st["scale"], st["maxdev"], st["rmsdev"] = (f32(float(v)) for v in log["KERNEL"][1:4])
    st["method"], st["units"] = "Deconvolution", "ergs s-1 cm-2 A-1"
    n_ = st["std"]["spec"].shape[2]
    tell = np.zeros((len(st["objorders"]), 4, n_), f32)
    vega = np.zeros_like(tell)
    for n, o in enumerate(st["objorders"]):
        s_ = st["std"]["spec"][st["stdorders"].index(o)]
        tc, te, sv = T.mktellspec(s_[0], s_[1], s_[2], st["vmag"], st["bmag"] - st["vmag"],
                                  st["kernels"][o], np.full(n_, st["scale"], f32), st["wvega"],
                                  st["fvega"], st["cfvega"], st["cf2vega"], st["vshift"])
        tell[n], vega[n] = s_, s_
        tell[n, 1], tell[n, 2] = tc, te
        vega[n, 1], vega[n, 2] = sv, 1.0
    st["tellspec"], st["vegaspec"] = tell, vega
    for n, o in enumerate(st["objorders"]):
        try:
            T.getshift(st, o, srange=SR[n])
        except ValueError:                    # minimum at the edge: IDL types 0
            st["shift"][n, 0] = 0.0
    ok = np.array_equal(st["shift"][:, 0].astype(f32), np.array([float(v) for v in log["SHIFTS"]], f32))
    print(f"  {'shifts':40s} {'identical' if ok else 'DIFFERENT'}")
    if not ok:
        failures.append("shifts")
    T.corrected(st)
    p, tp, vp = T.output_names(os.path.join(out, "tc_NGC134.fits"))
    T.write(p, st, tellpath=tp, model=True)
    for q in (p, tp, vp):
        check(os.path.basename(q))

    print("merged orders, from our corrected file:")
    cs = E.caaspecflags(rd(p))                # XMERGEORDERS_LOADSPEC
    m, used, done = M.merge_orders(cs, st["objorders"], anchor=3, sequence=[4, 5, 6, 7])
    pm = os.path.join(out, "merged_NGC134.fits")
    M.write(pm, m, fits.getheader(p), used, done, source="tc_NGC134.fits")
    check("merged_NGC134.fits")

    print("cleaned spectra, from our merged file (the panel's default fields):")
    for method, name in (("Savitzky-Golay", "clean_sg_NGC134.fits"), ("Gaussian", "clean_g_NGC134.fits")):
        S = CL.Session(pm)
        S.smooth(method=method)
        S.write(os.path.join(out, name))
        check(name)
    return m


def calibrations():
    """Flat and wavelength calibration from the raw frames, checked against IDL's."""
    print("flat (ref_flat.pro):")
    on = [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in range(2, 52)]
    off = [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in range(52, 104)]
    F.write_flat(os.path.join(OUT, "flat.fits"), F.make_flat(on, off), version="v4.1")
    for e, what in ((1, "flat"), (2, "variance"), (3, "flags")):
        same_data(f"flat.fits [{what}]", rd(os.path.join(OUT, "flat.fits"), e), rd(os.path.join(REF, "flat.fits"), e))
    same_header("flat.fits", os.path.join(OUT, "flat.fits"), os.path.join(REF, "flat.fits"))
    fl = F.read_flat(os.path.join(OUT, "flat.fits"))

    print("wavelength calibration on sky (ref_wavecal.pro), from our flat file:")
    sky = [os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits") for n in range(120, 132)]
    W.write_wavecal(os.path.join(OUT, "wavecal_NGC134.fits"), W.make_wavecal(sky, fl))
    for e, what in ((1, "wavelengths"), (2, "slit positions")):
        same_data(f"wavecal_NGC134.fits [{what}]", rd(os.path.join(OUT, "wavecal_NGC134.fits"), e),
                  rd(os.path.join(REF, "wavecal_NGC134.fits"), e))
    same_header("wavecal_NGC134.fits", os.path.join(OUT, "wavecal_NGC134.fits"),
                os.path.join(REF, "wavecal_NGC134.fits"))
    wc = W.read_wavecal(os.path.join(OUT, "wavecal_NGC134.fits"), rotation=fl["rotation"])

    return fl, wc


if __name__ == "__main__":
    fl, wc = calibrations()
    merged = {case: chain(case, fl, wc) for case in (sys.argv[1:] or ["gui", "sum"])}
    if len(merged) == 2:
        g, s_ = merged["gui"], merged["sum"]
        z = np.isfinite(g[1]) & np.isfinite(s_[1]) & (s_[1] != 0)
        r = g[1][z] / s_[1][z]
        print(f"\nmerged NGC134, GUI-default standard / sum standard: median {np.median(r):.4f}, "
              f"16-84% {np.percentile(r, 16):.4f}-{np.percentile(r, 84):.4f}, "
              f"range {r.min():.4f}-{r.max():.4f}")
    print("\nall identical" if not failures else f"\n{len(failures)} differences: {failures}")
