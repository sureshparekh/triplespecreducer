"""
Stages D-G (xcombspec, xtellcor, xmergeorders, xcleanspec), each compared
with IDL (validation/idl/ref_dg.pro) on IDL's own input to that stage, so a
difference belongs to the stage that shows it. The extraction of the four
pairs is compared too, as written files.
"""

import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import (flat as F, wavecal as W, extract as E, combine as C,   # noqa: E402
                                telluric as T, merge as M, clean as CL)
from triplespec_reducer.instrument import read_atran                                    # noqa: E402

RAW = os.path.join(os.path.dirname(HERE), "data", "raw", "2025-11-02")
REF = os.path.join(HERE, "idl_out")
OUT = os.path.join(HERE, "py_out")
os.makedirs(OUT, exist_ok=True)


def stat(name, a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.shape != b.shape:
        print(f"  {name:34s} SHAPE {a.shape} vs IDL {b.shape}")
        return
    m = np.isfinite(a) & np.isfinite(b)
    d = np.abs(a - b)[m]
    rel = d / np.maximum(np.abs(b[m]), 1e-30)
    fin = np.array_equal(np.isfinite(a), np.isfinite(b))
    print(f"  {name:34s} rel. diff median {np.median(rel):.1e}  99% {np.percentile(rel, 99):.1e}"
          f"  max {rel.max():.1e}   finite-mask {'same' if fin else 'DIFFERS (%d vs %d)' % (np.isfinite(a).sum(), np.isfinite(b).sum())}")


def idlcube(path):
    """An IDL [nwave, ncol, n] cube read by astropy is already [n, ncol, nwave]."""
    return fits.getdata(path).astype(float)


fl = F.read_flat(os.path.join(REF, "flat.fits"))
wc = W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"), rotation=fl["rotation"])
atran = read_atran(fl["rp"])

print("Extraction to files (xspextool A-B):")
pairs = [((116, 117), "point"), ((118, 119), "point"), ((120, 121), "extended"),
         ((126, 127), "extended")]
for (a, b), kind in pairs:
    files = [os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits") for n in (a, b)]
    op = [os.path.join(OUT, f"spectra{a:04d}.fits"), os.path.join(OUT, f"spectra{b:04d}.fits")]
    r = E.reduce_pair(files, fl, wc, op, source=kind, atran=atran, xs_positions=[20.1],
                      flatname="flat.fits", wavecalname="wavecal_NGC134.fits",
                      optimal=False)          # ref_dg.pro: sum with background
    for p in r["written"]:
        ref = os.path.join(REF, os.path.basename(p))
        if os.path.exists(ref):
            stat(os.path.basename(p), idlcube(p), idlcube(ref))
        else:
            print("  ", os.path.basename(p), "has no IDL counterpart")
    for n in (a, b):
        if os.path.exists(os.path.join(REF, f"spectra{n:04d}.fits")) and \
                not os.path.exists(os.path.join(OUT, f"spectra{n:04d}.fits")):
            print(f"   IDL wrote spectra{n:04d} and python did not")

# the numbers ref_dg.pro printed (validation/idl_out/ref_dg.log)
IDLLOG = {}
for _l in open(os.path.join(REF, "ref_dg.log")):
    _t = _l.split()
    if _t and _t[0].isupper() and _t[0].isalpha():
        IDLLOG.setdefault(_t[0], _t[1:])

print("xcombspec (on IDL's spectra):")
for files, oname in ((["spectra0116", "spectra0117", "spectra0118", "spectra0119"], "comb_HD225200"),
                     (["spectra0120", "spectra0126"], "comb_NGC134")):
    L = C.load([os.path.join(REF, f + ".fits") for f in files])
    sc, _ = C.scale(L, 5)
    print("   scales", np.round(sc, 6))
    comb = C.combine(L)
    ref = idlcube(os.path.join(REF, oname + ".fits"))
    stat(oname + " (w, f, e)", comb[:, :3], ref[:, :3])
    C.write(os.path.join(OUT, oname + ".fits"), L)

print("xtellcor (on IDL's combined spectra):")
st = T.load(os.path.join(REF, "comb_HD225200.fits"), os.path.join(REF, "comb_NGC134.fits"), 6.33, 6.31)
print(f"   airmass std {st['stdairmass']:.4f} obj {st['objairmass']:.4f}")
T.conkernel(st, order=6, normreg=[(0.990, 0.997), (1.013, 1.020)], normdegree=3,
            linereg=(0.997, 1.013))
print(f"   vshift {st['vshift']:.6f} (IDL {IDLLOG['VSHIFT'][0]})  scale {st['scale']:.8f} (IDL {IDLLOG['KERNEL'][1]})"
      f"  maxdev {st['maxdev']:.8f} (IDL {IDLLOG['KERNEL'][2]})  rmsdev {st['rmsdev']:.8f} (IDL {IDLLOG['KERNEL'][3]})")
for o in st["stdorders"]:
    stat(f"kernel resampled, order {o}", st["kernels"][o], fits.getdata(os.path.join(REF, f"tc_kern_{o}.fits")))
# The resampling alone, from IDL's own (FFT-made) kernel: exact.
ks = T.resample_kernel(fits.getdata(os.path.join(REF, "tc_kernel.fits")).astype(np.float32), st["stddisp"],
                       st["stdorders"].index(6), np.float32(1.49774e-5), st["stdorders"])
for o in st["stdorders"]:
    stat(f"IDL kernel resampled, order {o}", ks[o], fits.getdata(os.path.join(REF, f"tc_kern_{o}.fits")))
# Each stage below starts from IDL's own result of the stage before: the kernel
# is the one input IDL's FFT makes (see docs/validation.md), so the telluric
# spectrum is built from IDL's kernels, scale and velocity shift.
st["kernels"] = {o: fits.getdata(os.path.join(REF, f"tc_kern_{o}.fits")).astype(np.float32)
                 for o in st["stdorders"]}
st["vshift"] = np.float32(float(IDLLOG["VSHIFT"][0]))
st["scale"] = np.float32(float(IDLLOG["KERNEL"][1]))
reft = idlcube(os.path.join(REF, "tc_tellspec.fits"))
for n, o in enumerate(st["objorders"]):
    i = st["stdorders"].index(o)
    s_ = st["std"]["spec"][i]
    sc = np.full(s_.shape[1], st["scale"], np.float32)
    tc, te, _ = T.mktellspec(s_[0], s_[1], s_[2], st["vmag"], st["bmag"] - st["vmag"], st["kernels"][o],
                             sc, st["wvega"], st["fvega"], st["cfvega"], st["cf2vega"], st["vshift"])
    stat(f"telluric spectrum order {o} (f, e)", np.array([tc, te]), reft[n, 1:3])
st["tellspec"] = reft.astype(np.float32)
sr = [(2.00, 2.06), (1.43, 1.48), (1.11, 1.15), (0.93, 0.97), (0.94, 1.00)]
for n, o in enumerate(st["objorders"]):
    try:
        T.getshift(st, o, srange=sr[n])
    except ValueError:
        print(f"   order {o}: minimum at the edge (IDL too)")
        st["shift"][n, 0] = 0.0
idlsh = np.array([float(v) for v in IDLLOG["SHIFTS"]], dtype=np.float32)
print("   shifts identical:", np.array_equal(st["shift"][:, 0].astype(np.float32), idlsh),
      "\n   shifts python", st["shift"][:, 0], "\n   shifts IDL   ", idlsh)
st["shift"][:, 0] = idlsh
cor = T.corrected(st)
refc = idlcube(os.path.join(REF, "tc_NGC134.fits"))
for n, o in enumerate(st["objorders"]):
    stat(f"corrected order {o} (f, e)", cor[n, 1:3], refc[n, 1:3])

print("xmergeorders (on IDL's corrected spectrum):")
m, used, done = M.merge_orders(refc, st["objorders"], anchor=3, sequence=[4, 5, 6, 7])
refm = fits.getdata(os.path.join(REF, "m_NGC134.fits")).astype(float)
stat("merged (w, f, e)", m[:3], refm[:3])

print("xcleanspec (on IDL's merged spectrum):")
z = np.isfinite(refm[1])
cf, ce = CL.convolvespec(np.arange(refm.shape[1], dtype=np.float32), refm[1][z], 2.4, error=refm[2][z])
refg = fits.getdata(os.path.join(REF, "cl_gauss.fits")).astype(float)
stat("Gaussian FWHM 2.4 (f, e)", np.array([cf, ce]), refg)
sf, se = CL.savitzky_golay(refm[1][z], 7, degree=2, error=refm[2][z])
refs = fits.getdata(os.path.join(REF, "cl_sg.fits")).astype(float)
stat("Savitzky-Golay 7/2 (f, e)", np.array([sf, se]), refs)
