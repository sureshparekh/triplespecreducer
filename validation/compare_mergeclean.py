"""
xmergeorders, xcleanspec, spex2text and xgetoffset options, compared with
IDL (validation/idl/ref_mergeclean.pro, which runs xmergeorders' own data
procedures on a state structure and makes the widget events by hand):

    merge   auto scale on the untrimmed spectra (the region reaching into
            both cuts), a green and a blue cut on the same order, a typed
            scale, a gap, the 'Inside' case; the file and its header
    clean   defaults from SLTW_PIX, two uncommitted overlapping fixes, an
            'Error' removal, an undone fix, an uncommitted fix lost to the
            smoothing, the smoothing history; MC_CRANGE
    text    spex2text on four files (FLOAT/DOUBLE, 3/4 columns, 1/5/10 spectra)
    offset  xgetoffset's line fits on three sky lines of the home order

Every comparison is exact (bit for bit, text character for character).
"""

import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from triplespec_reducer import merge as M, clean as CL, textexport as TX, getoffset as G   # noqa: E402
from triplespec_reducer import wavecal as W                                               # noqa: E402
from triplespec_reducer.instrument import INSTRUMENT_DIR                                  # noqa: E402

REF = os.path.join(HERE, "idl_out")
OUT = os.path.join(HERE, "py_out")
os.makedirs(os.path.join(OUT, "s2t"), exist_ok=True)


def log():
    out = {}
    for line in open(os.path.join(REF, "ref_mergeclean.log")):
        t = line.replace("|", " ").split()
        if t and t[0].isupper():
            out.setdefault(t[0], []).append(line.rstrip("\n"))
    return out


def same(name, a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        print(f"  {name:44s} SHAPE {a.shape} vs IDL {b.shape}")
        return False
    ok = (a == b) | (np.isnan(a.astype(float)) & np.isnan(b.astype(float)))
    print(f"  {name:44s} {'identical' if ok.all() else f'{int((~ok).sum())} values DIFFER'}")
    return bool(ok.all())


STRUCT = ("SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND")


def header_same(name, p_ours, p_idl):
    """Card for card, apart from the structural cards and MC_ADDHISTLABEL's label
    line, which IDL writes from column 2 (not a valid FITS card; idl.addhistlabel
    writes a blank-keyword card instead, see its docstring)."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = [c.image for c in fits.getheader(p_ours).cards]
        b = [c.image for c in fits.getheader(p_idl).cards]
    keep = lambda L: [c for c in L if c[:8].strip() not in STRUCT and "#####" not in c and c.strip()]
    ok = keep(a) == keep(b)
    print(f"  {name:44s}", "identical (label line apart)" if ok else "DIFFER")
    return ok


def native(a):
    return a.astype(a.dtype.newbyteorder("="))


L = log()
cs = native(fits.getdata(os.path.join(REF, "tc_NGC134.fits")))
orders = [3, 4, 5, 6, 7]
trims = {5: [("Green", "Right", 1.47), ("Blue", "Left", 1.44)],
         6: [("Green", "Right", 1.20), ("Blue", "Left", 1.16)],
         7: ("Green", "Left", 0.97)}
scales = {5: "auto", 3: 1.05, 6: "auto"}
regions = {5: (1.42, 1.47), 6: (1.22, 1.15)}

# ---- merge: orders 4, 5, 3, 6 written as XMERGEORDERS_WRITEFILE writes them
m, used, done = M.merge_orders(cs, orders, anchor=4, sequence=[5, 3, 6], scales=scales,
                               regions=regions, trims=trims)
for o, key in ((5, "SCALE5"), (6, "SCALE6")):
    idl = np.float32(float(L[key][0].split()[1]))
    print(f"  auto scale order {o}: {used[o]!r} (IDL {float(idl)!r})",
          "identical" if np.float32(used[o]) == idl else "DIFFERS")
hdr = fits.getheader(os.path.join(REF, "comb_NGC134.fits"))
p = os.path.join(OUT, "mc_merged.fits")
M.write(p, m, hdr, used, done, aperture=1, source="tc_NGC134.fits")
ref = fits.getdata(os.path.join(REF, "mc_merged.fits"))
same("merged (cuts, auto scales, gap)", fits.getdata(p), ref)
header_same("merged header cards", p, os.path.join(REF, "mc_merged.fits"))

# without the IDL rule (scale on the cut spectra) order 6 would differ
blue = M.merge_orders(cs, orders, anchor=4, sequence=[5, 3], scales=scales, regions=regions, trims=trims)[0]
add = cs[3][:, np.isfinite(cs[3, 1])]
gcut = add[:, add[0] < np.float32(1.20)]
bcut = blue[:, (blue[0] > np.float32(1.16)) | ~np.isfinite(blue[0])]
print("   (scale from the cut spectra would be", float(M.autoscale(bcut, gcut, regions[6])), ")")

# ---- the 'Inside' case: order 7 cut to lie inside the merged range
m7, _, _ = M.merge_orders(cs, orders, anchor=4, sequence=[5, 3, 6, 7], scales=scales,
                          regions=regions, trims=trims)
same("merged, order 7 'Inside'", m7, fits.getdata(os.path.join(REF, "mc_inside.fits")).astype(float))
b = native(fits.getdata(os.path.join(REF, "mc_inside_blue.fits")))
g = native(fits.getdata(os.path.join(REF, "mc_inside_green.fits")))
r = M.mergespec(b[0], b[1], g[0], g[1], b[2], g[2], b[3], g[3])
same("  MC_MERGESPEC 'Inside' alone", np.array(r[:4], float),
     fits.getdata(os.path.join(REF, "mc_inside.fits")).astype(float))

# ---- clean: the session on IDL's merged file
S = CL.Session(os.path.join(REF, "mc_merged.fits"))
d = S.defaults
dl = L["DEFAULTS"][0].split()
print(f"  defaults sgwin {float(d['sgwin']):g} fwhm {float(d['fwhm']):#.6g} (IDL {dl[1]} {dl[2]})",
      "identical" if (np.float32(dl[1]) == d["sgwin"] and np.float32(dl[2]) == d["fwhm"]) else "DIFFER")
i = S.index()
S.fix(i, 1.27, 1.29)
S.fix(i, 1.30, 1.285)
S.commit()
S.remove(i, 1.82, 1.80, spectype="Error")
S.commit()
S.fix(i, 1.55, 1.56)
S.undo()
S.fix(i, 2.10, 2.12, spectype="Error")
same("committed edits (espec)", S.espec[0], fits.getdata(os.path.join(REF, "cl_edits.fits")))
S.smooth()
pc = os.path.join(OUT, "cl_session.fits")
S.write(pc)
same("smoothed and written", fits.getdata(pc), fits.getdata(os.path.join(REF, "cl_session.fits")))
hist_idl = L["SGHIST"][0].split("|")[1]
print(f"  {'smoothing history':44s}", "identical" if S.history == hist_idl else f"DIFFERS: {S.history!r}")
header_same("cleaned header cards", pc, os.path.join(REF, "cl_session.fits"))
S2 = CL.Session(os.path.join(REF, "mc_merged.fits"))
S2.smooth(method="Gaussian")
same("Gaussian on a fresh load", S2.espec[0], fits.getdata(os.path.join(REF, "cl_session_g.fits")))
print(f"  {'Gaussian history':44s}", "identical" if S2.history == L["GHIST"][0][6:] else f"DIFFERS {S2.history!r}")
try:
    S2.smooth()
    print("  second smoothing: allowed (WRONG)")
except RuntimeError:
    print("  second smoothing: refused, as XCLEANSPEC greys it out")
crange = L["CRANGE"][0].split()[1:]
res = []
for deg in (4., 3.):
    try:
        CL.check_degree(3.6, deg)
        res.append("0")
    except ValueError:
        res.append("1")
print(f"  MC_CRANGE cancel for degree 4, 3: {res} (IDL {crange})", "identical" if res == crange else "DIFFER")

# ---- spex2text
for f in ("comb_NGC134", "spectra0116", "mc_merged"):
    TX.spex2text(os.path.join(REF, f + ".fits"), outdir=os.path.join(OUT, "s2t"))
TX.spex2text(os.path.join(OUT, "HD225200_ci_spec.fits"), outdir=os.path.join(OUT, "s2t"))
for f in ("comb_NGC134", "spectra0116", "mc_merged", "HD225200_ci_spec"):
    a = open(os.path.join(OUT, "s2t", f + ".txt")).read()
    b_ = open(os.path.join(REF, "s2t", f + ".txt")).read()
    print(f"  {'spex2text ' + f:44s}", "identical" if a == b_ else "DIFFERS")

# ---- xgetoffset on three sky lines (template and observed home order 4)
cal = W.read_linecal(os.path.join(INSTRUMENT_DIR, "OHXD_LineCal.fits"))
k = W.nantrim(cal["xcorspec"][0], 2)
wh, fh = cal["xcorspec"][0][k], cal["xcorspec"][1][k]
s1 = native(fits.getdata(os.path.join(REF, "wc_spec1d_4.fits")))
clicks = [(784.3, 783.6), (1182.2, 1181.4), (920.7, 920.2)]
for n, line in enumerate(L["GOFF"]):
    v = [float(t) for t in line.split()[2:]]
    off, det = G.xgetoffset(wh, fh, s1[0], s1[1], *clicks[n])
    f1, f2 = det["template"][1], det["observed"][1]
    ok = (np.array_equal(np.float32(f1), np.float32(v[0:5])) and np.array_equal(np.float32(f2), np.float32(v[5:10]))
          and off == np.float32(v[12]))
    print(f"  xgetoffset line {n}: offset {float(off)!r} (IDL {v[12]!r})", "identical" if ok else "DIFFERS")
