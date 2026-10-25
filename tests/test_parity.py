"""
Every step against IDL's own output for the same frames.

The IDL outputs live in validation/idl_out (made by validation/idl/*.pro
with IDL 7.1, Spextool 4.1 and TS4_SOAR); the raw frames in
data/raw/2025-11-02. Without either, the tests are skipped. Every product
is required to be bit-identical to IDL's, except the telluric kernel, which
IDL makes with its own FFT (see docs/validation.md): there the stage after
it is fed IDL's kernel, and the kernel numbers themselves get a tolerance.

    pytest tests/test_parity.py -v          (about ten minutes)
"""

import os
import sys

import numpy as np
import pytest
from astropy.io import fits

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
REF = os.path.join(ROOT, "validation", "idl_out")
RAW = os.path.join(ROOT, "data", "raw", "2025-11-02")
pytestmark = pytest.mark.skipif(not (os.path.isdir(REF) and os.path.isdir(RAW)),
                                reason="IDL reference outputs or raw frames not present")

from triplespec_reducer import (idl as I, flat as F, wavecal as W, extract as E,   # noqa: E402
                                combine as C, telluric as T, merge as M, clean as CL)
from triplespec_reducer.instrument import read_atran                               # noqa: E402


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    assert a.shape == b.shape, f"shape {a.shape} vs IDL {b.shape}"
    assert np.array_equal(np.isfinite(a), np.isfinite(b)), "NaNs in different places"
    m = np.isfinite(a) & (np.abs(b) > 0)
    r = np.abs(a - b)[m] / np.abs(b[m])
    return np.median(r), np.percentile(r, 99)


def same(a, b, what=""):
    """Bit-identical, NaNs in the same places."""
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape, f"{what}: shape {a.shape} vs IDL {b.shape}"
    ok = (a == b) | (np.isnan(a.astype(float)) & np.isnan(b.astype(float)))
    assert ok.all(), f"{what}: {int((~ok).sum())} values differ from IDL"


def idllog():
    out = {}
    for line in open(os.path.join(REF, "ref_dg.log")):
        t = line.split()
        if t and t[0].isupper() and t[0].isalpha():
            out.setdefault(t[0], t[1:])
    return out


def raw(n):
    return os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits")


@pytest.fixture(scope="module")
def flatinfo():
    return F.read_flat(os.path.join(REF, "flat.fits"))


@pytest.fixture(scope="module")
def wcinfo(flatinfo):
    return W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"), rotation=flatinfo["rotation"])


def test_flat():
    on = [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in range(2, 52)]
    off = [os.path.join(RAW, f"SPEC_DFLAT_{i:04d}.fits") for i in range(52, 104)]
    r = F.make_flat(on, off)
    same(r["combined"], fits.getdata(os.path.join(REF, "flat_mean.fits")), "scaled flat mean")
    ref = F.read_flat(os.path.join(REF, "flat.fits"))
    x = np.arange(2048.)
    for i in range(5):
        for j in range(2):
            # the reference edges come back from FITS header text, to ~1e-9 px
            d = np.abs(I.poly(x, r["edgecoeffs"][i, j]) - I.poly(x, ref["edgecoeffs"][i, j])).max()
            assert d < 1e-8, f"order edge {i},{j} off by {d} px"
    assert np.array_equal(r["xranges"], ref["xranges"])
    m = ref["omask"] > 0
    same(np.where(m, r["flat"], np.nan), np.where(m, ref["flat"], np.nan), "normalized flat")
    same(np.where(m, r["model"], np.nan),
         np.where(m, fits.getdata(os.path.join(REF, "flat_model.fits")), np.nan), "normalization model")
    # the header keeps the RMS to seven decimals
    assert np.allclose(r["rms"], ref["rms"], rtol=0, atol=5.01e-8)


def test_wavecal(flatinfo):
    files = [raw(n) for n in range(120, 132)]
    r = W.make_wavecal(files, flatinfo)
    assert r["sol1d"]["nlines"] == [31, 64, 63, 76, 32]
    for i, o in enumerate(r["orders"]):
        ref = fits.getdata(os.path.join(REF, f"wc_lines2d_{o}.fits"))
        L = r["lines2d"][i]
        same(np.vstack([L[k] for k in ("sx", "sy", "ss", "sl", "sw")]), ref, f"traced lines, order {o}")
    ref = W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"), rotation=3)
    same(r["wavecal"], ref["wavecal"], "wavelength image")
    same(r["spatcal"], ref["spatcal"], "slit-position image")


def test_extract_point_and_extended(flatinfo, wcinfo):
    aw, at = read_atran(flatinfo["rp"])
    work, var, _, _ = E.load_image([raw(116), raw(117)], flatinfo)
    same(work, fits.getdata(os.path.join(REF, "ps_work.fits")), "work image")
    prof = E.mkspatprof2d(work, flatinfo["omask"], wcinfo["wavecal"], wcinfo["spatcal"],
                          flatinfo["slith_arc"] / flatinfo["slith_pix"], 1, aw, at, 0.4)
    for o in flatinfo["orders"]:
        same(np.array(prof[o]), fits.getdata(os.path.join(REF, f"ps_prof_{o}.fits")), f"profile {o}")
    pos, _ = E.findpeaks(prof, 2, mode="auto")
    same(pos, fits.getdata(os.path.join(REF, "ps_positions.fits")), "aperture positions")
    tc, _ = E.tracespec2d(work, wcinfo["wavecal"], wcinfo["spatcal"], flatinfo["edgecoeffs"],
                          flatinfo["xranges"], pos, 5, 5, 3.0, 0.1, 1)
    same(tc, fits.getdata(os.path.join(REF, "ps_trace.fits")), "trace coefficients")
    sign = E.psdefine_signs(prof, pos)
    sp = E.extspec2d(work, var, flatinfo["omask"], flatinfo["orders"], wcinfo["wavecal"],
                     wcinfo["spatcal"], flatinfo["ps"], sign, tc, np.ones((5, 2)),
                     psbginfo=[2.5, 5.0], bgdeg=1)
    k = 0
    for o in flatinfo["orders"]:
        for a in (1, 2):
            same(sp[(o, a)], fits.getdata(os.path.join(REF, f"ps_spec_{k}.fits")),
                 f"point source order {o} ap {a}")
            k += 1
    work, var, _, _ = E.load_image([raw(120), raw(121)], flatinfo)
    sp = E.extspec2d(work, var, flatinfo["omask"], flatinfo["orders"], wcinfo["wavecal"],
                     wcinfo["spatcal"], flatinfo["ps"], np.array([1]), np.array([[20.1, 0.0]] * 5),
                     np.ones((5, 1)))
    for k, o in enumerate(flatinfo["orders"]):
        same(sp[(o, 1)], fits.getdata(os.path.join(REF, f"xs_spec_{k}.fits")), f"extended order {o}")


def test_combine():
    for files, name in (((116, 117, 118, 119), "comb_HD225200"), ((120, 126), "comb_NGC134")):
        L = C.load([os.path.join(REF, f"spectra{n:04d}.fits") for n in files])
        C.scale(L, 5)
        comb = C.combine(L)
        same(comb[:, :3], fits.getdata(os.path.join(REF, name + ".fits"))[:, :3], name)


def test_telluric_merge_clean():
    log = idllog()
    st = T.load(os.path.join(REF, "comb_HD225200.fits"), os.path.join(REF, "comb_NGC134.fits"),
                6.33, 6.31)
    T.conkernel(st, order=6, normreg=[(0.990, 0.997), (1.013, 1.020)], normdegree=3,
                linereg=(0.997, 1.013))
    assert st["vshift"] == np.float32(float(log["VSHIFT"][0]))
    assert st["scale"] == np.float32(float(log["KERNEL"][1]))
    assert np.float32(st["maxdev"]) == np.float32(float(log["KERNEL"][2]))
    # the kernel comes from IDL's own FFT: agreement to its rounding
    for o in st["stdorders"]:
        k = fits.getdata(os.path.join(REF, f"tc_kern_{o}.fits"))
        assert np.max(np.abs(st["kernels"][o] - k)) / np.max(np.abs(k)) < 1e-5
    # its resampling to each order, from IDL's kernel: exact
    ks = T.resample_kernel(fits.getdata(os.path.join(REF, "tc_kernel.fits")).astype(np.float32),
                           st["stddisp"], st["stdorders"].index(6), np.float32(1.49774e-5),
                           st["stdorders"])
    for o in st["stdorders"]:
        same(ks[o], fits.getdata(os.path.join(REF, f"tc_kern_{o}.fits")), f"resampled kernel {o}")
    st["kernels"] = {o: fits.getdata(os.path.join(REF, f"tc_kern_{o}.fits")).astype(np.float32)
                     for o in st["stdorders"]}
    reft = fits.getdata(os.path.join(REF, "tc_tellspec.fits"))
    for n, o in enumerate(st["objorders"]):
        s_ = st["std"]["spec"][st["stdorders"].index(o)]
        tc, te, _ = T.mktellspec(s_[0], s_[1], s_[2], st["vmag"], st["bmag"] - st["vmag"],
                                 st["kernels"][o], np.full(s_.shape[1], st["scale"], np.float32),
                                 st["wvega"], st["fvega"], st["cfvega"], st["cf2vega"], st["vshift"])
        same(np.array([tc, te]), reft[n, 1:3], f"telluric spectrum order {o}")
    st["tellspec"] = reft.astype(np.float32)
    sr = [(2.00, 2.06), (1.43, 1.48), (1.11, 1.15), (0.93, 0.97), (0.94, 1.00)]
    for n, o in enumerate(st["objorders"]):
        try:
            T.getshift(st, o, srange=sr[n])
        except ValueError:                    # order 6: IDL's minimum is at the edge too
            st["shift"][n, 0] = 0.0
    same(st["shift"][:, 0], np.array([float(v) for v in log["SHIFTS"]], np.float32), "shifts")
    cor = T.corrected(st)
    refc = fits.getdata(os.path.join(REF, "tc_NGC134.fits"))
    same(cor[:, 1:3], refc[:, 1:3], "corrected spectrum")
    m, _, _ = M.merge_orders(refc, st["objorders"], anchor=3, sequence=[4, 5, 6, 7])
    refm = fits.getdata(os.path.join(REF, "m_NGC134.fits"))
    same(m[:3], refm[:3], "merged spectrum")
    z = np.isfinite(refm[1])
    cf, ce = CL.convolvespec(np.arange(refm.shape[1], dtype=np.float32), refm[1][z].astype(float),
                             2.4, error=refm[2][z].astype(float))
    same(np.array([cf, ce]), fits.getdata(os.path.join(REF, "cl_gauss.fits")), "Gaussian smoothing")
    sf, se = CL.savitzky_golay(refm[1][z].astype(float), 7, degree=2, error=refm[2][z].astype(float))
    same(np.array([sf, se]), fits.getdata(os.path.join(REF, "cl_sg.fits")), "Savitzky-Golay")


def test_wavecal_arc(flatinfo, tmp_path):
    """The CuHeAr lamp solution (ref_wavecal_arc.pro), file and header."""
    files = [os.path.join(RAW, f"SPEC_ARC_{n:04d}.fits") for n in range(104, 116)]
    r = W.make_wavecal(files, flatinfo, arctype="CuHeAr Arc")
    assert r["sol1d"]["nlines"] == [14, 16, 39, 23, 5]
    for i, o in enumerate(r["orders"]):
        ref = fits.getdata(os.path.join(REF, f"wca_lines2d_{o}.fits"))
        L = r["lines2d"][i]
        same(np.vstack([L[k] for k in ("sx", "sy", "ss", "sl", "sw")]), ref, f"arc lines, order {o}")
    p = str(tmp_path / "wavecal_arc.fits")
    W.write_wavecal(p, r)
    ours, idl = fits.getheader(p), fits.getheader(os.path.join(REF, "wavecal_arc.fits"))
    for k in idl:
        if k.startswith(("1DWC", "OR0", "DISPO", "WDEG", "ODEG", "WXDEG", "SXDEG", "HOMEORDR")):
            assert ours.cards[k].image == idl.cards[k].image, k
    a = W.read_wavecal(p, rotation=3)
    b = W.read_wavecal(os.path.join(REF, "wavecal_arc.fits"), rotation=3)
    same(a["wavecal"], b["wavecal"], "arc wavelength image")
    same(a["spatcal"], b["spatcal"], "arc slit-position image")
    same(a["disp"], b["disp"], "dispersions read back")


def test_header_round_trip():
    """FXADDPAR text and FXPAR typing (probed in IDL 7.1): what the next step reads."""
    h = fits.Header()
    I.fxaddpar(h, "F", np.float32(38.0), " c")
    I.fxaddpar(h, "D", 1.39090757573154, " c")
    I.fxaddpar(h, "G", 38.0, " c", fmt="G0.5")
    I.fxaddpar(h, "S", np.float32(2.4), " c", fmt="G0.5")
    assert h.cards["F"].image.startswith("F       =              38.0000 / c")
    assert h.cards["D"].image.startswith("D       =        1.39090757573 / c")
    assert isinstance(I.fxpar(h, "F"), np.float64)          # '38.0000 ' is 8 characters
    assert I.fxpar(h, "D") == 1.39090757573
    assert isinstance(I.fxpar(h, "G"), int)                  # '38': a LONG
    assert isinstance(I.fxpar(h, "S"), np.float32)           # '2.4 ': a FLOAT


def _mcl_log():
    out = {}
    for line in open(os.path.join(REF, "ref_mergeclean.log")):
        t = line.replace("|", " ").split()
        if t and t[0].isupper():
            out.setdefault(t[0], []).append(line.rstrip("\n"))
    return out


needs_mcl = pytest.mark.skipif(not os.path.exists(os.path.join(REF, "ref_mergeclean.log")),
                               reason="ref_mergeclean.pro outputs not present")
MCL_TRIMS = {5: [("Green", "Right", 1.47), ("Blue", "Left", 1.44)],
             6: [("Green", "Right", 1.20), ("Blue", "Left", 1.16)],
             7: ("Green", "Left", 0.97)}
MCL_SCALES = {5: "auto", 3: 1.05, 6: "auto"}
MCL_REGIONS = {5: (1.42, 1.47), 6: (1.22, 1.15)}


def _native(a):
    return a.astype(a.dtype.newbyteorder("="))


@needs_mcl
def test_merge_options(tmp_path):
    """ref_mergeclean.pro: auto scale on the uncut spectra, green + blue cuts, 'Inside'."""
    log = _mcl_log()
    cs = _native(fits.getdata(os.path.join(REF, "tc_NGC134.fits")))
    kw = dict(anchor=4, scales=MCL_SCALES, regions=MCL_REGIONS, trims=MCL_TRIMS)
    m, used, done = M.merge_orders(cs, [3, 4, 5, 6, 7], sequence=[5, 3, 6], **kw)
    assert np.float32(used[5]) == np.float32(float(log["SCALE5"][0].split()[1]))
    assert np.float32(used[6]) == np.float32(float(log["SCALE6"][0].split()[1]))
    p = str(tmp_path / "mc_merged.fits")
    M.write(p, m, fits.getheader(os.path.join(REF, "comb_NGC134.fits")), used, done,
            source="tc_NGC134.fits")
    same(fits.getdata(p), fits.getdata(os.path.join(REF, "mc_merged.fits")), "merged with cuts")
    m7, _, _ = M.merge_orders(cs, [3, 4, 5, 6, 7], sequence=[5, 3, 6, 7], **kw)
    same(m7, fits.getdata(os.path.join(REF, "mc_inside.fits")).astype(float), "'Inside' merge")


def test_merge_rules():
    """The add-order list, and 'auto' without a region refused as XMERGEORDERS refuses it."""
    assert M.addable_orders([3, 4, 5, 6, 7], [4]) == [3, 5]
    assert M.addable_orders([3, 4, 5, 6, 7], [4, 3]) == [5]
    assert M.addable_orders([3, 4, 5, 6, 7], [3, 4, 5, 6, 7]) == []
    w = np.linspace(1.0, 2.0, 50)
    spec = np.zeros((3, 4, 50), np.float32)
    for k, off in enumerate((0.0, 0.8, 1.6)):
        spec[k, 0], spec[k, 1], spec[k, 2] = w + off, 1.0, 0.1
    with pytest.raises(ValueError):
        M.merge_orders(spec, [3, 4, 5], anchor=3, sequence=[5, 4])
    with pytest.raises(ValueError):
        M.merge_orders(spec, [3, 4, 5], anchor=3, scales={4: "auto"})
    m, used, _ = M.merge_orders(spec, [3, 4, 5], anchor=3, scales={4: "auto"}, regions={4: "overlap"})
    assert used[4] == np.float32(1.0)


@needs_mcl
def test_clean_session(tmp_path):
    """ref_mergeclean.pro: xcleanspec's edit buffer, 'Error' type, defaults, one smoothing."""
    log = _mcl_log()
    S = CL.Session(os.path.join(REF, "mc_merged.fits"))
    d = log["DEFAULTS"][0].split()
    assert S.defaults["sgwin"] == np.float32(d[1]) and S.defaults["fwhm"] == np.float32(d[2])
    i = S.index()
    S.fix(i, 1.27, 1.29)
    S.fix(i, 1.30, 1.285)
    S.commit()
    S.remove(i, 1.82, 1.80, spectype="Error")
    S.commit()
    S.fix(i, 1.55, 1.56)
    S.undo()
    S.fix(i, 2.10, 2.12, spectype="Error")
    same(S.espec[0], fits.getdata(os.path.join(REF, "cl_edits.fits")), "committed edits")
    S.smooth()
    p = str(tmp_path / "cl.fits")
    S.write(p)
    same(fits.getdata(p), fits.getdata(os.path.join(REF, "cl_session.fits")), "smoothed session")
    assert S.history == log["SGHIST"][0].split("|")[1]
    with pytest.raises(RuntimeError):
        S.smooth()
    G = CL.Session(os.path.join(REF, "mc_merged.fits"))
    G.smooth(method="Gaussian")
    same(G.espec[0], fits.getdata(os.path.join(REF, "cl_session_g.fits")), "Gaussian")
    assert G.history == log["GHIST"][0][6:]
    with pytest.raises(ValueError):
        CL.check_degree(3.6, 4.)
    CL.check_degree(3.6, 3.)


@needs_mcl
def test_spex2text(tmp_path):
    from triplespec_reducer import textexport as TX
    srcs = [os.path.join(REF, f + ".fits") for f in ("comb_NGC134", "spectra0116", "mc_merged")]
    srcs.append(os.path.join(ROOT, "validation", "py_out", "HD225200_ci_spec.fits"))
    for s in srcs:
        if not os.path.exists(s):
            continue
        out, = TX.spex2text(s, outdir=str(tmp_path))
        ref = os.path.join(REF, "s2t", os.path.basename(out))
        assert open(out).read() == open(ref).read(), os.path.basename(out)


@needs_mcl
def test_xgetoffset():
    from triplespec_reducer import getoffset as GO
    from triplespec_reducer.instrument import INSTRUMENT_DIR
    cal = W.read_linecal(os.path.join(INSTRUMENT_DIR, "OHXD_LineCal.fits"))
    k = W.nantrim(cal["xcorspec"][0], 2)
    wh, fh = cal["xcorspec"][0][k], cal["xcorspec"][1][k]
    s1 = _native(fits.getdata(os.path.join(REF, "wc_spec1d_4.fits")))
    clicks = [(784.3, 783.6), (1182.2, 1181.4), (920.7, 920.2)]
    for n, line in enumerate(_mcl_log()["GOFF"]):
        v = [float(t) for t in line.split()[2:]]
        off, det = GO.xgetoffset(wh, fh, s1[0], s1[1], *clicks[n])
        assert np.array_equal(np.float32(det["template"][1]), np.float32(v[0:5]))
        assert np.array_equal(np.float32(det["observed"][1]), np.float32(v[5:10]))
        assert off == np.float32(v[12])
