"""
The xcombspec and xtellcor options beyond the default path, against IDL.

The IDL outputs (validation/idl_out/ct_*.fits, ref_combtell.log) come from
validation/idl/ref_combtell.pro, which runs Spextool's own XCOMBSPEC_COMBINESPEC,
_CORRECTSPEC, _WRITEFILE, XTELLCOR_WRITEFILE, XTELLCOR_BASIC_TELLCOR and
XTELLCOR_FINISH_TELLCOR. Data are required bit for bit, headers card for card
(MC_ADDHISTLABEL's banner is compared by its text: IDL writes it from column 2,
which is not a FITS record, and it is written here as blank-keyword cards). The
spectral-shape correction goes through IDL's FFT: there the corrected stacks
get a tolerance, and the combination of IDL's own corrected stacks must be exact.
"""

import os
import sys
import warnings

import numpy as np
import pytest
from astropy.io import fits

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
REF = os.path.join(ROOT, "validation", "idl_out")
pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(REF, "ct_fin_0122.fits")),
                                reason="IDL reference outputs (ref_combtell.pro) not present")

from triplespec_reducer import combine as C, telluric as T            # noqa: E402

warnings.filterwarnings("ignore", module="astropy")


def f(n):
    return os.path.join(REF, n)


def same(a, b, what=""):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape, f"{what}: shape {a.shape} vs IDL {b.shape}"
    ok = (a == b) | (np.isnan(a.astype(float)) & np.isnan(b.astype(float)))
    assert ok.all(), f"{what}: {int((~ok).sum())} values differ from IDL"


def records(path):
    with open(path, "rb") as fh:
        b = fh.read(2880 * 40)
    out = []
    for i in range(0, len(b), 80):
        r = b[i:i + 80].decode("ascii", "replace")
        if r.startswith("END     "):
            break
        out.append(r)
    return out


def split_header(path):
    """The cards (structure cards left out) and the banners' texts, in order."""
    cards, labels, prev = [], [], False
    for r in records(path):
        if r.startswith(("SIMPLE", "BITPIX", "NAXIS", "EXTEND")):
            continue
        if "#" in r[:8] or (r[:8].strip() == "" and "#" in r):
            text = r.strip()
            if prev:
                labels[-1] += text
            else:
                labels.append(text)
            prev = True
            continue
        prev = False
        if r.strip():
            cards.append(r.rstrip())
    return cards, labels


def same_file(ours, idl, data=True):
    if data:
        same(fits.getdata(ours), fits.getdata(idl), os.path.basename(idl))
    a, b = split_header(ours), split_header(idl)
    assert a[0] == b[0], f"{os.path.basename(idl)}: header cards differ"
    assert a[1] == b[1], f"{os.path.basename(idl)}: history labels differ"


def logged(tag):
    for line in open(f("ref_combtell.log")):
        if line.startswith(tag + " "):
            return np.array([float(v) for v in line.split()[1:]], np.float32)
    raise KeyError(tag)


FLAGS = [f"ct_flag{i}.fits" for i in range(4)]
TWO = [f"ct_2ap{i}.fits" for i in range(3)]


def prune_xc1(L):
    C.prune(L, 5, remove=[1], regions=[(1.20, 1.205, [0, 2]), (1.25, 1.252, [0, 2, 3])])
    C.prune(L, 3, remove=[3])
    C.prune(L, 4, remove=[0], regions=[(1.60, 1.62, [1, 2, 3])])


def test_pruned_statistics_and_masks():
    """Removed spectra and masked ranges: all nine statistics, IDL's full-stack mask."""
    L = C.load([f(n) for n in FLAGS])
    prune_xc1(L)
    for s in range(9):
        same(C.combine(L, statistic=C.STATISTICS[s]), fits.getdata(f(f"ct_xc1_s{s}.fits")),
             C.STATISTICS[s])
    C.combine(L)
    C.scale(L, 5)
    same(C.combine(L), fits.getdata(f("ct_xc1_scaled.fits")), "scaled after pruning")
    with pytest.raises(ValueError):
        C.scale(L, 5)                                 # XCOMBSPEC refuses a second scaling
    # MC_SPECCOR's robust mean with the full-stack mask
    for j in range(5):
        st = L["spc"][0][j][:, 1]
        z = np.where(L["spcmask"][0][j] == 1)[0]
        m, _, g = C.meancomb_stackmask(st[z], L["pixmask"][0][j].astype(np.float32), robust=3)
        r = fits.getdata(f(f"ct_xc1_sc_{j}.fits"))
        same(m, r[0], f"speccor mean {j}")
        same(g, r[1:], f"speccor mask {j}")


@pytest.mark.parametrize("name", ["ct_xc1", "ct_xc5"])
def test_corrected_after_pruning(name, tmp_path):
    if name == "ct_xc1":
        L = C.load([f(n) for n in FLAGS])
        prune_xc1(L)
        stat = 0
    else:
        L = C.load([f(n) for n in TWO], combine_apertures=True)
        C.prune(L, 5, remove=[5])
        stat = 4
    C.scale(L, 5)
    assert C.correct_shapes(L)
    with pytest.raises(ValueError):
        C.correct_shapes(L)
    for j in range(5):
        ref = fits.getdata(f(f"{name}_cor_{j}.fits")).astype(np.float32)
        r = np.abs(L["spc"][0][j][:, 1:3] - ref[:, 1:3]) / np.abs(ref[:, 1:3])
        assert np.nanmax(r) < 1e-5                    # IDL's FFT
        L["spc"][0][j] = ref
    C.combine(L, statistic=C.STATISTICS[stat])
    p = str(tmp_path / (name + ".fits"))
    C.write(p, L)
    same_file(p, f(name + ".fits"))


def test_scale_modes_basic_header(tmp_path):
    L = C.load([f(n) for n in FLAGS])
    C.scale(L, 5, mode="Spectrum", spectrum=2, srange=(1.21, 1.31))
    C.combine(L, statistic=C.STATISTICS[5])
    p = str(tmp_path / "ct_xc2.fits")
    C.write(p, L, basic=True)
    same_file(p, f("ct_xc2.fits"))
    L = C.load([f(n) for n in FLAGS])
    C.scale(L, 4, mode="None")
    C.combine(L, statistic=C.STATISTICS[6])
    p = str(tmp_path / "ct_xc3.fits")
    C.write(p, L)
    same_file(p, f("ct_xc3.fits"))
    L = C.load([f(n) for n in FLAGS])
    C.correct_shapes(L)
    with pytest.raises(ValueError):
        C.scale(L, 4)                                 # no scaling after the correction


def test_two_apertures(tmp_path):
    """One scale order per aperture; MC_SPECCOR refuses XCOMBSPEC's 3-D pixel mask."""
    L = C.load([f(n) for n in TWO])
    C.scale(L, 5, aperture=1)
    C.scale(L, 4, aperture=2)
    before = [s_.copy() for s_ in L["spc"][0]]
    assert C.correct_shapes(L, aperture=1) is False
    for a, b in zip(before, L["spc"][0]):
        same(a, b, "uncorrected")
    C.prune(L, 6, aperture=2, remove=[0], regions=[(1.05, 1.06, [1, 2])])
    C.combine(L)
    p = str(tmp_path / "ct_xc4.fits")
    C.write(p, L)
    same_file(p, f("ct_xc4.fits"))


def test_intended_masks_option():
    """The opt-in: only the kept spectra's mask rows (not what IDL computes)."""
    data = np.array([[1.], [2.], [3.]], np.float32)
    mask = np.ones((4, 1), int)
    assert C.meancomb_stackmask(data, mask)[0][0] == np.float32(1.5)       # IDL: 6/4
    assert C.medcomb_stackmask(data, mask)[0][0] == np.float32(2.5)        # IDL: [1,2,3,3]
    L = C.load([f(n) for n in FLAGS])
    C.prune(L, 5, remove=[1])
    a = C.combine(L, statistic="Mean (RMS)", intended_masks=True)[2, 1]
    st = L["spc"][0][2][[0, 2, 3], 1]
    np.testing.assert_allclose(a, st.astype(float).mean(axis=0), rtol=1e-6)


def tellcor_state():
    log = {}
    for line in open(f("ref_dg.log")):
        t = line.split()
        if t and t[0].isupper() and t[0].isalpha():
            log.setdefault(t[0], t[1:])
    st = T.load(f("comb_HD225200.fits"), f("comb_NGC134.fits"), 6.33, 6.31)
    st["vshift"] = np.float32(float(log["VSHIFT"][0]))
    st["scale"] = np.float32(float(log["KERNEL"][1]))
    st["maxdev"] = np.float32(float(log["KERNEL"][2]))
    st["rmsdev"] = np.float32(float(log["KERNEL"][3]))
    st["method"] = "Deconvolution"
    n_ = st["std"]["spec"].shape[2]
    tell = np.zeros((5, 4, n_), np.float32)
    vega = np.zeros_like(tell)
    for n, o in enumerate(st["objorders"]):
        s_ = st["std"]["spec"][st["stdorders"].index(o)]
        k = fits.getdata(f(f"tc_kern_{o}.fits")).astype(np.float32)
        tc, te, sv = T.mktellspec(s_[0], s_[1], s_[2], st["vmag"], st["bmag"] - st["vmag"], k,
                                  np.full(n_, st["scale"], np.float32), st["wvega"], st["fvega"],
                                  st["cfvega"], st["cf2vega"], st["vshift"])
        tell[n], vega[n] = s_, s_
        tell[n, 1], tell[n, 2] = tc, te
        vega[n, 1], vega[n, 2] = sv, 1.0
    st["tellspec"], st["vegaspec"], st["units"] = tell, vega, "ergs s-1 cm-2 A-1"
    return st, [float(v) for v in log["SHIFTS"]]


def test_xtellcor_three_files(tmp_path):
    """XTELLCOR_WRITEFILE: corrected, telluric and Model (modVega) files."""
    st, sh = tellcor_state()
    T.corrected(st, shifts={(o, 1): v for o, v in zip(st["objorders"], sh)})
    p, tp, vp = T.output_names(str(tmp_path / "ct_tc_NGC134.fits"))
    T.write(p, st, tellpath=tp, model=True)
    for ours, idl in ((p, "ct_tc_NGC134.fits"), (tp, "ct_tc_NGC134_tellspec.fits"),
                      (vp, "ct_tc_NGC134_modVega.fits")):
        same_file(ours, f(idl))


SR = [(2.00, 2.06), (1.43, 1.48), (1.11, 1.15), (0.93, 0.97), (0.94, 1.00)]


def autoshift(fn, st, order, sr):
    try:
        fn(st, order, 1, srange=sr)
    except ValueError:                               # IDL stops at the edge; 0 is typed
        T.set_shift(st, order, 1, 0.0)


def test_xtellcor_basic(tmp_path):
    b = T.load_basic(f("comb_HD225200.fits"), f("comb_NGC134.fits"))
    for o in b["objorders"]:
        autoshift(T.getshift_basic, b, o, None)
    same(b["shift"][:, 0], logged("BSHIFTS1"), "shifts, whole orders")
    T.tellcor_basic(b, restore=True, temp=9480, vmag=6.31)
    T.write_basic(str(tmp_path / "ct_bas1.fits"), b)
    same_file(str(tmp_path / "ct_bas1.fits"), f("ct_bas1.fits"))
    b = T.load_basic(f("ct_std_no5.fits"), f("comb_NGC134.fits"))
    autoshift(T.getshift_basic, b, 4, SR[1])
    same(b["shift"][:, 0], logged("BSHIFTS2"), "shift, chosen range")
    T.tellcor_basic(b, restore=False, units="Jy")
    T.write_basic(str(tmp_path / "ct_bas2.fits"), b)
    same_file(str(tmp_path / "ct_bas2.fits"), f("ct_bas2.fits"))
    b = T.load_basic(f("comb_HD225200.fits"), f("spectra0122.fits"))
    for n, o in enumerate(b["objorders"]):
        autoshift(T.getshift_basic, b, o, SR[n])
    same(b["shift"][:, 0], logged("BSHIFTS3"), "shifts, ranges")
    T.tellcor_basic(b, restore=True, temp=5800., vmag=7.5, units="W m-2 Hz-1")
    T.write_basic(str(tmp_path / "ct_bas3.fits"), b)
    same_file(str(tmp_path / "ct_bas3.fits"), f("ct_bas3.fits"))


def test_xtellcor_finish(tmp_path):
    fs = T.load_finish(f("ct_tc_NGC134_tellspec.fits"), f("spectra0122.fits"))
    for n, o in enumerate(fs["objorders"]):
        autoshift(T.getshift_finish, fs, o, SR[1] if n == 1 else None)
    same(fs["shift"][:, 0], logged("FSHIFTS"), "shifts")
    T.tellcor_finish(fs)
    p = str(tmp_path / "ct_fin_0122.fits")
    T.write_finish(p, fs)
    same_file(p, f("ct_fin_0122.fits"))
