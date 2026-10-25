"""
The TS4_SOAR instrument definition, read from the files the IDL uses.

XD.dat is the mode file (orders, plate scale, slit, order-finding and
wavelength-calibration parameters) and ts4_soar.dat the instrument file
(reduction defaults). Both are parsed as MC_READMODEFILE and
MC_READINSTRFILE parse them, so changing a value there changes it here.
"""

import os
import re
from dataclasses import dataclass, field

import numpy as np

from triplespec_reducer._paths import INSTRUMENT_DIR, SPEXTOOL_DIR


def _pairs(path):
    """KEY=VALUE lines, in order, skipping '#' comments and blank lines."""
    out = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("%"):
                continue
            if "=" not in s:
                continue
            k, v = s.split("=", 1)
            out.append((k.strip(), v.strip()))
    return out


def fsextract_index(text):
    """
    MC_FSEXTRACT(/INDEX): '3-7' or '1,4-6' to a list of integers. A
    descending range runs downwards ('7-3' is 7, 6, 5, 4, 3); '5-' is 5.
    """
    out = []
    for part in str(text).split(","):
        if not part.strip():
            continue
        try:
            d = [int(t) for t in part.split("-") if t.strip()]
        except ValueError:
            raise ValueError("Type conversion error.  Probably an INDEX versus FILENAME "
                             "File Read Mode problem.") from None
        if len(d) == 1:
            d = [d[0], d[0]]
        lo, hi = min(d[0], d[1]), max(d[0], d[1])
        rng = list(range(lo, hi + 1))
        out.extend(rng[::-1] if d[1] < d[0] else rng)
    return out


def fsextract_filename(text):
    """MC_FSEXTRACT(/FILENAME): comma-separated names, all blanks removed."""
    return ["".join(t.split()) for t in str(text).split(",")]


@dataclass
class ModeInfo:
    """MC_READMODEFILE's structure, one field per entry in XD.dat."""
    modename: str
    rotation: int
    orders: list
    ps: float
    slith_pix: int
    slith_pix_range: list
    slith_arc: float
    rpppix: float
    step: int
    flatfrac: float
    comwin: int
    edgedeg: int
    norm_nxg: int
    norm_nyg: int
    oversamp: float
    norm_ybuffer: int
    fixed: str
    guesspos: np.ndarray        # [norders, 2]  (x, y)
    findrange: np.ndarray       # [norders, 2]  (start, stop)
    type: str
    homeorder: int
    dispdeg: int
    ordrdeg: int
    linedeg: int
    linereg: int
    ystep: int
    ysum: int
    wxdeg: int
    wydeg: int
    sxdeg: int
    sydeg: int


def read_modefile(path=None):
    """Positional, as MC_READMODEFILE is: the order of entries is the format."""
    path = path or os.path.join(INSTRUMENT_DIR, "XD.dat")
    val = [v for _, v in _pairs(path)]
    i = 0

    def nxt():
        nonlocal i
        i += 1
        return val[i - 1]

    modename = nxt()
    rotation = int(nxt())
    orders = fsextract_index(nxt())
    n = len(orders)
    ps = float(nxt())
    slith_pix = int(nxt())
    slith_pix_range = [int(t) for t in nxt().split()]
    slith_arc = float(nxt())
    rpppix = float(nxt())
    step = int(nxt())
    flatfrac = float(nxt())
    comwin = int(nxt())
    edgedeg = int(nxt())
    norm_nxg = int(nxt())
    norm_nyg = int(nxt())
    oversamp = float(nxt())
    norm_ybuffer = int(nxt())
    fixed = nxt()
    guesspos = np.array([[int(t) for t in nxt().split()] for _ in range(n)])
    findrange = np.array([[int(t) for t in nxt().split()] for _ in range(n)])
    return ModeInfo(modename, rotation, orders, ps, slith_pix, slith_pix_range,
                    slith_arc, rpppix, step, flatfrac, comwin, edgedeg,
                    norm_nxg, norm_nyg, oversamp, norm_ybuffer, fixed,
                    guesspos, findrange, nxt(), int(nxt()), int(nxt()),
                    int(nxt()), int(nxt()), int(nxt()), int(nxt()),
                    int(nxt()), int(nxt()), int(nxt()), int(nxt()), int(nxt()))


def _idl_int(text):
    """FIX() / LONG() of a string: the number read, truncated toward zero ('0.1' is 0)."""
    t = str(text).strip()
    try:
        return int(t)
    except ValueError:
        return int(float(t))


def _idl_float(text):
    """FLOAT() of a string: the single-precision value (returned as a Python float)."""
    return float(np.float32(float(str(text).strip())))


def _idl_ints(text):
    """FIX(STRSPLIT(val, ' ', /EXTRACT)): the on/off switches, [value, sensitive]."""
    return [_idl_int(t) for t in str(text).split()]


def _fix(v):
    """IDL FIX() / LONG() of a string: the number, truncated toward zero."""
    return int(float(v.split()[0]))


def _f32(v):
    """IDL FLOAT() of a string, or MC_CFLD(...,4) of a widget showing it."""
    return np.float32(float(v.split()[0]))


@dataclass
class InstrInfo:
    """
    The parts of ts4_soar.dat the reduction reads, each in the type the IDL
    uses. MC_READINSTRFILE reads some entries with FIX() or LONG() and
    keeps the others as text, which become a widget's initial value and are
    then read back by MC_CFLD as FLOAT (type 4) or LONG (type 3). So
    TRACEsIGTHRESH=0.1 goes through FIX() and is 0, PSPSFRAD=2.2 is the
    FLOAT 2.2, PSBGDEG is the FLOAT 1.0 and LINCORMAX=80000.0 the LONG 80000.
    """
    ncols: int = 2048
    nrows: int = 1024
    nint: int = 4
    suffix: str = ".fits"
    reductionmode: str = "A-B"
    combstat: str = "Robust Mean (Std Error)"
    combthresh: np.float32 = np.float32(8.0)        # text -> MC_CFLD(...,4)
    ybuffer: int = 1                                # text -> MC_CFLD(...,3)
    oversamp: str = "2"                             # text
    atmosthresh: np.float32 = np.float32(0.4)       # text -> MC_CFLD(...,4)
    psnaps: int = 2                                 # FIX
    pspsfrad: np.float32 = np.float32(2.2)          # text -> MC_CFLD(...,4)
    psaprad: np.float32 = np.float32(1.0)           # text -> MC_CFLD(...,4)
    psbgsub: int = 1                                # FIX
    psbgstart: np.float32 = np.float32(2.5)         # text -> MC_CFLD(...,4)
    psbgwidth: np.float32 = np.float32(5.0)         # text -> MC_CFLD(...,4)
    psbgdeg: np.float32 = np.float32(1.0)           # text -> MC_CFLD(...,4)
    xsbgsub: int = 1                                # FIX
    xsbg: str = "0-2,13-15"                         # text
    xsbgdeg: int = 0                                # text -> MC_CFLD(...,3)
    tracedeg: int = 1                               # FIX
    tracestep: int = 5                              # FIX
    tracesumap: int = 5                             # FIX
    tracesigthresh: int = 0                         # FIX('0.1')
    tracewinthresh: int = 3                         # FIX
    badpixelthresh: np.float32 = np.float32(7.0)    # FLOAT
    lincormax: int = 80000                          # LONG('80000.0')
    flatfield: int = 1
    errorprop: int = 1
    fixbadpixels: int = 0
    optimalextraction: int = 0
    keywords: list = field(default_factory=list)
    # The rest of MC_READINSTRFILE's structure: text kept as text (STRTRIM),
    # FIX/LONG as int, FLOAT as float32; the two-value switches in full.
    instr: str = "TS4_SOAR"
    stdimage: int = 1024
    plotwinsize: list = field(default_factory=lambda: [1024, 512])
    irafname: str = "FILENAME"
    exptime: str = "EXPTIME"
    time: str = "UT"
    posangle: str = "ROTANG"
    ha: str = "HA"
    airmass: str = "AIRMASS"
    nsuffix: int = 6
    bdpxmk: str = "fakebmask.fits"                  # BADPIXMASK, may be 'None'
    calmodule: str = "kna_ts41cals2d"
    filereadmode: str = "Index"
    iprefix: str = "spc"
    oprefix: str = "spectra"
    fitsreadprogram: str = "kna_readts4fits"
    combmode: str = "A"
    combodir: str = "proc/"
    skystat: str = "Robust Weighted Mean"
    skythresh: np.float32 = np.float32(8.0)         # text -> MC_CFLD(...,4)
    checkseeing: int = 0
    seeingthresh: np.float32 = np.float32(3.0)
    switches: dict = field(default_factory=lambda: {
        "AMPCOR": [0, 1], "LINCORRECT": [0, 0], "ERRORPROPAGATION": [1, 1],
        "FLATFIELD": [1, 1], "FIXBADPIXELS": [0, 0], "OPTIMALEXTRACTION": [0, 0]})


def read_instrfile(path=None):
    """MC_READINSTRFILE, each value typed as the IDL ends up using it (see InstrInfo)."""
    path = path or os.path.join(INSTRUMENT_DIR, "ts4_soar.dat")
    info = InstrInfo()
    fix = {"NCOLS": "ncols", "NROWS": "nrows", "NINT": "nint", "PSNAPS": "psnaps",
           "PSBGSUB": "psbgsub", "XSBGSUB": "xsbgsub", "TRACEDEG": "tracedeg",
           "TRACESTEP": "tracestep", "TRACESUMAP": "tracesumap",
           "TRACESIGTHRESH": "tracesigthresh", "TRACEWINTHRESH": "tracewinthresh",
           "LINCORMAX": "lincormax", "YBUFFER": "ybuffer", "XSBGDEG": "xsbgdeg",
           "FLATFIELD": "flatfield", "ERRORPROPAGATION": "errorprop",
           "FIXBADPIXELS": "fixbadpixels", "OPTIMALEXTRACTION": "optimalextraction"}
    flt = {"COMBTHRESH": "combthresh", "ATMOSTHRESH": "atmosthresh", "PSPSFRAD": "pspsfrad",
           "PSAPRAD": "psaprad", "PSBGSTART": "psbgstart", "PSBGWIDTH": "psbgwidth",
           "PSBGDEG": "psbgdeg", "BADPIXELTHRESH": "badpixelthresh"}
    txt = {"XSBG": "xsbg", "COMBSTAT": "combstat", "OVERSAMP": "oversamp",
           "REDUCTIONMODE": "reductionmode", "INSTRUMENT": "instr", "FILENAME": "irafname",
           "EXPTIME": "exptime", "TIME": "time", "POSANGLE": "posangle", "HA": "ha",
           "AIRMASS": "airmass", "BADPIXMASK": "bdpxmk", "CALMODULE": "calmodule",
           "FILEREADMODE": "filereadmode", "IPREFIX": "iprefix", "OPREFIX": "oprefix",
           "FITSREADPROGRAM": "fitsreadprogram", "COMBMODE": "combmode",
           "COMBODIR": "combodir", "SKYSTAT": "skystat"}
    fix.update({"STDIMAGE": "stdimage", "NSUFFIX": "nsuffix", "CHECKSEEING": "checkseeing"})
    flt.update({"SKYTHRESH": "skythresh", "SEEINGTHRESH": "seeingthresh"})
    for k, v in _pairs(path):
        K = k.upper()
        if K == "KEYWORD":
            info.keywords.append(v)
        elif K == "SUFFIX":
            info.suffix = "" if v.strip() == "None" else v.strip()
        elif K == "PLOTWINSIZE":
            info.plotwinsize = _idl_ints(v)
        elif K in fix:
            if K in info.switches:
                info.switches[K] = _idl_ints(v)
            setattr(info, fix[K], _fix(v))          # 'x y' pairs: the first is the value
        elif K in info.switches:
            info.switches[K] = _idl_ints(v)
        elif K in flt:
            setattr(info, flt[K], _f32(v))
        elif K in txt:
            setattr(info, txt[K], v)
    return info


def spextool_version():
    """The Spextool version string xspextool reads from data/Version.dat."""
    with open(os.path.join(SPEXTOOL_DIR, "Version.dat")) as fh:
        return fh.read().split()[0]


def read_badpixmask(info=None):
    """
    The bad pixel mask as xspextool loads it: BYTE() of the BADPIXMASK
    image (so a value like 0.5 becomes 0, bad), or, for BADPIXMASK=None,
    INTARR(ncols, nrows)+1. Unrotated, as stored.
    """
    from astropy.io import fits
    info = info or read_instrfile()
    if info.bdpxmk == "None":
        return np.ones((info.nrows, info.ncols), dtype=np.int16)
    m = fits.getdata(os.path.join(INSTRUMENT_DIR, info.bdpxmk))
    return (np.trunc(np.asarray(m, dtype=float)).astype(np.int64) & 255).astype(np.uint8)


@dataclass
class Defaults:
    """
    Everything xspextool's panels start with for TS4_SOAR, from the
    instrument file, as the reduction uses it (see defaults()).
    `instr` is the InstrInfo itself (every value in its IDL type);
    the switches are (on, the user may change it).
    """
    instr: InstrInfo
    ampcor: tuple
    lincorrect: tuple
    errorprop: tuple
    flatfield: tuple
    fixbadpixels: tuple
    optimalextraction: tuple
    badpixmask: str             # path of the mask, or None (BADPIXMASK=None)
    suffix: str                 # file glob after prefix + index ('*.fits*')
    combodir: str               # where combined images go ('proc/' or 'cal/')
    oversamp: float             # MC_CFLD(...,4) of the OVERSAMP text
    version: str                # Spextool version written into products ('4.1')


def defaults(path=None):
    """
    The single place the reduction's instrument defaults come from: edit
    ts4_soar.dat (or pass another `path`) and they change, as in xspextool.
    """
    i = read_instrfile(path)
    sw = lambda k: (bool(i.switches[k][0]), bool(i.switches[k][1]))
    return Defaults(instr=i, ampcor=sw("AMPCOR"), lincorrect=sw("LINCORRECT"),
                    errorprop=sw("ERRORPROPAGATION"), flatfield=sw("FLATFIELD"),
                    fixbadpixels=sw("FIXBADPIXELS"), optimalextraction=sw("OPTIMALEXTRACTION"),
                    badpixmask=None if i.bdpxmk == "None" else os.path.join(INSTRUMENT_DIR, i.bdpxmk),
                    suffix=i.suffix, combodir=i.combodir, oversamp=_f32(i.oversamp),
                    version=spextool_version())


# Values the TS4.1 calibration code sets for itself, not read from a file
# (ts41cals_mkflat in kna_ts41cals2d.pro).
SLITW_ARC = 1.1
SLITW_PIX = 2.4
RESOLVINGPOWER = 3500


def atran_file(rp):
    """
    The ATRAN model xspextool picks for a resolving power.

    Rounded as XSPEXTOOL_LOADIMAGES rounds it: RP 3500 (what the TS4.1 flat
    records) becomes 4000.
    """
    from triplespec_reducer.idl import iround
    if rp >= 10000.0:
        trp = iround(rp / 10000.) * 10000
    elif rp >= 1000.0:
        trp = iround(rp / 1000.) * 1000
    elif rp >= 100.0:
        trp = iround(rp / 100.) * 100
    elif rp >= 10.0:
        trp = 100
    else:
        trp = 2000
    return os.path.join(SPEXTOOL_DIR, f"atran{trp}.fits")


def read_atran(rp):
    from astropy.io import fits
    spec = fits.getdata(atran_file(rp)).astype(float)
    return spec[0], spec[1]
