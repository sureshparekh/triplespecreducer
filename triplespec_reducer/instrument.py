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
    """MC_FSEXTRACT(/INDEX): '3-7' or '1,4-6' to a list of integers."""
    out = []
    for part in str(text).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


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


@dataclass
class InstrInfo:
    """The parts of ts4_soar.dat the reduction reads."""
    ncols: int = 2048
    nrows: int = 1024
    nint: int = 4
    suffix: str = ".fits"
    reductionmode: str = "A-B"
    combstat: str = "Robust Mean (Std Error)"
    combthresh: float = 8.0
    ybuffer: int = 1
    oversamp: int = 2
    atmosthresh: float = 0.4
    psnaps: int = 2
    pspsfrad: float = 2.2
    psaprad: float = 1.0
    psbgsub: int = 1
    psbgstart: float = 2.5
    psbgwidth: float = 5.0
    psbgdeg: int = 1
    xsbgsub: int = 1
    xsbg: str = "0-2,13-15"
    xsbgdeg: int = 0
    tracedeg: int = 1
    tracestep: int = 5
    tracesumap: int = 5
    tracesigthresh: float = 0.1
    tracewinthresh: float = 3.0
    badpixelthresh: float = 7.0
    lincormax: float = 80000.0
    flatfield: int = 1
    errorprop: int = 1
    fixbadpixels: int = 0
    optimalextraction: int = 0
    keywords: list = field(default_factory=list)


def read_instrfile(path=None):
    path = path or os.path.join(INSTRUMENT_DIR, "ts4_soar.dat")
    info = InstrInfo()
    first = lambda v: v.split()[0]
    for k, v in _pairs(path):
        K = k.upper()
        if K == "KEYWORD":
            info.keywords.append(v)
        elif K == "NCOLS":
            info.ncols = int(v)
        elif K == "NROWS":
            info.nrows = int(v)
        elif K == "NINT":
            info.nint = int(first(v))
        elif K == "COMBTHRESH":
            info.combthresh = float(v)
        elif K == "YBUFFER":
            info.ybuffer = int(v)
        elif K == "ATMOSTHRESH":
            info.atmosthresh = float(v)
        elif K == "PSNAPS":
            info.psnaps = int(v)
        elif K == "PSPSFRAD":
            info.pspsfrad = float(v)
        elif K == "PSAPRAD":
            info.psaprad = float(v)
        elif K == "PSBGSUB":
            info.psbgsub = int(v)
        elif K == "PSBGSTART":
            info.psbgstart = float(v)
        elif K == "PSBGWIDTH":
            info.psbgwidth = float(v)
        elif K == "PSBGDEG":
            info.psbgdeg = int(v)
        elif K == "XSBGSUB":
            info.xsbgsub = int(v)
        elif K == "XSBG":
            info.xsbg = v
        elif K == "XSBGDEG":
            info.xsbgdeg = int(v)
        elif K == "TRACEDEG":
            info.tracedeg = int(v)
        elif K == "TRACESTEP":
            info.tracestep = int(v)
        elif K == "TRACESUMAP":
            info.tracesumap = int(v)
        elif K == "TRACESIGTHRESH":
            info.tracesigthresh = float(v)
        elif K == "TRACEWINTHRESH":
            info.tracewinthresh = float(v)
        elif K == "BADPIXELTHRESH":
            info.badpixelthresh = float(v)
        elif K == "LINCORMAX":
            info.lincormax = float(first(v))
        elif K == "COMBSTAT":
            info.combstat = v
    return info


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
