"""
The reduction from the raw pairs to the cleaned spectrum as one chain,
against IDL's chain (validation/idl/ref_chain.pro, run by
validation/compare_chain.py): extraction -> combination -> telluric
correction -> merged orders -> cleaned spectra, every stage reading the
files of the stage before, every file bit for bit, card for card and history
line for line.

    gui   the Python defaults = the TS4_SOAR panel's (optimal extraction of
          the point sources)
    sum   optimal=False: the aperture sum with the background

The flat and wavecal are IDL's own (their parity is test_parity's). The
telluric kernel is IDL's FFT's: ours must agree to its rounding, and the
chain goes on from IDL's.
"""

import os
import sys
import warnings

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "validation"))
REF = os.path.join(ROOT, "validation", "idl_out")
RAW = os.path.join(ROOT, "data", "raw", "2025-11-02")
pytestmark = pytest.mark.skipif(
    not (os.path.isfile(os.path.join(REF, "chain_gui", "clean_sg_NGC134.fits")) and os.path.isdir(RAW)),
    reason="IDL chain references (ref_chain.pro) or raw frames not present")

import compare_chain as CC                                       # noqa: E402
from triplespec_reducer import flat as F, wavecal as W           # noqa: E402

warnings.filterwarnings("ignore", module="astropy")


@pytest.mark.parametrize("case", ["gui", "sum"])
def test_chain(case, tmp_path):
    fl = F.read_flat(os.path.join(REF, "flat.fits"))
    wc = W.read_wavecal(os.path.join(REF, "wavecal_NGC134.fits"), rotation=fl["rotation"])
    del CC.failures[:]
    CC.chain(case, fl, wc, outroot=str(tmp_path))
    assert not CC.failures, f"{case}: {CC.failures}"
