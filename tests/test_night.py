"""The night planner on the test night: every decision it makes, checked."""

import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
RAW = os.path.join(ROOT, "data", "raw", "2025-11-02")
pytestmark = pytest.mark.skipif(not os.path.isdir(RAW), reason="raw frames not present")

from triplespec_reducer import night as N   # noqa: E402


@pytest.fixture(scope="module")
def n(tmp_path_factory):
    return N.organize(RAW, str(tmp_path_factory.mktemp("out")), lookup=False)


def test_frames(n):
    assert len(n.flats_on) == 50 and n.flats_on[0] == "SPEC_DFLAT_0002.fits"
    assert len(n.flats_off) == 52 and len(n.arcs) == 12
    ign = {f: n.frame(f).note for f in n.ignored}
    assert "SPEC_141Tau_0001.fits" in ign
    for t in ("t20-21.fits", "t133-137.fits", "t168-169.fits", "t198-197.fits"):
        assert "not a raw frame" in ign[t]


def test_sequences(n):
    names = [s.name for s in n.sequences]
    assert names == ["HD225200_1", "NGC134", "NGC134_pos2", "HD1777",
                     "kks_03-gc1_two_objects_-_fainter", "HD225200_2",
                     "NGC_253_gc_73NGC_253_gc_73", "HD225200_3", "NGC_253_gc_213", "HD225200_4"]
    roles = {s.name: s.role for s in n.sequences}
    assert roles["NGC134"] == "science" and roles["HD1777"] == "standard"


def test_pairs(n):
    ngc = n.seq("NGC134")
    assert ngc.pattern == "object-sky" and ngc.source == "extended"
    short = [(a[-9:-5], b[-9:-5]) for a, b in ngc.pairs]
    assert short == [("0120", "0121"), ("0122", "0121"), ("0123", "0124"), ("0126", "0127"),
                     ("0128", "0127"), ("0129", "0130"), ("0131", "0130")]
    assert [(a[-9:-5], b[-9:-5]) for a, b in n.seq("NGC134_pos2").pairs] == [("0125", "0124")]
    hd = n.seq("HD225200_1")
    assert hd.pattern == "ABBA" and [(a[-9:-5], b[-9:-5]) for a, b in hd.pairs] == [("0116", "0117"),
                                                                                      ("0118", "0119")]
    gc = n.seq("NGC_253_gc_73NGC_253_gc_73")
    assert len(gc.pairs) == 16 and not gc.skipped


def test_calibration_choices(n):
    wc = {s.name: s.wavecal_from for s in n.sequences}
    assert wc["NGC134_pos2"] == "NGC134" and wc["HD225200_1"] == "NGC134"
    assert wc["HD225200_2"] == "kks_03-gc1_two_objects_-_fainter"
    std = {s.name: s.standard for s in n.sequences if s.role == "science"}
    assert std["NGC134"] == "HD225200_1"
    assert std["NGC_253_gc_73NGC_253_gc_73"] == "HD225200_2"
    assert std["NGC_253_gc_213"] == "HD225200_3"


def test_naming(n):
    p = n.spectrum_path("NGC134", "SPEC_SO2025D-006-0120.fits")
    assert p.endswith(os.path.join("proc", "NGC134_spec0120.fits"))
    assert n.path("NGC134", "tc").endswith("NGC134_tc.fits")
    assert n.wavecal_path("NGC134").endswith(os.path.join("cal", "NGC134_wavecal.fits"))


def test_roundtrip(n, tmp_path):
    p = n.save(str(tmp_path / "night.json"))
    m = N.Night.load(p)
    assert [s.pairs for s in m.sequences] == [s.pairs for s in n.sequences]
