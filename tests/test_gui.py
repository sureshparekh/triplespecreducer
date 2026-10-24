"""The interface opens, takes a plan, walks every step in both themes, and does not fail."""

import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
PLAN = os.path.join(ROOT, "reduced", "2025-11-02", "night.json")


@pytest.mark.skipif(not os.path.exists(PLAN), reason="no organized night")
def test_walk_the_steps():
    from PySide6.QtWidgets import QApplication
    from triplespec_reducer.gui import theme as TH
    from triplespec_reducer.gui.app import TSR
    from triplespec_reducer.night import Night
    app = QApplication.instance() or QApplication([])
    w = TSR()
    w._night_ready(Night.load(PLAN))
    for C in (TH.DARK, TH.LIGHT):
        w.apply_theme(C)
        for i in range(7):
            w.ws.go(i)
            app.processEvents()
    assert w.tbl_seq.rowCount() == len(w.night.sequences)
    w.close()
