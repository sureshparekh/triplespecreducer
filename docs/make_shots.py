"""
The manual's screenshots, made by driving the real interface offscreen
through the tutorial's reduction (NGC 134 with the A0 V star HD 225200).

    python docs/make_shots.py            # light theme, into docs/_static/gui/
    python docs/make_shots.py dark       # the dark theme too (light_NAME.png, dark_NAME.png)

It works in a copy of a reduced night (reduced/tutorial), so the products
already there are reused and the validated folder is never touched.
"""

import json
import os
import shutil
import sys
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from PySide6.QtWidgets import QApplication              # noqa: E402

from triplespec_reducer.gui import theme as TH          # noqa: E402
from triplespec_reducer.gui.app import TSR              # noqa: E402
from triplespec_reducer.night import Night              # noqa: E402

SRC = os.path.join(ROOT, "reduced", "2025-11-02-parallel")
WORK = os.path.join(ROOT, "reduced", "tutorial")
OUT = os.path.join(ROOT, "docs", "_static", "gui")
SCI, STD = "NGC134", "HD225200_1"


def prepare():
    """A copy of the reduced night whose plan points at the copy."""
    if not os.path.isdir(WORK):
        shutil.copytree(SRC, WORK)
    p = os.path.join(WORK, "night.json")
    d = json.load(open(p))
    d["out_dir"] = WORK
    json.dump(d, open(p, "w"), indent=1)
    return Night.load(p)


app = QApplication([])
TH.load_fonts()
themes = ["light"] + (["dark"] if "dark" in sys.argv[1:] else [])
night = prepare()


def settle(w, t=0.3):
    """Let Qt draw, and wait for a background step to finish."""
    end = time.time() + t
    while time.time() < end or (w.job is not None and w.job.isRunning()):
        app.processEvents()
        time.sleep(0.02)
    for _ in range(5):
        app.processEvents()


REGIONS = {}
WATCH = ("tbl_seq", "lbl_notes", "tbl_wc", "lbl_flat", "seg_mode", "sp_naps", "lst_ex", "lst_cb",
         "cb_stat", "lbl_kern", "tbl_cp", "tbl_shift", "tbl_mg", "seg_sm", "ed_raw", "lbl_tc_am",
         "btn_read", "lst_prod", "txt_log")


def regions(w):
    """Where the named widgets are in the window, clipped to what is visible."""
    from PySide6.QtCore import QPoint, QRect
    out = {}
    for name in WATCH:
        wd = getattr(w, name, None)
        if wd is None or not wd.isVisible():
            continue
        r = QRect(wd.mapTo(w, QPoint(0, 0)), wd.size())
        vis = wd.visibleRegion().boundingRect()
        if vis.isEmpty():
            continue
        r = QRect(wd.mapTo(w, vis.topLeft()), vis.size())
        out[name] = [r.x(), r.y(), r.width(), r.height()]
    return out


def shot(w, name, theme):
    settle(w)
    os.makedirs(OUT, exist_ok=True)
    w.grab().save(os.path.join(OUT, f"{theme}_{name}.png"))
    REGIONS[name] = regions(w)
    print("  ", f"{theme}_{name}.png")


def open_about(w, i, which=0):
    w.about_panels[i][which].head.setChecked(True)


def scroll_to(w, i, widget):
    """Scroll step i's panel so `widget` is in view."""
    from PySide6.QtCore import QPoint
    from PySide6.QtWidgets import QScrollArea
    page = w.ws.page(i)
    for sa in page.findChildren(QScrollArea):
        if sa.widget() is not None and sa.widget().isAncestorOf(widget):
            # the widget's panel near the top of the view, not just barely in it
            y = widget.mapTo(sa.widget(), QPoint(0, 0)).y()
            sa.verticalScrollBar().setValue(max(0, y - 120))
    app.processEvents()


for theme in themes:
    print(theme)
    w = TSR()
    w.apply_theme(TH.LIGHT if theme == "light" else TH.DARK)
    w.resize(1600, 1000)
    w.show()
    settle(w)
    shot(w, "00-start", theme)

    # 1. the night
    w.ed_raw.setText(night.raw_dir)
    w.ed_out.setText(night.out_dir)
    w._night_ready(night)
    w.tbl_seq.selectRow([s.name for s in night.sequences].index(SCI))
    shot(w, "00-folders", theme)
    scroll_to(w, 0, w.tbl_seq)
    shot(w, "01-night", theme)
    scroll_to(w, 0, w.ed_raw)
    open_about(w, 0)
    from PySide6.QtWidgets import QScrollArea
    for sa in w.ws.page(0).findChildren(QScrollArea):
        sa.verticalScrollBar().setValue(0)
    shot(w, "01-night-about", theme)

    # 2. calibrate
    w.ws.go(1)
    w._show_flat()
    scroll_to(w, 1, w.lbl_flat)
    shot(w, "02-flat", theme)
    names = [w.tbl_wc.item(r, 0).text() for r in range(w.tbl_wc.rowCount())]
    if SCI in names:
        w.tbl_wc.selectRow(names.index(SCI))
    scroll_to(w, 1, w.tbl_wc)
    shot(w, "02-wavecal", theme)

    # 3. extract: the point-source standard, then the galaxy
    w.ws.go(2)
    w.cb_ex_src.setCurrentText(STD)
    w._ex_preview()
    scroll_to(w, 2, w.sp_naps)
    shot(w, "03-extract-preview-std", theme)
    w.cb_ex_src.setCurrentText(SCI)
    w._ex_preview()
    shot(w, "03-extract-preview-sci", theme)
    if w.lst_ex.count():
        w.lst_ex.setCurrentRow(0)
    scroll_to(w, 2, w.lst_ex)
    shot(w, "03-extract-spectrum", theme)

    # 4. combine
    w.ws.go(3)
    w.cb_cb_src.setCurrentText(STD)
    w._cb_source_changed()
    w._cb_run()
    scroll_to(w, 3, w.lst_cb)
    shot(w, "04-combine-std", theme)
    w.cb_cb_src.setCurrentText(SCI)
    w._cb_source_changed()
    w._cb_run()
    scroll_to(w, 3, w.lst_cb)
    shot(w, "04-combine-sci", theme)

    # 5. telluric
    w.ws.go(4)
    w.cb_tc_src.setCurrentText(SCI)
    w.cb_tc_std.setCurrentText(STD)
    w._tc_load()
    scroll_to(w, 4, w.lbl_tc_am)
    shot(w, "05-telluric-loaded", theme)
    w._tc_kernel()
    scroll_to(w, 4, w.lbl_kern)
    shot(w, "05-telluric-kernel", theme)
    w.cb_hl_order.setCurrentText("4")
    w._tc_show_lines()
    scroll_to(w, 4, w.tbl_cp)
    shot(w, "05-telluric-lines", theme)
    w._tc_shifts()
    scroll_to(w, 4, w.tbl_shift)
    shot(w, "05-telluric-shifts", theme)
    w._tc_write()
    shot(w, "05-telluric-corrected", theme)

    # 6. merge
    w.ws.go(5)
    w.cb_mg_src.setCurrentText(SCI)
    w._mg_run()
    scroll_to(w, 5, w.tbl_mg)
    shot(w, "06-merge", theme)

    # 7. clean
    w.ws.go(6)
    w.cb_cl_src.setCurrentText(SCI)
    w._cl_show()
    scroll_to(w, 6, w.seg_sm)
    shot(w, "07-clean", theme)

    # the products
    w.seg_page.setValue("Products")
    w.pages.setCurrentIndex(1)
    w.ed_prod.setText(SCI)
    settle(w)
    if w.lst_prod.count():
        w.lst_prod.setCurrentRow(w.lst_prod.count() - 1)
    shot(w, "08-products", theme)
    w.close()
    if theme == "light":
        json.dump(REGIONS, open(os.path.join(OUT, "regions.json"), "w"), indent=1)
print("done")
