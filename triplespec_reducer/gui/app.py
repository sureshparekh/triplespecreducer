#!/usr/bin/env python3
"""
TripleSpec Reducer -- the desktop interface.

    triplespec-reducer            (or python -m triplespec_reducer)

The SOAR TripleSpec 4.1 reduction, in the order the IDL Spextool tool does
it: flat, wavelength calibration, extraction, combination, telluric
correction, merging, cleaning. The night is read and organized first, so
nothing has to be typed as a frame number; every step can then be run as
proposed, adjusted, or skipped, and every product is written as
<source>_<step>.fits under the output folder.

The look is BRAIN's: one step bar over the step's controls and a stage
that shows what the step is working on.
"""

import json
import os
import sys

import numpy as np

from PySide6.QtCore import Qt, QTimer, QSettings
from PySide6.QtGui import QIcon, QPixmap
from PySide6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                               QLabel, QPushButton, QLineEdit, QComboBox, QCheckBox, QSpinBox,
                               QDoubleSpinBox, QFileDialog, QTableWidget, QTableWidgetItem,
                               QHeaderView, QAbstractItemView, QStackedWidget, QPlainTextEdit,
                               QListWidget, QListWidgetItem, QMessageBox, QSizePolicy)

from triplespec_reducer.gui import theme as TH
from triplespec_reducer.gui import help as H
from triplespec_reducer.gui.widgets import (Segmented, Disclosure, Plate, Section, label, row,
                                            vdivider, failure_dialog)
from triplespec_reducer.gui.shell import Workspace, READY, TODO, ATTENTION, OFF
from triplespec_reducer.gui.stage import Stage
from triplespec_reducer.gui.worker import Job

from triplespec_reducer import night as NT
from triplespec_reducer.pipeline import Reducer, Options

APP_NAME = "TripleSpec Reducer"
ORG = "SOAR TripleSpec Reducer"
LOGO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "logo")


def saved_theme():
    """'dark' or 'light': the theme last used (dark the first time)."""
    return str(QSettings(ORG, APP_NAME).value("theme", "dark"))


def logo_pixmap(name, size, widget=None):
    """A logo image at `size` logical pixels, sharp on high-density screens."""
    dpr = widget.devicePixelRatioF() if widget is not None else 1.0
    pm = QPixmap(os.path.join(LOGO_DIR, name)).scaled(
        int(size * dpr), int(size * dpr), Qt.KeepAspectRatio, Qt.SmoothTransformation)
    pm.setDevicePixelRatio(dpr)
    return pm
ORDER_COLOURS = ["p_star", "p_gas", "p_fit", "p_neb", "p_spec_lines"]


# ------------------------------------------------------------- small builders

def spin(lo, hi, value, step=0.1, decimals=3, tip=""):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setSingleStep(step)
    s.setValue(value)
    s.setToolTip(tip)
    return s


def ispin(lo, hi, value, tip=""):
    s = QSpinBox()
    s.setRange(lo, hi)
    s.setValue(value)
    s.setToolTip(tip)
    return s


def button(text, kind=None, tip=""):
    b = QPushButton(text)
    if kind:
        b.setObjectName(kind)
    b.setToolTip(tip)
    return b


def ranges_text(r):
    return ", ".join(f"{a:g}-{b:g}" for a, b in r)


def parse_ranges(text):
    out = []
    for part in str(text).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        a, b = part.split("-") if part.count("-") == 1 else part.rsplit("-", 1)
        out.append([float(a), float(b)])
    return out


def wrapped(text, width=74):
    import textwrap
    return "\n".join("\n".join(textwrap.wrap(p, width) or [""]) for p in str(text).split("\n"))


def caption(text, warn=False, width=74):
    """
    A paragraph under a control. Wrapped here, to a fixed measure, because
    a word-wrapping QLabel in a form reports a height for a width it has
    not been given yet, and the paragraph is clipped by the row above.
    """
    import textwrap
    lines = []
    for para in str(text).split("\n"):
        lines += textwrap.wrap(para, width) or [""]
    l = QLabel("\n".join(lines))
    l.setObjectName("CaptionWarn" if warn else "Caption")
    l.setWordWrap(False)
    return l


class TSR(QMainWindow):

    def __init__(self):
        super().__init__()
        self.C = TH.DARK
        self.night = None
        self.opts = Options()
        self.reducer = None
        self.job = None
        self.tell = {}            # science source -> telluric state (xtellcor session)
        self.merged = {}          # source -> merged spectrum
        self.cl = None            # the open XCLEANSPEC session
        self.setWindowTitle(f"{APP_NAME}  SOAR TripleSpec 4.1")
        self.resize(1560, 980)
        self.setMinimumSize(1200, 760)
        self._build()
        self.apply_theme(TH.LIGHT if saved_theme() == "light" else TH.DARK)
        self._update_steps()
        self.centralWidget().setFocusPolicy(Qt.ClickFocus)
        QTimer.singleShot(0, self.centralWidget().setFocus)

    # ================================================================ frame
    def _build(self):
        root = QWidget()
        self.setCentralWidget(root)
        v = QVBoxLayout(root)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)
        v.addWidget(self._header())
        self.stage = Stage()
        self.stage.span.connect(self._span)
        self.stage.hover.connect(lambda t: self.lbl_cursor.setText(t))
        steps = [
            ("Night", "Point at the night's raw frames. They are read, sorted and paired "
                      "for you; check the plan and change anything that is wrong.", "folder_open"),
            ("Calibrate", "The flat from the lamp frames, and a wavelength solution from each "
                          "target's sky lines (or the arc lamp).", "sliders_horizontal"),
            ("Extract", "One spectrum from every pair of frames: A minus B, divided by the "
                        "flat, apertures found, traced and summed.", "scan_eye"),
            ("Combine", "Every spectrum of a target scaled to one another and combined.",
             "layers"),
            ("Telluric", "The A0 V standard removes the atmosphere and calibrates the flux.",
             "wand_sparkles"),
            ("Merge", "The five orders joined into one spectrum.", "chart_area"),
            ("Clean & save", "Optional smoothing and repairs, then the final file.", "save"),
        ]
        self.ws = Workspace(steps, self.stage)
        self.ws.step_changed.connect(self._step_changed)
        for i, fn in enumerate((self._page_night, self._page_cal, self._page_extract,
                                self._page_combine, self._page_telluric, self._page_merge,
                                self._page_clean)):
            self._about(self.ws.page(i), i)
            fn(self.ws.page(i))
            self.ws.page(i).finish()
        for i, kind in enumerate(("required", "required", "required", "required", "optional",
                                  "required", "optional")):
            self.ws.page(i).set_tag(kind)
        self.pages = QStackedWidget()
        self.pages.addWidget(self.ws)
        self.pages.addWidget(self._products_page())
        v.addWidget(self.pages, 1)
        v.addWidget(self._statusbar())

    # ================================================================ explanations
    def _h(self, form, key, widget):
        """
        A control with its explanation from help.py: the label (or a check
        box's text), the full explanation as a tooltip on both, and a
        one-line caption under it.
        """
        d = H.field(key)
        tip = H.tooltip(key)
        widget.setToolTip(tip)
        if isinstance(widget, QCheckBox):
            widget.setText(d["label"])
            form.addRow("", widget)
        else:
            lab = label(d["label"])
            lab.setToolTip(tip)
            form.addRow(lab, widget)
        if d["caption"]:
            form.addRow(caption(d["caption"]))

    def _about(self, page, i):
        """'About this step': why it is done, what it does, what to check."""
        st = H.STEPS[i]
        p = Plate("About this step")

        def full(text):
            # prose across the whole panel, not indented to the field column
            holder = QWidget()
            lay = QHBoxLayout(holder)
            lay.setContentsMargins(0, 0, 0, 0)
            lay.addWidget(caption(text, width=92))
            return holder
        p.form.addRow(full(st["why"]))
        how = Disclosure("What it does, in order", len(st["how"]))
        for n, t in enumerate(st["how"], 1):
            how.add(caption(f"{n}.  {t}", width=92))
        p.form.addRow(how)
        chk = Disclosure("What to check before moving on", len(st["check"]))
        for t in st["check"]:
            chk.add(caption(f"-  {t}", width=92))
        p.form.addRow(chk)
        p.form.addRow(full("The same step in Spextool: " + st["idl"]))
        page.add(p)
        self.about_panels = getattr(self, "about_panels", []) + [(how, chk)]
        self.about_plates = getattr(self, "about_plates", []) + [p]

    def _header(self):
        h = QWidget()
        h.setObjectName("Header")
        h.setAttribute(Qt.WA_StyledBackground, True)
        h.setFixedHeight(58)
        lay = QHBoxLayout(h)
        lay.setContentsMargins(16, 0, 16, 0)
        lay.setSpacing(10)
        mark = QLabel()
        mark.setPixmap(logo_pixmap("app_icon_round_64.png", 36, self))
        mark.setToolTip("SOAR TripleSpec Reducer")
        lay.addWidget(mark)
        word = QLabel(APP_NAME)
        word.setObjectName("Wordmark")
        word.setFont(TH.wordmark())
        lay.addWidget(word)
        lay.addSpacing(6)
        lay.addWidget(vdivider(20))
        lay.addSpacing(6)
        self.seg_page = Segmented(["Reduce", "Products"], init=0)
        self.seg_page.setFixedWidth(260)
        self.seg_page.setFixedHeight(38)
        self.seg_page.changed.connect(lambda n: self.pages.setCurrentIndex(0 if n == "Reduce" else 1))
        lay.addWidget(self.seg_page)
        lay.addStretch(1)
        self.lbl_status = QLabel("ready")
        self.lbl_status.setObjectName("Hint")
        self.lbl_status.setFont(TH.mono(TH.SIZE["label"]))
        lay.addWidget(self.lbl_status)
        self.btn_theme = QPushButton()
        self.btn_theme.setObjectName("Tool")
        self.btn_theme.setFixedSize(34, 32)
        self.btn_theme.setToolTip("Light / dark")
        self.btn_theme.clicked.connect(lambda: self.apply_theme(TH.LIGHT if self.C["name"] == "dark" else TH.DARK))
        lay.addWidget(self.btn_theme)
        self.btn_run = QPushButton("  Reduce everything")
        self.btn_run.setObjectName("Primary")
        self.btn_run.setToolTip("Every step, for every science target, as the plan stands")
        self.btn_run.clicked.connect(self._run_all)
        lay.addWidget(self.btn_run)
        return h

    def _statusbar(self):
        b = QWidget()
        b.setObjectName("StatusBar")
        b.setAttribute(Qt.WA_StyledBackground, True)
        b.setFixedHeight(34)
        h = QHBoxLayout(b)
        h.setContentsMargins(14, 0, 14, 0)
        self.lbl_info = QLabel("")
        self.lbl_info.setObjectName("Hint")
        self.lbl_info.setFont(TH.mono(TH.SIZE["label"]))
        h.addWidget(self.lbl_info)
        h.addStretch(1)
        self.lbl_cursor = QLabel("")
        self.lbl_cursor.setObjectName("Value")
        self.lbl_cursor.setFont(TH.mono(TH.SIZE["label"]))
        h.addWidget(self.lbl_cursor)
        return b

    def apply_theme(self, C):
        self.C = C
        QSettings(ORG, APP_NAME).setValue("theme", C["name"])     # the next splash matches
        app = QApplication.instance()
        app.setPalette(TH.palette(C))
        app.setStyleSheet(TH.stylesheet(C))
        TH.apply_matplotlib(C)
        self.ws.set_palette(C)
        self.stage.set_palette(C)
        for seg in self.findChildren(Segmented):
            seg.restyle()
        for d in self.findChildren(Disclosure):
            d.set_palette(C)
        self.btn_theme.setIcon(TH.icon("sun" if C["name"] == "dark" else "moon", C["text_dim"], 17))
        self.btn_run.setIcon(TH.icon("play", C["ink"] if C["name"] == "dark" else C["void"], 12))
        self._redraw()

    # ================================================================ jobs
    def _busy(self, on, text=""):
        self.btn_run.setEnabled(not on)
        self.lbl_status.setText(text or ("working" if on else "ready"))
        QApplication.setOverrideCursor(Qt.BusyCursor) if on else QApplication.restoreOverrideCursor()

    def run_job(self, fn, done, *args, label=""):
        if self.job is not None and self.job.isRunning():
            QMessageBox.information(self, APP_NAME, "A step is still running; wait for it to finish.")
            return
        self.job = Job(fn, *args)
        self.job.progress.connect(lambda t: self.lbl_status.setText(t))
        self.job.done.connect(lambda r: (self._busy(False), done(r), self._update_steps()))
        self.job.failed.connect(lambda t: (self._busy(False), failure_dialog(
            self, t, folder=self.night.out_dir if self.night else None, title=APP_NAME)))
        self._busy(True, label)
        self.job.start()

    def _reducer(self):
        self.reducer = Reducer(self.night, self.opts,
                               progress=lambda step, msg: self.job.say(step, msg) if self.job else None)
        return self.reducer

    # ================================================================ 1. night
    def _page_night(self, page):
        sec = Plate("Folders")
        f = sec.form
        self.ed_raw = QLineEdit()
        self.ed_raw.setPlaceholderText("the folder of raw SPEC_*.fits frames")
        b1 = button("Browse", "Ghost")
        b1.clicked.connect(lambda: self._browse(self.ed_raw, out=False))
        self._h(f, "night.raw", row(self.ed_raw, b1, stretch={0: 1}))
        self.ed_out = QLineEdit()
        self.ed_out.setPlaceholderText("where cal/, proc/ and night.json are written")
        b2 = button("Browse", "Ghost")
        b2.clicked.connect(lambda: self._browse(self.ed_out, out=True))
        self._h(f, "night.out", row(self.ed_out, b2, stretch={0: 1}))
        self.cb_simbad = QCheckBox("Look up the standards' B, V and type on SIMBAD")
        self.cb_simbad.setChecked(True)
        self._h(f, "night.simbad", self.cb_simbad)
        self.btn_read = button("Read the night", tip="Read every header and propose the plan")
        self.btn_read.clicked.connect(self._read_night)
        self.btn_loadplan = button("Open a saved plan", "Ghost")
        self.btn_loadplan.clicked.connect(self._open_plan)
        f.addRow("", row(self.btn_read, self.btn_loadplan, None))
        f.addRow(caption("The raw folder is only read. Copies are never needed: nothing in it is "
                         "written or moved."))
        page.add(sec)

        sq = Plate("Targets and their plan")
        self.tbl_seq = QTableWidget(0, 11)
        self.tbl_seq.setHorizontalHeaderLabels(["Source", "Role", "Source type", "Pattern", "Frames",
                                                "Pairs", "Airmass", "Wavecal from", "Standard", "B", "V"])
        self.tbl_seq.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.tbl_seq.setHorizontalScrollMode(QAbstractItemView.ScrollPerPixel)
        self.tbl_seq.verticalHeader().setVisible(False)
        self.tbl_seq.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tbl_seq.setMinimumHeight(260)
        self.tbl_seq.itemSelectionChanged.connect(self._seq_selected)
        self.tbl_seq.itemChanged.connect(self._seq_edited)
        sq.form.addRow(self.tbl_seq)
        self.lbl_notes = caption("")
        sq.form.addRow(self.lbl_notes)
        self.btn_saveplan = button("Save the plan", "Ghost")
        self.btn_saveplan.clicked.connect(self._save_plan)
        sq.form.addRow(row(None, self.btn_saveplan))
        page.add(sq)

        d = Disclosure("Every frame and what it was taken to be")
        self.tbl_frames = QTableWidget(0, 7)
        self.tbl_frames.setHorizontalHeaderLabels(["Frame", "Kind", "OBJECT", "EXPTIME", "Airmass",
                                                   "UT start", "Note"])
        self.tbl_frames.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.tbl_frames.horizontalHeader().setStretchLastSection(True)
        self.tbl_frames.verticalHeader().setVisible(False)
        self.tbl_frames.setMinimumHeight(320)
        self.tbl_frames.itemSelectionChanged.connect(self._frame_selected)
        d.add(self.tbl_frames)
        d.add(caption("Change a frame's kind in the Kind column: flat_on, flat_off, arc, object "
                      "or ignored. Read the night again to rebuild the targets from it."))
        page.add(d)

    def _browse(self, field, out):
        d = QFileDialog.getExistingDirectory(self, "Choose a folder", field.text() or os.getcwd())
        if d:
            field.setText(d)
            if not out and not self.ed_out.text():
                base = os.path.basename(os.path.normpath(d))
                self.ed_out.setText(os.path.join(os.path.dirname(os.path.dirname(os.path.normpath(d))),
                                                 "reduced", base))

    def _read_night(self):
        raw, out = self.ed_raw.text().strip(), self.ed_out.text().strip()
        if not os.path.isdir(raw):
            QMessageBox.warning(self, APP_NAME, "Choose the folder that holds the raw frames.")
            return
        if not out:
            out = os.path.join(os.path.dirname(os.path.normpath(raw)), "reduced",
                               os.path.basename(os.path.normpath(raw)))
            self.ed_out.setText(out)
        lookup = self.cb_simbad.isChecked()

        def work(say):
            return NT.organize(raw, out, lookup=lookup, progress=say)

        self.run_job(work, self._night_ready, label="reading the headers")

    def _open_plan(self):
        p, _ = QFileDialog.getOpenFileName(self, "Open a night plan", self.ed_out.text() or os.getcwd(),
                                           "night.json (night.json)")
        if p:
            self._night_ready(NT.Night.load(p))

    def _night_ready(self, night):
        self.night = night
        self.ed_raw.setText(night.raw_dir)
        self.ed_out.setText(night.out_dir)
        night.save()
        self.tell.clear()
        self.merged.clear()
        self._fill_night()
        self._fill_sources()
        self._update_steps()
        self.lbl_info.setText(f"{len(night.frames)} frames, {len(night.sequences)} targets  "
                              f"plan saved to {os.path.join(night.out_dir, 'night.json')}")

    def _fill_night(self):
        n = self.night
        t = self.tbl_seq
        t.blockSignals(True)
        t.setRowCount(len(n.sequences))
        names = [s.name for s in n.sequences]
        stds = [""] + [s.name for s in n.sequences if s.role == "standard"]
        wcs = [s.name for s in n.sequences if s.exptime >= 60] + ["arcs"]
        for r, s in enumerate(n.sequences):
            it = QTableWidgetItem(s.name)
            it.setToolTip(f"OBJECT = {s.object}\nEvery product is {s.name}_<step>.fits")
            t.setItem(r, 0, it)
            for c, (opts, val) in ((1, (["science", "standard"], s.role)),
                                   (2, (["point", "extended"], s.source)),
                                   (7, (wcs, s.wavecal_from)),
                                   (8, (stds, s.standard))):
                cb = QComboBox()
                cb.addItems(opts)
                if val in opts:
                    cb.setCurrentText(val)
                cb.currentTextChanged.connect(lambda v, r=r, c=c: self._seq_combo(r, c, v))
                t.setCellWidget(r, c, cb)
            for c, val in ((3, s.pattern), (4, str(len(s.frames))), (5, str(len(s.pairs))),
                           (6, f"{s.airmass:.3f}")):
                it = QTableWidgetItem(val)
                it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                t.setItem(r, c, it)
            for c, val in ((9, s.bmag), (10, s.vmag)):
                t.setItem(r, c, QTableWidgetItem("" if not np.isfinite(val) else f"{val:.3f}"))
        t.blockSignals(False)
        notes = list(n.warnings)
        for s in n.sequences:
            notes += [f"{s.name}: {x}" for x in s.notes]
        import textwrap
        self.lbl_notes.setText("\n".join("\n".join(textwrap.wrap(x, 74)) for x in notes)
                               if notes else "No warnings.")
        ft = self.tbl_frames
        ft.setRowCount(len(n.frames))
        for r, f in enumerate(n.frames):
            vals = [f.name, f.kind, f.object, f"{f.exptime:g}", f"{f.airmass:.3f}", f.time[11:19], f.note]
            for c, v in enumerate(vals):
                it = QTableWidgetItem(v)
                if c != 1:
                    it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                ft.setItem(r, c, it)

    def _seq_combo(self, r, c, v):
        s = self.night.sequences[r]
        if c == 1:
            s.role = v
        elif c == 2:
            s.source = v
        elif c == 7:
            s.wavecal_from = v
        elif c == 8:
            s.standard = v
        self._fill_sources()

    def _seq_edited(self, item):
        r, c = item.row(), item.column()
        s = self.night.sequences[r]
        if c == 0:
            s.name = NT.safe_name(item.text())
        elif c in (9, 10):
            try:
                val = float(item.text())
            except ValueError:
                val = float("nan")
            if c == 9:
                s.bmag = val
            else:
                s.vmag = val

    def _save_plan(self):
        if self.night:
            for f, r in zip(self.night.frames, range(self.tbl_frames.rowCount())):
                k = self.tbl_frames.item(r, 1).text()
                if k in ("flat_on", "flat_off", "arc", "object", "ignored"):
                    f.kind = k
            self.night.flats_on = [f.name for f in self.night.frames if f.kind == "flat_on"]
            self.night.flats_off = [f.name for f in self.night.frames if f.kind == "flat_off"]
            self.night.arcs = [f.name for f in self.night.frames if f.kind == "arc"]
            p = self.night.save()
            self.lbl_info.setText(f"plan saved to {p}")

    def _seq_selected(self):
        rows = {i.row() for i in self.tbl_seq.selectedIndexes()}
        if not rows or not self.night:
            return
        s = self.night.sequences[min(rows)]
        from triplespec_reducer.frames import read_ts4
        if s.pairs:
            a, b = s.pairs[0]
            d, _, _, _ = read_ts4([self.night.frame(a).path, self.night.frame(b).path], pair=True)
            self.stage.image(d[0], f"{s.name}", f"{a} minus {b}  ({s.pattern}, {s.source} source)")
        else:
            d, _, _, _ = read_ts4([self.night.frame(s.frames[0]).path])
            self.stage.image(d[0], s.name, s.frames[0])

    def _frame_selected(self):
        rows = {i.row() for i in self.tbl_frames.selectedIndexes()}
        if not rows or not self.night:
            return
        f = self.night.frames[min(rows)]
        from triplespec_reducer.frames import read_ts4
        d, _, _, _ = read_ts4([f.path])
        self.stage.image(d[0], f.name, f"{f.kind}  {f.object}  {f.exptime:g} s  rotated: orders run "
                                        "left to right, order 3 at the bottom")

    # ================================================================ 2. calibrate
    def _page_cal(self, page):
        p = Plate("Flat field")
        self.lbl_flat = caption("Read the night first.")
        p.form.addRow(self.lbl_flat)
        self.cb_flatfield = QCheckBox("Divide the images by the flat")
        self.cb_flatfield.setChecked(True)
        self.cb_flatfield.setToolTip("ts4_soar.dat FLATFIELD=1")
        self._h(p.form, "cal.flatfield", self.cb_flatfield)
        self.ed_flatname = QLineEdit("flat")
        self._h(p.form, "cal.flatname", self.ed_flatname)
        b = button("Build the flat")
        b.clicked.connect(self._make_flat)
        bs = button("Show it", "Ghost")
        bs.clicked.connect(self._show_flat)
        p.form.addRow("", row(b, bs, None))
        p.form.addRow(caption("TS41CALS_MKFLAT: the lamps-off frames are averaged and subtracted, "
                              "the orders found from their outrigger peaks, the frames scaled to "
                              "one another, averaged, and each order divided by a smooth model of "
                              "itself."))
        page.add(p)
        w = Plate("Wavelength calibration")
        self.tbl_wc = QTableWidget(0, 5)
        self.tbl_wc.setHorizontalHeaderLabels(["Solution", "Lines", "Frames", "RMS (A)", "Status"])
        self.tbl_wc.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.tbl_wc.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tbl_wc.verticalHeader().setVisible(False)
        self.tbl_wc.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tbl_wc.setMinimumHeight(180)
        self.tbl_wc.itemSelectionChanged.connect(self._show_wavecal)
        w.form.addRow(self.tbl_wc)
        b1 = button("Calibrate the selected")
        b1.clicked.connect(lambda: self._make_wavecal(selected=True))
        b2 = button("Calibrate all", "Ghost")
        b2.clicked.connect(lambda: self._make_wavecal(selected=False))
        self.cb_autox = QCheckBox("Find the offset from the template automatically (Auto X-Correlate)")
        self.cb_autox.setChecked(True)
        self.sp_offset = spin(-200, 200, 0.0, 0.1, 3, "The pixel offset of the home order (4) from "
                                                     "the template, when not found automatically")
        self.sp_offset.setEnabled(False)
        self.cb_autox.toggled.connect(lambda on: self.sp_offset.setEnabled(not on))
        self._h(w.form, "cal.autox", self.cb_autox)
        bgo = button("Measure on a line", "Ghost", "XGETOFFSET: mark one line in the template, then "
                                                   "the same line in the observed spectrum")
        bgo.clicked.connect(self._xgetoffset)
        self._h(w.form, "cal.offset", row(self.sp_offset, bgo, None))
        self.seg_wctype = Segmented(["As planned", "Sky", "CuHeAr Arc"], init=0)
        self._h(w.form, "cal.wctype", self.seg_wctype)
        self.ed_wcframes = QLineEdit()
        self.ed_wcframes.setPlaceholderText("blank: the target's own frames (or the arcs); or numbers, e.g. 120-131")
        self._h(w.form, "cal.wcframes", self.ed_wcframes)
        self.cb_shortlines = QCheckBox("Leave out lines traced on too few points (not IDL: it stops)")
        self._h(w.form, "cal.shortlines", self.cb_shortlines)
        bqa = button("Quality plots", "Ghost", "The _Wavecal and _Distortion plots and the rejected lines")
        bqa.clicked.connect(self._wavecal_qa)
        w.form.addRow(row(b1, b2, bqa, None))
        w.form.addRow(caption("Sky: a target's own OH lines, median-combined over its frames, as the "
                              "tutorial recommends (the arcs were taken hours earlier). The "
                              "standards use the solution of the target nearest in time. 'arcs' "
                              "uses the CuHeAr lamp."))
        page.add(w)

    def _make_flat(self):
        if not self.night:
            return
        self.opts.flatfield = self.cb_flatfield.isChecked()
        self.night.flat_name = self.ed_flatname.text().strip() or "flat"

        def work(say):
            r = self._reducer()
            r.force = True
            return r.flat()

        self.run_job(work, lambda p: (self._fill_cal(), self._show_flat()), label="building the flat")

    def _show_flat(self):
        if not self.night or not os.path.exists(self.night.flat_path()):
            return
        from triplespec_reducer.flat import read_flat
        fl = read_flat(self.night.flat_path())
        x = np.arange(fl["flat"].shape[1])
        ov = []
        from triplespec_reducer import idl as I
        for i, o in enumerate(fl["orders"]):
            xx = np.arange(fl["xranges"][i][0], fl["xranges"][i][1] + 1)
            for j in (0, 1):
                ov.append(dict(x=xx, y=I.poly(xx, fl["edgecoeffs"][i, j]), color="accent", lw=0.9))
        self.stage.image(fl["flat"], "Normalized flat",
                         "order edges in blue  rms " + "  ".join(f"{o}:{r:.3f}" for o, r in zip(fl["orders"], fl["rms"])),
                         overlays=ov, lo=1, hi=99)

    def _wavecal_sources(self):
        n = self.night
        return list(dict.fromkeys(s.wavecal_from for s in n.sequences if s.wavecal_from))

    def _fill_cal(self):
        n = self.night
        if not n:
            return
        made = os.path.exists(n.flat_path())
        self.lbl_flat.setText(wrapped(f"{len(n.flats_on)} lamps-on frames ({n.flats_on[0] if n.flats_on else '-'} to "
                              f"{n.flats_on[-1] if n.flats_on else '-'}), {len(n.flats_off)} lamps-off. "
                              + (f"Built: cal/{os.path.basename(n.flat_path())}" if made else "Not built yet.")))
        t = self.tbl_wc
        srcs = self._wavecal_sources()
        t.setRowCount(len(srcs))
        for r, s in enumerate(srcs):
            frames = len(n.arcs) if s == "arcs" else len(n.seq(s).frames)
            path = n.wavecal_path(s)
            rms, status = "", "to do"
            logp = n.path(s, "log", "json")
            if os.path.exists(path):
                status = "done"
                if os.path.exists(logp):
                    with open(logp) as fh:
                        lg = json.load(fh).get("wavecal", {})
                    rms = f"{lg.get('rms_angstrom', float('nan')):.3f}"
            for c, v in enumerate((s, "CuHeAr" if s == "arcs" else "OH sky", str(frames), rms, status)):
                it = QTableWidgetItem(v)
                it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                t.setItem(r, c, it)

    def _make_wavecal(self, selected):
        if not self.night:
            return
        srcs = self._wavecal_sources()
        if selected:
            rows = sorted({i.row() for i in self.tbl_wc.selectedIndexes()})
            srcs = [srcs[r] for r in rows]
        if not srcs:
            return

        self.opts.wavecal_offset = None if self.cb_autox.isChecked() else self.sp_offset.value()
        self.opts.short_lines = "skip" if self.cb_shortlines.isChecked() else "abort"
        kind, nums = self.seg_wctype.value(), self.ed_wcframes.text().strip()
        if selected:
            from triplespec_reducer.instrument import fsextract_index
            for s in srcs:
                if kind and kind != "As planned":
                    self.opts.wavecal_arctype = dict(self.opts.wavecal_arctype or {}, **{s: kind})
                if nums:
                    want = set(fsextract_index(nums))
                    names = [f.name for f in self.night.frames if f.number in want]
                    self.opts.wavecal_frames = dict(self.opts.wavecal_frames or {}, **{s: names})

        def work(say):
            r = self._reducer()
            r.force = True
            if not os.path.exists(self.night.flat_path()):
                r.flat()
            for s in srcs:
                r.wavecal(s)
            return srcs

        self.run_job(work, lambda _: (self._fill_cal(), self._show_wavecal()), label="wavelength calibrating")

    def _selected_wavecal(self):
        rows = sorted({i.row() for i in self.tbl_wc.selectedIndexes()})
        return self._wavecal_sources()[rows[0]] if self.night and rows else None

    def _wavecal_qa(self):
        """The QA TS41CALS_MKWAVECAL leaves: the two plots, and the lines whose points were rejected."""
        s = self._selected_wavecal()
        if s is None:
            QMessageBox.information(self, APP_NAME, "Select a solution first.")
            return
        base = self.night.wavecal_path(s)[:-5]
        files = [base + "_Wavecal.pdf", base + "_Distortion.pdf"]
        if not os.path.exists(files[0]):
            QMessageBox.information(self, APP_NAME, "Calibrate it first (the plots are written with it).")
            return
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices
        for f in files:
            QDesktopServices.openUrl(QUrl.fromLocalFile(f))
        txt = base + "_Distortion.txt"
        if os.path.exists(txt):
            with open(txt) as fh:
                QMessageBox.information(self, APP_NAME, fh.read()[:4000] or "No line had points rejected.")

    def _xgetoffset(self):
        """XGETOFFSET: one line marked in the template, then in the observed home order."""
        s = self._selected_wavecal()
        if s is None or not os.path.exists(self.night.flat_path()):
            QMessageBox.information(self, APP_NAME, "Build the flat and select a solution first.")
            return
        from triplespec_reducer import getoffset as GO
        r = self._reducer()
        if s == "arcs":
            files, kind = r._raw(self.night.arcs), "CuHeAr Arc"
        else:
            files, kind = r._raw(self.night.seq(s).frames), "Sky"
        wh, fh, wo, fo, home = GO.home_spectra(files, r.flatinfo(), arctype=kind)
        self._xgo = dict(spectra=(wh, fh, wo, fo), clicks=[])
        self.stage.spectra([dict(series=[dict(x=wh, y=fh, color="p_fit")], ylabel="template"),
                            dict(series=[dict(x=wo, y=fo, color="p_obs")], ylabel="observed")],
                           f"Offset of order {home}: {s}", "drag across a line in the template, then the "
                           "same line in the observed spectrum", xlabel="column")
        self.stage.select_span(axis_index=0)

    def _xgetoffset_span(self, lo, hi):
        st = self._xgo
        st["clicks"].append((lo + hi) / 2)
        if len(st["clicks"]) == 1:
            self.stage.select_span(axis_index=1)
            return
        from triplespec_reducer import getoffset as GO
        wh, fh, wo, fo = st["spectra"]
        off, _ = GO.xgetoffset(wh, fh, wo, fo, st["clicks"][0], st["clicks"][1])
        self._xgo = None
        self.cb_autox.setChecked(False)
        self.sp_offset.setValue(float(off))
        self.lbl_info.setText(f"offset {float(off):.4f} px, from the line at {st['clicks'][0]:.1f}")

    def _show_wavecal(self):
        rows = sorted({i.row() for i in self.tbl_wc.selectedIndexes()})
        if not self.night or not rows:
            return
        s = self._wavecal_sources()[rows[0]]
        p = self.night.wavecal_path(s)
        if not os.path.exists(p):
            self.stage.message("Not calibrated yet.", "Press Calibrate the selected.")
            return
        from triplespec_reducer.wavecal import read_wavecal
        w = read_wavecal(p)
        self.stage.image(w["wavecal"], f"Wavelength image: {s}",
                         f"{w['header'].get('TSRTYPE', '')}  rms {w['header'].get('TSRRMS1D', float('nan')):.3f} A  "
                         f"offset {w['header'].get('TSROFFS', float('nan')):.3f} px", cmap="viridis", lo=0, hi=100)

    # ================================================================ 3. extract
    def _source_combo(self, roles=("science", "standard")):
        cb = QComboBox()
        cb.setMinimumWidth(260)
        cb.setProperty("roles", roles)
        return cb

    def _fill_sources(self):
        if not self.night:
            return
        for cb in (self.cb_ex_src, self.cb_cb_src):
            cur = cb.currentText()
            cb.blockSignals(True)
            cb.clear()
            cb.addItems([s.name for s in self.night.sequences])
            if cur:
                cb.setCurrentText(cur)
            cb.blockSignals(False)
        for cb in (self.cb_tc_src, self.cb_mg_src, self.cb_cl_src):
            cur = cb.currentText()
            cb.blockSignals(True)
            cb.clear()
            cb.addItems([s.name for s in self.night.sequences if s.role == "science"])
            if cur:
                cb.setCurrentText(cur)
            cb.blockSignals(False)
        self._ex_source_changed()
        self._tc_source_changed()

    def _page_extract(self, page):
        top = Plate("Target")
        self.cb_ex_src = self._source_combo()
        self.cb_ex_src.currentTextChanged.connect(self._ex_source_changed)
        self._h(top.form, "ex.source", self.cb_ex_src)
        self.lbl_ex = caption("")
        top.form.addRow(self.lbl_ex)
        page.add(top)

        md = Plate("How the frames are reduced")
        f = md.form
        self.seg_mode = Segmented(["A-B pairs", "Combined", "A", "A-Sky"], init=0)
        self.seg_mode.setToolTip("A-B: every pair differenced (the tutorial). Combined: all pairs "
                                 "averaged into one image first (XSPEXTOOL_COMBINEIMGS), one "
                                 "spectrum out. A: each frame alone. A-Sky: each object frame minus "
                                 "a sky made from the sky frames.")
        self._h(f, "ex.mode", self.seg_mode)
        from triplespec_reducer.combineimgs import COMBSTATS
        self.cb_cistat = QComboBox()
        self.cb_cistat.addItems(COMBSTATS)
        self.cb_cistat.setCurrentText(self.opts.combimg_stat)
        self.sp_cithr = spin(0.5, 50, self.opts.combimg_thresh, 0.5, 1, "COMBTHRESH")
        self._h(f, "ex.combine_with", row(self.cb_cistat, self.sp_cithr, None))
        self.cb_abab = QCheckBox("ABBA to ABAB, so the pairs add up")
        self.cb_abab.setChecked(True)
        self.cb_cisky = QCheckBox("Subtract the sky level column by column first")
        self._h(f, "ex.abab", self.cb_abab)
        self._h(f, "ex.cisky", self.cb_cisky)
        page.add(md)

        ps = Plate("Point source")
        f = ps.form
        self.sp_naps = ispin(1, 4, self.opts.psnaps, "Apertures to find: 2 for an A-B pair "
                                                    "(one positive, one negative); 4 for two objects")
        self._h(f, "ex.naps", self.sp_naps)
        self.seg_find = Segmented(["Auto", "Guess", "Fix"], init=0)
        self.seg_find.setToolTip("MC_FINDPEAKS: find the peaks, refine typed guesses, or use "
                                 "the typed positions as they are")
        self._h(f, "ex.find", self.seg_find)
        self.ed_pos = QLineEdit()
        self.ed_pos.setPlaceholderText("arcsec along the slit, e.g. 9.8, 28.1 (Guess and Fix)")
        self._h(f, "ex.positions", self.ed_pos)
        self.sp_aprad = spin(0.1, 19, self.opts.apradius, 0.1, 2, "Aperture radius, arcsec (PSAPRAD)")
        self._h(f, "ex.radius", self.sp_aprad)
        self.sp_bgs = spin(0.0, 19, self.opts.psbgstart, 0.1, 2, "PSBGSTART")
        self.sp_bgw = spin(0.0, 38, self.opts.psbgwidth, 0.1, 2, "PSBGWIDTH")
        self.sp_bgd = ispin(0, 7, self.opts.psbgdeg, "PSBGDEG")
        self._h(f, "ex.background", row(label("start"), self.sp_bgs, label("width"), self.sp_bgw,
                                        label("degree"), self.sp_bgd, None))
        self.cb_psbg = QCheckBox("Subtract the background and extract optimally")
        self.cb_psbg.setChecked(self.opts.psbgsub)
        self._h(f, "ex.psbg", self.cb_psbg)
        self.sp_psf = spin(0.1, 19, self.opts.psfradius, 0.1, 2, "PSF radius, arcsec (PSPSFRAD)")
        self._h(f, "ex.psf", self.sp_psf)
        self.cb_medprof = QCheckBox("One median profile for the whole order")
        self._h(f, "ex.medprofile", self.cb_medprof)
        self.cb_sumbg = QCheckBox("Sum with the background instead (not in the IDL tool for TS4)")
        self._h(f, "ex.sumbg", self.cb_sumbg)
        self.ed_psorders = QLineEdit()
        self.ed_psorders.setPlaceholderText("blank: all orders; or e.g. 3,5,6")
        self._h(f, "ex.orders", self.ed_psorders)
        self.cb_trace = QCheckBox("Trace the apertures (otherwise keep them at the found position)")
        self.cb_trace.setChecked(True)
        self._h(f, "ex.trace", self.cb_trace)
        self.ed_signs = QLineEdit()
        self.ed_signs.setPlaceholderText("blank: from the profile; or e.g. +, -  (Override Signs)")
        self._h(f, "ex.signs", self.ed_signs)
        self.cb_store = QCheckBox("Store the trace in cal/<source>_trace.dat")
        self._h(f, "ex.store", self.cb_store)
        self.ed_itrace = QLineEdit()
        self.ed_itrace.setPlaceholderText("Import Trace: a stored trace file instead of tracing")
        bi = button("Browse", "Ghost")
        bi.clicked.connect(lambda: self.ed_itrace.setText(QFileDialog.getOpenFileName(
            self, "A trace file", self.night.cal_dir if self.night else "", "*.dat")[0]))
        self._h(f, "ex.import", row(self.ed_itrace, bi, stretch={0: 1}))
        dt = Disclosure("Trace settings", 5)
        self.sp_tstep = ispin(1, 20, self.opts.tracestep, "TRACESTEP")
        self.sp_tsum = ispin(1, 39, self.opts.tracesumap, "TRACESUMAP (odd)")
        self.sp_tdeg = ispin(0, 4, self.opts.tracedeg, "TRACEDEG")
        self.sp_twin = spin(0.1, 20, self.opts.tracewinthresh, 0.5, 2, "TRACEWINTHRESH (pixels)")
        self.sp_tsig = spin(0.0, 10, self.opts.tracesigthresh, 0.05, 3, "TRACESIGTHRESH")
        dt.add(row(label("step"), self.sp_tstep, label("columns"), self.sp_tsum,
                   label("degree"), self.sp_tdeg, None))
        dt.add(row(label("window"), self.sp_twin, label("min sigma"), self.sp_tsig, None))
        f.addRow(dt)
        f.addRow(caption(H.field("ex.trace_settings")["caption"]))
        dt.setToolTip(H.tooltip("ex.trace_settings"))
        page.add(ps)

        xs = Plate("Extended source")
        g = xs.form
        self.ed_xspos = QLineEdit()
        self.ed_xspos.setPlaceholderText("blank: the peak of each pair's profile")
        self.ed_xspos.setToolTip("Aperture centres in arcsec; several separated by commas give "
                                 "spatially resolved spectra")
        self._h(g, "ex.xs_positions", self.ed_xspos)
        self.ed_xsrad = QLineEdit()
        self.ed_xsrad.setPlaceholderText("1.0  (one per position)")
        self._h(g, "ex.xs_radii", self.ed_xsrad)
        self.ed_xsbg = QLineEdit(self.opts.xs_bg)
        self.ed_xsbg.setToolTip("Background regions in arcsec, e.g. 2-6, 32-36")
        self._h(g, "ex.xs_bg", self.ed_xsbg)
        self.sp_xsbgd = ispin(0, 7, self.opts.xs_bgdeg, "XSBGDEG")
        self._h(g, "ex.xs_bgdeg", self.sp_xsbgd)
        self.ed_xsorders = QLineEdit()
        self.ed_xsorders.setPlaceholderText("blank: all orders; or e.g. 3,5,6")
        self._h(g, "ex.orders", self.ed_xsorders)
        self.cb_xsbgsub = QCheckBox("XS BG Subtraction (used only without a wavecal file)")
        self.cb_xsbgsub.setChecked(self.opts.xsbgsub)
        self._h(g, "ex.xsbgsub", self.cb_xsbgsub)
        self.cb_xsbg = QCheckBox("Subtract the background (an option; the IDL tool never does)")
        self._h(g, "ex.xs_background", self.cb_xsbg)
        g.addRow(caption("In the IDL tool's 2D mode the extended-source background is drawn but "
                         "never subtracted: MC_EXTSPEC2D passes XSBGINFO=xsbgino and the extract "
                         "step passes an undefined degree. Left off, the result matches the IDL "
                         "tool; on, the regions are fitted and removed as for point sources."))
        page.add(xs)

        act = Plate("Run")
        self.cb_nowc = QCheckBox("Extract in pixels, without a wavecal file")
        self._h(act.form, "ex.nowavecal", self.cb_nowc)
        self.sp_bdpix = spin(1, 50, self.opts.bdpixthresh, 0.5, 1, "Bad Pix Thresh (sigma)")
        self._h(act.form, "ex.bdpix", self.sp_bdpix)
        b1 = button("Preview the first pair", "Ghost", "Profiles and apertures, nothing written")
        b1.clicked.connect(self._ex_preview)
        b2 = button("Extract every pair")
        b2.clicked.connect(self._ex_run)
        act.form.addRow(row(b1, b2, None))
        self.lst_ex = QListWidget()
        self.lst_ex.setMinimumHeight(140)
        self.lst_ex.itemSelectionChanged.connect(self._ex_show_file)
        act.form.addRow(self.lst_ex)
        page.add(act)

    def _ex_opts(self):
        o = self.opts
        o.psnaps = self.sp_naps.value()
        o.findmode = {"Auto": "auto", "Guess": "guess", "Fix": "fixed"}[self.seg_find.value() or "Auto"]
        o.apradius = self.sp_aprad.value()
        o.psbgstart, o.psbgwidth, o.psbgdeg = self.sp_bgs.value(), self.sp_bgw.value(), self.sp_bgd.value()
        o.tracestep, o.tracesumap, o.tracedeg = self.sp_tstep.value(), self.sp_tsum.value(), self.sp_tdeg.value()
        o.tracewinthresh, o.tracesigthresh = self.sp_twin.value(), self.sp_tsig.value()
        t = self.ed_xspos.text().strip()
        o.xs_positions = [float(v) for v in t.split(",")] if t else None
        t = self.ed_xsrad.text().strip()
        o.xs_radii = [float(v) for v in t.split(",")] if t else None
        o.xs_bg = self.ed_xsbg.text().strip() or "0-2,13-15"
        o.xs_bgdeg = self.sp_xsbgd.value()
        o.xs_background = self.cb_xsbg.isChecked()
        m = self.seg_mode.value() or "A-B pairs"
        o.combine_images = m == "Combined"
        o.mode = {"A-B pairs": "A-B", "Combined": "A-B", "A": "A", "A-Sky": "A-Sky"}[m]
        o.combimg_stat, o.combimg_thresh = self.cb_cistat.currentText(), self.sp_cithr.value()
        o.abba2abab, o.combimg_subtract_sky = self.cb_abab.isChecked(), self.cb_cisky.isChecked()
        t = self.ed_signs.text().replace(" ", "")
        o.signs = [1 if c == "+" else -1 for c in t.split(",") if c] if t else None
        o.store_trace = self.cb_store.isChecked()
        o.import_trace = self.ed_itrace.text().strip() or None
        t = self.ed_pos.text().replace(",", " ").split()
        o.guess = [float(v) for v in t] if t else None
        if o.findmode != "auto" and not o.guess:
            raise ValueError("Guess and Fix need positions: type them in arcsec along the slit.")
        o.trace = self.cb_trace.isChecked() and o.findmode != "fixed"
        o.psbgsub = self.cb_psbg.isChecked()
        o.optimal = False if (o.psbgsub and self.cb_sumbg.isChecked()) else None
        o.psfradius = self.sp_psf.value()
        o.medprofile = self.cb_medprof.isChecked()
        orders = lambda e: [int(v) for v in e.text().replace(",", " ").split()] or None
        o.ps_orders, o.xs_orders = orders(self.ed_psorders), orders(self.ed_xsorders)
        o.xsbgsub = self.cb_xsbgsub.isChecked()
        o.no_wavecal = self.cb_nowc.isChecked()
        o.bdpixthresh = self.sp_bdpix.value()
        return o

    def _ex_opts_checked(self):
        """The settings, or None (with the reason shown) when they cannot be used."""
        try:
            return self._ex_opts()
        except ValueError as e:
            QMessageBox.information(self, APP_NAME, str(e))
            return None

    def _ex_source_changed(self, *_):
        if not self.night or not self.cb_ex_src.currentText():
            return
        s = self.night.seq(self.cb_ex_src.currentText())
        self.lbl_ex.setText(wrapped(f"{s.source} source, {s.pattern}; {len(s.pairs)} pairs: "
                            + ", ".join(f"{a[-9:-5]}-{b[-9:-5]}" for a, b in s.pairs[:10])
                            + (" ..." if len(s.pairs) > 10 else "")
                            + f".  Wavelength solution: {s.wavecal_from}."))
        self.lst_ex.clear()
        for f in self.night.frames:
            p = self.night.spectrum_path(s.name, f.name) if f.name in s.frames else None
            if p and os.path.exists(p):
                self.lst_ex.addItem(os.path.basename(p))

    def _ex_preview(self):
        if not self.night:
            return
        s = self.night.seq(self.cb_ex_src.currentText())
        o = self._ex_opts_checked()
        if o is None:
            return

        def work(say):
            from triplespec_reducer import extract as E
            r = self._reducer()
            r.flat()
            r.wavecal(s.wavecal_from)
            fl, wc = r.flatinfo(), r.wcinfo(s.wavecal_from)
            work_img, var, _, _ = E.load_image(r._raw(s.pairs[0]), fl)
            prof = E.mkspatprof2d(work_img, fl["omask"], wc["wavecal"], wc["spatcal"],
                                  fl["slith_arc"] / fl["slith_pix"], o.ybuffer, *r.atran(), o.atmosthresh)
            npk = (len(o.guess) if o.guess and o.findmode != "auto" else o.psnaps) if s.source == "point" else 1
            typed = [o.guess] * len(prof) if o.guess and s.source == "point" else None
            mode = o.findmode if s.source == "point" else "auto"
            pos, _ = E.findpeaks(prof, npk, mode=mode, guess=typed if mode == "guess" else None,
                                 fixed=typed if mode == "fixed" else None)
            return s, prof, pos

        self.run_job(work, self._ex_preview_done, label="making the profiles")

    def _ex_preview_done(self, res):
        s, prof, pos = res
        o = self.opts
        panels = []
        for i, (order, (x, p)) in enumerate(prof.items()):
            reg = []
            if s.source == "point":
                for c in pos[i]:
                    reg.append((c - o.apradius, c + o.apradius, "p_gas"))
                    reg.append((c - o.psbgstart - o.psbgwidth, c - o.psbgstart, "p_mask"))
                    reg.append((c + o.psbgstart, c + o.psbgstart + o.psbgwidth, "p_mask"))
            else:
                cs = o.xs_positions or [float(np.round(np.nanmedian(pos[:, 0]), 1))]
                rr = o.xs_radii or [o.apradius] * len(cs)
                for c, rad in zip(cs, rr):
                    reg.append((c - rad, c + rad, "p_gas"))
                for a, b in parse_ranges(o.xs_bg):
                    reg.append((a, b, "p_mask"))
            panels.append(dict(series=[dict(x=x, y=p, color=ORDER_COLOURS[i % 5], drawstyle="steps-mid")],
                               regions=reg, note=f"order {order}   peaks " +
                               ", ".join(f"{v:.2f}" for v in pos[i]), hline=0.0))
        self.stage.spectra(panels, f"Spatial profiles: {s.name}",
                           "apertures in green, background in red", xlabel="slit position (arcsec)")

    def _ex_run(self):
        if not self.night:
            return
        name = self.cb_ex_src.currentText()
        if self._ex_opts_checked() is None:
            return

        def work(say):
            r = self._reducer()
            r.force = True
            r.flat() if not os.path.exists(self.night.flat_path()) else None
            r.force = True
            return r.extract(name)

        self.run_job(work, lambda w: (self._ex_source_changed(), self.lst_ex.setCurrentRow(0)),
                     label=f"extracting {name}")

    def _ex_show_file(self):
        it = self.lst_ex.currentItem()
        if not it:
            return
        from triplespec_reducer.extract import read_spectrum
        d = read_spectrum(os.path.join(self.night.proc_dir, it.text()))
        series = []
        for k in range(d["spec"].shape[0]):
            o = d["orders"][k // d["naps"]]
            series.append(dict(x=d["spec"][k, 0], y=d["spec"][k, 1], color=ORDER_COLOURS[(k // d["naps"]) % 5],
                               label=f"order {o}" + (f" ap {k % d['naps'] + 1}" if d["naps"] > 1 else "")))
        self.stage.spectra([dict(series=series, ylabel="DN / s")], it.text(),
                           f"{d['header'].get('OBJECT', '')}  AIMAGE {d['header'].get('AIMAGE', '')}  "
                           f"SKY {d['header'].get('SKY', '')}")

    # ================================================================ 4. combine
    def _page_combine(self, page):
        p = Plate("Combine the spectra of a target")
        f = p.form
        self.cb_cb_src = self._source_combo()
        self.cb_cb_src.currentTextChanged.connect(self._cb_source_changed)
        self._h(f, "cb.source", self.cb_cb_src)
        self.cb_scaleorder = QComboBox()
        self.cb_scaleorder.addItems(["3", "4", "5", "6", "7"])
        self.cb_scaleorder.setCurrentText(str(self.opts.scaleorder))
        self.cb_scaleorder.setToolTip("XMC_SCALESPEC 'Median' over the middle 80% of this order; "
                                      "the tutorial uses order 4 or 5")
        self._h(f, "cb.scaleorder", self.cb_scaleorder)
        self.cb_scale = QCheckBox("Scale the spectra to a common level")
        self.cb_scale.setChecked(self.opts.scale_spectra)
        self._h(f, "cb.scale", self.cb_scale)
        self.seg_scalemode = Segmented(["Median", "Spectrum"], init=0)
        self.sp_scaleref = ispin(1, 99, 1, "The spectrum the others are scaled to, for 'Spectrum'")
        self._h(f, "cb.scalemode", row(self.seg_scalemode, self.sp_scaleref, None))
        self.ed_scalerange = QLineEdit()
        self.ed_scalerange.setPlaceholderText("blank: the middle 80% of the order, e.g. 1.20-1.26")
        self._h(f, "cb.scalerange", self.ed_scalerange)
        self.ed_prune = QLineEdit()
        self.ed_prune.setPlaceholderText("order: spectra, e.g. 5: 2 3   (several: 5: 2; 6: 4)")
        self._h(f, "cb.prune", self.ed_prune)
        self.ed_pmask = QLineEdit()
        self.ed_pmask.setPlaceholderText("order: ranges, e.g. 5: 1.30-1.32")
        self._h(f, "cb.pixmask", self.ed_pmask)
        from triplespec_reducer.combine import STATISTICS
        self.cb_stat = QComboBox()
        self.cb_stat.addItems(STATISTICS)
        self.cb_stat.setCurrentText(self.opts.statistic)
        self._h(f, "cb.statistic", self.cb_stat)
        self.sp_thresh = spin(0.5, 50, self.opts.combthresh, 0.5, 1, "Robust threshold (sigma)")
        self._h(f, "cb.thresh", self.sp_thresh)
        self.cb_shape = QCheckBox("Correct the spectral shapes (MC_SPECCOR)")
        self._h(f, "cb.shape", self.cb_shape)
        self.cb_cbaps = QCheckBox("Combine the apertures into one (Combine apertures: Yes)")
        self._h(f, "cb.apertures", self.cb_cbaps)
        self.cb_cbbasic = QCheckBox("Basic header (xcombspec /BASIC)")
        self._h(f, "cb.basic", self.cb_cbbasic)
        self.lst_cb = QListWidget()
        self.lst_cb.setMinimumHeight(160)
        self._h(f, "cb.spectra", self.lst_cb)
        b = button("Combine")
        b.clicked.connect(self._cb_run)
        f.addRow("", row(b, None))
        page.add(p)

    def _cb_source_changed(self, *_):
        if not self.night or not self.cb_cb_src.currentText():
            return
        name = self.cb_cb_src.currentText()
        self.lst_cb.clear()
        r = Reducer(self.night, self.opts)
        for pth in r.spectra(name):
            it = QListWidgetItem(os.path.basename(pth))
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked)
            self.lst_cb.addItem(it)
        if os.path.exists(self.night.path(name, "comb")):
            self._cb_show(name)

    def _cb_prune(self):
        """The two prune fields as XCOMBSPEC_PRUNESPEC sessions, one per order (and aperture: 5/2)."""
        out = {}
        for text, key in ((self.ed_prune.text(), "remove"), (self.ed_pmask.text(), "regions")):
            for part in text.split(";"):
                if ":" not in part:
                    continue
                head, body = part.split(":", 1)
                o, _, ap = head.strip().partition("/")
                d = out.setdefault((int(o), int(ap or 1)), dict(order=int(o), aperture=int(ap or 1)))
                d[key] = ([int(t) - 1 for t in body.replace(",", " ").split()] if key == "remove"
                          else [tuple(r) for r in parse_ranges(body)])
        return list(out.values())

    def _cb_run(self):
        if not self.night:
            return
        name = self.cb_cb_src.currentText()
        o = self.opts
        o.scaleorder = int(self.cb_scaleorder.currentText())
        o.statistic = self.cb_stat.currentText()
        o.combthresh = self.sp_thresh.value()
        o.correct_shapes = self.cb_shape.isChecked()
        o.combine_apertures = self.cb_cbaps.isChecked()
        o.combine_basic = self.cb_cbbasic.isChecked()
        o.scale_spectra = self.cb_scale.isChecked()
        o.scale_mode = self.seg_scalemode.value() or "Median"
        o.scale_spectrum = self.sp_scaleref.value() - 1 if o.scale_mode == "Spectrum" else None
        rg = parse_ranges(self.ed_scalerange.text())
        o.scale_range = tuple(rg[0]) if rg else None
        try:
            o.prune = dict(o.prune or {}, **{name: self._cb_prune()})
        except ValueError:
            QMessageBox.information(self, APP_NAME, "Prune: 'order: spectrum numbers' and 'order: ranges'.")
            return
        keep = [self.lst_cb.item(i).text() for i in range(self.lst_cb.count())
                if self.lst_cb.item(i).checkState() == Qt.Checked]

        def work(say):
            r = Reducer(self.night, o, force=True)
            r.combine(name, files=[os.path.join(self.night.proc_dir, k) for k in keep])
            with open(self.night.path(name, "log", "json")) as fh:
                notes = json.load(fh)["combine"].get("notes", [])
            return name, r.loaded, notes

        def done(res):
            if res[2]:
                QMessageBox.information(self, APP_NAME, "\n".join(res[2]))
            self._update_steps()
            self._cb_show(res[0], res[1])

        self.run_job(work, done, label=f"combining {name}")

    def _cb_show(self, name, L=None):
        from triplespec_reducer.extract import read_spectrum
        d = read_spectrum(self.night.path(name, "comb"))
        panels = []
        for i, o in enumerate(d["orders"]):
            series = []
            if L is not None:
                st = L["spc"][0][i]
                for k in range(st.shape[0]):
                    series.append(dict(x=st[k, 0], y=st[k, 1], color="p_faint", lw=0.6))
            series.append(dict(x=d["spec"][i * d["naps"], 0], y=d["spec"][i * d["naps"], 1],
                               color=ORDER_COLOURS[i % 5], lw=0.9))
            panels.append(dict(series=series, note=f"order {o}"))
        self.stage.spectra(panels, f"Combined: {name}", "scaled spectra in grey, combination in colour",
                           xlabel="wavelength (um)")

    # ================================================================ 5. telluric
    def _page_telluric(self, page):
        p = Plate("Target and standard")
        f = p.form
        self.cb_tc_src = self._source_combo(("science",))
        self.cb_tc_src.currentTextChanged.connect(self._tc_source_changed)
        self._h(f, "tc.target", self.cb_tc_src)
        self.cb_tc_std = QComboBox()
        self._h(f, "tc.standard", self.cb_tc_std)
        self.sp_b = spin(-5, 25, 0.0, 0.01, 3, "B magnitude of the standard")
        self.sp_v = spin(-5, 25, 0.0, 0.01, 3, "V magnitude of the standard")
        bl = button("Look up", "Ghost", "SIMBAD")
        bl.clicked.connect(self._tc_lookup)
        self._h(f, "tc.bv", row(self.sp_b, self.sp_v, bl, None))
        self.lbl_tc_am = caption("")
        f.addRow(self.lbl_tc_am)
        self.seg_route = Segmented(["A0 V standard", "Other standard", "Telluric file"], init=0)
        self._h(f, "tc.route", self.seg_route)
        self.sp_bbtemp = spin(1000, 50000, 9500, 100, 0, "Blackbody temperature of the standard (K)")
        self.cb_restore = QCheckBox("Restore the continuum with a blackbody")
        self.cb_restore.setChecked(True)
        self._h(f, "tc.basic", row(self.cb_restore, self.sp_bbtemp, None))
        self.ed_tellfile = QLineEdit()
        self.ed_tellfile.setPlaceholderText("a <source>_tellspec.fits written before")
        btf = button("Browse", "Ghost")
        btf.clicked.connect(lambda: self.ed_tellfile.setText(QFileDialog.getOpenFileName(
            self, "A telluric spectrum", self.night.proc_dir if self.night else "", "FITS (*.fits)")[0]
            or self.ed_tellfile.text()))
        self._h(f, "tc.tellfile", row(self.ed_tellfile, btf))
        bload = button("Load the two spectra")
        bload.clicked.connect(self._tc_load)
        bdirect = button("Correct (other standard or file)", "Ghost",
                         "XTELLCOR_BASIC or XTELLCOR_FINISH, with shifts found on ATRAN regions")
        bdirect.clicked.connect(self._tc_direct)
        f.addRow("", row(bload, bdirect, None))
        page.add(p)

        k = Plate("Kernel")
        f = k.form
        self.seg_meth = Segmented(["Deconvolution", "IP (not in IDL for TS4)"], init=0)
        self._h(f, "tc.method", self.seg_meth)
        self.sp_korder = ispin(3, 7, self.opts.kernel_order, "The order holding the H I line the "
                                                            "kernel is built on (Pa delta in 6)")
        self._h(f, "tc.order", self.sp_korder)
        self.ed_norm = QLineEdit(ranges_text(self.opts.normreg))
        self.ed_norm.setToolTip("Continuum windows either side of the line (the IDL 'n' key). "
                                "Drag on the plot with this field focused to add one.")
        self._h(f, "tc.continuum", self.ed_norm)
        self.sp_ndeg = ispin(1, 15, self.opts.normdegree, "Fit degree of the continuum (the IDL "
                                                         "default is 1; the tutorial uses 3)")
        self._h(f, "tc.degree", self.sp_ndeg)
        self.ed_line = QLineEdit(ranges_text([self.opts.linereg]))
        self.ed_line.setToolTip("The line window (the IDL 's' key)")
        self._h(f, "tc.line", self.ed_line)
        from triplespec_reducer.telluric import UNITS
        self.cb_units = QComboBox()
        self.cb_units.addItems(list(UNITS))
        self._h(f, "tc.units", self.cb_units)
        self.cb_prop = QCheckBox("Propagate the standard's uncertainty (always on in IDL)")
        self.cb_prop.setChecked(True)
        self._h(f, "tc.prop", self.cb_prop)
        self.sp_tension = spin(0.001, 100, 10.0, 1, 3, "The spline tension of the line scales "
                                                       "(fixed at 10 in IDL, which has no control for it)")
        self._h(f, "tc.tension", self.sp_tension)
        self.ed_vshift = QLineEdit()
        self.ed_vshift.setPlaceholderText("blank: as measured (km/s)")
        self._h(f, "tc.vshift", self.ed_vshift)
        bk = button("Build the kernel")
        bk.clicked.connect(self._tc_kernel)
        f.addRow("", row(bk, None))
        self.lbl_kern = caption("")
        f.addRow(self.lbl_kern)
        page.add(k)

        hl = Plate("Hydrogen lines")
        f = hl.form
        self.cb_hl_order = QComboBox()
        self.cb_hl_order.addItems(["3", "4", "5", "6", "7"])
        self.cb_hl_order.currentTextChanged.connect(self._tc_show_lines)
        self._h(f, "tc.hl_order", self.cb_hl_order)
        self.tbl_cp = QTableWidget(0, 3)
        self.tbl_cp.setHorizontalHeaderLabels(["wavelength (um)", "line", "scale"])
        self.tbl_cp.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tbl_cp.verticalHeader().setVisible(False)
        self.tbl_cp.setMinimumHeight(170)
        self.tbl_cp.itemChanged.connect(self._tc_cp_edited)
        f.addRow(self.tbl_cp)
        self.cb_estimate = QCheckBox("A drag on the plot estimates the scale there (XMC_SCALELINES 'e')")
        self._h(f, "tc.estimate", self.cb_estimate)
        self.ed_cut = QLineEdit()
        self.ed_cut.setPlaceholderText("cut regions bridged by a straight line, e.g. 1.81-1.83")
        self._h(f, "tc.cut", self.ed_cut)
        f.addRow(caption("Each control point scales Vega's line depth at that line; a spline under "
                         "tension (10) joins them. Push a point down if the corrected standard "
                         "shows a peak there, up for a dip -- the aim is a flat telluric spectrum."))
        page.add(hl)

        sh = Plate("Shift")
        f = sh.form
        self.tbl_shift = QTableWidget(0, 4)
        self.tbl_shift.setHorizontalHeaderLabels(["order", "aperture", "region (um)", "shift (pixels)"])
        self.tbl_shift.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tbl_shift.verticalHeader().setVisible(False)
        self.tbl_shift.setMinimumHeight(170)
        f.addRow(self.tbl_shift)
        b1 = button("Find every shift")
        b1.clicked.connect(self._tc_shifts)
        f.addRow("", row(b1, None))
        f.addRow(caption("XMC_FINDSHIFT: the telluric spectrum is shifted from -1.5 to +1.5 pixels "
                         "and the shift that leaves least structure in the region is kept. Regions "
                         "are chosen from ATRAN unless you type one."))
        self.cb_tcmodel = QCheckBox("Also write the model (<source>_tc_modVega.fits)")
        self._h(f, "tc.model", self.cb_tcmodel)
        bw = button("Write the corrected spectrum", "Primary")
        bw.clicked.connect(self._tc_write)
        f.addRow("", row(bw, None))
        page.add(sh)

    def _tc_source_changed(self, *_):
        if not self.night or not self.cb_tc_src.currentText():
            return
        s = self.night.seq(self.cb_tc_src.currentText())
        self.cb_tc_std.blockSignals(True)
        self.cb_tc_std.clear()
        stds = [x for x in self.night.sequences if x.role == "standard"]
        for x in stds:
            self.cb_tc_std.addItem(f"{x.name}", x.name)
        if s.standard:
            self.cb_tc_std.setCurrentText(s.standard)
        self.cb_tc_std.blockSignals(False)
        try:
            self.cb_tc_std.currentTextChanged.disconnect()
        except (RuntimeError, TypeError):
            pass
        self.cb_tc_std.currentTextChanged.connect(self._tc_std_changed)
        self._tc_std_changed()

    def _tc_std_changed(self, *_):
        if not self.night or not self.cb_tc_std.currentText():
            return
        std = self.night.seq(self.cb_tc_std.currentText())
        sci = self.night.seq(self.cb_tc_src.currentText())
        sci.standard = std.name
        if np.isfinite(std.bmag):
            self.sp_b.setValue(std.bmag)
        if np.isfinite(std.vmag):
            self.sp_v.setValue(std.vmag)
        d = std.airmass - sci.airmass
        self.lbl_tc_am.setText(wrapped(f"std airmass {std.airmass:.3f}, target {sci.airmass:.3f}, "
                               f"difference {d:+.3f}" + ("  (more than 0.1: xtellcor would warn)" if abs(d) > 0.1 else "")
                               + (f"   {std.sptype}" if std.sptype else "")))

    def _tc_lookup(self):
        if not self.night:
            return
        std = self.night.seq(self.cb_tc_std.currentText())
        r = NT.lookup_standard(std.object)
        if r:
            std.bmag, std.vmag, std.sptype = r["B"], r["V"], r["sptype"]
            self._tc_std_changed()
            self._fill_night()
        else:
            QMessageBox.information(self, APP_NAME, f"SIMBAD did not answer for {std.object}; type B and V.")

    def _tc_load(self):
        if not self.night:
            return
        sci = self.cb_tc_src.currentText()
        std = self.cb_tc_std.currentText()
        s = self.night.seq(std)
        s.bmag, s.vmag = self.sp_b.value(), self.sp_v.value()

        def work(say):
            from triplespec_reducer import telluric as T
            r = self._reducer()
            sp, tp = r.combine(std), r.combine(sci)
            return T.load(sp, tp, s.bmag, s.vmag)

        def done(st):
            self.tell[sci] = st
            self._tc_show_std()

        self.run_job(work, done, label="loading the standard and the target")

    def _tc_direct(self):
        """The two XTELLCOR routes without the kernel: XTELLCOR_BASIC or XTELLCOR_FINISH."""
        if not self.night:
            return
        route = self.seg_route.value()
        sci, stdname = self.cb_tc_src.currentText(), self.cb_tc_std.currentText()
        o = self.opts
        if route == "Other standard":
            s = self.night.seq(stdname)
            s.vmag = self.sp_v.value()
            o.basic_standards = dict(o.basic_standards or {}, **{stdname: dict(
                temp=self.sp_bbtemp.value(), restore=self.cb_restore.isChecked())})
            o.telluric_from = {k: v for k, v in (o.telluric_from or {}).items() if k != sci}
        elif route == "Telluric file":
            if not os.path.exists(self.ed_tellfile.text()):
                QMessageBox.information(self, APP_NAME, "Choose a telluric spectrum written before.")
                return
            o.telluric_from = dict(o.telluric_from or {}, **{sci: self.ed_tellfile.text()})
        else:
            QMessageBox.information(self, APP_NAME, "With an A0 V standard: load the two spectra and "
                                                    "build the kernel below.")
            return
        o.units = self.cb_units.currentText()

        def work(say):
            out = self.night.path(sci, "tc")
            if os.path.exists(out):             # only the correction is redone, not the combined spectra
                os.remove(out)
            return Reducer(self.night, o).telluric(sci)

        def done(path):
            from triplespec_reducer.extract import read_spectrum
            self._update_steps()
            d = read_spectrum(path)
            series = [dict(x=d["spec"][k, 0], y=d["spec"][k, 1], color=ORDER_COLOURS[(k // d["naps"]) % 5])
                      for k in range(d["spec"].shape[0])]
            self.stage.spectra([dict(series=series, ylabel=o.units)], f"Telluric corrected: {sci}",
                               os.path.basename(path))

        self.run_job(work, done, label=f"correcting {sci}")

    def _tc_state(self):
        st = self.tell.get(self.cb_tc_src.currentText())
        if st is None:
            QMessageBox.information(self, APP_NAME, "Load the two spectra first.")
        return st

    def _tc_show_std(self):
        st = self._tc_state()
        if st is None:
            return
        o = self.sp_korder.value()
        i = st["stdorders"].index(o)
        w, f = st["std"]["spec"][i, 0], st["std"]["spec"][i, 1]
        hv = [(h, "p_neb", n) for h, n in zip(st["hlines"], st["hnames"]) if np.nanmin(w) < h < np.nanmax(w)]
        reg = [(a, b, "p_gas") for a, b in parse_ranges(self.ed_norm.text())] + \
              [(a, b, "accent") for a, b in parse_ranges(self.ed_line.text())]
        self.stage.spectra([dict(series=[dict(x=w, y=f, color="p_obs")], regions=reg, vlines=hv,
                                 ylabel="DN / s")],
                           f"Standard, order {o}", "continuum windows in green, the line in blue; "
                           "drag to set a window into the focused field")
        self.stage.select_span()
        self._span_target = "kernel"

    def _span(self, lo, hi):
        """A drag on the plot: put the range where the current step wants one."""
        step = self.ws.current()
        if getattr(self, "_xgo", None):
            self._xgetoffset_span(lo, hi)
            return
        if step == 4:
            if self.ed_line.hasFocus():
                self.ed_line.setText(f"{lo:.4f}-{hi:.4f}")
            elif self.ed_norm.hasFocus():
                t = self.ed_norm.text().strip()
                self.ed_norm.setText((t + ", " if t else "") + f"{lo:.4f}-{hi:.4f}")
            elif self.ed_cut.hasFocus():
                t = self.ed_cut.text().strip()
                self.ed_cut.setText((t + ", " if t else "") + f"{lo:.4f}-{hi:.4f}")
            elif self.cb_estimate.isChecked():
                self._tc_estimate(lo, hi)
            else:
                rows = sorted({i.row() for i in self.tbl_shift.selectedIndexes()})
                if rows:
                    self.tbl_shift.item(rows[0], 1).setText(f"{lo:.4f}-{hi:.4f}")
            self._tc_redraw_current()
        elif step == 5:
            self.ed_mg_region.setText(f"{lo:.4f}-{hi:.4f}")
        elif step == 6:
            self.ed_cl_region.setText(f"{lo:.4f}-{hi:.4f}")

    def _tc_redraw_current(self):
        if getattr(self, "_tc_view", "std") == "std":
            self._tc_show_std()
        else:
            self._tc_show_lines()

    def _tc_kernel(self):
        st = self._tc_state()
        if st is None:
            return
        o = self.opts
        o.kernel_order = self.sp_korder.value()
        o.normreg = parse_ranges(self.ed_norm.text())
        o.normdegree = self.sp_ndeg.value()
        o.linereg = parse_ranges(self.ed_line.text())[0]
        o.tcmethod = (self.seg_meth.value() or "Deconvolution").split(" (")[0]
        units, tension, vs = self.cb_units.currentText(), self.sp_tension.value(), self.ed_vshift.text().strip()
        o.propstdunc = self.cb_prop.isChecked()
        ip = None                                  # chosen by slit width in conkernel

        def work(say):
            from triplespec_reducer import telluric as T
            T.conkernel(st, order=o.kernel_order, normreg=o.normreg, normdegree=o.normdegree,
                        linereg=o.linereg, method=o.tcmethod, ipcoeffs=ip)
            if vs:
                st["vshift"] = np.float32(float(vs))        # FLOAT fields in XTELLCOR
            st["tension"] = np.float32(tension)
            T.initial_points(st)
            T.contellspec(st, units=units)
            return st

        self.run_job(work, lambda _: self._tc_kernel_done(), label="building the kernel")

    def _tc_kernel_done(self):
        st = self._tc_state()
        if st["method"] == "Deconvolution":
            k = st["kernel_plot"]
            self.lbl_kern.setText(wrapped(f"Vega velocity {st['vshift']:.3f} km/s   line EW scale {st['scale']:.5f}   "
                                  f"max deviation {100 * st['maxdev']:.2f}%   rms {100 * st['rmsdev']:.2f}%"))
            self.stage.spectra([dict(series=[dict(x=k["w"], y=k["data"], color="p_obs", label="A0 V"),
                                             dict(x=k["w"], y=k["vega"], color="p_mask", label="Vega"),
                                             dict(x=k["w"], y=k["conv"], color="p_gas", label="convolved Vega"),
                                             dict(x=k["w"], y=k["resid"], color="p_fit", label="residual")],
                                     hline=0.0, ylabel="residual flux")],
                               "Kernel", "the residual should stay inside +-1%")
        else:
            self.lbl_kern.setText(wrapped("IP kernels built from IP_coefficients.dat"))
        self._fill_cp()
        self._fill_shift_table()

    def _fill_cp(self):
        st = self._tc_state()
        if st is None or "cpoints" not in st:
            return
        o = int(self.cb_hl_order.currentText())
        cp = st["cpoints"][o]
        self.tbl_cp.blockSignals(True)
        self.tbl_cp.setRowCount(cp.shape[0])
        for r, (w, sc) in enumerate(cp):
            name = ""
            m = np.abs(st["hlines"] - w) < 1e-6
            if m.any():
                name = st["hnames"][int(np.argmax(m))]
            for c, v in enumerate((f"{w:.6f}", name, f"{sc:.6f}")):
                it = QTableWidgetItem(v)
                if c != 2:
                    it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                self.tbl_cp.setItem(r, c, it)
        self.tbl_cp.blockSignals(False)

    def _tc_cp_edited(self, item):
        st = self._tc_state()
        if st is None or item.column() != 2:
            return
        o = int(self.cb_hl_order.currentText())
        try:
            st["cpoints"][o][item.row(), 1] = float(item.text())
        except ValueError:
            return
        self._tc_rebuild(o)

    def _tc_rebuild(self, o):
        from triplespec_reducer import telluric as T
        st = self._tc_state()
        st.setdefault("cutreg", {})[o] = parse_ranges(self.ed_cut.text())
        T.contellspec(st, units=st.get("units", "ergs s-1 cm-2 A-1"))
        self._tc_show_lines()

    def _tc_estimate(self, lo, hi):
        from triplespec_reducer import telluric as T
        st = self._tc_state()
        o = int(self.cb_hl_order.currentText())
        r = T.estimate_scale(st, o, (lo, hi))
        if r is None:
            QMessageBox.information(self, APP_NAME, "No H I line falls in that region.")
            return
        self._fill_cp()
        self._tc_rebuild(o)

    def _tc_show_lines(self, *_):
        st = self.tell.get(self.cb_tc_src.currentText()) if self.night else None
        if st is None or "tellspec" not in st:
            return
        self._tc_view = "lines"
        o = int(self.cb_hl_order.currentText())
        self._fill_cp()
        n = st["objorders"].index(o)
        t = st["tellspec"][n]
        w = t[0]
        hv = [(h, "p_neb", nm) for h, nm in zip(st["hlines"], st["hnames"]) if np.nanmin(w) < h < np.nanmax(w)]
        aw, at = st["awave"], st["atrans"]
        z = (aw > np.nanmin(w)) & (aw < np.nanmax(w))
        inv = 1. / t[1]
        sc = np.nanmedian(inv)
        self.stage.spectra([dict(series=[dict(x=aw[z], y=at[z] * sc, color="p_faint", label="ATRAN"),
                                         dict(x=w, y=inv, color="p_obs", label="1 / telluric")],
                                 vlines=hv, ylabel="arbitrary"),
                            dict(series=[dict(x=w, y=self._tc_obj_corrected(st, n), color="p_fit",
                                              label="target x telluric")], ylabel="flux")],
                           f"Telluric spectrum, order {o}", "a residual H I line shows as a peak or a dip "
                           "in the top panel", heights=[1, 1])
        self.stage.select_span()

    def _tc_obj_corrected(self, st, n):
        from triplespec_reducer.combine import interpspec
        t = st["tellspec"][n]
        k = n * st["objnaps"]
        ob = st["obj"]["spec"][k]
        y, _ = interpspec(t[0], t[1], ob[0])
        return np.interp(t[0], ob[0], ob[1] * y, left=np.nan, right=np.nan)

    def _fill_shift_table(self):
        st = self._tc_state()
        if st is None:
            return
        r = Reducer(self.night, self.opts)
        r._atran = (st["awave"], st["atrans"])
        cands = r.shift_regions(st)
        naps = st["objnaps"]
        self.tbl_shift.setRowCount(len(st["objorders"]) * naps)
        for i, o in enumerate(st["objorders"]):
            reg = (self.opts.shift_regions or {}).get(o) or (cands[o][0] if cands.get(o) else None)
            for ap in range(1, naps + 1):
                vals = (str(o), str(ap), f"{reg[0]:.4f}-{reg[1]:.4f}" if reg else "",
                        f"{st['shift'][i, ap - 1]:.4f}")
                for c, v in enumerate(vals):
                    it = QTableWidgetItem(v)
                    if c < 2:
                        it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                    self.tbl_shift.setItem(i * naps + ap - 1, c, it)
        self._shift_cands = cands

    def _tc_shifts(self):
        st = self._tc_state()
        if st is None:
            return
        from triplespec_reducer import telluric as T
        notes = []
        for r in range(self.tbl_shift.rowCount()):
            o, ap = int(self.tbl_shift.item(r, 0).text()), int(self.tbl_shift.item(r, 1).text())
            i = st["objorders"].index(o)
            txt = self.tbl_shift.item(r, 2).text().strip()
            regs = [parse_ranges(txt)[0]] if txt else []
            regs += [c for c in self._shift_cands.get(o, []) if c not in regs]
            ok = False
            for reg in regs:
                try:
                    T.getshift(st, o, ap, srange=reg)
                    self.tbl_shift.item(r, 2).setText(f"{reg[0]:.4f}-{reg[1]:.4f}")
                    ok = True
                    break
                except ValueError:
                    T.set_shift(st, o, ap, 0.0)
            if not ok:
                notes.append(f"order {o}, aperture {ap}: no region gave a shift inside +-1.5 px; 0 kept")
            self.tbl_shift.item(r, 3).setText(f"{st['shift'][i, ap - 1]:.4f}")
        if notes:
            QMessageBox.information(self, APP_NAME, "\n".join(notes))

    def _tc_write(self):
        st = self._tc_state()
        if st is None:
            return
        if "tellspec" not in st:
            QMessageBox.information(self, APP_NAME, "Build the kernel first.")
            return
        from triplespec_reducer import telluric as T
        name = self.cb_tc_src.currentText()
        std = self.night.seq(self.cb_tc_std.currentText())
        for r in range(self.tbl_shift.rowCount()):
            try:
                T.set_shift(st, int(self.tbl_shift.item(r, 0).text()), int(self.tbl_shift.item(r, 1).text()),
                            float(self.tbl_shift.item(r, 3).text()))
            except ValueError:
                pass
        T.corrected(st, propstdunc=self.opts.propstdunc)
        model = self.cb_tcmodel.isChecked()
        T.write(self.night.path(name, "tc"), st, tellpath=self.night.path(name, "tellspec"),
                std_name=std.object, model=model, vegapath=self.night.path(name, "tc_modVega") if model else None)
        Reducer(self.night, self.opts).log(name, "telluric", standard=std.name, bmag=st["bmag"],
                                           vmag=st["vmag"], vshift=st["vshift"], scale=st["scale"],
                                           maxdev=st["maxdev"], rmsdev=st["rmsdev"],
                                           shifts=st["shift"].tolist(),
                                           cpoints={str(k): v.tolist() for k, v in st["cpoints"].items()})
        self._update_steps()
        sp = st["corrected"]
        series = [dict(x=sp[k, 0], y=sp[k, 1], color=ORDER_COLOURS[(k // st["objnaps"]) % 5])
                  for k in range(sp.shape[0])]
        self.stage.spectra([dict(series=series, ylabel="erg s-1 cm-2 A-1")], f"Telluric corrected: {name}",
                           os.path.basename(self.night.path(name, "tc")))

    # ================================================================ 6. merge
    def _page_merge(self, page):
        p = Plate("Merge the orders")
        f = p.form
        self.cb_mg_src = self._source_combo(("science",))
        self._h(f, "mg.source", self.cb_mg_src)
        self.sp_anchor = ispin(3, 7, self.opts.anchor, "The order the others are joined to")
        self._h(f, "mg.anchor", self.sp_anchor)
        self.tbl_mg = QTableWidget(4, 4)
        self.tbl_mg.setHorizontalHeaderLabels(["add order", "scale", "region (um) for auto",
                                               "trim: Green|Blue Left|Right um"])
        self.tbl_mg.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tbl_mg.verticalHeader().setVisible(False)
        self.tbl_mg.setMinimumHeight(160)
        for r, o in enumerate((4, 5, 6, 7)):
            for c, v in enumerate((str(o), "1.0", "", "")):
                self.tbl_mg.setItem(r, c, QTableWidgetItem(v))
        f.addRow(self.tbl_mg)
        self.sp_mg_ap = ispin(1, 4, 1, "The aperture whose orders are merged")
        self._h(f, "mg.aperture", self.sp_mg_ap)
        self.ed_mg_region = QLineEdit()
        self.ed_mg_region.setPlaceholderText("drag on the plot: the region for the selected row's auto scale")
        self._h(f, "mg.region", self.ed_mg_region)
        self.cb_mgw = QCheckBox("Weight the overlap by the errors")
        self.cb_mgw.setChecked(True)
        self.cb_mgsum = QCheckBox("Sum the overlap instead of averaging")
        self._h(f, "mg.weighted", self.cb_mgw)
        self._h(f, "mg.sum", self.cb_mgsum)
        f.addRow(caption("A scale of 1 keeps the order as it is, as the tutorial does. 'auto' scales it "
                         "to the merged spectrum over the region (XMERGEORDERS_AUTOSCALE). Rows are "
                         "added top to bottom; trims as 'green right 1.47; blue left 1.44'."))
        b = button("Merge")
        b.clicked.connect(self._mg_run)
        f.addRow("", row(b, None))
        page.add(p)

    def _mg_run(self):
        if not self.night:
            return
        name = self.cb_mg_src.currentText()
        if not os.path.exists(self.night.path(name, "tc")):
            QMessageBox.information(self, APP_NAME, "Write the telluric-corrected spectrum first.")
            return
        rows = sorted({i.row() for i in self.tbl_mg.selectedIndexes()})
        if rows and self.ed_mg_region.text():
            self.tbl_mg.item(rows[0], 2).setText(self.ed_mg_region.text())
        seq, scales, trims, regions = [], {}, {}, {}
        for r in range(self.tbl_mg.rowCount()):
            o = int(self.tbl_mg.item(r, 0).text())
            seq.append(o)
            s = self.tbl_mg.item(r, 1).text().strip()
            scales[o] = "auto" if s.lower() == "auto" else float(s or 1)
            cuts = []
            for part in self.tbl_mg.item(r, 3).text().split(";"):
                tt = part.split()
                if len(tt) == 3:
                    cuts.append((tt[0].capitalize(), tt[1].capitalize(), float(tt[2])))
            if cuts:
                trims[o] = cuts
            if scales[o] == "auto":
                regions[o] = tuple(parse_ranges(self.tbl_mg.item(r, 2).text())[0])
        from astropy.io import fits
        from triplespec_reducer import merge as M
        with fits.open(self.night.path(name, "tc")) as h:
            d = h[0].data
            spec = d.astype(d.dtype.newbyteorder("="))         # MC_READSPEC: the stored type, as the pipeline
            hdr = h[0].header
        orders = [int(t) for t in str(hdr["ORDERS"]).split(",")]
        naps, ap = int(hdr["NAPS"]), self.sp_mg_ap.value()
        if ap > naps:
            QMessageBox.information(self, APP_NAME, f"{name} has {naps} aperture(s).")
            return
        try:
            m, used, done = M.merge_orders(spec, orders, naps=naps, aperture=ap, anchor=self.sp_anchor.value(),
                                           sequence=seq, scales=scales, trims=trims, regions=regions,
                                           nonidl_unweighted=not self.cb_mgw.isChecked(),
                                           nonidl_sum_flux=self.cb_mgsum.isChecked())
        except ValueError as e:                     # what XMERGEORDERS refuses, as it says it
            QMessageBox.information(self, APP_NAME, str(e))
            return
        out = self.night.path(name, "merged" if naps == 1 else f"merged_ap{ap}")
        M.write(out, m, hdr, used, done, aperture=ap, source=os.path.basename(self.night.path(name, "tc")))
        Reducer(self.night, self.opts).log(name, "merge" if naps == 1 else f"merge_ap{ap}",
                                           anchor=self.sp_anchor.value(), sequence=done, scales=used)
        for r in range(self.tbl_mg.rowCount()):
            o = int(self.tbl_mg.item(r, 0).text())
            self.tbl_mg.item(r, 1).setText(f"{used.get(o, 1.0):g}")
        self.merged[name] = m
        series = [dict(x=spec[k, 0], y=spec[k, 1] * used.get(orders[k], 1.0), color=ORDER_COLOURS[k % 5],
                       lw=0.6, alpha=0.6) for k in range(spec.shape[0])]
        self._update_steps()
        self.stage.spectra([dict(series=series + [dict(x=m[0], y=m[1], color="p_obs", lw=0.8, label="merged")],
                                 ylabel="erg s-1 cm-2 A-1")], f"Merged: {name}",
                           "orders in colour, the merged spectrum on top")
        self.stage.select_span()

    # ================================================================ 7. clean
    def _page_clean(self, page):
        p = Plate("Smooth and repair")
        f = p.form
        self.cb_cl_src = self._source_combo(("science",))
        self._h(f, "cl.source", self.cb_cl_src)
        self.cb_cl_file = QComboBox()
        self._h(f, "cl.file", self.cb_cl_file)
        self.sp_cl_ord = ispin(0, 7, 0, "Order (a merged spectrum has one)")
        self.sp_cl_ap = ispin(1, 4, 1, "Aperture")
        self._h(f, "cl.order", row(self.sp_cl_ord, self.sp_cl_ap, None))
        self.cb_cl_src.currentTextChanged.connect(self._cl_files)
        self.cb_cl_file.currentTextChanged.connect(lambda _: self._cl_open())
        self.sp_cl_ord.valueChanged.connect(lambda _: self._cl_show())
        self.sp_cl_ap.valueChanged.connect(lambda _: self._cl_show())
        self.seg_sm = Segmented(["None", "Savitzky-Golay", "Gaussian"], init=0)
        self._h(f, "cl.smoothing", self.seg_sm)
        self.sp_sgw = spin(3, 101, 3.6, 1, 1, "SG width in pixels (the IDL default, 1.5 slit widths, "
                                              "truncates to 3: no smoothing)")
        self.sp_sgd = ispin(1, 10, 2, "SG degree")
        self.sp_fwhm = spin(0.5, 50, 2.4, 0.1, 2, "Gaussian FWHM in pixels (default: the slit width)")
        self._h(f, "cl.sg", row(self.sp_sgw, self.sp_sgd, None))
        self._h(f, "cl.fwhm", self.sp_fwhm)
        self.seg_cl_type = Segmented(["Flux", "Error"], init=0)
        self._h(f, "cl.spectype", self.seg_cl_type)
        self.ed_cl_region = QLineEdit()
        self.ed_cl_region.setPlaceholderText("drag on the plot, then Fix or Remove")
        self._h(f, "cl.region", self.ed_cl_region)
        bf = button("Fix", "Ghost", "A straight line between the region's ends ('f')")
        bf.clicked.connect(lambda: self._cl_edit("fix"))
        br = button("Remove", "Ghost", "Blank the region ('r')")
        br.clicked.connect(lambda: self._cl_edit("remove"))
        bk = button("Keep", "Ghost", "Keep the edits so far ('s')")
        bk.clicked.connect(lambda: self._cl_do("commit"))
        bu = button("Undo", "Ghost", "Back to the last kept state ('u')")
        bu.clicked.connect(lambda: self._cl_do("undo"))
        f.addRow("", row(bf, br, bk, bu, None))
        self.btn_cl_smooth = button("Smooth", "Ghost", "Once per file opened, from the kept spectrum")
        self.btn_cl_smooth.clicked.connect(self._cl_smooth)
        bs = button("Save", "Primary", "Write the kept spectrum")
        bs.clicked.connect(self._cl_save)
        bt = button("Also write a text table", "Ghost", "SPEX2TEXT: every order, with the header")
        bt.clicked.connect(self._cl_text)
        f.addRow("", row(self.btn_cl_smooth, bs, bt, None))
        page.add(p)

    def _cl_files(self, name=None):
        """The source's spectrum files XCLEANSPEC could open."""
        name = name or self.cb_cl_src.currentText()
        if not self.night or not name:
            return
        kinds = ["merged"] + [f"merged_ap{k}" for k in range(1, 5)] + ["tc", "comb"]
        self.cb_cl_file.blockSignals(True)
        self.cb_cl_file.clear()
        self.cb_cl_file.addItems([k for k in kinds if os.path.exists(self.night.path(name, k))])
        self.cb_cl_file.blockSignals(False)
        self._cl_open()

    def _cl_open(self):
        """Open the chosen file: a new XCLEANSPEC session, its panel defaults from SLTW_PIX."""
        from triplespec_reducer import clean as CL
        self.cl = None
        kind = self.cb_cl_file.currentText()
        if self.night and kind:
            self.cl = CL.Session(self.night.path(self.cb_cl_src.currentText(), kind))
            d = self.cl.defaults
            if d:
                self.sp_sgw.setValue(float(d["sgwin"]))
                self.sp_sgd.setValue(int(d["sgdeg"]))
                self.sp_fwhm.setValue(float(d["fwhm"]))
            self.sp_cl_ord.setRange(min(self.cl.orders), max(self.cl.orders))
            self.sp_cl_ord.setValue(self.cl.orders[0])
            self.sp_cl_ap.setRange(1, self.cl.naps)
            self.btn_cl_smooth.setEnabled(True)
        self._cl_show()

    def _cl_idx(self):
        return self.cl.index(self.sp_cl_ord.value() if len(self.cl.orders) > 1 else None, self.sp_cl_ap.value())

    def _cl_show(self):
        if getattr(self, "cl", None) is None:
            self.stage.message("Nothing to clean yet.", "Merge the orders first, or write a corrected spectrum.")
            return
        s = self.cl.pspec[self._cl_idx()]
        with np.errstate(divide="ignore", invalid="ignore"):
            sn = s[1] / s[2]
        self.stage.spectra([dict(series=[dict(x=s[0], y=s[1], color="p_obs")], ylabel="flux"),
                            dict(series=[dict(x=s[0], y=sn, color="p_gas")], ylabel="S/N")],
                           f"Clean: {self.cb_cl_src.currentText()}", os.path.basename(self.cl.path),
                           heights=[3, 1])
        self.stage.select_span()

    def _cl_edit(self, what):
        rg = parse_ranges(self.ed_cl_region.text())
        if getattr(self, "cl", None) is None or not rg:
            return
        fn = self.cl.fix if what == "fix" else self.cl.remove
        fn(self._cl_idx(), rg[0][0], rg[0][1], spectype=self.seg_cl_type.value())
        self._cl_show()

    def _cl_do(self, what):
        if getattr(self, "cl", None) is not None:
            getattr(self.cl, what)()
            self._cl_show()

    def _cl_smooth(self):
        method = self.seg_sm.value()
        if getattr(self, "cl", None) is None or not method or method == "None":
            return
        try:
            self.cl.smooth(method=method, sgwin=self.sp_sgw.value(), sgdeg=self.sp_sgd.value(),
                           fwhm=self.sp_fwhm.value())
        except (ValueError, RuntimeError) as e:     # what XCLEANSPEC refuses
            QMessageBox.information(self, APP_NAME, str(e))
            return
        self.btn_cl_smooth.setEnabled(False)
        self._cl_show()

    def _cl_save(self):
        if getattr(self, "cl", None) is None:
            return
        name = self.cb_cl_src.currentText()
        out = self.night.path(name, "clean")
        self.cl.write(out)
        self._update_steps()
        self.lbl_info.setText(f"wrote {out}")

    def _cl_text(self):
        from triplespec_reducer import textexport as TX
        if getattr(self, "cl", None) is None:
            return
        name = self.cb_cl_src.currentText()
        src = self.night.path(name, "clean")
        src = src if os.path.exists(src) else self.cl.path
        written = TX.spex2text(src)
        self.lbl_info.setText("wrote " + ", ".join(written))

    # ================================================================ products
    def _products_page(self):
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(24, 16, 24, 16)
        left = QVBoxLayout()
        t = QLabel("Products")
        t.setObjectName("StepTitle")
        left.addWidget(t)
        left.addWidget(caption("Every file is <source>_<step>.fits: type a source to find all of it."))
        self.ed_prod = QLineEdit()
        self.ed_prod.setPlaceholderText("source name")
        self.ed_prod.textChanged.connect(self._fill_products)
        left.addWidget(self.ed_prod)
        self.lst_prod = QListWidget()
        left.addWidget(self.lst_prod, 1)
        h.addLayout(left, 1)
        self.txt_log = QPlainTextEdit()
        self.txt_log.setReadOnly(True)
        h.addWidget(self.txt_log, 2)
        self.lst_prod.itemSelectionChanged.connect(self._show_product)
        return w

    def _fill_products(self, *_):
        self.lst_prod.clear()
        if not self.night:
            return
        src = self.ed_prod.text().strip()
        names = [s.name for s in self.night.sequences if src.lower() in s.name.lower()] if src else \
            [s.name for s in self.night.sequences]
        for n in names:
            for step, p in self.night.products(n).items():
                self.lst_prod.addItem(f"{n}  {step}  ->  {os.path.relpath(p, self.night.out_dir)}")

    def _show_product(self):
        it = self.lst_prod.currentItem()
        if not it:
            return
        rel = it.text().split("->")[1].strip()
        p = os.path.join(self.night.out_dir, rel)
        if p.endswith(".json"):
            with open(p) as fh:
                self.txt_log.setPlainText(fh.read())
        elif p.endswith(".fits"):
            from astropy.io import fits
            self.txt_log.setPlainText(repr(fits.getheader(p)))

    # ================================================================ state
    def _step_changed(self, i):
        self._update_steps()
        if not self.night:
            return
        if i == 1:
            self._fill_cal()
        elif i == 2:
            if self.lst_ex.count():
                self.lst_ex.setCurrentRow(0)
                self._ex_show_file()
            else:
                self.stage.message("Preview the first pair, or extract every pair.",
                                   "The profiles show where the apertures and the background\n"
                                   "will be; each extracted spectrum can then be plotted here.")
        elif i == 3:
            self._cb_source_changed()
        elif i == 4:
            st = self.tell.get(self.cb_tc_src.currentText())
            if st is not None:
                self._tc_show_std()
            else:
                self.stage.message("Load the two spectra.",
                                   "The target's and the standard's combined spectra are read;\n"
                                   "the kernel is then built on the standard's H I line.")
        elif i == 5:
            self.stage.message("Merge the orders.", "Choose the scales and press Merge.")
        elif i == 6:
            self._cl_files() if self.cl is None else self._cl_show()

    def _redraw(self):
        if hasattr(self, "ws"):
            self._step_changed(self.ws.current())

    def _update_steps(self):
        n = self.night
        ws = self.ws
        if n is None:
            ws.set_state(0, TODO, "no night read")
            for i in range(1, 7):
                ws.set_state(i, TODO, "")
            return
        ws.set_state(0, READY, f"{len(n.sequences)} targets")
        flat = os.path.exists(n.flat_path())
        wcs = self._wavecal_sources()
        nwc = sum(os.path.exists(n.wavecal_path(s)) for s in wcs)
        ws.set_state(1, READY if flat and nwc == len(wcs) else TODO,
                     f"flat {'done' if flat else 'to do'}, wavecal {nwc}/{len(wcs)}")
        sci = [s for s in n.sequences if s.role == "science"]
        allseq = n.sequences
        ex = sum(1 for s in allseq if n.products(s.name).get(f"spec{n.frame(s.frames[0]).number:04d}")
                 or any(k.startswith("spec") for k in n.products(s.name)))
        ws.set_state(2, READY if ex == len(allseq) else TODO, f"{ex}/{len(allseq)} targets")
        cb = sum(os.path.exists(n.path(s.name, "comb")) for s in allseq)
        ws.set_state(3, READY if cb == len(allseq) else TODO, f"{cb}/{len(allseq)} combined")
        tc = sum(os.path.exists(n.path(s.name, "tc")) for s in sci)
        ws.set_state(4, READY if tc == len(sci) else TODO, f"{tc}/{len(sci)} corrected")
        mg = sum(os.path.exists(n.path(s.name, "merged")) for s in sci)
        ws.set_state(5, READY if mg == len(sci) else TODO, f"{mg}/{len(sci)} merged")
        cl = sum(os.path.exists(n.path(s.name, "clean")) for s in sci)
        ws.set_state(6, READY if cl else OFF, f"{cl}/{len(sci)} saved")

    def _run_all(self):
        if not self.night:
            QMessageBox.information(self, APP_NAME, "Read a night first.")
            return
        if self._ex_opts_checked() is None:
            return
        missing = [s.name for s in self.night.sequences if s.role == "standard"
                   and not (np.isfinite(s.bmag) and np.isfinite(s.vmag))]
        if missing:
            QMessageBox.information(self, APP_NAME, "Enter B and V for: " + ", ".join(missing))
            return

        def work(say):
            r = self._reducer()
            return r.run()

        self.run_job(work, lambda prods: (self.lbl_info.setText(
            f"done: {len(prods)} targets reduced into {self.night.proc_dir}"), self._fill_products()),
            label="reducing everything")


def main(argv=None):
    os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")
    app = QApplication.instance() or QApplication(argv or sys.argv)
    app.setApplicationName(APP_NAME)
    app.setWindowIcon(QIcon(os.path.join(LOGO_DIR, "app_icon.png")))
    TH.load_fonts()
    w = TSR()
    if argv and len(argv) > 1 and os.path.isdir(argv[1]):
        w.ed_raw.setText(argv[1])
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
