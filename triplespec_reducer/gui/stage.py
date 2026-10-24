"""
The stage: the one plot area every step draws into.

It follows BRAIN's: a strip naming what is shown, home / pan / zoom tools,
no box, grid or legend panel inside the data, numbers in the mono face,
and the theme's colours throughout. A step can ask for a detector image,
one spectrum panel, or a column of panels sharing a wavelength axis, and
can switch on a drag-to-select that reports the span it was given --
the same gesture as the IDL tool's 'select a region' keys.
"""

import numpy as np

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QSizePolicy

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.widgets import SpanSelector

from triplespec_reducer.gui import theme as TH
from triplespec_reducer.gui.widgets import Strip


class Stage(QWidget):
    span = Signal(float, float)
    hover = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("Stage")
        self.setAttribute(Qt.WA_StyledBackground, True)
        self.C = TH.DARK
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        bar = Strip()
        bar.setFixedHeight(40)
        h = QHBoxLayout(bar)
        h.setContentsMargins(16, 0, 12, 0)
        h.setSpacing(10)
        self.title = QLabel("")
        self.title.setObjectName("SectionTitle")
        h.addWidget(self.title)
        self.subtitle = QLabel("")
        self.subtitle.setObjectName("Hint")
        h.addWidget(bar.optional(self.subtitle, 700))
        h.addStretch(1)
        self.fig = Figure(figsize=(9, 7), dpi=100)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.setFocusPolicy(Qt.ClickFocus)
        self._nav = NavigationToolbar2QT(self.canvas, self)
        self._nav.hide()
        self._tools = []
        for kind, tip, slot, check in (("home", "Reset the view", self.home, False),
                                       ("pan", "Drag to pan", self._pan, True),
                                       ("zoom", "Drag a box to zoom in", self._zoom, True)):
            b = QPushButton()
            b.setObjectName("Tool")
            b.setFixedSize(32, 32)
            b.setCheckable(check)
            b.setToolTip(tip)
            b.clicked.connect(slot)
            h.addWidget(b)
            self._tools.append((b, kind))
        self.btn_pan, self.btn_zoom = self._tools[1][0], self._tools[2][0]

        well = QWidget()
        well.setObjectName("Well")
        well.setAttribute(Qt.WA_StyledBackground, True)
        wv = QVBoxLayout(well)
        wv.setContentsMargins(0, 0, 0, 0)
        wv.addWidget(self.canvas)
        v.addWidget(bar)
        v.addWidget(well, 1)

        self.axes = []
        self._span = None
        self._span_color = None
        self._image = None
        self.canvas.mpl_connect("motion_notify_event", self._motion)
        self.canvas.mpl_connect("scroll_event", self._scroll)
        self.message("Nothing to show yet.", "Each step draws what it is working on here.")

    # ------------------------------------------------------------- tools
    def set_palette(self, C):
        self.C = C
        for b, k in self._tools:
            b.setIcon(TH.icon(k, C["text_dim"], 16))
        self.fig.set_facecolor(C["surface"])
        for a in self.axes:
            self._style(a)
        self.canvas.draw_idle()

    def home(self):
        self._nav.home()

    def _pan(self):
        self._nav.pan()
        self._sync()

    def _zoom(self):
        self._nav.zoom()
        self._sync()

    def _sync(self):
        mode = str(self._nav.mode)
        self.btn_pan.setChecked("pan" in mode)
        self.btn_zoom.setChecked("zoom" in mode)
        if self._span is not None:
            self._span.set_active(not mode)

    def _scroll(self, ev):
        """The wheel zooms the horizontal axis about the pointer; Shift, the vertical."""
        if ev.inaxes is None or self._image is not None and ev.inaxes is self.axes[0]:
            return
        f = 1 / 1.25 if ev.button == "up" else 1.25
        a = ev.inaxes
        if "shift" in (ev.modifiers or ()):
            lo, hi = a.get_ylim()
            a.set_ylim(ev.ydata - (ev.ydata - lo) * f, ev.ydata + (hi - ev.ydata) * f)
        else:
            lo, hi = a.get_xlim()
            a.set_xlim(ev.xdata - (ev.xdata - lo) * f, ev.xdata + (hi - ev.xdata) * f)
        self.canvas.draw_idle()

    def _motion(self, ev):
        if ev.inaxes is None or ev.xdata is None:
            self.hover.emit("")
            return
        if self._image is not None and ev.inaxes is self.axes[0]:
            img = self._image
            c, r = int(round(ev.xdata)), int(round(ev.ydata))
            if 0 <= r < img.shape[0] and 0 <= c < img.shape[1]:
                self.hover.emit(f"x {c:5d}  y {r:5d}  {img[r, c]:12.5g}")
            return
        self.hover.emit(f"{ev.xdata:10.5f}  {ev.ydata:12.5g}")

    # ------------------------------------------------------------- layout
    def _style(self, a):
        C = self.C
        a.set_facecolor(C["void"])
        for name, sp in a.spines.items():
            sp.set_visible(name in ("left", "bottom"))
            sp.set_color(C["border_2"])
        a.tick_params(colors=C["text_mute"], labelsize=9, length=3, width=1.0)
        for lab in a.get_xticklabels() + a.get_yticklabels():
            lab.set_fontfamily("monospace")
        a.xaxis.label.set_color(C["text_dim"])
        a.yaxis.label.set_color(C["text_dim"])
        a.title.set_color(C["text_dim"])

    def _reset(self, n=1, sharex=True, heights=None):
        self.fig.clear()
        self.fig.set_facecolor(self.C["surface"])
        self._image = None
        self._span = None
        if n == 1:
            self.axes = [self.fig.add_subplot(111)]
        else:
            gs = self.fig.add_gridspec(n, 1, height_ratios=heights, hspace=0.12)
            self.axes = []
            for i in range(n):
                a = self.fig.add_subplot(gs[i], sharex=self.axes[0] if (sharex and self.axes) else None)
                self.axes.append(a)
        self.fig.subplots_adjust(left=0.08, right=0.98, top=0.97, bottom=0.07)
        for a in self.axes:
            self._style(a)
        return self.axes

    def set_title(self, title, subtitle=""):
        self.title.setText(title)
        self.subtitle.setText(subtitle)

    def message(self, head, text="", title=""):
        """The empty state says what will appear, rather than that nothing has."""
        self.set_title(title, "")
        a = self._reset()[0]
        a.set_xticks([])
        a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
        a.text(0.04, 0.56, head, transform=a.transAxes, color=self.C["text_dim"], fontsize=13)
        a.text(0.04, 0.49, text, transform=a.transAxes, color=self.C["text_mute"], fontsize=10,
               va="top")
        self.canvas.draw_idle()

    # ------------------------------------------------------------- drawing
    def image(self, img, title="", subtitle="", overlays=(), lo=2, hi=98, cmap="gray"):
        a = self._reset()[0]
        self.set_title(title, subtitle)
        good = np.isfinite(img)
        v0, v1 = (np.percentile(img[good], [lo, hi]) if good.any() else (0, 1))
        a.imshow(img, origin="lower", cmap=cmap, vmin=v0, vmax=v1, aspect="auto",
                 interpolation="nearest")
        for ov in overlays:
            a.plot(ov["x"], ov["y"], ov.get("style", "-"), color=self.C.get(ov.get("color", "p_fit"), ov.get("color")),
                   lw=ov.get("lw", 1.0), ms=ov.get("ms", 3), alpha=ov.get("alpha", 0.9))
        a.set_xlim(0, img.shape[1] - 1)
        a.set_ylim(0, img.shape[0] - 1)
        self._image = img
        self.canvas.draw_idle()

    def spectra(self, panels, title="", subtitle="", xlabel="wavelength (um)", heights=None):
        """
        panels: a list of dicts, one per row, each with 'series' (x, y, colour
        key, label, lw, alpha, style), optional 'regions' ((x0, x1, colour key)),
        'vlines' ((x, colour key, label)), 'ylabel', 'ylim', 'hline'.
        """
        axes = self._reset(len(panels), heights=heights)
        self.set_title(title, subtitle)
        for a, p in zip(axes, panels):
            for s in p.get("series", []):
                col = self.C.get(s.get("color", "p_obs"), s.get("color"))
                a.plot(s["x"], s["y"], s.get("style", "-"), color=col, lw=s.get("lw", 0.9),
                       alpha=s.get("alpha", 1.0), ms=s.get("ms", 3), label=s.get("label"),
                       drawstyle=s.get("drawstyle", "default"))
            for (x0, x1, ck) in p.get("regions", []):
                a.axvspan(x0, x1, color=self.C.get(ck, ck), alpha=0.16, lw=0)
            for (x, ck, lab) in p.get("vlines", []):
                a.axvline(x, color=self.C.get(ck, ck), lw=0.8, ls=":")
                if lab:
                    a.text(x, 0.98, lab, transform=a.get_xaxis_transform(), rotation=90,
                           va="top", ha="right", fontsize=8, color=self.C["text_mute"])
            if "hline" in p:
                a.axhline(p["hline"], color=self.C["border_2"], lw=0.8)
            if p.get("ylim") is not None:
                a.set_ylim(*p["ylim"])
            elif p.get("robust", True):
                # Opaque telluric bands blow up to many times the spectrum;
                # the limits follow the bulk of the data instead.
                ys = np.concatenate([np.asarray(s["y"], float).ravel() for s in p.get("series", [])
                                     if np.size(s["y"])] or [np.array([0.0, 1.0])])
                ys = ys[np.isfinite(ys)]
                if ys.size > 10:
                    lo, hi = np.percentile(ys, [0.5, 99.5])
                    pad = 0.08 * (hi - lo if hi > lo else abs(hi) + 1e-30)
                    a.set_ylim(lo - pad, hi + pad)
            if p.get("ylabel"):
                a.set_ylabel(p["ylabel"], fontsize=9)
            if p.get("note"):
                a.text(0.005, 0.96, p["note"], transform=a.transAxes, va="top", ha="left",
                       fontsize=9, color=self.C["text_mute"])
            labs = [s.get("label") for s in p.get("series", []) if s.get("label")]
            if labs:
                leg = a.legend(loc="upper right", fontsize=8, frameon=False)
                for t in leg.get_texts():
                    t.set_color(self.C["text_dim"])
        for a in axes[:-1]:
            a.tick_params(labelbottom=False)
        axes[-1].set_xlabel(xlabel, fontsize=9)
        self.canvas.draw_idle()
        return axes

    def select_span(self, colour="accent", axis_index=0):
        """Let a drag across the plot pick a wavelength range."""
        a = self.axes[axis_index]
        # The selector adds an empty rectangle at x = 0 that autoscaling
        # would otherwise stretch the axes to include.
        xl, yl = a.get_xlim(), a.get_ylim()
        self._span = SpanSelector(a, lambda lo, hi: self.span.emit(float(min(lo, hi)), float(max(lo, hi))),
                                  "horizontal", useblit=True, minspan=1e-6, interactive=False,
                                  props=dict(facecolor=self.C["accent"], alpha=0.22))
        a.set_xlim(xl)
        a.set_ylim(yl)
        self._sync()
