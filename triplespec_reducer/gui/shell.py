"""
TripleSpec Reducer: the workspace shell (from BRAIN).

Every workspace (a spectrum, a datacube) is a step bar over two columns:

    step bar: 1 Spectrum > 2 Corrections > 3 Masks > 4 Model > 5 Fit
    the current step's controls | the stage

The bar is the workflow in the order it is done. Each step shows an icon,
its name and one line of state, and whether it is ready, so the state of a
whole fit can be read across the top without opening anything. Only the
current step's controls are on screen, in a column as wide as the window
allows; the stage is the thing being worked on (a spectrum and its fit, a
cube and its maps) and stays put whichever step is open.

A step is numbered because the steps are a sequence: corrections are applied
before the fit, and the model is chosen before it runs. A step that is off
(an optional stage nobody switched on) is drawn quietly rather than hidden,
because a stage that disappears cannot be read to find out what it would
have done.
"""

from PySide6.QtCore import Qt, Signal, QRectF, QSize
from PySide6.QtGui import QPainter, QColor, QPen, QFontMetrics, QFont, QIcon
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel,
                               QPushButton, QScrollArea, QStackedWidget,
                               QSplitter, QSizePolicy, QFrame)

from triplespec_reducer.gui import theme as TH

READY, TODO, ATTENTION, OFF = "ready", "todo", "attention", "off"


class StepItem(QWidget):
    """One step in the bar: an icon in a ring, its name and a line of state."""

    clicked = Signal()

    HEIGHT = 60

    def __init__(self, number, title, icon=None, parent=None):
        super().__init__(parent)
        self.number = number
        self.title = title
        self.icon_name = icon
        self.summary = ""
        self.state = TODO
        self.current = False
        self.last = False
        self._hover = False
        self.C = TH.DARK
        self.setFixedHeight(self.HEIGHT)
        self.setMinimumWidth(150)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setAccessibleName(title)

    def set_palette(self, colours):
        self.C = colours
        self.update()

    def set_state(self, state, summary=None):
        self.state = state
        if summary is not None:
            self.summary = summary
            self.setToolTip(f"{self.title}: {summary}")
        self.update()

    def set_current(self, on):
        self.current = bool(on)
        self.update()

    def enterEvent(self, e):
        self._hover = True
        self.update()

    def leaveEvent(self, e):
        self._hover = False
        self.update()

    def mousePressEvent(self, e):
        if e.button() == Qt.LeftButton:
            self.clicked.emit()

    def keyPressEvent(self, e):
        if e.key() in (Qt.Key_Return, Qt.Key_Enter, Qt.Key_Space):
            self.clicked.emit()
        else:
            super().keyPressEvent(e)

    def sizeHint(self):
        return QSize(200, self.HEIGHT)

    def paintEvent(self, e):
        C = self.C
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        r = QRectF(self.rect())
        card = r.adjusted(4, 8, -4, -8)
        if self.current or self._hover:
            # The current step is a raised neutral card, as the current page
            # is in the website's sidebar; hovering shows a fainter one.
            p.setPen(QPen(QColor(C["border_2"] if self.current
                                 else C["border"]), 1) if self.current
                     else Qt.NoPen)
            p.setBrush(QColor(C["surface_3"] if self.current
                              else C["surface_2"]))
            p.drawRoundedRect(card.adjusted(0.5, 0.5, -0.5, -0.5), 12, 12)
        if self.hasFocus():
            p.setPen(QPen(QColor(C["accent"]), 1.4))
            p.setBrush(Qt.NoBrush)
            p.drawRoundedRect(card.adjusted(0.5, 0.5, -0.5, -0.5), 12, 12)

        badge = QRectF(card.left() + 10, card.center().y() - 15, 30, 30)
        if self.current:
            # the strongest neutral, like the primary action
            ink = C["ink"] if C["name"] == "dark" else C["void"]
            ring, fill = C["text"], C["text"]
        elif self.state == ATTENTION:
            ring, fill, ink = C["crit"], None, C["crit"]
        elif self.state == READY:
            ring, fill, ink = C["good"], None, C["good"]
        else:
            ring, fill, ink = C["border_2"], None, C["text_mute"]
        p.setPen(QPen(QColor(ring), 1.4))
        p.setBrush(QColor(fill) if fill else Qt.NoBrush)
        p.drawEllipse(badge)

        glyph = ("check" if self.state == READY and not self.current
                 else self.icon_name)
        if self.state == ATTENTION and not self.current:
            p.setPen(QColor(ink))
            f = TH.sans(13, 700)
            f.setPixelSize(14)
            p.setFont(f)
            p.drawText(badge, Qt.AlignCenter, "!")
        elif glyph:
            TH.icon(glyph, ink, 16).paint(
                p, badge.adjusted(7, 7, -7, -7).toRect())
        else:
            p.setPen(QColor(ink))
            f = TH.mono(TH.SIZE["label"])
            f.setPixelSize(TH.SIZE["label"] + 1)
            p.setFont(f)
            p.drawText(badge, Qt.AlignCenter, str(self.number))

        x0 = badge.right() + 10
        width = card.right() - x0 - 8
        # The theme's helpers size in points and the stylesheet in pixels;
        # pixel sizes are set here so painted text matches the rest.
        title_font = TH.sans(TH.SIZE["body"] + 1)
        title_font.setPixelSize(TH.SIZE["body"] + 2)
        title_font.setWeight(QFont.Weight.DemiBold if self.current
                             else QFont.Weight.Medium)
        p.setFont(title_font)
        title_colour = (C["text"] if self.current else
                        C["text_mute"] if self.state == OFF else C["text_dim"])
        p.setPen(QColor(title_colour))
        cy = card.center().y()
        p.drawText(QRectF(x0, cy - 17, width, 18),
                   Qt.AlignLeft | Qt.AlignVCenter,
                   QFontMetrics(title_font).elidedText(
                       f"{self.number}  {self.title}", Qt.ElideRight,
                       int(width)))

        if self.summary:
            small = TH.sans(TH.SIZE["caption"])
            small.setPixelSize(TH.SIZE["caption"])
            p.setFont(small)
            p.setPen(QColor(C["crit"] if self.state == ATTENTION
                            else C["text_mute"]))
            text = QFontMetrics(small).elidedText(self.summary, Qt.ElideRight,
                                                  int(width))
            p.drawText(QRectF(x0, cy + 1, width, 16),
                       Qt.AlignLeft | Qt.AlignVCenter, text)
        p.end()


class StepBar(QWidget):
    """The workflow, left to right, one StepItem per step."""

    selected = Signal(int)

    def __init__(self, titles, icons=None, parent=None):
        super().__init__(parent)
        self.setObjectName("StepBar")
        self.setAttribute(Qt.WA_StyledBackground, True)
        h = QHBoxLayout(self)
        h.setContentsMargins(2 * TH.UNIT, 0, 2 * TH.UNIT, 0)
        h.setSpacing(0)
        self.items = []
        self.arrows = []
        icons = list(icons or [None] * len(titles))
        for i, title in enumerate(titles):
            if i:
                arrow = QLabel()
                arrow.setFixedWidth(18)
                arrow.setAlignment(Qt.AlignCenter)
                h.addWidget(arrow)
                self.arrows.append(arrow)
            item = StepItem(i + 1, title, icons[i])
            item.clicked.connect(lambda i=i: self.selected.emit(i))
            h.addWidget(item, 1)
            self.items.append(item)

    def set_current(self, index):
        for i, item in enumerate(self.items):
            item.set_current(i == index)

    def set_palette(self, colours):
        for item in self.items:
            item.set_palette(colours)
        for arrow in self.arrows:
            arrow.setPixmap(TH.icon("chevron_right", colours["border_2"],
                                    14).pixmap(14, 14))


# The old name, for anything that still asks for a rail.
StepRail = StepBar


class StepPage(QWidget):
    """
    One step's controls: a heading, a sentence saying what the step is for,
    the controls themselves, and a way on to the next step.
    """

    back = Signal()
    forward = Signal()

    def __init__(self, title, description, next_title=None, parent=None):
        super().__init__(parent)
        self.setObjectName("StepPage")
        self.setAttribute(Qt.WA_StyledBackground, True)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        head = QWidget()
        hv = QVBoxLayout(head)
        hv.setContentsMargins(4 * TH.UNIT, 3 * TH.UNIT, 4 * TH.UNIT,
                              2 * TH.UNIT)
        hv.setSpacing(TH.UNIT)
        t = QLabel(title)
        t.setObjectName("StepTitle")
        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.setSpacing(TH.UNIT)
        top.addWidget(t)
        from triplespec_reducer.gui.widgets import make_tag
        self.tag = make_tag()
        top.addWidget(self.tag, 0, Qt.AlignVCenter)
        top.addStretch(1)
        hv.addLayout(top)
        if description:
            d = QLabel(description)
            d.setObjectName("StepLede")
            d.setWordWrap(True)
            hv.addWidget(d)
        outer.addWidget(head)

        self.scroll = QScrollArea()
        self.scroll.setObjectName("StepScroll")
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.NoFrame)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        body = QWidget()
        body.setMinimumWidth(60 * TH.UNIT)
        self.body = QVBoxLayout(body)
        self.body.setContentsMargins(4 * TH.UNIT, TH.UNIT, 4 * TH.UNIT,
                                     3 * TH.UNIT)
        self.body.setSpacing(2 * TH.UNIT)
        self.scroll.setWidget(body)
        outer.addWidget(self.scroll, 1)

        foot = QWidget()
        foot.setObjectName("StepFoot")
        foot.setAttribute(Qt.WA_StyledBackground, True)
        fh = QHBoxLayout(foot)
        fh.setContentsMargins(4 * TH.UNIT, 12,
                              4 * TH.UNIT, 12)
        fh.setSpacing(TH.UNIT)
        self.btn_back = QPushButton("Back")
        self.btn_back.setObjectName("Ghost")
        self.btn_back.clicked.connect(self.back)
        fh.addWidget(self.btn_back)
        fh.addStretch(1)
        self.foot_extra = QHBoxLayout()
        self.foot_extra.setSpacing(TH.UNIT)
        fh.addLayout(self.foot_extra)
        # "&" in a button label is Qt's shortcut marker; doubled, it is an "&".
        nxt = (next_title or "").replace("&", "&&")
        self.btn_next = QPushButton(f"Next: {nxt}" if next_title else "")
        self.btn_next.setObjectName("Next")
        self.btn_next.setLayoutDirection(Qt.RightToLeft)
        self.btn_next.clicked.connect(self.forward)
        self.btn_next.setVisible(bool(next_title))
        fh.addWidget(self.btn_next)
        outer.addWidget(foot)

    def set_tag(self, kind, text=None):
        """Required, Optional or Not needed, beside the step's title."""
        from triplespec_reducer.gui.widgets import set_tag
        set_tag(self.tag, kind, text)

    def add(self, widget):
        self.body.addWidget(widget)
        return widget

    def finish(self):
        self.body.addStretch(1)


class Workspace(QWidget):
    """
    Step bar over step and stage, wired together.

    `steps` is a list of (title, description) pairs, or (title, description,
    icon) triples naming a glyph in icons/. Pages are created here and filled
    by the caller through `page(i)`.
    """

    step_changed = Signal(int)

    def __init__(self, steps, stage, parent=None):
        super().__init__(parent)
        steps = [tuple(s) + (None,) * (3 - len(s)) for s in steps]
        self.titles = [t for t, _, _ in steps]
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(0)

        self.rail = StepBar(self.titles, [i for _, _, i in steps])
        self.rail.selected.connect(self.go)
        v.addWidget(self.rail)

        self.pages = QStackedWidget()
        self.pages.setObjectName("Pages")
        self._pages = []
        for i, (title, description, _) in enumerate(steps):
            nxt = self.titles[i + 1] if i + 1 < len(steps) else None
            page = StepPage(title, description, nxt)
            page.btn_back.setVisible(i > 0)
            page.back.connect(lambda i=i: self.go(i - 1))
            page.forward.connect(lambda i=i: self.go(i + 1))
            self.pages.addWidget(page)
            self._pages.append(page)

        split = QSplitter(Qt.Horizontal)
        split.setChildrenCollapsible(False)
        split.setHandleWidth(1)
        split.addWidget(self.pages)
        split.addWidget(stage)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([660, 840])
        self.splitter = split
        v.addWidget(split, 1)
        self.go(0)

    def page(self, index):
        return self._pages[index]

    def go(self, index):
        index = max(0, min(int(index), len(self._pages) - 1))
        self.pages.setCurrentIndex(index)
        self.rail.set_current(index)
        self.step_changed.emit(index)

    def current(self):
        return self.pages.currentIndex()

    def set_state(self, index, state, summary=None):
        self.rail.items[index].set_state(state, summary)

    def set_palette(self, colours):
        self.rail.set_palette(colours)
        for page in self._pages:
            page.btn_next.setIcon(TH.icon("chevron_right", colours["ink"], 14))
