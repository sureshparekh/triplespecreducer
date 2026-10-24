"""
TripleSpec Reducer — custom Qt widgets (from BRAIN).

Qt has no dual-handle range slider and no segmented control, so both are
drawn here. Everything takes its colours from the active palette, so the
light/dark switch reaches every pixel.
"""

import math
import os

from PySide6.QtCore import Qt, Signal, QRectF, QPointF
from PySide6.QtGui import QPainter, QColor, QPen, QBrush
from PySide6.QtWidgets import (QWidget, QSizePolicy, QHBoxLayout, QPushButton,
                               QButtonGroup, QLabel, QVBoxLayout, QFrame,
                               QCheckBox, QFormLayout)

from triplespec_reducer.gui import theme as TH


class RangeSlider(QWidget):
    """
    Two handles selecting a low and a high value.

    Continuous when given `lo`/`hi`; discrete when given `values`, in which
    case the handles snap to list entries — which is what template ages and
    metallicities need, since only catalogue values exist.
    """

    changed = Signal(float, float)

    def __init__(self, lo=0.0, hi=1.0, values=None, init=None,
                 fmt=lambda v: f"{v:g}", log=False, parent=None):
        super().__init__(parent)
        self.C = TH.LIGHT
        self.values = list(values) if values else None
        self.log = log and self.values is None
        self.lo_b, self.hi_b = ((0, len(self.values) - 1) if self.values
                                else (float(lo), float(hi)))
        self.fmt = fmt
        self.a, self.b = (self.lo_b, self.hi_b)
        if init:
            self.a, self.b = self._pos(init[0]), self._pos(init[1])
        self._drag = None
        self._hot = None
        self._key_handle = "a"
        self.setFixedHeight(56)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

    def set_palette(self, C):
        self.C = C
        self.update()

    def _pos(self, v):
        if self.values:
            return min(range(len(self.values)),
                       key=lambda i: abs(self.values[i] - v))
        return max(self.lo_b, min(self.hi_b, float(v)))

    def _val(self, p):
        return self.values[int(round(p))] if self.values else p

    def value(self):
        a, b = sorted((self.a, self.b))
        return self._val(a), self._val(b)

    def setValue(self, lo, hi, emit=True):
        self.a, self.b = self._pos(lo), self._pos(hi)
        self.update()
        if emit:
            self.changed.emit(*self.value())

    def setBounds(self, lo, hi, values=None, keep=True):
        """Re-range the slider, clamping the current selection into it."""
        cur = self.value()
        self.values = list(values) if values else None
        self.lo_b, self.hi_b = ((0, len(self.values) - 1) if self.values
                                else (float(lo), float(hi)))
        if keep:
            self.a = self._pos(max(min(cur[0], hi), lo))
            self.b = self._pos(max(min(cur[1], hi), lo))
        else:
            self.a, self.b = self.lo_b, self.hi_b
        self.update()
        self.changed.emit(*self.value())

    def _frac(self, p):
        if self.values:
            return (p - self.lo_b) / max(self.hi_b - self.lo_b, 1)
        if self.log:
            l0, l1 = math.log10(max(self.lo_b, 1e-30)), math.log10(max(self.hi_b, 1e-30))
            return (math.log10(max(p, 1e-30)) - l0) / max(l1 - l0, 1e-30)
        return (p - self.lo_b) / max(self.hi_b - self.lo_b, 1e-30)

    def _unfrac(self, f):
        f = max(0.0, min(1.0, f))
        if self.values:
            return round(self.lo_b + f * (self.hi_b - self.lo_b))
        if self.log:
            l0, l1 = math.log10(max(self.lo_b, 1e-30)), math.log10(max(self.hi_b, 1e-30))
            return 10 ** (l0 + f * (l1 - l0))
        return self.lo_b + f * (self.hi_b - self.lo_b)

    @property
    def _pad(self):
        return 13

    TRACK_Y = 44
    GRAB = 18

    def _x(self, p):
        w = max(self.width(), 40)
        return self._pad + self._frac(p) * (w - 2 * self._pad)

    def _near(self, x):
        return "a" if abs(x - self._x(self.a)) <= abs(x - self._x(self.b)) else "b"

    def _set(self, name, value):
        """
        Move one handle, and do not let it pass the other.

        Two handles that can cross still report a sorted pair, so nothing
        downstream notices, but the one under the cursor swaps identity
        mid-drag and the handle appears to stick while the other one runs
        away. Clamping each against the other is what makes the drag feel
        like it is attached to the thing being dragged.
        """
        if name == "a":
            self.a = min(value, self.b)
        else:
            self.b = max(value, self.a)

    def mousePressEvent(self, e):
        self.setFocus(Qt.MouseFocusReason)
        self._drag = self._near(e.position().x())
        self.mouseMoveEvent(e)

    def mouseMoveEvent(self, e):
        x = e.position().x()
        if self._drag:
            w = max(self.width(), 40)
            self._set(self._drag,
                      self._unfrac((x - self._pad) / max(w - 2 * self._pad, 1)))
            self.update()
            self.changed.emit(*self.value())
        else:
            near = abs(e.position().y() - self.TRACK_Y) < self.GRAB
            h = self._near(x) if near else None
            if h != self._hot:
                self._hot = h
                self.update()

    def mouseReleaseEvent(self, e):
        self._drag = None
        self.update()

    def leaveEvent(self, e):
        self._hot = None
        self.update()

    def keyPressEvent(self, e):
        """
        Arrow keys move the handle the keyboard is on.

        A custom-painted control that never takes focus cannot be reached
        without a mouse at all, which is not a styling shortcoming but a
        control that half the ways of using it cannot operate.
        """
        step = {Qt.Key_Left: -1, Qt.Key_Down: -1,
                Qt.Key_Right: 1, Qt.Key_Up: 1}.get(e.key())
        if step is None:
            if e.key() == Qt.Key_Tab:
                self._key_handle = "b" if self._key_handle == "a" else "a"
                self.update()
                return
            return super().keyPressEvent(e)
        if e.modifiers() & Qt.ShiftModifier:
            step *= 10
        name = self._key_handle
        current = getattr(self, name)
        if self.values:
            self._set(name, max(self.lo_b, min(self.hi_b, current + step)))
        else:
            span = self.hi_b - self.lo_b
            self._set(name, max(self.lo_b, min(
                self.hi_b, self._unfrac(self._frac(current) + step * 0.01))))
        self.update()
        self.changed.emit(*self.value())

    def paintEvent(self, _):
        """
        A track, the span between the handles, and two handles.

        The accent used to paint the whole track, both handles and every tick,
        which is most of a panel's width given to the one colour that is
        supposed to mean something. It now marks the span that is selected and
        the handle being touched, and the rest of the control is drawn in the
        neutral ramp.

        The values sit under their own handle rather than pinned to the ends
        of the widget, so reading one does not mean working out which end it
        belongs to, and they are clamped so neither runs off an edge.

        Catalogue ticks are drawn above the track rather than below it. Below,
        they are the lowest thing the control paints and they collided with
        the label of whatever the layout put next; above, they sit between the
        numbers and the track they belong to and nothing can reach them.
        """
        C = self.C
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        w = max(self.width(), 40)
        y = self.TRACK_Y
        xa, xb = sorted((self._x(self.a), self._x(self.b)))

        p.setPen(Qt.NoPen)
        p.setBrush(QBrush(QColor(C["surface_3"])))
        p.drawRoundedRect(QRectF(self._pad, y - 1.5, w - 2 * self._pad, 3),
                          1.5, 1.5)
        p.setBrush(QBrush(QColor(C["accent"])))
        p.drawRoundedRect(QRectF(xa, y - 1.5, max(xb - xa, 3), 3), 1.5, 1.5)

        if self.values and len(self.values) <= 40:
            for i in range(len(self.values)):
                x = self._x(i)
                inside = xa - 0.5 <= x <= xb + 0.5
                p.setPen(QPen(QColor(C["text_mute"] if inside
                                     else C["border_2"]), 1))
                p.drawLine(QPointF(x, y - 10), QPointF(x, y - 7))

        focused = self.hasFocus()
        for pos, name in ((self.a, "a"), (self.b, "b")):
            x = round(self._x(pos))
            hot = self._hot == name or self._drag == name
            ring = (focused and self._key_handle == name)
            if hot or ring:
                p.setPen(Qt.NoPen)
                p.setBrush(QBrush(QColor(C["accent_soft"])))
                p.drawEllipse(QPointF(x, y), 11, 11)
            p.setPen(QPen(QColor(C["accent"] if (hot or ring)
                                 else C["text_mute"]), 1.6))
            p.setBrush(QBrush(QColor(C["surface_2"] if not hot
                                     else C["accent"])))
            p.drawEllipse(QPointF(x, y), 6, 6)

        lo, hi = self.value()
        p.setFont(TH.mono(TH.SIZE["label"]))
        for value, pos in ((lo, min(self.a, self.b)),
                           (hi, max(self.a, self.b))):
            text = self.fmt(value)
            width = p.fontMetrics().horizontalAdvance(text) + 8
            x = self._x(pos) - width / 2
            x = max(0.0, min(float(w) - width, x))
            p.setPen(QColor(C["text"]))
            p.drawText(QRectF(x, 6, width, 18),
                       Qt.AlignHCenter | Qt.AlignVCenter, text)
        p.end()

class Segmented(QWidget):
    """
    One row of mutually exclusive buttons.

    An option can be marked unavailable, which disables its button, strikes
    its label through and puts the reason on it as a tooltip. A choice the
    run would silently ignore is worse than no choice at all: the fit
    finishes, the numbers look ordinary, and nothing says the setting did
    not apply. The strikethrough carries the same message as the greying,
    so the state does not depend on seeing colour.
    """

    changed = Signal(str)

    def __init__(self, options, init=0, parent=None):
        super().__init__(parent)
        self.reasons = {}
        # The row is a track the chosen option sits in as a pill; the track
        # paints itself, so it takes its fill from the stylesheet.
        self.setObjectName("SegTrack")
        self.setAttribute(Qt.WA_StyledBackground, True)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(3, 3, 3, 3)
        lay.setSpacing(2)
        self.group = QButtonGroup(self)
        self.group.setExclusive(True)
        self.buttons = []
        for i, name in enumerate(options):
            b = QPushButton(name)
            b.setCheckable(True)
            # Focus from the keyboard only: a click would otherwise leave a
            # focus ring on the option just left, so two would look chosen.
            b.setFocusPolicy(Qt.TabFocus)
            b.setObjectName("SegL" if i == 0 else
                            "SegR" if i == len(options) - 1 else "Seg")
            if i not in (0, len(options) - 1):
                b.setProperty("class", "Seg")
            b.setChecked(i == init)
            self.group.addButton(b, i)
            lay.addWidget(b, 1)
            self.buttons.append(b)
        self.options = list(options)
        self.group.idClicked.connect(
            lambda i: self.changed.emit(self.options[i]))

    def setProperty(self, name, value):
        """
        Carry `error` down to the buttons that are actually painted.

        A Segmented is a layout holding push buttons and paints nothing
        itself, so a stylesheet rule keyed on a property of the container
        reaches no pixels. An unanswered control that has been told it is in
        error and shows nothing is worse than one that was never asked: the
        caller believes the user has been warned. The flag is pushed onto the
        children, which are the things with a border.
        """
        ok = super().setProperty(name, value)
        if name == "error":
            for b in getattr(self, "buttons", []):
                b.setProperty("error", bool(value))
                b.style().unpolish(b)
                b.style().polish(b)
        return ok

    def value(self):
        """
        The chosen option, or None when nothing has been chosen.

        QButtonGroup answers -1 when no button is checked, and indexing the
        options with that returns the LAST one. A control built with no
        initial selection would therefore report a confident answer nobody
        gave, which is this project's characteristic way of being wrong. The
        caller that asks a question with no default is entitled to hear that
        it has not been answered.
        """
        index = self.group.checkedId()
        if index < 0:
            return None
        return self.options[index]

    def setValue(self, name):
        if name in self.options:
            b = self.buttons[self.options.index(name)]
            if b.isEnabled():
                b.setChecked(True)

    def set_available(self, name, available, reason=""):
        """Enable or disable one option, carrying the reason on the button."""
        if name not in self.options:
            return
        b = self.buttons[self.options.index(name)]
        self.reasons[name] = "" if available else reason
        f = b.font()
        f.setStrikeOut(not available)
        b.setFont(f)
        b.setToolTip(reason)
        if not available and b.isChecked():
            self.group.setExclusive(False)
            b.setChecked(False)
            self.group.setExclusive(True)
        b.setEnabled(available)
        b.style().unpolish(b)
        b.style().polish(b)

    def available(self):
        """The options a user can still choose, in panel order."""
        return [n for n, b in zip(self.options, self.buttons) if b.isEnabled()]

    def value(self):
        i = self.group.checkedId()
        return self.options[i] if i >= 0 else None

    def clear(self):
        """
        Put the control back to unanswered, the state `init=None` also means.

        setValue cannot express it: a name that is not an option is ignored,
        so handing it None leaves whatever was chosen before still chosen.
        That matters when a saved configuration is loaded over a page a user
        has already touched -- without this, a file that declares no answer
        inherits the last one and the page shows a setting the file does not
        contain. Exclusivity is lifted for the unchecking because a Qt button
        group in exclusive mode refuses to let its last checked button go.
        """
        self.group.setExclusive(False)
        for b in self.buttons:
            b.setChecked(False)
        self.group.setExclusive(True)

    def restyle(self):
        for b in self.buttons:
            b.style().unpolish(b)
            b.style().polish(b)


class Disclosure(QWidget):
    """
    A heading that folds a group of controls away, shut by default.

    Every setting the reducer accepts has to be reachable, but a panel that shows
    all of them at once is a wall of spin boxes in which the handful that
    change a result are no easier to find than the solver tolerances. The
    ones that decide the physics stay in the open; the rest live behind one
    of these, which states how many settings it holds so their existence is
    visible even when their controls are not.
    """

    toggled = Signal(bool)

    def __init__(self, title, count=None, parent=None):
        super().__init__(parent)
        self.title = title
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(8)
        self.head = QPushButton()
        self.head.setObjectName("Disclosure")
        self.head.setCheckable(True)
        self.head.setText(f"{title}"
                          + (f"   ({count})" if count is not None else ""))
        self.head.toggled.connect(self._show)
        v.addWidget(self.head)
        self.body = QWidget()
        self.bodyv = QVBoxLayout(self.body)
        self.bodyv.setContentsMargins(0, 2, 0, 6)
        self.bodyv.setSpacing(9)
        self.body.setVisible(False)
        v.addWidget(self.body)

    def _show(self, on):
        self.body.setVisible(on)
        self.restyle()
        self.toggled.emit(on)

    def add(self, w):
        self.bodyv.addWidget(w)
        return w

    def restyle(self):
        self.head.setIcon(TH.icon("chevron_down" if self.head.isChecked()
                                  else "chevron_right", self._arrow, 11))

    def set_palette(self, C):
        self._arrow = C["text_mute"]
        self.restyle()

    _arrow = "#8b95a3"


TAG_TEXT = {"required": "Required", "optional": "Optional",
            "unused": "Not needed"}


def make_tag():
    """A small label beside a heading saying whether the panel must be
    answered. Hidden until set_tag gives it a kind."""
    tag = QLabel("")
    tag.setObjectName("Tag")
    tag.setVisible(False)
    return tag


def set_tag(tag, kind, text=None):
    """
    Mark a panel Required, Optional or Not needed.

    Required: the run is refused, or silently wrong, without an answer here.
    Optional: the defaults are a sound fit; change them only with a reason.
    Not needed: the current choices never read this panel.
    """
    if not kind:
        tag.setVisible(False)
        return
    tag.setText(text or TAG_TEXT[kind])
    tag.setProperty("kind", kind)
    tag.style().unpolish(tag)
    tag.style().polish(tag)
    tag.setVisible(True)


class Plate(QWidget):
    """
    A panel with a heading, and no bevel anywhere on it.

    A QGroupBox draws a sunken frame with its title notched into the top
    edge. That frame is the most recognisable mark of an interface nobody
    has looked at, and the title inside it repeats a heading the sidebar has
    already said, so the same four words arrived three times on one screen.
    What is left here is a surface one value above the page, a heading in
    the interface grotesque, and the form beneath it. Depth comes from value,
    as it does everywhere else in this application.

    An optional stage carries a checkbox in place of its heading rather than
    beside it, because the heading IS the thing being switched on, and the
    body is disabled rather than hidden when it is off: a stage that
    disappears when it is not running cannot be read to find out what it
    would have done.
    """

    toggled = Signal(bool)

    def __init__(self, title, checkable=False, parent=None):
        super().__init__(parent)
        self.setObjectName("Plate")
        self.setAttribute(Qt.WA_StyledBackground, True)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(TH.SPACE["wide"], TH.SPACE["base"],
                                 TH.SPACE["wide"], TH.SPACE["wide"])
        outer.setSpacing(TH.SPACE["base"])

        self.switch = None
        top = QHBoxLayout()
        top.setContentsMargins(0, 0, 0, 0)
        top.setSpacing(TH.SPACE["snug"])
        if checkable:
            self.switch = QCheckBox(title)
            self.switch.setObjectName("PlateSwitch")
            self.switch.setFont(TH.heading())
            self.switch.toggled.connect(self._switched)
            top.addWidget(self.switch)
        else:
            head = QLabel(title)
            head.setObjectName("SectionTitle")
            head.setFont(TH.heading())
            top.addWidget(head)
        top.addStretch(1)
        self.tag = make_tag()
        top.addWidget(self.tag)
        outer.addLayout(top)

        self.body = QWidget()
        self.form = FieldForm(self.body)
        self.form.setContentsMargins(0, 0, 0, 0)
        self.form.setSpacing(TH.SPACE["base"])
        self.form.setLabelAlignment(Qt.AlignLeft)
        outer.addWidget(self.body)
        if checkable:
            self.body.setEnabled(False)

    def _switched(self, on):
        self.body.setEnabled(bool(on))
        self.toggled.emit(bool(on))

    def isChecked(self):
        """True when the stage is on, and true always when it cannot be off."""
        return True if self.switch is None else self.switch.isChecked()

    def setChecked(self, on):
        if self.switch is not None:
            self.switch.setChecked(bool(on))

    def set_tag(self, kind, text=None):
        set_tag(self.tag, kind, text)

    def set_error(self, on):
        """Mark the panel as holding something the run would refuse."""
        self.setProperty("error", "true" if on else "false")
        self.style().unpolish(self)
        self.style().polish(self)


def vdivider(height=18):
    """A short vertical rule, for separating items inside a row."""
    f = QFrame()
    f.setObjectName("VDivider")
    f.setFrameShape(QFrame.VLine)
    f.setFixedWidth(1)
    f.setFixedHeight(height)
    return f


def divider():
    f = QFrame()
    f.setObjectName("Divider")
    f.setFrameShape(QFrame.HLine)
    return f


class Section(QWidget):
    """
    A titled group under a hairline.

    The title is set in true small capitals where the face carries them,
    which is a heading rather than the tracked-out uppercase micro-label that
    every generated interface stacks above every section. Where the face is
    missing, the fallback capitalises instead.
    """

    def __init__(self, title, hint=None, parent=None):
        super().__init__(parent)
        self.title = title
        self.v = QVBoxLayout(self)
        self.v.setContentsMargins(0, 0, 0, 0)
        self.v.setSpacing(8)
        self.v.addWidget(divider())
        self.v.addSpacing(2)
        head = QHBoxLayout()
        head.setContentsMargins(0, 0, 0, 0)
        t = QLabel(title)
        t.setObjectName("SectionTitle")
        t.setFont(TH.smallcaps(TH.SIZE["heading"]))
        head.addWidget(t)
        head.addStretch(1)
        if hint:
            h = QLabel(hint)
            h.setObjectName("Hint")
            head.addWidget(h)
        self.tag = make_tag()
        head.addWidget(self.tag)
        self.v.addLayout(head)
        self._title_label = t
        self._chrome = [self.v.itemAt(0).widget(), t, self.tag]

    def set_title(self, title):
        self.title = title
        self._title_label.setText(title)

    def set_tag(self, kind, text=None):
        set_tag(self.tag, kind, text)

    def hide_title(self):
        """For a section that is the only thing in a step of the same name."""
        for w in self._chrome:
            if w is not None:
                w.setVisible(False)

    def add(self, w):
        self.v.addWidget(w)
        return w

    def addLayout(self, l):
        self.v.addLayout(l)
        return l


WRAPPING = ("Caption", "CaptionWarn", "Hint", "Note")


def label(text, kind="FieldLabel"):
    """
    A label, wrapping only where wrapping is what it is for.

    Every label used to word-wrap. A word-wrapped QLabel reports its height
    as a function of its width, and a QScrollArea that is resizable asks its
    content for a height before the content knows how wide it will be. The
    answer comes back short, the layout takes the difference out of whichever
    items can give it, and a one word field label is squeezed to no height at
    all and paints on top of the control above it. The prose labels still
    wrap, because a paragraph has to; `Age` and `Metallicity` never did.
    """
    l = QLabel(text)
    l.setObjectName(kind)
    l.setWordWrap(kind in WRAPPING)
    if kind not in WRAPPING:
        l.setMinimumHeight(l.fontMetrics().height())
    return l


LABEL_COL = 17 * TH.UNIT
"""
The width of every panel's label column, in pixels.

A QFormLayout sizes its label column to its own longest label, so each
panel found a different one and the fields started at a different x in
every panel down the page. Fixing the column makes one vertical edge that
every control on the page is aligned to, which is the thing a column of
settings is read down.
"""


class FieldForm(QFormLayout):
    """
    A form whose field column starts in the same place in every panel, and
    whose prose starts there too.

    Two alignment faults came from the plain QFormLayout. The field column
    moved from panel to panel, because each form measured only its own
    labels. And a caption added as a single widget spanned both columns, so
    the paragraph explaining a control began at the panel's left edge while
    the control it explains began a hundred-odd pixels to the right, with
    the text running the full width of the plate and touching both edges.

    Here a label is given the fixed column, and a wrapping caption is
    indented to where the fields begin. Labels sit to the left of that edge;
    everything else -- field, button, table, paragraph -- begins on it.
    """

    def addRow(self, first, second=None):
        if second is None:
            if isinstance(first, QWidget) and first.objectName() in WRAPPING:
                gap = self.horizontalSpacing()
                holder = QWidget()
                inner = QHBoxLayout(holder)
                inner.setContentsMargins(
                    LABEL_COL + (gap if gap > 0 else TH.SPACE["base"]), 0, 0, 0)
                inner.setSpacing(0)
                inner.addWidget(first)
                return super().addRow(holder)
            return super().addRow(first)
        if isinstance(first, str):
            first = label(first)
        if isinstance(first, QWidget):
            first.setFixedWidth(LABEL_COL)
        return super().addRow(first, second)


def row(*widgets, stretch=None):
    """
    Widgets side by side, with `None` standing for a spring.

    A horizontal layout gives space it has no instruction for to whichever
    children will take it, and a QLabel will. Two call sites asked for a
    trailing spring by naming a `stretch` index one past the last widget,
    which matches nothing, so the slack went to the labels instead and
    opened arbitrary gaps in the middle of a row of steppers. Passing None
    puts the spring where it is meant to go and is visible in the call.
    """
    w = QWidget()
    h = QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    h.setSpacing(8)
    position = 0
    for x in widgets:
        if x is None:
            h.addStretch(1)
            continue
        if isinstance(x, str):
            x = label(x)
        h.addWidget(x, (stretch or {}).get(position, 0))
        position += 1
    return w


class SeriesKey(QWidget):
    """
    Which lines are on the plot, and in what colour.

    A key belongs in the chrome, not floating inside the data area on a
    filled panel of its own. Each entry is a short rule in the series colour
    with its name beside it, in the order the series are drawn.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.items = []
        self.setFont(TH.sans(TH.SIZE["label"]))
        self.setFixedHeight(20)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.C = TH.DARK

    def set_palette(self, C):
        self.C = C
        self.update()

    def set_items(self, items):
        self.items = list(items)
        self.setFixedWidth(self._measure())
        self.updateGeometry()
        self.update()

    def _measure(self):
        if not self.items:
            return 0
        fm = self.fontMetrics()
        return sum(14 + 5 + fm.horizontalAdvance(n) + 14
                   for n, _ in self.items) - 14

    def paintEvent(self, _):
        if not self.items:
            return
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setFont(self.font())
        fm = p.fontMetrics()
        y = self.height() / 2.0
        x = 0.0
        for name, colour in self.items:
            p.setPen(Qt.NoPen)
            p.setBrush(QBrush(QColor(colour)))
            p.drawRect(QRectF(x, y - 1.0, 14.0, 2.0))
            x += 19.0
            p.setPen(QColor(self.C["text_dim"]))
            p.drawText(QRectF(x, 0, fm.horizontalAdvance(name) + 2,
                              self.height()),
                       Qt.AlignLeft | Qt.AlignVCenter, name)
            x += fm.horizontalAdvance(name) + 14.0
        p.end()


class Strip(QWidget):
    """
    A row that drops its least important parts before anything overlaps.

    Qt will happily let two labels in one layout paint over each other rather
    than admit the row does not fit, which is how a strip that reads cleanly
    at 1500 pixels turns into overlapping text at 1000. Each optional child
    is registered with the width below which it goes, widest-first, and the
    row hides them in that order on every resize.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._optional = []

    def optional(self, widget, below):
        self._optional.append((widget, below))
        return widget

    def resizeEvent(self, e):
        super().resizeEvent(e)
        w = self.width()
        for widget, below in self._optional:
            widget.setVisible(w >= below)


def failure_headline(message):
    """
    The sentence a failed run's message should open with.

    A run that raised comes back as the exception and a Python traceback
    from the worker thread that runs each step. The
    traceback is for a bug report; the person at the window needs the
    exception's own words, which are its last line.
    """
    text = str(message or "").strip()
    if "Traceback (most recent call last)" in text:
        lines = [l.strip() for l in text.splitlines() if l.strip()]
        last = lines[-1] if lines else text
        head = text.split("Traceback (most recent call last)")[0].strip()
        return head.splitlines()[0] if head else last
    return text.split("\n\n")[0].strip() or "The run failed."


def failure_dialog(parent, message, folder=None, title="TripleSpec Reducer"):
    """
    A failed run, said in one readable sentence, with the whole message one
    click away under "Show Details" and kept in a log file.
    """
    import datetime
    from PySide6.QtWidgets import QMessageBox
    text = str(message or "")
    saved = ""
    if folder:
        try:
            os.makedirs(folder, exist_ok=True)
            path = os.path.join(os.path.abspath(folder), "tsr_error.log")
            with open(path, "a") as fh:
                fh.write(f"\n=== {datetime.datetime.now():%Y-%m-%d %H:%M:%S}"
                         f" ===\n{text}\n")
            saved = f"\n\nThe full message is in {path}."
        except OSError:
            pass
    box = QMessageBox(parent)
    box.setIcon(QMessageBox.Critical)
    box.setWindowTitle(title)
    box.setText("The run stopped with an error.")
    box.setInformativeText(failure_headline(text)[:600] + saved)
    box.setDetailedText(text)
    box.exec()
