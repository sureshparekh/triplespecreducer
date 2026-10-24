"""
TripleSpec Reducer — palette, typography and stylesheet (BRAIN's design system).

One place for every colour, font and piece of chrome in the application, so
the two themes stay in step and nothing is styled ad hoc.

The interface wears the documentation website's clothes, so the two read as
one thing. The neutrals are warm greys, never blue-grey: an off-white ground
and white wells in light mode, near-black in dark mode with panels one step
off it. Colour is spent sparingly. The primary action is the strongest
neutral there is -- near-black in light mode, near-white in dark -- rather
than a blue slab, and cobalt, the logo's blue, is kept for small signals: a
checked box, a focus ring, progress, and BRAIN's own fit drawn over the data.
The colours that carry data sit apart from it: the stellar component is
amber, gas green, nebular violet, so the fit is the only blue line on a plot
of the observed spectrum.

Type is the website's three faces: Hanken Grotesk for every control and
sentence, Source Serif 4 for the wordmark and the title of each step, and
JetBrains Mono for every number, with tabular figures so a column of values
does not shift between 1 and 8.

Geometry follows one rule. Anything pressed is a capsule (buttons, chips,
the segmented choices, whose chosen option is lifted out of its track as a
pill); anything typed into is a 10 px rounded field; panels are rounder,
14 px, than anything inside them. Depth comes from value and a hairline rim,
not from shadow: in dark mode a panel's top edge is a shade lighter than its
sides, where light would catch a sheet of glass.

Every control has a visible focus ring and a visible disabled state, and a
control holding something the run would refuse carries `error` as a Qt
property, which paints its border in the critical colour, so a validation
message sits next to the field it is about instead of only in a dialog.

`palette(C)` gives Qt's own palette the same colours. The stylesheet does not
reach every widget Qt builds itself -- the file chooser's views, message boxes
-- and those used the desktop's palette, so on a dark desktop the light theme
put dark text on a black file dialog.
"""

import math
import os
import tempfile

from PySide6.QtSvg import QSvgRenderer

from PySide6.QtGui import (QFont, QFontDatabase, QFontInfo, QIcon, QPixmap,
                           QPainter, QColor, QPen, QPainterPath)
from PySide6.QtCore import Qt, QRectF, QPointF, QByteArray

UNIT = 8

SPACE = dict(hair=2, tight=4, snug=8, base=12, wide=16, room=24, air=32)

SIZE = dict(wordmark=19, heading=13, body=12, small=12, label=11, caption=11,
            micro=10)

RADIUS = 10

# Panels are rounder than the fields inside them, so a field never looks like
# a panel: 10 px on a field, 14 px on a panel, and a full capsule on anything
# pressed (half a row's height).
RADIUS_PANEL = 14

ROW = 32
ROW_TALL = 36


DARK = dict(
    name="dark",
    void="#0A0A09",
    bg="#0D0D0C",
    surface="#141413",
    surface_2="#1B1B1A",
    surface_3="#262624",
    border="#242422",
    border_2="#3A3A37",
    rim="#34332F",
    text="#F0EFEC",
    text_dim="#CFCEC4",
    text_mute="#A3A19A",
    accent="#7896FF",
    accent_hi="#9CB2FF",
    accent_dim="#2C4290",
    accent_soft="#172347",
    ink="#0D0D0C",
    good="#4CCF8A",
    warn="#F4C97E",
    crit="#FF7B6E",
    p_grid="#1D1D1B",
    p_obs="#E6E4DD",
    p_faint="#33322F",
    p_fit="#7896FF",
    p_star="#F0A04B",
    p_gas="#3FB984",
    p_neb="#B49AF0",
    p_resid="#8C8A83",
    p_mask="#FF7B6E",
    p_weight="#3FB984",
    p_scrim="#0D0D0C",
    # the Results spectrum: data, total fit, SSP, emission lines
    p_spec_obs="#B8B5AB",
    p_spec_total="#F0A04B",
    p_spec_ssp="#7896FF",
    p_spec_lines="#45C4B0",
)

LIGHT = dict(
    name="light",
    void="#FFFFFF",
    bg="#F7F7F5",
    surface="#FCFCFB",
    surface_2="#F2F1EE",
    surface_3="#E9E8E3",
    border="#E5E4DE",
    border_2="#D0CEC6",
    rim="#E5E4DE",
    text="#141413",
    text_dim="#46453F",
    text_mute="#63615B",
    accent="#2A4FC4",
    accent_hi="#1F43B0",
    accent_dim="#B9C8F0",
    accent_soft="#EDF1FC",
    ink="#FFFFFF",
    good="#0D6B3C",
    warn="#8A5A00",
    crit="#B42318",
    p_grid="#EFEEE9",
    p_obs="#2B2A27",
    p_faint="#D9D7CF",
    p_fit="#2A4FC4",
    p_star="#D9822B",
    p_gas="#1C7C54",
    p_neb="#7A4FB5",
    p_resid="#77756E",
    p_mask="#B42318",
    p_weight="#1C7C54",
    p_scrim="#FFFFFF",
    p_spec_obs="#5A5850",
    p_spec_total="#D9822B",
    p_spec_ssp="#2A4FC4",
    p_spec_lines="#0F8A7A",
)


from triplespec_reducer._paths import FONT_DIR as _FONT_DIR, ICON_DIR as _ICON_DIR

WORDMARK_STACK = ["Source Serif 4", "DejaVu Serif"]

_WORDMARK = WORDMARK_STACK
_DISPLAY = ["Source Serif 4", "DejaVu Serif"]
_SANS = ["Hanken Grotesk", "DejaVu Sans"]
_CAPS = ["Hanken Grotesk", "DejaVu Sans"]
_MONO = ["JetBrains Mono", "DejaVu Sans Mono", "Liberation Mono",
         "Courier New"]
# Prose is set in the interface face; the serif is kept for titles.
_SERIF = _SANS

_cache = {}
_loaded = False


def load_fonts():
    """
    Register the faces that ship with the package, before anything asks for one.

    A font that is not installed falls back without saying so, and an
    interface whose type silently changes on somebody else's machine is one
    nobody designed. The faces in `fonts/` are loaded from there, so the
    stacks below resolve the same way everywhere. Calling this more than
    once is harmless and does nothing after the first time.
    """
    global _loaded
    if _loaded:
        return
    _loaded = True
    if not os.path.isdir(_FONT_DIR):
        return
    for name in sorted(os.listdir(_FONT_DIR)):
        if name.lower().endswith((".ttf", ".otf")):
            QFontDatabase.addApplicationFont(os.path.join(_FONT_DIR, name))
    _cache.clear()


def _pick(prefs, key):
    load_fonts()
    if key in _cache:
        return _cache[key]
    have = set(QFontDatabase.families())
    for p in prefs:
        if p in have:
            _cache[key] = p
            return p
    _cache[key] = QFontInfo(QFont(prefs[-1])).family()
    return _cache[key]


def _weight(w):
    """Accept a plain CSS-style number as well as a QFont.Weight."""
    if isinstance(w, QFont.Weight):
        return w
    return min(QFont.Weight, key=lambda k: abs(int(k) - int(w)))


def _tabular(font):
    """
    Ask for tabular figures, on the Qt versions that can.

    A proportional 1 is narrower than a 0 in most grotesques, so a column of
    numbers that changes by one digit visibly shifts. The feature is a
    request: a face that does not carry it is left as it is rather than
    swapped for one that does.
    """
    try:
        font.setFeature(QFont.Tag("tnum"), 1)
    except (AttributeError, TypeError, ValueError):
        pass
    return font


def display(size=SIZE["wordmark"], weight=500, spacing=-0.4):
    f = QFont(_pick(_DISPLAY, "display"), size)
    f.setWeight(_weight(weight))
    if spacing:
        f.setLetterSpacing(QFont.AbsoluteSpacing, spacing)
    return f



def wordmark(size=SIZE["wordmark"], weight=600, spacing=-0.3):
    """
    The logotype, and only the logotype.

    A wordmark is display type and is allowed a face the interface does not
    otherwise use; what it is not allowed to do is leak into a heading or a
    button. It is named once, in WORDMARK, so supplying a licensed copy of
    that family is the whole of the change needed to use it.
    """
    f = QFont(_pick(_WORDMARK, "wordmark"), size)
    f.setWeight(_weight(weight))
    if spacing:
        f.setLetterSpacing(QFont.AbsoluteSpacing, spacing)
    return f


def wordmark_family():
    return _pick(_WORDMARK, "wordmark")


def serif(size=SIZE["body"], weight=QFont.Normal, italic=False):
    f = QFont(_pick(_SERIF, "serif"), size)
    f.setWeight(_weight(weight))
    f.setItalic(italic)
    return f


def sans(size=SIZE["body"], weight=QFont.Normal, spacing=0.0):
    f = QFont(_pick(_SANS, "sans"), size)
    f.setWeight(_weight(weight))
    if spacing:
        f.setLetterSpacing(QFont.AbsoluteSpacing, spacing)
    return _tabular(f)


def heading(size=SIZE["heading"], weight=600, spacing=0.0):
    """
    A section heading, in the interface grotesque and in sentence case.

    Headings used to be set in small capitals, and to fall back to
    `title.upper()` wherever that face was missing. A tracked-out uppercase
    micro-label is the most recognisable signature of a generated interface,
    and a fallback that produces one on every machine without one specific
    font is not a fallback, it is the default. The heading is now the words
    as they were written, in the same face as the controls beneath it, and
    it earns its rank from weight and colour rather than from shouting.
    """
    f = QFont(_pick(_CAPS, "caps"), size)
    f.setWeight(_weight(weight))
    if spacing:
        f.setLetterSpacing(QFont.AbsoluteSpacing, spacing)
    return _tabular(f)


smallcaps = heading


def mono(size=SIZE["label"], weight=QFont.Normal):
    f = QFont(_pick(_MONO, "mono"), size)
    f.setWeight(_weight(weight))
    return _tabular(f)


def font_report():
    return (f"display={_pick(_DISPLAY,'display')}  "
            f"sans={_pick(_SANS,'sans')}  caps={_pick(_CAPS,'caps')}  "
            f"mono={_pick(_MONO,'mono')}")


_icon_cache = {}


def icon(kind, colour, size=18):
    """
    One icon set, on one grid, recoloured at paint time.

    The files in `icons/` are Lucide and carry `stroke="currentColor"` with
    no colour of their own, so the palette decides what an icon looks like
    rather than the asset does. The same file is text_dim while it is idle
    and the accent once it is active, without a second copy of it existing,
    and a palette change reaches every glyph without touching any of them.
    That property is most of the difference between a set that belongs to an
    interface and a set that was pasted into one.

    Rendered at three times the requested size and scaled down: a 1.75 unit
    stroke on a 24 unit grid falls between pixels at 16 px, and a renderer
    asked for the final size directly rounds it heavier on one edge than the
    other, which is visible as a glyph that leans.

    A missing file returns an empty icon rather than raising. A button whose
    glyph did not load still has its label and its tooltip; a window that
    refuses to open has neither.
    """
    key = (kind, str(colour), int(size))
    if key in _icon_cache:
        return _icon_cache[key]
    px = int(size * 3)
    pm = QPixmap(px, px)
    pm.fill(Qt.transparent)
    path = os.path.join(_ICON_DIR, f"{kind}.svg")
    if os.path.exists(path):
        with open(path) as handle:
            body = handle.read()
        body = body.replace('stroke="currentColor"', f'stroke="{colour}"')
        body = body.replace('stroke-width="2"', 'stroke-width="1.75"')
        renderer = QSvgRenderer(QByteArray(body.encode("utf-8")))
        painter = QPainter(pm)
        painter.setRenderHint(QPainter.Antialiasing)
        renderer.render(painter, QRectF(0, 0, px, px))
        painter.end()
    # One pixmap per screen density, so the glyph stays sharp on a high-DPI
    # display (and in the documentation's 2x screenshots) instead of being
    # drawn at 1x and stretched.
    out = QIcon()
    for ratio in (1, 2, 3):
        scaled = pm.scaled(int(size * ratio), int(size * ratio),
                           Qt.KeepAspectRatio, Qt.SmoothTransformation)
        scaled.setDevicePixelRatio(ratio)
        out.addPixmap(scaled)
    _icon_cache[key] = out
    return out


_STYLE_ICON_DIR = os.path.join(
    tempfile.gettempdir(), f"tsr-style-icons-{os.getuid()}")


def style_icon(kind, colour):
    """
    One glyph with its colour resolved, written where a stylesheet can read it.

    `image: url(...)` hands the file to Qt's SVG renderer with nothing to
    inherit from, so the `stroke="currentColor"` every Lucide file carries
    resolves to nothing and the glyph is drawn empty. Every check mark and
    every drop-down chevron in the application was therefore invisible: the
    control was styled, the state was set, and the one mark that says which
    state it is in was not drawn. The copies are written once per colour,
    outside the repository, and the originals are left as they are so
    `icon()` keeps recolouring them at paint time.

    A glyph that is not in the set, or a cache that cannot be written,
    returns an empty string: a stylesheet rule with no image loses its mark
    and keeps its control, which is what the old behaviour already was.
    """
    source = os.path.join(_ICON_DIR, f"{kind}.svg")
    if not os.path.exists(source):
        return ""
    target = os.path.join(_STYLE_ICON_DIR,
                          f"{kind}-{str(colour).lstrip('#')}.svg")
    if os.path.exists(target):
        return target
    try:
        os.makedirs(_STYLE_ICON_DIR, exist_ok=True)
        with open(source) as handle:
            body = handle.read()
        with open(target, "w") as handle:
            handle.write(body.replace('stroke="currentColor"',
                                      f'stroke="{colour}"'))
    except OSError:
        return ""
    return target


def icon_names():
    """Every glyph the set ships, for a caller that wants to check one exists."""
    if not os.path.isdir(_ICON_DIR):
        return ()
    return tuple(sorted(n[:-4] for n in os.listdir(_ICON_DIR)
                        if n.endswith(".svg")))



def logo(path, height):
    """
    Load a logo and trim the transparent margin before scaling.

    The source artwork is a 2000x2000 square with a lot of empty space
    around the mark; scaled naively it lands in the header looking half the
    size it should be. Cropping to the visible pixels first fixes that.
    """
    pm = QPixmap(path)
    if pm.isNull():
        return pm
    img = pm.toImage()
    if img.hasAlphaChannel():
        rect = img.createAlphaMask().rect()
        try:
            from PySide6.QtGui import QRegion, QBitmap
            region = QRegion(QBitmap.fromImage(img.createAlphaMask()))
            if not region.isEmpty():
                rect = region.boundingRect()
        except Exception:
            pass
        if rect.width() > 4 and rect.height() > 4:
            pm = pm.copy(rect)
    # Scaled for the screen's density, so the logo is sharp on a high-DPI
    # display rather than drawn at `height` pixels and stretched.
    from PySide6.QtGui import QGuiApplication
    app = QGuiApplication.instance()
    ratio = app.devicePixelRatio() if app is not None else 1.0
    out = pm.scaledToHeight(int(round(height * ratio)), Qt.SmoothTransformation)
    out.setDevicePixelRatio(ratio)
    return out



def apply_matplotlib(C):
    """Point matplotlib at the same faces the interface uses."""
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [_pick(_SANS, "sans"), "DejaVu Sans"],
        "font.monospace": [_pick(_MONO, "mono"), "DejaVu Sans Mono"],
        "axes.titlesize": SIZE["label"],
        "axes.labelsize": SIZE["micro"],
        "axes.titleweight": "medium",
        "mathtext.fontset": "dejavusans",
        "figure.facecolor": C["surface"],
        "savefig.facecolor": C["surface"],
        "path.simplify": True,
        "path.simplify_threshold": 1.0,
        "agg.path.chunksize": 20000,
        "lines.solid_joinstyle": "round",
        "lines.solid_capstyle": "round",
    })


def palette(C):
    """
    Qt's own palette in this theme's colours.

    Set on the application alongside the stylesheet. Widgets Qt builds for
    itself -- the file chooser's list and tree, its path bar, message boxes --
    take their colours from here, and without it they took the desktop's.
    """
    from PySide6.QtGui import QPalette
    pal = QPalette()
    roles = {
        QPalette.Window: C["bg"], QPalette.WindowText: C["text"],
        QPalette.Base: C["void"], QPalette.AlternateBase: C["surface_2"],
        QPalette.Text: C["text"], QPalette.Button: C["surface_2"],
        QPalette.ButtonText: C["text"], QPalette.BrightText: C["text"],
        QPalette.ToolTipBase: C["surface_3"], QPalette.ToolTipText: C["text"],
        QPalette.Highlight: C["accent"], QPalette.HighlightedText: C["ink"],
        QPalette.Link: C["accent"], QPalette.PlaceholderText: C["text_mute"],
        QPalette.Light: C["surface_3"], QPalette.Midlight: C["surface_2"],
        QPalette.Mid: C["border_2"], QPalette.Dark: C["border"],
        QPalette.Shadow: C["void"],
    }
    for role, colour in roles.items():
        pal.setColor(role, QColor(colour))
    for role in (QPalette.WindowText, QPalette.Text, QPalette.ButtonText):
        pal.setColor(QPalette.Disabled, role, QColor(C["text_mute"]))
    return pal


def stylesheet(C):
    """
    Qt stylesheet for a palette dict.

    Sizes come from SIZE, spacings from SPACE and corners from RADIUS rather
    than from numbers written here, so the scale is one decision made once.
    The four blocks at the end are the ones that make the interface
    answerable: focus, so the keyboard is usable; disabled, so a control that
    cannot be used says so; error, so a field the run would refuse is marked
    where it is; and hover, so a thing that can be clicked looks like one.

    An accent fill carries `ink`: white on the light theme's cobalt (6.3:1),
    near-black on the dark theme's lighter blue (7:1).
    """
    f_sans = _pick(_SANS, "sans")
    f_mono = _pick(_MONO, "mono")
    f_caps = _pick(_CAPS, "caps")
    f_serif = _pick(_SERIF, "serif")
    f_word = _pick(_WORDMARK, "wordmark")
    body, label, micro = SIZE["body"], SIZE["label"], SIZE["micro"]
    r = RADIUS
    rp = RADIUS_PANEL
    return f"""
/* A plain container paints nothing and shows whatever it sits on. The
   ground is painted once, by the window and by the named surfaces below.
   Giving every QWidget the page background instead filled each layout
   holder inside a panel with the page's own colour, so a black rectangle
   appeared behind the form in every lighter plate -- one per stage, the
   width of the fields and as tall as the rows. The holders are structure,
   not surfaces, and structure has no colour of its own. */
QWidget {{
    background: transparent;
    color: {C['text']};
    font-family: "{f_sans}";
    font-size: {body}px;
}}
QMainWindow > QWidget, QScrollArea, QScrollArea > QWidget > QWidget {{
    background: transparent;
}}
QSplitter::handle {{ background: {C['border']}; }}
QMainWindow, QDialog {{ background: {C['bg']}; }}

#Panel, #Header, #StatusBar, #Well {{ background: {C['surface']}; }}
#Header  {{ border-bottom: 1px solid {C['border']}; }}
#StatusBar {{ border-top: 1px solid {C['border']}; }}
#Panel   {{ border-right: 1px solid {C['border']}; }}
#Stage   {{ background: {C['bg']}; }}
#Card {{
    background: {C['surface']};
    border: 1px solid {C['border']};
    border-radius: {r}px;
}}
/* Each stage is a ruled box. With the ground at true black a panel can no
   longer be told apart by its fill -- surface sits one step off black and
   reads as the same colour -- so the edge is what separates one stage from
   the next. The fill stays for the light theme, where it still carries. */
/* ---- workspace shell: rail | step | stage ---- */
#StepBar  {{ background: {C['surface']}; border-bottom: 1px solid {C['border']}; }}
#Pages, #StepPage, #StepScroll, #StepScroll > QWidget > QWidget {{
    background: {C['bg']};
}}
#StepFoot {{ background: {C['bg']}; border-top: 1px solid {C['border']}; }}
QLabel#StepTitle {{
    color: {C['text']};
    font-family: "{f_sans}";
    font-size: {SIZE['heading'] + 7}px;
    font-weight: 500;
}}
QLabel#StepLede {{
    color: {C['text_dim']};
    font-family: "{f_serif}";
    font-size: {SIZE['body'] + 1}px;
}}

#Plate {{
    background: {C['surface']};
    border: 1px solid {C['border']};
    border-radius: {r}px;
}}
#Plate[error="true"] {{ border: 1px solid {C['crit']}; }}
QCheckBox#PlateSwitch {{
    color: {C['text']};
    font-family: "{f_caps}";
    font-size: {SIZE['heading']}px;
    font-weight: 600;
    spacing: 9px;
}}
QCheckBox#PlateSwitch::indicator {{
    width: 15px; height: 15px;
    border: 1px solid {C['border_2']};
    border-radius: 2px;
    background: {C['void']};
}}
QCheckBox#PlateSwitch::indicator:checked {{
    background: {C['accent_soft']}; border-color: {C['accent']};
    image: url({style_icon('check', C['accent'])});
}}
QCheckBox#PlateSwitch::indicator:hover {{ border-color: {C['accent']}; }}
QCheckBox#PlateSwitch:!checked {{ color: {C['text_dim']}; }}

/* ---- headings and labels ---- */
#SectionTitle {{
    color: {C['text']};
    font-family: "{f_caps}";
    font-size: {SIZE['heading']}px;
    font-weight: 500;
    letter-spacing: 0.4px;
}}
#Wordmark {{
    color: {C['text']};
    font-family: "{f_word}";
    font-size: {SIZE['wordmark']}px;
    font-weight: 700;
    letter-spacing: -0.4px;
}}
#Caption   {{ color: {C['text_mute']}; font-size: {label}px;
             font-family: "{f_serif}"; }}
#CaptionWarn {{ color: {C['warn']}; font-size: {label}px;
               font-family: "{f_serif}"; }}
#FieldLabel {{ color: {C['text_dim']}; font-size: {label}px; }}
#Value     {{ color: {C['text']}; font-family: "{f_mono}";
              font-size: {label}px; }}
#Readout   {{ color: {C['text']}; font-family: "{f_mono}";
              font-size: {SIZE['heading']}px; }}
#Hint      {{ color: {C['text_mute']}; font-size: {label}px; }}
#Error     {{ color: {C['crit']}; font-size: {label}px; }}
#Good      {{ color: {C['good']}; font-size: {label}px; }}
#Accent    {{ color: {C['accent']}; }}

/* ---- primary action ---- */
QPushButton#Primary {{
    background: {C['accent']};
    color: {C['ink']};
    border: 1px solid {C['accent']};
    border-radius: {r}px;
    padding: 0 18px;
    min-height: {ROW_TALL}px;
    font-size: {SIZE['heading']}px;
    font-weight: 700;
}}
QPushButton#Primary:hover    {{ background: {C['accent_hi']};
                                border-color: {C['accent_hi']}; }}
QPushButton#Primary:pressed  {{ background: {C['accent_dim']};
                                border-color: {C['accent_dim']}; }}
QPushButton#Primary:disabled {{ background: {C['surface_2']};
                                border-color: {C['border']};
                                color: {C['text_mute']}; }}

/* ---- secondary / ghost ---- */
QPushButton {{
    background: {C['surface_2']};
    color: {C['text_dim']};
    border: 1px solid {C['border_2']};
    border-radius: {r}px;
    padding: 0 12px;
    min-height: {ROW}px;
    font-size: {label}px;
}}
QPushButton:hover   {{ background: {C['surface_3']}; color: {C['text']};
                       border-color: {C['text_mute']}; }}
QPushButton:pressed {{ background: {C['surface_2']}; }}
QPushButton#Tool {{
    background: transparent;
    border: 1px solid transparent;
    border-radius: {r}px;
    padding: 0;
    color: {C['text_dim']};
}}
QPushButton#Tool:hover {{
    background: {C['surface_2']}; border-color: {C['border_2']};
}}
QPushButton#Tool:checked {{
    background: {C['accent_soft']}; border-color: {C['accent']};
}}
QPushButton#Chip {{
    background: transparent;
    border: 1px solid {C['border_2']};
    border-radius: {r}px;
    padding: 0 12px; min-height: {ROW}px;
    color: {C['text_dim']};
    font-family: "{f_mono}";
}}
QPushButton#Chip:hover   {{ color: {C['text']}; border-color: {C['text_mute']}; }}
QPushButton#Chip:checked {{ background: {C['accent_soft']};
                            border-color: {C['accent']}; color: {C['text']}; }}

QPushButton#Step {{
    background: {C['surface_2']};
    border: 1px solid {C['border']};
    border-radius: {r}px;
    padding: 0;
}}
QPushButton#Step:hover  {{ background: {C['surface_3']};
                           border-color: {C['border_2']}; }}
QPushButton#Step:pressed {{ background: {C['surface']}; }}

QPushButton#Ghost {{
    background: transparent; border: 1px solid {C['border']};
    border-radius: {r}px;
    padding: 0 10px; min-height: {ROW}px;
    color: {C['text_dim']};
}}
QPushButton#Ghost:hover {{ color: {C['text']};
                           border-color: {C['border_2']};
                           background: {C['surface_2']}; }}
QPushButton#Ghost:pressed {{ background: {C['surface']}; }}
QPushButton#Ghost:disabled {{ color: {C['text_mute']};
                              border-color: {C['border']}; }}
QPushButton:focus, QPushButton#Ghost:focus, QPushButton#Tool:focus {{
    outline: none;
    border: 1px solid {C['accent']};
}}
QPushButton#Primary:focus {{ border: 1px solid {C['text']}; }}
QPushButton#Tool:hover   {{ color: {C['text']}; }}
QPushButton#Tool:checked {{ color: {C['accent']}; }}
QPushButton[error="true"], QPushButton#Seg[error="true"],
QPushButton#SegL[error="true"], QPushButton#SegR[error="true"] {{
    border: 1px solid {C['crit']};
    color: {C['crit']};
}}
QPushButton:disabled {{
    background: {C['surface']}; color: {C['text_mute']};
    border-color: {C['border']};
}}
QPushButton#Primary:disabled {{
    background: {C['accent_dim']}; color: {C['text_mute']};
}}

/* ---- inputs ---- */
QLineEdit, QSpinBox, QDoubleSpinBox {{
    background: {C['void']};
    border: 1px solid {C['border_2']};
    border-radius: {r}px;
    padding: 0 8px;
    min-height: {ROW}px;
    selection-background-color: {C['accent']};
    selection-color: {C['ink']};
    font-family: "{f_mono}";
    font-size: {label}px;
}}
QLineEdit:hover, QSpinBox:hover, QDoubleSpinBox:hover {{
    border-color: {C['text_mute']};
}}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
    border-color: {C['accent']};
}}
QLineEdit#Search {{ font-family: "{f_sans}"; }}

QComboBox {{
    background: {C['void']};
    border: 1px solid {C['border_2']};
    border-radius: {r}px;
    padding: 0 8px;
    min-height: {ROW}px;
    font-size: {label}px;
}}
QComboBox:hover {{ border-color: {C['text_mute']}; }}
QComboBox::drop-down {{ border: none; width: 22px; }}
QComboBox::down-arrow {{
    image: url({style_icon('chevron_down', C['text_dim'])});
    width: 12px; height: 12px;
}}
QComboBox::down-arrow:hover, QComboBox::down-arrow:on {{
    image: url({style_icon('chevron_down', C['text'])});
}}
QComboBox::down-arrow:disabled {{
    image: url({style_icon('chevron_down', C['text_mute'])});
}}
QComboBox QAbstractItemView::item {{
    min-height: {ROW}px; padding: 0 8px; border: none;
}}
QComboBox QAbstractItemView::item:hover {{
    background: {C['surface_3']}; color: {C['text']};
}}
QComboBox QAbstractItemView::item:selected {{
    background: {C['accent_soft']}; color: {C['text']};
    border-left: 2px solid {C['accent']};
}}
QComboBox QAbstractItemView {{
    background: {C['surface_2']};
    border: 1px solid {C['border_2']};
    selection-background-color: {C['accent']};
    selection-color: {C['ink']};
    outline: none;
    padding: 2px;
}}

/* ---- lists: a well, not a card ---- */
QListWidget {{
    background: {C['void']};
    border: 1px solid {C['border']};
    border-radius: 0;
    font-family: "{f_mono}";
    font-size: {label}px;
    outline: none;
    padding: 0;
}}
QListWidget::item {{ padding: 4px 8px; border-radius: 0;
                     border-left: 2px solid transparent; }}
QListWidget::item:selected {{ background: {C['accent_soft']};
                              color: {C['text']};
                              border-left: 2px solid {C['accent']}; }}
QListWidget::item:hover:!selected {{ background: {C['surface_2']};
                                     color: {C['text']}; }}

/* ---- checkboxes ---- */
QCheckBox {{ color: {C['text_dim']}; font-size: {label}px; spacing: 8px; }}
QCheckBox::indicator {{
    width: 15px; height: 15px;
    border: 1px solid {C['border_2']};
    border-radius: 2px;
    background: {C['void']};
}}
QCheckBox::indicator:checked {{
    background: {C['accent_soft']}; border-color: {C['accent']};
    image: url({style_icon('check', C['accent'])});
}}
QCheckBox::indicator:hover {{ border-color: {C['accent']}; }}

/* ---- scrollbars: theirs by default, ours here ---- */
QScrollArea {{ border: none; background: {C['surface']}; }}
QScrollBar:vertical {{
    background: transparent; width: 9px; margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {C['border_2']}; border-radius: 0; min-height: 30px;
}}
QScrollBar::handle:vertical:hover {{ background: {C['text_mute']}; }}
QScrollBar:horizontal {{ background: transparent; height: 9px; margin: 0; }}
QScrollBar::handle:horizontal {{
    background: {C['border_2']}; border-radius: 0; min-width: 30px;
}}
QScrollBar::handle:horizontal:hover {{ background: {C['text_mute']}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}

/* ---- segmented control ---- */
QPushButton#Seg, QPushButton#SegL, QPushButton#SegR {{
    background: {C['surface_2']};
    border: 1px solid {C['border_2']};
    border-radius: 0;
    padding: 0 6px;
    min-height: {ROW}px;
    font-size: {label}px;
    color: {C['text_dim']};
}}
QPushButton#Seg:checked, QPushButton#SegL:checked,
QPushButton#SegR:checked {{
    background: {C['accent']}; color: {C['ink']};
    border-color: {C['accent']}; font-weight: 700;
}}
QPushButton#Seg:hover:!checked, QPushButton#SegL:hover:!checked,
QPushButton#SegR:hover:!checked {{ background: {C['surface_3']};
                                   color: {C['text']}; }}
QPushButton#Seg:disabled, QPushButton#SegL:disabled,
QPushButton#SegR:disabled {{
    background: {C['surface']}; color: {C['text_mute']};
    border-color: {C['border']};
}}
QPushButton#SegL {{ border-top-left-radius: {r}px;
                    border-bottom-left-radius: {r}px; }}
QPushButton#SegR {{ border-top-right-radius: {r}px;
                    border-bottom-right-radius: {r}px; }}

/* ---- disclosure: settings that exist but are not the first question ---- */
QPushButton#Disclosure {{
    background: transparent;
    border: none;
    border-top: 1px solid {C['border']};
    border-radius: 0;
    padding: 10px 0 8px 0;
    /* a heading line; at 0 a form squeezes the text out from under the chevron */
    min-height: {SIZE['heading'] + 8}px;
    text-align: left;
    color: {C['text_dim']};
    font-family: "{f_caps}";
    font-size: {SIZE['heading']}px;
    font-weight: 500;
    letter-spacing: 0.4px;
}}
QPushButton#Disclosure:hover {{ color: {C['text']}; }}
QPushButton#Disclosure:checked {{ color: {C['text']}; }}
QPushButton#Disclosure:focus {{ border: none;
                                border-top: 1px solid {C['accent']}; }}

QToolTip {{
    background: {C['surface_3']}; color: {C['text']};
    border: 1px solid {C['border_2']}; padding: 4px 7px;
    border-radius: {r}px; font-size: {label}px;
}}
#Divider {{ background: {C['border']}; max-height: 1px; min-height: 1px;
            border: none; }}
#VDivider {{ background: {C['border']}; border: none; }}

/* ---- optional stages, tables, log and progress ---- */
QGroupBox {{
    background: {C['surface']};
    border: 1px solid {C['border']};
    border-radius: {r}px;
    margin-top: 14px;
    padding: 14px 12px 12px 12px;
    font-size: {label}px;
    font-weight: 600;
    color: {C['text']};
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 10px;
    padding: 0 5px;
    color: {C['text']};
}}
QGroupBox::indicator {{
    width: 13px; height: 13px;
    border: 1px solid {C['border_2']};
    border-radius: 2px;
    background: {C['void']};
}}
QGroupBox::indicator:checked {{
    background: {C['accent']}; border-color: {C['accent']};
}}
QGroupBox:!enabled {{ color: {C['text_mute']}; }}

QTableWidget {{
    background: {C['void']};
    alternate-background-color: {C['surface']};
    border: 1px solid {C['border']};
    border-radius: 0;
    gridline-color: {C['border']};
    font-family: "{f_mono}";
    font-size: {label}px;
    selection-background-color: {C['accent']};
    selection-color: {C['ink']};
}}
QHeaderView::section {{
    background: {C['surface_2']};
    color: {C['text_dim']};
    border: none;
    border-right: 1px solid {C['border']};
    border-bottom: 1px solid {C['border']};
    padding: 4px 7px;
    font-family: "{f_sans}";
    font-size: {micro}px;
    font-weight: 600;
}}
QTableCornerButton::section {{ background: {C['surface_2']}; border: none; }}

QTextEdit, QPlainTextEdit {{
    background: {C['void']};
    color: {C['text_dim']};
    border: 1px solid {C['border']};
    border-radius: 0;
    font-family: "{f_mono}";
    font-size: {micro}px;
    padding: 6px;
    selection-background-color: {C['accent']};
    selection-color: {C['ink']};
}}

QProgressBar {{
    background: {C['surface_2']};
    border: 1px solid {C['border']};
    border-radius: 0;
    height: 6px;
    text-align: center;
    font-size: {micro}px;
    color: {C['text_dim']};
}}
QProgressBar::chunk {{ background: {C['accent']}; border-radius: 0; }}


/* ---- the rest of the widget set, so none of it arrives from the platform ---- */

QRadioButton {{ spacing: 7px; font-size: {label}px; color: {C['text_dim']}; }}
QRadioButton::indicator {{
    width: 15px; height: 15px;
    border: 1px solid {C['border_2']};
    border-radius: 8px;
    background: {C['void']};
}}
QRadioButton::indicator:hover {{ border-color: {C['accent']}; }}
QRadioButton::indicator:checked {{
    border: 5px solid {C['accent']};
    background: {C['void']};
}}
QRadioButton:disabled {{ color: {C['text_mute']}; }}

QCheckBox {{ spacing: 7px; font-size: {label}px; color: {C['text_dim']}; }}
QCheckBox:disabled {{ color: {C['text_mute']}; }}
QCheckBox::indicator:disabled {{
    border-color: {C['border']}; background: {C['surface']};
}}

QMenuBar {{ background: {C['surface']}; color: {C['text_dim']};
            border-bottom: 1px solid {C['border']}; }}
QMenuBar::item {{ padding: 5px 10px; background: transparent;
                  border-radius: {r}px; }}
QMenuBar::item:selected {{ background: {C['surface_3']}; color: {C['text']}; }}
QMenu {{
    background: {C['surface_2']};
    border: 1px solid {C['border_2']};
    border-radius: {r}px;
    padding: 4px;
    font-size: {label}px;
}}
QMenu::item {{ padding: 6px 22px 6px 12px; border-radius: {r}px;
               color: {C['text_dim']}; }}
QMenu::item:selected {{ background: {C['accent_soft']}; color: {C['text']}; }}
QMenu::item:disabled {{ color: {C['text_mute']}; }}
QMenu::separator {{ height: 1px; background: {C['border']}; margin: 4px 6px; }}
QMenu::indicator {{ width: 14px; height: 14px; margin-left: 6px; }}

QSplitter::handle {{ background: {C['border']}; }}
QSplitter::handle:horizontal {{ width: 1px; }}
QSplitter::handle:vertical {{ height: 1px; }}
QSplitter::handle:hover {{ background: {C['accent_dim']}; }}
QSplitter::handle:pressed {{ background: {C['accent']}; }}

QTabWidget::pane {{
    border: 1px solid {C['border']}; border-radius: {r}px;
    background: {C['surface']}; top: -1px;
}}
QTabBar::tab {{
    background: transparent; color: {C['text_mute']};
    padding: 7px 14px; margin-right: 2px;
    border: 1px solid transparent;
    border-top-left-radius: {r}px; border-top-right-radius: {r}px;
    font-size: {label}px;
}}
QTabBar::tab:hover {{ color: {C['text']}; background: {C['surface_2']}; }}
QTabBar::tab:selected {{
    color: {C['text']}; background: {C['surface']};
    border-color: {C['border']}; border-bottom-color: {C['surface']};
}}
QTabBar::tab:focus {{ border-color: {C['accent']}; }}

QScrollBar::handle:vertical:pressed,
QScrollBar::handle:horizontal:pressed {{ background: {C['accent']}; }}

QHeaderView::section:hover {{ background: {C['surface_3']};
                              color: {C['text']}; }}

QToolTip {{ opacity: 245; }}

QTextEdit, QPlainTextEdit {{
    background: {C['void']};
    border: 1px solid {C['border_2']};
    border-radius: {r}px;
    selection-background-color: {C['accent']};
    selection-color: {C['ink']};
    font-family: "{f_mono}";
    font-size: {label}px;
}}

QSpinBox::up-button, QDoubleSpinBox::up-button,
QSpinBox::down-button, QDoubleSpinBox::down-button {{
    background: {C['surface_2']}; border: none; width: 16px;
}}
QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover,
QSpinBox::down-button:hover, QDoubleSpinBox::down-button:hover {{
    background: {C['surface_3']};
}}
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {{
    image: url({style_icon('chevron_up', C['text_dim'])});
    width: 9px; height: 9px;
}}
QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {{
    image: url({style_icon('chevron_down', C['text_dim'])});
    width: 9px; height: 9px;
}}

QDialogButtonBox QPushButton {{ min-width: 84px; }}

QPushButton#NavOn {{
    background: {C['surface_3']};
    border: 1px solid transparent;
    border-left: 2px solid {C['accent']};
    color: {C['text']};
    text-align: left;
    padding: 0 12px;
}}
QPushButton#Nav {{
    background: transparent;
    border: 1px solid transparent;
    border-left: 2px solid transparent;
    color: {C['text_dim']};
    text-align: left;
    padding: 0 12px;
}}
QPushButton#Nav:hover {{ background: {C['surface_2']}; color: {C['text']}; }}

/* ---- focus: the keyboard has to be able to find where it is ---- */
QPushButton:focus, QComboBox:focus, QCheckBox:focus, QGroupBox:focus,
QTableWidget:focus, QListWidget:focus, QTextEdit:focus {{
    border: 1px solid {C['accent']};
    outline: none;
}}
QPushButton#Primary:focus {{ border: 1px solid {C['accent_hi']}; }}
QPushButton#Nav:focus {{ border-left: 2px solid {C['accent']}; }}

/* ---- disabled: a control that cannot be used should not look ready ---- */
QLineEdit:disabled, QComboBox:disabled, QSpinBox:disabled,
QDoubleSpinBox:disabled {{
    background: {C['surface']};
    color: {C['text_mute']};
    border-color: {C['border']};
}}
QLabel:disabled, QCheckBox:disabled {{ color: {C['text_mute']}; }}
QPushButton:disabled {{
    background: {C['surface']};
    color: {C['text_mute']};
    border-color: {C['border']};
}}

/* ---- error: marked where it is, not only in a dialog ---- */
QLineEdit[error="true"], QComboBox[error="true"], QSpinBox[error="true"],
QDoubleSpinBox[error="true"] {{
    border: 1px solid {C['crit']};
    background: {C['void']};
}}
QGroupBox[error="true"] {{ border: 1px solid {C['crit']}; }}
QGroupBox[error="true"]::title {{ color: {C['crit']}; }}

/* ---- Required / Optional beside a panel's heading ---- */
QLabel#Tag {{
    border-radius: {r}px;
    padding: 1px 7px;
    font-size: {micro}px;
    font-weight: 600;
    color: {C['text_mute']};
    border: 1px solid {C['border_2']};
    background: transparent;
}}
QLabel#Tag[kind="required"] {{
    color: {C['accent']};
    border: 1px solid {C['accent']};
    background: {C['accent_soft']};
}}
QLabel#Tag[kind="unused"] {{
    color: {C['text_mute']};
    border: 1px dashed {C['border_2']};
}}

/* ---- views Qt builds itself: the file chooser, mostly ---- */
QTreeView, QListView {{
    background: {C['void']};
    color: {C['text']};
    border: 1px solid {C['border']};
    alternate-background-color: {C['surface']};
    selection-background-color: {C['accent']};
    selection-color: {C['ink']};
    outline: none;
}}
QTreeView::item, QListView::item {{ padding: 3px 4px; }}
QTreeView::item:hover:!selected, QListView::item:hover:!selected {{
    background: {C['surface_2']};
}}
QTreeView::item:selected, QListView::item:selected {{
    background: {C['accent']}; color: {C['ink']};
}}
QFileDialog QToolButton {{
    background: {C['surface_2']}; border: 1px solid {C['border_2']};
    border-radius: {r}px; padding: 3px;
}}
QFileDialog QToolButton:hover {{ background: {C['surface_3']}; }}
QMessageBox QLabel {{ color: {C['text']}; }}

QTableWidget::item:hover {{ background: {C['surface_2']}; }}
QTableWidget::item:selected {{ background: {C['accent']};
                               color: {C['ink']}; }}

/* ================= the current look =================
   Later rules win over the earlier ones of the same selector; these are the
   interface as it is drawn now, in the website's language: warm neutrals,
   capsules for anything pressed, 10 px fields, 14 px panels, and cobalt
   kept for small signals. */

#Plate, #Card {{
    background: {C['surface']}; border: 1px solid {C['border']};
    border-top-color: {C['rim']}; border-radius: {rp}px;
}}
#Header, #StepBar {{ background: {C['surface']}; border-bottom: 1px solid {C['border']}; }}
#StatusBar {{ background: {C['surface']}; border-top: 1px solid {C['border']}; }}
QLabel#StepTitle {{
    font-family: "{_pick(_DISPLAY, 'display')}"; font-size: {SIZE['heading'] + 11}px;
    font-weight: 600; letter-spacing: -0.2px; color: {C['text']};
}}
QLabel#StepLede {{ font-size: {SIZE['body'] + 1}px; color: {C['text_dim']}; }}
#SectionTitle {{ letter-spacing: 0; font-weight: 600; }}

/* buttons: capsules */
QPushButton {{
    background: {C['void']}; color: {C['text']};
    border: 1px solid {C['border_2']}; border-top-color: {C['rim'] if C['name'] == 'dark' else C['border_2']};
    border-radius: {ROW // 2}px; padding: 0 16px; min-height: {ROW}px;
    font-size: {label + 1}px; font-weight: 500;
}}
QPushButton:hover   {{ background: {C['surface_2']}; border-color: {C['text_mute']}; }}
QPushButton:pressed {{ background: {C['surface_3']}; }}
QPushButton:disabled {{ background: {C['surface']}; color: {C['text_mute']};
                        border-color: {C['border']}; }}
QPushButton#Ghost {{ background: transparent; border: 1px solid {C['border']};
                     color: {C['text_dim']}; border-radius: {ROW // 2}px; }}
QPushButton#Ghost:hover {{ background: {C['surface_2']}; color: {C['text']};
                           border-color: {C['border_2']}; }}
/* the primary action: the strongest neutral, not a blue slab */
QPushButton#Primary, QPushButton#Next {{
    background: {C['text']}; color: {C['ink'] if C['name'] == 'dark' else C['void']};
    border: 1px solid {C['text']}; font-weight: 600;
}}
QPushButton#Primary {{ border-radius: {ROW_TALL // 2}px; min-height: {ROW_TALL}px;
                       padding: 0 20px; font-size: {SIZE['heading']}px; }}
QPushButton#Next {{ border-radius: {ROW // 2}px; padding: 0 14px 0 16px; }}
QPushButton#Primary:hover, QPushButton#Next:hover {{
    background: {C['text_dim']}; border-color: {C['text_dim']};
}}
QPushButton#Primary:pressed, QPushButton#Next:pressed {{ background: {C['text_mute']}; }}
QPushButton#Primary:disabled, QPushButton#Next:disabled {{
    background: {C['surface_3']}; border-color: {C['surface_3']}; color: {C['text_mute']};
}}
QPushButton#Primary:focus, QPushButton#Next:focus {{ border: 2px solid {C['accent']}; }}
QPushButton#Tool {{ border-radius: {r}px; background: transparent;
                    border: 1px solid transparent; padding: 0; color: {C['text_dim']}; }}
QPushButton#Tool:hover {{ background: {C['surface_2']}; border-color: {C['border']};
                          color: {C['text']}; }}
QPushButton#Tool:checked {{ background: {C['surface_3']}; border-color: {C['border_2']};
                            color: {C['text']}; }}
QPushButton#Chip {{ border-radius: 14px; min-height: 28px; padding: 0 14px;
                    background: {C['void']}; color: {C['text_dim']}; }}
QPushButton#Chip:checked {{ background: {C['text']}; border-color: {C['text']};
                            color: {C['ink'] if C['name'] == 'dark' else C['void']}; }}
QPushButton#Step {{ border-radius: {r}px; }}
QPushButton#NavOn {{ background: {C['surface_3']}; border: 1px solid transparent;
                     border-radius: {r}px; color: {C['text']}; }}
QPushButton#Nav {{ border: 1px solid transparent; border-radius: {r}px; }}
QPushButton#Nav:focus {{ border: 1px solid {C['accent']}; }}

/* segmented controls: a soft track with the choice lifted out of it as a pill */
#SegTrack {{
    background: {C['surface_2']}; border: 1px solid {C['border']};
    border-radius: {ROW // 2 + 1}px;
}}
QPushButton#Seg, QPushButton#SegL, QPushButton#SegR {{
    background: transparent; color: {C['text_dim']};
    border: 1px solid transparent; border-radius: {ROW // 2 - 3}px;
    min-height: {ROW - 6}px; padding: 0 12px; font-weight: 500;
}}
QPushButton#Seg:hover:!checked, QPushButton#SegL:hover:!checked,
QPushButton#SegR:hover:!checked {{ background: {C['surface_3']}; color: {C['text']}; }}
QPushButton#Seg:checked, QPushButton#SegL:checked, QPushButton#SegR:checked {{
    background: {C['void'] if C['name'] == 'light' else C['surface_3']};
    color: {C['text']}; font-weight: 600;
    border: 1px solid {C['border_2'] if C['name'] == 'light' else C['rim']};
}}
QPushButton#Seg:disabled, QPushButton#SegL:disabled, QPushButton#SegR:disabled {{
    background: transparent; color: {C['text_mute']}; border-color: transparent;
}}
QPushButton#Seg:focus, QPushButton#SegL:focus, QPushButton#SegR:focus {{
    border: 1px solid {C['text_mute']};
}}

/* fields */
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
    background: {C['void']}; border: 1px solid {C['border_2']};
    border-radius: {r}px; min-height: {ROW}px; padding: 0 10px;
}}
QLineEdit:hover, QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover {{
    border-color: {C['text_mute']};
}}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus,
QComboBox:on {{ border: 1px solid {C['accent']}; }}
QSpinBox::up-button, QDoubleSpinBox::up-button {{
    subcontrol-origin: border; subcontrol-position: top right;
    width: 22px; border: none; border-left: 1px solid {C['border']};
    border-top-right-radius: {r}px; background: transparent;
}}
QSpinBox::down-button, QDoubleSpinBox::down-button {{
    subcontrol-origin: border; subcontrol-position: bottom right;
    width: 22px; border: none; border-left: 1px solid {C['border']};
    border-bottom-right-radius: {r}px; background: transparent;
}}
QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover,
QSpinBox::down-button:hover, QDoubleSpinBox::down-button:hover {{
    background: {C['surface_2']};
}}

/* drop-downs: a list that opens under the box, not over it */
QComboBox {{ combobox-popup: 0; padding-right: 28px; }}
QComboBox::drop-down {{ border: none; width: 28px; }}
QComboBox QAbstractItemView {{
    background: {C['void']}; border: 1px solid {C['border_2']};
    border-radius: {rp - 2}px; padding: 4px; outline: none;
    selection-background-color: {C['surface_3']};
    selection-color: {C['text']};
}}
QComboBox QAbstractItemView::item {{
    min-height: {ROW - 2}px; padding: 0 10px; border-radius: {r - 3}px;
    border: none; color: {C['text']};
}}
QComboBox QAbstractItemView::item:hover {{ background: {C['surface_2']}; }}
QComboBox QAbstractItemView::item:selected {{
    background: {C['surface_3']}; color: {C['text']}; border: none;
}}

/* check boxes and radios: cobalt when on, a small signal */
QCheckBox::indicator, QCheckBox#PlateSwitch::indicator {{
    width: 17px; height: 17px; border-radius: 5px;
    border: 1px solid {C['border_2']}; background: {C['void']};
}}
QCheckBox::indicator:hover {{ border-color: {C['text_mute']}; }}
QCheckBox::indicator:checked, QCheckBox#PlateSwitch::indicator:checked {{
    background: {C['accent']}; border-color: {C['accent']};
    image: url({style_icon('check', C['ink'] if C['name'] == 'dark' else C['void'])});
}}
QRadioButton::indicator {{ width: 17px; height: 17px; border-radius: 9px; }}
QRadioButton::indicator:checked {{ border: 5px solid {C['accent']};
                                   background: {C['void']}; }}

/* scroll bars: slim capsules, out of the way until used */
QScrollBar:vertical {{ background: transparent; width: 10px; margin: 2px; }}
QScrollBar:horizontal {{ background: transparent; height: 10px; margin: 2px; }}
QScrollBar::handle:vertical {{ background: {C['border_2']}; border-radius: 3px;
                               min-height: 36px; }}
QScrollBar::handle:horizontal {{ background: {C['border_2']}; border-radius: 3px;
                                 min-width: 36px; }}
QScrollBar::handle:hover {{ background: {C['text_mute']}; }}
QScrollBar::handle:pressed {{ background: {C['text_dim']}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ width: 0; height: 0; border: none; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}

/* lists and tables: wells, the selection a neutral highlight */
QListWidget {{ border-radius: {r}px; border: 1px solid {C['border']}; padding: 4px; }}
QListWidget::item {{ border-radius: {r - 3}px; padding: 6px 8px; border: none; }}
QListWidget::item:selected {{ background: {C['surface_3']}; color: {C['text']};
                              border: none; }}
QListWidget::item:hover:!selected {{ background: {C['surface_2']}; }}
QTableWidget {{ border-radius: {r}px; selection-background-color: {C['surface_3']};
                selection-color: {C['text']}; }}
QTableWidget::item:selected {{ background: {C['surface_3']}; color: {C['text']}; }}
QHeaderView::section {{ padding: 6px 8px; font-size: {label}px; background: {C['surface_2']};
                        color: {C['text_dim']}; }}
QTreeView::item:selected, QListView::item:selected {{
    background: {C['surface_3']}; color: {C['text']};
}}
QTextEdit, QPlainTextEdit {{ border-radius: {r}px; padding: 8px; }}

QProgressBar {{ border: none; border-radius: 4px; background: {C['surface_3']};
                min-height: 8px; max-height: 8px; color: transparent; }}
QProgressBar::chunk {{ background: {C['accent']}; border-radius: 4px; }}

/* menus and tooltips */
QMenu {{ background: {C['void']}; border: 1px solid {C['border_2']};
         border-radius: {rp - 2}px; padding: 4px; }}
QMenu::item {{ border-radius: {r - 3}px; padding: 7px 24px 7px 12px; color: {C['text']}; }}
QMenu::item:selected {{ background: {C['surface_3']}; color: {C['text']}; }}
QToolTip {{ background: {C['text']}; color: {C['bg']}; border: none;
            border-radius: {r - 2}px; padding: 5px 9px; }}

/* tabs: an underline, not a folder */
QTabWidget::pane {{ border: none; border-top: 1px solid {C['border']}; }}
QTabBar::tab {{ background: transparent; border: none; padding: 8px 14px;
                border-bottom: 2px solid transparent; color: {C['text_dim']}; }}
QTabBar::tab:selected {{ color: {C['text']}; border-bottom: 2px solid {C['text']}; }}
QTabBar::tab:hover {{ color: {C['text']}; }}

/* tags beside a panel's heading: quiet capsules */
QLabel#Tag {{ border-radius: 9px; padding: 1px 8px; background: transparent;
              color: {C['text_mute']}; border: 1px solid {C['border_2']}; }}
QLabel#Tag[kind="required"] {{ color: {C['text']}; background: {C['surface_3']};
                               border: 1px solid {C['surface_3']}; }}
QLabel#Tag[kind="unused"] {{ color: {C['text_mute']}; border: 1px dashed {C['border_2']}; }}

QSplitter::handle {{ background: {C['border']}; }}
QSplitter::handle:hover {{ background: {C['border_2']}; }}
QSplitter::handle:pressed {{ background: {C['accent']}; }}
"""
