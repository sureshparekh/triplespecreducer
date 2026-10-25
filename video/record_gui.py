"""
Film 2's screen recording: the real interface, driven through the tutorial's
beats (video/tutorial_script.py) on a fresh output folder, captured frame by
frame offscreen.

    python video/record_gui.py            (from the repository's root, after video/voice.py tutorial)

Everything the program computes is computed: reading the night, the flat,
the wavelength solution, the extractions, the telluric correction. The
waiting is recorded as it happens.

Writes video/build/tutorial/rec/: the frames (JPEG, written only when the
window changed), and rec.json: every sample (its recording time, frame,
cursor), the clicks, the camera's zoom targets, the beats and the waits.
"""

import hashlib
import importlib.util
import json
import os
import shutil
import sys
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

from PySide6.QtCore import QPoint, QRect                                    # noqa: E402
from PySide6.QtGui import QPixmap                                           # noqa: E402
from PySide6.QtWidgets import (QApplication, QMessageBox, QPushButton, QScrollArea,   # noqa: E402
                               QCheckBox, QComboBox, QLineEdit, QListWidget, QTableWidget)

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tutorial_script", os.path.join(HERE, "tutorial_script.py"))
SCRIPT = importlib.util.module_from_spec(spec)
spec.loader.exec_module(SCRIPT)
VOICE = os.path.join(HERE, "build", "tutorial", "voice")
REC = os.path.join(HERE, "build", "tutorial", "rec")
FPS = 30
WAIT_FPS = 10                       # while a step computes, samples are taken in real time
W, H = 1600, 1000

app = QApplication([])
# a dialog would block an offscreen window: report it instead
for name in ("information", "warning", "critical", "question"):
    setattr(QMessageBox, name, staticmethod(lambda *a, **k: print("dialog:", a[2] if len(a) > 2 else a) or QMessageBox.Ok))

from triplespec_reducer.gui import app as A                                 # noqa: E402
from triplespec_reducer.gui import theme as TH                              # noqa: E402
A.failure_dialog = lambda parent, text, **k: print("FAILED:", text[-2000:])

shutil.rmtree(REC, ignore_errors=True)
os.makedirs(REC)
shutil.rmtree(os.path.join(ROOT, SCRIPT.OUT), ignore_errors=True)

samples, clicks, zooms, beats, waits, keys = [], [], [], [], [], []
state = dict(t=0.0, cur=(W * 0.62, H * 0.55), last=None, nframes=0, img=None)


def frame_file(qimg):
    """Save the window's image if it changed; return its file name."""
    h = hashlib.md5(bytes(qimg.constBits())).hexdigest()
    if h != state["last"]:
        state["last"] = h
        state["img"] = f"f{state['nframes']:06d}.jpg"
        qimg.save(os.path.join(REC, state["img"]), "JPG", 90)
        state["nframes"] += 1
    return state["img"]


def capture(dt=1.0 / FPS, kind="win"):
    app.processEvents()
    if kind == "win":
        img = frame_file(w.grab().toImage())
    else:
        img = kind
    samples.append([round(state["t"], 4), img, round(state["cur"][0], 1), round(state["cur"][1], 1)])
    state["t"] += dt


def hold(sec):
    for _ in range(max(1, int(round(sec * FPS)))):
        capture()


def rect_of(wd):
    p = wd.mapTo(w, QPoint(0, 0))
    vis = wd.visibleRegion().boundingRect()
    if not vis.isEmpty():
        p = wd.mapTo(w, vis.topLeft())
        return QRect(p, vis.size())
    return QRect(p, wd.size())


def resolve(tg):
    """A target -> (widget or None, rect in window pixels, what clicking it does)."""
    if tg is None:
        return None, QRect(0, 0, W, H), None
    if isinstance(tg, str):
        if tg == "STAGE":
            return w.stage, rect_of(w.stage), None
        if tg == "BAR":
            return w.ws.rail, rect_of(w.ws.rail), None
        if tg == "ABOUT":
            p = w.about_plates[w.ws.current()]
            return p, rect_of(p), None
        if tg == "TOPRIGHT":
            return None, QRect(W - 420, 0, 420, 64), None
        if tg == "panel":
            pg = w.ws.page(w.ws.current())
            return pg, rect_of(pg), None
        wd = getattr(w, tg)
        return wd, rect_of(wd), None
    kind = tg[0]
    if kind == "button":
        want = tg[1]
        found = [b for b in w.findChildren(QPushButton)
                 if b.isVisible() and b.text().replace("&&", "&").strip() == want]
        if not found:
            raise SystemExit(f"no button {want!r} on the visible page")
        b = found[0]
        return b, rect_of(b), b.click
    if kind == "rail":
        it = w.ws.rail.items[tg[1]]
        return it, rect_of(it), it.clicked.emit
    if kind == "about":
        h = w.about_panels[tg[1]][tg[2]].head
        return h, rect_of(h), h.click
    if kind == "seg":
        for b in w.seg_page.findChildren(QPushButton):
            if b.text() == tg[1]:
                return b, rect_of(b), lambda: w.seg_page.setValue(tg[1]) or w.pages.setCurrentIndex(
                    0 if tg[1] == "Reduce" else 1)
        return w.seg_page, rect_of(w.seg_page), lambda: w.seg_page.setValue(tg[1])
    if kind == "row":
        tbl = getattr(w, tg[1])
        for r in range(tbl.rowCount()):
            if tbl.item(r, 0) is not None and tbl.item(r, 0).text().strip() == tg[2]:
                ir = tbl.visualItemRect(tbl.item(r, 0))
                p = tbl.viewport().mapTo(w, ir.topLeft())
                return tbl, QRect(p, ir.size()), lambda r=r: (tbl.selectRow(r), tbl.setCurrentCell(r, 0))
        raise SystemExit(f"no row {tg[2]!r} in {tg[1]}")
    if kind == "lastitem":
        lst = getattr(w, tg[1])
        n = lst.count() - 1
        ir = lst.visualItemRect(lst.item(n))
        p = lst.viewport().mapTo(w, ir.topLeft())
        return lst, QRect(p, ir.size()), lambda: lst.setCurrentRow(n)
    raise SystemExit(f"unknown target {tg!r}")


def ensure_visible(wd):
    """If wd lies in a scrolled panel outside the view, scroll the panel to it, as a person would."""
    if wd is None:
        return
    for sa in w.findChildren(QScrollArea):
        inner = sa.widget()
        if inner is None or not sa.isVisible() or not inner.isAncestorOf(wd):
            continue
        top = wd.mapTo(inner, QPoint(0, 0)).y()
        bar, view = sa.verticalScrollBar(), sa.viewport().height()
        if top >= bar.value() + 20 and top + wd.height() <= bar.value() + view - 20:
            return
        v0, v1 = bar.value(), max(0, min(bar.maximum(), top - view // 3))
        n = int(0.55 * FPS)
        for i in range(1, n + 1):
            u = i / n
            bar.setValue(int(v0 + (v1 - v0) * u * u * (3 - 2 * u)))
            capture()
        return


def glide(x1, y1, speed=1.0):
    x0, y0 = state["cur"]
    d = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5
    if d < 2:
        return
    n = max(8, int(round((0.45 + d / 1800) / speed * FPS)))
    for i in range(1, n + 1):
        u = i / n
        e = u * u * (3 - 2 * u)
        # a slight arc, as a hand moves
        bow = 0.06 * d * 4 * u * (1 - u)
        state["cur"] = (x0 + (x1 - x0) * e, y0 + (y1 - y0) * e - bow)
        capture()


def point_in(rc):
    return rc.x() + min(rc.width() * 0.5, 60 + rc.width() * 0.08) if rc.width() > 260 else rc.center().x(), \
        rc.center().y()


def click(tg):
    ensure_visible(resolve(tg)[0])
    wd, rc, act = resolve(tg)
    x, y = point_in(rc)
    glide(x, y)
    hold(0.12)
    clicks.append(round(state["t"], 4))
    if act is None and wd is not None:
        if isinstance(wd, QPushButton):
            act = wd.click
        elif isinstance(wd, QCheckBox):
            act = lambda: wd.setChecked(not wd.isChecked())
    if act:
        act()
    hold(0.25)


def type_into(tg, text):
    ensure_visible(resolve(tg)[0])
    wd, rc, _ = resolve(tg)
    glide(rc.x() + 40, rc.center().y())
    hold(0.1)
    clicks.append(round(state["t"], 4))
    wd.setFocus()
    wd.selectAll()
    hold(0.25)
    wd.clear()
    for i in range(1, len(text) + 1):
        wd.setText(text[:i])
        keys.append(round(state["t"], 4))
        capture(1.0 / 16)              # about sixteen characters a second
    hold(0.3)


def select(tg, text):
    ensure_visible(resolve(tg)[0])
    wd, rc, _ = resolve(tg)
    glide(*point_in(rc))
    hold(0.12)
    clicks.append(round(state["t"], 4))
    hold(0.25)
    wd.setCurrentText(text)
    hold(0.4)


def scroll_to(step, tg):
    wd, _, _ = resolve(tg)
    for sa in w.ws.page(step).findChildren(QScrollArea):
        if sa.widget() is not None and sa.widget().isAncestorOf(wd):
            bar = sa.verticalScrollBar()
            v0, v1 = bar.value(), max(0, min(bar.maximum(), wd.mapTo(sa.widget(), QPoint(0, 0)).y() - 110))
            n = int(0.6 * FPS)
            for i in range(1, n + 1):
                u = i / n
                bar.setValue(int(v0 + (v1 - v0) * u * u * (3 - 2 * u)))
                capture()


def wait_job(beat):
    t0 = state["t"]
    real0 = time.time()
    last = time.time()
    while w.job is not None and w.job.isRunning():
        app.processEvents()
        now = time.time()
        if now - last >= 1.0 / WAIT_FPS:
            capture(dt=now - last)
            last = now
        else:
            time.sleep(0.005)
    for _ in range(10):
        app.processEvents()
    capture(dt=max(time.time() - last, 1.0 / FPS))
    waits.append(dict(beat=beat, t0=round(t0, 3), t1=round(state["t"], 3), real=round(time.time() - real0, 1)))
    print(f"    waited {time.time() - real0:5.1f} s", flush=True)


def zoom(tg):
    _, rc, _ = resolve(tg)
    zooms.append(dict(t=round(state["t"], 4), rect=[rc.x(), rc.y(), rc.width(), rc.height()],
                      full=tg is None))


def launch():
    """The splash, as the program shows it, then the window."""
    import triplespec_reducer.gui.launch as L                                # noqa: F401
    pm = QPixmap(os.path.join(A.LOGO_DIR, "splash_dark.png"))
    splash = os.path.join(REC, "splash.png")
    pm.save(splash)
    for _ in range(int(3.0 * FPS)):
        samples.append([round(state["t"], 4), "@splash", -1, -1])
        state["t"] += 1.0 / FPS
    w.show()
    hold(0.6)


def word_time(beat_voice, word):
    """Seconds into the beat's narration when `word` (or 'word#n') is spoken."""
    if word == "end":
        return beat_voice["duration"]
    wd, _, nth = word.partition("#")
    n = int(nth or 1)
    hits = [x for l in beat_voice["lines"] for x in l["words"]
            if x[0].strip(".,;:!?'\"").lower() == wd.lower()]
    if len(hits) < n:
        print(f"    (word {word!r} not spoken; acting now)")
        return 0.0
    return hits[n - 1][1]


TH.load_fonts()
w = A.TSR()
w.apply_theme(TH.LIGHT)
w.resize(W, H)
LEAD = 0.5

for bi, b in enumerate(SCRIPT.BEATS):
    v = json.load(open(os.path.join(VOICE, b["id"] + ".json")))
    start = state["t"]
    if not w.isVisible() and b["actions"] and b["actions"][0][0] != "launch":
        pass
    beats.append(dict(id=b["id"], step=b["step"], t0=round(start, 4), voice=v["duration"], lead=LEAD))
    print(f"{b['id']:12s} t={start:7.1f}  narration {v['duration']:5.1f} s", flush=True)
    if not samples:                                # before the program opens: the title card
        for _ in range(int(LEAD * FPS)):
            samples.append([round(state["t"], 4), "@title", -1, -1])
            state["t"] += 1.0 / FPS
    else:
        hold(LEAD) if w.isVisible() else None
    for act in b["actions"]:
        op = act[0]
        at = None
        if op in ("move", "click", "zoom") and len(act) > 2:
            at = act[2]
        elif op in ("type", "select", "scroll") and len(act) > 3:
            at = act[3]
        if at is not None:
            due = start + LEAD + word_time(v, at)
            while state["t"] < due:
                if w.isVisible():
                    capture()
                else:
                    samples.append([round(state["t"], 4), "@title", -1, -1])
                    state["t"] += 1.0 / FPS
        if op == "launch":
            due = start + LEAD + word_time(v, "logo")
            while state["t"] < due:
                samples.append([round(state["t"], 4), "@title", -1, -1])
                state["t"] += 1.0 / FPS
            launch()
        elif op == "move":
            if act[1] not in ("STAGE", "BAR", "panel"):
                ensure_visible(resolve(act[1])[0])
            _, rc, _ = resolve(act[1])
            glide(*point_in(rc))
        elif op == "click":
            click(act[1])
        elif op == "type":
            type_into(act[1], act[2])
        elif op == "select":
            select(act[1], act[2])
        elif op == "scroll":
            scroll_to(act[1], act[2])
        elif op == "wait":
            wait_job(b["id"])
        elif op == "pause":
            hold(act[1])
        elif op == "zoom":
            zoom(act[1])
        else:
            raise SystemExit(f"unknown action {op}")
    # the beat lasts as long as its narration at least
    due = start + LEAD + v["duration"] + 0.35
    while state["t"] < due:
        if w.isVisible():
            capture()
        else:
            samples.append([round(state["t"], 4), "@title", -1, -1])
            state["t"] += 1.0 / FPS
    beats[-1]["t1"] = round(state["t"], 4)

hold(1.5)
json.dump(dict(fps=FPS, size=[W, H], samples=samples, clicks=clicks, keys=keys, zooms=zooms, beats=beats,
               waits=waits, duration=round(state["t"], 3)),
          open(os.path.join(REC, "rec.json"), "w"))
print(f"recorded {state['t'] / 60:.1f} min, {state['nframes']} distinct frames, {len(samples)} samples")
