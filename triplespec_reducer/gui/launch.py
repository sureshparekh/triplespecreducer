"""
Start the desktop interface behind a splash screen.

The splash is shown before the reduction code is imported (about two
seconds of numerical libraries), in the theme last used, and stays at least
SPLASH_SECONDS so the logo is seen even on a fast machine. Only Qt is
imported until the splash is on screen.
"""

import os
import sys
import time

SPLASH_SECONDS = 3.0


def main(argv=None):
    argv = argv or sys.argv
    os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")
    from PySide6.QtCore import Qt, QSettings, QTimer
    from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap, QFont
    from PySide6.QtWidgets import QApplication, QSplashScreen

    app = QApplication.instance() or QApplication(argv)
    app.setApplicationName("TripleSpec Reducer")
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    logo = os.path.join(here, "data", "logo")
    app.setWindowIcon(QIcon(os.path.join(logo, "app_icon.png")))
    theme = str(QSettings("SOAR TripleSpec Reducer", "TripleSpec Reducer").value("theme", "dark"))
    light = theme == "light"
    bg = QColor(245, 247, 250) if light else QColor(10, 14, 20)
    fg = QColor(70, 80, 96) if light else QColor(150, 162, 180)

    # the logo, with a band beneath it for the loading line
    dpr = app.primaryScreen().devicePixelRatio() if app.primaryScreen() else 1.0
    side, band = 460, 54
    art = QPixmap(os.path.join(logo, f"splash_{'light' if light else 'dark'}.png")).scaled(
        int(side * dpr), int(side * dpr), Qt.KeepAspectRatio, Qt.SmoothTransformation)
    pm = QPixmap(int(side * dpr), int((side + band) * dpr))
    pm.fill(bg)
    p = QPainter(pm)
    p.drawPixmap(0, 0, art)
    p.end()
    pm.setDevicePixelRatio(dpr)
    splash = QSplashScreen(pm, Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint)
    f = QFont()
    f.setPointSizeF(10.5)
    splash.setFont(f)

    def say(text):
        splash.showMessage(text, Qt.AlignHCenter | Qt.AlignBottom, fg)
        app.processEvents()

    t0 = time.monotonic()
    splash.show()
    say("Loading the reduction modules ...")
    from triplespec_reducer.gui import theme as TH
    from triplespec_reducer.gui.app import TSR
    say("Preparing the interface ...")
    TH.load_fonts()
    w = TSR()
    if len(argv) > 1 and os.path.isdir(argv[1]):
        w.ed_raw.setText(argv[1])
    say("Ready")

    def show():
        w.show()
        splash.finish(w)

    QTimer.singleShot(max(0, int((SPLASH_SECONDS - (time.monotonic() - t0)) * 1000)), show)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
