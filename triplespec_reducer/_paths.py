"""Where the package's own files are, wherever it is installed."""

import os

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(PACKAGE_DIR, "data")
INSTRUMENT_DIR = os.path.join(DATA_DIR, "instrument")
SPEXTOOL_DIR = os.path.join(DATA_DIR, "spextool")
FONT_DIR = os.path.join(DATA_DIR, "fonts")
ICON_DIR = os.path.join(DATA_DIR, "icons")
LOGO_DIR = os.path.join(DATA_DIR, "logo")
