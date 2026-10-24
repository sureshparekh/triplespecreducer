"""Sphinx configuration for the TripleSpec Reducer documentation (BRAIN's design)."""

import os
import shutil
import sys
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(HERE, "_ext"))

project = "TripleSpec Reducer"
author = "Suresh Parekh"
copyright = f"{date.today().year}, Suresh Parekh (UFRGS)"
release = "1.0.0"
version = "1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
    "tsr_ext",
]

source_suffix = {".md": "markdown"}
master_doc = "index"
exclude_patterns = ["_build", "_ext", "_themes", "README.md"]

myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist",
                          "attrs_inline", "attrs_block", "fieldlist"]
myst_heading_anchors = 3
suppress_warnings = ["myst.header"]

autodoc_member_order = "bysource"
autodoc_typehints = "none"
autodoc_mock_imports = ["PySide6", "numba"]
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

html_theme = "tsr"
html_theme_path = ["_themes"]
html_title = "TripleSpec Reducer documentation"
html_static_path = ["_static"]
html_show_sourcelink = False
html_copy_source = False
html_show_sphinx = False
html_last_updated_fmt = "%d %B %Y"
html_permalinks_icon = "#"

# The top bar's tabs, and in each the sidebar's groups of pages, in reading
# order (which is also the order of Previous / Next).
tsr_nav = [
    {"title": "Get started", "slug": "start", "index": "start/index", "prefixes": ["start/"], "groups": [
        ("Introduction", ["start/index", "start/install"]),
        ("Before you reduce", ["start/data", "start/quickstart"]),
    ]},
    {"title": "Tutorial", "slug": "tutorial", "index": "tutorial/index", "prefixes": ["tutorial/"], "groups": [
        ("One galaxy, start to finish", ["tutorial/index", "tutorial/1-night", "tutorial/2-calibrate",
                                         "tutorial/3-extract", "tutorial/4-combine", "tutorial/5-telluric",
                                         "tutorial/6-merge", "tutorial/7-clean", "tutorial/8-products"]),
        ("Watch it", ["tutorial/video"]),
    ]},
    {"title": "Interface", "slug": "gui", "index": "gui/index", "prefixes": ["gui/"], "groups": [
        ("The program", ["gui/index", "gui/settings"]),
        ("The night and its files", ["gui/automation", "gui/naming"]),
        ("Without the interface", ["gui/cli", "gui/options"]),
    ]},
    {"title": "Science", "slug": "science", "index": "science/index", "prefixes": ["science/"], "groups": [
        ("The data", ["science/index", "science/instrument"]),
        ("Calibrations", ["science/flat", "science/wavecal"]),
        ("Spectra", ["science/extraction", "science/combine"]),
        ("Flux and atmosphere", ["science/telluric", "science/merge"]),
    ]},
    {"title": "Validation", "slug": "validation", "index": "validation/index", "prefixes": ["validation/"],
     "groups": [
        ("Against IDL", ["validation/index", "validation/idl_differences"]),
    ]},
    {"title": "Reference", "slug": "reference", "index": "reference/index", "prefixes": ["reference/"], "groups": [
        ("Overview", ["reference/index"]),
        ("Python API", ["reference/api"]),
        ("Help", ["reference/troubleshooting", "reference/citing", "reference/sources"]),
    ]},
]
tsr_parents = {}

STEP_OF = {"night": 0, "cal": 1, "ex": 2, "cb": 3, "tc": 4, "mg": 5, "cl": 6}


def _write_settings():
    """gui/settings.md, from the interface's own explanations (gui/help.py)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "tsr_help", os.path.join(ROOT, "triplespec_reducer", "gui", "help.py"))
    H = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(H)
    out = ["# Every setting", "",
           "Every control of the desktop interface, step by step: what it is, what it does, and where "
           "the IDL tool starts it. The same text appears in the interface, as the caption under each "
           "control and as its tooltip; this page is generated from it "
           "(`triplespec_reducer/gui/help.py`), so the two always agree.", ""]
    for i, st in enumerate(H.STEPS):
        out += [f"## {i + 1}. {st['title']}", "", st["why"], ""]
        keys = [k for k in H.FIELDS if STEP_OF[k.split(".")[0]] == i]
        if keys:
            out += ["| Setting | What it does | Default |", "|---|---|---|"]
            for k in keys:
                d = H.field(k)
                text = d["caption"]
                if d["details"] and d["details"] != d["caption"]:
                    text += " " + d["details"]
                text = text.replace("|", "\\|").replace("\n", " ")
                out.append(f"| **{d['label'] or '(positions)'}** | {text} | {d['default'] or '—'} |")
            out.append("")
        out += [f"*Reproduces:* {st['idl']}", ""]
    path = os.path.join(HERE, "gui", "settings.md")
    text = "\n".join(out) + "\n"
    if not os.path.exists(path) or open(path).read() != text:
        open(path, "w").write(text)


def _copy_logos():
    src = os.path.join(ROOT, "triplespec_reducer", "data", "logo")
    dst = os.path.join(HERE, "_static", "logo")
    os.makedirs(dst, exist_ok=True)
    for f in os.listdir(src):
        shutil.copy2(os.path.join(src, f), os.path.join(dst, f))


_write_settings()
_copy_logos()


# The modules' docstrings are written as plain text, for a reader of the
# source: aligned tables of steps and formulas. A module docstring is shown
# as written; elsewhere an indented block after a paragraph is shown as a
# literal block, as BRAIN's documentation does.
def _module_docstring_as_text(app, what, name, obj, options, lines):
    if what == "module" and lines:
        lines[:] = ["::", ""] + ["    " + l for l in lines] + [""]


def _literal_blocks(app, what, name, obj, options, lines):
    if any(set(l.strip()) == {"-"} and len(l.strip()) >= 3 for l in lines):
        return
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        indented = line.startswith((" ", "\t")) and line.strip()
        prev = out[-1] if out else ""
        text_before = [l for l in out if l.strip()]
        after_text = bool(text_before) and \
            not text_before[-1].startswith((" ", "\t")) and \
            not text_before[-1].rstrip().endswith("::")
        if indented and after_text:
            out += ["", "::", ""]
            while i < len(lines) and (lines[i].startswith((" ", "\t"))
                                      or not lines[i].strip()):
                if not lines[i].strip() and i + 1 < len(lines) and \
                        lines[i + 1].strip() and \
                        not lines[i + 1].startswith((" ", "\t")):
                    break
                out.append(lines[i])
                i += 1
            out.append("")
            continue
        out.append(line)
        i += 1
    lines[:] = out


def setup(app):
    app.connect("autodoc-process-docstring", _module_docstring_as_text, priority=300)
    app.connect("autodoc-process-docstring", _literal_blocks, priority=400)
