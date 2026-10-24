"""
The TripleSpec Reducer documentation's own Sphinx pieces (after BRAIN's).

* Navigation: NAV in conf.py (tabs, groups, pages) becomes the top bar, the
  sidebar, the breadcrumbs and the previous/next links of every page, with
  each page's title read from the page itself.
* ``cards``/``linkcard``: a grid of linked cards with a line icon.
* ``diagram``: an SVG from _diagrams/, inlined, so it takes the page's
  colours in the light and the dark theme.
* ``shots``: a screenshot of the interface in the light theme and the dark
  one; the reader sees the one matching the page's theme.
"""

import hashlib
import html
import os
import re

from docutils import nodes
from docutils.parsers.rst import Directive, directives
from sphinx.util.docutils import SphinxDirective

# --------------------------------------------------------------------------
# icons: 24 x 24, stroked with currentColor (drawn for these pages)
# --------------------------------------------------------------------------
ICONS = {
    "window": '<rect x="3" y="4.5" width="18" height="15" rx="2.5"/><path d="M3 8.5h18M6 6.5h.01M8.5 6.5h.01"/><path d="M6.5 15.5l2.5-3 2 2 3-4 3.5 5"/>',
    "code": '<path d="M8.5 8l-4.5 4 4.5 4M15.5 8l4.5 4-4.5 4M13.5 5.5l-3 13"/>',
    "terminal": '<rect x="3" y="4.5" width="18" height="15" rx="2.5"/><path d="M7 10l3 2.5L7 15M12.5 15.5H17"/>',
    "spectrum": '<path d="M3 17c2 0 2.5-3.5 4-3.5s1.6 4.5 3 4.5 1.8-12 3.4-12 1.7 9 3.2 9 1.7-3 4.4-3"/><path d="M3 20.5h18"/>',
    "cube": '<path d="M12 3l8 4.5v9L12 21l-8-4.5v-9z"/><path d="M4 7.5l8 4.5 8-4.5M12 12v9"/>',
    "layers": '<path d="M12 3.5l9 4.5-9 4.5-9-4.5z"/><path d="M3 12.5l9 4.5 9-4.5M3 16.5l9 4.5 9-4.5"/>',
    "book": '<path d="M4.5 5.5c2.5-1.3 5-1.3 7.5 0v14c-2.5-1.3-5-1.3-7.5 0zM19.5 5.5c-2.5-1.3-5-1.3-7.5 0v14c2.5-1.3 5-1.3 7.5 0z"/>',
    "atom": '<circle cx="12" cy="12" r="1.6"/><ellipse cx="12" cy="12" rx="9" ry="3.6"/><ellipse cx="12" cy="12" rx="9" ry="3.6" transform="rotate(60 12 12)"/><ellipse cx="12" cy="12" rx="9" ry="3.6" transform="rotate(120 12 12)"/>',
    "gauge": '<path d="M4.5 17.5a8.5 8.5 0 1 1 15 0"/><path d="M12 13l4-4.5M12 13h.01"/>',
    "settings": '<circle cx="12" cy="12" r="3"/><path d="M12 2.8v2.4M12 18.8v2.4M2.8 12h2.4M18.8 12h2.4M5.5 5.5l1.7 1.7M16.8 16.8l1.7 1.7M5.5 18.5l1.7-1.7M16.8 7.2l1.7-1.7"/>',
    "list": '<path d="M9 6.5h11M9 12h11M9 17.5h11M4.5 6.5h.01M4.5 12h.01M4.5 17.5h.01"/>',
    "download": '<path d="M12 4v11M7.5 10.5L12 15l4.5-4.5M4.5 19.5h15"/>',
    "check": '<circle cx="12" cy="12" r="8.5"/><path d="M8 12.3l2.6 2.6L16 9.5"/>',
    "play": '<circle cx="12" cy="12" r="8.5"/><path d="M10 8.5l5.5 3.5-5.5 3.5z"/>',
    "map": '<rect x="3.5" y="3.5" width="17" height="17" rx="2.5"/><path d="M3.5 13.5h4v-4h4v6h4v-4h5"/>',
    "chart": '<path d="M4 4v16h16"/><path d="M7.5 15l3.5-4.5 3 2.5 5-6.5"/>',
    "flask": '<path d="M9.5 3.5h5M10.5 3.5v5.5L5 18.5a1.5 1.5 0 0 0 1.3 2.2h11.4a1.5 1.5 0 0 0 1.3-2.2L13.5 9V3.5"/><path d="M7.5 14.5h9"/>',
    "dust": '<circle cx="7" cy="8" r="1.2"/><circle cx="14" cy="6" r="1.6"/><circle cx="17.5" cy="12" r="1"/><circle cx="10" cy="13" r="2"/><circle cx="6" cy="17" r="1"/><circle cx="15" cy="17.5" r="1.4"/>',
    "lines": '<path d="M3 19h18"/><path d="M5 19c1.5 0 2-1 2.5-1s.7-12 1.5-12 .8 12 1.5 12h3c.7 0 .8-7 1.5-7s.8 7 1.5 7H19"/>',
    "sigma": '<path d="M17.5 5.5h-11l6 6.5-6 6.5h11"/>',
    "grid": '<rect x="4" y="4" width="6.5" height="6.5" rx="1.5"/><rect x="13.5" y="4" width="6.5" height="6.5" rx="1.5"/><rect x="4" y="13.5" width="6.5" height="6.5" rx="1.5"/><rect x="13.5" y="13.5" width="6.5" height="6.5" rx="1.5"/>',
    "cpu": '<rect x="6" y="6" width="12" height="12" rx="2"/><rect x="9.5" y="9.5" width="5" height="5" rx="1"/><path d="M9.5 3v3M14.5 3v3M9.5 18v3M14.5 18v3M3 9.5h3M3 14.5h3M18 9.5h3M18 14.5h3"/>',
    "mask": '<path d="M3 18h18"/><path d="M5 18V8h4v10M15 18V5h4v13"/><path d="M5 11l4 3M15 8l4 3M15 12l4 3"/>',
    "path": '<circle cx="6" cy="6" r="2"/><circle cx="18" cy="18" r="2"/><path d="M8 6h6.5a3.5 3.5 0 0 1 0 7h-5a3.5 3.5 0 0 0 0 7H16"/>',
    "help": '<circle cx="12" cy="12" r="8.5"/><path d="M9.6 9.5a2.5 2.5 0 1 1 3.4 2.3c-.6.3-1 .8-1 1.5v.7M12 16.8v.2"/>',
    "quote": '<path d="M5 17c0-4 1-7 4.5-9M13 17c0-4 1-7 4.5-9"/><circle cx="7" cy="15.5" r="2"/><circle cx="15" cy="15.5" r="2"/>',
    "stars": '<path d="M12 4l1.7 4.3L18 10l-4.3 1.7L12 16l-1.7-4.3L6 10l4.3-1.7z"/><path d="M18.5 15.5l.8 1.7 1.7.8-1.7.8-.8 1.7-.8-1.7-1.7-.8 1.7-.8z"/>',
    "file": '<path d="M6.5 3.5h7.5l4.5 4.5v12.5h-12z"/><path d="M14 3.5V8h4.5M9 13h6M9 16.5h6"/>',
    "arrow": '<path d="M5 12h14M13 6l6 6-6 6"/>',
    "folder": '<path d="M3.5 7.5a2 2 0 0 1 2-2h4l2 2.5h7a2 2 0 0 1 2 2v7.5a2 2 0 0 1-2 2h-13a2 2 0 0 1-2-2z"/>',
}


def icon(name, cls=""):
    body = ICONS.get(name, ICONS["arrow"])
    return (f'<svg viewBox="0 0 24 24" aria-hidden="true" class="{cls}">'
            f'{body}</svg>')


# --------------------------------------------------------------------------
# cards
# --------------------------------------------------------------------------

class cards_node(nodes.General, nodes.Element):
    pass


class card_node(nodes.General, nodes.Element):
    pass


class Cards(SphinxDirective):
    has_content = True
    option_spec = {"columns": directives.positive_int}

    def run(self):
        node = cards_node()
        node["columns"] = self.options.get("columns", 2)
        self.state.nested_parse(self.content, self.content_offset, node)
        return [node]


class Card(SphinxDirective):
    """``card`` TITLE, with :link: (a document, or a URL) and :icon:."""
    required_arguments = 1
    final_argument_whitespace = True
    has_content = True
    option_spec = {"link": directives.unchanged_required,
                   "icon": directives.unchanged, "go": directives.unchanged}

    def run(self):
        node = card_node()
        node["title"] = self.arguments[0]
        node["link"] = self.options.get("link", "")
        node["icon"] = self.options.get("icon", "")
        node["go"] = self.options.get("go", "")
        node["docname"] = self.env.docname
        if node["link"] and "://" not in node["link"]:
            self.env.note_dependency(self.env.doc2path(self.env.docname))
        self.state.nested_parse(self.content, self.content_offset, node)
        return [node]


def _href(translator, link):
    if not link or "://" in link or link.startswith("#"):
        return link
    import posixpath
    target, _, anchor = link.partition("#")
    builder = translator.builder
    here = builder.current_docname
    # a link is relative to the page it is written in, as in a toctree;
    # one starting with / is from the top of the site
    doc = (target.lstrip("/") if target.startswith("/") else
           posixpath.normpath(posixpath.join(posixpath.dirname(here), target)))
    uri = builder.get_relative_uri(here, doc)
    return uri + ("#" + anchor if anchor else "")


def visit_cards(self, node):
    self.body.append(f'<div class="cards" style="--cols:{node["columns"]}">')


def depart_cards(self, node):
    self.body.append("</div>")


def visit_card(self, node):
    href = _href(self, node["link"])
    self.body.append(f'<a class="card" href="{html.escape(href)}">')
    if node["icon"]:
        self.body.append(f'<span class="card-icon">{icon(node["icon"])}</span>')
    self.body.append(f'<span class="card-title">{html.escape(node["title"])}</span>'
                     f'<span class="card-text">')
    # the description, as inline text: block elements are not allowed in a link
    text = " ".join(p.astext() for p in node.children)
    self.body.append(html.escape(text))
    self.body.append("</span>")
    if node["go"]:
        self.body.append(f'<span class="card-go">{html.escape(node["go"])} →</span>')
    self.body.append("</a>")
    raise nodes.SkipNode


# --------------------------------------------------------------------------
# diagrams and screenshots
# --------------------------------------------------------------------------

_BOX = re.compile(r'<rect class="(d-box|d-accent)"([^>]*?)/>')


def _glass_boxes(svg, key):
    """Give every box of a diagram the glass of the rest of the site: over
    each box, a sheen fading down from its top and a rim brighter along the
    top edge. The gradients' colours are set by tsr.css for each theme."""
    sheen, rim = f"gs-{key}", f"gr-{key}"
    defs = (f'<linearGradient id="{sheen}" x1="0" y1="0" x2="0" y2="1">'
            '<stop class="g-sheen-top" offset="0"/><stop class="g-sheen-end" offset=".55"/></linearGradient>'
            f'<linearGradient id="{rim}" x1="0" y1="0" x2="0" y2="1">'
            '<stop class="g-rim-top" offset="0"/><stop class="g-rim-end" offset="1"/></linearGradient>')
    if "<defs>" in svg:
        svg = svg.replace("<defs>", "<defs>" + defs, 1)
    else:
        svg = re.sub(r"(<svg[^>]*>)", r"\1<defs>" + defs + "</defs>", svg, count=1)
    return _BOX.sub(lambda m: m.group(0) + f'<rect class="d-sheen"{m.group(2)} '
                    f'fill="url(#{sheen})" stroke="url(#{rim})"/>', svg)


class Diagram(SphinxDirective):
    required_arguments = 1
    option_spec = {"caption": directives.unchanged, "alt": directives.unchanged}

    def run(self):
        path = os.path.join(self.env.srcdir, "_diagrams", self.arguments[0])
        self.env.note_dependency(path)
        with open(path, encoding="utf-8") as handle:
            svg = handle.read()
        if svg.startswith("<?xml"):
            svg = svg.split("?>", 1)[1]
        alt = self.options.get("alt", self.options.get("caption", ""))
        svg = svg.replace("<svg ", f'<svg role="img" aria-label="{html.escape(alt)}" ', 1)
        svg = _glass_boxes(svg, os.path.splitext(self.arguments[0])[0])
        out = [nodes.raw("", f'<div class="diagram">{svg}</div>', format="html")]
        if self.options.get("caption"):
            out.append(nodes.raw("", '<p class="diagram-caption">'
                                 f'{html.escape(self.options["caption"])}</p>',
                                 format="html"))
        return out


class Shots(SphinxDirective):
    """``shots`` NAME: _static/gui/light_NAME.png and dark_NAME.png, with
    the content as caption."""
    required_arguments = 1
    has_content = True
    option_spec = {"alt": directives.unchanged, "width": directives.unchanged}

    def run(self):
        name = self.arguments[0]
        depth = self.env.docname.count("/")
        up = "../" * depth
        alt = html.escape(self.options.get("alt", " ".join(self.content)))
        width = self.options.get("width", "100%")
        fig = nodes.container(classes=["shots-figure"])
        fig += nodes.raw("", (
            f'<figure class="shots" style="max-width:{width}">'
            f'<img class="shot-light" src="{up}_static/gui/light_{name}.png" alt="{alt}" loading="lazy">'
            f'<img class="shot-dark" src="{up}_static/gui/dark_{name}.png" alt="{alt}" loading="lazy">'),
            format="html")
        if self.content:
            cap = nodes.container(classes=["figcaption-wrap"])
            self.state.nested_parse(self.content, self.content_offset, cap)
            fig += nodes.raw("", "<figcaption>", format="html")
            fig += cap
            fig += nodes.raw("", "</figcaption>", format="html")
        fig += nodes.raw("", "</figure>", format="html")
        return [fig]


class Themed(SphinxDirective):
    """``themed`` STEM: _static/STEM-light.png and _static/STEM-dark.png (a
    figure in each theme), with the content as caption."""
    required_arguments = 1
    has_content = True
    option_spec = {"alt": directives.unchanged, "width": directives.unchanged}

    def run(self):
        stem = self.arguments[0]
        up = "../" * self.env.docname.count("/")
        alt = html.escape(self.options.get("alt", " ".join(self.content)))
        width = self.options.get("width", "100%")
        fig = nodes.container(classes=["themed-figure"])
        fig += nodes.raw("", (
            f'<figure class="themed" style="max-width:{width}">'
            f'<img class="shot-light" src="{up}_static/{stem}-light.png" alt="{alt}" loading="lazy">'
            f'<img class="shot-dark" src="{up}_static/{stem}-dark.png" alt="{alt}" loading="lazy">'),
            format="html")
        if self.content:
            cap = nodes.container(classes=["figcaption-wrap"])
            self.state.nested_parse(self.content, self.content_offset, cap)
            fig += nodes.raw("", "<figcaption>", format="html")
            fig += cap
            fig += nodes.raw("", "</figcaption>", format="html")
        fig += nodes.raw("", "</figure>", format="html")
        return [fig]


# --------------------------------------------------------------------------
# navigation
# --------------------------------------------------------------------------

def _title(env, doc):
    t = env.titles.get(doc)
    return t.astext() if t is not None else doc


def _nav(app):
    env = app.env
    out = []
    for tab in app.config.tsr_nav:
        groups = []
        for gtitle, docs in tab["groups"]:
            groups.append({"title": gtitle, "pages": [
                {"doc": d, "title": _title(env, d)} for d in docs
                if d in env.found_docs]})
        out.append({"title": tab["title"], "slug": tab["slug"],
                    "index": tab["index"], "prefixes": tab.get("prefixes", []),
                    "groups": groups})
    return out


_TABLE = re.compile(r'(<table class="[^"]*docutils.*?</table>)', re.S)


def _wrap_tables(body):
    """Put each table in its own scrolling box, so a wide table scrolls
    inside the page and its header row always spans the whole table."""
    return _TABLE.sub(r'<div class="table-wrap">\1</div>', body)


def _asset_version(app, name):
    """A short hash of a theme file, added to its URL so that a browser
    fetches the file again whenever it changes."""
    path = os.path.join(app.confdir, "_themes", "tsr", "static", name)
    try:
        with open(path, "rb") as f:
            return hashlib.md5(f.read()).hexdigest()[:8]
    except OSError:
        return "0"


def page_context(app, pagename, templatename, context, doctree):
    context["tsr_js_version"] = _asset_version(app, "tsr.js")
    if "body" in context:
        context["body"] = _wrap_tables(context["body"])
    nav = _nav(app)
    context["nav"] = nav
    context["is_home"] = pagename == app.config.master_doc
    tab, group = None, None
    for t in nav:
        for g in t["groups"]:
            if any(p["doc"] == pagename for p in g["pages"]):
                tab, group = t, g
        if tab:
            break
    if tab is None:
        for t in nav:
            if any(pagename.startswith(p) for p in t["prefixes"]):
                tab = t
                break
    context["tab"] = tab
    crumbs = []
    if tab:
        crumbs.append({"title": tab["title"], "doc": tab["index"]})
        if group:
            crumbs.append({"title": group["title"], "doc": None})
        elif pagename != tab["index"]:
            parent = app.config.tsr_parents.get(
                next((p for p in tab["prefixes"] if pagename.startswith(p)), ""))
            if parent:
                crumbs.append({"title": _title(app.env, parent), "doc": parent})
    context["crumbs"] = crumbs
    # the sidebar entry to highlight: the page, or for a page not listed
    # (a generated API page) the page it belongs under
    active = pagename
    if group is None and tab is not None:
        prefix = next((p for p in tab["prefixes"] if pagename.startswith(p)),
                      "")
        active = app.config.tsr_parents.get(prefix, pagename)
    context["active_doc"] = active
    flat = [p for g in (tab["groups"] if tab else []) for p in g["pages"]]
    idx = next((i for i, p in enumerate(flat) if p["doc"] == pagename), None)
    context["prev_page"] = flat[idx - 1] if idx else None
    context["next_page"] = (flat[idx + 1] if idx is not None
                            and idx + 1 < len(flat) else None)


def no_em_dashes(app, doctree):
    """The house style has no em dashes. Text the pages draw from the code
    (docstrings, the settings' documentation) has some: each becomes a
    comma, or a hyphen between two words written without spaces."""
    for node in list(doctree.traverse(nodes.Text)):
        text = node.astext()
        if "\u2014" not in text:
            continue
        new = text.replace(" \u2014 ", ", ").replace("\u2014 ", ", ")
        new = new.replace(" \u2014", ",").replace("\u2014", "-")
        node.parent.replace(node, nodes.Text(new))


def setup(app):
    app.add_config_value("tsr_nav", [], "html")
    app.add_config_value("tsr_parents", {}, "html")
    app.add_node(cards_node, html=(visit_cards, depart_cards))
    app.add_node(card_node, html=(visit_card, None))
    app.add_directive("cards", Cards)
    app.add_directive("linkcard", Card)
    app.add_directive("diagram", Diagram)
    app.add_directive("shots", Shots)
    app.add_directive("themed", Themed)
    app.connect("html-page-context", page_context)
    app.connect("doctree-read", no_em_dashes)
    return {"version": "1.0", "parallel_read_safe": True,
            "parallel_write_safe": True}
