"""
A night of TripleSpec data, organized: what each frame is, how the frames
pair, and which calibrations and standard each target uses.

This is the bookkeeping the IDL tool leaves to the observer -- reading the
headers, looking at quick-look images, typing frame numbers into fields.
Every decision here only chooses inputs; the reduction steps that consume
them are the Spextool ones, unchanged. Every decision is also recorded in
the night's plan (night.json) and can be edited before anything runs.

Rules, taken from the tutorial and the headers:

  frames      OBSTYPE/OBJECT: Dflat + DFLAT_ON / DFLAT_OFF, Calibration +
              arc, Object. Anything else (focus, tests) is set aside.
  sequences   consecutive Object frames with one OBJECT name.
  nod         offsets of each frame from the first, from RA/DEC. Positions
              within the slit length of each other are an along-slit nod
              (ABBA): a point source, paired frame by frame; a frame that
              repeats the previous position is left out, as the tutorial
              leaves out frame 432. A position farther than the slit is
              sky: an extended source, each object frame paired with the
              nearest sky frame in time (433-434, 435-434 in the tutorial).
  standards   HD/HIP/HR/BD names or exposures of 30 s or less.
  wavecal     a science sequence of long exposures is calibrated on its own
              OH lines ("Sky"); a standard uses the nearest science
              target's; with no long exposures at all, the arcs.
  telluric    for each target, the standard nearest in airmass, then time.
"""

import datetime as _dt
import glob
import json
import os
import re
from dataclasses import dataclass, asdict, field

import numpy as np
from astropy.io import fits

SLIT_LENGTH = 38.0          # arcsec, XD.dat SLITH_ARC
SAME_POSITION = 2.0         # arcsec: two pointings closer than this are one
STANDARD_NAME = re.compile(r"^(HD|HIP|HR|BD|SAO|TYC|GSC|2MASS)[\s_-]*\d", re.I)


def safe_name(name):
    """An OBJECT string as a file-name prefix: 'kks_03-gc1 (two ob' -> 'kks_03-gc1_two_ob'."""
    s = re.sub(r"[^A-Za-z0-9+\-.]+", "_", str(name).strip())
    return re.sub(r"_+", "_", s).strip("_") or "unknown"


def _sexa(text, hours):
    t = str(text).strip()
    neg = t.startswith("-")
    p = [float(v) for v in t.replace("+", "").replace("-", "").split(":")]
    v = p[0] + p[1] / 60. + p[2] / 3600.
    v = -v if neg else v
    return v * 15.0 if hours else v


@dataclass
class Frame:
    path: str
    name: str
    number: int
    prefix: str
    object: str
    obstype: str
    exptime: float
    airmass: float
    ra: float
    dec: float
    time: str
    kind: str = ""              # flat_on, flat_off, arc, object, ignored
    note: str = ""

    @property
    def mjd(self):
        t = _dt.datetime.fromisoformat(self.time)
        return (t - _dt.datetime(1858, 11, 17)).total_seconds() / 86400.


def read_frame(path):
    h = fits.getheader(path)
    base = os.path.basename(path)
    m = re.match(r"^(.*?)(\d+)\.fits(\.fz|\.gz)?$", base)
    prefix, number = (m.group(1), int(m.group(2))) if m else (base, -1)
    try:
        ra = _sexa(h.get("RA", "0:0:0"), True)
        dec = _sexa(h.get("DEC", "0:0:0"), False)
    except ValueError:
        ra = dec = np.nan
    fr = Frame(path=os.path.abspath(path), name=base, number=number, prefix=prefix,
               object=str(h.get("OBJECT", "")).strip(), obstype=str(h.get("OBSTYPE", "")).strip(),
               exptime=float(h.get("EXPTIME", 0) or 0), airmass=float(h.get("AIRMASS", 0) or 0),
               ra=ra, dec=dec, time=str(h.get("DATE-OBS", "1970-01-01T00:00:00")))
    # A raw frame is the file the instrument wrote: its FILENAME is its own
    # name. A renamed copy or a derived image (t20-21.fits, an A-B made at
    # the telescope, still carries the first frame's FILENAME) is not.
    recorded = str(h.get("FILENAME", base)).strip()
    if recorded and recorded != base and recorded + ".fz" != base and recorded + ".gz" != base:
        fr.kind = "ignored"
        fr.note = f"not a raw frame: its header names it {recorded}"
    return fr


def classify(fr):
    o, t = fr.object.upper(), fr.obstype.upper()
    if "DFLAT_ON" in o or (("FLAT" in t) and "ON" in o):
        return "flat_on"
    if "DFLAT_OFF" in o or (("FLAT" in t) and "OFF" in o):
        return "flat_off"
    if t.startswith("CALIB") and ("ARC" in o or "LAMP" in o):
        return "arc"
    if t == "OBJECT":
        return "object"
    return "ignored"


def offset_arcsec(a, b):
    """Angular separation of two frames' pointings, in arcsec."""
    d = np.radians
    c = (np.sin(d(a.dec)) * np.sin(d(b.dec)) +
         np.cos(d(a.dec)) * np.cos(d(b.dec)) * np.cos(d(a.ra - b.ra)))
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600.)


@dataclass
class Sequence:
    name: str                   # the file-name prefix for every product
    object: str
    frames: list                # Frame.name, in time order
    role: str                   # science | standard
    source: str                 # point | extended
    pattern: str                # ABBA | ABAB | A-B | object-sky | single
    pairs: list                 # [[a, b], ...]  a minus b
    beams: dict                 # frame -> beam label
    skipped: list = field(default_factory=list)
    airmass: float = 0.0
    mjd: float = 0.0
    ra: float = 0.0
    dec: float = 0.0
    exptime: float = 0.0
    wavecal_from: str = ""      # a sequence name, or 'arcs'
    standard: str = ""          # for science: the standard sequence
    bmag: float = float("nan")
    vmag: float = float("nan")
    sptype: str = ""
    notes: list = field(default_factory=list)


def _cluster(frames):
    """Pointing clusters: the label of each frame and the cluster centres."""
    labels, centres = [], []
    for f in frames:
        for k, c in enumerate(centres):
            if offset_arcsec(f, c) < SAME_POSITION:
                labels.append(k)
                break
        else:
            centres.append(f)
            labels.append(len(centres) - 1)
    return labels, centres


def analyse_sequence(frames):
    """Nod pattern, source type and frame pairs for one target's frames."""
    labels, centres = _cluster(frames)
    names = [f.name for f in frames]
    counts = np.bincount(labels)
    main = int(np.argmax(counts))
    far = [k for k, c in enumerate(centres) if offset_arcsec(c, centres[main]) > SLIT_LENGTH]
    beams = {}
    pairs, skipped, notes = [], [], []
    if len(frames) == 1:
        return "single", "point", [], {names[0]: "A"}, [], ["A single frame: reduce in A mode."]
    if far:
        sky = [i for i, l in enumerate(labels) if l in far]
        obj = [i for i, l in enumerate(labels) if l not in far]
        for i in sky:
            beams[names[i]] = "sky"
        for i in obj:
            beams[names[i]] = "object" if labels[i] == main else f"object{labels[i]}"
            j = min(sky, key=lambda s: (abs(frames[s].mjd - frames[i].mjd), s))
            pairs.append([names[i], names[j]])
        return "object-sky", "extended", pairs, beams, [], notes
    # An along-slit nod has two positions; pointing jitter can split one of
    # them into small clusters, so each frame goes to the nearer of the two
    # most used positions.
    top = sorted(range(len(centres)), key=lambda k: -counts[k])[:2]
    letters = {}
    for i, f in enumerate(frames):
        l = min(top, key=lambda k: offset_arcsec(f, centres[k]))
        labels[i] = l
        if l not in letters:
            letters[l] = "AB"[len(letters)]
        beams[names[i]] = letters[l]
    i = 0
    while i < len(frames) - 1:
        if labels[i] != labels[i + 1]:
            pairs.append([names[i], names[i + 1]])
            i += 2
        else:
            skipped.append(names[i])
            notes.append(f"{names[i]} is at the same position as {names[i + 1]}; left out, "
                         "as frame 432 is in the tutorial.")
            i += 1
    if i == len(frames) - 1:
        skipped.append(names[i])
        notes.append(f"{names[i]} has no partner; left out.")
    # The pattern from the pairs themselves: ABBA when every second pair is
    # B-A (the case xspextool's 'ABBA to ABAB' switch is for), ABAB when all
    # are A-B (the switch must then stay off).
    first = [beams[a] == beams[pairs[0][0]] for a, _ in pairs] if pairs else []
    if all(first):
        pattern = "ABAB"
    elif all(f == (k % 2 == 0) for k, f in enumerate(first)):
        pattern = "ABBA"
    else:
        pattern = "A-B"
        notes.append("The pairs do not alternate regularly, so combining their images first would "
                     "cancel some of them; extract them pair by pair.")
    return pattern, "point", pairs, beams, skipped, notes


@dataclass
class Night:
    raw_dir: str
    out_dir: str
    frames: list
    flats_on: list
    flats_off: list
    arcs: list
    ignored: list
    sequences: list
    warnings: list = field(default_factory=list)
    xs_background: bool = False         # the IDL default: see extract.py
    flat_name: str = "flat"             # xspextool's Output Flat Name: cal/<flat_name>.fits

    # -------------------------------------------------- names and places
    @property
    def cal_dir(self):
        return os.path.join(self.out_dir, "cal")

    @property
    def proc_dir(self):
        return os.path.join(self.out_dir, "proc")

    @property
    def qa_dir(self):
        return os.path.join(self.out_dir, "qa")

    def frame(self, name):
        return next(f for f in self.frames if f.name == name)

    def seq(self, name):
        return next(s for s in self.sequences if s.name == name)

    def path(self, source, step, ext="fits"):
        """Every product is <source>_<step>.<ext>; calibrations go in cal/."""
        d = self.cal_dir if step in ("flat", "wavecal") else self.proc_dir
        return os.path.join(d, f"{source}_{step}.{ext}")

    def flat_path(self):
        from triplespec_reducer.flat import flat_output_path
        return flat_output_path(self.cal_dir, self.flat_name)

    def wavecal_path(self, source):
        return os.path.join(self.cal_dir, f"{source}_wavecal.fits")

    def spectrum_path(self, source, frame):
        n = self.frame(frame).number
        return os.path.join(self.proc_dir, f"{source}_spec{n:04d}.fits")

    def products(self, source):
        """Everything already made for `source`, by step: the 'give a name, get the files' lookup."""
        found = {}
        longer = [s.name + "_" for s in self.sequences
                  if s.name != source and s.name.startswith(source + "_")]
        for d in (self.cal_dir, self.proc_dir, self.qa_dir):
            for p in sorted(glob.glob(os.path.join(d, f"{source}_*"))):
                if any(os.path.basename(p).startswith(x) for x in longer):
                    continue
                step = os.path.basename(p)[len(source) + 1:].split(".")[0]
                found.setdefault(step, p)
        return found

    # -------------------------------------------------- persistence
    def save(self, path=None):
        path = path or os.path.join(self.out_dir, "night.json")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        d = asdict(self)
        with open(path, "w") as fh:
            json.dump(d, fh, indent=1, default=str)
        return path

    @classmethod
    def load(cls, path):
        with open(path) as fh:
            d = json.load(fh)
        d["frames"] = [Frame(**f) for f in d["frames"]]
        d["sequences"] = [Sequence(**s) for s in d["sequences"]]
        return cls(**d)


def lookup_standard(name, cache_path=None):
    """B, V and spectral type from SIMBAD, remembered between runs."""
    cache_path = cache_path or os.path.join(os.path.expanduser("~"), ".triplespec_reducer",
                                            "simbad.json")
    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path) as fh:
            cache = json.load(fh)
    key = re.sub(r"\s+", "", name.upper())
    if key in cache:
        return cache[key]
    try:
        import warnings
        from astroquery.simbad import Simbad
        s = Simbad()
        s.TIMEOUT = 20
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s.add_votable_fields("B", "V", "sp_type")
            t = s.query_object(re.sub(r"^(HD|HIP|HR)[_-]?", r"\1 ", name, flags=re.I))
        if t is None or len(t) == 0:
            return None
        r = dict(B=float(t["B"][0]), V=float(t["V"][0]), sptype=str(t["sp_type"][0]),
                 simbad=str(t["main_id"][0]))
    except Exception:
        return None
    cache[key] = r
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "w") as fh:
        json.dump(cache, fh, indent=1)
    return r


def flat_sanity(on, off):
    """The quick-look check of the tutorial: lamps-on frames are much brighter."""
    def level(f):
        d = fits.getdata(f.path).astype(float)
        return float(np.median(d[300:700, 200:900]))
    if not on:
        return ["No lamps-on flats were found; the flat step cannot run."]
    warn = []
    lon = level(on[len(on) // 2])
    if off:
        loff = level(off[len(off) // 2])
        if not lon > 5 * max(loff, 1.0):
            warn.append(f"Lamps-on and lamps-off flats look alike (median {lon:.0f} vs {loff:.0f} DN); "
                        "check their labels.")
    return warn


def organize(raw_dir, out_dir, lookup=True, progress=None):
    """Read every frame in raw_dir and propose the whole night's plan."""
    say = progress or (lambda *_: None)
    from triplespec_reducer.frames import find_raw_frames
    from triplespec_reducer.instrument import read_instrfile
    paths = find_raw_frames(raw_dir, read_instrfile().suffix or "*.fits*")
    frames = []
    for p in paths:
        try:
            f = read_frame(p)
        except (OSError, KeyError, ValueError):
            continue
        if f.kind != "ignored":
            f.kind = classify(f)
        frames.append(f)
    frames.sort(key=lambda f: (f.time, f.name))
    say(f"{len(frames)} frames read")
    on = [f for f in frames if f.kind == "flat_on"]
    off = [f for f in frames if f.kind == "flat_off"]
    # The tutorial's flats share one file-name prefix; a lone lamps-on frame
    # with another prefix (SPEC_141Tau_0001 here) is set aside with a note.
    if on:
        main = max(set(f.prefix for f in on), key=lambda p: sum(f.prefix == p for f in on))
        for f in on:
            if f.prefix != main:
                f.kind, f.note = "ignored", "a lamps-on flat with a different file prefix"
        on = [f for f in on if f.kind == "flat_on"]
    arcs = [f for f in frames if f.kind == "arc"]
    objs = [f for f in frames if f.kind == "object"]
    ignored = [f for f in frames if f.kind == "ignored"]
    warnings = flat_sanity(on, off)

    groups = []
    for f in objs:
        if groups and groups[-1][-1].object == f.object:
            groups[-1].append(f)
        else:
            groups.append([f])
    used = {}
    seqs = []
    for g in groups:
        base = safe_name(g[0].object)
        used[base] = used.get(base, 0) + 1
        seqs.append((base, g))
    totals = {b: sum(1 for x, _ in seqs if x == b) for b, _ in seqs}
    count = {}
    out = []
    for base, g in seqs:
        count[base] = count.get(base, 0) + 1
        name = base if totals[base] == 1 else f"{base}_{count[base]}"
        pattern, source, pairs, beams, skipped, notes = analyse_sequence(g)
        std = bool(STANDARD_NAME.match(g[0].object)) or max(f.exptime for f in g) <= 30
        # Object frames at a second pointing (an offset along a galaxy, say)
        # are another spectrum of another place: they become their own
        # target rather than being averaged into the first.
        groups_by_beam = {}
        for a, b in pairs:
            groups_by_beam.setdefault(beams.get(a, "A") if pattern == "object-sky" else "main", []).append([a, b])
        for k, (beam, bp) in enumerate(sorted(groups_by_beam.items(), key=lambda kv: kv[0] != "object")):
            sub = name if k == 0 else f"{name}_pos{k + 1}"
            fr = [f for f in g if any(f.name in p for p in bp)] if pattern == "object-sky" else g
            first_obj = next(f for f in fr if beams.get(f.name) != "sky")
            main_obj = next(f for f in g if beams.get(f.name) == "object") if pattern == "object-sky" else g[0]
            subnotes = list(notes) if k == 0 else [
                f"Its object frames point {offset_arcsec(first_obj, main_obj):.1f} arcsec from "
                f"{name}'s; kept apart as a separate target."]
            out.append(Sequence(name=sub, object=g[0].object, frames=[f.name for f in fr],
                                role="standard" if std else "science", source="point" if std else source,
                                pattern=pattern, pairs=bp, beams={f.name: beams[f.name] for f in fr},
                                skipped=skipped if k == 0 else [],
                                airmass=float(np.mean([f.airmass for f in fr])),
                                mjd=float(np.mean([f.mjd for f in fr])),
                                ra=float(np.mean([f.ra for f in fr])), dec=float(np.mean([f.dec for f in fr])),
                                exptime=float(np.max([f.exptime for f in fr])), notes=subnotes))
        if not pairs:
            out.append(Sequence(name=name, object=g[0].object, frames=[f.name for f in g],
                                role="standard" if std else "science", source="point" if std else source,
                                pattern=pattern, pairs=pairs, beams=beams, skipped=skipped,
                                airmass=float(np.mean([f.airmass for f in g])),
                                mjd=float(np.mean([f.mjd for f in g])),
                                ra=float(np.mean([f.ra for f in g])), dec=float(np.mean([f.dec for f in g])),
                                exptime=float(np.max([f.exptime for f in g])), notes=notes))
    sci = [s for s in out if s.role == "science"]
    stds = [s for s in out if s.role == "standard"]
    # Science targets first, so a standard can take the solution its
    # nearest target actually uses.
    for s in out:
        if s.role == "science" and s.exptime >= 60:
            # An offset pointing shares its parent's sky lines.
            parent = s.name.split("_pos")[0]
            s.wavecal_from = parent if "_pos" in s.name and any(x.name == parent for x in out) else s.name
    for s in out:
        if s.wavecal_from:
            continue
        near = min((x for x in sci if x.exptime >= 60), key=lambda x: abs(x.mjd - s.mjd), default=None)
        s.wavecal_from = near.wavecal_from if near else "arcs"
    for s in sci:
        if stds:
            # Within xtellcor's 0.1 airmass tolerance the nearest in time,
            # otherwise the nearest in airmass.
            ok = [t for t in stds if abs(t.airmass - s.airmass) <= 0.1]
            best = (min(ok, key=lambda t: abs(t.mjd - s.mjd)) if ok else
                    min(stds, key=lambda t: abs(t.airmass - s.airmass)))
            s.standard = best.name
            if abs(best.airmass - s.airmass) > 0.1:
                s.notes.append(f"The closest standard in airmass ({best.name}) differs by "
                               f"{abs(best.airmass - s.airmass):.2f}; xtellcor warns above 0.1.")
    if lookup:
        for s in stds:
            r = lookup_standard(s.object)
            if r:
                s.bmag, s.vmag, s.sptype = r["B"], r["V"], r["sptype"]
                if not s.sptype.upper().startswith("A0"):
                    s.notes.append(f"SIMBAD lists {s.object} as {s.sptype}, not A0 V; the Vega "
                                   "model assumes an A0 V star.")
            else:
                s.notes.append("B and V were not found on SIMBAD; enter them before the telluric step.")
    night = Night(raw_dir=os.path.abspath(raw_dir), out_dir=os.path.abspath(out_dir), frames=frames,
                  flats_on=[f.name for f in on], flats_off=[f.name for f in off],
                  arcs=[f.name for f in arcs], ignored=[f.name for f in ignored],
                  sequences=out, warnings=warnings)
    return night
