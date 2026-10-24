"""
The reduction of a night, step by step, on the plan night.organize() made.

    flat       cal/flat.fits                      TS41CALS_MKFLAT
    wavecal    cal/<source>_wavecal.fits          TS41CALS_MKWAVECAL (Sky; CuHeAr for 'arcs')
    extract    proc/<source>_specNNNN.fits         XSPEXTOOL, A-B, every pair
    combine    proc/<source>_comb.fits             XCOMBSPEC
    telluric   proc/<source>_tellspec.fits         XTELLCOR (the standard's file)
               proc/<source>_tc.fits               (the corrected target)
    merge      proc/<source>_merged.fits           XMERGEORDERS
    clean      proc/<source>_clean.fits            XCLEANSPEC (only when asked)

Every setting defaults to the IDL tool's own default. Where the IDL asks the
user to pick something -- an aperture position, a continuum window, a
region to measure the telluric shift on -- a choice is made here, written
into <source>_log.json with its value, and can be overridden.
"""

import json
import os
import time
from dataclasses import dataclass, field, asdict

import numpy as np

from triplespec_reducer import (idl as I, flat as F, wavecal as W, extract as E, combine as C,
                                telluric as T, merge as M, clean as CL)
from triplespec_reducer.instrument import read_instrfile, read_atran, read_modefile
from triplespec_reducer.night import Night


@dataclass
class Options:
    """Every knob, at the IDL default unless noted."""
    # xspextool
    flatfield: bool = True
    errorprop: bool = True
    ybuffer: int = 1
    atmosthresh: float = 0.4
    psnaps: int = 2                     # an A-B pair or combined pairs: a positive and a negative beam
    psnaps_single: int = 1              # A and A-Sky: the source appears once
    findmode: str = "auto"              # auto | guess | fixed
    apradius: float = 1.0
    psbgstart: float = 2.5
    psbgwidth: float = 5.0
    psbgdeg: int = 1
    tracestep: int = 5
    tracesumap: int = 5
    tracewinthresh: float = 3.0
    tracesigthresh: float = 0.1
    tracedeg: int = 1
    xs_positions: list = None           # None: the profile peak of each pair
    xs_radii: list = None
    xs_bg: str = "0-2,13-15"            # ts4_soar.dat XSBG
    xs_bgdeg: int = 0
    xs_background: bool = False         # IDL default: not subtracted (see extract.py)
    wavecal_offset: float = None        # None: automatic cross-correlation (the default)
    mode: str = "A-B"                   # A-B | A | A-Sky
    combine_images: bool = False        # XSPEXTOOL_COMBINEIMGS before extracting
    combimg_stat: str = "Robust Mean (Std Error)"   # ts4_soar.dat COMBSTAT
    combimg_thresh: float = 8.0
    abba2abab: bool = True              # so the combined pairs add up
    combimg_subtract_sky: bool = False
    sky_frames: list = None             # A-Sky: frames combined into the sky
    signs: list = None                  # override the aperture signs
    store_trace: bool = False
    import_trace: str = None
    # xcombspec
    scaleorder: int = 5                 # the tutorial: "usually order 4 or 5"
    statistic: str = "Robust Weighted Mean"
    combthresh: float = 8.0
    correct_shapes: bool = False
    combine_apertures: bool = False
    # xtellcor
    kernel_order: int = 6               # Pa delta, as the TS4 README advises
    normreg: list = field(default_factory=lambda: [[0.990, 0.997], [1.013, 1.020]])
    normdegree: int = 3                 # the tutorial's choice; XMC_CONKERN's own field starts at 1
    linereg: list = field(default_factory=lambda: [0.997, 1.013])
    tcmethod: str = "Deconvolution"
    propstdunc: bool = True
    units: str = "ergs s-1 cm-2 A-1"
    tension: float = 10.
    vshift: float = None                # override the measured Vega velocity
    shift_regions: dict = None          # order -> [w0, w1]; None: chosen from ATRAN
    find_shift: bool = True
    # xmergeorders
    anchor: int = 3
    merge_scales: dict = None
    merge_trims: dict = None
    merge_weighted: bool = True
    merge_sum: bool = False
    # xcleanspec
    clean: bool = False
    clean_method: str = "Savitzky-Golay"
    workers: int = 0                    # parallel processes for wavecal and extraction; 0: automatic
    sgwin: float = None
    sgdeg: int = 2
    fwhm: float = None

    @classmethod
    def tutorial(cls):
        """The choices made in the SOAR TripleSpec tutorial (Dahmer-Hahn)."""
        return cls(normdegree=3, xs_bg="2-6,32-36", scaleorder=5)

    @classmethod
    def idl_defaults(cls):
        """Exactly what the IDL panels start with, where that differs from the above."""
        return cls(normdegree=1)


def _bg_regions(text):
    out = []
    for part in str(text).split(","):
        a, b = part.split("-")
        out.append((float(a), float(b)))
    return out


def E_reduce(files, *a, **k):
    from triplespec_reducer.extract import reduce_pair
    return reduce_pair(files or [], *a, **k)


class Reducer:
    """Runs the steps for a Night, writing every product and a log per source."""

    def __init__(self, night, options=None, progress=None, force=False):
        self.n = night
        self.o = options or Options()
        self.say = progress or (lambda *a: None)
        self.force = force
        os.makedirs(night.cal_dir, exist_ok=True)
        os.makedirs(night.proc_dir, exist_ok=True)
        os.makedirs(night.qa_dir, exist_ok=True)
        self._flat = None
        self._wc = {}
        self._atran = None

    # ------------------------------------------------------------- helpers
    def log(self, source, step, **entries):
        path = self.n.path(source, "log", "json")
        d = {}
        if os.path.exists(path):
            with open(path) as fh:
                d = json.load(fh)
        d[step] = dict(entries, time=time.strftime("%Y-%m-%d %H:%M:%S"))
        with open(path, "w") as fh:
            json.dump(d, fh, indent=1, default=lambda v: v.tolist() if hasattr(v, "tolist") else str(v))

    def _raw(self, names):
        return [self.n.frame(x).path for x in names]

    def flatinfo(self):
        if self._flat is None:
            self._flat = F.read_flat(self.n.flat_path())
        return self._flat

    def atran(self):
        if self._atran is None:
            self._atran = read_atran(self.flatinfo()["rp"])
        return self._atran

    def wcinfo(self, source):
        if source not in self._wc:
            self._wc[source] = W.read_wavecal(self.n.wavecal_path(source),
                                              rotation=self.flatinfo()["rotation"])
        return self._wc[source]

    # ------------------------------------------------------------- steps
    def flat(self):
        out = self.n.flat_path()
        if os.path.exists(out) and not self.force:
            return out
        self.say("flat", "Building the flat")
        from astropy.io import fits as _f
        from triplespec_reducer.instrument import INSTRUMENT_DIR
        bpm = _f.getdata(os.path.join(INSTRUMENT_DIR, "fakebmask.fits"))
        r = F.make_flat(self._raw(self.n.flats_on), self._raw(self.n.flats_off) or None,
                        badpixmask=bpm, progress=lambda s: self.say("flat", s))
        F.write_flat(out, r)
        self._flat = None
        self.log("night", "flat", on=self.n.flats_on, off=self.n.flats_off,
                 rms=r["rms"], xranges=r["xranges"])
        return out

    def wavecal(self, source):
        out = self.n.wavecal_path(source)
        if os.path.exists(out) and not self.force:
            return out
        if source == "arcs":
            files, kind = self._raw(self.n.arcs), "CuHeAr Arc"
        else:
            files, kind = self._raw(self.n.seq(source).frames), "Sky"
        self.say("wavecal", f"Calibrating {source} on {kind}")
        r = W.make_wavecal(files, self.flatinfo(), arctype=kind,
                           progress=lambda s: self.say("wavecal", s), offset=self.o.wavecal_offset)
        W.write_wavecal(out, r)
        self._wc.pop(source, None)
        self.last_log = dict(type=kind, frames=[os.path.basename(f) for f in files],
                             offset=r["offset"], rms_angstrom=r["sol1d"]["rms"],
                             nlines=r["sol1d"]["nlines"],
                             stats=[{k: v for k, v in s.items() if k in ("order", "nlines", "spatial_rms",
                                                                         "wave_rms_pix")} for s in r["stats"]])
        self.log(source, "wavecal", **self.last_log)
        return out

    def _extract_kwargs(self, source, wcname, single=False):
        o = self.o
        return dict(atran=self.atran(), apradius=o.apradius,
                    psbg=(o.psbgstart, o.psbgwidth, o.psbgdeg),
                    npeaks=o.psnaps_single if single else o.psnaps,
                    findmode=o.findmode, tracestep=o.tracestep, sumap=o.tracesumap,
                    winthresh=o.tracewinthresh, sigthresh=o.tracesigthresh,
                    tracedeg=o.tracedeg, ybuffer=o.ybuffer, atmosthresh=o.atmosthresh,
                    flatname=os.path.basename(self.n.flat_path()),
                    wavecalname=os.path.basename(self.n.wavecal_path(wcname)),
                    signs=o.signs, flatfield=o.flatfield, errorprop=o.errorprop,
                    trace_in=o.import_trace,
                    trace_out=self.n.path(source, "trace", "dat") if o.store_trace else None)

    def _xs_kwargs(self, files, mode, combined=None, sky=None):
        """The extended-source settings; positions from the profile peak unless typed."""
        from triplespec_reducer import extract as E
        o = self.o
        pos = o.xs_positions
        if pos is None:
            fl, wc = self.flatinfo(), self._wc_current
            work, _, _, _ = E.load_image(files, fl, mode=mode, combined=combined, sky=sky)
            prof = E.mkspatprof2d(work, fl["omask"], wc["wavecal"], wc["spatcal"],
                                  fl["slith_arc"] / fl["slith_pix"], o.ybuffer,
                                  *self.atran(), o.atmosthresh)
            p, _ = E.findpeaks(prof, 1, mode="auto")
            pos = [float(np.round(np.nanmedian(p[:, 0]), 1))]
        return dict(xs_positions=pos, xs_radii=o.xs_radii or [o.apradius] * len(pos),
                    xsbg=_bg_regions(o.xs_bg), xsbgdeg=o.xs_bgdeg, xs_background=o.xs_background)

    def units(self, source):
        """What one extraction is made from: pairs, single frames, or one combined image."""
        s = self.n.seq(source)
        o = self.o
        if o.combine_images:
            return [("combined", [f for pr in s.pairs for f in pr])]
        if o.mode == "A-Sky":
            return [("A-Sky", [f]) for f, b in s.beams.items() if b != "sky"]
        if o.mode == "A":
            return [("A", [f]) for f in s.frames]
        if s.pattern == "single":
            return [("A", [s.frames[0]])]
        return [("A-B", list(p)) for p in s.pairs]

    def _sky(self, source):
        from triplespec_reducer import combineimgs as CI
        s = self.n.seq(source)
        o = self.o
        skyframes = o.sky_frames or [f for f, b in s.beams.items() if b == "sky"]
        sky = os.path.join(self.n.cal_dir, f"{source}_sky.fits")
        if not os.path.exists(sky) or self.force:
            CI.combine_images(self._raw(skyframes), sky, mode="A",
                              statistic=o.combimg_stat, thresh=o.combimg_thresh)
        return sky

    def extract_unit(self, source, k):
        """One unit of units(source): extract it, write its spectra, return a log entry."""
        from triplespec_reducer import combineimgs as CI
        s = self.n.seq(source)
        o = self.o
        fl = self.flatinfo()
        wcname = s.wavecal_from
        wc = self._wc_current = self.wcinfo(wcname)
        mode, frames = self.units(source)[k]
        if mode == "combined":
            out = self.n.path(source, "specall")
            if os.path.exists(out) and not self.force:
                return None, [out]
            cimg = self.n.path(source, "combimg")
            self.say("extract", f"{source}: combining {len(frames)} frames into one image")
            CI.combine_images(self._raw(frames), cimg, mode="A-B",
                              statistic=o.combimg_stat, thresh=o.combimg_thresh,
                              abba2abab=o.abba2abab and s.pattern == "ABBA",
                              subtract_sky=o.combimg_subtract_sky, flatinfo=fl)
            kw = self._extract_kwargs(source, wcname)
            if s.source == "extended":
                kw.update(self._xs_kwargs([], "A", combined=cimg))
            r = E_reduce(None, fl, wc, [out], source=s.source, mode="A", combined=cimg, **kw)
            return (dict(combined=os.path.basename(cimg), frames=frames,
                         positions_found=np.round(r["positions"], 3).tolist(),
                         signs=r["signs"].tolist(), written=[os.path.basename(out)]), r["written"])
        sky = self._sky(source) if mode == "A-Sky" else None
        ops = [self.n.spectrum_path(source, f) for f in frames] + [None]
        if os.path.exists(ops[0]) and not self.force:
            return None, [p for p in ops[:2] if p and os.path.exists(p)]
        self.say("extract", f"{source}: {mode} {k + 1} of {len(self.units(source))}")
        kw = self._extract_kwargs(source, wcname, single=mode != "A-B")
        if s.source == "extended":
            kw.update(self._xs_kwargs(self._raw(frames), mode, sky=sky))
        r = E_reduce(self._raw(frames), fl, wc, ops, source=s.source, mode=mode, sky=sky, **kw)
        return (dict(frames=frames, mode=mode, positions=kw.get("xs_positions"),
                     positions_found=np.round(r["positions"], 3).tolist(),
                     signs=r["signs"].tolist(),
                     written=[os.path.basename(p) for p in r["written"]]), r["written"])

    def _log_extract(self, source, entries):
        s = self.n.seq(source)
        o = self.o
        if entries:
            self.log(source, "extract", source_type=s.source, pattern=s.pattern, wavecal=s.wavecal_from,
                     mode="combined image" if o.combine_images else o.mode, pairs=entries,
                     xs_background=o.xs_background if s.source == "extended" else None)

    def extract(self, source):
        """
        Every pair (A-B), every frame (A, A-Sky), or one combined image,
        through xspextool's extraction. Returns the spectra written.
        """
        s = self.n.seq(source)
        if not os.path.exists(self.n.wavecal_path(s.wavecal_from)):
            self.wavecal(s.wavecal_from)
        written, entries = [], []
        for k in range(len(self.units(source))):
            e, w = self.extract_unit(source, k)
            written += w
            if e:
                entries.append(e)
        self._log_extract(source, entries)
        return written

    def spectra(self, source):
        s = self.n.seq(source)
        if self.o.combine_images and os.path.exists(self.n.path(source, "specall")):
            return [self.n.path(source, "specall")]
        out = []
        for f in s.frames:
            p = self.n.spectrum_path(source, f)
            if os.path.exists(p):
                out.append(p)
        return out

    def combine(self, source):
        out = self.n.path(source, "comb")
        if os.path.exists(out) and not self.force:
            return out
        files = self.spectra(source)
        if not files:
            files = self.extract(source)
        self.say("combine", f"{source}: combining {len(files)} spectra")
        L = C.load(files, combine_apertures=self.o.combine_apertures)
        scales = {}
        for ap in range(1, L["naps"] + 1):
            sc, mr = C.scale(L, self.o.scaleorder, aperture=ap)
            scales[ap] = sc
            if self.o.correct_shapes:
                C.correct_shapes(L, aperture=ap)
        C.combine(L, statistic=self.o.statistic, thresh=self.o.combthresh)
        C.write(out, L)
        notes = []
        for ap, sc in scales.items():
            med = np.median(sc)
            for f, v in zip(files, sc):
                if v <= 0 or v > 3 * med or v < med / 3:
                    notes.append(f"{os.path.basename(f)} needed a scale of {v:.3g} (median {med:.3g}): "
                                 "look at it, and leave it out if it is not like the others.")
        self.log(source, "combine", files=[os.path.basename(f) for f in files],
                 scaleorder=self.o.scaleorder, scales={k: v.tolist() for k, v in scales.items()},
                 statistic=self.o.statistic, thresh=self.o.combthresh, notes=notes)
        return out

    def shift_regions(self, st):
        """
        A region per order to measure the telluric shift on, from ATRAN.

        The IDL user drags across a stretch where the corrected spectrum
        wiggles. Here it is the 0.03 um stretch of the order where the
        atmosphere has the most structure without being opaque, tried in
        order of that structure until the best shift is not at the edge of
        XMC_FINDSHIFT's search.
        """
        aw, at = self.atran()
        cands = {}
        for n, o in enumerate(st["objorders"]):
            w = st["obj"]["spec"][n * st["objnaps"], 0]
            w = w[np.isfinite(w)]
            lo, hi = w.min() + 0.01, w.max() - 0.01
            step = 0.005
            res = []
            for a in np.arange(lo, hi - 0.03, step):
                z = (aw >= a) & (aw <= a + 0.03)
                if z.sum() < 10:
                    continue
                t = at[z]
                if np.min(t) < 0.2:
                    continue
                res.append((float(np.std(t)), [round(float(a), 4), round(float(a + 0.03), 4)]))
            res.sort(key=lambda r: -r[0])
            cands[o] = [r[1] for r in res[:12]]
        return cands

    def telluric(self, source):
        s = self.n.seq(source)
        out = self.n.path(source, "tc")
        tell = self.n.path(source, "tellspec")
        if os.path.exists(out) and not self.force:
            return out
        if not s.standard:
            raise ValueError(f"{source} has no telluric standard assigned.")
        std = self.n.seq(s.standard)
        if not (np.isfinite(std.bmag) and np.isfinite(std.vmag)):
            raise ValueError(f"Enter B and V for {std.object} before the telluric step.")
        sp, tp = self.combine(std.name), self.combine(source)
        o = self.o
        self.say("telluric", f"{source}: with {std.name}")
        st = T.load(sp, tp, std.bmag, std.vmag)
        T.conkernel(st, order=o.kernel_order, normreg=o.normreg, normdegree=o.normdegree,
                    linereg=o.linereg, method=o.tcmethod)
        if o.vshift is not None:
            st["vshift"] = np.float32(o.vshift)                # FLOAT fields in XTELLCOR
        st["tension"] = np.float32(o.tension)
        T.initial_points(st)
        T.contellspec(st, units=o.units)
        chosen, notes = {}, []
        if o.find_shift:
            cands = {k: [v] for k, v in (o.shift_regions or {}).items()} or self.shift_regions(st)
            for n, order in enumerate(st["objorders"]):
                for ap in range(1, st["objnaps"] + 1):
                    done = False
                    for reg in cands.get(order, []):
                        try:
                            sh, _ = T.getshift(st, order, ap, srange=reg)
                            chosen[f"{order}/{ap}"] = dict(region=reg, shift=sh)
                            done = True
                            break
                        except ValueError:
                            st["shift"][n, ap - 1] = 0.0
                    if not done:
                        notes.append(f"order {order}: no region gave a shift inside +-1.5 px; 0 used")
        T.corrected(st, propstdunc=o.propstdunc)
        T.write(out, st, tellpath=tell, std_name=std.object)
        self.log(source, "telluric", standard=std.name, bmag=std.bmag, vmag=std.vmag,
                 sptype=std.sptype, dairmass=st["dairmass"], vshift=st["vshift"],
                 scale=st["scale"], maxdev=st["maxdev"], rmsdev=st["rmsdev"],
                 kernel_order=o.kernel_order, normreg=o.normreg, normdegree=o.normdegree,
                 linereg=o.linereg, shifts=chosen, notes=notes)
        return out

    def merge(self, source):
        out = self.n.path(source, "merged")
        if os.path.exists(out) and not self.force:
            return out
        tp = self.telluric(source)
        from astropy.io import fits as _f
        with _f.open(tp) as h:
            d = h[0].data
            spec = d.astype(d.dtype.newbyteorder("="))     # as MC_READSPEC reads it
            hdr = h[0].header
        orders = [int(t) for t in str(hdr["ORDERS"]).split(",")]
        naps = int(hdr["NAPS"])
        results = []
        for ap in range(1, naps + 1):
            m, used, done = M.merge_orders(spec, orders, naps=naps, aperture=ap,
                                           anchor=self.o.anchor, scales=self.o.merge_scales,
                                           trims=self.o.merge_trims, weighted=self.o.merge_weighted,
                                           sum_=self.o.merge_sum)
            p = out if naps == 1 else self.n.path(source, f"merged_ap{ap}")
            M.write(p, m, hdr, used, done, aperture=ap, source=os.path.basename(tp))
            results.append(p)
            self.log(source, "merge" if naps == 1 else f"merge_ap{ap}", anchor=self.o.anchor,
                     sequence=done, scales=used)
        return results[0]

    def clean(self, source):
        out = self.n.path(source, "clean")
        mp = self.merge(source)
        cube, hdr, two_d = CL.read(mp)
        sm, hist = CL.smooth(cube, method=self.o.clean_method, sgwin=self.o.sgwin,
                             sgdeg=self.o.sgdeg, fwhm=self.o.fwhm,
                             slitw_pix=I.fxpar(hdr, "SLTW_PIX"))
        CL.write(out, sm, hdr, hist, two_d)
        self.log(source, "clean", method=self.o.clean_method, history=hist)
        return out

    # ------------------------------------------------------------- everything
    def _pool(self):
        n = self.o.workers or max(1, min(16, (os.cpu_count() or 2) - 2))
        return n

    def _fail(self, source, step, message):
        if source not in self.failures:
            self.failures[source] = (step, message)
            self.say(step, f"FAILED for {source}: {message}")

    def _guard(self, source, step, fn, *args, record=None):
        """Run one step for one target; a failure is noted and the night goes on."""
        try:
            return fn(*args)
        except Exception as e:                  # noqa: BLE001 - reported, not hidden
            msg = f"{type(e).__name__}: {e}"
            if record is not None:
                record[source] = msg
            else:
                self._fail(source, step, msg)

    def run(self, sources=None, steps=("flat", "wavecal", "extract", "combine", "telluric", "merge")):
        """
        The steps for `sources` and their standards. Wavelength solutions and
        extractions are independent of one another and run in parallel
        processes; each process runs the same code on the same inputs, so
        the products do not depend on how many there are.
        """
        sci = [s.name for s in self.n.sequences if s.role == "science"]
        sources = sources or sci
        self.failures = {}
        if "flat" in steps:
            self.flat()                      # nothing can go on without it
        needed_stds = {self.n.seq(x).standard for x in sources if self.n.seq(x).standard}
        everyone = list(dict.fromkeys(list(needed_stds) + list(sources)))
        ok = lambda x: x not in self.failures
        nw = self._pool()
        if "wavecal" in steps:
            todo = [w for w in dict.fromkeys(self.n.seq(x).wavecal_from for x in everyone)
                    if self.force or not os.path.exists(self.n.wavecal_path(w))]
            bad = {}
            if nw > 1 and len(todo) > 1:
                self.say("wavecal", f"{len(todo)} solutions on {min(nw, len(todo))} processes")
                for (w,), entry in _parallel(self, "wavecal", [(w,) for w in todo], nw):
                    if isinstance(entry, _Failed):
                        bad[w] = entry.message
                    else:
                        self.log(w, "wavecal", **entry)
            else:
                for w in todo:
                    self._guard(w, "wavecal", self.wavecal, w, record=bad)
            for x in everyone:
                w = self.n.seq(x).wavecal_from
                if w in bad:
                    self._fail(x, "wavecal", f"wavelength solution from {w}: {bad[w]}")
        if "extract" in steps:
            todo = [x for x in everyone if ok(x)]
            jobs = [(x, k) for x in todo for k in range(len(self.units(x)))]
            if nw > 1 and len(jobs) > 1:
                self.say("extract", f"{len(jobs)} extractions on {min(nw, len(jobs))} processes")
                entries = {x: [] for x in todo}
                for (x, k), res in _parallel(self, "extract_unit", jobs, nw):
                    if isinstance(res, _Failed):
                        self._fail(x, "extract", f"unit {k + 1}: {res.message}")
                    elif res[0]:
                        entries[x].append((k, res[0]))
                for x in todo:
                    if ok(x):
                        self._log_extract(x, [e for _, e in sorted(entries[x], key=lambda t: t[0])])
            else:
                for x in todo:
                    self._guard(x, "extract", self.extract, x)
        if "combine" in steps:
            for x in everyone:
                if ok(x):
                    self._guard(x, "combine", self.combine, x)
        for x in sources:
            std = self.n.seq(x).standard
            if ok(x) and std and not ok(std) and ("telluric" in steps or "merge" in steps):
                self._fail(x, "telluric", f"its standard {std} failed")
            if ok(x) and "telluric" in steps:
                self._guard(x, "telluric", self.telluric, x)
            if ok(x) and "merge" in steps:
                self._guard(x, "merge", self.merge, x)
            if ok(x) and (self.o.clean or "clean" in steps):
                self._guard(x, "clean", self.clean, x)
        if self.failures:
            self.say("summary", f"{len(self.failures)} target(s) not finished:")
            for x, (step, msg) in self.failures.items():
                self.say("summary", f"  {x}: {step}: {msg}")
        return {x: self.n.products(x) for x in sources}


# ------------------------------------------------------------- parallel workers

_WORKER = {}


def _worker_reducer(night, options, force):
    r = _WORKER.get("r")
    if r is None or _WORKER.get("key") != (night.out_dir, repr(options), force):
        r = Reducer(night, options, force=force)
        _WORKER.update(r=r, key=(night.out_dir, repr(options), force))
    return r


class _Failed:
    """A worker's exception, carried back to the parent as text."""
    def __init__(self, message):
        self.message = message


def _run_job(args):
    night, options, force, method, job = args
    r = _worker_reducer(night, options, force)
    try:
        if method == "wavecal":
            # The worker returns what the log needs; only the parent writes logs.
            r.log = lambda *a, **k: None
            r.wavecal(job[0])
            return job, r.last_log
        return job, getattr(r, method)(*job)
    except Exception as e:                      # noqa: BLE001 - reported by the parent
        return job, _Failed(f"{type(e).__name__}: {e}")


def _parallel(red, method, jobs, nworkers):
    """Run red.<method>(*job) for every job in a pool of processes, yielding (job, result)."""
    import concurrent.futures as cf
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    args = [(red.n, red.o, red.force, method, j) for j in jobs]
    with cf.ProcessPoolExecutor(max_workers=min(nworkers, len(jobs)), mp_context=ctx) as ex:
        futs = {ex.submit(_run_job, a): a[4] for a in args}
        done = 0
        for f in cf.as_completed(futs):
            done += 1
            red.say(method.split("_")[0], f"{done} of {len(jobs)} done")
            yield f.result()
