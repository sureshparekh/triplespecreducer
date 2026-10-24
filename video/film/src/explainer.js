// Film 1: the science of the reduction, and how the program does each step.
// One builder per scene of video/explainer_script.py; update(lt, D) draws the
// scene at lt seconds into its D seconds.
import { W, H, C, ORDER, clamp, lerp, ease, easeOut, spring, win, vis, el, place, show, plot, svgEl, shotFrame } from "./lib.js";

const SHOT = n => `/gui/light_${n}.png`;
const IMG = n => `/explainer/${n}`;

function card(root, x, y, w, h, cls = "glass") { return place(el("div", cls, root), x, y, w, h); }
function text(root, x, y, w, html, cls = "body") { return place(el("div", cls, root, html), x, y, w); }

export function scenes(D) {
  const S = {};

  S.open = {
    build(r) {
      // the card first, the logo over it: glass blurs whatever is behind it
      this.ring = card(r, 742, 132, 436, 436); this.ring.style.borderRadius = "40px";
      this.logo = place(el("img", "", r), 760, 150, 400, 400); this.logo.src = "/logo/logo_dark.png"; this.logo.style.borderRadius = "30px";
      this.t = text(r, 360, 610, 1200, `From raw frames<br>to a <em>calibrated spectrum</em>`, "title"); this.t.style.textAlign = "center";
      this.s = text(r, 360, 820, 1200, "SOAR TripleSpec 4.1 · 0.94 – 2.47 µm · R ≈ 3500", "label"); this.s.style.textAlign = "center"; this.s.style.fontSize = "24px";
    },
    update(t) {
      show(this.ring, easeOut(win(t, .1, 1.2)), 0, lerp(.9, 1, spring(win(t, .1, 1.4))));
      show(this.logo, easeOut(win(t, .3, 1.2)), 0, lerp(.92, 1, spring(win(t, .3, 1.4))));
      show(this.t, easeOut(win(t, 1.4, 1)), 24); show(this.s, easeOut(win(t, 2.2, 1)), 16);
    },
  };

  S.method = {
    build(r) {
      const cards = [["Spextool", "Cushing, Vacca & Rayner 2004, PASP 116, 362", C.accent],
                     ["Telluric correction", "Vacca, Cushing & Rayner 2003, PASP 115, 389", C.sky],
                     ["TS4_SOAR module", "Spextool for TripleSpec 4.1 at SOAR", C.violet]];
      this.cards = cards.map(([a, b, c], i) => {
        const g = card(r, 140, 190 + i * 210, 760, 170);
        el("div", "eyebrow", g, "the method").style.cssText = `position:absolute;left:36px;top:30px;color:${c}`;
        el("div", "h2", g, a).style.cssText = "position:absolute;left:36px;top:62px";
        el("div", "label", g, b).style.cssText = "position:absolute;left:38px;top:122px;font-size:20px";
        return g;
      });
      this.arrow = svgEl("svg", { width: 260, height: 600 }); Object.assign(this.arrow.style, { position: "absolute", left: "910px", top: "200px" }); r.appendChild(this.arrow);
      this.paths = [0, 1, 2].map(i => svgEl("path", { d: `M10 ${75 + i * 210} C 140 ${75 + i * 210}, 120 300, 250 300`, stroke: "rgba(255,255,255,.45)", "stroke-width": 2.5, fill: "none" }, this.arrow));
      this.tsr = card(r, 1180, 360, 620, 270, "glass");
      const lg = el("img", "", this.tsr); lg.src = "/logo/app_icon_round_256.png"; lg.style.cssText = "position:absolute;left:36px;top:46px;width:120px;border-radius:26px";
      el("div", "h2", this.tsr, "TripleSpec<br>Reducer").style.cssText = "position:absolute;left:184px;top:52px";
      el("div", "label", this.tsr, "the same steps and science").style.cssText = "position:absolute;left:186px;top:176px;font-size:20px";
      this.badge = place(el("div", "glass pill chip", r, `<span class="dot" style="background:${C.green};box-shadow:0 0 14px ${C.green}"></span>results identical to the IDL tool's, bit for bit`), 1180, 680);
    },
    update(t, D) {
      this.cards.forEach((g, i) => show(g, easeOut(win(t, .2 + i * .9, .8)), 30));
      this.paths.forEach((p, i) => { const L = 360; p.style.strokeDasharray = L; p.style.strokeDashoffset = L * (1 - ease(win(t, 3.2 + i * .25, 1))); });
      show(this.tsr, easeOut(win(t, 4.2, .9)), 20, lerp(.96, 1, spring(win(t, 4.2, 1.2))));
      show(this.badge, easeOut(win(t, D * .62, .8)), 16);
    },
  };

  S.instrument = {
    build(r) {
      const s = this.s = svgEl("svg", { width: W, height: 700 }); Object.assign(s.style, { position: "absolute", left: 0, top: 0 }); r.appendChild(s);
      this.slit = svgEl("rect", { x: 170, y: 170, width: 14, height: 240, rx: 5, fill: C.text }, s);
      this.beam = [svgEl("path", { d: "M184 290 L 520 290", stroke: C.accent, "stroke-width": 4, fill: "none" }, s),
                   svgEl("path", { d: "M560 290 L 860 290", stroke: C.accent, "stroke-width": 4, fill: "none" }, s)];
      this.grat = svgEl("g", {}, s);
      for (let i = 0; i < 14; i++) svgEl("line", { x1: 520 + i * 3, y1: 190, x2: 528 + i * 3, y2: 390, stroke: C.sky, "stroke-width": 1.4, opacity: .8 }, this.grat);
      this.prism = svgEl("path", { d: "M860 400 L 960 400 L 910 190 Z", fill: "rgba(200,162,255,.16)", stroke: C.violet, "stroke-width": 2.5 }, s);
      this.fan = [3, 4, 5, 6, 7].map((o, i) => svgEl("path", { d: `M935 300 L 1240 ${525 - i * 95}`, stroke: ORDER[o], "stroke-width": 3, fill: "none", opacity: .9 }, s));
      this.det = svgEl("rect", { x: 1240, y: 80, width: 520, height: 500, rx: 14, fill: "rgba(255,255,255,.03)", stroke: "rgba(255,255,255,.3)", "stroke-width": 2 }, s);
      this.strips = [3, 4, 5, 6, 7].map((o, i) => svgEl("path", { d: `M1265 ${545 - i * 95} Q 1500 ${510 - i * 95} 1735 ${492 - i * 95}`, stroke: ORDER[o], "stroke-width": 26, "stroke-linecap": "round", fill: "none", opacity: .75 }, s));
      this.lab = ["slit", "grating", "prism", "detector"].map((n, i) => { const e = svgEl("text", { x: [177, 540, 910, 1500][i], y: [450, 450, 450, 620][i], fill: C.muted, "font-size": 24, "text-anchor": "middle" }, s); e.textContent = n; return e; });
      this.olab = [3, 4, 5, 6, 7].map((o, i) => { const e = svgEl("text", { x: 1780, y: 500 - i * 95, fill: ORDER[o], "font-size": 24 }, s); e.textContent = `order ${o}`; return e; });
      // wavelength coverage of the orders
      const P = this.cov = plot(r, 260, 636, 1400, 124, [0.9, 2.5], [0, 1], { xticks: [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4], xlabel: "wavelength (µm)", pad: [40, 18, 10, 46] });
      this.bars = Object.entries(D.orders).map(([o, [a, b]], i) => P.band(a, b, ORDER[o], .85));
      this.bars.forEach((b, i) => { b.setAttribute("y", 20 + (i % 2) * 40); b.setAttribute("height", 26); b.setAttribute("rx", 8); });
    },
    update(t, Dur) {
      const a = k => easeOut(win(t, k, .7));
      this.slit.style.opacity = a(.3); this.lab[0].style.opacity = a(.3);
      this.beam[0].style.strokeDasharray = 340; this.beam[0].style.strokeDashoffset = 340 * (1 - ease(win(t, 1.0, .9)));
      this.grat.style.opacity = a(1.2); this.lab[1].style.opacity = a(1.2);
      this.beam[1].style.strokeDasharray = 300; this.beam[1].style.strokeDashoffset = 300 * (1 - ease(win(t, 2.0, .9)));
      this.prism.style.opacity = a(2.4); this.lab[2].style.opacity = a(2.4);
      this.fan.forEach((f, i) => { f.style.strokeDasharray = 330; f.style.strokeDashoffset = 330 * (1 - ease(win(t, 3.2 + i * .15, .8))); });
      this.det.style.opacity = a(3.4); this.lab[3].style.opacity = a(3.4);
      this.strips.forEach((s_, i) => { s_.style.strokeDasharray = 480; s_.style.strokeDashoffset = 480 * (1 - ease(win(t, 4.0 + i * .2, .9))); });
      this.olab.forEach((s_, i) => s_.style.opacity = a(4.6 + i * .2));
      this.cov.svg.style.opacity = a(Dur * .55);
      this.bars.forEach((b, i) => b.style.opacity = .85 * a(Dur * .58 + i * .25));
    },
  };

  S.frame = {
    build(r) {
      const R = D.raw, w = 1400, h = Math.round(w * R.size[1] / R.size[0]);
      this.g = card(r, 260 - 14, 120 - 14, w + 28, h + 28, "glass dark");
      this.im = place(el("img", "", r), 260, 120, w, h); this.im.src = IMG("raw.png"); this.im.style.borderRadius = "16px";
      const s = this.s = svgEl("svg", { width: w, height: h }); Object.assign(s.style, { position: "absolute", left: "260px", top: "120px" }); r.appendChild(s);
      const X = x => x / R.nx * w, Y = y => h - y / R.ny * h;
      this.edges = Object.entries(R.edges).map(([o, ed]) => ed.map(row => svgEl("path", { d: row.map((y, i) => `${i ? "L" : "M"}${X(i * 32).toFixed(1)},${Y(y).toFixed(1)}`).join(""), stroke: ORDER[o], "stroke-width": 2.5, fill: "none" }, s)));
      this.labs = Object.entries(R.edges).map(([o, ed]) => { const e = svgEl("text", { x: w - 190, y: Y(ed[1][Math.floor(ed[1].length * .85)]) - 10, fill: ORDER[o], "font-size": 26, "font-weight": 600 }, s); e.textContent = `order ${o}`; return e; });
      this.note = place(el("div", "glass pill chip", r, `<span class="dot"></span>the star, twice: at A and at B`), 760, 960);
    },
    update(t, Dur) {
      show(this.g, easeOut(win(t, 0, .8))); show(this.im, easeOut(win(t, .1, .8)));
      this.edges.forEach((pr, i) => pr.forEach(p => { const L = 1500; p.style.strokeDasharray = L; p.style.strokeDashoffset = L * (1 - ease(win(t, 1 + i * .35, 1.2))); }));
      this.labs.forEach((l, i) => l.style.opacity = easeOut(win(t, 1.8 + i * .35, .6)));
      show(this.note, easeOut(win(t, Dur * .6, .7)), 14);
    },
  };

  S.atmosphere = {
    build(r) {
      const P = this.P = plot(r, 140, 190, 1640, 560, [0.92, 2.5], [-0.02, 1.05], { xticks: [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4], xlabel: "wavelength (µm)" });
      this.gaps = [[1.34, 1.46], [1.80, 1.96]].map(([a, b]) => P.band(a, b, C.red, .14));
      this.curve = P.line(D.atran, { stroke: "rgba(200,210,230,.85)", "stroke-width": 1.3 });
      this.bands = [["J", 1.25], ["H", 1.63], ["K", 2.2]].map(([n, x]) => P.text(x, 0.45, n, { "font-size": 64, "font-family": "Instrument Serif", fill: C.accent, "text-anchor": "middle", "font-style": "italic" }));
      this.t = text(r, 140, 108, 1600, "What the atmosphere lets through", "h2");
      this.gl = text(r, 1360, 790, 420, "opaque between the bands", "label"); this.gl.style.color = C.red;
    },
    update(t, Dur) {
      show(this.t, easeOut(win(t, 0, .8)), 16);
      this.curve.draw(ease(win(t, .5, 3.5)));
      this.bands.forEach((b, i) => b.style.opacity = easeOut(win(t, 3.5 + i * .4, .6)));
      this.gaps.forEach(g => g.style.opacity = .09 * easeOut(win(t, Dur * .45, .8)));
      show(this.gl, easeOut(win(t, Dur * .48, .8)), 10);
    },
  };

  S.sky = {
    build(r) {
      const ys = D.oh.map(p => p[1]), top = Math.max(...ys);
      const P = this.P = plot(r, 140, 200, 1640, 560, [D.oh[0][0], D.oh[D.oh.length - 1][0]], [-top * .03, top * 1.08], { xticks: [1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8], xlabel: "wavelength (µm), order 4" });
      this.curve = P.line(D.oh, { stroke: C.green, "stroke-width": 1.6 });
      // the strongest lines, marked
      const peaks = D.oh.filter((p, i) => i > 2 && i < D.oh.length - 3 && p[1] > top * .35 && p[1] >= D.oh[i - 1][1] && p[1] >= D.oh[i + 1][1]);
      this.marks = peaks.slice(0, 16).map(p => svgEl("circle", { cx: P.X(p[0]), cy: P.Y(p[1]) - 14, r: 6, fill: "none", stroke: C.accent, "stroke-width": 2.5 }, P.svg));
      this.t = text(r, 140, 108, 1600, "The night sky: hydroxyl emission lines", "h2");
      this.k = place(el("div", "glass pill chip", r, `<span class="dot"></span>known to a fraction of an ångström: they calibrate the wavelengths`), 140, 790);
    },
    update(t, Dur) {
      show(this.t, easeOut(win(t, 0, .8)), 16);
      this.curve.draw(ease(win(t, .4, 3.2)));
      this.marks.forEach((m, i) => { m.style.opacity = easeOut(win(t, Dur * .55 + i * .08, .4)); m.setAttribute("r", (6 * spring(win(t, Dur * .55 + i * .08, .6))).toFixed(2)); });
      show(this.k, easeOut(win(t, Dur * .55, .7)), 14);
    },
  };

  S.nodding = {
    build(r) {
      const mk = (x, lab) => {
        const g = card(r, x, 170, 190, 560, "glass dark");
        const s = svgEl("svg", { width: 190, height: 560 }); Object.assign(s.style, { position: "absolute", left: 0, top: 0 }); g.appendChild(s);
        for (let k = 0; k < 9; k++) svgEl("line", { x1: 25, y1: 60 + k * 55, x2: 165, y2: 60 + k * 55, stroke: C.blue, "stroke-width": 3, opacity: .45 }, s);
        const star = svgEl("circle", { cx: 95, cy: 160, r: 17, fill: C.accent }, s);
        el("div", "h2", g, lab).style.cssText = "position:absolute;left:0;right:0;top:-70px;text-align:center";
        return { g, s, star };
      };
      this.A = mk(260, "A"); this.B = mk(620, "B"); this.B.star.setAttribute("cy", 400);
      this.minus = text(r, 487, 400, 100, "−", "title"); this.eq = text(r, 845, 400, 100, "=", "title");
      this.R = mk(980, "A − B"); this.R.s.querySelectorAll("line").forEach(l => l.setAttribute("opacity", 0));
      this.R.star.setAttribute("fill", C.green); this.neg = svgEl("circle", { cx: 95, cy: 400, r: 17, fill: C.red }, this.R.s);
      this.note = card(r, 1260, 280, 520, 330, "glass");
      this.note.innerHTML = `<div class="h3" style="position:absolute;left:36px;top:34px">The sky cancels</div><div class="body" style="position:absolute;left:36px;top:96px;width:450px;font-size:24px">Sky lines, the telescope's warm glow and the dark current are the same in A and B. The star is left twice: <span style="color:${C.green}">positive</span> and <span style="color:${C.red}">negative</span>.</div>`;
      this.real = ["A.png", "B.png", "AB.png"].map((n, i) => { const im = place(el("img", "", r), 150 + i * 560, 330, 520, 248); im.src = IMG(n); im.style.borderRadius = "14px"; im.style.objectFit = "cover"; return im; });
      this.rl = ["A", "B", "A − B: the sky gone"].map((n, i) => text(r, 150 + i * 560, 600, 520, n, "label"));
    },
    update(t, Dur) {
      // the drawings, then the real frames in their place
      const out = 1 - easeOut(win(t, Dur * .42, .6));
      show(this.A.g, easeOut(win(t, .2, .7)) * out, 30);
      // the star moves from A to B along the slit
      const m = ease(win(t, 1.6, 1.2));
      show(this.B.g, easeOut(win(t, 1.4, .7)) * out, 30);
      this.B.star.setAttribute("cy", lerp(160, 400, m).toFixed(1));
      show(this.minus, easeOut(win(t, 3.0, .5)) * out); show(this.eq, easeOut(win(t, 3.6, .5)) * out);
      show(this.R.g, easeOut(win(t, 3.8, .8)) * out, 30);
      show(this.note, easeOut(win(t, 4.6, .8)) * out, 20);
      this.real.forEach((im, i) => show(im, easeOut(win(t, Dur * .48 + i * .4, .7)), 20));
      this.rl.forEach((l, i) => show(l, easeOut(win(t, Dur * .5 + i * .4, .6))));
    },
  };

  S.night = {
    build(r) {
      const items = [["Flats", "lamp on and off", C.accent], ["Arcs", "lamp lines (optional)", C.sky], ["Targets", "nodded A and B", C.green], ["A0 V standards", "near each target", C.red]];
      this.cards = items.map(([a, b, c], i) => { const g = card(r, 110, 140 + i * 175, 560, 145); g.innerHTML = `<div class="h3" style="position:absolute;left:34px;top:30px;color:${c}">${a}</div><div class="label" style="position:absolute;left:36px;top:84px;font-size:21px">${b}</div>`; return g; });
      this.F = shotFrame(r, SHOT("01-night"), 760, 150, 1040, 650);
      this.chip = place(el("div", "glass pill chip", r, `<span class="dot"></span>one folder → the headers → the plan`), 760, 860);
    },
    update(t, Dur) {
      this.cards.forEach((g, i) => show(g, easeOut(win(t, .2 + i * .5, .7)), 24));
      const a = easeOut(win(t, Dur * .45, .8)); this.F.all.forEach(e => show(e, a, 24));
      const z = ease(win(t, Dur * .62, 1.6)); this.F.set(lerp(0, 40, z), lerp(0, 340, z), lerp(1600, 580, z), lerp(1000, 360, z));
      show(this.chip, easeOut(win(t, Dur * .5, .8)), 14);
    },
  };

  S.flat = {
    build(r) {
      this.lamp = place(el("img", "", r), 140, 200, 520, 474); this.lamp.src = IMG("lamp.png");
      this.flatim = place(el("img", "", r), 760, 200, 520, 474); this.flatim.src = IMG("flat.png");
      [this.lamp, this.flatim].forEach(i => { i.style.borderRadius = "16px"; i.style.objectFit = "cover"; });
      this.l1 = text(r, 140, 690, 520, "a lamps-on flat", "label"); this.l2 = text(r, 760, 690, 520, "the normalised flat: 1 ± a few per cent", "label");
      this.eq = card(r, 1340, 260, 470, 300, "glass");
      this.eq.innerHTML = `<div class="eyebrow" style="position:absolute;left:36px;top:34px">the flat</div><div class="eq" style="position:absolute;left:36px;top:96px"><span class="frac"><span>on − off</span><span>smooth model</span></span></div>`;
      this.F = shotFrame(r, SHOT("02-flat"), 560, 600, 760, 380);
    },
    update(t, Dur) {
      show(this.lamp, easeOut(win(t, .2, .8)), 20); show(this.l1, easeOut(win(t, .6, .6)));
      show(this.eq, easeOut(win(t, 2.2, .8)), 20);
      const f = easeOut(win(t, Dur * .4, .9)); show(this.flatim, f, 20); show(this.l2, f);
      const g = vis(t, Dur * .7, Dur + 1, .8, .1);
      this.F.all.forEach(e => show(e, g, 20)); this.F.set(660, 120, 940, 470);
      [this.lamp, this.flatim, this.l1, this.l2].forEach(e => e.style.filter = `blur(${(g * 6).toFixed(1)}px)`);
    },
  };

  S.wave = {
    build(r) {
      const ys = D.oh.map(p => p[1]), top = Math.max(...ys);
      this.P = plot(r, 120, 120, 1100, 420, [D.oh[0][0], D.oh[D.oh.length - 1][0]], [0, top * 1.05], { xticks: [1.5, 1.6, 1.7, 1.8], xlabel: "µm" });
      this.c = this.P.line(D.oh, { stroke: C.green, "stroke-width": 1.4 });
      this.eq = card(r, 1290, 150, 520, 330, "glass");
      this.eq.innerHTML = `<div class="eyebrow" style="position:absolute;left:36px;top:34px">one polynomial, all five orders</div><div class="eq" style="position:absolute;left:36px;top:104px">λ m / 4 = Σ a<span class="dn">ij</span> x<span class="up">i</span> m<span class="up">j</span></div><div class="label" style="position:absolute;left:38px;top:220px;font-size:21px">x: column · m: order · rms ≈ 0.3 Å</div>`;
      this.wimg = place(el("img", "", r), 120, 190, 800, 520); this.wimg.src = IMG("wavecal.png"); this.wimg.style.borderRadius = "16px"; this.wimg.style.objectFit = "cover";
      this.wl = text(r, 120, 730, 800, "a wavelength for every pixel", "label");
      this.F = shotFrame(r, SHOT("02-wavecal"), 1000, 190, 800, 500);
    },
    update(t, Dur) {
      const out = 1 - easeOut(win(t, Dur * .48, .6));
      this.c.draw(ease(win(t, .3, 2.5))); this.P.svg.style.opacity = out;
      show(this.eq, easeOut(win(t, Dur * .2, .8)) * out, 20);
      show(this.wimg, easeOut(win(t, Dur * .52, .8)), 20); show(this.wl, easeOut(win(t, Dur * .55, .6)));
      const a = easeOut(win(t, Dur * .68, .8)); this.F.all.forEach(e => show(e, a, 20));
      const z = ease(win(t, Dur * .75, 1.5)); this.F.set(lerp(0, 40, z), lerp(0, 440, z), lerp(1600, 600, z), lerp(1000, 375, z));
    },
  };

  S.extract = {
    build(r) {
      this.ab = place(el("img", "", r), 120, 130, 760, 362); this.ab.src = IMG("AB.png"); this.ab.style.borderRadius = "16px";
      this.abl = text(r, 120, 500, 760, "A − B, divided by the flat", "label");
      const pr = D.profile, ys = pr.curve.map(p => p[1]);
      const P = this.P = plot(r, 960, 120, 860, 420, [0, 38], [Math.min(...ys) * 1.15, Math.max(...ys) * 1.25], { xticks: [0, 10, 20, 30], xlabel: "arcsec along the slit" });
      this.ap = pr.peaks.map(c => P.band(c - pr.radius, c + pr.radius, C.green, .28));
      this.bg = []; pr.peaks.forEach(c => [-1, 1].forEach(sg => { const a = c + sg * pr.bg[0], b = c + sg * pr.bg[1]; this.bg.push(P.band(Math.min(a, b), Math.max(a, b), C.red, .15)); }));
      this.curve = P.line(pr.curve, { stroke: C.text, "stroke-width": 2.6 });
      this.eq = card(r, 300, 640, 1320, 210, "glass");
      this.eq.innerHTML = `<div class="eyebrow" style="position:absolute;left:40px;top:32px">the extracted flux, at every wavelength</div><div class="eq" style="position:absolute;left:40px;top:84px">f(λ) = Σ<span class="dn">y</span> w(y) [ (A − B)(λ, y) − background(λ, y) ]</div>`;
      this.F = shotFrame(r, SHOT("03-extract-preview-std"), 560, 600, 800, 400);
    },
    update(t, Dur) {
      show(this.ab, easeOut(win(t, .2, .8)), 20); show(this.abl, easeOut(win(t, .5, .6)));
      this.curve.draw(ease(win(t, Dur * .25, 1.8)));
      this.ap.forEach(b => b.style.opacity = .28 * easeOut(win(t, Dur * .4, .7)));
      this.bg.forEach(b => b.style.opacity = .16 * easeOut(win(t, Dur * .5, .7)));
      const e = vis(t, Dur * .6, Dur * .82, .8, .6); show(this.eq, e, 20);
      const g = easeOut(win(t, Dur * .82, .8)); this.F.all.forEach(x => show(x, g, 20)); this.F.set(660, 120, 940, 500);
    },
  };

  S.combine = {
    build(r) {
      const all = D.combine.spectra.flat().concat(D.combine.combined), ys = all.map(p => p[1]);
      const P = this.P = plot(r, 140, 130, 1100, 600, [1.20, 1.27], [Math.min(...ys) * .97, Math.max(...ys) * 1.02], { xticks: [1.20, 1.22, 1.24, 1.26], xlabel: "µm, order 5" });
      const cols = ["#8b95a3", "#b7a27c", "#8fb59a", "#b79ab8"];
      this.ss = D.combine.spectra.map((s, i) => P.line(s, { stroke: cols[i % 4], "stroke-width": 1.4, opacity: .85 }));
      this.cb = P.line(D.combine.combined, { stroke: C.text, "stroke-width": 3 });
      this.eq = card(r, 1300, 160, 500, 300, "glass");
      this.eq.innerHTML = `<div class="eyebrow" style="position:absolute;left:36px;top:34px">robust weighted mean</div><div class="eq" style="position:absolute;left:36px;top:94px"><span class="frac"><span>Σ f<span class="dn">i</span>/σ<span class="dn">i</span><span class="up">2</span></span><span>Σ 1/σ<span class="dn">i</span><span class="up">2</span></span></span></div><div class="label" style="position:absolute;left:38px;top:236px;font-size:21px">outliers beyond 8σ rejected</div>`;
      this.sn = place(el("div", "glass pill chip", r, `<span class="dot"></span>signal to noise ∝ √N`), 1300, 520);
      this.F = shotFrame(r, SHOT("04-combine-sci"), 1240, 640, 560, 350);
    },
    update(t, Dur) {
      this.ss.forEach((s, i) => s.draw(ease(win(t, .3 + i * .5, 1.6))));
      this.cb.draw(ease(win(t, Dur * .45, 2.0)));
      show(this.eq, easeOut(win(t, Dur * .42, .8)), 20); show(this.sn, easeOut(win(t, Dur * .7, .7)), 14);
      const g = easeOut(win(t, Dur * .6, .8)); this.F.all.forEach(x => show(x, g, 20)); this.F.set(660, 120, 940, 520);
    },
  };

  S.telluric = {
    build(r) {
      this.e1 = card(r, 360, 210, 1200, 180, "glass");
      this.e1.innerHTML = `<div class="eq" style="position:absolute;left:0;right:0;top:56px;text-align:center">O<span class="dn">galaxy</span> = F<span class="dn">galaxy</span> × T<span class="dn">atmosphere</span> × R<span class="dn">instrument</span></div>`;
      this.e2 = card(r, 360, 470, 1200, 230, "glass");
      this.e2.innerHTML = `<div class="eq" style="position:absolute;left:0;right:0;top:64px;text-align:center;color:${C.accent}">F<span class="dn">galaxy</span> = O<span class="dn">galaxy</span> × <span class="frac"><span>F<span class="dn">Vega model</span></span><span>O<span class="dn">A0 V star</span></span></span></div>`;
      this.l = place(el("div", "glass pill chip", r, `<span class="dot"></span>the star went through the same atmosphere and instrument, minutes earlier`), 430, 780);
    },
    update(t, Dur) {
      show(this.e1, easeOut(win(t, .3, .9)), 20);
      show(this.e2, easeOut(win(t, Dur * .45, .9)), 20, lerp(.96, 1, spring(win(t, Dur * .45, 1.2))));
      show(this.l, easeOut(win(t, Dur * .35, .8)), 14);
    },
  };

  S.kernel = {
    build(r) {
      const K = D.kernel, ys = K.vega.concat(K.star).map(p => p[1]);
      const P = this.P = plot(r, 140, 120, 1120, 640, [K.star[0][0], K.star[K.star.length - 1][0]], [Math.min(...ys) * 1.1, 0.04], { xticks: [0.998, 1.002, 1.006, 1.010], xlabel: "µm · Paschen δ, order 6" });
      this.star = P.line(K.star, { stroke: C.text, "stroke-width": 5, opacity: .55 });
      this.vega = P.line(K.vega, { stroke: C.red, "stroke-width": 2.4 });
      this.conv = P.line(K.vega, { stroke: C.green, "stroke-width": 2.6, "stroke-dasharray": "none" });
      this.legend = card(r, 1320, 150, 480, 250, "glass");
      this.legend.innerHTML = `<div class="label" style="position:absolute;left:34px;top:36px;font-size:22px;color:#eef0f5">━ the star</div><div class="label" style="position:absolute;left:34px;top:86px;font-size:22px;color:${C.red}">━ Vega's model</div><div class="label" style="position:absolute;left:34px;top:136px;font-size:22px;color:${C.green}">━ Vega ∗ kernel</div><div class="label" style="position:absolute;left:34px;top:186px;font-size:20px">shifted ${K.vshift.toFixed(1)} km/s</div>`;
      this.eq = card(r, 1320, 450, 480, 200, "glass");
      this.eq.innerHTML = `<div class="eyebrow" style="position:absolute;left:34px;top:32px">the kernel</div><div class="eq" style="position:absolute;left:34px;top:82px;font-size:42px">K̂ = Ô<span class="dn">star</span> / V̂<span class="dn">Vega</span></div>`;
      this.dev = place(el("div", "glass pill chip", r, `<span class="dot" style="background:${C.green};box-shadow:0 0 14px ${C.green}"></span>max ${(100 * K.maxdev).toFixed(1)} % · rms ${(100 * K.rmsdev).toFixed(1)} %`), 1320, 700);
    },
    update(t, Dur) {
      this.star.draw(ease(win(t, .2, 1.4))); this.vega.draw(ease(win(t, 1.0, 1.4)));
      show(this.legend, easeOut(win(t, 1.2, .8)), 16);
      show(this.eq, easeOut(win(t, Dur * .45, .8)), 16);
      const u = ease(win(t, Dur * .55, 2.4));
      this.P.morph(this.conv, D.kernel.vega, D.kernel.conv, u); this.conv.style.opacity = easeOut(win(t, Dur * .55, .4));
      show(this.dev, easeOut(win(t, Dur * .8, .7)), 14);
    },
  };

  S.tellspec = {
    build(r) {
      const T = D.tell, mk = (pts, y, h, col, lab) => {
        const ys = pts.map(p => p[1]).sort((a, b) => a - b), lo = ys[Math.floor(ys.length * .01)], hi = ys[Math.floor(ys.length * .99)];
        const P = plot(r, 140, y, 1100, h, [1.45, 1.80], [Math.min(0, lo), hi * 1.15], { xticks: y > 600 ? [1.5, 1.6, 1.7, 1.8] : [], xlabel: y > 600 ? "µm, order 4" : "", pad: [60, 18, 14, y > 600 ? 46 : 10] });
        const l = P.line(pts, { stroke: col, "stroke-width": 1.6 });
        const tx = text(r, 220, y + 6, 900, lab, "label"); tx.style.color = col;
        return { P, l, tx };
      };
      this.a = mk(T.std, 100, 260, C.text, "the A0 V standard, as observed");
      this.b = mk(T.tell, 380, 250, C.accent, "the telluric spectrum = Vega model ÷ observed star");
      this.c = mk(T.before, 650, 330, "rgba(160,170,190,.8)", "NGC 134 as observed (grey), and corrected (blue)");
      this.d = this.c.P.line(T.after, { stroke: C.blue, "stroke-width": 2 });
      this.F = shotFrame(r, SHOT("05-telluric-corrected"), 1290, 300, 520, 330);
    },
    update(t, Dur) {
      [this.a, this.b, this.c].forEach((p, i) => { p.l.draw(ease(win(t, .3 + i * Dur * .18, 1.8))); show(p.tx, easeOut(win(t, .3 + i * Dur * .18, .6))); });
      this.d.draw(ease(win(t, Dur * .62, 2)));
      const g = easeOut(win(t, Dur * .7, .8)); this.F.all.forEach(x => show(x, g, 20)); this.F.set(660, 120, 940, 600);
    },
  };

  S.merge = {
    build(r) {
      const ys = D.merge.merged.map(p => p[1]).sort((a, b) => a - b), hi = ys[Math.floor(ys.length * .985)] * 1.35;
      this.P = plot(r, 140, 120, 1640, 560, [0.92, 2.5], [0, hi], { xticks: [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4], xlabel: "wavelength (µm)" });
      this.ol = Object.entries(D.merge.orders).map(([o, pts]) => ({ o, e: this.P.line(pts, { stroke: ORDER[o], "stroke-width": 1.5 }) }));
      this.m = this.P.line(D.merge.merged, { stroke: C.text, "stroke-width": 1.6 });
      this.eq = card(r, 300, 720, 1320, 130, "glass");
      this.eq.innerHTML = `<div class="eyebrow" style="position:absolute;left:40px;top:22px">in each overlap</div><div class="eq" style="position:absolute;left:40px;top:54px;font-size:40px">f = ( f<span class="dn">1</span>/σ<span class="dn">1</span><span class="up">2</span> + f<span class="dn">2</span>/σ<span class="dn">2</span><span class="up">2</span> ) / ( 1/σ<span class="dn">1</span><span class="up">2</span> + 1/σ<span class="dn">2</span><span class="up">2</span> )</div>`;
    },
    update(t, Dur) {
      this.ol.forEach(({ e }, i) => { e.draw(ease(win(t, .2 + i * .4, 1.2))); e.style.transform = `translateY(${(60 * (i - 2) * (1 - ease(win(t, Dur * .45, 1.4)))).toFixed(1)}px)`; e.style.opacity = 1 - .75 * ease(win(t, Dur * .7, .8)); });
      this.m.draw(ease(win(t, Dur * .62, 1.6)));
      show(this.eq, easeOut(win(t, Dur * .3, .8)), 16);
    },
  };

  S.products = {
    build(r) {
      this.tree = card(r, 140, 160, 700, 640, "glass");
      const rows = ["night.json", "cal/flat.fits", "cal/NGC134_wavecal.fits", "proc/NGC134_spec0120.fits  …", "proc/NGC134_comb.fits", "proc/NGC134_tellspec.fits", "proc/NGC134_tc.fits", "proc/NGC134_merged.fits", "proc/NGC134_log.json"];
      this.rows = rows.map((s, i) => { const e = el("div", "mono", this.tree, s); e.style.cssText = `position:absolute;left:40px;top:${44 + i * 62}px;font-size:24px;color:${s.includes("merged") ? C.accent : "#dfe3ec"}`; return e; });
      this.F = shotFrame(r, SHOT("08-products"), 920, 180, 900, 560);
      this.chip = place(el("div", "glass pill chip", r, `<span class="dot"></span>&lt;source&gt;_&lt;step&gt;.fits · a log of every step`), 920, 800);
    },
    update(t, Dur) {
      show(this.tree, easeOut(win(t, .2, .8)), 20);
      this.rows.forEach((e, i) => show(e, easeOut(win(t, .6 + i * .3, .5)), 12));
      const g = easeOut(win(t, Dur * .35, .8)); this.F.all.forEach(x => show(x, g, 20));
      show(this.chip, easeOut(win(t, Dur * .5, .7)), 14);
    },
  };

  S.validation = {
    build(r) {
      this.tab = card(r, 140, 120, 980, 820, "glass");
      const rows = [["flat field", "identical"], ["wavelength solution", "identical"], ["extracted spectra", "identical"], ["combined spectra", "identical"], ["telluric spectrum", "identical *"], ["corrected spectrum", "identical *"], ["merged and cleaned", "identical"], ["every header value", "identical"], ["telluric kernel (IDL's FFT)", "10⁻⁶"]];
      this.rows = rows.map(([a, b], i) => { const e = el("div", "", this.tab, `<span class="body" style="font-size:26px">${a}</span><span class="mono" style="position:absolute;right:40px;color:${b.startsWith("identical") ? C.green : C.accent}">${b.startsWith("identical") ? "✓ " : ""}${b}</span>`); e.style.cssText = `position:absolute;left:40px;right:0;top:${46 + i * 78}px`; return e; });
      this.note = text(r, 180, 880, 900, "* given IDL's kernel", "label");
      this.logo = place(el("img", "", r), 1260, 200, 520, 520); this.logo.src = "/logo/logo_dark.png"; this.logo.style.borderRadius = "36px";
      this.doc = place(el("div", "glass pill chip", r, `<span class="dot"></span>the documentation: tutorial · science · sources`), 1250, 790);
    },
    update(t, Dur) {
      show(this.tab, easeOut(win(t, .2, .8)), 20);
      this.rows.forEach((e, i) => show(e, easeOut(win(t, .8 + i * .55, .5)), 10));
      show(this.note, easeOut(win(t, 6, .6)));
      show(this.logo, easeOut(win(t, Dur * .7, 1)), 20, lerp(.95, 1, spring(win(t, Dur * .7, 1.3))));
      show(this.doc, easeOut(win(t, Dur * .78, .8)), 14);
    },
  };
  return S;
}
