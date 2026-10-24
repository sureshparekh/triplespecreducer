// Shared pieces of both films: elements, easing, time windows, the ground's
// glows, the caption pill and SVG plots that draw themselves.

export const W = 1920, H = 1080;
export const C = { bg: "#06070b", text: "#eef0f5", muted: "#a3aaba", dim: "#6b7385", accent: "#f4a259",
  orange: "#f39b4a", red: "#ff6a3d", blue: "#5b8cff", sky: "#8fb6ff", violet: "#c8a2ff", green: "#7fd0a8" };
export const ORDER = { 3: "#f4a259", 4: "#7fd0a8", 5: "#5b8cff", 6: "#c8a2ff", 7: "#ff6a3d" };

export const clamp = (x, a = 0, b = 1) => Math.max(a, Math.min(b, x));
export const lerp = (a, b, u) => a + (b - a) * u;
export const ease = u => { u = clamp(u); return u < .5 ? 4 * u * u * u : 1 - Math.pow(-2 * u + 2, 3) / 2; };
export const easeOut = u => 1 - Math.pow(1 - clamp(u), 3);
export const spring = u => { u = clamp(u); return 1 - Math.pow(2, -10 * u) * Math.cos(u * 9.42); };
// progress of t through [a, a+d]
export const win = (t, a, d) => clamp((t - a) / d);
// fade in at a over fi, out at b over fo
export const vis = (t, a, b, fi = .6, fo = .6) => Math.min(easeOut(win(t, a, fi)), 1 - ease(win(t, b - fo, fo)));

export function el(tag, cls, parent, html) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html != null) e.innerHTML = html;
  parent?.appendChild(e);
  return e;
}
export function place(e, x, y, w, h) {
  Object.assign(e.style, { position: "absolute", left: `${x}px`, top: `${y}px` });
  if (w != null) e.style.width = `${w}px`;
  if (h != null) e.style.height = `${h}px`;
  return e;
}
export function show(e, a, dy = 0, s = 1) {
  e.style.opacity = a.toFixed(3);
  e.style.visibility = a > 0.002 ? "visible" : "hidden";
  e.style.transform = `translate3d(0, ${(dy * (1 - a)).toFixed(2)}px, 0) scale(${s.toFixed(4)})`;
}

// the ground: three soft glows in the site's orange and blue, drifting slowly
export function sky(root) {
  const g = [["#f39b4a", 900, .16], ["#5b8cff", 980, .15], ["#c8a2ff", 700, .08], ["#ff6a3d", 600, .07]]
    .map(([c, s, o]) => { const d = el("div", "glow", root); Object.assign(d.style, { width: `${s}px`, height: `${s}px`, background: c, opacity: o }); return d; });
  return t => {
    const p = [[-.08, .78, .11], [.86, .05, .07], [.55, .55, .05], [.15, .1, .09]];
    g.forEach((d, i) => {
      const [x, y, sp] = p[i];
      const dx = Math.sin(t * sp + i) * 120, dy = Math.cos(t * sp * .8 + i * 2) * 90;
      d.style.transform = `translate3d(${x * W + dx - 450}px, ${y * H + dy - 450}px, 0)`;
    });
  };
}

// the caption: the line being spoken, its words lit as they are said
export function captions(root, lines) {
  const box = el("div", "glass dark cap pill", root);
  let cur = -1, spans = [];
  return t => {
    let k = -1;
    for (let i = 0; i < lines.length; i++) if (t >= lines[i].t0 - .15 && t <= lines[i].t1 + .5) k = i;
    if (k !== cur) {
      cur = k; box.innerHTML = "";
      spans = k < 0 ? [] : lines[k].shown.split(/\s+/).map((w, j) => {
        const s = el("span", "w", box); s.textContent = (j ? " " : "") + w; return s; });
    }
    if (k < 0) { box.style.opacity = 0; return; }
    const L = lines[k], a = vis(t, L.t0 - .15, L.t1 + .5, .25, .3);
    box.style.opacity = a.toFixed(3);
    // light the shown words in proportion to the spoken ones
    const spoken = L.words.filter(w => w[1] <= t).length / Math.max(1, L.words.length);
    const lit = Math.round(spoken * spans.length);
    spans.forEach((s, j) => { s.style.color = j < lit ? "#ffffff" : "rgba(233,236,243,.55)"; });
  };
}

// an SVG plot: axes in data units, paths that draw on or morph
const NS = "http://www.w3.org/2000/svg";
export function svgEl(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  parent?.appendChild(e);
  return e;
}
export function plot(parent, x, y, w, h, xlim, ylim, opt = {}) {
  const svg = svgEl("svg", { width: w, height: h, viewBox: `0 0 ${w} ${h}` });
  Object.assign(svg.style, { position: "absolute", left: `${x}px`, top: `${y}px`, overflow: "visible" });
  parent.appendChild(svg);
  const pad = opt.pad ?? [60, 18, 20, 46];      // left, right, top, bottom
  const X = v => pad[0] + (v - xlim[0]) / (xlim[1] - xlim[0]) * (w - pad[0] - pad[1]);
  const Y = v => h - pad[3] - (v - ylim[0]) / (ylim[1] - ylim[0]) * (h - pad[2] - pad[3]);
  const axes = svgEl("g", { opacity: 1 }, svg);
  svgEl("line", { x1: pad[0], y1: h - pad[3], x2: w - pad[1], y2: h - pad[3], stroke: "rgba(255,255,255,.22)", "stroke-width": 1.2 }, axes);
  if (opt.xticks) for (const v of opt.xticks) {
    svgEl("line", { x1: X(v), y1: h - pad[3], x2: X(v), y2: h - pad[3] + 7, stroke: "rgba(255,255,255,.3)" }, axes);
    const tx = svgEl("text", { x: X(v), y: h - pad[3] + 30, fill: C.muted, "font-size": 18, "text-anchor": "middle" }, axes); tx.textContent = v;
  }
  if (opt.xlabel) { const tx = svgEl("text", { x: (pad[0] + w - pad[1]) / 2, y: h + 8, fill: C.dim, "font-size": 18, "text-anchor": "middle" }, axes); tx.textContent = opt.xlabel; }
  const P = {
    svg, X, Y, axes,
    line(pts, attrs) {
      // inside the axes only: points beyond the x range dropped, values held to the y range
      pts = pts.filter(p => p[0] >= xlim[0] && p[0] <= xlim[1]).map(p => [p[0], clamp(p[1], ylim[0], ylim[1])]);
      // a jump in x (a masked band) breaks the line
      const dx = pts.slice(1).map((p, i) => p[0] - pts[i][0]).sort((u, v) => u - v)[pts.length >> 1] || 0;
      const d = pts.map((p, i) => `${i && p[0] - pts[i - 1][0] < 25 * dx ? "L" : "M"}${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join("");
      const e = svgEl("path", Object.assign({ d, fill: "none", stroke: C.text, "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round" }, attrs), svg);
      e._len = null;
      e.draw = u => { if (e._len == null) e._len = e.getTotalLength() || 1; e.style.strokeDasharray = e._len; e.style.strokeDashoffset = (e._len * (1 - clamp(u))).toFixed(1); };
      e.pts = pts;
      return e;
    },
    morph(e, a, b, u) {      // a path between two series on the same x
      const n = Math.min(a.length, b.length);
      const d = []; for (let i = 0; i < n; i++) d.push(`${i ? "L" : "M"}${X(a[i][0]).toFixed(1)},${Y(lerp(a[i][1], b[i][1], u)).toFixed(1)}`);
      e.setAttribute("d", d.join(""));
    },
    band(x0, x1, color, op = .2) {
      return svgEl("rect", { x: X(x0), y: pad[2], width: Math.max(0, X(x1) - X(x0)), height: h - pad[2] - pad[3], fill: color, opacity: op }, svg);
    },
    text(xv, yv, s, attrs) { const e = svgEl("text", Object.assign({ x: X(xv), y: Y(yv), fill: C.text, "font-size": 20 }, attrs), svg); e.textContent = s; return e; },
  };
  return P;
}

// a screenshot in a glass frame, with a camera that can zoom onto a region
export function shotFrame(parent, src, x, y, w, h, imgW = 1600, imgH = 1000) {
  const g = place(el("div", "glass dark", parent), x - 12, y - 12, w + 24, h + 24);
  g.style.borderRadius = "28px";
  const s = place(el("div", "shot", parent), x, y, w, h);
  const im = el("img", "", s); im.src = src; im.width = imgW; im.height = imgH;
  const base = Math.min(w / imgW, h / imgH);
  const set = (rx, ry, rw, rh) => {      // show the image region (rx, ry, rw, rh) in the frame
    const k = Math.min(w / rw, h / rh);
    const ox = w / 2 - (rx + rw / 2) * k, oy = h / 2 - (ry + rh / 2) * k;
    im.style.transform = `translate3d(${ox.toFixed(2)}px, ${oy.toFixed(2)}px, 0) scale(${k.toFixed(5)})`;
  };
  set(0, 0, imgW, imgH);
  return { g, s, im, base, set, all: [g, s] };
}
