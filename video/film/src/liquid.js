const NS = "http://www.w3.org/2000/svg";
const SEL = ".glass";
const ua = navigator.userAgent;
const SUPPORTED = /Chrome\/\d+/.test(ua) && !/CriOS|FxiOS|EdgiOS/.test(ua) && !matchMedia("(prefers-reduced-transparency: reduce)").matches;

let defs = null, n = 0;
const filters = new Map();
const owner = new WeakMap();

function radiusOf(el, w, h) {
  const r = parseFloat(getComputedStyle(el).borderTopLeftRadius) || 0;
  return Math.min(r, w / 2, h / 2);
}

function displacementMap(w, h, r, rim) {
  const k = Math.max(w, h) > 600 ? 0.25 : Math.max(w, h) > 200 ? 0.5 : 1;
  const W = Math.max(2, Math.round(w * k)), H = Math.max(2, Math.round(h * k));
  const c = document.createElement("canvas"); c.width = W; c.height = H;
  const g = c.getContext("2d"), img = g.createImageData(W, H), d = img.data;
  const hx = w / 2, hy = h / 2;
  for (let j = 0; j < H; j++) {
    for (let i = 0; i < W; i++) {
      const px = (i + 0.5) / k - hx, py = (j + 0.5) / k - hy;
      const qx = Math.abs(px) - (hx - r), qy = Math.abs(py) - (hy - r);
      let nx, ny, sd;
      if (qx > 0 && qy > 0) {
        const L = Math.hypot(qx, qy) || 1;
        nx = qx / L; ny = qy / L; sd = L - r;
      } else if (qx > qy) { nx = 1; ny = 0; sd = qx - r; }
      else { nx = 0; ny = 1; sd = qy - r; }
      nx *= Math.sign(px) || 1; ny *= Math.sign(py) || 1;
      const t = Math.min(1, Math.max(0, -sd / rim));
      const m = Math.pow(1 - t, 2.4);
      const o = (j * W + i) * 4;
      d[o] = 128 - nx * m * 127;
      d[o + 1] = 128 - ny * m * 127;
      d[o + 2] = 128; d[o + 3] = 255;
    }
  }
  g.putImageData(img, 0, 0);
  return c.toDataURL();
}

function el(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  parent?.appendChild(e);
  return e;
}

function filterFor(w, h, r, blur) {
  const key = `${w}x${h}x${r}x${blur}`;
  let f = filters.get(key);
  if (f) { f.uses++; return f; }
  const rim = Math.max(10, Math.min(Math.min(w, h) * 0.5, 56));
  const scale = Math.max(10, Math.min(rim * 1.5, 84));
  const id = `lg${++n}`;
  const node = el("filter", { id, x: "0", y: "0", width: "100%", height: "100%", "color-interpolation-filters": "sRGB" }, defs);
  el("feGaussianBlur", { in: "SourceGraphic", stdDeviation: blur, result: "b" }, node);
  el("feImage", { href: displacementMap(w, h, r, rim), x: 0, y: 0, width: w, height: h, preserveAspectRatio: "none", result: "m" }, node);
  if (w * h > 250000) {
    el("feDisplacementMap", { in: "b", in2: "m", scale: scale.toFixed(1), xChannelSelector: "R", yChannelSelector: "G" }, node);
  } else {
    const ch = [["r", scale, "1 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 1 0"], ["g", scale * 0.93, "0 0 0 0 0  0 1 0 0 0  0 0 0 0 0  0 0 0 1 0"], ["bl", scale * 0.86, "0 0 0 0 0  0 0 0 0 0  0 0 1 0 0  0 0 0 1 0"]];
    for (const [nm, s, mat] of ch) {
      el("feDisplacementMap", { in: "b", in2: "m", scale: s.toFixed(1), xChannelSelector: "R", yChannelSelector: "G", result: `d${nm}` }, node);
      el("feColorMatrix", { in: `d${nm}`, type: "matrix", values: mat, result: nm }, node);
    }
    el("feBlend", { in: "r", in2: "g", mode: "screen", result: "rg" }, node);
    el("feBlend", { in: "rg", in2: "bl", mode: "screen" }, node);
  }
  f = { id, key, node, uses: 1 };
  filters.set(key, f);
  return f;
}

function release(f) {
  if (!f || --f.uses > 0) return;
  f.node.remove(); filters.delete(f.key);
}

function apply(target) {
  const w = Math.round(target.offsetWidth / 2) * 2, h = Math.round(target.offsetHeight / 2) * 2;
  const prev = owner.get(target);
  if (w < 4 || h < 4) { if (prev) { release(prev); owner.delete(target); target.style.removeProperty("--lgf"); } return; }
  const r = Math.round(radiusOf(target, w, h));
  const blur = target.classList.contains("dark") ? 9 : Math.min(w, h) < 80 ? 4 : 7;
  const key = `${w}x${h}x${r}x${blur}`;
  if (prev?.key === key) return;
  const f = filterFor(w, h, r, blur);
  release(prev);
  owner.set(target, f);
  target.style.setProperty("--lgf", `url(#${f.id})`);
}

export function initLiquidGlass() {
  if (!SUPPORTED) return;
  const svg = el("svg", { width: "0", height: "0", "aria-hidden": "true", style: "position:absolute;width:0;height:0;overflow:hidden" });
  defs = el("defs", {}, svg);
  document.body.prepend(svg);
  document.documentElement.classList.add("liquid");
  let queue = new Set(), raf = 0;
  const flush = () => { raf = 0; for (const t of queue) if (t.isConnected) apply(t); queue = new Set(); };
  const schedule = t => { queue.add(t); if (!raf) raf = requestAnimationFrame(flush); };
  const ro = new ResizeObserver(es => es.forEach(e => schedule(e.target)));
  const watch = t => { ro.observe(t); schedule(t); };
  document.querySelectorAll(SEL).forEach(watch);
  new MutationObserver(ms => {
    for (const m of ms) {
      m.addedNodes.forEach(nd => { if (nd.nodeType !== 1) return; if (nd.matches(SEL)) watch(nd); nd.querySelectorAll?.(SEL).forEach(watch); });
      m.removedNodes.forEach(nd => {
        if (nd.nodeType !== 1) return;
        const all = nd.matches(SEL) ? [nd, ...nd.querySelectorAll(SEL)] : [...(nd.querySelectorAll?.(SEL) || [])];
        all.forEach(t => { ro.unobserve(t); release(owner.get(t)); owner.delete(t); });
      });
    }
  }).observe(document.body, { childList: true, subtree: true });
}
