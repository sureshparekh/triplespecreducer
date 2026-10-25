// The film runner: ?film=explainer|tutorial. window.FILM.setup() builds the
// scenes once; window.FILM.frame(t) draws the frame at t seconds, so a
// renderer can step through the film exactly, frame by frame.
import { sky, el, place, vis } from "./lib.js";
import { initLiquidGlass } from "./liquid.js";

const film = new URLSearchParams(location.search).get("film") || "explainer";
const root = document.getElementById("world"), hud = document.getElementById("hud");
let TL, frameFns = [];

// each scene, centred and enlarged into the frame below the brand chip:
// its content's box is measured once, as built
function fit(layer) {
  let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
  for (const c of layer.children) {
    const r = c.getBoundingClientRect();
    if (r.width < 2 || r.height < 2) continue;
    x0 = Math.min(x0, r.left); y0 = Math.min(y0, r.top); x1 = Math.max(x1, r.right); y1 = Math.max(y1, r.bottom);
  }
  if (x1 < x0) return;
  const k = Math.min(1.15, 1820 / (x1 - x0), 930 / (y1 - y0));
  const tx = 960 - k * (x0 + x1) / 2, ty = 578 - k * (y0 + y1) / 2;
  Object.assign(layer.style, { transformOrigin: "0 0", transform: `translate(${tx.toFixed(1)}px, ${ty.toFixed(1)}px) scale(${k.toFixed(4)})` });
}

async function setup() {
  TL = await (await fetch(`/${film}/timeline.json`)).json();
  await document.fonts.ready;
  await Promise.all(["400 40px Geist", "600 40px Geist", "italic 40px 'Instrument Serif'", "400 40px 'Instrument Serif'", "500 20px 'Geist Mono'"].map(f => document.fonts.load(f)));
  const drift = sky(document.getElementById("sky"));
  frameFns.push(drift);
  const brand = place(el("div", "glass pill chip brand", hud,
    `<img src="/logo/app_icon_round_64.png" style="width:30px;height:30px;border-radius:8px">TripleSpec Reducer <span class="muted" style="font-size:19px">· ${TL.title}</span>`), 40, 34);
  brand.style.padding = "9px 20px 9px 10px";
  if (film === "explainer") {
    const D = await (await fetch("/explainer/data.json")).json();
    const { scenes } = await import("./explainer.js");
    const S = scenes(D);
    for (const sc of TL.scenes) {
      const s = S[sc.id];
      if (!s) continue;
      const layer = el("div", "", root); Object.assign(layer.style, { position: "absolute", inset: 0 });
      s.build(layer);
      fit(layer);
      frameFns.push(t => {
        const a = vis(t, sc.t0 - .3, sc.t1 + .55, .7, .7);
        layer.style.opacity = a.toFixed(3);
        layer.style.visibility = a > .002 ? "visible" : "hidden";
        if (a > .002) s.update(t - sc.t0, sc.t1 - sc.t0);
      });
    }
  } else {
    const { tutorial } = await import("./tutorial.js");
    frameFns.push(await tutorial(root, hud, TL));
  }
  initLiquidGlass();
  await new Promise(r => setTimeout(r, 400));
  return { duration: TL.duration };
}

async function frame(t) {
  for (const f of frameFns) await f(t);
  await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));
}

window.FILM = { setup, frame };
