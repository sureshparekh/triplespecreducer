// Film 2: the screen recording of the interface, composited: a glass window
// on the ground, a camera that zooms to what the narration is about, the
// cursor and its clicks, the step being worked on, and a badge when a long
// computation is shown faster than it ran.
import { W, H, C, clamp, lerp, ease, easeOut, spring, win, vis, el, place, show } from "./lib.js";

const STEPS = ["Night", "Calibrate", "Extract", "Combine", "Telluric", "Merge", "Clean & save"];

export async function tutorial(root, hud, TL) {
  const [SW, SH] = TL.size;                       // the recorded window, 1600 x 1000
  const k0 = 0.72, ww = SW * k0, wh = SH * k0;    // the window on screen, clear of the brand and the captions
  const wx = (W - ww) / 2, wy = 136;

  // the title card, before the program opens
  const title = place(el("div", "glass", root), 360, 250, 1200, 520);
  title.innerHTML = `<div class="eyebrow" style="position:absolute;left:70px;top:70px">hands-on tutorial</div>
    <div class="title" style="position:absolute;left:66px;top:120px">Reducing NGC 134<br><em>from scratch</em></div>
    <div class="body muted" style="position:absolute;left:70px;top:350px;width:1060px">SOAR TripleSpec 4.1 · the night of 2025 November 2 · every step in the interface, from the raw frames to the final spectrum</div>`;
  // the splash, as the program shows it while it loads
  const splash = place(el("img", "", root), (W - 460) / 2, (H - 514) / 2 - 20, 460, 514);
  splash.src = "/tutorial/rec/splash.png";
  splash.style.cssText += ";border-radius:18px;box-shadow:0 30px 80px rgba(0,0,0,.6)";

  // the window: a glass frame, a title bar, and the screen inside
  const frame = place(el("div", "glass dark", root), wx - 14, wy - 50, ww + 28, wh + 64);
  frame.style.borderRadius = "24px";
  frame.innerHTML = `<div style="position:absolute;left:24px;top:17px;display:flex;gap:9px">${["#ff5f57", "#febc2e", "#28c840"].map(c => `<span style="width:13px;height:13px;border-radius:50%;background:${c};opacity:.85"></span>`).join("")}</div>
    <div class="label" style="position:absolute;left:0;right:0;top:14px;text-align:center;font-size:18px">TripleSpec Reducer — SOAR TripleSpec 4.1</div>`;
  const screen = place(el("div", "", root), wx, wy, ww, wh);
  Object.assign(screen.style, { overflow: "hidden", borderRadius: "10px", background: "#fff" });
  const cam = el("div", "", screen);
  Object.assign(cam.style, { position: "absolute", left: 0, top: 0, width: `${SW}px`, height: `${SH}px`, transformOrigin: "0 0" });
  const img = el("img", "", cam);
  Object.assign(img.style, { position: "absolute", left: 0, top: 0, width: `${SW}px`, height: `${SH}px` });
  // the cursor, a ripple per click
  const cur = el("div", "", cam);
  cur.innerHTML = `<svg width="34" height="40" viewBox="0 0 17 20"><path d="M1 1 L1 15.5 L4.6 12.2 L7.1 18.2 L9.6 17.1 L7.1 11.3 L12 11.3 Z" fill="#111" stroke="#fff" stroke-width="1.3" stroke-linejoin="round"/></svg>`;
  Object.assign(cur.style, { position: "absolute", left: 0, top: 0, filter: "drop-shadow(0 2px 3px rgba(0,0,0,.35))" });
  const ripples = [0, 1, 2].map(() => { const r = el("div", "", cam); Object.assign(r.style, { position: "absolute", width: "60px", height: "60px", marginLeft: "-30px", marginTop: "-30px", borderRadius: "50%", border: `3px solid ${C.accent}`, opacity: 0 }); return r; });

  // the step chip and the speed badge
  const step = place(el("div", "glass pill chip", hud, ""), 0, 34);
  const badge = place(el("div", "glass pill chip", hud, ""), 0, 34);

  // the camera: the zoom keyframes, eased from one to the next
  const Z = TL.zooms.slice().sort((a, b) => a.t - b.t);
  const fit = z => {
    if (!z || z.full) return [0, 0, SW, SH];
    let [x, y, w, h] = z.rect;
    // some air around the region, the window's shape kept, never closer than 2.3x
    const cx = x + w / 2, cy = y + h / 2;
    w = Math.max(w * 1.45, SW / 2.3); h = Math.max(h * 1.45, SH / 2.3);
    if (w / h > SW / SH) h = w * SH / SW; else w = h * SW / SH;
    x = clamp(cx - w / 2, 0, SW - w); y = clamp(cy - h / 2, 0, SH - h);
    return [x, y, w, h];
  };
  const camAt = t => {
    let i = -1;
    for (let k = 0; k < Z.length; k++) if (Z[k].t <= t) i = k;
    const now = fit(i >= 0 ? Z[i] : null), prev = fit(i > 0 ? Z[i - 1] : null);
    const u = i >= 0 ? ease(win(t, Z[i].t, 1.1)) : 1;
    return now.map((v, j) => lerp(prev[j], v, u));
  };

  let lastSrc = "";
  return async t => {
    const fi = Math.min(TL.frames.length - 1, Math.max(0, Math.floor(t * TL.fps)));
    const name = TL.images[TL.frames[fi]];
    const onTitle = name === "@title", onSplash = name === "@splash";
    // the title card, then the splash, then the window
    const firstWin = TL.frames.findIndex(f => !TL.images[f].startsWith("@")) / TL.fps;
    const splashT = TL.frames.findIndex(f => TL.images[f] === "@splash") / TL.fps;
    show(title, vis(t, 0, splashT, .8, .7), 26);
    show(splash, vis(t, splashT, firstWin + .3, .5, .5), 0, lerp(.92, 1, spring(win(t, splashT, .9))));
    const wa = easeOut(win(t, firstWin - .15, .9));
    [frame, screen].forEach(e => show(e, wa, 30, lerp(.975, 1, wa)));
    if (!onTitle && !onSplash && name !== lastSrc) {
      lastSrc = name;
      img.src = `/tutorial/rec/${name}`;
      try { await img.decode(); } catch (e) { /* the previous frame stays */ }
    }
    // the camera
    const [cx, cy, cw, ch] = camAt(t);
    const k = ww / cw;
    cam.style.transform = `translate3d(${(-cx * k).toFixed(2)}px, ${(-cy * k).toFixed(2)}px, 0) scale(${k.toFixed(5)})`;
    // the cursor, scaled back so it keeps its size on screen
    const [px, py] = TL.cursor[fi];
    cur.style.display = px < 0 ? "none" : "block";
    cur.style.transform = `translate3d(${(px - 2).toFixed(1)}px, ${(py - 2).toFixed(1)}px, 0) scale(${(k0 / k).toFixed(4)})`;
    cur.style.transformOrigin = "2px 2px";
    const recent = TL.clicks.filter(c => c <= t && t - c < .6).slice(-3);
    ripples.forEach((r, j) => {
      const c = recent[j];
      if (c == null) { r.style.opacity = 0; return; }
      const u = (t - c) / .6, ci = Math.min(TL.cursor.length - 1, Math.floor(c * TL.fps)), [qx, qy] = TL.cursor[ci];
      r.style.left = `${qx}px`; r.style.top = `${qy}px`;
      r.style.opacity = (.9 * (1 - u)).toFixed(3);
      r.style.transform = `scale(${(lerp(.3, 1.4, easeOut(u)) * k0 / k).toFixed(3)})`;
    });
    // the step being worked on
    const sc = TL.scenes.find(s => t >= s.t0 && t < s.t1) || TL.scenes[TL.scenes.length - 1];
    const label = sc.step ? `<span class="dot"></span>Step ${sc.step} · ${STEPS[sc.step - 1]}` : "";
    if (step._l !== label) { step._l = label; step.innerHTML = label; }
    step.style.left = `${W - 40 - step.offsetWidth}px`;
    show(step, label ? wa : 0, -8);
    // a computation shown faster than it ran
    const b = TL.badges.find(x => t >= x.t0 - .2 && t <= x.t1 + .2);
    if (b) { const s = `<span class="dot" style="background:${C.sky};box-shadow:0 0 14px ${C.sky}"></span>computing · shown ${b.factor}× faster`; if (badge._l !== s) { badge._l = s; badge.innerHTML = s; } }
    badge.style.left = `${(W - badge.offsetWidth) / 2}px`;
    show(badge, b ? vis(t, b.t0 - .2, b.t1 + .2, .3, .3) : 0, -8);
  };
}
