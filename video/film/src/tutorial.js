// Film 2: the screen recording of the interface, composited: a short title
// card, the splash, then the window, large and centred, with the cursor and
// its clicks, and a badge when a long computation is shown faster than it
// ran. The window never moves: no zooming in or out.
import { W, H, C, lerp, easeOut, spring, win, vis, el, place, show } from "./lib.js";

export async function tutorial(root, hud, TL) {
  const [SW, SH] = TL.size;                          // the recorded window, 1600 x 1000
  const k = Math.min((W - 64) / SW, (H - 50) / SH);  // as large as the frame allows
  const ww = SW * k, wh = SH * k, wx = (W - ww) / 2, wy = (H - wh) / 2;

  // the title card, before the program opens
  const title = place(el("div", "glass", root), 460, 330, 1000, 380);
  title.innerHTML = `<div class="eyebrow" style="position:absolute;left:64px;top:64px">hands-on tutorial</div>
    <div class="title" style="position:absolute;left:60px;top:112px">Reducing NGC 134<br><em>from scratch</em></div>`;
  // the splash, as the program shows it while it loads
  const splash = place(el("img", "", root), (W - 460) / 2, (H - 514) / 2, 460, 514);
  splash.src = "/tutorial/rec/splash.png";
  splash.style.cssText += ";border-radius:18px;box-shadow:0 30px 80px rgba(0,0,0,.6)";

  // the window: a thin glass rim, and the screen inside
  const frame = place(el("div", "glass dark", root), wx - 8, wy - 8, ww + 16, wh + 16);
  frame.style.borderRadius = "18px";
  const screen = place(el("div", "", root), wx, wy, ww, wh);
  Object.assign(screen.style, { overflow: "hidden", borderRadius: "11px", background: "#fff" });
  const cam = el("div", "", screen);
  Object.assign(cam.style, { position: "absolute", left: 0, top: 0, width: `${SW}px`, height: `${SH}px`, transformOrigin: "0 0",
                             transform: `scale(${k.toFixed(5)})` });
  const img = el("img", "", cam);
  Object.assign(img.style, { position: "absolute", left: 0, top: 0, width: `${SW}px`, height: `${SH}px` });
  // the cursor, a ripple per click
  const cur = el("div", "", cam);
  cur.innerHTML = `<svg width="26" height="31" viewBox="0 0 17 20"><path d="M1 1 L1 15.5 L4.6 12.2 L7.1 18.2 L9.6 17.1 L7.1 11.3 L12 11.3 Z" fill="#111" stroke="#fff" stroke-width="1.3" stroke-linejoin="round"/></svg>`;
  Object.assign(cur.style, { position: "absolute", left: 0, top: 0, filter: "drop-shadow(0 2px 3px rgba(0,0,0,.35))" });
  const ripples = [0, 1, 2].map(() => { const r = el("div", "", cam); Object.assign(r.style, { position: "absolute", width: "48px", height: "48px", marginLeft: "-24px", marginTop: "-24px", borderRadius: "50%", border: `3px solid ${C.accent}`, opacity: 0 }); return r; });

  // the speed badge, over the bottom of the window
  const badge = place(el("div", "glass pill chip", hud, ""), 0, wy + wh - 76);
  const brand = hud.querySelector(".brand");

  const firstWin = TL.frames.findIndex(f => !TL.images[f].startsWith("@")) / TL.fps;
  const splashT = TL.frames.findIndex(f => TL.images[f] === "@splash") / TL.fps;
  let lastSrc = "";
  return async t => {
    const fi = Math.min(TL.frames.length - 1, Math.max(0, Math.floor(t * TL.fps)));
    const name = TL.images[TL.frames[fi]];
    // the title card, then the splash, then the window
    show(title, vis(t, 0, splashT, .6, .6), 26);
    show(splash, vis(t, splashT, firstWin + .3, .5, .5), 0, lerp(.92, 1, spring(win(t, splashT, .9))));
    const wa = easeOut(win(t, firstWin - .15, .9));
    [frame, screen].forEach(e => show(e, wa, 30, lerp(.975, 1, wa)));
    // the brand chip belongs to the title card; the window shows its own name
    if (brand) brand.style.opacity = (1 - wa).toFixed(3);
    if (!name.startsWith("@") && name !== lastSrc) {
      lastSrc = name;
      img.src = `/tutorial/rec/${name}`;
      try { await img.decode(); } catch (e) { /* the previous frame stays */ }
    }
    const [px, py] = TL.cursor[fi];
    cur.style.display = px < 0 ? "none" : "block";
    cur.style.transform = `translate3d(${(px - 2).toFixed(1)}px, ${(py - 2).toFixed(1)}px, 0)`;
    const recent = TL.clicks.filter(c => c <= t && t - c < .6).slice(-3);
    ripples.forEach((r, j) => {
      const c = recent[j];
      if (c == null) { r.style.opacity = 0; return; }
      const u = (t - c) / .6, ci = Math.min(TL.cursor.length - 1, Math.floor(c * TL.fps)), [qx, qy] = TL.cursor[ci];
      r.style.left = `${qx}px`; r.style.top = `${qy}px`;
      r.style.opacity = (.9 * (1 - u)).toFixed(3);
      r.style.transform = `scale(${lerp(.3, 1.4, easeOut(u)).toFixed(3)})`;
    });
    // a computation shown faster than it ran
    const b = TL.badges.find(x => t >= x.t0 - .2 && t <= x.t1 + .2);
    if (b) { const s = `<span class="dot" style="background:${C.sky};box-shadow:0 0 14px ${C.sky}"></span>computing · shown ${b.factor}× faster`; if (badge._l !== s) { badge._l = s; badge.innerHTML = s; } }
    badge.style.left = `${(W - badge.offsetWidth) / 2}px`;
    show(badge, b ? vis(t, b.t0 - .2, b.t1 + .2, .3, .3) : 0, 8);
  };
}
