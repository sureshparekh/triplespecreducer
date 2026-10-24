"""
Render a film frame by frame, as the reels are rendered: headless Chromium
(Playwright) steps window.FILM.frame(t) and ffmpeg encodes the screenshots
with the mixed audio.

    python video/render.py explainer                 # the whole film, with build/explainer/mix.wav
    python video/render.py tutorial
    python video/render.py explainer preview 3 40.5  # stills at those times, to build/<film>/previews
    python video/render.py explainer parallel 8      # the same film in 8 parts at once, then joined
    python video/render.py explainer part 2 8        # part 2 of 8 again; then: render.py explainer join 8

Every frame depends only on t, so the parts are independent.

Needs the film's page served: video/serve.sh (Vite, port 5241).
"""

import asyncio
import os
import subprocess
import sys
import time
from pathlib import Path

from playwright.async_api import async_playwright

HERE = Path(__file__).parent
GPU = ["--use-angle=vulkan", "--enable-features=Vulkan", "--ignore-gpu-blocklist", "--enable-gpu", "--headless=new"]
CHROME = os.environ.get("CHROME", str(Path.home() / ".cache/ms-playwright/chromium-1208/chrome-linux64/chrome"))
FPS = 30


async def main(film, mode="video", times=(), part=None):
    out = HERE / "build" / film
    async with async_playwright() as p:
        b = await p.chromium.launch(args=GPU, executable_path=CHROME)
        pg = await b.new_page(viewport={"width": 1920, "height": 1080}, device_scale_factor=1)
        pg.on("pageerror", lambda e: print("pageerror", e, flush=True))
        pg.on("response", lambda r: print("missing", r.url, flush=True) if r.status >= 400 else None)
        pg.on("console", lambda m: print("console", m.text, flush=True) if m.type == "error" else None)
        await pg.goto(f"http://localhost:5241/?film={film}", wait_until="networkidle")
        info = await pg.evaluate("window.FILM.setup()")
        # a reload would leave a page that was never set up, and blank frames: stop instead
        moved = []
        pg.on("framenavigated", lambda f: moved.append(f.url) if f == pg.main_frame else None)
        if mode == "preview":
            (out / "previews").mkdir(parents=True, exist_ok=True)
            for t in times:
                await pg.evaluate(f"window.FILM.frame({t})")
                await pg.screenshot(path=str(out / "previews" / f"t{float(t):07.2f}.png"))
                print("preview", t, flush=True)
        else:
            n = int(round(info["duration"] * FPS))
            a, z = 0, n
            enc = ["-c:v", "libx264", "-preset", "slow", "-crf", "18", "-pix_fmt", "yuv420p", "-profile:v", "high", "-r", str(FPS)]
            if part:                     # frames [a, z) only, picture only
                k, m = part
                a, z = n * k // m, n * (k + 1) // m
                (out / "parts").mkdir(parents=True, exist_ok=True)
                dst = out / "parts" / f"p{k:02d}.mp4"
                cmd = ["-i", "-"] + enc + [str(dst)]
            else:
                dst = out / f"TripleSpec_Reducer_{film}.mp4"
                cmd = (["-i", "-", "-i", str(out / "mix.wav"), "-map", "0:v", "-map", "1:a"] + enc +
                       ["-c:a", "aac", "-b:a", "192k", "-af", "apad", "-ar", "48000", "-t", f"{n / FPS:.3f}", "-movflags", "+faststart", str(dst)])
            ff = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "image2pipe", "-framerate", str(FPS), "-c:v", "mjpeg"] + cmd,
                                  stdin=subprocess.PIPE)
            st = time.time()
            for i in range(a, z):
                await pg.evaluate(f"window.FILM.frame({i / FPS:.5f})")
                if moved:
                    ff.stdin.close(); ff.wait(); dst.unlink(missing_ok=True)
                    sys.exit(f"the page was reloaded at frame {i}; nothing written")
                ff.stdin.write(await pg.screenshot(type="jpeg", quality=94))
                if (i - a) % 300 == 0 and i > a:
                    el = time.time() - st
                    print(f"{i - a}/{z - a} t={i / FPS:.0f}s  {el / 60:.1f} min, eta {(el / (i - a)) * (z - i) / 60:.1f} min", flush=True)
            ff.stdin.close()
            ff.wait()
            print("written", dst, flush=True)
        await b.close()


def parallel(film, m):
    """m renderers at once, then the parts joined and the soundtrack laid under them"""
    out = HERE / "build" / film
    procs = [subprocess.Popen([sys.executable, __file__, film, "part", str(k), str(m)],
                              stdout=open(out / f"part{k:02d}.log", "w"), stderr=subprocess.STDOUT) for k in range(m)]
    bad = [k for k, p in enumerate(procs) if p.wait() != 0]
    if bad:
        sys.exit(f"parts {bad} failed; see build/{film}/part*.log")
    join(film, m)


def join(film, m):
    """the m parts joined, the soundtrack under them (after re-rendering some parts alone)"""
    out = HERE / "build" / film
    (out / "parts" / "list.txt").write_text("".join(f"file 'p{k:02d}.mp4'\n" for k in range(m)))
    dst = out / f"TripleSpec_Reducer_{film}.mp4"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", str(out / "parts" / "list.txt"),
                    "-i", str(out / "mix.wav"), "-map", "0:v", "-map", "1:a", "-c:v", "copy", "-c:a", "aac", "-b:a", "192k",
                    "-ar", "48000", "-shortest", "-movflags", "+faststart", str(dst)], check=True)
    print("written", dst, flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) > 1 and a[1] == "parallel":
        parallel(a[0], int(a[2]))
    elif len(a) > 1 and a[1] == "join":
        join(a[0], int(a[2]))
    elif len(a) > 1 and a[1] == "part":
        asyncio.run(main(a[0], "video", part=(int(a[2]), int(a[3]))))
    else:
        asyncio.run(main(a[0], a[1] if len(a) > 1 else "video", [float(x) for x in a[2:]]))
