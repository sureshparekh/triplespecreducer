"""
Put the rendered films on the documentation site: a web encode of each
(GitHub keeps files under 100 MB), WebVTT subtitles from the timeline's
caption lines, a poster frame, and the chapter list for the page.

    python video/publish.py explainer 12.0    # the poster at 12 s
    python video/publish.py tutorial 160.0
"""

import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
DOCS = HERE.parent / "docs" / "_static" / "video"
NAME = {"explainer": "science", "tutorial": "tutorial"}


def stamp(t):
    h, r = divmod(max(0.0, t), 3600)
    m, s = divmod(r, 60)
    return f"{int(h):02d}:{int(m):02d}:{s:06.3f}"


def main(film, poster_t):
    tl = json.load(open(HERE / "film" / "public" / film / "timeline.json"))
    src = HERE / "build" / film / f"TripleSpec_Reducer_{film}.mp4"
    DOCS.mkdir(parents=True, exist_ok=True)
    name = NAME[film]
    # the subtitles: one cue per caption line, as shown on screen
    cues = ["WEBVTT", ""]
    for s in tl["scenes"]:
        for l in s["lines"]:
            cues += [f"{stamp(l['t0'])} --> {stamp(l['t1'] + .3)}", l["shown"], ""]
    (DOCS / f"{name}.vtt").write_text("\n".join(cues))
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(src), "-c:v", "libx264", "-preset", "slow",
                    "-crf", "27", "-tune", "stillimage" if film == "tutorial" else "animation", "-pix_fmt", "yuv420p",
                    "-c:a", "aac", "-b:a", "128k", "-movflags", "+faststart", str(DOCS / f"{name}.mp4")], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-ss", str(poster_t), "-i", str(src), "-frames:v", "1",
                    "-vf", "scale=1280:-2", "-q:v", "3", str(DOCS / f"{name}.jpg")], check=True)
    mb = (DOCS / f"{name}.mp4").stat().st_size / 1e6
    print(f"{name}.mp4 {mb:.0f} MB, {tl['duration'] / 60:.1f} min")
    for s in tl["scenes"]:
        m, sec = divmod(int(s["t0"]), 60)
        print(f"| {m}:{sec:02d} | {s['id']} |")


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]))
