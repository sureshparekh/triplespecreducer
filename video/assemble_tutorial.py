"""
Film 2's timeline: the recording (video/build/tutorial/rec/rec.json) and the
narration (video/build/tutorial/voice/) on one clock.

While a step computes, the recording plays in real time as long as the
narration about it lasts (and 1.5 s more); beyond that, the rest of the wait
is played faster, at most 6 s on screen, with a badge saying by how much.

Writes video/film/public/tutorial/timeline.json (frames, cursor, clicks,
zooms, captions, badges, steps) and video/build/tutorial/placements.json
(where each narration clip starts, every click and keystroke) for audio.py.
"""

import bisect
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, "build", "tutorial")
rec = json.load(open(os.path.join(B, "rec", "rec.json")))
FPS = 30
MAX_WAIT = 6.0

# ---------------------------------------------------------------- the warp
beats = {b["id"]: b for b in rec["beats"]}
knots = [(0.0, 0.0)]            # (recording time, output time)
badges = []
segs = []
for w in rec["waits"]:
    b = beats[w["beat"]]
    keep_until = max(w["t0"], b["t0"] + b["lead"] + b["voice"] + 1.5)
    rest = w["t1"] - keep_until
    if rest > MAX_WAIT:
        segs.append((keep_until, w["t1"], MAX_WAIT / rest))
for a, b_, k in sorted(segs):
    r0, o0 = knots[-1]
    knots.append((a, o0 + (a - r0)))
    knots.append((b_, o0 + (a - r0) + (b_ - a) * k))
    badges.append(dict(t0=round(o0 + (a - r0), 3), t1=round(o0 + (a - r0) + (b_ - a) * k, 3), factor=round(1 / k)))
r0, o0 = knots[-1]
knots.append((rec["duration"], o0 + rec["duration"] - r0))
KR = [k[0] for k in knots]
KO = [k[1] for k in knots]


def out_t(r):
    i = max(0, min(len(KR) - 2, bisect.bisect_right(KR, r) - 1))
    u = (r - KR[i]) / max(1e-9, KR[i + 1] - KR[i])
    return KO[i] + u * (KO[i + 1] - KO[i])


def rec_t(o):
    i = max(0, min(len(KO) - 2, bisect.bisect_right(KO, o) - 1))
    u = (o - KO[i]) / max(1e-9, KO[i + 1] - KO[i])
    return KR[i] + u * (KR[i + 1] - KR[i])


# ---------------------------------------------------------------- frames
S = rec["samples"]
ST = [s[0] for s in S]
dur = KO[-1]
nf = int(dur * FPS)
imgs, img_id, frames, cursor = [], {}, [], []
for i in range(nf):
    r = rec_t(i / FPS)
    j = max(0, bisect.bisect_right(ST, r) - 1)
    s = S[j]
    if s[1] not in img_id:
        img_id[s[1]] = len(imgs)
        imgs.append(s[1])
    frames.append(img_id[s[1]])
    # the cursor between this sample and the next
    if j + 1 < len(S) and S[j + 1][2] >= 0 and s[2] >= 0:
        u = (r - s[0]) / max(1e-9, S[j + 1][0] - s[0])
        u = max(0.0, min(1.0, u))
        cursor.append([round(s[2] + u * (S[j + 1][2] - s[2]), 1), round(s[3] + u * (S[j + 1][3] - s[3]), 1)])
    else:
        cursor.append([s[2], s[3]])

# ---------------------------------------------------------------- captions and steps
scenes, placements = [], []
for b in rec["beats"]:
    v = json.load(open(os.path.join(B, "voice", b["id"] + ".json")))
    start = out_t(b["t0"] + b["lead"])
    lines = [dict(shown=l["text"], t0=round(start + l["t0"], 3), t1=round(start + l["t1"], 3),
                  words=[[x[0], round(start + x[1], 3), round(start + x[2], 3)] for x in l["words"]])
             for l in v["lines"]]
    scenes.append(dict(id=b["id"], step=b["step"], t0=round(out_t(b["t0"]), 3), t1=round(out_t(b["t1"]), 3), lines=lines))
    placements.append(dict(id=b["id"], t=round(start, 3)))

tl = dict(title="Hands-on tutorial: NGC 134 from scratch", duration=round(dur, 3), fps=FPS, size=rec["size"],
          images=imgs, frames=frames, cursor=cursor,
          clicks=[round(out_t(c), 3) for c in rec["clicks"]],
          zooms=[dict(t=round(out_t(z["t"]), 3), rect=z["rect"], full=z["full"]) for z in rec["zooms"]],
          badges=badges, scenes=scenes)
os.makedirs(os.path.join(HERE, "film", "public", "tutorial"), exist_ok=True)
json.dump(tl, open(os.path.join(HERE, "film", "public", "tutorial", "timeline.json"), "w"))
json.dump(dict(duration=round(dur, 3), voice=placements, clicks=tl["clicks"],
               keys=[round(out_t(k), 3) for k in rec["keys"]]),
          open(os.path.join(B, "placements.json"), "w"), indent=1)
print(f"output {dur / 60:.1f} min ({nf} frames), recording {rec['duration'] / 60:.1f} min, "
      f"{len(badges)} sped-up waits: {[b['factor'] for b in badges]}")
