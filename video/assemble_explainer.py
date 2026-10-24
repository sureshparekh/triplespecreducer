"""
Film 1's timeline: its scenes back to back, each lasting its narration plus a
lead-in and a moment for the visuals to settle.

Writes video/film/public/explainer/timeline.json and
video/build/explainer/placements.json (for audio.py).
"""

import importlib.util
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, "build", "explainer")
spec = importlib.util.spec_from_file_location("explainer_script", os.path.join(HERE, "explainer_script.py"))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
LEAD, TAIL = 0.9, 1.4

t, scenes, placements = 0.0, [], []
for sc in S.SCENES:
    v = json.load(open(os.path.join(B, "voice", sc["id"] + ".json")))
    start = t + LEAD
    lines = [dict(shown=l["text"], t0=round(start + l["t0"], 3), t1=round(start + l["t1"], 3),
                  words=[[x[0], round(start + x[1], 3), round(start + x[2], 3)] for x in l["words"]]) for l in v["lines"]]
    end = start + v["duration"] + TAIL
    scenes.append(dict(id=sc["id"], t0=round(t, 3), t1=round(end, 3), lines=lines))
    placements.append(dict(id=sc["id"], t=round(start, 3)))
    t = end
json.dump(dict(title="The science and the steps", duration=round(t + 1.0, 3), scenes=scenes),
          open(os.path.join(HERE, "film", "public", "explainer", "timeline.json"), "w"))
json.dump(dict(duration=round(t + 1.0, 3), voice=placements, clicks=[], keys=[]),
          open(os.path.join(B, "placements.json"), "w"), indent=1)
print(f"{t / 60:.1f} min")
