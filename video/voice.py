"""
The narration, in the voice of the reels: Chatterbox Turbo cloned from the
Kokoro am_michael reference, one generation per scene, the pauses tightened,
then MMS forced alignment for every word's time (captions and action cues).

    conda run -n tts python video/voice.py tutorial      # or: explainer
    ... [scene ids]                                       # regenerate only these

Writes video/build/<film>/voice/<id>.wav (24 kHz) and <id>.json: the lines
as displayed, each with its words' start and end times in seconds from the
start of that scene's audio.
"""

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
import torch
import torchaudio
from chatterbox.tts_turbo import ChatterboxTurboTTS

HERE = Path(__file__).parent
FILM = sys.argv[1]
ONLY = set(sys.argv[2:])
REF = "/home/suresh/PhD/content/voice_samples/ref_am_michael.wav"
OUT = HERE / "build" / FILM / "voice"
OUT.mkdir(parents=True, exist_ok=True)

# how the voice should say what the captions write
SAY = [(r"\bNGC 134\b", "N G C one thirty-four"), (r"\bHD 225200\b", "H D two two five two zero zero"),
       (r"\bHD225200\b", "H D two two five two zero zero"), (r"\bTSR\b", "T S R"), (r"\bA0 V\b", "A zero V"),
       (r"\bA1 V\b", "A one V"), (r"\bRMS\b", "R M S"), (r"\bSIMBAD\b", "Simbad"), (r"\bIDL\b", "I D L"),
       (r"\bFITS\b", "fits"), (r"\bOH\b", "O H"), (r"\bA − B\b", "A minus B"), (r"\bµm\b", "microns"),
       (r"\b2025\b", "twenty twenty-five"), (r"\bTripleSpec\b", "Triple Spec"), (r"\bSpextool\b", "Spex tool"),
       (r"\bxtellcor\b", "X tell cor"), (r"\bxspextool\b", "X spex tool"), (r"\bSOAR\b", "soar")]


def spoken(text):
    for pat, rep in SAY:
        text = re.sub(pat, rep, text)
    return text


def norm(w):
    return re.sub(r"[^a-z']", "", w.lower().replace("-", ""))


def load_scenes():
    name = {"tutorial": "tutorial_script", "explainer": "explainer_script"}[FILM]
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.BEATS if FILM == "tutorial" else m.SCENES


tts = ChatterboxTurboTTS.from_pretrained(device="cuda")
SR = tts.sr
B = torchaudio.pipelines.MMS_FA
fa, tok, aligner = B.get_model().cuda(), B.get_tokenizer(), B.get_aligner()


def trim(a):
    nz = np.flatnonzero(np.abs(a) > 0.01 * np.abs(a).max())
    return a[max(0, nz[0] - int(0.03 * SR)): nz[-1] + int(0.1 * SR)] if len(nz) else a


def tighten(a, cap=0.32):
    """Shorten every internal pause to at most `cap` s: a tutorial's pace, unhurried but without dead air."""
    n = int(0.02 * SR)
    e = np.array([np.sqrt((a[i:i + n] ** 2).mean()) for i in range(0, len(a) - n, n)])
    quiet = e < 0.03 * e.max()
    out, i, k = [], 0, 0
    while k < len(quiet):
        if quiet[k]:
            j = k
            while j < len(quiet) and quiet[j]:
                j += 1
            if (j - k) * 0.02 > cap and k > 0 and j < len(quiet):
                keep = int(cap / 2 * SR)
                out.append(a[i:k * n + keep])
                i = j * n - keep
            k = j
        else:
            k += 1
    out.append(a[i:])
    return np.concatenate(out)


def align(a, words):
    x = torchaudio.functional.resample(torch.tensor(a, dtype=torch.float32)[None], SR, B.sample_rate).cuda()
    keep = [i for i, w in enumerate(words) if norm(w)]
    with torch.inference_mode():
        em, _ = fa(x)
        spans = aligner(em[0], tok([norm(words[i]) for i in keep]))
    r = x.shape[1] / em.shape[1] / B.sample_rate
    out, score = [None] * len(words), []
    for i, sp in zip(keep, spans):
        out[i] = (sp[0].start * r, sp[-1].end * r)
        score.append(np.mean([s.score for s in sp]))
    for i in range(len(words)):
        if out[i] is None:
            p = out[i - 1] if i else (0, 0)
            out[i] = (p[1], p[1])
    return out, float(np.mean(score)) if score else 0.0


total = 0.0
for sc in load_scenes():
    shown = sc["lines"]
    said = [spoken(l) for l in shown]
    text = " ".join(said)
    words = text.split()
    wav, js = OUT / f"{sc['id']}.wav", OUT / f"{sc['id']}.json"
    if wav.exists() and js.exists() and json.loads(js.read_text()).get("text") == text and sc["id"] not in ONLY:
        a, _ = sf.read(wav)
    else:
        best = None
        for seed in (11, 23, 37):            # retry when the alignment says words were mangled or skipped
            torch.manual_seed(seed)
            a = tighten(trim(tts.generate(text, audio_prompt_path=REF, temperature=0.7).squeeze().cpu().numpy()))
            _, s_ = align(a, words)
            cps = len(text) / (len(a) / SR)
            print(f"  {sc['id']} seed {seed}: {len(a) / SR:5.1f}s align {s_:.2f} chars/s {cps:.1f}", flush=True)
            if best is None or s_ > best[1]:
                best = (a, s_)
            if s_ > 0.55 and 9 < cps < 19:
                break
        a = best[0]
        sf.write(wav, a.astype(np.float32), SR)
    spans, score = align(a, words)
    lines, wi = [], 0
    for disp, sp in zip(shown, said):
        n = len(sp.split())
        ws = [[words[wi + k], round(spans[wi + k][0], 3), round(spans[wi + k][1], 3)] for k in range(n)]
        lines.append(dict(text=disp, t0=ws[0][1], t1=ws[-1][2], words=ws))
        wi += n
    js.write_text(json.dumps(dict(id=sc["id"], text=text, sr=SR, duration=round(len(a) / SR, 3),
                                  align=round(score, 3), lines=lines), indent=1))
    total += len(a) / SR
    print(f"{sc['id']:12s} {len(a) / SR:6.1f} s  align {score:.2f}", flush=True)
print(f"total {total / 60:.1f} min")
