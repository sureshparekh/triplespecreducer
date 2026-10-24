"""
A film's soundtrack, as the reels are mixed: the narration through the
broadcast voice chain, a quiet pad under it (ducked while the voice speaks),
a soft air at each scene change and, in the tutorial, a tick for every click
and keystroke of the recording.

    python video/audio.py explainer      # build/explainer/placements.json -> build/explainer/mix.wav
    python video/audio.py tutorial

The bed has no drums: these are long films to listen to while reading.
"""

import json
import sys
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import butter, fftconvolve, resample_poly, sosfilt

HERE = Path(__file__).parent
SR = 48000
rng = np.random.default_rng(7)


def lp(x, f, o=2): return sosfilt(butter(o, f, "low", fs=SR, output="sos"), x, axis=0)
def hp(x, f, o=2): return sosfilt(butter(o, f, "high", fs=SR, output="sos"), x, axis=0)
def bp(x, a, b, o=2): return sosfilt(butter(o, [a, b], "band", fs=SR, output="sos"), x, axis=0)
def ss(x): x = np.clip(x, 0, 1); return x * x * (3 - 2 * x)
def mtof(m): return 440.0 * 2 ** ((m - 69) / 12)


def add(buf, sig, t, pan=0.0, gain=1.0):
    i = int(t * SR)
    if sig.ndim == 1:
        l, r = np.cos((pan + 1) * np.pi / 4), np.sin((pan + 1) * np.pi / 4)
        sig = np.stack([sig * l * 1.414, sig * r * 1.414], 1)
    if i < 0: sig, i = sig[-i:], 0
    n = min(len(sig), len(buf) - i)
    if n > 0: buf[i:i + n] += sig[:n] * gain


def reverb_ir(sec=4.0, decay=1.1, damp=5000, seed=1):
    r = np.random.default_rng(seed); n = int(sec * SR); t = np.arange(n) / SR
    ir = np.stack([lp(r.standard_normal(n), damp) * np.exp(-t / decay * 3.0) for _ in range(2)], 1)
    ir[: int(0.018 * SR)] = 0
    return ir / np.sqrt((ir ** 2).sum(0, keepdims=True))


IR = reverb_ir()


def verb(x, wet):
    y = np.stack([fftconvolve(x[:, c], IR[:, c])[: len(x)] for c in range(2)], 1)
    return x * (1 - wet) + y * wet


def env_follow(x, att=0.03, rel=0.35):
    a = np.abs(x).astype(np.float64)
    win = int(0.04 * SR); a = np.convolve(a, np.ones(win) / win, mode="same")
    a = resample_poly(a, 1, 480); out = np.zeros_like(a); r = SR / 480
    ka, kr = np.exp(-1 / (att * r)), np.exp(-1 / (rel * r)); v = 0
    for i, s in enumerate(a):
        v = ka * v + (1 - ka) * s if s > v else kr * v + (1 - kr) * s; out[i] = v
    return np.interp(np.arange(len(x)) / 480, np.arange(len(out)), out)


def voice_fx(v):
    """broadcast chain: low cut, presence lift, gentle compression, soft limiter"""
    v = hp(v, 85, 2)
    v = v + 0.32 * bp(v, 2200, 5500, 2) + 0.12 * bp(v, 140, 260, 2)
    e = env_follow(np.abs(v), att=0.005, rel=0.12) + 1e-6
    thr = np.percentile(e[e > 1e-3], 70)
    g = np.where(e > thr, (e / thr) ** (-0.45), 1.0)
    v = v * g
    v = np.tanh(1.4 * v / (np.abs(v).max() + 1e-9)) / np.tanh(1.4)
    return v * 0.9


# ---------------------------------------------------------------- the bed
CHORDS = {"D": [50, 57, 61, 64, 66], "Bm": [47, 54, 57, 61, 62], "G": [43, 50, 54, 57, 59], "A": [45, 52, 57, 59, 61],
          "Em": [40, 52, 55, 59, 62]}
PROG = ["D", "Bm", "G", "A", "D", "Em", "G", "A"]


def pad_note(f, dur, att=2.6, rel=3.2):
    n = int((dur + rel) * SR); t = np.arange(n) / SR
    env = np.minimum(1, t / att) ** 1.6 * np.where(t > dur, np.exp(-(t - dur) / rel * 3), 1)
    out = np.zeros((n, 2))
    for k, det in enumerate([-0.07, 0.0, 0.07]):
        ff = f * 2 ** (det / 12)
        for h in range(1, 7):
            if ff * h > 6000: break
            ph = rng.uniform(0, 2 * np.pi); vib = 0.0015 * np.sin(2 * np.pi * (0.13 + 0.05 * k) * t + ph)
            s = np.sin(2 * np.pi * ff * h * t * (1 + vib) + ph) / h ** 1.7
            pan = (k - 1) * 0.55
            out[:, 0] += s * np.cos((pan + 1) * np.pi / 4); out[:, 1] += s * np.sin((pan + 1) * np.pi / 4)
    return out * env[:, None]


def pluck(f, dur=4.0):
    n = int(dur * SR); t = np.arange(n) / SR; s = np.zeros(n)
    for h, a in [(1, 1), (2, 0.42), (3, 0.16), (4, 0.09)]:
        s += a * 0.9 ** (h - 1) * np.sin(2 * np.pi * f * h * t) * np.exp(-t * (1.1 + 0.9 * h))
    return s * np.minimum(1, t / 0.006)


def bed(N, dur):
    """slow pads, a chord every 8 s, a sparse pluck now and then"""
    pad = np.zeros((N, 2)); arp = np.zeros((N, 2))
    cl = 8.0
    for c in range(int(np.ceil(dur / cl)) + 1):
        notes = CHORDS[PROG[c % len(PROG)]]; ct = c * cl
        for m in notes: add(pad, pad_note(mtof(m), cl + 0.6, att=2.4, rel=3.0), ct - 0.6, gain=0.11 / np.sqrt(len(notes)))
        tones = [m + 12 for m in notes[1:]]
        for k in range(4):
            if rng.random() < 0.6:
                add(arp, pluck(mtof(tones[rng.integers(len(tones))])), ct + k * 2.0 + rng.normal(0, .01), pan=rng.uniform(-.5, .5), gain=0.04)
    mus = verb(lp(pad, 3000), 0.45) + verb(hp(arp, 250), 0.4)
    t = np.arange(N) / SR
    return mus * (ss(t / 2.0) * (1 - ss((t - (dur - 3.0)) / 3.5)))[:, None]


# ---------------------------------------------------------------- effects
def pink(n):
    w = rng.standard_normal((n, 2))
    return lp(w, 900, 1) * 1.6 + lp(w, 4000, 1) * 0.6 + w * 0.15


def whoosh_soft(dur=1.4, pan0=-0.6, pan1=0.6):
    n = int(dur * SR); t = np.arange(n) / SR; a = t / dur
    env = np.sin(np.pi * np.clip(a, 0, 1)) ** 2.4
    x = pink(n)
    y = (bp(x, 200, 1200, 2) * (1 - a)[:, None] + bp(x, 1200, 5500, 2) * a[:, None]) * env[:, None]
    pan = pan0 + (pan1 - pan0) * a; L, R = np.cos((pan + 1) * np.pi / 4), np.sin((pan + 1) * np.pi / 4)
    y = np.stack([y[:, 0] * L, y[:, 1] * R], 1)
    return y / (np.abs(y).max() + 1e-9)


def click():
    n = int(0.05 * SR); t = np.arange(n) / SR
    s = hp(rng.standard_normal(n), 2500, 2) * np.exp(-t * 320) + np.sin(2 * np.pi * 2100 * t) * np.exp(-t * 220) * 0.5
    return s / np.abs(s).max()


def key():
    n = int(0.03 * SR); t = np.arange(n) / SR
    s = bp(rng.standard_normal(n), 1800, 7000, 2) * np.exp(-t * 500)
    return s / np.abs(s).max()


def main(film):
    B = HERE / "build" / film
    P = json.load(open(B / "placements.json"))
    dur = P["duration"] + 0.3
    N = int((dur + 2) * SR)
    v = np.zeros(N)
    for c in P["voice"]:
        w, sr = sf.read(B / "voice" / f"{c['id']}.wav")
        if w.ndim > 1: w = w.mean(1)
        if sr != SR: w = resample_poly(w, SR // 1000, sr // 1000)
        i = int(c["t"] * SR); n = min(len(w), N - i)
        v[i:i + n] += w[:n]
    v = voice_fx(v)
    print("bed...", flush=True)
    mus = bed(N, dur)
    fx = np.zeros((N, 2))
    for c in P["voice"][1:]:       # a breath of air as each scene comes in
        add(fx, whoosh_soft(1.4, rng.uniform(-.7, -.3), rng.uniform(.3, .7)), c["t"] - 1.4, gain=0.03)
    for t in P.get("clicks", []): add(fx, click(), t, pan=rng.uniform(-.15, .15), gain=0.05)
    for t in P.get("keys", []): add(fx, key(), t, pan=rng.uniform(-.1, .1), gain=0.018)
    fx = verb(fx, 0.12)
    e = env_follow(v); act = ss((e / (np.percentile(e[e > 1e-4], 60) + 1e-9) - 0.15) / 0.5)
    mus *= 0.045 / np.sqrt((mus ** 2).mean())
    mix = np.stack([v, v], 1) * 0.95 + mus * (1 - 0.45 * act)[:, None] + fx
    mix = mix[: int(dur * SR)]
    sf.write(B / "mix.wav", (mix / max(1.0, np.abs(mix).max() / 0.98)).astype(np.float32), SR)
    print("written", B / "mix.wav", f"{len(mix) / SR / 60:.1f} min, peak {np.abs(mix).max():.2f}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
