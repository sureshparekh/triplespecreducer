"""
The science pages' figures, made from the tutorial night with the package's
own routines (so every figure shows what the reduction actually does), in a
light and a dark version.

    python docs/make_figures.py          (after docs/make_shots.py has made reduced/tutorial)
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402
from astropy.io import fits                            # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from triplespec_reducer import flat as F, wavecal as W, extract as E, telluric as T   # noqa: E402
from triplespec_reducer.frames import read_ts4                                           # noqa: E402
from triplespec_reducer.instrument import read_atran                                     # noqa: E402

RED = os.path.join(ROOT, "reduced", "tutorial")
RAW = os.path.join(ROOT, "data", "raw", "2025-11-02")
OUT = os.path.join(ROOT, "docs", "_static", "science")
os.makedirs(OUT, exist_ok=True)
SCI, STD = "NGC134", "HD225200_1"
ORDER_COL = {3: "#C8892A", 4: "#2E8B57", 5: "#3A6FD8", 6: "#8E5BD0", 7: "#D0556F"}
STYLES = {
    "light": dict(bg="white", fg="#1B2230", grid="#D9DEE6", mute="#6B7585", img="gray"),
    "dark": dict(bg="#0E131B", fg="#D5DCE6", grid="#2A3342", mute="#8B95A3", img="gray"),
}


def raw(n):
    return os.path.join(RAW, f"SPEC_SO2025D-006-{n:04d}.fits")


def style(name):
    s = STYLES[name]
    plt.rcParams.update({
        "figure.facecolor": s["bg"], "axes.facecolor": s["bg"], "savefig.facecolor": s["bg"],
        "text.color": s["fg"], "axes.labelcolor": s["fg"], "axes.edgecolor": s["mute"],
        "xtick.color": s["mute"], "ytick.color": s["mute"], "grid.color": s["grid"],
        "axes.grid": True, "grid.linewidth": 0.6, "font.size": 10.5, "axes.titlesize": 11.5,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "figure.dpi": 110, "savefig.dpi": 170, "savefig.bbox": "tight"})
    return s


AXES = {}


def save(fig, name, theme):
    """Save, and record where each axes lies in the saved image (fractions, origin top left)."""
    fig.savefig(os.path.join(OUT, f"{name}-{theme}.png"))
    r = fig.canvas.get_renderer()
    tb = fig.get_tightbbox(r).padded(plt.rcParams["savefig.pad_inches"])
    W, H = tb.width, tb.height
    axes = []
    for ax in fig.axes:
        b = ax.get_window_extent(r).transformed(fig.dpi_scale_trans.inverted())
        axes.append(dict(x0=(b.x0 - tb.x0) / W, x1=(b.x1 - tb.x0) / W,
                         y0=1 - (b.y1 - tb.y0) / H, y1=1 - (b.y0 - tb.y0) / H,
                         xlim=list(ax.get_xlim()), ylim=list(ax.get_ylim())))
    AXES[name] = axes
    plt.close(fig)
    print("  ", f"{name}-{theme}.png")


# ----------------------------------------------------------------- the data, read once
fl = F.read_flat(os.path.join(RED, "cal", "flat.fits"))
wc = W.read_wavecal(os.path.join(RED, "cal", f"{SCI}_wavecal.fits"), rotation=fl["rotation"])
aw, at = read_atran(fl["rp"])
frames_ab, var_ab, _, _ = read_ts4([raw(116), raw(117)], rotation=fl["rotation"])
diff = (frames_ab[0] - frames_ab[1]) / fl["flat"]
sky_raw, _, _, _ = read_ts4([raw(121)], rotation=fl["rotation"])
sky = sky_raw[0] / fl["flat"]
prof = E.mkspatprof2d(diff, fl["omask"], wc["wavecal"], wc["spatcal"],
                      fl["slith_arc"] / fl["slith_pix"], 1, aw, at, 0.4)
pos, signs = E.findpeaks(prof, 2)


def spec_file(name):
    with fits.open(os.path.join(RED, "proc", name)) as h:
        d = h[0].data
        return d.astype(d.dtype.newbyteorder("=")), h[0].header


comb_std, h_std = spec_file(f"{STD}_comb.fits")
comb_sci, h_sci = spec_file(f"{SCI}_comb.fits")
tell, _ = spec_file(f"{SCI}_tellspec.fits")
tc, _ = spec_file(f"{SCI}_tc.fits")
merged, _ = spec_file(f"{SCI}_merged.fits")
orders = [int(t) for t in str(h_sci["ORDERS"]).split(",")]
st = T.load(os.path.join(RED, "proc", f"{STD}_comb.fits"), os.path.join(RED, "proc", f"{SCI}_comb.fits"),
            float(h_std.get("A0VBMAG", 6.40)) if "A0VBMAG" in h_std else 6.40, 6.39)
T.conkernel(st, order=6, normreg=[(0.990, 0.997), (1.013, 1.020)], normdegree=3, linereg=(0.997, 1.013))
specs = sorted(f for f in os.listdir(os.path.join(RED, "proc")) if f.startswith(f"{STD}_spec"))

for theme in ("light", "dark"):
    s = style(theme)
    print(theme)

    # 1. a raw frame: the five orders
    fig, ax = plt.subplots(figsize=(9, 4.6))
    img = frames_ab[0]
    lo, hi = np.nanpercentile(img, [5, 99.5])
    ax.imshow(img, origin="lower", cmap=s["img"], vmin=lo, vmax=hi, aspect="auto")
    x = np.arange(img.shape[1])
    for i, o in enumerate(fl["orders"]):
        for j in (0, 1):
            ax.plot(x, np.polyval(fl["edgecoeffs"][i, j][::-1], x), color=ORDER_COL[o], lw=0.9)
        yc = np.polyval(fl["edgecoeffs"][i, 1][::-1], 1700) + 8
        ax.text(1700, yc, f"order {o}", color=ORDER_COL[o], fontsize=10, fontweight="bold")
    ax.set(xlabel="column (dispersion)", ylabel="row", title="One raw frame of HD 225200, with the order edges")
    ax.grid(False)
    save(fig, "raw-orders", theme)

    # 2. the atmosphere and the five orders
    fig, ax = plt.subplots(figsize=(9, 3.4))
    ax.plot(aw, at, color=s["mute"], lw=0.5)
    for k, o in enumerate(orders):
        wv = comb_sci[k, 0]
        g = np.isfinite(wv)
        y0 = -0.12 - 0.14 * (k % 2)
        ax.plot([np.nanmin(wv[g]), np.nanmax(wv[g])], [y0, y0], color=ORDER_COL[o], lw=6,
                solid_capstyle="butt")
        ax.text(np.nanmean(wv[g]), y0 + 0.05, f"order {o}", color=ORDER_COL[o], ha="center", fontsize=8.5)
    for name, (a, b) in (("J", (1.17, 1.33)), ("H", (1.49, 1.78)), ("K", (2.03, 2.37))):
        ax.text((a + b) / 2, 0.45, name, ha="center", fontsize=13, fontweight="bold", color=s["fg"], alpha=0.6)
    ax.set(xlim=(0.92, 2.5), ylim=(-0.32, 1.05), xlabel="wavelength (µm)", ylabel="transmission",
           title="What the atmosphere lets through (ATRAN, R = 4000), and the five TripleSpec orders")
    save(fig, "atmosphere", theme)

    # 3. the flat: raw lamp frame and the normalised flat
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    with fits.open(os.path.join(RED, "cal", "flat.fits")) as h:
        nf = fl["flat"]
    lamp, _, _, _ = read_ts4([os.path.join(RAW, "SPEC_DFLAT_0002.fits")], rotation=fl["rotation"])
    sl = (slice(150, 560), slice(800, 1250))
    axs[0].imshow(lamp[0][sl], origin="lower", cmap=s["img"], aspect="auto")
    axs[0].set_title("a lamps-on flat (detail)")
    axs[1].imshow(nf[sl], origin="lower", cmap=s["img"], vmin=0.9, vmax=1.1, aspect="auto")
    axs[1].set_title("the normalised flat: 1 ± pixel-to-pixel variations")
    for a in axs:
        a.grid(False)
        a.set_xticks([])
        a.set_yticks([])
    save(fig, "flat", theme)

    # 4. sky lines: the OH spectrum of order 4, from the sky frame and the wavelength image
    fig, ax = plt.subplots(figsize=(9, 3.4))
    o = 4
    m = (fl["omask"] == o) & (wc["spatcal"] > 14) & (wc["spatcal"] < 24)
    cols = np.unique(np.nonzero(m)[1])
    wl, fx = [], []
    for c in cols:
        r = m[:, c]
        wl.append(np.nanmedian(wc["wavecal"][r, c]))
        fx.append(np.nanmedian(sky[r, c]))
    ax.plot(wl, fx, color=ORDER_COL[o], lw=0.7)
    ax.set(xlabel="wavelength (µm)", ylabel="DN s⁻¹",
           title="The night sky in order 4 (H band): OH emission lines at known wavelengths")
    save(fig, "oh-lines", theme)

    # 5. the wavelength image
    fig, ax = plt.subplots(figsize=(9, 4.2))
    im = ax.imshow(wc["wavecal"], origin="lower", cmap="viridis", aspect="auto")
    cb = fig.colorbar(im, ax=ax, pad=0.01)
    cb.set_label("wavelength (µm)")
    ax.set(xlabel="column", ylabel="row", title="The wavelength of every pixel (the 2D solution)")
    ax.grid(False)
    save(fig, "wavecal-image", theme)

    # 6. A, B and A-B
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    sl = (slice(370, 560), slice(900, 1300))
    for a, (im_, t) in zip(axs, ((frames_ab[0], "A"), (frames_ab[1], "B"), (diff, "A − B, divided by the flat"))):
        v = np.nanpercentile(np.abs(diff[sl]), 99)
        if t == "A" or t == "B":
            lo, hi = np.nanpercentile(im_[sl], [2, 99.5])
            a.imshow(im_[sl], origin="lower", cmap=s["img"], vmin=lo, vmax=hi, aspect="auto")
        else:
            a.imshow(im_[sl], origin="lower", cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
        a.set_title(t)
        a.grid(False)
        a.set_xticks([])
        a.set_yticks([])
    save(fig, "a-minus-b", theme)

    # 7. the spatial profile of order 5 with apertures and background
    fig, ax = plt.subplots(figsize=(8, 3.4))
    y, p = prof[5]
    ax.plot(y, p, color=s["fg"], lw=1.1)
    i5 = fl["orders"].index(5)
    for j in range(2):
        c = pos[i5, j]
        ax.axvspan(c - 1.0, c + 1.0, color="#2E8B57", alpha=0.25)
        for side in (-1, 1):
            a0, a1 = sorted((c + side * 2.5, c + side * 7.5))
            ax.axvspan(a0, a1, color="#D0556F", alpha=0.12)
        ax.text(c, np.nanmax(p) * 1.05, "+" if signs[j] > 0 else "−", ha="center", fontsize=14, color="#2E8B57")
    ax.set(xlabel="position along the slit (arcsec)", ylabel="normalised flux",
           title="Spatial profile, order 5: apertures (green) and background regions (red)")
    save(fig, "profile", theme)

    # 8. combining: every scaled spectrum of order 5 and the combination
    fig, ax = plt.subplots(figsize=(9, 3.4))
    k5 = orders.index(5)
    for n, f in enumerate(specs):
        d, _ = spec_file(f)
        ax.plot(d[k5, 0], d[k5, 1], lw=0.9, alpha=0.8, label=f.split("_")[-1].replace(".fits", ""),
                color=["#9AA6B8", "#B7A27C", "#8FB59A", "#B79AB8"][n % 4])
    ax.plot(comb_std[k5, 0], comb_std[k5, 1], color=s["fg"], lw=1.4, label="combined")
    ax.set_xlim(1.20, 1.27)
    z = (comb_std[k5, 0] > 1.20) & (comb_std[k5, 0] < 1.27)
    ax.set_ylim(np.nanmin(comb_std[k5, 1][z]) * 0.9, np.nanmax(comb_std[k5, 1][z]) * 1.08)
    ax.legend(loc="lower left", ncol=5, fontsize=8.5)
    ax.set(xlabel="wavelength (µm)", ylabel="DN s⁻¹",
           title=f"HD 225200, order 5 (detail): the {len(specs)} spectra of its two pairs and their combination")
    save(fig, "combine", theme)

    # 9. the kernel fit at Pa delta
    kp = st["kernel_plot"]
    fig, axs = plt.subplots(2, 1, figsize=(8, 4.6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
    axs[0].plot(kp["w"], kp["data"], color=s["fg"], lw=2.2, alpha=0.6, label="the star (normalised)")
    axs[0].plot(kp["w"], kp["vega"], color="#D0556F", lw=1.0, label="Vega's model")
    axs[0].plot(kp["w"], kp["conv"], color="#2E8B57", lw=1.0, ls="--", label="Vega convolved with the kernel")
    axs[0].legend(loc="lower left", fontsize=9)
    axs[0].set(ylabel="flux − 1", title="Building the kernel on Pa δ (1.005 µm), order 6")
    axs[1].plot(kp["w"], 100 * kp["resid"], color=s["mute"], lw=0.8)
    axs[1].set(xlabel="wavelength (µm)", ylabel="residual (%)")
    save(fig, "kernel", theme)

    # 10. the telluric correction of order 4: standard, telluric spectrum, target before and after
    fig, axs = plt.subplots(3, 1, figsize=(9, 6.4), sharex=True)
    k4 = orders.index(4)
    axs[0].plot(comb_std[k4, 0], comb_std[k4, 1], color=s["fg"], lw=0.7)
    axs[0].set(ylabel="DN s⁻¹", title="Order 4 (H band)")
    axs[0].text(0.01, 0.85, "the A0 V standard, as observed", transform=axs[0].transAxes)
    axs[1].plot(tell[k4, 0], tell[k4, 1], color="#C8892A", lw=0.7)
    axs[1].text(0.01, 0.85, "the telluric spectrum = Vega model / observed standard", transform=axs[1].transAxes)
    axs[1].set(ylabel="flux / (DN s⁻¹)")
    for a_, w_, y_ in ((axs[0], comb_std[k4, 0], comb_std[k4, 1]), (axs[1], tell[k4, 0], tell[k4, 1])):
        inner = (w_ > 1.45) & (w_ < 1.78)                 # away from the order's ragged ends
        a_.set_ylim(0, np.nanmax(y_[inner]) * 1.25)
    ax2 = axs[2]
    ax2.plot(comb_sci[k4, 0], comb_sci[k4, 1] / np.nanmedian(comb_sci[k4, 1]), color=s["mute"], lw=0.6,
             label="NGC 134 as observed (scaled)")
    ax2.plot(tc[k4, 0], tc[k4, 1] / np.nanmedian(tc[k4, 1]), color="#3A6FD8", lw=0.8,
             label="corrected and flux calibrated (scaled)")
    ax2.legend(loc="upper left", fontsize=9)
    ax2.set_ylim(0, 2.2)
    axs[0].set_xlim(1.42, 1.83)
    ax2.set(xlabel="wavelength (µm)", ylabel="relative flux")
    save(fig, "telluric", theme)

    # 11. merging
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(10, 5.6), sharex=True)
    mm = merged if merged.ndim == 2 else merged[0]
    top = np.nanpercentile(mm[1], 98) * 1.4
    for k, o in enumerate(orders):
        ax.plot(tc[k, 0], tc[k, 1], color=ORDER_COL[o], lw=0.6, alpha=0.9, label=f"order {o}")
    ax.legend(ncol=5, fontsize=8.5, loc="upper right")
    ax.set(ylim=(0, top), ylabel="erg s⁻¹ cm⁻² Å⁻¹",
           title="NGC 134: the five corrected orders (top), merged into one spectrum (bottom)")
    bx.plot(mm[0], mm[1], color=s["fg"], lw=0.6)
    bx.set(ylim=(0, top), xlabel="wavelength (µm)", ylabel="erg s⁻¹ cm⁻² Å⁻¹")
    save(fig, "merge", theme)
import json                                               # noqa: E402
json.dump(AXES, open(os.path.join(OUT, "axes.json"), "w"), indent=1)
print("done")
