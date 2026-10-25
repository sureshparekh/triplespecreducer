# 7. Cleaning and saving

:::{admonition} Why
:class: note
Optional last touches. Smoothing lowers the noise at the cost of resolution; a bad region
— the noise in the opaque gaps between the bands, or a residual from a bad pixel — can be
bridged by a straight line or blanked.
:::

## What you do

1. Go to **7 Clean & save** and choose **NGC134**. The **File** is its merged spectrum; any
    other of its spectra (corrected or combined, any order and aperture) can be opened too.

```{shots} 07-clean
:alt: The final spectrum and its signal-to-noise

The final spectrum and its signal-to-noise
```

The stage shows the spectrum (top) and its **signal-to-noise ratio** (bottom): about
100–200 per pixel in the J, H and K windows, dropping in the opaque gaps.

2. **Fix** or **Remove** a region: drag across it on the plot, then press *Fix* (a straight
    line between its ends) or *Remove* (blanked), for the flux and its error or, with *Edit:
    Error*, for the error only. As in the IDL tool, edits are drawn on a working copy:
    **Keep** makes them stick, **Undo** goes back to the last Keep.
3. **Smoothing**: choose *Savitzky-Golay* (a moving polynomial that keeps line shapes) or
    *Gaussian*, with a width in pixels, and press **Smooth**. It starts from the kept
    spectrum, and is allowed once per file opened, as in the IDL tool.
4. Press **Save**: the kept spectrum is written to `proc/NGC134_clean.fits`. **Also write a
    text table** writes it as the IDL tool's SPEX2TEXT does, header and every order included.

:::{admonition} Smoothing and resolution
:class: warning
Smoothing widens every line. Keep the width below the width of the lines you will measure
— or measure on the unsmoothed `NGC134_merged.fits`. The IDL panel's default
Savitzky–Golay width, 1.5 slit widths (3.6 pixels), truncates to 3 pixels, which with
degree 2 leaves the spectrum unchanged; the default here does the same.
:::

**Next:** [8. The products](8-products.md) · *The science:* [cleaning](../science/merge.md#cleaning)
