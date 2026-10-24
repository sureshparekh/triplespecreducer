# 7. Cleaning and saving

:::{admonition} Why
:class: note
Optional last touches. Smoothing lowers the noise at the cost of resolution; a bad region
— the noise in the opaque gaps between the bands, or a residual from a bad pixel — can be
bridged by a straight line or blanked.
:::

## What you do

1. Go to **7 Clean & save** and choose **NGC134**.

```{shots} 07-clean
:alt: The final spectrum and its signal-to-noise

The final spectrum and its signal-to-noise
```

The stage shows the merged spectrum (top) and its **signal-to-noise ratio** (bottom): about
100–200 per pixel in the J, H and K windows, dropping in the opaque gaps.

2. **Smoothing** — choose *Savitzky-Golay* (a moving polynomial that keeps line shapes) or
    *Gaussian*, with a width in pixels, or *None*.
3. **Fix** or **Remove** a region — drag across it on the plot, then press *Fix* (a straight
    line between its ends) or *Remove* (blanked). *Undo all* reverts the repairs.
4. Press **Smooth and save**. The result is `proc/NGC134_clean.fits`; **Also write a text table**
    saves it as columns of wavelength, flux, error and flag.

:::{admonition} Smoothing and resolution
:class: warning
Smoothing widens every line. Keep the width below the width of the lines you will measure
— or measure on the unsmoothed `NGC134_merged.fits`. The IDL panel's default
Savitzky–Golay width, 1.5 slit widths (3.6 pixels), truncates to 3 pixels, which with
degree 2 leaves the spectrum unchanged; the default here does the same.
:::

**Next:** [8. The products](8-products.md) · *The science:* [cleaning](../science/merge.md#cleaning)
