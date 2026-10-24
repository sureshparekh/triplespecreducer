# Wavelength calibration

*The same step in Spextool: xspextool's wavelength calibration, as modified for TripleSpec in
the TS4_SOAR module; the cross-dispersed solution is described in Cushing, Vacca & Rayner
(2004).*

## Why

Each order is curved on the detector, and the image of the slit is tilted within it, by a
different amount in each order: a pixel's wavelength depends on both its column and its row.
The calibration gives **every pixel two numbers** — its wavelength and its position along the
slit — so that extraction can work on straight lines of constant wavelength.

## Sky lines or arcs

The upper atmosphere emits hundreds of narrow lines of the hydroxyl radical (OH), at
wavelengths known to a fraction of an ångström. They are in every target frame, recorded at
the same moment and through the same optics as the target, so they calibrate it exactly —
including any flexure since the start of the night. The CuHeAr arc lamp gives fewer lines,
and its frames are usually taken hours earlier. The tutorial, and the default here, use the
**sky**; the arcs remain an option (and are [checked against IDL](../validation/index.md) too).

```{themed} science/oh-lines
:alt: The OH sky lines in order 4

The OH sky lines in order 4
```

## The algorithm

1. **The image.** The target's frames are median-combined (the target itself is rejected by
    the median, since it moves between A and B) and divided by the flat.

2. **One-dimensional spectra.** In every order, column by column, the median of the central
    ±5″ of the slit gives a sky spectrum. (Spextool takes it from the *first* frame only — it
    indexes the stack of frames with two subscripts — and so does this program; the 2D
    tracing below uses the median image.)

3. **The offset.** The spectrum of order 4 is cross-correlated with a stored template, made on
    a reference night; a Gaussian fit to the correlation peak gives the shift, a fraction of a
    pixel to a few pixels.

4. **The lines.** Every line of the line list is fitted with a Gaussian where the shifted
    template predicts it; a line is kept if its centre lands within a pixel of the
    prediction.

5. **The 1D solution.** All the lines of all the orders are fitted together. The grating
    equation makes $m\lambda$ a smooth function of the column $x$ alone, so one polynomial
    serves all five orders:

    $$
    \frac{\lambda\,m}{m_{\rm home}} \;=\; \sum_{i=0}^{3}\sum_{j=0}^{3} a_{ij}\, x^i\, m^j ,
    $$
    
    with $m$ the order number and $m_{\rm home}=4$; points further than 3.5σ are rejected
    and the fit repeated. A typical RMS is **0.3 Å**.

6. **Along the slit.** Every isolated line is followed up and down the slit, every 3 rows,
    in 5-row medians; its centre at each row shows how the line is tilted and curved. Two
    surfaces per order are fitted to these points: the column each pixel would have at
    mid-slit (degree 2 × 2) and its position along the slit (degree 4 × 4). Evaluating the 1D
    solution at the mid-slit column gives each pixel's wavelength.

```{themed} science/wavecal-image
:alt: The wavelength image

The wavelength image
```

## The file

`cal/<source>_wavecal.fits`: the wavelength image and the slit-position image, with the 1D
and 2D coefficients and each order's dispersion in the header — the dispersions go on into
every spectrum file and are used by the telluric step to resample its kernel.

## Which solution for which target

Each science target with long exposures is calibrated on its own sky; a second pointing uses
its parent's solution; each standard (seconds-long exposures, too short for sky lines) uses
the solution of the science target nearest in time.
