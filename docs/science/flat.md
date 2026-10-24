# The flat field

*The same step in Spextool: xspextool's flat-field calibration, as modified for TripleSpec in
the TS4_SOAR module; described in Cushing, Vacca & Rayner (2004).*

## Why

No two pixels of an infrared array respond identically: their gains differ by a few per cent,
some are dead, some hot. If every frame is divided by an image of the detector's response to
even illumination — a **flat field** — those differences cancel. The same frames also show
where the five orders lie, which every later step needs.

## The frames

A quartz lamp lights the slit evenly. Two sets are taken: **lamps on** and **lamps off**. The
lamps-off frames contain everything except the lamp — the warm optics' thermal glow, which is
strong in the K band — and are subtracted.

## The algorithm

1. **Reading.** Every frame becomes a count rate (DN s⁻¹), with its variance:

    $$
    r = \frac{\rm DN}{{\rm COADDS}\cdot t_{\rm exp}},\qquad
    \sigma_r^2 = \frac{|{\rm DN}|\,\big(1-\frac{t_{\rm read}(n^2-1)}{3\,t\,n}\big)}{{\rm COADDS}^2\,t^2\,g}
    + \frac{2\,{\rm RON}^2}{n\;{\rm COADDS}\;t^2\,g^2},
    $$
    
    with $n$ the number of Fowler samples, $g$ the gain and RON the read noise
    (the noise model of Cushing, Vacca & Rayner (2004).

2. **Background.** The lamps-off frames are averaged and subtracted from every lamps-on frame.
3. **The orders.** On the average of the lamps-on frames, the program walks along each
    order from the guess position in the instrument file, every 3 columns, fitting a Gaussian
    to each of the two narrow bright *outrigger* peaks that border every TripleSpec order.
    Where a pair is found 98–122 pixels apart it is kept; each edge is fitted with a robust
    polynomial of degree 4, and moved 13 % of the slit height inward.

4. **Scaling.** Each lamps-on frame is scaled to the median level inside the orders, so lamp
    variations do not matter, and they are averaged.

5. **Normalisation.** Each order is divided by a smooth model of itself, so that only the
    pixel-to-pixel structure is left: the order is resampled onto a grid of (column,
    arcsec along the slit), median-filtered over 5 pixels, fitted with a surface of bicubic
    patches on a 5 × 3 grid, and the original pixels are divided by that surface
    (the normalisation of Cushing, Vacca & Rayner (2004).

```{themed} science/flat
:alt: A lamps-on flat and the normalised flat

A lamps-on flat and the normalised flat
```

The normalised flat is 1 on average inside every order, with the detector's structure;
outside the orders it is 1 exactly, so nothing there is changed.

## What it does not remove

The lamp's own spectrum and the instrument's large-scale response are divided out of the
flat by the normalisation, so they are **not** corrected here: they are left in every frame,
and removed together with the atmosphere by the [telluric correction](telluric.md), because
the standard star carries exactly the same ones.

## The file

`cal/flat.fits`, in Spextool's format: the flat, its variance and a flag image, with the
order edges, extraction ranges, slit sizes and plate scale in the header. Spextool passes
these numbers to the next steps through the header text; this program writes and reads them
exactly as IDL does (see [Behaviour of the IDL tool](../validation/idl_differences.md#header-values-text-then-read-back)).
