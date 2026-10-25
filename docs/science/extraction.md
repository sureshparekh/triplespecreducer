# Extraction

*The same step in Spextool: xspextool's extraction ("Do All Steps"), described in Cushing,
Vacca & Rayner (2004) and the Spextool manual.*

## Sky subtraction by nodding

```{themed} science/a-minus-b
:alt: A, B and A − B

A, B and A − B
```

The target is observed at two positions, A and B, and the frames are subtracted. Everything
that is the same in both — the sky's lines and continuum, the telescope's thermal glow, the
detector's dark current — cancels:

$$
{\rm work} = \frac{A - B}{\rm flat}, \qquad \sigma^2_{\rm work} = \frac{\sigma_A^2 + \sigma_B^2}{{\rm flat}^2} .
$$

* **Point sources: ABBA along the slit.** Both positions are on the slit, about 18″ apart, so
  A − B holds the target twice: positive at A, negative at B. Both are extracted; each pair
  gives two spectra, which doubles the data at no cost.
* **Extended sources: object–sky.** A galaxy fills much of the slit, so the telescope moves
  off it to blank sky (≈ 50″ away for NGC 134). A − B holds the target once. A sky frame can
  serve the object frames on either side of it.

Other modes exist: **A** (one frame, no sky subtraction), **A − Sky** (each object frame minus
a sky made from all the sky frames), and **Combined** (all pairs averaged into one image first,
then extracted once).

## The spatial profile

Each order is resampled onto a grid of (wavelength, arcsec along the slit), using the
wavelength and slit-position images. At each wavelength the median along the slit is
subtracted; wavelengths where the atmosphere transmits less than 40 % are discarded (they are
noise); the **profile** is the median over wavelength, normalised.

```{themed} science/profile
:alt: The spatial profile of HD 225200, order 5

The spatial profile of HD 225200, order 5
```

## Apertures

**Point sources.** The strongest peaks of |profile| are found and fitted with Gaussians — two
for an A − B pair. Their signs (+ for A, − for B) come from the profile. Positions can instead
be typed as guesses (refined by the fit) or fixed.

**Extended sources.** The aperture is placed where you choose; by default on the profile's
peak (the nucleus), pair by pair. Several positions give spatially resolved spectra.

## The trace

A point source's spectrum is not exactly parallel to the order: atmospheric dispersion and
the optics move it by a fraction of a pixel along the order. Every 5 columns the profile of 5
summed columns is fitted with a Gaussian; the centres, as a function of wavelength, are fitted
with a polynomial (degree 1), points more than 3 pixels off rejected.
The aperture follows this trace.

## Measuring the light

For every wavelength of the resampled grid, the cell's value is the exact-area sum of the
detector pixels under it. Then, for a point source:

1. **Background.** The difference image still holds a residual: the sky changed a little
    between A and B. At each wavelength it is fitted across the slit, between 2.5″ and 7.5″
    from each aperture, with a polynomial of degree 1, and subtracted; its uncertainty is
    added to the error.

2. **Optimal extraction.** With the background subtracted, as the IDL tool does by default,
    the flux is not a plain sum but a profile-weighted estimate (Horne 1986). The spatial
    profile of the order, measured from the image itself, is fitted along the order, so that
    at each wavelength a profile $P(y)$ is known. Rebuilt at the wavelength and interpolated
    onto the slit, it is scaled to the data with a robust straight-line fit (points beyond
    6σ left out), forced to the aperture's sign, and normalized to unit sum within the PSF
    radius (2.2″). Every pixel in the aperture then gives an estimate of the whole flux,
    ${\rm data}(y)/P(y)$, with variance $\sigma^2(y)/P(y)^2$, and the flux is their
    variance-weighted mean:

    $$
    f(\lambda) = \frac{\sum_y {\rm data}(y)\,P(y)/\sigma^2(y)}{\sum_y P(y)^2/\sigma^2(y)},\qquad
    \sigma_f^2 = \frac{1}{\sum_y P(y)^2/\sigma^2(y)} .
    $$

    Because the profile is normalized within the PSF radius, this is the flux of the whole
    star inside 2.2″, not only of the 1″ aperture: on this night's standard it is 17–22%
    higher than the 1″ sum, with the same signal-to-noise.

3. **Without the background** (the switch off), the flux inside the aperture is summed, the
    pixels at the aperture's edges counted by the fraction inside it:

    $$
    f(\lambda) = \sum_{y} w(y)\,{\rm work}(\lambda,y),\qquad
    \sigma_f^2 = \sum_y w(y)^2\,\sigma^2_{\rm work} ,
    $$

    with $w(y)\in[0,1]$ the fraction of each slit pixel inside the aperture.

An extended source is always summed this way, over the apertures you place.

## The files

`proc/<source>_spec<NNNN>.fits`, one for each image of the target in each pair (A and B for
a point source), named after the frame it came from: for every order, wavelength (µm),
flux (DN s⁻¹), error and flags, in Spextool's format.

:::{admonition} One difference from the documentation of the IDL tool
:class: note
For extended sources the IDL tool *draws* the background regions but never subtracts
them. The default here gives the IDL tool's
results; *Subtract the background* applies what was intended. See
[Behaviour of the IDL tool](../validation/idl_differences.md).
:::
