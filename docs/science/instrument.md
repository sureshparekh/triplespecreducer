# Near-infrared spectroscopy with TripleSpec

## The instrument

**TripleSpec 4.1** is a cross-dispersed near-infrared spectrograph on the 4.1 m SOAR
telescope. In one exposure it records the whole range from **0.94 to 2.47 µm** — the z, J, H
and K bands — through a slit 1.1″ wide and 38″ long, at a resolving power
$R = \lambda/\Delta\lambda \approx 3500$ (about 85 km s⁻¹).

A grating disperses the light; its orders 3 to 7 overlap in angle, so a prism disperses them
again at right angles. On the 2048 × 1024 detector they appear as five curved strips, one
above the other, each the image of the slit spread out in wavelength:

```{themed} science/raw-orders
:alt: A raw frame with the five orders

A raw frame with the five orders
```

| order | covers (µm) | band | dispersion (Å pixel⁻¹) |
|---|---|---|---|
| 3 | 1.88–2.47 | K | 2.9 |
| 4 | 1.41–1.85 | H | 2.2 |
| 5 | 1.13–1.48 | J | 1.7 |
| 6 | 0.95–1.24 | z/J | 1.4 |
| 7 | 0.94–1.06 | z | 1.2 |

The plate scale along the slit is 0.417″ per pixel; the slit is 2.4 pixels wide, which sets
the spectral resolution.

## What stands between the detector and the spectrum

What a pixel records is

$$
{\rm DN}(x, y) \;=\; g(x,y)\Big[\,T_\oplus(\lambda)\,R(\lambda)\,F^{\rm target}_\lambda \;+\; S_{\rm sky}(\lambda, t) \;+\; B_{\rm thermal}\Big] \;+\;\text{noise},
$$

and every step of the reduction removes one of these terms:

* **$g(x,y)$**, the gain of each pixel, which differs by a few per cent from pixel to pixel
  → removed by the [flat field](flat.md).
* **$\lambda(x,y)$**: which wavelength falls on which pixel is a curved, tilted map
  → measured by the [wavelength calibration](wavecal.md).
* **$S_{\rm sky}$**, the sky's emission — mostly OH lines from the upper atmosphere, far
  brighter than faint targets and varying by the minute — and **$B_{\rm thermal}$**, the
  glow of the warm telescope → removed by nodding and subtracting, then
  [extraction](extraction.md).
* **$T_\oplus(\lambda)$**, the atmosphere's transmission, and **$R(\lambda)$**, the
  instrument's response → removed together by the [telluric correction](telluric.md),
  which also gives the flux its units.

## The atmosphere

```{themed} science/atmosphere
:alt: Atmospheric transmission and the TripleSpec orders

Atmospheric transmission and the TripleSpec orders
```

Water vapour and carbon dioxide absorb strongly between the bands: around 1.35–1.45 µm and
1.8–1.95 µm almost nothing gets through, which defines the J, H and K *windows*. Within the
windows hundreds of narrow absorption lines remain, changing with airmass and humidity —
which is why a standard star must be observed close in time and airmass to every target.

## Why nodding

The sky in the H band is about 14 mag arcsec⁻², a galaxy's outskirts 18 or fainter: the sky
dominates every pixel. It cannot be modelled well enough to subtract, because the OH lines
vary in strength by tens of per cent in minutes. Instead the telescope **nods**: it observes
the target at position A, then moves it to B, and the two frames are subtracted. The sky
(and the thermal background, and the detector's dark current) is the same in both and
cancels; the target does not. See [extraction](extraction.md) for the two nodding patterns.
