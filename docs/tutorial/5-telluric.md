# 5. Telluric correction

:::{admonition} Why
:class: note
What reaches the detector is the galaxy's light **times** the transmission of Earth's
atmosphere **times** the response of the telescope and instrument. The A0 V star was
observed through nearly the same atmosphere minutes earlier, and its true spectrum is
known — it is like Vega's, which is modelled in detail. So *model ÷ observed star* is
exactly the factor that turns the galaxy's counts into flux and removes the atmosphere.
The difficulty is that the star has its own deep hydrogen lines, which must not end up in
the correction: Vega's model is first made to match the star's lines exactly.
:::

## 1. Load the two spectra

Go to **5 Telluric**. The target **NGC134** and its standard **HD225200_1** are chosen from
the plan, with B = 6.40 and V = 6.39. The caption reports the airmasses — 1.088 and 1.079,
a difference of 0.008 — and the spectral type SIMBAD gives. Press **Load the two spectra**.

```{shots} 05-telluric-loaded
:alt: The two spectra loaded

The two spectra loaded
```

The stage shows the standard's order 6 around **Paschen δ** (1.005 µm), with the continuum
windows (green) either side of the line and the line region (blue).

## 2. Build the kernel

Leave the method at **Deconvolution**, the order at 6, the continuum windows at 0.990–0.997
and 1.013–1.020 µm with degree 3, and the line at 0.997–1.013 µm. Press **Build the kernel**.

```{shots} 05-telluric-kernel
:alt: The kernel fit

The kernel fit
```

What happened:

1. The star was divided by a polynomial through the continuum windows.
2. Vega's model was cross-correlated with the star: Vega must be shifted by
    **34.0 km s⁻¹** to line up.
3. The ratio of the Fourier transforms of the two lines gives the **kernel**: the blurring
    that turns Vega's sharp line into the star's, as TripleSpec sees it (the instrument's
    resolution, and any difference in the stars' rotation).
4. Vega's line, convolved with that kernel, is drawn over the star's.

:::{admonition} Check
:class: tip
The convolved Vega line should lie on the star's. The caption gives the **maximum and
RMS deviation**: here **1.96 %** and **0.62 %** — good. Several per cent means the
continuum windows or the degree need changing (try degree 3 if you used 1, or narrower
windows).
:::

## 3. The hydrogen lines

The correction for every order is now computed. The **Hydrogen lines** panel lists, for the
selected order, a **control point** at each H I line of the star, with the scale applied to
Vega's line depth there.

```{shots} 05-telluric-lines
:alt: The telluric spectrum of order 4 with its control points

The telluric spectrum of order 4 with its control points
```

The stage shows the **telluric spectrum** of the order (top) — Vega's model divided by the
observed star — with the control points marked, and the corrected target beneath it.

:::{admonition} Check
:class: tip
The telluric spectrum should be **smooth through every H I line**. A bump at a line means
Vega's line is too deep there: lower its scale; a dip means raise it. Type a new scale in
the table, or tick *Estimate the scale by dragging* and drag across the line. Here the
defaults leave it smooth, and nothing needs changing.
:::

## 4. The shifts

Small flexures between the standard's and the target's observations shift the absorption
features by a fraction of a pixel. Press **Find every shift**:

```{shots} 05-telluric-shifts
:alt: The shifts found per order

The shifts found per order
```

For each order the telluric spectrum was shifted from −1.5 to +1.5 pixels and the shift that
leaves the least structure in the corrected target was kept: −0.18, −0.27, −0.25, −0.25 and
−0.33 pixels for orders 3 to 7. A region can be typed for an order where the automatic one
does not suit.

## 5. Write the corrected spectrum

Press **Write the corrected spectrum**.

```{shots} 05-telluric-corrected
:alt: NGC 134 corrected

NGC 134 corrected
```

The galaxy is now in **erg s⁻¹ cm⁻² Å⁻¹**, with the atmosphere removed: the J, H and K
windows join smoothly, and only the opaque gaps (around 1.35–1.45 and 1.8–1.95 µm) remain
noisy.

:::{admonition} Check
:class: tip
No residual absorption bands outside the opaque gaps, and no bumps or dips at the
hydrogen lines (Paβ 1.282, Brγ 2.166 µm and the Brackett series in H).
:::

**Next:** [6. Merging the orders](6-merge.md) · *The science:* [telluric correction](../science/telluric.md)
