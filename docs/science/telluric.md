# Telluric correction and flux calibration

*The same step in Spextool: xtellcor, the method of Vacca, Cushing & Rayner
(2003, PASP 115, 389).*

## The idea

What was recorded of the target is its true spectrum times the atmosphere's transmission
$T_\oplus$ times the instrument's response $R$:

$$
O_{\rm target}(\lambda) = F_{\rm target}(\lambda)\; T_\oplus(\lambda)\, R(\lambda).
$$

An **A0 V star** observed minutes away, at a similar airmass, went through the same
$T_\oplus R$. Its true spectrum is known: A0 V stars are almost featureless apart from
hydrogen lines, and they are all like **Vega**, whose spectrum has been modelled in detail
and calibrated in absolute flux. So

$$
\frac{F_{\rm Vega\ model}(\lambda)}{O_{\rm star}(\lambda)} = \frac{1}{T_\oplus R}
\qquad\Rightarrow\qquad
F_{\rm target} = O_{\rm target}\times\frac{F_{\rm Vega\ model}}{O_{\rm star}} .
$$

The ratio is the **telluric spectrum**; it removes the atmosphere and the instrument, and —
because Vega's model is in physical units, scaled to the star's magnitude — gives the target
its flux.

## The difficulty: the star's own lines

The A0 V star has deep, broad hydrogen lines (Paschen and Brackett series). If Vega's model
lines are not *exactly* like the star's as TripleSpec saw them — same depth, same width, same
velocity — the mismatch is copied into every corrected target as false emission or absorption
at those wavelengths. Vega's lines must be made to match:

1. **Velocity.** Vega's model is cross-correlated with the star around a hydrogen line, on a
    logarithmic wavelength grid (where a velocity is a constant shift); a fit to the
    correlation peak gives the shift, here 34 km s⁻¹.

2. **The kernel.** The star's line is the model's line blurred by the instrument (and by any
    difference in rotation). That blurring, the **kernel** $K$, is found by deconvolution: in
    Fourier space a convolution is a product, so

    $$
    \hat K = \frac{\widehat{O_{\rm line}}}{\widehat{V_{\rm line}}}\;\Phi,
    $$
    
    where $O_{\rm line}$ and $V_{\rm line}$ are the star's and Vega's normalised line
    profiles (1 − the continuum-normalised flux) around **Paδ** (1.005 µm) in order 6, and
    $\Phi$ is a smooth window that cuts the frequencies where noise dominates (at 10 times the
    narrower line's width). The equivalent widths are matched with a scale factor.

3. **Its quality.** Vega's line convolved with $K$ is compared with the star's: the maximum
    and RMS deviation, 1–2 % being typical.

4. **Every order.** The kernel, made on Vega's velocity grid, is resampled to each order's
    pixel scale (its dispersion).

```{themed} science/kernel
:alt: The kernel built on Paδ

The kernel built on Paδ
```

## The telluric spectrum

For every order, Vega's model is rebuilt as the star would look through
TripleSpec:

1. Vega's model comes as three parts: the full spectrum $F_V$, its continuum $C_V$, and a
    smoother continuum $C'_V$. The **line spectrum** $F_V/C_V - 1$ and the **residual
    continuum** $C_V/C'_V - 1$ are each multiplied by the **line scales** and convolved with
    the kernel; the model is then $(1 + \ldots)(1 + \ldots)\,C'_V$.

2. The line scales: at each H I line a **control point** holds a scale (1 to start), and a
    spline under tension joins them, so each line's depth can be adjusted.

3. The model is shifted to the star's velocity and resampled to its wavelengths.
4. It is **reddened** by the star's colour excess $E(B-V) = (B-V)_\star - (B-V)_{\rm Vega}$
    (Vega's being 0; never negative) with the Cardelli–Clayton–Mathis law, $R_V = 3.1$,
    keeping the V-band flux fixed, and **scaled** to the star's V magnitude,
    $10^{-0.4\,(V_\star - 0.03)}$ (Vega: V = 0.03).

5. It is divided by the observed star; the error follows from the star's.

The control points are where the user's judgement enters: if the telluric spectrum shows a
bump at a hydrogen line, Vega's line is too deep there and its scale is lowered; a dip, and
it is raised. An automatic estimate fits a Lorentzian to the line.

```{themed} science/telluric
:alt: The standard, the telluric spectrum, and the target before and after

The standard, the telluric spectrum, and the target before and after
```

## The shift

The target and the standard were observed with the instrument pointing differently, and it
flexes: the absorption features may be offset by a fraction of a pixel between them. For each
order the telluric spectrum is shifted by −1.5 to +1.5 pixels in steps of 0.02; the shift
that leaves the least structure (lowest RMS) in the corrected target over a region of strong
absorption is kept, refined with a parabola.

## The corrected target

The target is multiplied by the shifted telluric spectrum, order by order. Its error includes
the standard's (unless that is switched off):

$$
\sigma_F^2 = T^2\sigma_O^2 + O^2\sigma_T^2 .
$$

The result is in erg s⁻¹ cm⁻² Å⁻¹, or in W m⁻² µm⁻¹, erg s⁻¹ cm⁻² Hz⁻¹, W m⁻² Hz⁻¹ or Jy.

## The IP method

Instead of deriving the kernel from the star, it can be built from a stored instrument
profile. Spextool holds these only for slits of 0.3, 0.5, 0.8, 1.6 and 3.0″ (and TS4_SOAR's
file is a copy of another instrument's); for TripleSpec's 1.1″ slit the IDL tool stops, and so
does this program, pointing to the deconvolution method.
