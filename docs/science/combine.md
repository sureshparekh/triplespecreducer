# Combining spectra

*The same step in Spextool: xcombspec, described in Cushing, Vacca & Rayner (2004).*

## Why

N spectra of the same target, combined, give a signal-to-noise ratio √N times higher, and
anything that is in only one of them — a cosmic ray, a bad pixel, a satellite trail — can be
rejected. Their flux levels differ by a few per cent, because the seeing and the guiding
change during the sequence and the slit catches a different fraction of the light; so they
are put on a common level first.

```{themed} science/combine
:alt: Four spectra of HD 225200 and their combination

Four spectra of HD 225200 and their combination
```

## The algorithm

1. **A common grid.** Every spectrum is interpolated onto the wavelengths of the first.
2. **Scale factors.** In one order (5 by default), over its middle 80 %, each spectrum is
    compared with the median spectrum; the median of the ratio is its scale factor $s_i$,
    applied to all its orders.

3. **Combination**, at every wavelength, with one of nine statistics. The default, the
    **robust weighted mean**:

    $$
    \bar f = \frac{\sum_i f_i/\sigma_i^2}{\sum_i 1/\sigma_i^2},\qquad
    \sigma_{\bar f}^2 = \frac{1}{\sum_i 1/\sigma_i^2},
    $$
    
    over the values within 8 robust standard deviations (1.4826 × the median absolute
    deviation) of the median. The others are the robust mean (with its RMS or its standard
    error), the plain mean and weighted mean, the median (with the MAD or the median error),
    and the sum.

4. **Shape correction** (optional). When slit losses depend on wavelength, the continua of
    the spectra differ in shape; each is then multiplied by the smooth ratio of the mean
    spectrum to it before combining.

## Pruning

A spectrum taken through cloud, or with the target drifting off the slit, shows a scale
factor well below the others. Untick it and combine again; the log keeps the scale factors of
every combination.
