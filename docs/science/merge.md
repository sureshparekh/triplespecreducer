# Merging and cleaning

## Merging the orders

*The same step in Spextool: xmergeorders (Cushing, Vacca & Rayner (2004).*

The five orders each cover part of the range and overlap their neighbours at the ends. One
order is the **anchor** (3 by default); the others are added to it one at a time. Each can be

* **scaled** — by 1 (the tutorial's choice: after the telluric correction the orders are on
  one flux scale already), by a typed factor, or *automatically*: the factor that best matches
  the merged spectrum over a chosen stretch of the overlap, weighted by the errors;
* **trimmed** — its end, or the merged spectrum's, cut at a wavelength, to drop a noisy edge.

Where two orders overlap, the result is the error-weighted mean,

$$
f = \frac{f_1/\sigma_1^2 + f_2/\sigma_2^2}{1/\sigma_1^2 + 1/\sigma_2^2},\qquad
\sigma^2 = \frac{1}{1/\sigma_1^2 + 1/\sigma_2^2},
$$

or the plain mean (unweighted), or the sum. Flags are combined.

```{themed} science/merge
:alt: The five orders and the merged spectrum

The five orders and the merged spectrum
```

The merged spectrum is the final product: `proc/<source>_merged.fits`, one array of
wavelength, flux, error and flag.

## Cleaning

*The same step in Spextool: xcleanspec, with its 'f' (fix) and 'r' (remove) keys; the filter is
Savitzky & Golay's (1964).*

Optional. **Smoothing** lowers the noise at the cost of resolution:

* **Savitzky–Golay**: at each pixel, a polynomial of the given degree is fitted to the pixels
  within the window and evaluated there. It preserves the heights and widths of features
  better than a moving average.
* **Gaussian**: convolution with a Gaussian of the given full width at half maximum.

The errors are propagated through either (the smoothing weights, squared, applied to the
variances).

**Repairs**: *fix* replaces a region by a straight line between its ends (and flags it); *remove*
sets it to NaN — for the opaque gaps between the bands, or a residual from a bad pixel.
