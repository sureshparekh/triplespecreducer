# Validation against IDL

Every step was run twice on the same frames of the night of 2025-11-02 — once by the IDL tool
(IDL 7.1.1, Spextool 4.1 with the TS4_SOAR module, its steps called without its widgets by the
scripts in `validation/idl/`), once by this program — and the outputs compared value by value.
Each stage was fed IDL's own output from the stage before, so a difference belongs to the stage
that shows it.

**The result: every product is bit-identical to IDL's** — the same 32- or 64-bit number in every
pixel — with one exception, the telluric kernel, explained below.

## The steps

| step | compared | agreement |
|---|---|---|
| flat | combined (scaled) flat, 2 097 152 pixels | **identical** |
| | order edges and x ranges | **identical** |
| | normalization model, normalized flat, its variance | **identical** |
| | per-order RMS | **identical** |
| wavelength | 1D sky spectra, cross-correlation offset | **identical** |
| | 318 OH line centroids (1D) and the 1D solution, RMS | **identical** |
| | 7 046 line points traced along the slit | **identical** |
| | wavelength image, slit-position image | **identical** |
| | CuHeAr arc lamp: offset, lines per order, 1D solution, traced lines, both images | **identical** |
| extraction | A−B work image | **identical** |
| | spatial profiles, aperture positions, trace coefficients | **identical** |
| | point-source spectra (2 apertures × 5 orders, PS background) | **identical** |
| | extended-source spectra (5 orders) | **identical** |
| | the six spectrum files written by the A−B extraction | **identical** |
| combine | scale factors, combined spectra (point and extended source) | **identical** |
| telluric | Vega velocity, EW scale, maximum deviation, line centre | **identical** |
| | kernel resampled to each order | 10⁻⁶ (IDL's FFT, see below); from IDL's kernel **identical** |
| | telluric spectrum (given IDL's kernel) | **identical** |
| | shifts found on the object (autofind) | **identical** |
| | corrected object spectrum | **identical** |
| merge | merged spectrum; with automatic scale, trims and a gap | **identical** |
| clean | Gaussian, Savitzky–Golay | **identical** |
| options | combined image (A−B, ABBA→ABAB) and its A-mode extraction | **identical** |
| | combined image with the sky subtracted (`mc_subspecsky`) | **identical** |
| | `xcleanspec` fix and remove of a region | **identical** |
| | units (Jy, W m⁻² µm⁻¹) | **identical** |
| | edited H I control points, tension spline, telluric from them | **identical** |
| | 'estimate scale' around Br 11: Lorentzian fit and scale | **identical** |

| files | flat, wavelength-calibration, spectrum and combined-spectrum headers, card by card | **identical** |

## The chain

The table compares each stage on IDL's input. `validation/compare_chain.py` runs the stages
as a chain instead — flat, wavelength calibration, extraction of four pairs, combination —
each Python stage reading only the files the Python stage before it wrote, and compares every
file with IDL's chain: **all data and all header cards are identical**. This matters because
Spextool passes numbers between steps as header text (see [Behaviour of the IDL tool](idl_differences.md)): the order
edges the wavelength calibration uses are the flat file's 12-digit values, not the fit's full
precision, and the dispersions the telluric step uses are the wavelength file's 6-digit ones.
Writing and reading headers exactly as IDL does moves a whole-night reduction by about 10⁻⁷
(median) from one that keeps full precision — and puts it on IDL's numbers.

`tests/test_parity.py` asserts these identities (and the kernel's tolerance), so any change
that moves a single bit fails.

## The one exception: IDL's FFT

The Vega kernel (`mc_vegaconv`) divides two Fourier transforms. IDL computes them with its own
single-precision mixed-radix FFT, whose rounding is not that of any available FFT (numpy's in
double, scipy's in single: both agree with IDL to 10⁻⁷ but match exactly only a few values in
a hundred). The kernel therefore agrees to its last bit or two, about 10⁻⁶ of its peak, and
everything computed from it — the telluric spectrum, then the corrected object — inherits
differences of that size (a few parts in 10⁶). Fed IDL's own kernel, every later step is
identical, as the table shows. The optional shape correction of `xcombspec` (also FFT-based)
is in the same situation.

## How exact agreement was reached

Getting the right algorithm is not enough for identical numbers: IDL's arithmetic has to be
reproduced too. Each rule below was established by running probes through IDL 7.1 — thousands
of values each — and is applied wherever the IDL applies it.

* **Types.** Spextool keeps most arrays in single precision (`fltarr`, `replicate(!values.f_nan)`,
  structure fields); assigning into them rounds. Each Python quantity carries the IDL type of
  its counterpart, including where IDL mixes them (a `FLOAT` literal such as `2.354` or `0.1`
  enters a double expression as its single-precision value).
* **Running sums.** `TOTAL` adds in order; above 100 000 elements IDL splits the sum over its
  thread pool (64 chunks on the reference machine) and adds the partial sums in double. Matrix
  products (`#`) are running sums in the result's type. Neither is what numpy does.
* **Library functions.** `EXP`, `ALOG`, `SIN`, `COS` of single-precision data are the C
  library's `expf` and friends; `^` with a non-integer exponent is `pow`; integer powers are
  repeated squaring; `x^2` is `x*x`. numpy's own single-precision routines differ in up to a
  third of the values.
* **INTERPOLATE** rounds its subscripts to single precision even when they are double, and has
  its own forms for linear and cubic interpolation; **INVERT** is a Gauss–Jordan elimination
  that does not stop at a singular matrix; **CURVEFIT**, **MPFIT**, **POLY_FIT**, **SPLINE**,
  **LUDC/LUSOL**, **CONVOL**, **C_CORRELATE**, **MOMENT** follow the same published algorithms,
  computed in those types.
* **The exact-area pixel clipping** of the extraction was checked against the IDL tool's own
  compiled clipping library on 21 819 pixel areas.

## The building blocks

Thirty-odd IDL routines are checked on random inputs, edge cases included (NaNs, even-length
medians, kernel edges, points off an image) by `validation/primitives_parity.py`: all agree
exactly, except bilinear (not cubic) `INTERPOLATE` of a double image, which no step uses, and
the values `TRIGRID` interpolates, which agree to 10⁻¹⁵: they are only ever used as pixel
coordinates, which `INTERPOLATE` rounds to single precision, and after that rounding all
218 673 coordinates of a test order are identical.

## Running it yourself

```bash
pytest tests/test_parity.py -v                 # every step against the stored IDL outputs
python validation/primitives_parity.py         # the building blocks (needs IDL)
idl -e "@validation/idl/ref_flat"              # regenerate an IDL reference (needs IDL)
python validation/compare_flat.py              # ... and compare
```

`compare_wavecal.py`, `compare_extract.py`, `compare_dg.py` (combine, telluric, merge, clean)
and `compare_extras.py` (the options) do the same for the other steps.

```{toctree}
:hidden:

idl_differences
```
