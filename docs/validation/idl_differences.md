# Behaviour of the IDL tool

TripleSpec Reducer follows the steps and the science of Spextool and its TripleSpec module, and
its results are compared with the IDL tool's on the same frames (see
[Agreement with IDL](index.md)). That comparison brought out places where the IDL tool does
something other than what its interface or its documentation suggests. The rule followed here:
**the IDL tool's behaviour is the default**, so results match it; where a behaviour is plainly
unintended, the intended one is available **as an option, off by default**.

## Matched, with an option for the intended behaviour

**Extended sources: the background is never subtracted.** In the IDL tool's two-dimensional
mode, the background regions chosen for an extended source are drawn on the profile, and are
recorded in the output header, but no background is fitted or subtracted from the spectrum.
For a galaxy observed object–sky, the only sky removal is therefore the A − sky subtraction.
*Option*: `xs_background` (Extract → "Subtract the background").

**Combine images with "scale orders" stops** with an error in the IDL tool when every scale
factor is finite. *Option*: `scale_orders="intended"` in `combineimgs.combine_images`.

## Matched as they are

* **The one-dimensional sky spectra for the wavelength solution come from the first frame
  only**, not from the median of the frames; the tracing of the lines along the slit does use
  the median image.
* **In that one-dimensional extraction the aperture's edge pixels count for nothing**: the
  fractional weights of the pixels at the aperture's edges are truncated to zero.
* **The kernel step uses Vega's first continuum model** where the second would be expected,
  both for the velocity and for the normalised Vega spectrum.
* **The aperture signs of a point source are taken from the last order.**
* **Extended sources in A−B mode write only the A file**: every aperture has a positive sign.
  This is the "second file does not exist" message seen in the tutorial.
* **Applying an existing telluric spectrum to another object** computes the error after
  correcting the flux, so its second term uses the corrected flux; an object order with no
  telluric counterpart ends the loop, leaving the later orders uncorrected.
* **Merging without flags** clears one spectrum's flags twice and the other's never (the
  merge panel always passes flags, so it does not arise there).
* **Combine images, "subtract sky"** abandons the rest of an order when one column falls off
  the array.
* **Combine images, "Weighted Mean"** is a plain mean: no variance reaches the combination.
* **The flat is a plain mean** of the scaled lamp frames, although its history says "Robust
  Mean": no rejection threshold reaches the combination.
* **The Savitzky–Golay default of the cleaning panel** (1.5 slit widths = 3.6 pixels)
  truncates to 3 pixels, which with degree 2 leaves the spectrum unchanged.
* **The amplifier-correction switch** of the TripleSpec reader is read and never used.
* **Bridging a cut region of the telluric spectrum** replaces the flux by a straight line, but
  the error is bridged only where the original array is already undefined.
* **Each order's kernel is chosen by the order's position among the object's orders**, while
  the kernels are made in the standard's order. The two agree whenever the standard and the
  object have the same orders, as they always do for TripleSpec.
* **"Estimate scale" in the hydrogen-line panel** starts the Lorentzian fit from the
  continuum's slope where it would be expected to start from the line's depth.
* **Merging fails** ("subscript out of range") when the order being added extends past both
  ends of the merged spectrum — reachable only by trimming the merged spectrum. The same
  happens here (an error).
* **The combined spectrum's header is the first spectrum's**: the "average" UT and airmass
  written are the first spectrum's values.
* **The telluric step takes the slit width from the object**, not the standard, and keeps the
  slit width, the airmasses, the magnitudes, the velocity shift and the deviations in single
  precision.
* **The IP kernel method cannot run for TripleSpec's 1.1″ slit.** Instrument-profile
  coefficients exist only for 0.3, 0.5, 0.8, 1.6 and 3.0″ slits, and the IDL tool stops for any
  other width (TS4_SOAR's coefficient file is in any case another instrument's). Here the IP
  method raises the same condition as a clear error and points to the deconvolution method.
* **The modified-Vega file's `GRAT` keyword is 0**: it is read from the standard's header,
  which has no such keyword.

## Header values: text, then read back

Spextool passes numbers from one step to the next through FITS headers. It writes each number as
text — double precision with 12 significant digits, single precision with 6, the slit sizes in
the spectrum files with 5 — and the next step reads that text back, with a type that depends on
the text: a real is double precision when its text, up to the `/`, is eight characters or
longer, single precision when shorter, and an integer when it has no decimal point. So the order
edges a wavelength calibration uses are the flat file's 12-digit values, and a slit height
written as `38` is read back as an integer. Every header written here holds the same text, and
every reader applies the same typing, so a reduction run entirely in this program carries
exactly the numbers the IDL chain carries (`validation/compare_chain.py`). One cosmetic
exception: the `#### Xspextool History ####` banner line of the IDL files is not a valid FITS
record; it is written here as a standard blank-keyword card with the same text.

## Not reproducible, so done right

**Combine images without scaling or sky subtraction writes a scrambled image** in the IDL tool:
the pixel values are right, but stored in the wrong orientation, so the file cannot then be
extracted. Here the same values are written in the right orientation.

**Merging orders unweighted stops in the IDL tool** with an undefined-variable error. Here the
unweighted merge is done as the panel intends: every error is set to 1, so the overlap is the
plain mean.

## Choices made differently from a panel's starting value

**The kernel's continuum degree is 3, not the panel's 1** — the tutorial's choice. On the test
night degree 1 leaves a 3.5 % maximum and 1.4 % rms deviation in the kernel fit, degree 3 2.0 %
and 0.6 %. `Options.idl_defaults()` restores 1.

## Where agreement stops being exact

Everywhere except one place, the numbers are IDL's to the last bit (see [Agreement with IDL](index.md)).
The exception is the telluric kernel: IDL makes it with its own single-precision FFT, which no
available FFT reproduces bit for bit, so the kernel — and the telluric correction made from
it — agree to a few parts in 10⁶.

## The reference runs

The IDL reference outputs were made by scripts that call the IDL tool's steps without its
widgets (`validation/idl/`). Two reviews of those scripts against the tools' own behaviour found
places where the scripts had departed from what the tools do — the keywords copied from the raw
headers, the header of the combined spectrum, the type of the linearity maximum, the axis titles,
and the type of two stored values — and they were corrected and every reference regenerated;
the data did not change.
