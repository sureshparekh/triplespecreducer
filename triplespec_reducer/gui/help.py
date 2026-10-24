"""
Every explanation the interface shows, in one place.

The desktop interface draws its captions, tooltips and step notes from here,
and the manual's settings reference is generated from the same text
(docs/gen_reference.py), so the two cannot drift apart.

STEPS[i]  : one entry per step: `why` (what the step is for), `how` (what
            it does, in order), `check` (what to look at before moving on),
            `idl` (the Spextool tool that does the same step, and where it is described).
FIELDS[k] : one entry per control: `label` (beside the control), `caption`
            (one line under it), `details` (the full explanation: tooltip
            and manual), `default` (the IDL tool's starting value).
"""

STEPS = [
    dict(
        title="Night",
        why=(
            "A TripleSpec night is a few hundred frames: lamp flats with the lamp on and off, "
            "arc lamps, science targets nodded along the slit, and A0 V standard stars taken "
            "close in time and airmass to each target. Everything after this step needs to "
            "know which frames belong together: which are flats, which form the A-B pairs of "
            "each target, which standard corrects which target, and where each target's "
            "wavelength solution comes from."),
        how=[
            "Every header in the raw folder is read (nothing is copied, written or moved there).",
            "Frames are sorted by type: lamps-on and lamps-off flats, arcs, and objects. Focus "
            "frames and tests are set aside.",
            "Objects are grouped into targets by name and pointing; a pointing more than a few "
            "arcseconds away becomes a separate target (for example NGC134_pos2).",
            "Each target's nod pattern (ABBA, or object-sky for extended sources) is read from "
            "the telescope offsets and its frames are paired A with B.",
            "Targets are classed as science or A0 V standard (the standards are looked up on "
            "SIMBAD for B, V and spectral type), and each science target is given the standard "
            "closest to it in time and airmass.",
            "The plan is written to night.json in the output folder. Every choice can be "
            "changed in the table and saved again."],
        check=[
            "Each target has the role, source type and pattern you expect.",
            "Each science target has a standard, with an airmass difference below about 0.1.",
            "The notes under the table: they report anything unusual (a standard that is not "
            "A0 V, a lone frame, a pointing kept apart).",
            "B and V of each standard: they set the flux scale of the final spectrum."],
        idl="xspextool's file lists, typed by hand in the IDL tool; here read from the headers.",
    ),
    dict(
        title="Calibrate",
        why=(
            "Two calibrations are needed before any spectrum can be extracted. The flat field "
            "corrects the pixel-to-pixel sensitivity of the detector and tells the software "
            "where each of the five cross-dispersed orders falls on it. The wavelength "
            "calibration maps every pixel to a wavelength and to a position along the slit, "
            "because TripleSpec's orders are curved and its slit image is tilted."),
        how=[
            "Flat: the lamps-off frames are averaged and subtracted from the lamps-on frames, "
            "which removes the thermal background.",
            "The order edges are found by following the edges of the flat's illuminated bands "
            "up the detector and fitting a polynomial to each (degree 4).",
            "The flats are scaled to a common level and averaged; each order is then divided "
            "by a smooth two-dimensional model of itself, so the flat keeps only the "
            "pixel-to-pixel variations, normalised to 1.",
            "Wavelengths: a one-dimensional spectrum is extracted down the middle of the slit "
            "in every order, from the target's own sky lines (OH emission, median of its "
            "frames) or from the CuHeAr arc lamp.",
            "That spectrum is cross-correlated with a stored template to find the pixel offset, "
            "the known lines are located, and one polynomial in column and order number is "
            "fitted to all of them (lambda x order / home order = P(column, order)).",
            "Each line is then traced along the slit, and a two-dimensional fit gives the "
            "wavelength and the slit position of every pixel. These two images are the "
            "wavelength calibration file."],
        check=[
            "Show the flat: the five orders should be outlined by their edges, with the "
            "normalised flat close to 1 inside them.",
            "The table: every solution should have tens of lines per order and an RMS of a "
            "fraction of an angstrom.",
            "Select a solution: the traced lines should follow the tilt of the slit smoothly."],
        idl="xspextool's calibration step (flat field and 2D wavelength calibration; "
            "Cushing, Vacca & Rayner 2004).",
    ),
    dict(
        title="Extract",
        why=(
            "The infrared sky is bright and variable; it is removed by nodding the target "
            "between two positions on the slit (A and B) and subtracting one frame from the "
            "other. In A-B the target appears twice, once positive and once negative, and the "
            "sky cancels. Extraction turns each pair into a spectrum per order by adding up the "
            "light inside an aperture around each image of the target."),
        how=[
            "The pair is read, converted to counts per second, differenced (A-B) and divided by "
            "the flat. The variance of every pixel is carried along.",
            "A spatial profile is made for every order: the light at each slit position, "
            "averaged along wavelength (with strong telluric absorption excluded).",
            "The apertures are found on the profile: the positive and the negative image of a "
            "point source, or the positions you give for an extended source.",
            "Point sources are traced: the aperture centre is followed along the order and "
            "fitted with a polynomial, so the aperture follows any tilt of the spectrum.",
            "The background left in the difference is fitted across the slit between two "
            "radii and subtracted, and the flux inside the aperture radius is summed, with "
            "fractional pixels at the edges.",
            "Each pair writes two spectra (A and B) for a point source, one for an extended "
            "source; every file holds wavelength, flux, error and flags for each order."],
        check=[
            "Preview the first pair: the profile should show the positive and negative peaks, "
            "with the apertures (shaded) on them and the background regions clear of them.",
            "The extracted spectra should look alike from pair to pair, apart from noise.",
            "A spectrum that stands out (a cloud, a guiding error) can be left out in the "
            "next step."],
        idl="xspextool's extraction ('Do All Steps': profiles, apertures, trace, extraction; "
            "Cushing, Vacca & Rayner 2004).",
    ),
    dict(
        title="Combine",
        why=(
            "A target observed in many pairs gives many spectra; they are combined into one "
            "with a higher signal-to-noise ratio. Their flux levels differ slightly because "
            "seeing, guiding and transparency change during the sequence, so they are first "
            "scaled to a common level."),
        how=[
            "Every spectrum is interpolated onto the wavelength grid of the first.",
            "Scale factors are computed from one order (the median over its middle 80 %) so "
            "that every spectrum matches the median level of all of them.",
            "The scaled spectra are combined pixel by pixel. The default, a robust weighted "
            "mean, rejects values more than the threshold (in sigma) away and weights the rest "
            "by their inverse variance; the error is that of the weighted mean.",
            "Optionally, the spectral shapes are matched first (shape correction), for spectra "
            "whose continua differ because of slit losses."],
        check=[
            "The scale factors should be close to 1 (a few per cent); a much lower factor marks "
            "a spectrum taken through cloud or with the target off the slit.",
            "Untick a spectrum that looks wrong and combine again."],
        idl="xcombspec (Cushing, Vacca & Rayner 2004).",
    ),
    dict(
        title="Telluric",
        why=(
            "Earth's atmosphere absorbs much of the near infrared (water and carbon dioxide), "
            "and the instrument's response varies strongly with wavelength. Both are measured "
            "with an A0 V star observed near the target: its true spectrum is known (it looks "
            "like Vega, whose spectrum is modelled), so dividing the model by the observed star "
            "gives the correction that turns the target's counts into flux and removes the "
            "atmosphere. This is Vacca, Cushing and Rayner's (2003) method."),
        how=[
            "Vega's model is shifted to the star's velocity, found by cross-correlating the "
            "star with Vega around a hydrogen line (Pa delta in order 6).",
            "A convolution kernel is derived that turns Vega's line into the star's line: the "
            "deconvolution method divides the Fourier transforms of the two line profiles. It "
            "carries the instrumental resolution and any difference in line width.",
            "The kernel is resampled to every order. Vega's model is convolved with it, scaled "
            "to the star's B and V magnitudes and reddened by their difference.",
            "The telluric spectrum is the convolved model divided by the observed star. The "
            "depth of every hydrogen line can be adjusted (control points joined by a tension "
            "spline) so the star's own lines do not imprint on the correction.",
            "The correction is shifted by a fraction of a pixel to best match the target's "
            "absorption features, then multiplied into the target's spectrum."],
        check=[
            "The kernel fit: the convolved Vega line over the star's should leave residuals of "
            "a few per cent at most (maximum and RMS deviation are reported).",
            "The telluric spectrum should be smooth through the hydrogen lines; adjust the "
            "scales where it shows a peak or a dip.",
            "The corrected target should not show residual absorption bands outside the "
            "opaque regions (around 1.35-1.45 and 1.8-1.95 um)."],
        idl="xtellcor (Vacca, Cushing & Rayner 2003; Cushing, Vacca & Rayner 2004).",
    ),
    dict(
        title="Merge",
        why=(
            "TripleSpec records five orders (3 to 7), each covering part of 0.94-2.47 um, "
            "with overlaps between neighbours. Merging joins them into one continuous spectrum."),
        how=[
            "One order is the anchor; the others are added to it one at a time.",
            "Each added order can be scaled (by hand, or automatically to match the merged "
            "spectrum over a region of the overlap) and trimmed at either end.",
            "In the overlap the two orders are averaged, weighted by their errors."],
        check=[
            "The orders should meet smoothly in the overlaps. A step there means a scale is "
            "needed (the tutorial keeps 1, since the telluric correction already put the "
            "orders on one flux scale)."],
        idl="xmergeorders (Cushing, Vacca & Rayner 2004).",
    ),
    dict(
        title="Clean & save",
        why=(
            "Optional last touches: smoothing, and repairing or removing bad regions, such as "
            "the opaque gaps between the atmospheric windows or a residual from a bad pixel."),
        how=[
            "Smoothing with a Savitzky-Golay filter (keeps line shapes) or a Gaussian of a "
            "given width in pixels; the errors are propagated.",
            "Fix: a region is replaced by a straight line between its ends. Remove: a region "
            "is set to NaN.",
            "The result is written next to the merged spectrum, and can be saved as text."],
        check=["Smoothing lowers the resolution: keep it narrower than the lines you measure."],
        idl="xcleanspec (Cushing, Vacca & Rayner 2004); Savitzky & Golay (1964).",
    ),
]


FIELDS = {
    # ------------------------------------------------------------ 1. night
    "night.raw": dict(
        label="Raw frames",
        caption="The folder of the night's raw FITS files. It is only read.",
        details="Every SPEC_*.fits file in this folder is read for its header. Nothing in it is "
                "written, moved or renamed, so there is no need to make a copy first."),
    "night.out": dict(
        label="Output",
        caption="Where everything is written: cal/ for calibrations, proc/ for spectra.",
        details="The plan (night.json), the calibration files (cal/), every spectrum (proc/) and "
                "the log of every step (qa/) go here. Each file is named after its target "
                "(sourcename_step.fits), so everything made for one object sorts together."),
    "night.simbad": dict(
        label="Look up the standards' B, V and type on SIMBAD",
        caption="Needs the internet; without it, type B and V in the table.",
        details="The B and V magnitudes of the A0 V standard set the absolute flux scale of the "
                "corrected spectrum, and the spectral type is checked (the method assumes an "
                "A0 V star)."),
    # ------------------------------------------------------------ 2. calibrate
    "cal.flatfield": dict(
        label="Divide the images by the flat",
        caption="On in the IDL tool; turn off only to inspect the raw response.",
        default="on (ts4_soar.dat FLATFIELD=1)",
        details="Every frame is divided by the normalised flat before extraction, which removes "
                "the pixel-to-pixel sensitivity differences of the detector."),
    "cal.autox": dict(
        label="Find the offset automatically (Auto X-Correlate)",
        caption="The spectrum is cross-correlated with the stored template.",
        default="on",
        details="The one-dimensional sky (or arc) spectrum of order 4 is cross-correlated with "
                "the template the line list was made on; the peak gives the pixel offset that "
                "places every known line. Turn it off to type the offset yourself when the "
                "correlation fails (very few lines, a shifted grating)."),
    "cal.offset": dict(
        label="Offset (px)",
        caption="Used only when the automatic offset is off.",
        details="The shift, in pixels, between this night's spectrum and the template, in "
                "order 4. Find it by comparing one line in both."),
    # ------------------------------------------------------------ 3. extract
    "ex.source": dict(
        label="Source",
        caption="The target whose frames are extracted.",
        details="Standards and science targets are both extracted the same way."),
    "ex.mode": dict(
        label="Mode",
        caption="How the frames are turned into images before extraction.",
        default="A-B pairs",
        details="A-B pairs: every pair is differenced and gives its own spectra (the tutorial's "
                "way). Combined: all pairs are averaged into one image first, giving one "
                "spectrum. A: each frame on its own, without sky subtraction. A-Sky: each object "
                "frame minus a sky image made from the sky frames (object-sky nodding)."),
    "ex.combine_with": dict(
        label="Combine with",
        caption="For the Combined mode: the statistic and its rejection threshold.",
        details="The images are combined pixel by pixel with this statistic; the robust ones "
                "reject values more than the threshold (in sigma) from the mean."),
    "ex.abab": dict(
        label="ABBA to ABAB, so the pairs add up",
        caption="Reorders ABBA frames so every pair is A then B.",
        details="With ABBA nodding the second pair is B-A; reordering makes every difference "
                "A-B so they can be averaged."),
    "ex.cisky": dict(
        label="Subtract the sky level column by column first",
        caption="Removes any sky left in each order before combining.",
        details="The median of each column inside each order is subtracted."),
    "ex.naps": dict(
        label="Apertures",
        caption="2 for one point source in an A-B pair (its positive and negative image).",
        default="2",
        details="The number of peaks to find in the spatial profile. Two objects on the slit in "
                "an A-B pair give 4."),
    "ex.find": dict(
        label="Positions",
        caption="Auto finds the peaks; Guess refines your positions; Fix uses them as typed.",
        default="Auto",
        details="Auto takes the strongest peaks of the profile and fits each. Guess starts the "
                "fit from the positions you type (arcsec along the slit). Fix uses the typed "
                "positions without fitting."),
    "ex.positions": dict(
        label="",
        caption="Positions in arcsec along the slit, for Guess and Fix.",
        details="Read them off the profile plot: the horizontal axis is the slit position."),
    "ex.radius": dict(
        label="Radius",
        caption="Half-width of the aperture, in arcsec.",
        default="1.0 (PSAPRAD)",
        details="Light within this distance of the aperture centre is summed, partial pixels "
                "included. Comparable to the seeing: wider collects more light but more noise."),
    "ex.background": dict(
        label="Background",
        caption="Fitted between start and start+width arcsec from each aperture.",
        default="start 2.5, width 5.0, degree 1",
        details="The difference image still holds a little sky. It is fitted across the slit "
                "with a polynomial of this degree in the background region on either side of "
                "the aperture, and subtracted before summing."),
    "ex.trace": dict(
        label="Trace the apertures",
        caption="Follows the spectrum along the order; off keeps a fixed position.",
        default="on",
        details="The aperture centre is measured every few columns along the order and fitted "
                "with a polynomial, so the aperture follows any tilt or curvature of the "
                "spectrum on the detector."),
    "ex.signs": dict(
        label="Signs",
        caption="Blank: from the profile. Or type them, e.g. +, -.",
        details="Whether each aperture is the positive (A) or negative (B) image. Found from "
                "the profile; typed signs override (the IDL 'Override Signs')."),
    "ex.store": dict(
        label="Store the trace",
        caption="Saves the trace in cal/<source>_trace.dat for reuse.",
        details="A stored trace can be imported for a faint target, using a bright star's trace."),
    "ex.import": dict(
        label="Import",
        caption="A stored trace file, used instead of tracing.",
        details="For targets too faint to trace on their own."),
    "ex.trace_settings": dict(
        label="Trace settings",
        caption="Step and columns summed, fit degree, window and minimum significance.",
        default="step 5, columns 5, degree 1, window 3.0, min sigma 0.1",
        details="The trace samples the profile every `step` columns, summing `columns` columns, "
                "fits each peak, rejects points further than `window` pixels or fainter than "
                "`min sigma`, and fits the rest with a polynomial of `degree`."),
    "ex.xs_positions": dict(
        label="Positions",
        caption="Aperture centres in arcsec; blank uses the profile's peak.",
        details="For an extended source, apertures are placed where you choose. Several "
                "positions, separated by commas, give spatially resolved spectra."),
    "ex.xs_radii": dict(
        label="Radii",
        caption="One radius per position, in arcsec.",
        details="Half-widths of the apertures."),
    "ex.xs_bg": dict(
        label="Background",
        caption="Regions in arcsec, e.g. 2-6, 32-36.",
        details="Background regions for an extended source (see the option below)."),
    "ex.xs_bgdeg": dict(
        label="Degree",
        caption="Polynomial degree of the background fit.",
        details="Used only when the background is subtracted."),
    "ex.xs_background": dict(
        label="Subtract the background",
        caption="Off gives the IDL tool's results, which never subtract it.",
        default="off",
        details="In the IDL tool the extended-source background regions are drawn but the "
                "background is never subtracted. Off matches the IDL results; on fits and "
                "removes it as for point sources. For object-sky nodding the sky is already "
                "removed by the subtraction."),
    # ------------------------------------------------------------ 4. combine
    "cb.source": dict(
        label="Source",
        caption="The target whose spectra are combined.",
        details="All the spectra extracted for this target are listed below."),
    "cb.scaleorder": dict(
        label="Scale on order",
        caption="The order whose median level sets the scale factors.",
        default="5",
        details="Pick an order with good signal and little telluric absorption; the tutorial "
                "uses 4 or 5."),
    "cb.statistic": dict(
        label="Statistic",
        caption="How the spectra are combined at each wavelength.",
        default="Robust Weighted Mean",
        details="Robust Weighted Mean rejects outliers and weights by inverse variance (best "
                "with good errors). The robust means reject outliers without weights; Median is "
                "the most resistant; Sum adds the spectra."),
    "cb.thresh": dict(
        label="Threshold",
        caption="Rejection threshold, in sigma, for the robust statistics.",
        default="8",
        details="Values further than this from the mean are rejected before the final average."),
    "cb.shape": dict(
        label="Correct the spectral shapes",
        caption="Matches each spectrum's continuum shape to the mean's first.",
        default="off",
        details="For spectra whose shapes differ (slit losses vary with wavelength). Each is "
                "multiplied by the smooth ratio of the mean spectrum to it."),
    "cb.apertures": dict(
        label="Combine the apertures into one",
        caption="A and B images of the same source combined together.",
        default="off",
        details="Every aperture of every file goes into one stack, giving one spectrum."),
    "cb.spectra": dict(
        label="Spectra",
        caption="Untick a spectrum to leave it out (the prune step).",
        details="The scale factor of each is shown once combined."),
    # ------------------------------------------------------------ 5. telluric
    "tc.target": dict(
        label="Target",
        caption="The science target to correct.",
        details="Its combined spectrum is corrected."),
    "tc.standard": dict(
        label="A0 V standard",
        caption="The star used to measure the atmosphere and response.",
        details="Chosen in the plan as the closest in time and airmass; change it here."),
    "tc.bv": dict(
        label="B, V",
        caption="The standard's magnitudes: they set the flux scale.",
        details="Vega's model is scaled to the V magnitude and reddened by B-V. Look up fills "
                "them from SIMBAD."),
    "tc.method": dict(
        label="Method",
        caption="Deconvolution builds the kernel from the star's own line.",
        default="Deconvolution",
        details="Deconvolution derives the convolution kernel from a hydrogen line of the star "
                "(recommended). IP uses a stored instrument profile, available only for 0.3, "
                "0.5, 0.8, 1.6 and 3.0 arcsec slits; for TripleSpec's 1.1 arcsec slit the IDL "
                "tool stops, and so does this one."),
    "tc.order": dict(
        label="Order",
        caption="The order holding the H I line the kernel is built on.",
        default="6 (Pa delta, 1.005 um)",
        details="Pa delta in order 6 is isolated and free of strong telluric absorption."),
    "tc.continuum": dict(
        label="Continuum",
        caption="Windows either side of the line, in um, e.g. 0.99-0.997, 1.013-1.02.",
        details="The star is normalised by a polynomial fitted in these windows before the line "
                "is compared with Vega's. Drag on the plot with the field focused to set one."),
    "tc.degree": dict(
        label="Degree",
        caption="Degree of the continuum polynomial.",
        default="IDL panel 1; the tutorial uses 3",
        details="Degree 3 follows the continuum better around Pa delta (smaller residuals)."),
    "tc.line": dict(
        label="Line",
        caption="The line window in um, e.g. 0.997-1.013.",
        details="The region over which the kernel is derived."),
    "tc.units": dict(
        label="Units",
        caption="Flux units of the corrected spectrum.",
        default="ergs s-1 cm-2 A-1",
        details="The conversion is applied to the telluric spectrum."),
    "tc.prop": dict(
        label="Propagate the standard's uncertainty",
        caption="The standard's noise is added to the target's errors.",
        default="on",
        details="Off keeps only the target's own noise."),
    "tc.tension": dict(
        label="Tension",
        caption="Stiffness of the spline joining the line scales.",
        default="10",
        details="Higher tension gives straighter segments between the control points."),
    "tc.vshift": dict(
        label="Vega velocity",
        caption="Blank: as measured. Or type a velocity in km/s.",
        details="The radial velocity applied to Vega's model so its lines match the star's."),
    "tc.hl_order": dict(
        label="Order",
        caption="The order whose hydrogen-line scales are shown.",
        details="Each order has its own control points."),
    "tc.estimate": dict(
        label="Estimate the scale by dragging on the plot",
        caption="Fits the line under the drag and sets its scale (the IDL 'e' key).",
        details="A Lorentzian fit to the line gives the scale that removes it."),
    "tc.cut": dict(
        label="Cut",
        caption="Regions bridged by a straight line, e.g. 1.81-1.83.",
        details="For a feature in the telluric spectrum that should not be there."),
    "tc.shift": dict(
        label="Shift",
        caption="Fraction-of-a-pixel shifts that best align the correction.",
        details="The telluric spectrum is shifted from -1.5 to +1.5 pixels and the shift that "
                "leaves least structure in the region is kept."),
    # ------------------------------------------------------------ 6. merge
    "mg.source": dict(
        label="Source",
        caption="The corrected target whose orders are merged.",
        details=""),
    "mg.anchor": dict(
        label="Anchor",
        caption="The order the others are joined to.",
        default="3",
        details="Its flux level is kept; the others are added to it in turn."),
    "mg.table": dict(
        label="Orders",
        caption="Order, scale (a number or 'auto'), region for auto, trims.",
        details="A scale of 1 keeps the order as it is. 'auto' scales it to the merged spectrum "
                "over the region. A trim, e.g. 'green right 1.85', cuts the order or the merged "
                "spectrum at a wavelength."),
    "mg.region": dict(
        label="Region",
        caption="Drag on the plot: the region for the selected row's auto scale.",
        details=""),
    "mg.weighted": dict(
        label="Weight the overlap by the errors",
        caption="The overlap is the error-weighted mean of the two orders.",
        default="on",
        details="Off gives the plain mean."),
    "mg.sum": dict(
        label="Sum the overlap instead of averaging",
        caption="For spectra whose orders split the flux.",
        default="off",
        details=""),
    # ------------------------------------------------------------ 7. clean
    "cl.source": dict(label="Source", caption="The merged target to clean.", details=""),
    "cl.smoothing": dict(
        label="Smoothing",
        caption="None, Savitzky-Golay (keeps line shapes) or Gaussian.",
        default="None",
        details="Smoothing lowers noise at the cost of resolution."),
    "cl.sg": dict(
        label="Width, degree",
        caption="Savitzky-Golay window in pixels and polynomial degree.",
        default="1.5 slit widths (3.6 px, which truncates to 3: no change), degree 2",
        details="The window must be wider than the degree."),
    "cl.fwhm": dict(
        label="FWHM",
        caption="Gaussian full width at half maximum, in pixels.",
        default="the slit width (2.4 px)",
        details=""),
    "cl.region": dict(
        label="Region",
        caption="Drag on the plot, then Fix (straight line) or Remove (NaN).",
        details=""),
}


def field(key):
    """The entry for `key`, with empty strings for anything not given."""
    d = dict(label="", caption="", details="", default="")
    d.update(FIELDS[key])
    return d


def tooltip(key):
    """The full explanation as a tooltip: details, then the IDL default."""
    d = field(key)
    t = d["details"] or d["caption"]
    if d["default"]:
        t += f"\n\nDefault: {d['default']}"
    return t
