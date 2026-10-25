# Every setting

Every control of the desktop interface, step by step: what it is, what it does, and where the IDL tool starts it. The same text appears in the interface, as the caption under each control and as its tooltip; this page is generated from it (`triplespec_reducer/gui/help.py`), so the two always agree.

## 1. Night

A TripleSpec night is a few hundred frames: lamp flats with the lamp on and off, arc lamps, science targets nodded along the slit, and A0 V standard stars taken close in time and airmass to each target. Everything after this step needs to know which frames belong together: which are flats, which form the A-B pairs of each target, which standard corrects which target, and where each target's wavelength solution comes from.

| Setting | What it does | Default |
|---|---|---|
| **Raw frames** | The folder of the night's raw FITS files. It is only read. Every SPEC_*.fits file in this folder is read for its header. Nothing in it is written, moved or renamed, so there is no need to make a copy first. | — |
| **Output** | Where everything is written: cal/ for calibrations, proc/ for spectra. The plan (night.json), the calibration files (cal/), every spectrum (proc/) and the log of every step (qa/) go here. Each file is named after its target (sourcename_step.fits), so everything made for one object sorts together. | — |
| **Look up the standards' B, V and type on SIMBAD** | Needs the internet; without it, type B and V in the table. The B and V magnitudes of the A0 V standard set the absolute flux scale of the corrected spectrum, and the spectral type is checked (the method assumes an A0 V star). | — |

*Reproduces:* xspextool's file lists, typed by hand in the IDL tool; here read from the headers.

## 2. Calibrate

Two calibrations are needed before any spectrum can be extracted. The flat field corrects the pixel-to-pixel sensitivity of the detector and tells the software where each of the five cross-dispersed orders falls on it. The wavelength calibration maps every pixel to a wavelength and to a position along the slit, because TripleSpec's orders are curved and its slit image is tilted.

| Setting | What it does | Default |
|---|---|---|
| **Divide the images by the flat** | On in the IDL tool; turn off only to inspect the raw response. Every frame is divided by the normalised flat before extraction, which removes the pixel-to-pixel sensitivity differences of the detector. | on (ts4_soar.dat FLATFIELD=1) |
| **Flat name** | The flat's file name in cal/ (xspextool's Output Flat Name). A different name keeps more than one flat for the night, e.g. one per flat set. | flat |
| **Calibrate on** | For the selected solutions: as planned, the sky's OH lines, or the CuHeAr arc lines. xspextool leaves the choice of frames and line list free; the plan uses each target's own sky, and the arcs only when no target has long enough exposures. | As planned |
| **Frames** | For the selected solutions: the frames to calibrate on, as frame numbers. Blank uses all of the target's frames (or all the arcs). Any list may be given, e.g. part of a sequence, as in xspextool's Wavecal Images field. | — |
| **Short lines** | Leave out a line traced on too few points instead of stopping. Spextool stops the whole calibration when a line along the slit has fewer than two points (its POLY_FIT fails). On, such a line is left out and listed: this is not what the IDL tool does. | off |
| **Find the offset automatically (Auto X-Correlate)** | The spectrum is cross-correlated with the stored template. The one-dimensional sky (or arc) spectrum of order 4 is cross-correlated with the template the line list was made on; the peak gives the pixel offset that places every known line. Turn it off to type the offset yourself when the correlation fails (very few lines, a shifted grating). | on |
| **Offset (px)** | Used only when the automatic offset is off. The shift, in pixels, between this night's spectrum and the template, in order 4. 'Measure on a line' finds it as XGETOFFSET does: drag across a line in the template, then across the same line in the observed spectrum; each is fitted with a Gaussian on a straight line, and the offset is the difference of the two centres. | — |

*Reproduces:* xspextool's calibration step (flat field and 2D wavelength calibration; Cushing, Vacca & Rayner 2004).

## 3. Extract

The infrared sky is bright and variable; it is removed by nodding the target between two positions on the slit (A and B) and subtracting one frame from the other. In A-B the target appears twice, once positive and once negative, and the sky cancels. Extraction turns each pair into a spectrum per order by adding up the light inside an aperture around each image of the target.

| Setting | What it does | Default |
|---|---|---|
| **Source** | The target whose frames are extracted. Standards and science targets are both extracted the same way. | — |
| **Mode** | How the frames are turned into images before extraction. A-B pairs: every pair is differenced and gives its own spectra (the tutorial's way). Combined: all pairs are averaged into one image first, giving one spectrum. A: each frame on its own, without sky subtraction. A-Sky: each object frame minus a sky image made from the sky frames (object-sky nodding). | A-B pairs |
| **Combine with** | For the Combined mode: the statistic and its rejection threshold. The images are combined pixel by pixel with this statistic; the robust ones reject values more than the threshold (in sigma) from the mean. | — |
| **ABBA to ABAB, so the pairs add up** | Reorders ABBA frames so every pair is A then B. With ABBA nodding the second pair is B-A; reordering makes every difference A-B so they can be averaged. | — |
| **Subtract the sky level column by column first** | Removes any sky left in each order before combining. In each column of each order, the robust weighted mean of the pixels (outliers beyond 3 sigma left out) is subtracted (MC_SUBSPECSKY). | — |
| **Apertures** | 2 for one point source in an A-B pair (its positive and negative image). The number of peaks to find in the spatial profile. Two objects on the slit in an A-B pair give 4. | 2 |
| **Positions** | Auto finds the peaks; Guess refines your positions; Fix uses them as typed. Auto takes the strongest peaks of the profile and fits each. Guess starts the fit from the positions you type (arcsec along the slit). Fix uses the typed positions without fitting. | Auto |
| **(positions)** | Positions in arcsec along the slit, for Guess and Fix. Read them off the profile plot: the horizontal axis is the slit position. | — |
| **Radius** | Half-width of the aperture, in arcsec. Light within this distance of the aperture centre is summed, partial pixels included. Comparable to the seeing: wider collects more light but more noise. | 1.0 (PSAPRAD) |
| **Background** | Fitted between start and start+width arcsec from each aperture. The difference image still holds a little sky. It is fitted across the slit with a polynomial of this degree in the background region on either side of the aperture, and subtracted before summing. | start 2.5, width 5.0, degree 1 |
| **Background** | Subtract the background around each aperture, and extract optimally. As in the IDL tool, where PSBGSUB=1 also switches optimal extraction on: each wavelength's spatial profile, fitted along the order, weights the pixels within the PSF radius (Horne 1986), and the flux of the whole profile is returned. Off gives a plain sum inside the aperture, with no background. | on |
| **PSF radius** | Radius in arcsec over which the profile is fitted for optimal extraction. The aperture radius must not be larger than it. | 2.2 (PSPSFRAD) |
| **Median profile** | Use one median profile for the whole order instead of a profile that changes along it. xspextool's Median Profile. | off |
| **Sum with background** | Keep the background but sum the aperture instead of extracting optimally. Not available in the IDL tool for TripleSpec, where background subtraction always comes with optimal extraction; results have no IDL counterpart. | off |
| **Orders** | Extract only these orders, e.g. 3,5,6; blank for all. xspextool's PS Orders and XS Orders. | all |
| **No wavelength calibration** | Extract in pixels, without a wavecal file (xspextool with a blank Full Wavecal Name). The spectra are then in pixels, the bad-pixel mask and the threshold below are used, and the extended-source background is subtracted when its switch is on. | off |
| **Bad pixel threshold** | Sigma for replacing bad pixels, when extracting without a wavecal file. | 7 (BADPIXELTHRESH) |
| **XS background** | Extended-source background subtraction, as xspextool's XS BG Subtraction. Only takes effect without a wavecal file: with one, the IDL tool never subtracts it (see the option below). | on (XSBGSUB=1) |
| **Trace the apertures** | Follows the spectrum along the order; off keeps a fixed position. The aperture centre is measured every few columns along the order and fitted with a polynomial, so the aperture follows any tilt or curvature of the spectrum on the detector. | on |
| **Signs** | Blank: from the profile. Or type them, e.g. +, -. Whether each aperture is the positive (A) or negative (B) image. Found from the profile; typed signs override (the IDL 'Override Signs'). | — |
| **Store the trace** | Saves the trace in cal/<source>_trace.dat for reuse. A stored trace can be imported for a faint target, using a bright star's trace. | — |
| **Import** | A stored trace file, used instead of tracing. For targets too faint to trace on their own. | — |
| **Trace settings** | Step and columns summed, fit degree, window and minimum significance. The trace samples the profile every `step` columns, summing `columns` columns, fits each peak, rejects points further than `window` pixels or fainter than `min sigma`, and fits the rest with a polynomial of `degree`. The instrument file's min sigma, 0.1, is read by the IDL tool as an integer: 0. | step 5, columns 5, degree 1, window 3.0, min sigma 0 |
| **Positions** | Aperture centres in arcsec; blank uses the profile's peak. For an extended source, apertures are placed where you choose. Several positions, separated by commas, give spatially resolved spectra. | — |
| **Radii** | One radius per position, in arcsec. Half-widths of the apertures. | — |
| **Background** | Regions in arcsec, e.g. 2-6, 32-36. Background regions for an extended source (see the option below). | — |
| **Degree** | Polynomial degree of the background fit. Used only when the background is subtracted. | — |
| **Subtract the background** | Off gives the IDL tool's results, which never subtract it. In the IDL tool the extended-source background regions are drawn but the background is never subtracted. Off matches the IDL results; on fits and removes it as for point sources. For object-sky nodding the sky is already removed by the subtraction. | off |

*Reproduces:* xspextool's extraction ('Do All Steps': profiles, apertures, trace, extraction; Cushing, Vacca & Rayner 2004).

## 4. Combine

A target observed in many pairs gives many spectra; they are combined into one with a higher signal-to-noise ratio. Their flux levels differ slightly because seeing, guiding and transparency change during the sequence, so they are first scaled to a common level.

| Setting | What it does | Default |
|---|---|---|
| **Source** | The target whose spectra are combined. All the spectra extracted for this target are listed below. | — |
| **Scale on order** | The order whose median level sets the scale factors. Pick an order with good signal and little telluric absorption; the tutorial uses 4 or 5. | 5 |
| **Scale the spectra** | Bring the spectra to a common level before combining (XCOMBSPEC 'Scale Spectra'). Off combines them as extracted, as XCOMBSPEC does when Scale is never pressed; the history then says nothing about scaling. | on, as in the tutorial |
| **Scale to** | The median spectrum, or one chosen spectrum (its number in the list). XMC_SCALESPEC's 'Scale to'. The reference spectrum must not be left out. | Median |
| **Scale range** | The wavelengths the scale factors are measured over. XMC_SCALESPEC's 's' selection. | the middle 80% of the order |
| **Leave out of an order** | Spectra to drop from one order only, e.g. '5: 2 3' (numbers as listed). XCOMBSPEC's prune step (XMC_MASKSPEC). For another aperture: '5/2: 1'. As in XCOMBSPEC, the unweighted means then still divide by the full number of spectra, and the median repeats the last one; the program reproduces this. | — |
| **Mask ranges** | Wavelength ranges masked in one order, e.g. '5: 1.30-1.32'. XMC_MASKSPEC's region mask: those pixels of every kept spectrum are ignored. | — |
| **Basic header** | Write the shorter header of xcombspec /BASIC. | off |
| **Statistic** | How the spectra are combined at each wavelength. Robust Weighted Mean rejects outliers and weights by inverse variance (best with good errors). The robust means reject outliers without weights; Median is the most resistant; Sum adds the spectra. | Robust Weighted Mean |
| **Threshold** | Rejection threshold, in sigma, for the robust statistics. Values further than this from the mean are rejected before the final average. | 8 |
| **Correct the spectral shapes** | Matches each spectrum's continuum shape to the mean's first. For spectra whose shapes differ (slit losses vary with wavelength). Each is multiplied by the smooth ratio of the mean spectrum to it. | off |
| **Combine the apertures into one** | A and B images of the same source combined together. Every aperture of every file goes into one stack, giving one spectrum. | off |
| **Spectra** | Untick a spectrum to leave it out (the prune step). The scale factor of each is shown once combined. | — |

*Reproduces:* xcombspec (Cushing, Vacca & Rayner 2004).

## 5. Telluric

Earth's atmosphere absorbs much of the near infrared (water and carbon dioxide), and the instrument's response varies strongly with wavelength. Both are measured with an A0 V star observed near the target: its true spectrum is known (it looks like Vega, whose spectrum is modelled), so dividing the model by the observed star gives the correction that turns the target's counts into flux and removes the atmosphere. This is Vacca, Cushing and Rayner's (2003) method.

| Setting | What it does | Default |
|---|---|---|
| **Target** | The science target to correct. Its combined spectrum is corrected. | — |
| **Standard** | The star used to measure the atmosphere and response. Chosen in the plan as the closest in time and airmass; change it here. | — |
| **B, V** | The standard's magnitudes: they set the flux scale. Vega's model is scaled to the V magnitude and reddened by B-V. Look up fills them from SIMBAD. | — |
| **Correct with** | An A0 V star (xtellcor), a star of another type (xtellcor_basic), or a telluric spectrum written before (xtellcor_finish). A0 V: Vega's model, convolved with a kernel from the star's own line, gives the atmosphere and the flux scale. Other standard: the target is divided by the star and, optionally, multiplied by a blackbody at its temperature scaled to its V. Telluric file: an existing <source>_tellspec.fits is applied, with its own shift. | A0 V standard |
| **Blackbody** | For another standard: restore its continuum with a blackbody of this temperature. Scaled to the V magnitude above. Off leaves the plain ratio of target to star. | on |
| **Telluric file** | For the telluric-file route: the telluric spectrum to apply. Its standard's name, magnitudes and airmass are carried into the header; a difference of more than 0.1 in airmass is reported. | — |
| **Model** | Also write the convolved, scaled Vega model (XTELLCOR 'Model'). | off |
| **Method** | Deconvolution builds the kernel from the star's own line. Deconvolution derives the convolution kernel from a hydrogen line of the star (recommended). IP uses a stored instrument profile, available only for 0.3, 0.5, 0.8, 1.6 and 3.0 arcsec slits; for TripleSpec's 1.1 arcsec slit the IDL tool stops, and so does this one. The IDL tool does not offer IP for TS4 at all: its Method buttons are only shown for SpeX. | Deconvolution |
| **Order** | The order holding the H I line the kernel is built on. Pa delta in order 6 is isolated and free of strong telluric absorption. | 6 (Pa delta, 1.005 um) |
| **Continuum** | Windows either side of the line, in um, e.g. 0.99-0.997, 1.013-1.02. The star is normalised by a polynomial fitted in these windows before the line is compared with Vega's. Drag on the plot with the field focused to set one. | — |
| **Degree** | Degree of the continuum polynomial. Degree 3 follows the continuum better around Pa delta (smaller residuals). | IDL panel 1; the tutorial uses 3 |
| **Line** | The line window in um, e.g. 0.997-1.013. The region over which the kernel is derived. | — |
| **Units** | Flux units of the corrected spectrum. The conversion is applied to the telluric spectrum. | ergs s-1 cm-2 A-1 |
| **Propagate the standard's uncertainty** | The standard's noise is added to the target's errors. Off keeps only the target's own noise. IDL has no such switch: it always propagates. | on |
| **Tension** | Stiffness of the spline joining the line scales. Higher tension gives straighter segments between the control points. IDL has no control for it: it is always 10. | 10 |
| **Vega velocity** | Blank: as measured. Or type a velocity in km/s. The radial velocity applied to Vega's model so its lines match the star's. | — |
| **Order** | The order whose hydrogen-line scales are shown. Each order has its own control points. | — |
| **Estimate the scale by dragging on the plot** | Fits the line under the drag and sets its scale (the IDL 'e' key). A Lorentzian fit to the line gives the scale that removes it. | — |
| **Cut** | Regions bridged by a straight line, e.g. 1.81-1.83. For a feature in the telluric spectrum that should not be there. | — |
| **Shift** | Fraction-of-a-pixel shifts that best align the correction, per order and aperture. The telluric spectrum is shifted from -1.5 to +1.5 pixels and the shift that leaves least structure in the region is kept. | — |

*Reproduces:* xtellcor (Vacca, Cushing & Rayner 2003; Cushing, Vacca & Rayner 2004).

## 6. Merge

TripleSpec records five orders (3 to 7), each covering part of 0.94-2.47 um, with overlaps between neighbours. Merging joins them into one continuous spectrum.

| Setting | What it does | Default |
|---|---|---|
| **Source** | The corrected target whose orders are merged. | — |
| **Anchor** | The order the others are joined to. Its flux level is kept; the others are added to it in turn. | 3 |
| **Orders** | Order, scale (a number or 'auto'), region for auto, trims; added top to bottom. A scale of 1 keeps the order as it is. 'auto' scales it to the merged spectrum over the region, found before any trim, as XMERGEORDERS does. A trim cuts the order being added (Green) or the merged spectrum (Blue) at a wavelength, e.g. 'green right 1.47'; one of each, separated by ';', may be given. Each order must be next to the orders merged so far, as XMERGEORDERS offers only those. | — |
| **Region** | Drag on the plot: the region for the selected row's auto scale. XMERGEORDERS has no default region: 'auto' without one is refused. | — |
| **Aperture** | Which aperture's orders are merged. Each aperture is merged on its own. With one aperture the file is <source>_merged.fits; with several, <source>_merged_ap<n>.fits. | 1 |
| **Weight the overlap by the errors** | The overlap is the error-weighted mean of the two orders. Off is not available in Spextool 4.1: its XMERGEORDERS stops on an undefined variable, so the plain mean it gives here has no IDL counterpart. | on |
| **Sum the overlap instead of averaging** | For spectra whose orders split the flux. Not in Spextool 4.1. XMERGEORDERS 4.1 switches weighting off for a sum and then stops on the same undefined variable, so a sum has no IDL counterpart. | off |

*Reproduces:* xmergeorders (Cushing, Vacca & Rayner 2004).

## 7. Clean & save

Optional last touches: smoothing, and repairing or removing bad regions, such as the opaque gaps between the atmospheric windows or a residual from a bad pixel.

| Setting | What it does | Default |
|---|---|---|
| **Source** | The target whose spectrum is cleaned. | — |
| **File** | Which of its spectra: merged, one aperture's merged, corrected, or combined. XCLEANSPEC opens any Spextool spectrum file. | merged |
| **Order, aperture** | The spectrum of the file that is shown and edited. A merged file has one spectrum; a corrected or combined one has one per order and aperture. | — |
| **Edit** | Fix or remove the flux and its error together, or the error only. XCLEANSPEC's spectrum type: 'Error' changes only the uncertainties. | Flux |
| **Smoothing** | None, Savitzky-Golay (keeps line shapes) or Gaussian. Smoothing lowers noise at the cost of resolution. | None |
| **Width, degree** | Savitzky-Golay window in pixels and polynomial degree. The degree must be smaller than the window (XCLEANSPEC refuses it otherwise). | 1.5 slit widths (from the file's SLTW_PIX: 3.6 px, which truncates to 3: no change), degree 2 |
| **FWHM** | Gaussian full width at half maximum, in pixels. | the slit width (the file's SLTW_PIX, 2.4 px) |
| **Region** | Drag on the plot, then Fix (straight line) or Remove (NaN); Keep makes the edits stick, Undo goes back to the last Keep. As in XCLEANSPEC ('f', 'r', 's', 'u'), edits are drawn on a working copy: only kept edits are smoothed and written. Smoothing is allowed once per file opened, and starts from the kept spectrum. | — |

*Reproduces:* xcleanspec (Cushing, Vacca & Rayner 2004); Savitzky & Golay (1964).

