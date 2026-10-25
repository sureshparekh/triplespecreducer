"""
Film 1: the science of the reduction, and how the program does each step.

Each scene: `id`, `lines` (narration as displayed; video/voice.py speaks
them), and the visual it plays over is built by video/film/src/explainer.js
under the same id.
"""

SCENES = [
    dict(id="open", lines=[
        "TripleSpec Reducer turns a night of near-infrared spectra from the SOAR TripleSpec "
        "spectrograph into calibrated, telluric-corrected spectra.",
        "In this film we follow the light, from the telescope to the final spectrum, and see what "
        "each step of the reduction does, why it is needed, and how the program does it."]),
    dict(id="method", lines=[
        "The reduction follows the steps and the science of Spextool, the package Michael Cushing, "
        "William Vacca and John Rayner wrote for NASA's Infrared Telescope Facility, with the "
        "telluric correction method they published in 2003, and its module for TripleSpec.",
        "Run on the same frames, the program's results and those of the IDL tool are identical, "
        "bit for bit."]),
    dict(id="instrument", lines=[
        "TripleSpec is a cross-dispersed spectrograph on the four point one metre SOAR telescope.",
        "Light from a slit, thirty-eight arcseconds long, is spread out by a grating, and a prism "
        "separates the grating's orders. So one exposure holds five spectra, one above the other: "
        "orders three to seven, from zero point nine four to two point four seven microns, at a "
        "resolving power of about three and a half thousand."]),
    dict(id="frame", lines=[
        "This is a real frame. Each curved strip is one order, the image of the slit spread out in "
        "wavelength, and inside each strip the star appears twice, at the two positions it was "
        "observed in."]),
    dict(id="atmosphere", lines=[
        "Two things stand between the detector and the spectrum we want. The first is the atmosphere.",
        "Water vapour and carbon dioxide absorb much of the near infrared: between the bands almost "
        "nothing gets through, and inside them hundreds of narrow lines remain, changing with the "
        "airmass and the humidity."]),
    dict(id="sky", lines=[
        "The second is the sky's own glow. Hydroxyl molecules high in the atmosphere emit bright, "
        "narrow lines, far brighter than a faint galaxy, and they change from minute to minute.",
        "Those same lines, though, are at wavelengths known to a fraction of an angstrom, and they "
        "will calibrate our wavelengths."]),
    dict(id="nodding", lines=[
        "The sky is removed by nodding. The telescope observes the target at one position, A, then "
        "moves it to another, B, and the two frames are subtracted.",
        "Everything that is the same in both, the sky lines, the telescope's warm glow, the "
        "detector's dark current, cancels. A star nodded along the slit appears twice, once "
        "positive and once negative; a galaxy is nodded off, to blank sky, and appears once."]),
    dict(id="night", lines=[
        "So a night holds four kinds of frames: lamp flats, arc lamps, the science targets in their "
        "A and B positions, and A0 V standard stars, observed close in time and airmass to every target.",
        "The program reads every header, sorts the frames, pairs them, and proposes the whole plan "
        "in one table, where every decision can be checked and changed."]),
    dict(id="flat", lines=[
        "The first calibration is the flat field. A lamp lights the slit evenly; the lamp-off frames, "
        "which hold only the thermal glow, are subtracted.",
        "The edges of the five orders are traced, and each order is divided by a smooth model of "
        "itself, leaving only the pixel-to-pixel sensitivity of the detector, close to one, which "
        "every frame is divided by."]),
    dict(id="wave", lines=[
        "Next, the wavelengths. A sky spectrum is extracted from each order and cross-correlated with "
        "a template; the hydroxyl lines are identified, and one polynomial in column and order number "
        "is fitted to all of them, with an RMS of about a third of an angstrom.",
        "Then each line is followed along the slit, to map how the slit's image is tilted. The result "
        "is a wavelength, and a slit position, for every pixel."]),
    dict(id="extract", lines=[
        "Now extraction. For every pair, B is subtracted from A and the result divided by the flat.",
        "A spatial profile shows where the light falls along the slit. Apertures are placed on the "
        "target, and the light left between them is fitted as background and removed.",
        "A star is then measured optimally: its own profile weights every pixel near it, giving the "
        "whole star's flux at the best signal to noise. An extended source is summed inside its apertures."]),
    dict(id="combine", lines=[
        "The spectra of one target are then combined. They are first scaled to a common level, because "
        "the seeing and the guiding change during a sequence.",
        "Then, at every wavelength, a robust weighted mean rejects the outliers and weights the rest by "
        "their errors. The signal to noise grows as the square root of the number of spectra."]),
    dict(id="telluric", lines=[
        "The telluric correction is the heart of the reduction. What we recorded of the galaxy is its "
        "true spectrum, times the atmosphere, times the instrument.",
        "An A0 V star went through the same atmosphere and instrument minutes earlier, and its true "
        "spectrum is known: it looks like Vega's, which has been modelled in detail. So the model "
        "divided by the observed star is exactly the correction."]),
    dict(id="kernel", lines=[
        "The catch is the star's own deep hydrogen lines. Vega's model lines must first match them "
        "exactly, as TripleSpec saw them.",
        "Vega is shifted to the star's velocity, and a convolution kernel is found by dividing the "
        "Fourier transforms of the two lines. Convolved with that kernel, Vega's line lies on the "
        "star's, to within a percent or two."]),
    dict(id="tellspec", lines=[
        "With the kernel, Vega's model is rebuilt for every order, scaled to the star's magnitude, and "
        "divided by the observed star: the telluric spectrum.",
        "Slid by a fraction of a pixel to line up with the galaxy's features, it multiplies the "
        "galaxy's spectrum, removing the atmosphere and giving it physical flux units at once."]),
    dict(id="merge", lines=[
        "Last, the five orders are merged. They overlap at their ends, and are joined one by one, the "
        "overlaps averaged with weights from their errors, into one continuous spectrum."]),
    dict(id="products", lines=[
        "Every file is named after its source and its step, in Spextool's own format, and every step "
        "keeps a log of its settings and its results.",
        "The whole night can be reduced with one button, or step by step, with every step explaining "
        "itself."]),
    dict(id="validation", lines=[
        "Every product was compared with the IDL tool's on the same frames: the flat, the wavelength "
        "solution, every extracted spectrum, the combined, corrected and merged spectra, even every "
        "number in every header. They are identical, bit for bit.",
        "The one exception is the telluric kernel, made with IDL's own Fourier transform, which agrees "
        "to a few parts in a million.",
        "The documentation has the tutorial, the science of every step, and every source it rests on."]),
]
