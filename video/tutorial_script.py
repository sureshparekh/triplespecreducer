"""
Film 2: a hands-on tutorial, reducing NGC 134 from scratch in the interface.

Each beat: `step` (shown in the film's step chip), `lines` (the narration,
as displayed; video/voice.py speaks them), and `actions`, done in order by
video/record_gui.py. An action may carry `at`: a word of the narration it
waits for (the first occurrence, or "word#2" for the second), or "end" for
after the narration. Targets are TSR attribute names, ("button", text),
("row", table, text-in-first-column), ("rail", step index), ("about", step,
0|1), or None.

Actions:
  move T          the cursor glides to T
  click T         move, press, and invoke T
  type T text     click into the field, select all, type the text
  select T text   a combo box set to text (with a click on it)
  step i          click step i in the step bar
  scroll i T      scroll step i's panel to T
  wait            the running job finishes (the waiting is recorded)
  pause s         hold s seconds
"""

RAW = "data/raw/2025-11-02"
OUT = "reduced/tutorial_film"
SCI, STD = "NGC134", "HD225200_1"

BEATS = [
    dict(id="intro", step=None, lines=[
        "Let's reduce one galaxy, NGC 134, from the raw frames to its final spectrum."],
        actions=[]),
    dict(id="launch", step=None, lines=[
        "Open the program with the command TSR."],
        actions=[("launch",)]),
    dict(id="layout", step=1, lines=[
        "Across the top is the step bar: the seven steps of the reduction, from Night to Clean and save.",
        "The panel on the left holds the current step, and the stage on the right shows whatever the "
        "step is working on."],
        actions=[("move", ("rail", 0)), ("move", ("rail", 6), "Clean"),
                 ("move", "ed_raw", "panel"), ("move", "STAGE", "stage")]),
    dict(id="about", step=1, lines=[
        "Every step opens with a panel called About this step. It says why the step is needed, "
        "what it does in order, and what to check before moving on.",
        "Every setting below it has a caption, and a tooltip with the full explanation."],
        actions=[("click", ("about", 0, 0), "what"), ("pause", 2.5),
                 ("click", ("about", 0, 0), "check"), ("pause", 0.6)]),
    dict(id="raw", step=1, lines=[
        "First, the raw frames: type, or browse to, the folder that holds the night's files.",
        "The program only reads this folder. Nothing in it is ever written, renamed or moved."],
        actions=[("type", "ed_raw", RAW, "type")]),
    dict(id="out", step=1, lines=[
        "Then the output folder, where everything will be written: calibrations in a folder called "
        "cal, and spectra in a folder called proc."],
        actions=[("type", "ed_out", OUT, "Then")]),
    dict(id="read", step=1, lines=[
        "Leave the SIMBAD box ticked, so the standard stars' magnitudes and spectral types are "
        "looked up online, and press Read the night."],
        actions=[("move", "cb_simbad", "Leave"), ("click", "btn_read", "press"), ("wait",), ("pause", 0.5)]),
    dict(id="plan", step=1, lines=[
        "Two hundred and thirty-five frames were read and sorted into ten targets.",
        "Each row is one target: whether it is science or a standard star, a point source or an "
        "extended one, its nod pattern, how many frames and pairs it has, its airmass, where its "
        "wavelength solution comes from, and the standard that will correct it."],
        actions=[("scroll", 0, "tbl_seq"), ("pause", 1.0)]),
    dict(id="row", step=1, lines=[
        "Our galaxy is NGC 134. It fills the slit, so it was nodded off the galaxy, to blank sky: an "
        "extended source, with eleven frames making seven pairs.",
        "Its standard was chosen automatically: of the standards within a tenth of its airmass, the "
        "one observed closest in time, here HD 225200, just before it. Any other can be picked in the table.",
        "Selecting the row shows its first A minus B image on the stage."],
        actions=[("click", ("row", "tbl_seq", SCI), "Our")]),
    dict(id="notes", step=1, lines=[
        "Under the table, the notes report every decision that needed a judgement.",
        "Here, SIMBAD lists the standard as A1 V rather than A0 V, close enough for the method, "
        "and a second pointing on the galaxy was kept apart as a target of its own."],
        actions=[("move", "lbl_notes", "Under")]),
    dict(id="save", step=1, lines=[
        "Any cell can be changed if a decision is wrong. This plan is right, so we save it and go "
        "to the next step."],
        actions=[("click", ("button", "Save the plan"), "save"), ("click", ("button", "Next: Calibrate"), "next")]),
    dict(id="cal", step=2, lines=[
        "Step two, Calibrate. Two calibrations come before any spectrum: the flat field, for the "
        "detector's pixel-to-pixel sensitivity, and the wavelength calibration, which gives every "
        "pixel a wavelength."],
        actions=[]),
    dict(id="flat", step=2, lines=[
        "Press Build the flat. Fifty lamp-on frames are combined and the lamp-off frames "
        "subtracted; the edges of the five orders are traced; and each order is divided by a "
        "smooth model of itself.",
        "This takes about half a minute."],
        actions=[("scroll", 1, "lbl_flat"), ("click", ("button", "Build the flat"), "Press"), ("wait",)]),
    dict(id="flat_seen", step=2, lines=[
        "Here is the normalised flat. It is close to one inside the orders, and the fitted edges "
        "follow every order across the detector.",
        "Order seven, at the top, is short, because it runs off the detector."],
        actions=[]),
    dict(id="wave", step=2, lines=[
        "Below are the wavelength solutions. Each science target is calibrated on its own sky lines, "
        "the hydroxyl emission recorded in its frames; the standards borrow the solution of the "
        "target nearest in time.",
        "Select NGC 134, and press Calibrate the selected."],
        actions=[("scroll", 1, "tbl_wc"), ("click", ("row", "tbl_wc", SCI), "Select"),
                 ("click", ("button", "Calibrate the selected"), "press")]),
    dict(id="wave_wait", step=2, wait_job=True, lines=[
        "While it runs: a sky spectrum is extracted down the middle of every order and matched to a "
        "template; the sky lines are identified, and fitted with one polynomial across all five "
        "orders; and then each line is followed along the slit, to map how the slit's image is tilted.",
        "This is the longest step, about a minute."],
        actions=[("wait",)]),
    dict(id="wave_done", step=2, lines=[
        "Done: between thirty and seventy-six lines per order, and an RMS of a third of an angstrom.",
        "The stage shows the wavelength of every pixel."],
        actions=[("pause", 2.0)]),
    dict(id="ext", step=3, lines=[
        "Step three, Extract. We start with the standard star, because a point source shows every "
        "part of extraction. The mode stays at A minus B pairs. Press Preview the first pair."],
        actions=[("click", ("button", "Next: Extract")), ("select", "cb_ex_src", STD, "standard"),
                 ("click", ("button", "Preview the first pair"), "Press"), ("wait",)]),
    dict(id="profiles", step=3, lines=[
        "These are the spatial profiles, one per order: the light along the slit.",
        "The star appears twice, positive where it sat in A, and negative in B.",
        "The green bands are the apertures, one arcsecond around each peak, and the red bands are "
        "where the leftover background is measured. They must not overlap."],
        actions=[]),
    dict(id="ext_std", step=3, lines=[
        "Press Extract every pair. Each pair is differenced, divided by the flat, and traced along "
        "the order. With the background subtracted, the star is extracted optimally: its own profile "
        "weights every pixel near it."],
        actions=[("click", ("button", "Extract every pair"), "Press"), ("wait",)]),
    dict(id="ext_sci", step=3, lines=[
        "Now the galaxy. Select NGC 134, and preview its first pair."],
        actions=[("select", "cb_ex_src", SCI, "Select"), ("click", ("button", "Preview the first pair"), "preview"), ("wait",)]),
    dict(id="nucleus", step=3, lines=[
        "One peak: the nucleus, nineteen arcseconds along the slit, and the aperture goes there.",
        "To extract other parts of the galaxy, type their positions in the extended-source box."],
        actions=[("pause", 1.5), ("scroll", 2, "ed_xspos", "To"), ("move", "ed_xspos", "type")]),
    dict(id="ext_sci_run", step=3, lines=[
        "Extract every pair: seven pairs, about twelve seconds each."],
        actions=[("click", ("button", "Extract every pair"), "Extract"), ("wait",)]),
    dict(id="spectrum", step=3, lines=[
        "Here is one extracted spectrum, each colour one order. It is still in counts per second, "
        "and it still carries the atmosphere: the deep gaps near one point four and one point nine "
        "microns are water vapour."],
        actions=[]),
    dict(id="comb", step=4, lines=[
        "Step four, Combine. The spectra of a target are scaled to a common level and averaged. "
        "The standard first."],
        actions=[("click", ("button", "Next: Combine")), ("select", "cb_cb_src", STD, "standard")]),
    dict(id="comb_set", step=4, lines=[
        "The scale factors come from order five, and the spectra are combined with a robust "
        "weighted mean, which rejects outliers beyond eight sigma and weights the rest by their errors.",
        "Every spectrum is listed with a tick; untick one to leave it out. Press Combine."],
        actions=[("move", "cb_scaleorder", "scale"), ("move", "cb_stat", "robust"), ("move", "lst_cb", "Every"),
                 ("click", ("button", "Combine"), "Press"), ("wait",)]),
    dict(id="comb_sci", step=4, lines=[
        "The grey lines are the individual spectra, and the coloured line their combination. "
        "And the same for the galaxy."],
        actions=[("select", "cb_cb_src", SCI, "galaxy"), ("click", ("button", "Combine"), "end"), ("wait",)]),
    dict(id="tell", step=5, lines=[
        "Step five, the telluric correction, the heart of the reduction.",
        "The A0 V star went through the same atmosphere minutes earlier, and its true spectrum is "
        "known. So a model of it, divided by what was observed of it, is the correction for the galaxy.",
        "The target and its standard are already chosen, with the star's magnitudes from SIMBAD. "
        "Press Load the two spectra."],
        actions=[("click", ("button", "Next: Telluric"), "Step"), ("click", ("button", "Load the two spectra"), "Press"), ("wait",)]),
    dict(id="pad", step=5, lines=[
        "The stage shows the star's Paschen delta line, in order six, with the continuum windows "
        "either side in green and the line in blue.",
        "Vega's model has to be blurred to match this line exactly. Press Build the kernel."],
        actions=[("pause", 2.0), ("scroll", 4, "sp_korder", "Vega"),
                 ("click", ("button", "Build the kernel"), "Press"), ("wait",)]),
    dict(id="kernel", step=5, lines=[
        "Vega's line, convolved with the kernel, now lies on the star's.",
        "The caption gives the quality of the match: a maximum deviation of about two percent, and "
        "an RMS of six tenths of a percent. That is a good kernel."],
        actions=[("pause", 2.0)]),
    dict(id="hlines", step=5, lines=[
        "The telluric spectrum is now computed for every order. Each hydrogen line has a control "
        "point that scales Vega's line there.",
        "If the correction showed a bump at a line, we would lower its scale; a dip, raise it. Here "
        "it is smooth, and nothing needs changing."],
        actions=[("scroll", 4, "tbl_cp"), ("pause", 1.0)]),
    dict(id="shift", step=5, lines=[
        "Press Find every shift. The correction is slid by fractions of a pixel, order by order, to "
        "line up with the galaxy's own absorption features."],
        actions=[("scroll", 4, "tbl_shift"), ("click", ("button", "Find every shift"), "Press"),
                 ("wait",)]),
    dict(id="write", step=5, lines=[
        "And press Write the corrected spectrum. The galaxy is now in physical flux units, with the "
        "atmosphere removed."],
        actions=[("click", ("button", "Write the corrected spectrum"), "press"), ("wait",)]),
    dict(id="merge", step=6, lines=[
        "Step six, Merge: the five orders are joined into one spectrum, with order three as the "
        "anchor. Press Merge."],
        actions=[("click", ("button", "Next: Merge")), ("click", ("button", "Merge"), "Press"), ("wait",)]),
    dict(id="merged", step=6, lines=[
        "The orders meet without a step, because the telluric correction has already put them on "
        "one flux scale."],
        actions=[]),
    dict(id="clean", step=7, lines=[
        "Step seven, Clean and save. The stage shows the final spectrum and its signal to noise: one "
        "to two hundred per pixel, in the bands.",
        "Smoothing is optional. We leave it off, and press Save to write the final file."],
        actions=[("click", ("button", "Next: Clean & save")), ("click", ("button", "Save"), "press"), ("pause", 1.0)]),
    dict(id="products", step=None, lines=[
        "Finally, the Products page. Type a source's name, and every file made for it is listed, "
        "with its header, and the log of every step: the settings used, and what came out."],
        actions=[("click", ("seg", "Products"), "Products"), ("type", "ed_prod", SCI, "Type"), ("click", ("lastitem", "lst_prod"), "listed")]),
    dict(id="all", step=None, lines=[
        "Everything we did step by step, the button Reduce everything does for the whole night at once.",
        "And that is one galaxy, reduced from the raw frames to a calibrated spectrum."],
        actions=[("click", ("seg", "Reduce"), "Everything"), ("move", ("button", "Reduce everything"), "button"),
                 ("pause", 2.0)]),
]
