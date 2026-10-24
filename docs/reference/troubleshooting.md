# Troubleshooting

### The window does not open, or opens without the logo

The interface needs PySide6 (`pip install -e .` installs it). On a remote machine run it with
X forwarding (`ssh -X`) or on a desktop session; `QT_QPA_PLATFORM=offscreen` only renders off
screen. The logo files are in the package (`triplespec_reducer/data/logo`); a reinstall
restores them.

### "No lamps-on flats were found"

The flats are recognised by `OBJECT = DFLAT_ON` / `DFLAT_OFF` and `OBSTYPE = Dflat`. If your
night names them differently, expand **Every frame and what it was taken to be** on the Night
step and set their kind by hand, then **Save the plan**.

### A target is split in two, or two targets are merged

Targets are grouped by `OBJECT` and pointing; a pointing more than the slit length away
becomes a sky (object–sky nodding) or a second target (`_pos2`). Rename a target, or change
its type, in the Night table.

### A science target has no standard, or the wrong one

The standard is the A0 V star nearest in time among those within 0.1 airmass. Choose another
in the table's **Standard** column. A note warns when none is within 0.1 airmass; the
correction is then less accurate in the deep absorption bands.

### B and V are empty

SIMBAD was not reached (no network) or does not know the star. Type the magnitudes in the
Night table (or on the Telluric step). They set the absolute flux scale.

### A wavelength solution has few lines or a large RMS

Usually the automatic offset failed (very short exposures, few sky lines). On the Calibrate
step, untick *Find the offset automatically*, type the shift in pixels of one recognisable
line between this night's spectrum and the template (order 4), and calibrate again. The arc
lamps are an alternative: set the target's **Wavecal from** to `arcs`.

### The profile shows no clear peaks

The target is faint or the pair is bad. Preview another pair; use **Guess** with typed
positions; or import the trace of a bright star (**Import**). For an extended source, type
the aperture position.

### The kernel deviation is several per cent

Change the continuum windows (they must avoid the line wings and telluric features) or the
degree (3 usually beats 1), and build the kernel again. A standard much fainter than usual,
or not A0 V, also shows here.

### Bumps or dips at the hydrogen lines of the corrected spectrum

The line scales need adjusting: on the Telluric step, select the order, and lower the scale
of a line where the telluric spectrum shows a bump, raise it where it dips — or tick
*Estimate the scale by dragging on the plot* and drag across the line.

### Residual absorption bands in the corrected spectrum

The standard and the target were observed at different airmasses or far apart in time; or
the shift was not found well. Check the airmass difference (Night step), and the shifts
(Telluric step): type a region of strong, unsaturated absorption for an order whose shift
looks wrong.

### A step between two orders in the merged spectrum

Set that order's scale to `auto`, drag over a clean stretch of the overlap, and merge again;
or trim a noisy order end.

### Where are the logs?

Every step's settings and results are in `proc/<source>_log.json`, shown on the Products
page. A failure dialog offers to open the output folder; the command line prints the step
and the reason.
