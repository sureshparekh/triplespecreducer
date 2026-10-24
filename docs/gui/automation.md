# How the night is organized

The IDL tool asks for frame numbers, prefixes, which star goes with which target and where the
apertures are. Here those inputs are proposed from the headers and the data. The rules are
below; every result is in `night.json` and the Night table, and every one can be changed.
None of them alters a computation — they only choose what the IDL user would type.

## Frames

| kind | rule |
|---|---|
| lamps-on flat | `OBSTYPE` Dflat and `OBJECT` DFLAT_ON, with the most common file prefix (a lone lamps-on frame under another prefix, `SPEC_141Tau_0001` in the test night, is set aside) |
| lamps-off flat | `OBJECT` DFLAT_OFF |
| arc | `OBSTYPE` Calibration and `OBJECT` arc |
| object | `OBSTYPE` Object |
| ignored | anything else; and any file whose header `FILENAME` is not its own name — a renamed copy or a difference image made at the telescope (`t20-21.fits`) |

A quick-look check compares the lamps-on and lamps-off levels, as the tutorial does by eye.

## Targets, nods and pairs

Consecutive object frames with one `OBJECT` form a target. Their pointings (`RA`, `DEC`) are
grouped (2″ tolerance):

* **All within the slit length (38″) of each other**: an along-slit nod, **ABBA** — a point
  source. Frames are paired in time order; a frame that repeats its successor's position is
  left out with a note, as the tutorial leaves out frame 432.
* **Some farther than the slit**: those are **sky** — an extended source. Each object frame is
  paired with the sky frame nearest in time, so a sky frame can serve two object frames
  (433−434 and 435−434 in the tutorial).
* **A second object pointing** (an offset along a galaxy): its frames become their own target,
  `<name>_pos2`, so two places are not averaged together.

**Standards** are targets whose names start HD, HIP, HR, BD, SAO, TYC, … or whose exposures
are 30 s or shorter. Their B and V magnitudes and spectral type come from SIMBAD (cached in
`~/.triplespec_reducer/simbad.json`); a note warns when the type is not A0 V.

## Calibrations

* **Wavelength**: a science target of exposures ≥ 60 s is calibrated on its own OH lines; a
  second pointing uses its parent's; a standard uses the solution of the science target
  nearest in time; with no long exposures at all, the arcs.
* **Telluric standard**: among the standards within 0.1 airmass of the target (xtellcor
  warns beyond that), the nearest in time; if none is that close, the nearest in airmass,
  with a note.

## Inputs chosen during the reduction

| IDL asks for | chosen here | recorded in |
|---|---|---|
| extended-source aperture position | the profile peak (one peak found automatically), per pair, to 0.1″ | `<source>_log.json` → extract |
| continuum windows and line region for the kernel | 0.990–0.997 and 1.013–1.020 µm around Paδ 0.997–1.013 µm, degree 3 | → telluric |
| the region for each order's shift | the 0.03 µm stretch with the most ATRAN structure that is nowhere below 20 % transmission; the next one if the best shift falls at the edge of the ±1.5 pixel search | → telluric.shifts |
| merge scale factors | 1, as in the tutorial | → merge |
| combine-scale order | 5 | → combine |

`<source>_log.json` also flags what deserves a look — a spectrum that needed a scale factor
far from the others, an order whose shift could not be found.
