# The data you need

A TripleSpec night produces a few hundred FITS files. All of them go in one folder; the
program reads their headers and works out what each is. This page says what each kind of
frame is, why it was taken, and how the program recognises it. The examples are from the
night of 2025 November 2, the night the [tutorial](../tutorial/index.md) reduces: 235 frames.

:::{admonition} Only read, never changed
:class: note
The raw folder is only ever read. Nothing in it is written, renamed or moved, so you can
point the program at your original data without making a copy.
:::

## What a frame looks like

```{themed} science/raw-orders
:alt: A raw frame with the five orders outlined

A raw frame with the five orders outlined
```

TripleSpec is a *cross-dispersed* spectrograph: a grating disperses the light along one axis
and prisms separate the grating's orders along the other, so one 2048 × 1024 frame holds five
spectra side by side — orders 3 to 7, from the K band at the bottom to the z and J bands at the
top. The program works on the frames turned so that wavelength runs along the columns, as
Spextool does. Each order is a curved strip 108 pixels (38″) tall: the image of the slit.
In the frame above, the bright streaks inside each strip are the star, twice — once at each
of the two positions it was nodded to.

## The frames, kind by kind

| kind | header | in the test night | what it is for |
|---|---|---|---|
| **Lamps-on flats** | `OBJECT = DFLAT_ON`, `OBSTYPE = Dflat` | 50 frames, 2 s | the flat field: the detector's pixel-to-pixel sensitivity, and where the orders are |
| **Lamps-off flats** | `OBJECT = DFLAT_OFF` | 52 frames | the thermal glow of the telescope and instrument, subtracted from the lamps-on frames |
| **Arcs** | `OBSTYPE = Calibration`, `OBJECT = arc` | 12 frames (CuHeAr lamp) | wavelength calibration from lamp lines; usually the sky's own lines are used instead |
| **Science targets** | `OBSTYPE = Object` | NGC 134 (11 frames, 180 s), three more | what you want to measure |
| **Telluric standards** | `OBSTYPE = Object`, names HD…, short exposures | HD 225200 (four visits), HD 1777 | A0 V stars near each target, to remove the atmosphere and calibrate the flux |

Anything else — focus frames, tests, a copy of a frame under another name, a difference
image made at the telescope — is set aside with a note saying why.

### Flats

A quartz lamp illuminates the slit evenly. In the lamps-on frames every order is a bright
band; the edges of those bands tell the program where each order is, and the small
variations inside them are the sensitivity differences of the detector's pixels. The warm
telescope and instrument also glow in the infrared, adding a background to every flat; the
lamps-off frames measure that glow so it can be subtracted. → [The flat field](../science/flat.md)

### Arcs

A copper–helium–argon lamp gives emission lines at known wavelengths. TripleSpec's arcs are
usually taken at the start of the night, hours before the targets, and the instrument flexes
in the meantime, so the tutorial — and this program, by default — calibrates every target
on the night sky's own OH emission lines, recorded in the target's frames at the same moment
as the target. The arcs remain available as an option. → [Wavelength calibration](../science/wavecal.md)

### Science targets: nodding

The near-infrared sky is bright — thousands of times brighter than a faint galaxy — and it
varies by the minute. It is removed by **nodding**: the telescope moves the target between
two places, A and B, and one frame is subtracted from the next.

* **Point sources** (stars, compact nuclei) are nodded *along the slit*, in an **ABBA**
  pattern: both A and B contain the target, 15–20″ apart on the slit, and A − B contains it
  twice — once positive, once negative — with the sky gone.
* **Extended sources** (a galaxy filling the slit) are nodded *off the target* to blank sky:
  **object–sky**. In NGC 134's sequence, frames 120, 122, 123, … are on the galaxy and 121,
  127, 130 are 50″ away on sky.

The program reads the pointing (RA, Dec) of every frame, finds the pattern, and pairs the
frames: A with B, or each object frame with the sky frame nearest in time. A second pointing
on the same object (along a galaxy's disc, for example) is kept as a separate target
(`NGC134_pos2`), so two places are never averaged together.

### Telluric standards

Earth's atmosphere absorbs much of the near infrared, and how much changes with the
airmass and the water vapour. An **A0 V star**, observed just before or after each target at
a similar airmass, measures it: the star's true spectrum is known (it is like Vega's), so
whatever differs is the atmosphere and the instrument. Its **B and V magnitudes** set the
flux scale of the final spectrum; the program looks them up on SIMBAD (or you type them).
Each science target is given the standard nearest to it in time among those within 0.1 in
airmass. → [Telluric correction](../science/telluric.md)

## What the header must contain

The program uses only standard TripleSpec keywords: `OBJECT`, `OBSTYPE`, `EXPTIME`,
`COADDS`, `FSAMPLE` (the number of Fowler samples, for the noise model), `RA`, `DEC`,
`AIRMASS`, `DATE-OBS` and `UT`. Files must be named as the instrument writes them
(`SPEC_<prefix>_<number>.fits`).

## The night in numbers (the tutorial's)

| target | role | type | pattern | pairs | airmass | standard |
|---|---|---|---|---|---|---|
| HD225200_1 | standard | point | ABBA | 2 | 1.09 | — |
| **NGC134** | science | extended | object–sky | 7 | 1.08 | HD225200_1 |
| NGC134_pos2 | science | extended | object–sky | 1 | 1.08 | HD225200_1 |
| HD1777 | standard | point | ABBA | 2 | 1.26 | — |
| kks_03-gc1 (two objects) | science | point | ABBA | 14 | 1.45 | HD225200_3 |
| HD225200_2 | standard | point | ABBA | 2 | 1.05 | — |
| NGC 253 gc 73 | science | point | ABBA | 16 | 1.11 | HD225200_2 |
| HD225200_3 | standard | point | ABBA | 2 | 1.44 | — |
| NGC 253 gc 213 | science | point | ABBA | 12 | 1.56 | HD225200_3 |
| HD225200_4 | standard | point | ABBA | 2 | 2.53 | — |

Every one of these decisions is shown in the Night table and can be changed there; see
[How the night is organized](../gui/automation.md) for the rules.
