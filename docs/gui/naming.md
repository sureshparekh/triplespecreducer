# Files and names

Everything is written under the output folder of the night; the raw folder is only read.

```
reduced/2025-11-02/
├── night.json                     the plan: every frame, every target, every decision
├── cal/
│   ├── flat.fits                  the normalized flat (Spextool's format)
│   ├── <source>_wavecal.fits      a wavelength solution (Spextool's 2D format)
│   ├── <source>_sky.fits          the combined sky of an A-Sky reduction
│   └── <source>_trace.dat         a stored trace
└── proc/
    ├── <source>_specNNNN.fits     the spectrum extracted from the pair (or frame) starting at NNNN
    ├── <source>_combimg.fits      all pairs averaged into one image (Combined mode)
    ├── <source>_specall.fits      the spectrum extracted from it
    ├── <source>_comb.fits         the combined spectrum               (xcombspec)
    ├── <source>_tellspec.fits     the telluric correction spectrum    (xtellcor)
    ├── <source>_tc.fits           the telluric-corrected, flux-calibrated spectrum
    ├── <source>_merged.fits       the orders merged                   (xmergeorders)
    ├── <source>_clean.fits        smoothed / repaired                 (xcleanspec)
    ├── <source>_clean.txt         the same as a text table
    └── <source>_log.json          every setting and every result of every step
```

**The rule is `<source>_<step>`.** The source name comes from the OBJECT keyword, with
characters a file name cannot hold replaced by `_` (`kks_03-gc1 (two ob…` becomes
`kks_03-gc1_two_ob…`); a target observed several times gets `_1`, `_2`, … in time order, and
a second pointing of the same target `_pos2`. A name can be changed in the Night table.

Because of the rule, giving a source name is enough to find every product:
`tsr products <out> <source>`, the Products page of the interface, or in Python

```python
from triplespec_reducer.night import Night
n = Night.load("reduced/2025-11-02/night.json")
n.products("NGC134")     # {'wavecal': ..., 'spec0120': ..., 'comb': ..., 'tc': ..., 'merged': ...}
```

## Formats

The FITS files are the ones Spextool writes, so they can be opened by its IDL tools (`xvspec`,
`xtellcor_finish`, …) and by anything that reads Spextool spectra:

* a spectrum is a cube `[n_orders × n_apertures, 3 or 4, n_pixels]` — wavelength (µm), flux,
  error and, where present, a flag — with `ORDERS`, `NAPS`, `DISPOnn`, `ITOT`, `AIRMASS` and the
  Spextool keywords in its header and the history of every step — every number written as
  Spextool writes it, so the files are interchangeable with IDL's;
* the merged and cleaned spectra are one order, written as xmergeorders writes them: a 2-D
  array `[4, n_pixels]` (IDL `[npix, 4]`);
* fluxes are DN s⁻¹ until the telluric step and erg s⁻¹ cm⁻² Å⁻¹ after it (or the unit chosen
  there); wavelengths are vacuum microns.
