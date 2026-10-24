# Quick start

## The whole night, in the interface

1. **Open it.** `tsr` (or `triplespec-reducer`). The logo shows for a few seconds while the
   reduction modules load; then the window opens on the **Night** step.
2. **Night.** Choose the folder of raw frames (`SPEC_*.fits`). The output folder is proposed
   next to it (`reduced/<night>`); change it if you like. Press **Read the night**.

   The table shows every target, what it is taken to be and how its frames pair:

   ```{shots} 01-night
:alt: The night

The night
```

   Check it. In particular:
   * **Role** — `science` or `standard`. Standards are recognized by their names (HD, HIP, HR, …)
     or by exposures of 30 s or less.
   * **Source type** — `point` (an along-slit ABBA nod) or `extended` (object and sky
     pointings). It follows the nod pattern; change it if you disagree.
   * **Standard** — the A0 V star used for each target: the one closest in time among those
     within 0.1 in airmass (xtellcor's own tolerance), else the closest in airmass.
   * **B, V** — the standard's magnitudes, looked up on SIMBAD. Type them if SIMBAD was not
     reached. A note warns when SIMBAD's spectral type is not A0 V.
   * The notes under the table say what was decided and why: a frame left out because it
     repeats its neighbour's position, a second pointing kept as its own target, a standard
     more than 0.1 airmass away.

3. **Reduce everything.** The button at the top right runs every step for every science
   target and its standard, with the defaults. It takes several minutes per target.

   If one target fails (a missing frame, a wavelength solution that cannot be found), the
   others carry on; a summary at the end lists what did not finish and why.

That is all a standard reduction needs. The result for each target is
`proc/<source>_merged.fits` — telluric corrected, flux calibrated, all five orders joined.

## The same, step by step

The [tutorial](../tutorial/index.md) walks through every step with screenshots. In short:

Each step in the bar along the top can be opened, adjusted and run on its own:

| step | press | look at |
|---|---|---|
| Calibrate | **Build the flat**, then **Calibrate all** | the normalized flat and its order edges; the RMS of each wavelength solution (≈ 0.3 Å is typical) |
| Extract | **Preview the first pair**, then **Extract every pair** | the spatial profiles with apertures (green) and background (red) |
| Combine | **Combine** | every scaled spectrum (grey) under their combination |
| Telluric | **Load the two spectra**, **Build the kernel**, adjust the H I lines, **Find every shift**, **Write the corrected spectrum** | the kernel residual (inside ±1 % is good), the telluric spectrum around each H I line |
| Merge | **Merge** | the overlaps between orders |
| Clean & save | **Smooth and save** | the final spectrum and its S/N |

## The whole night, from the command line

```bash
tsr organize data/raw/2025-11-02 reduced/2025-11-02     # writes reduced/2025-11-02/night.json
# (edit night.json if a decision is wrong, or open it in the interface)
tsr reduce reduced/2025-11-02                           # every science target
tsr reduce reduced/2025-11-02 NGC134                    # just one
tsr products reduced/2025-11-02 NGC134                  # every file made for it
```
