# Every option

`triplespec_reducer.pipeline.Options` holds every setting; the interface edits the same
object and `tsr reduce --options file.json` reads one. The IDL column gives the panel field or
instrument-file keyword it corresponds to. Defaults are the IDL tool's, with the TripleSpec
instrument file (`ts4_soar.dat`), unless the column says otherwise; options marked
**not in IDL** do something the IDL tool cannot, and their results have no IDL counterpart.

## Calibrations (`xspextool`, TS4_SOAR)

| option | default | IDL | meaning |
|---|---|---|---|
| `Night.flat_name` | `flat` | Output Flat Name | the flat's file in `cal/`; another name keeps several flats |
| `flatfield` | on | FLATFIELD | divide the images by the flat |
| `wavecal_frames` | the target's own frames (or the arcs) | Wavecal Images | per solution: any list of frames |
| `wavecal_arctype` | as planned | Sky / CuHeAr Arc | per solution: the line list |
| `wavecal_offset` | automatic | Auto X-Correlate | the home order's pixel offset; *Measure on a line* gives it as XGETOFFSET does |
| `short_lines` | `abort` | — | a line traced on fewer than two points stops the calibration, as in IDL; `skip` leaves it out (**not in IDL**) |

Each calibration writes `<name>_Wavecal.pdf`, `<name>_Distortion.pdf` and
`<name>_Distortion.txt` (the lines with rejected points) beside it, as the IDL tool writes its
`.eps`/`.ps` plots. Raw frames are found with the instrument file's `*.fits*`, so `.fits.gz`
and `.fits.fz` are read too. The bad-pixel mask is the instrument file's.

## Extraction (`xspextool`)

| option | default | IDL | meaning |
|---|---|---|---|
| `mode` | `A-B` | Reduction mode | `A-B`, `A`, `A-Sky` |
| `combine_images` | off | Combine Images | average the pairs into one image, then extract once |
| `combimg_stat`, `combimg_thresh` | Robust Mean (Std Error), 8 | COMBSTAT, COMBTHRESH | its statistic |
| `abba2abab` | on, for ABBA sequences only | ABBA to ABAB | swap every second pair so all are A−B; the plan decides from each sequence's own pairs, and an ABAB or irregular sequence is never swapped |
| `combimg_subtract_sky` | off | Subtract sky | the sky level of each column removed first |
| `flatfield`, `errorprop` | on, on | FLATFIELD, ERRORPROPAGATION | divide by the flat; propagate errors |
| `psnaps` / `psnaps_single` | 2 / 2 | PSNAPS | apertures per point source in A−B / in A and A−Sky (1–4) |
| `findmode`, `guess` | auto, none | Auto / Guess / Fix | how apertures are found; one typed list for every order, its length sets the number of apertures; Fix also stops the tracing |
| `trace` | on | (Use Ap Positions off) | trace the apertures |
| `apradius` | 1.0″ | PSAPRAD | aperture radius |
| `psbgsub` | on | PSBGSUB | subtract the background, and extract optimally, as the IDL tool does |
| `psfradius` | 2.2″ | PSPSFRAD | the radius of the profile used by optimal extraction |
| `medprofile` | off | Median Profile | one profile for the whole order |
| `optimal` | follows `psbgsub` | — | `False` with the background on sums the aperture (**not in IDL** for TripleSpec) |
| `psbgstart`, `psbgwidth`, `psbgdeg` | 2.5″, 5.0″, 1 | PSBGSTART, PSBGWIDTH, PSBGDEG | point-source background |
| `ps_orders`, `xs_orders` | all | PS Orders, XS Orders | extract only these orders |
| `tracestep`, `tracesumap`, `tracedeg` | 5, 5, 1 | TRACESTEP, TRACESUMAP, TRACEDEG | the trace |
| `tracewinthresh`, `tracesigthresh` | 3, 0 | TRACEWINTHRESH, TRACESIGTHRESH | acceptance of a trace point (the file says 0.1, which the IDL tool reads as an integer: 0) |
| `signs` | from the profile | Override Signs | e.g. `[1, -1]` |
| `store_trace`, `import_trace` | off, none | Store / Import Trace | an imported trace also fixes the orders extracted |
| `xs_positions`, `xs_radii` | profile peak, 1.0″ | Ap Positions, Ap Radii | extended-source apertures (several for resolved spectra) |
| `xs_bg`, `xs_bgdeg` | 0−2,13−15, 0 | XSBG, XSBGDEG | extended-source background |
| `xsbgsub` | on | XSBGSUB | used only without a wavecal file |
| `xs_background` | **off** | — | subtract the extended-source background with a wavecal file (IDL never does; see [Behaviour of the IDL tool](../validation/idl_differences.md)) |
| `no_wavecal`, `bdpixthresh` | off, 7 | blank Full Wavecal Name, Bad Pix Thresh | extract in pixels, with the bad-pixel mask |
| `ybuffer`, `atmosthresh` | 1, 0.4 | YBUFFER, ATMOSTHRESH | profile construction |

Bad-pixel fixing and the linearity correction are switched off and unavailable in the
TS4_SOAR instrument file ("0 0"); the frames are read as the TS4 module reads them.

## Combining (`xcombspec`)

| option | default | IDL | meaning |
|---|---|---|---|
| `scale_spectra` | on, as in the tutorial | Scale Spectra | bring the spectra to a common level |
| `scale_mode`, `scale_spectrum` | Median | Scale to: Median / Spectrum | the reference: the median spectrum, or one chosen spectrum |
| `scaleorder`, `scale_range` | 5, the middle 80% | Scale Order, 's' region | where the factors are measured |
| `prune` | none | Prune (XMC_MASKSPEC) | per source: spectra left out of one order, wavelength ranges masked |
| `statistic` | Robust Weighted Mean | Combination Statistic | one of nine |
| `combthresh` | 8 | Robust Threshold | |
| `correct_shapes` | off | Correct Spectral Shape | match the spectra's continuum shapes first (not done for apertures that are not combined, as in IDL) |
| `combine_apertures` | off | Combine Apertures | all apertures into one |
| `combine_basic` | off | `/BASIC` | the shorter header |

## Telluric correction (`xtellcor`, `xtellcor_basic`, `xtellcor_finish`)

| option | default | IDL | meaning |
|---|---|---|---|
| `tcmethod` | Deconvolution | Method | IP is offered by the IDL tool only for SpeX (**not in IDL** for TripleSpec) |
| `kernel_order` | 6 | Order | the order with the H I line |
| `normreg`, `normdegree` | see [How the night is organized](automation.md); **3** | 'n' regions; Fit Degree (1 in the panel; 1–15) | the continuum around the line. The tutorial uses 3, and on this night it halves the kernel deviation; `Options.idl_defaults()` starts at 1 as the panel does |
| `linereg` | 0.997–1.013 µm | 's' region | the line |
| `tension` | 10 | — | the line-scale spline; fixed at 10 in IDL, which has no control for it |
| `vshift` | measured | Vega dv | override Vega's velocity |
| `shift_regions`, `find_shift` | from ATRAN, on | 's' region | per order and aperture |
| `units` | ergs s⁻¹ cm⁻² Å⁻¹ | Units | or erg/s/cm²/Hz, W/m²/µm, W/m²/Hz, Jy |
| `propstdunc` | on | — | propagate the standard's uncertainty; IDL always does |
| `telluric_model` | off | Model | also write `<source>_tc_modVega.fits` |
| `basic_standards` | none | `xtellcor_basic` | per standard that is not an A0 V star: divide by it, and restore the continuum with a blackbody (`temp`, `restore`) scaled to its V |
| `telluric_from` | none | `xtellcor_finish` | per target: apply a telluric spectrum written before |

## Merging (`xmergeorders`) and cleaning (`xcleanspec`)

| option | default | IDL | meaning |
|---|---|---|---|
| `anchor` | 3 | Anchor Order | each order added must be next to those merged so far |
| `merge_scales` | 1 for every order | Scale | a number, or `auto` with a region |
| `merge_regions` | none | 's' region | per order, for `auto`: `(w0, w1)` or `overlap`; found before any cut |
| `merge_trims` | none | Cut | per order, one or two of (`Green`/`Blue`, `Left`/`Right`, µm) |
| `merge_unweighted`, `merge_sum` | off, off | Weighted / Sum | **not in IDL**: the IDL tool stops on an undefined variable with either |
| `clean`, `clean_method` | off, Savitzky-Golay | Smoothing | or Gaussian; once per file opened |
| `sgwin`, `sgdeg`, `fwhm` | from the file's slit width: 1.5 slit widths, 2, 1 slit width | SG Window, SG Degree, FWHM | the SG default truncates to 3 pixels, which with degree 2 does not smooth, as in IDL; the degree must be below the width |

In the interface, the clean step opens any spectrum file (merged, corrected or combined),
any order and aperture; *Fix* and *Remove* act on the flux and its error or on the error
only, *Keep* and *Undo* are the IDL tool's 's' and 'u', and the text export is SPEX2TEXT's.
