# Every option

`triplespec_reducer.pipeline.Options` holds every setting; the interface edits the same
object and `tsr reduce --options file.json` reads one. The IDL column gives the panel field or
instrument-file keyword it corresponds to.

## Extraction (`xspextool`)

| option | default | IDL | meaning |
|---|---|---|---|
| `mode` | `A-B` | Reduction mode | `A-B`, `A`, `A-Sky` |
| `combine_images` | off | Combine Images | average the pairs into one image, then extract once |
| `combimg_stat`, `combimg_thresh` | Robust Mean (Std Error), 8 | COMBSTAT, COMBTHRESH | its statistic |
| `abba2abab` | on | ABBA to ABAB | swap every second pair so all are A−B |
| `combimg_subtract_sky` | off | Subtract sky | column-by-column sky level removed first |
| `flatfield`, `errorprop` | on, on | FLATFIELD, ERRORPROPAGATION | divide by the flat; propagate errors |
| `psnaps` / `psnaps_single` | 2 / 1 | PSNAPS | apertures per point source in A−B / in A |
| `findmode` | auto | Auto / Guess / Fix | how apertures are found |
| `apradius` | 1.0″ | PSAPRAD | aperture radius |
| `psbgstart`, `psbgwidth`, `psbgdeg` | 2.5″, 5.0″, 1 | PSBGSTART, PSBGWIDTH, PSBGDEG | point-source background |
| `tracestep`, `tracesumap`, `tracedeg` | 5, 5, 1 | TRACESTEP, TRACESUMAP, TRACEDEG | the trace |
| `tracewinthresh`, `tracesigthresh` | 3, 0.1 | TRACEWINTHRESH, TRACESIGTHRESH | acceptance of a trace point |
| `signs` | from the profile | Override Signs | e.g. `[1, -1]` |
| `store_trace`, `import_trace` | off, none | Store / Import Trace | |
| `xs_positions`, `xs_radii` | profile peak, 1.0″ | Ap Positions, Ap Radii | extended-source apertures (several for resolved spectra) |
| `xs_bg`, `xs_bgdeg` | 0−2,13−15, 0 | XSBG, XSBGDEG | extended-source background |
| `xs_background` | **off** | — | subtract it (IDL never does; see [Behaviour of the IDL tool](../validation/idl_differences.md)) |
| `ybuffer`, `atmosthresh` | 1, 0.4 | YBUFFER, ATMOSTHRESH | profile construction |
| `wavecal_offset` | automatic | Auto X-Correlate | the pixel offset of the home order |

Optimal extraction, bad-pixel fixing and the linearity correction are switched off and
unavailable in the TS4_SOAR instrument file ("0 0"), and are absent here too.

## Combining (`xcombspec`)

| option | default | IDL | meaning |
|---|---|---|---|
| `scaleorder` | 5 | Scale Order | the order the scale factors come from |
| `statistic` | Robust Weighted Mean | Combination Statistic | one of nine |
| `combthresh` | 8 | Robust Threshold | |
| `correct_shapes` | off | Correct Spectral Shape | match the spectra's continuum shapes first |
| `combine_apertures` | off | Combine Apertures | all apertures into one |

## Telluric correction (`xtellcor`)

| option | default | IDL | meaning |
|---|---|---|---|
| `tcmethod` | Deconvolution | Method | or IP |
| `kernel_order` | 6 | Order | the order with the H I line |
| `normreg`, `normdegree` | see [How the night is organized](automation.md); **3** | 'n' regions; Fit Degree (1 in the panel) | the continuum around the line. The tutorial uses 3, and on this night it halves the kernel deviation; `Options.idl_defaults()` starts at 1 as the panel does |
| `linereg` | 0.997–1.013 µm | 's' region | the line |
| `tension` | 10 | Tension | the line-scale spline |
| `vshift` | measured | Vega dv | override Vega's velocity |
| `shift_regions` | from ATRAN | 's' region | per order |
| `units` | ergs s⁻¹ cm⁻² Å⁻¹ | Units | or erg/s/cm²/Hz, W/m²/µm, W/m²/Hz, Jy |
| `propstdunc` | on | Propagate Std Uncertainty | |

## Merging (`xmergeorders`) and cleaning (`xcleanspec`)

| option | default | IDL | meaning |
|---|---|---|---|
| `anchor` | 3 | Anchor Order | |
| `merge_scales` | 1 for every order | Scale | a number or `auto` (with a region) |
| `merge_trims` | none | Cut | (`Green`/`Blue`, `Left`/`Right`, µm) |
| `merge_weighted`, `merge_sum` | on, off | Weighted / Sum | |
| `clean_method` | Savitzky-Golay | Smoothing | or Gaussian |
| `sgwin`, `sgdeg`, `fwhm` | 1.5 slit widths, 2, 1 slit width | SG Window, SG Degree, FWHM | note: the SG default truncates to 3 pixels, which with degree 2 does not smooth — as in IDL |
