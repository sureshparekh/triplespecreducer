# 8. The products

Switch the header to **Products** and type a source name:

```{shots} 08-products
:alt: Every product of NGC 134

Every product of NGC 134
```

Every file made for the source is listed — calibrations, the spectrum of every pair, the
combined, corrected, merged and cleaned spectra — with its header, and the **log** of every
step: the settings used and what came out (the scale factors, the kernel's deviations, the
shifts). Select a file to plot it.

## What is where

```
reduced/2025-11-02/
├── night.json                  the plan
├── cal/
│   ├── flat.fits               the normalised flat
│   └── NGC134_wavecal.fits     NGC 134's wavelength solution
└── proc/
    ├── NGC134_spec0120.fits    the spectrum of pair 120−121 (… one per pair)
    ├── NGC134_comb.fits        combined
    ├── NGC134_tellspec.fits    the telluric correction spectrum
    ├── NGC134_tc.fits          telluric corrected and flux calibrated, five orders
    ├── NGC134_merged.fits      the orders merged: the final spectrum
    ├── NGC134_clean.fits       smoothed / repaired (if you did step 7)
    └── NGC134_log.json         every setting and result
```

**The rule is `<source>_<step>`**, so the source name finds everything. The files are in
Spextool's own format, readable by its IDL tools and by anything that reads Spextool spectra:
see [Files and names](../gui/naming.md).

## Reading the result in Python

```python
from astropy.io import fits
spec = fits.getdata("reduced/2025-11-02/proc/NGC134_merged.fits")   # [4, n]
wave, flux, error, flag = spec            # µm, erg s-1 cm-2 A-1, same, bit flags
```

## Doing it again, all at once

Everything in this tutorial is what **Reduce everything** does for every target, with the
settings you chose. From the command line: `tsr reduce reduced/2025-11-02` — see
[Command line and Python](../gui/cli.md).
