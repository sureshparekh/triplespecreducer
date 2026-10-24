# Installation

TripleSpec Reducer needs Python 3.10 or newer. It does **not** need IDL.

```bash
git clone <repository> triplespec-reducer      # or copy the folder
cd triplespec-reducer
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

This installs the package and two commands, `triplespec-reducer` and its short form `tsr`.
Everything the reduction reads ships with the package: the TS4_SOAR instrument files
(`XD.dat`, `ts4_soar.dat`, the OH and Ar line lists, the line-calibration templates), the
ATRAN atmospheric models and the Vega models of Spextool 4.1, and the fonts and icons of the
interface.

Dependencies: numpy, scipy, astropy, matplotlib, PySide6 (the interface), numba (the
exact-area extraction is compiled), astroquery (B, V and spectral types of the standards
from SIMBAD; optional — without a network they are typed in).

## Checking the installation

```bash
tsr organize /path/to/raw/night /path/to/reduced/night
```

prints the frames it found and the plan it proposes. `tsr` with no arguments opens the
interface: the logo appears while it loads, then the window.

## Re-running the IDL comparison (optional)

The `validation/` folder holds the scripts that compare every step with IDL. They need a
working IDL with Spextool 4.1 and the TS4_SOAR module, and are only for checking the
reduction against the IDL tool; reducing data never needs IDL. See [Agreement with IDL](../validation/index.md).
