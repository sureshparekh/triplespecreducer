# Command line and Python

## Commands

```text
tsr                                   the interface
tsr gui RAW                           the interface, with the raw folder filled in
tsr organize RAW OUT [--no-simbad]    read the night, write OUT/night.json, print the plan
tsr reduce OUT [SOURCE ...]           run the steps (all science targets by default)
      --steps flat,wavecal,extract,combine,telluric,merge[,clean]
      --force                         redo steps whose products exist
      --idl-defaults                  every setting as the IDL panels start
      --xs-background                 subtract extended-source backgrounds
      --combine-images                average the pairs into one image first
      --options FILE.json             any Options fields, e.g. {"apradius": 1.2}
tsr products OUT SOURCE               every file made for SOURCE
```

A step is skipped when its product exists, so `tsr reduce` can be stopped and started again;
`--force` redoes them.

A target that fails at some step does not stop the others: the run carries on, prints a
summary of what did not finish and why, and `tsr reduce` then exits with status 1 (0 when
everything finished). `Reducer.failures` holds the same in Python.

## In Python

```python
from triplespec_reducer.night import organize, Night
from triplespec_reducer.pipeline import Reducer, Options

night = organize("data/raw/2025-11-02", "reduced/2025-11-02")
night.seq("NGC134").standard = "HD225200_1"          # change any decision
night.save()

r = Reducer(night, Options(apradius=1.2), progress=print)
r.flat()
r.wavecal("NGC134")
r.extract("NGC134")
r.combine("NGC134")
r.telluric("NGC134")
r.merge("NGC134")
```

Every IDL routine is also available on its own, under the name of the module that holds it:
`flat.ts4_findorders`, `flat.normspecflat`, `wavecal.xdwavecal1d`, `wavecal.findlines2d`,
`extract.mkspatprof2d`, `extract.findpeaks`, `extract.tracespec2d`, `extract.extspec2d`,
`combine.combine`, `telluric.vegaconv`, `telluric.mktellspec`, `telluric.findshift`,
`telluric.apply_telluric` (xtellcor_finish), `merge.mergespec`, `clean.savitzky_golay`, … —
see [Python API](../reference/api.md).

## Speed

Wavelength solutions and extractions are independent of one another, and `tsr reduce` (and
**Reduce everything**) runs them in parallel processes — by default up to 16, two fewer than
the machine has cores; `{"workers": N}` in an options file sets it, 1 runs everything in one
process. The products do not depend on the number: the test night reduced in one process and
on 16 gives the same 132 files, byte for byte.

The resampling grid of each order — a triangulation of ~150 000 pixels — depends only on the
wavelength solution, so it is made once per solution and reused for every frame. On the
reference machine a pair takes about 12 s, a wavelength solution a little over a minute, and
the whole test night (56 pairs, 4 wavelength solutions) under five minutes on 16 processes.
