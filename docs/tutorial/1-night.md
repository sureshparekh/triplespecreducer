# 1. The night

:::{admonition} Why
:class: note
A night is a few hundred frames of different kinds. Before anything can be reduced, the
program must know which frames are flats, which are a target's A and B, which standard
star goes with which target, and where each target's wavelength calibration comes from.
In the IDL tool you type all of this as frame numbers; here it is read from the headers
and proposed to you as a plan.
:::

## What you do

1. Open the program: `tsr`. The logo shows while it loads, then the **Night** step opens.
2. **Raw frames** → *Browse*, and choose the folder with the `SPEC_*.fits` files.
3. **Output** is proposed next to it (`reduced/<night>`); change it if you like. Everything
    the reduction writes goes there — the raw folder is only read.
4. Leave **Look up the standards' B, V and type on SIMBAD** ticked if you are online.
5. Press **Read the night**. It takes a few seconds.

```{shots} 01-night
:alt: The Night step with the plan

The Night step with the plan
```

## What you see

The **Targets and their plan** table has one row per target:

| column | meaning | for NGC 134 |
|---|---|---|
| Source | the target's name, from `OBJECT` (made safe for file names) | `NGC134` |
| Role | science or standard | science |
| Source type | point (nodded along the slit) or extended (nodded to sky) | extended |
| Pattern | ABBA or object-sky | object-sky |
| Frames, Pairs | how many frames, and how many A−B pairs they make | 11 frames, 7 pairs (a sky frame serves the object frames either side of it) |
| Airmass | the mean airmass | 1.08 |
| Wavecal from | whose sky lines calibrate it | NGC134 (its own) |
| Standard | the A0 V star that corrects it | HD225200_1 |
| B, V | the standard's magnitudes | 6.40, 6.39 |

Selecting a row shows the target's first A − B image on the stage. Under the table, **notes**
report anything that needed a decision — here, that SIMBAD lists HD 225200 as A1 V rather
than A0 V (close enough; the method assumes A0 V), and that a second pointing on NGC 134 was
kept apart as `NGC134_pos2`.

Open **About this step** for the full account of how the plan is made:

```{shots} 01-night-about
:alt: About this step, opened

About this step, opened
```

:::{admonition} Check before going on
:class: tip
* **Every target has the right role and type.** A galaxy filling the slit must be
  *extended*; a star or compact source *point*. Change a cell by clicking it.
* **Every science target has a standard**, ideally within 0.1 in airmass (NGC 134: 1.079
  against 1.088). A note warns when none is that close.
* **B and V** of each standard are filled. They set the absolute flux scale; type them if
  SIMBAD was not reached.
* Expand **Every frame and what it was taken to be** if a frame seems missing: each frame
  is listed with its kind, and a note for anything set aside.

Press **Save the plan** after any change. The plan is `night.json` in the output folder.
:::

**Next:** [2. Calibrations](2-calibrate.md) · *The science:* [the data you need](../start/data.md),
[how the night is organized](../gui/automation.md)
