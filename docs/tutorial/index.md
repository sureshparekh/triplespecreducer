# Tutorial: one galaxy, start to finish

This tutorial reduces one target from the raw frames to its final spectrum, the way you would
on your own night: **NGC 134**, a spiral galaxy observed on 2025 November 2 with the slit on
its nucleus, corrected with the A0 V star **HD 225200** observed a few minutes earlier.

Every step has its own page. Each page says **why** the step is needed, **what you do** in
the interface, **what you see**, and **what to check** before going on; the science behind
it is linked at the bottom. The screenshots are the program running on these frames.

```{raw} html
<figure class="video">
<video controls preload="metadata" poster="../_static/video/tutorial.jpg">
  <source src="../_static/video/tutorial.mp4" type="video/mp4">
  <track kind="subtitles" srclang="en" label="English" src="../_static/video/tutorial.vtt">
</video>
<figcaption>The whole tutorial as a narrated screen recording (9 minutes). Chapters on the <a href="video.html">videos page</a>.</figcaption>
</figure>
```

::::{cards}
:columns: 3

:::{linkcard} 1. The night
:link: 1-night
:icon: folder

Point the program at the raw frames; check the plan it proposes.
:::

:::{linkcard} 2. Calibrations
:link: 2-calibrate
:icon: grid

The flat field and the wavelength solution.
:::

:::{linkcard} 3. Extraction
:link: 3-extract
:icon: chart

A spectrum from every pair of frames.
:::

:::{linkcard} 4. Combining
:link: 4-combine
:icon: layers

The spectra of a target averaged into one.
:::

:::{linkcard} 5. Telluric correction
:link: 5-telluric
:icon: atom

The atmosphere removed, the flux calibrated.
:::

:::{linkcard} 6. Merging the orders
:link: 6-merge
:icon: lines

Five orders, one spectrum.
:::

:::{linkcard} 7. Cleaning and saving
:link: 7-clean
:icon: stars

The final file.
:::

:::{linkcard} 8. The products
:link: 8-products
:icon: file

What was written, and where.
:::

:::{linkcard} The videos
:link: video
:icon: play

The science and the steps, and the whole reduction as a narrated screen recording.
:::
::::

## Before you start

* The program is [installed](../start/install.md), and opens with `tsr`.
* You have a night's raw frames in one folder. To follow along exactly, use the night of
  2025-11-02; any TripleSpec 4.1 night works the same way.
* About **15 minutes** for the whole tutorial, most of it reading. The computations take a
  few seconds each; the wavelength solution about a minute.

:::{admonition} The fast way
:class: tip
Everything in this tutorial is also what **Reduce everything** (top right) does with the
default settings, for every target at once. The tutorial goes step by step so that you
see — and can check — what happens at each one.
:::

## How the interface is laid out

```{shots} 01-night
:alt: The interface

The interface
```

* **The step bar** across the top: the seven steps, each with a tick when done and a line of
  progress (`8/10 combined`). Click any step to go to it.
* **The left panel**: the current step. It opens with **About this step** — why the step
  exists, what it does in order, and what to check — followed by its settings. Every setting
  has a caption under it; hover over it for the full explanation.
* **The stage**, on the right: whatever the step is working on — a frame, the flat, a
  profile, a spectrum. Its tools (top right) reset, pan and zoom; the mouse wheel zooms.
* **The header**: *Reduce* and *Products* pages, what is running, the light/dark switch,
  and **Reduce everything**.

```{toctree}
:hidden:

1-night
2-calibrate
3-extract
4-combine
5-telluric
6-merge
7-clean
8-products
video
```
