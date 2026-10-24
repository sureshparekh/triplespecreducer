# 3. Extraction

:::{admonition} Why
:class: note
The sky is removed by subtracting the B frame from the A frame of every pair. What is
left is the target — twice for a point source (once positive, once negative), once for an
extended source. Extraction finds the target on the slit and adds up its light at every
wavelength, in every order, giving one spectrum per image of the target.
:::

The standard is extracted first, because it is a point source and shows every part of the
step; then the galaxy.

## A point source: HD 225200

1. Go to **3 Extract** and choose **HD225200_1** as the source.
2. Leave the mode at **A-B pairs**, and the point-source settings at their defaults
    (2 apertures, found automatically, radius 1″, background from 2.5″ to 7.5″, traced).
3. Press **Preview the first pair**.

```{shots} 03-extract-preview-std
:alt: The spatial profiles of the standard

The spatial profiles of the standard
```

The stage shows the **spatial profile** of every order: the light along the slit, averaged
over wavelength. The two peaks are the star at A (positive, 9.9″) and at B (negative, 28.1″).
The **green bands** are the apertures (±1″ around each peak) and the **red bands** the
regions where the background is measured.

:::{admonition} Check
:class: tip
* Two clean peaks of opposite sign in every order, each inside its green band.
* No red band overlaps a peak. If the star is extended or the seeing poor, make the
  **Radius** larger and push **Background start** out.
:::

Press **Extract every pair**. Each pair writes two spectra (`HD225200_1_spec0116.fits` for
A, `…0117.fits` for B); click one in the list to plot it.

## An extended source: NGC 134

Choose **NGC134**. It is *extended*, nodded to blank sky, so A − B has the galaxy once and the
sky cancels. Press **Preview the first pair**:

```{shots} 03-extract-preview-sci
:alt: The spatial profile of NGC 134

The spatial profile of NGC 134
```

The galaxy's nucleus peaks at 19.2″ along the slit. With **Positions** blank the aperture is
centred on that peak in every pair; type positions (comma-separated, in arcsec) to extract
other parts of the galaxy — several positions give spatially resolved spectra.

Press **Extract every pair**, then click a spectrum in the list:

```{shots} 03-extract-spectrum
:alt: An extracted spectrum of NGC 134

An extracted spectrum of NGC 134
```

Each colour is one order, in counts per second (DN s⁻¹). The shape is not the galaxy's yet:
it is multiplied by the atmosphere's transmission (the deep gaps at 1.4 and 1.9 µm) and the
instrument's response. Those are removed in step 5.

:::{admonition} The extended-source background
:class: note
The IDL tool draws background regions for extended sources but never subtracts them (see
[Behaviour of the IDL tool](../validation/idl_differences.md)). The default here does the same; **Subtract the background** applies it. For object–sky
nodding the sky is already gone, so it rarely matters.
:::

**Next:** [4. Combining](4-combine.md) · *The science:* [extraction](../science/extraction.md)
