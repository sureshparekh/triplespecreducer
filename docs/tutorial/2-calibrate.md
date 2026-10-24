# 2. Calibrations

:::{admonition} Why
:class: note
Two calibrations come before any spectrum. The **flat field** corrects the different
sensitivities of the detector's pixels and tells the program where the five orders lie.
The **wavelength calibration** gives every pixel its wavelength and its position along
the slit: TripleSpec's orders are curved and the image of the slit is tilted, so neither
is a simple function of row and column.
:::

## The flat field

1. Go to **2 Calibrate**. **Divide the images by the flat** stays ticked.
2. Press **Build the flat**. It takes about half a minute: the 52 lamps-off frames are
    averaged and subtracted from the 50 lamps-on frames, the order edges are traced, and each
    order is divided by a smooth model of itself.
3. Press **Show it**.

```{shots} 02-flat
:alt: The normalised flat with its order edges

The normalised flat with its order edges
```

The stage shows the **normalised flat**: close to 1 inside the orders, with the pixel-to-pixel
pattern of the detector, and the fitted edges of every order drawn over it. The title gives
the RMS of the smooth model in each order — a few per cent.

:::{admonition} Check
:class: tip
The edges should follow the orders all the way along, with no order missing or cut short.
Order 7 (the top one) is short: it falls off the detector.
:::

## The wavelength solution

The table lists one solution per target that needs one — each science target is calibrated
on its own sky lines; each standard uses the solution of the target nearest in time.

1. Leave **Find the offset automatically** ticked.
2. Press **Calibrate all** (or select a row and **Calibrate the selected**). About a minute per
    solution.
3. Select **NGC134** in the table.

```{shots} 02-wavecal
:alt: The wavelength solution of NGC 134

The wavelength solution of NGC 134
```

What happened: the median of NGC 134's 11 frames was divided by the flat; a spectrum of the
sky was extracted down the middle of every order; it was cross-correlated with a stored
template (offset −0.66 pixels this night); the OH lines were found — 31, 64, 63, 76 and 32 in
orders 3 to 7; one polynomial in column and order was fitted to all 266 lines, with an RMS
of **0.34 Å**; and every line was traced along the slit to map the tilt. The stage shows the
result: the wavelength of every pixel.

:::{admonition} Check
:class: tip
* **RMS** a fraction of an ångström (0.2–0.5 Å is typical), and tens of lines per order.
* The wavelength image changes smoothly along each order, from blue (short) to yellow (long).
* If a solution fails or has few lines, untick *Find the offset automatically*, type the
  offset of one line between this night and the template, and calibrate again.
:::

**Next:** [3. Extraction](3-extract.md) · *The science:* [the flat field](../science/flat.md),
[wavelength calibration](../science/wavecal.md)
