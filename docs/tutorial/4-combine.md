# 4. Combining

:::{admonition} Why
:class: note
Every pair gave its own spectrum. Averaging them increases the signal-to-noise ratio — by
√N for N spectra — and rejects anything that appears in only one (a cosmic ray, a bad
pixel). Their levels differ a little, because seeing and guiding change during the
sequence, so they are scaled to one another first.
:::

## What you do

1. Go to **4 Combine** and choose **HD225200_1**.
2. Leave **Scale on order** at 5, **Statistic** at *Robust Weighted Mean* and **Threshold**
    at 8.
3. Every spectrum of the target is listed, ticked. Press **Combine**.

```{shots} 04-combine-std
:alt: The standard's spectra combined

The standard's spectra combined
```

The stage shows every order: each scaled spectrum in grey and the combination on top. Now
choose **NGC134** and **Combine** again:

```{shots} 04-combine-sci
:alt: NGC 134's spectra combined

NGC 134's spectra combined
```

## What happened

1. Every spectrum was interpolated onto the wavelengths of the first.
2. In order 5, the middle 80 % of each spectrum was compared with the median of all of them;
    the ratio's median is that spectrum's **scale factor**, applied to all its orders.
3. At every wavelength the scaled values were combined: values more than 8σ from the mean
    were rejected, and the rest averaged with weights 1/σ². The error of the result is that
    of the weighted mean.

:::{admonition} Check
:class: tip
* The scale factors are close to 1 (within a few per cent). One much lower than the rest
  marks a spectrum taken through cloud or with the target off the slit: **untick it** and
  combine again.
* The grey spectra agree with one another apart from noise.
:::

**Next:** [5. Telluric correction](5-telluric.md) · *The science:* [combining spectra](../science/combine.md)
