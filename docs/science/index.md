# The science

These pages explain what every step of the reduction does to the data and why — the physics,
the technique, and the exact algorithm, with the Spextool step each corresponds to. They go
with the [tutorial](../tutorial/index.md), which shows the same steps in the interface.

```{raw} html
<figure class="video">
<video controls preload="metadata" poster="../_static/video/science.jpg">
  <source src="../_static/video/science.mp4" type="video/mp4">
  <track kind="subtitles" srclang="en" label="English" src="../_static/video/science.vtt">
</video>
<figcaption>The science and the steps, in seven narrated minutes. Chapters on the <a href="../tutorial/video.html">videos page</a>.</figcaption>
</figure>
```

| page | the question it answers |
|---|---|
| [Near-infrared spectroscopy with TripleSpec](instrument.md) | what the instrument records, and what stands between it and the spectrum you want |
| [The flat field](flat.md) | how pixel sensitivities are measured, and where the orders are |
| [Wavelength calibration](wavecal.md) | how every pixel gets a wavelength and a slit position |
| [Extraction](extraction.md) | how sky is removed, the target found, and its light summed |
| [Combining spectra](combine.md) | how many spectra become one, robustly |
| [Telluric correction and flux calibration](telluric.md) | how an A0 V star removes the atmosphere and gives physical units |
| [Merging and cleaning](merge.md) | how five orders become one spectrum |

The method is that of **Spextool** (Cushing, Vacca & Rayner 2004) and its telluric correction
(Vacca, Cushing & Rayner 2003), with K. Allers' modifications for TripleSpec 4.1 at SOAR. The
workflow and its choices follow L. G. Dahmer-Hahn's TripleSpec reduction tutorial. Every paper,
manual and page the method comes from is listed on [Sources](../reference/sources.md).

## The whole reduction in one line

$$
F_\lambda^{\rm target} \;=\;
\frac{\big(A-B\big)_{\rm target}}{\rm flat}\Big|_{\rm extracted,\ combined}
\;\times\;
\underbrace{\frac{F_\lambda^{\rm Vega\ model} \ast K \;\big|_{V,\ B-V}}{\big(A-B\big)_{\rm A0V\ star}/\rm flat\,|_{\rm extracted,\ combined}}}_{\text{telluric spectrum}}
$$

The target's counts, divided by the flat, sky-subtracted, extracted and combined, are
multiplied by the *telluric spectrum*: a model of the A0 V star (Vega's spectrum, blurred by
a kernel $K$ to the star's resolution and scaled to its magnitudes) divided by what was
actually observed of it. Everything the atmosphere and the instrument did to the star cancels,
and so does what they did to the target — leaving its true spectrum in physical units.

```{toctree}
:hidden:

instrument
flat
wavecal
extraction
combine
telluric
merge
```
