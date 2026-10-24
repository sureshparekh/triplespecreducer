# 6. Merging the orders

:::{admonition} Why
:class: note
The corrected spectrum is still five separate orders, each covering part of the
0.94–2.47 µm range, overlapping their neighbours at the ends. Merging joins them into one
continuous spectrum. After the telluric correction the orders are already on one flux
scale, so they are joined as they are.
:::

## What you do

1. Go to **6 Merge** and choose **NGC134**.
2. Leave the **Anchor** at order 3, and every order's **scale** at 1.
3. Press **Merge**.

```{shots} 06-merge
:alt: The merged spectrum

The merged spectrum
```

The stage shows the five orders in colour and the merged spectrum on top. Where two orders
overlap, the merged spectrum is their error-weighted mean.

:::{admonition} Check
:class: tip
The orders should meet without a step. If one sits higher or lower than its neighbour,
type `auto` as its scale, drag over a stretch of the overlap to set its region, and merge
again — it is then scaled to match. A noisy order end can be cut with a trim, e.g.
`green right 1.85` (the order being added, its right end, at 1.85 µm).
:::

**Next:** [7. Cleaning and saving](7-clean.md) · *The science:* [merging](../science/merge.md)
