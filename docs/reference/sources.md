# Sources

Where the method comes from. TripleSpec Reducer follows the steps and the science of Spextool 4.1
and its TripleSpec 4.1 module for SOAR; everything it does is described in the papers, manuals and
pages below. Each entry says what it supports. Where a value or a file has no published source, the
entry says so.

## The method

**Cushing, Vacca & Rayner (2004).** Cushing, M. C., Vacca, W. D., & Rayner, J. T. 2004,
"Spextool: A Spectral Extraction Package for SpeX, a 0.8–5.5 Micron Cross-Dispersed Spectrograph",
*PASP*, 116, 362–376.
[doi:10.1086/382907](https://doi.org/10.1086/382907) ·
[2004PASP..116..362C](https://ui.adsabs.harvard.edu/abs/2004PASP..116..362C/abstract)

The reduction as a whole: flat fielding, locating the orders, wavelength calibration from sky lines, pair subtraction, aperture extraction with background subtraction, combining and merging.

**Vacca, Cushing & Rayner (2003).** Vacca, W. D., Cushing, M. C., & Rayner, J. T. 2003, "A Method
of Correcting Near-Infrared Spectra for Telluric Absorption", *PASP*, 115, 389–409.
[doi:10.1086/346193](https://doi.org/10.1086/346193) ·
[2003PASP..115..389V](https://ui.adsabs.harvard.edu/abs/2003PASP..115..389V/abstract) ·
[astro-ph/0211255](https://arxiv.org/abs/astro-ph/0211255)

The telluric correction and flux calibration with an A0 V star and a model of Vega: the convolution kernel from an H I line, the velocity shift, scaling to the star's V magnitude and reddening to its B−V, scaling of the H I lines, and the instrument-profile form of their eq. 13.

**Kirkpatrick et al. (2011).** Kirkpatrick, J. D., Cushing, M. C., Gelino, C. R., et al. 2011,
"The First Hundred Brown Dwarfs Discovered by the Wide-field Infrared Survey Explorer (WISE)",
*ApJS*, 197, 19, §4.4 (and §4.5).
[doi:10.1088/0067-0049/197/2/19](https://doi.org/10.1088/0067-0049/197/2/19) ·
[2011ApJS..197...19K](https://ui.adsabs.harvard.edu/abs/2011ApJS..197...19K/abstract) ·
[arXiv:1108.4677](https://arxiv.org/abs/1108.4677)

The changes Spextool needed for spectrographs whose slit is tilted on the detector. §4.4 (FIRE) describes the two-dimensional wavelength solution from OH lines, with the optical distortions mapped along both axes, and extraction through pseudo-rectangles cut out with polygon clipping; §4.5 says TripleSpec needed the same changes.

**Smith et al. (2007).** Smith, J. D. T., Armus, L., Dale, D. A., et al. 2007, "Spectral Mapping
Reconstruction of Extended Sources", *PASP*, 119, 1133–1144.
[doi:10.1086/522634](https://doi.org/10.1086/522634) ·
[2007PASP..119.1133S](https://ui.adsabs.harvard.edu/abs/2007PASP..119.1133S/abstract)

Extraction over curved, tilted orders through polygon-clipped pseudo-rectangles (cited by Kirkpatrick et al. 2011, §4.4).

**Sutherland & Hodgman (1974).** Sutherland, I. E., & Hodgman, G. W. 1974, "Reentrant polygon
clipping", *Communications of the ACM*, 17, 32–42.
[doi:10.1145/360767.360802](https://doi.org/10.1145/360767.360802)

The overlap of each pixel with the extraction polygons.

**Horne (1986).** Horne, K. 1986, "An Optimal Extraction Algorithm for CCD Spectroscopy", *PASP*,
98, 609–617. [doi:10.1086/131801](https://doi.org/10.1086/131801)

Listed only to say that optimal extraction is *not* used: the TripleSpec module supports aperture (summed) extraction only.

## The instrument

**Herter et al. (2020).** Herter, T. L., Henderson, C., Bonati, M., et al. 2020, "TripleSpec 4:
near-IR spectroscopy for the SOAR telescope", *Proc. SPIE*, 11447, 114476L.
[doi:10.1117/12.2563035](https://doi.org/10.1117/12.2563035) ·
[2020SPIE11447E..6LH](https://ui.adsabs.harvard.edu/abs/2020SPIE11447E..6LH/abstract)

TripleSpec 4.1, the instrument moved from the Blanco to SOAR.

**Schlawin et al. (2014).** Schlawin, E., Herter, T. L., Henderson, C., et al. 2014, "Design
updates and status of the fourth generation TripleSpec spectrograph", *Proc. SPIE*, 9147, 91472H.
[doi:10.1117/12.2055233](https://doi.org/10.1117/12.2055233)

TripleSpec 4 (ARCoIRIS on the Blanco), which became TripleSpec 4.1.

**Herter et al. (2008).** Herter, T. L., Henderson, C. P., Wilson, J. C., et al. 2008, "The
performance of TripleSpec at Palomar", *Proc. SPIE*, 7014, 70140X.
[doi:10.1117/12.789660](https://doi.org/10.1117/12.789660)

**Wilson et al. (2004).** Wilson, J. C., Henderson, C. P., Herter, T. L., et al. 2004, "Mass
producing an efficient NIR spectrograph", *Proc. SPIE*, 5492, 1295.
[doi:10.1117/12.550925](https://doi.org/10.1117/12.550925)

The original TripleSpec cross-dispersed design.

**Rayner et al. (2003).** Rayner, J. T., Toomey, D. W., Onaka, P. M., et al. 2003, "SpeX: A
Medium-Resolution 0.8–5.5 Micron Spectrograph and Imager for the NASA Infrared Telescope Facility",
*PASP*, 115, 362–382. [doi:10.1086/367745](https://doi.org/10.1086/367745)

The instrument Spextool was first written for; background only.

**NOIRLab, TripleSpec 4.1.**
[Instrument pages](https://noirlab.edu/science/programs/ctio/instruments/triplespec41-nir-imaging-spectrograph)
(slit 1.1″ × 28″, 0.94–2.46 µm, R ≈ 3500; its
[differences from the Blanco](https://noirlab.edu/science/programs/ctio/instruments/triplespec41-nir-imaging-spectrograph/differences-blanco),
[calibration](https://noirlab.edu/science/programs/ctio/instruments/triplespec41-nir-imaging-spectrograph/instrument-calibration)
and [history](https://noirlab.edu/science/programs/ctio/instruments/triplespec41-nir-imaging-spectrograph/history)),
and Points, S. 2020, "TripleSpec 4.1" (slides, NOIRLab),
[the step-by-step guide](https://noirlab.edu/science/programs/ctio/instruments/triplespec41-nir-imaging-spectrograph/TripleSpec4.1-step-by-step-guide).
The 0.80–2.47 µm quoted in some places is the instrument's range on the Blanco; at SOAR a dichroic
cuts the light below about 0.95 µm.

## Data and models used in the reduction

**The Vega model.** Kurucz, R. L., [Vega models](http://kurucz.harvard.edu/stars/vega/)
(T<sub>eff</sub> 9550 K, log g 3.95, v<sub>rot</sub> 25 km/s, R = 500 000, binned by 5), as used by
Vacca et al. (2003), who cite only the web page. The parameters are those published by
Castelli, F., & Kurucz, R. L. 1994, "Model atmospheres for Vega", *A&A*, 281, 817
([1994A&A...281..817C](https://ui.adsabs.harvard.edu/abs/1994A%26A...281..817C/abstract)).

**Mégessier (1995).** Mégessier, C. 1995, "Accuracy of the astrophysical absolute flux
calibrations: visible and near-infrared", *A&A*, 296, 771
([1995A&A...296..771M](https://ui.adsabs.harvard.edu/abs/1995A%26A...296..771M/abstract)).

The absolute flux of the Vega model, scaled to Vega at 5556 Å.

**Vega's V = 0.03 and B−V = 0.00.** Spextool takes these values without a stated source. V = 0.03
for Vega is the value commonly adopted, e.g. by Bessell, M. S., Castelli, F., & Plez, B. 1998,
*A&A*, 333, 231 ([1998A&A...333..231B](https://ui.adsabs.harvard.edu/abs/1998A%26A...333..231B/abstract)).

**SIMBAD.** Wenger, M., Ochsenbein, F., Egret, D., et al. 2000, "The SIMBAD astronomical
database", *A&AS*, 143, 9. [doi:10.1051/aas:2000332](https://doi.org/10.1051/aas:2000332)

The B and V magnitudes and spectral types of the A0 V standards, as Vacca et al. (2003) adopt.

**Cardelli, Clayton & Mathis (1989).** Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989,
"The Relationship between Infrared, Optical, and Ultraviolet Extinction", *ApJ*, 345, 245.
[doi:10.1086/167900](https://doi.org/10.1086/167900)

Reddening the Vega model to the standard's E(B−V), with R<sub>V</sub> = 3.1. Over TripleSpec's range only its infrared power law applies (their eqs. 2a, 2b), which was fitted to the data of Rieke, G. H., & Lebofsky, M. J. 1985, *ApJ*, 288, 618 ([doi:10.1086/162827](https://doi.org/10.1086/162827)), the curve Vacca et al. (2003) name. Below 0.909 µm the optical coefficients are those of O'Donnell, J. E. 1994, *ApJ*, 422, 158 ([doi:10.1086/173713](https://doi.org/10.1086/173713)).

**The OH sky-line list.** Made by M. Cushing (2008) and extended to 0.8 µm by K. Allers (2016),
in vacuum wavelengths, from:
Rousselot, P., Lidman, C., Cuby, J.-G., Moreels, G., & Monnet, G. 2000, "Night-sky spectral atlas
of OH emission lines in the near-infrared", *A&A*, 354, 1134
([2000A&A...354.1134R](https://ui.adsabs.harvard.edu/abs/2000A%26A...354.1134R/abstract)), and
Osterbrock, D. E., Fulbright, J. P., & Bida, T. A. 1997, "Night-Sky High-Resolution Spectral Atlas
of OH Emission Lines for Echelle Spectrograph Wavelength Calibration. II.", *PASP*, 109, 614
([doi:10.1086/133920](https://doi.org/10.1086/133920)), converted from air to vacuum with the
relation of Ciddor, P. E. 1996, *Applied Optics*, 35, 1566
([doi:10.1364/AO.35.001566](https://doi.org/10.1364/AO.35.001566)).

**The arc-lamp line list** (CuHeAr) comes with the TripleSpec module (K. Allers, 2015) without a
recorded source.

**H I wavelengths.** From the Atomic Line List, [linelist.pa.uky.edu](https://linelist.pa.uky.edu/atomic/);
van Hoof, P. A. M. 2018, "Recent Development of the Atomic Line List", *Galaxies*, 6, 63
([doi:10.3390/galaxies6020063](https://doi.org/10.3390/galaxies6020063)).

**ATRAN.** Lord, S. D. 1992, "A new software tool for computing Earth's atmospheric transmission
of near- and far-infrared radiation", NASA Technical Memorandum 103957
([NTRS](https://ntrs.nasa.gov/citations/19930010877)).

The atmospheric transmission drawn behind the telluric plots, as a guide; it does not enter the correction.

## Numerical methods

**Savitzky & Golay (1964).** Savitzky, A., & Golay, M. J. E. 1964, "Smoothing and Differentiation
of Data by Simplified Least Squares Procedures", *Analytical Chemistry*, 36, 1627.
[doi:10.1021/ac60214a047](https://doi.org/10.1021/ac60214a047)

The smoothing in the last step, computed as in Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. 1992, *Numerical Recipes*, 2nd ed. (Cambridge University Press), §14.8.

**Markwardt (2009).** Markwardt, C. B. 2009, "Non-linear Least Squares Fitting in IDL with MPFIT",
*ASP Conf. Ser.*, 411, 251 ([arXiv:0902.2850](https://arxiv.org/abs/0902.2850)), and
**Moré (1978)**, "The Levenberg-Marquardt algorithm: Implementation and theory", *Lecture Notes in
Mathematics*, 630, 105 ([doi:10.1007/BFb0067700](https://doi.org/10.1007/BFb0067700));
[MINPACK](https://www.netlib.org/minpack/).

The Levenberg–Marquardt fits.

**Bevington (1969).** Bevington, P. R. 1969, *Data Reduction and Error Analysis for the Physical
Sciences* (McGraw-Hill).

The gradient-expansion least-squares fit used for line profiles.

**Landsman (1993).** Landsman, W. B. 1993, "The IDL Astronomy User's Library", *ASP Conf. Ser.*,
52, 246 ([library](https://asd.gsfc.nasa.gov/archive/idlastro/)).

The reddening, air-to-vacuum and FITS routines Spextool relies on, whose behaviour this program matches.

**FITS.** Hanisch, R. J., Farris, A., Greisen, E. W., et al. 2001, "Definition of the Flexible
Image Transport System (FITS)", *A&A*, 376, 359.
[doi:10.1051/0004-6361:20010923](https://doi.org/10.1051/0004-6361:20010923)

## Manuals, tutorials and software

* Cushing, M., Vacca, W., & Rayner, J., *Spextool User's Manual*, version 4.1 (2016), with §5 on
  telluric correction (xtellcor), and the *Xtellcor general User's Manual*, version 4.0 (2015). Both
  come with the Spextool 4.1 distribution.
* [Spextool 4.1](https://irtfweb.ifa.hawaii.edu/~spex/Spextool/Spextool_v4.1.tar.gz), from the
  [IRTF data reduction resources](https://irtfweb.ifa.hawaii.edu/research/dr_resources/). Its
  successor for SpeX is [pySpextool](https://github.com/pyspextool/pyspextool).
* NOIRLab, [Reducing TSpec data](https://noirlab.edu/science/observing-noirlab/observing-ctio/observing-soar/data-reduction/triplespec-data):
  the Spextool version modified by K. Allers for SOAR, its download and four instructional videos.
* NOIRLab, [ARCoIRIS data reduction](https://noirlab.edu/science/programs/ctio/retired-instruments/ARCoIRIS/Data-Reduction),
  and K. Allers' [video series](https://www.youtube.com/playlist?list=PLxY2tf8cmN8Qep9tOgn_OWritFsTmkNyJ).
* The README of the TripleSpec module (K. Allers, 2016), which asks that reductions be cited as
  done with "a modified version of Spextool" (see [Citing](citing.md)).
* L. G. Dahmer-Hahn, *SOAR TripleSpec Data Reduction Tutorial* (video, 2025), whose workflow and
  choices this program's defaults follow.
