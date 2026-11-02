# Astro 320 Observing Project: CCD Characterization and Photometry

Code and analysis for my Astro 320 (observational astronomy) midterm report at Penn State. The project characterizes the CCD camera on the 24-inch telescope and then uses it for stellar photometry, color photometry, astrophotography and astrometry, split into five exercises (A–E).

Most of the image handling (calibration, alignment, aperture photometry, plot profiles) was done in **AstroImageJ**. The Python scripts here do the parts AstroImageJ can't, like fitting, custom calibration and the plots in the report.

**Requirements:** Python 3 with `numpy`, `scipy`, `astropy`, `matplotlib`, `pandas`

```bash
pip install numpy scipy astropy matplotlib pandas
```

> Scripts that read FITS files or tables have a SETTINGS block at the top with file paths from my computer (`/Users/santiago/...`). Change those to point at your own data before running. The plotting scripts for Exercise C have the measurements typed in and run as is.

---

## Exercise A: Read Noise and Dark Current

### Full-frame read noise script (provided by the course)
Measures the read noise over the whole chip from the 20 bias frames.
- Compares each pixel only with itself across all 20 frames, so fixed patterns in the bias (the gradient, dark columns, bright edge) cancel out and only the random read noise is left.
- **Result:** 9.637 DN (11.57 e⁻)

### ROI read noise script (provided by the course)
Does the same measurement in a 110 × 110 pixel defect-free region (ROI) centered at (x = 2300, y = 2030).
- **Result:** 9.636 DN, which agrees with the full frame to 0.001 DN and shows the read noise is uniform across the detector.

### `dark_current.py`
Measures the dark current from the dark frames and finds the activation energy.
- Subtracts the master bias from every dark and takes the median as the dark level. It trims the edges to avoid the bright strip in the bias.
- Reads the exposure time and CCD temperature from each FITS header.
- **Time series (~7 °C):** fits a line to dark level vs. exposure time. The slope × gain gives the dark current (0.436 e⁻/pixel/s), and the intercept soaks up the leftover offset from the master bias being taken at a colder temperature.
- **Temperature series (300 s darks):** removes the temperature-dependent offset using the coldest frames, then fits ln *I* vs. 1/*T* (the Arrhenius form *I* = *C* e<sup>−A/kT</sup>) to get the activation energy, A = 1.01 ± 0.02 (stat) ± 0.12 (sys) eV.
- Estimates the systematic error by trying a second offset model.
- **Outputs:** `dark_vs_time`, `dark_vs_temp`, `dark_arrhenius` (PDF + PNG)

---

## Exercise B: Flat Fields and Gain

No scripts here. The gain (*g* = 1.201 ± 0.003 e⁻/DN) came from a photon transfer curve of the dome flats. That script was not saved. The 1-D illumination profiles of the dome and sky flats were made with AstroImageJ's *Plot Profile*.

---

## Exercise C: Image Quality and Photometric Uncertainty

### `calibrate.py`
My full calibration pipeline, used for Exercises C and D:

**calibrated = (raw − bias − dark) / flat**

- **Bias:** the Exercise A master bias was taken at 1×1 binning, so the script block-averages it 2×2 to match the science frames.
- **Dark:** no 2×2 darks were taken, so the dark is computed from the Exercise A dark-current law at each frame's CCD temperature.
- **Flat:** normalizes the master flat and checks that its `FILTER` header matches the science frames.
- **Orientation:** automatically detects whether the calibration frames are flipped relative to the raw images.
- **Defects:** finds and repairs hot pixels and bad columns.
- **Inputs:** each calibration input can be a FITS file or a folder holding one.
- **Outputs:** `<name>_calibrated.fit` files plus `calibration_log.csv`, `hot_pixels.csv` and `calibration_summary.txt`

### `focus_curve.py`
Plots FWHM vs. focus position for *o* Aquilae, from 17,000 to 14,000 µm in 500 µm steps, measured with AstroImageJ's *Plot Profile*.
- Error bars of ±0.27 px come from the scatter between repeated images at the same focus.
- Adds a second y-axis in arcseconds using the plate scale from Exercise E (1.239″/px).
- The points form a V shape, with the sharpest image at 15,000 µm (7.9 px ≈ 9.7″).
- **Output:** `fwhm_vs_focus` (PDF + PNG)

### `photometry_tests.py`
Makes three plots from AstroImageJ Multi-Aperture photometry (r = 15/25/60 px).
- **Repeatability:** count rate of one unsaturated star in ten 10 s exposures, with the wind-smeared frame 7 dropped. Computes each frame's percent deviation from the mean, the frame-to-frame scatter and the photon-noise error. The scatter is larger than the photon noise, so the measurements are not photon-limited.
- **FWHM vs. exposure time:** exposure series of the TYC 4004-1801-1 field (5–180 s), comparing a faint star that never saturates with a brighter one that saturates from 60 s on. The unsaturated star's FWHM stays constant, while the saturated star's grows.
- **Count rate vs. exposure time:** same two stars. The unsaturated star's count rate stays constant, as expected for a linear detector, while the saturated star's drops.
- **Outputs:** `repeatability`, `fwhm_vs_exptime`, `countrate_vs_exptime` (PDF + PNG)

---

## Exercise D: Color Photometry and Astrophotography

### `color_color.py`
Turns the AstroImageJ Multi-Aperture tables for V, R and I (10 stars in NGC 6811) into a photometry table and diagrams.
- Computes instrumental magnitudes, *m* = −2.5 log<sub>10</sub>(*N*/*t*), with uncertainties σ<sub>m</sub> = 1.086 σ<sub>N</sub>/*N*.
- Handles both target (T1) and comparison (C2…C10) aperture labels.
- Puts V on a real magnitude scale with a zero point from the median offset of the non-variable comparison stars (21.40 ± 0.07 mag). The two δ Scuti variables are left out.
- **Outputs:** `ngc6811_photometry.csv`, a LaTeX table (`ngc6811_photometry.tex`), the color-color diagram (`color_color`, V−R vs. R−I) and the color-magnitude diagram (`cmd_V_VR`, V vs. V−R), with stars colored from hot (blue) to cool (red)

### `ngc6811_targets.txt`
RA/Dec (degrees) of the 10 stars measured in NGC 6811.

The Fireworks Galaxy (NGC 6946) color image was made entirely in AstroImageJ: calibrate each filter, align and median-combine, align the three filters, stretch, then *Merge Channels* (I → red, R → green, V → blue).

---

## Exercise E: Plate Scale

### `double_stars_grid.py`
Makes a 5 × 2 grid of cutouts of the 10 double stars. Each panel is centered on the brighter star, sized from the catalog separation, asinh-stretched and labeled with a 10″ scale bar.
- **Output:** `double_stars` (PDF + PNG)

Star centroids and FWHMs were measured in AstroImageJ (Multi-Aperture with centroiding on). Separations come from the Washington Double Star Catalog, accessed through Stelle Doppie.

---

## Key Results

| Quantity | Value |
|---|---|
| Read noise | 9.64 DN (11.6 e⁻) |
| Dark current at 7 °C | 0.436 ± 0.003 e⁻/pixel/s |
| Activation energy | 1.01 ± 0.02 (stat) ± 0.12 (sys) eV |
| Gain | 1.201 ± 0.003 e⁻/DN |
| Full well | ~44,000 DN (~53,000 e⁻) |
| Photometric repeatability | 0.68% (photon noise 0.28%) |
| Plate scale (2×2 binning) | 1.239 ± 0.006 ″/px |
| Seeing | 3.5 ± 0.3″ |
| Field of view | 31.5 ± 0.2′ |
