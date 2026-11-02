#!/usr/bin/env python3
"""
Exercise E: 5 x 2 grid of the 10 double stars

Fill in the FITS file name of each pair below.
Each panel is centered on the brightest star, sized from the catalog separation, and
stretched with asinh so both stars show. A scale bar marks 10" (using 1.239"/px).
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import median_filter


FOLDER = Path("/Users/santiago/Documents/Astronomy/Astro 320/Exercise E")   # where the FITS files are
PLATE = 1.239                                                               # "/px (binned 2x2)

# name shown on panel : (FITS file name, catalog separation in arcsec)
PAIRS = {
    r"HR 8384":            ("HR8384.fits",        91),
    r"61 Cygni":           ("61Cyg.fits",         32),
    r"$\zeta$ Lyrae":      ("zetaLyr.fits",       43),
    r"$\epsilon$ Lyrae":   ("epsLyr.fits",       208),
    r"$\beta$ Cygni":      ("betaCyg.fits",      34.3),
    r"$\nu$ Draconis":     ("nuDra.fits",        61.6),
    r"40/41 Draconis":     ("40_41Dra.fits",      19),
    r"$\delta$ Cephei":    ("deltaCep.fits",      41),
    r"8 Lacertae":         ("8Lac.fits",          22),
    r"100 Herculis":       ("100Her.fits",        14),
}

CENTER = {}


plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 10})
fig, axes = plt.subplots(2, 5, figsize=(12, 5.4))

for ax, (name, (fname, sep)) in zip(axes.flat, PAIRS.items()):
    img = fits.getdata(FOLDER / fname).astype(float)
    img = np.flipud(img)                                
    sky = np.median(img)
    sep_px = sep / PLATE

    if name in CENTER:
        x0, y0 = CENTER[name]
    else:                                              
        sm = median_filter(img, size=5)
        y0, x0 = np.unravel_index(np.argmax(sm), sm.shape)

    half = int(max(30, 1.6 * sep_px))                   
    ny, nx = img.shape
    ya, yb = max(0, y0 - half), min(ny, y0 + half)
    xa, xb = max(0, x0 - half), min(nx, x0 + half)
    cut = img[ya:yb, xa:xb] - sky

    noise = 1.4826 * np.median(np.abs(cut - np.median(cut)))
    stretched = np.arcsinh(np.clip(cut, 0, None) / (3 * noise))
    ax.imshow(stretched, cmap='gray', origin='upper',
              vmin=0, vmax=np.percentile(stretched, 99.9))

    # 10" scale bar
    bar = 10 / PLATE
    w = cut.shape[1]
    ax.plot([0.07 * w, 0.07 * w + bar], [0.92 * cut.shape[0]] * 2, color='w', lw=2)
    ax.text(0.07 * w + bar / 2, 0.88 * cut.shape[0], '10$^{\\prime\\prime}$', color='w',
            ha='center', va='bottom', fontsize=8)

    ax.set_title(f"{name}  ({sep}$^{{\\prime\\prime}}$)", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

plt.tight_layout()
plt.savefig(FOLDER / "double_stars.png", dpi=300)
plt.savefig(FOLDER / "double_stars.pdf")
plt.show()
print(f"Saved {FOLDER / 'double_stars.png'}")
