import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 10,
    'axes.linewidth': 0.8,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'legend.frameon': False,
})

# ================= DATA (from AstroImageJ Multi-Aperture, r = 15/25/60) =================
# Exposure series, night 1, Images 53-63, TYC 4004-1801-1 field
t = np.array([5, 10, 20, 30, 40, 50, 60, 80, 100, 120, 180])

# Star A: faint, never saturates (peak 1,248 -> 31,699 DN)
rate_A = np.array([1784.1, 1527.7, 1476.8, 1486.4, 1508.3, 1483.9, 1476.7, 1519.8, 1506.2, 1526.1, 1539.8])
err_A  = np.array([84.4, 59.9, 42.5, 34.7, 30.1, 26.7, 24.3, 21.1, 18.9, 17.3, 14.1])
fwhm_A = np.array([7.91, 3.19, 3.05, 2.61, 2.75, 2.66, 2.51, 2.71, 2.80, 2.90, 2.98])

# Star B: bright, saturates from 60 s on (peak 3,201 -> 49,683 DN)
rate_B = np.array([4561.8, 4685.0, 4803.3, 4794.4, 4882.8, 4889.4, 4876.8, 4807.7, 4675.6, 4592.7, 4287.5])
err_B  = np.array([87.5, 62.4, 44.4, 36.3, 31.5, 27.9, 25.5, 22.0, 19.7, 18.0, 14.7])
fwhm_B = np.array([2.99, 2.80, 2.69, 2.73, 2.72, 2.46, 2.73, 2.81, 3.04, 3.32, 4.22])
t_sat = 60   # first exposure where Star B's peak hits ~47,500 DN

# Repeatability, night 2, 10 x 10 s, Star C
rate_C = np.array([32257.1, 32303.0, 32323.7, 32091.4, 32704.6, 32486.2, 28002.0, 32320.9, 32713.4, 32605.8])
err_C  = np.array([91.4, 91.3, 91.4, 91.2, 91.5, 91.4, 88.9, 91.2, 91.5, 91.4])
fwhm_C = np.array([2.89, 4.05, 3.59, 4.41, 2.79, 3.28, 5.63, 2.92, 3.43, 4.06])
bad = 6      # frame 7: wind-smeared (width 12.7 px), excluded

# FWHM uncertainty: frame-to-frame scatter of Star C at fixed exposure (seeing changes)
good = np.arange(10) != bad
fwhm_err = fwhm_C[good].std(ddof=1)

# ================= PART 3: repeatability =================
frames = np.arange(1, 11)
mean_C = rate_C[good].mean()
dev = 100 * (rate_C - mean_C) / mean_C
dev_err = 100 * err_C / mean_C
scatter = dev[good].std(ddof=1)
print(f"Part 3: mean = {mean_C:.0f} +/- {rate_C[good].std(ddof=1)/np.sqrt(good.sum()):.0f} DN/s, "
      f"scatter = {scatter:.2f}%, photon error = {dev_err.mean():.2f}%")

fig, ax = plt.subplots(figsize=(5, 3.8))
ax.axhline(0, color='k', lw=0.8)
ax.axhspan(-scatter, scatter, color='gray', alpha=0.15, lw=0, label=rf'$\pm 1\sigma$ scatter ({scatter:.2f}%)')
ax.errorbar(frames[good], dev[good], yerr=dev_err[good], fmt='o', ms=4, color='tab:blue',
            ecolor='k', elinewidth=0.8, capsize=2, label='photon-noise error bars')
ax.set_xlabel('Frame number')
ax.set_ylabel('Deviation from mean (%)')
ax.set_title('Repeatability of Photometry (10 s exposures)')
ax.set_xticks(frames)
ax.legend(fontsize=8, loc='upper left')
plt.tight_layout()
plt.savefig('repeatability.pdf'); plt.savefig('repeatability.png', dpi=300)
plt.show()

# ================= PART 4: FWHM vs exposure time =================
fig, ax = plt.subplots(figsize=(5, 3.8))
ax.errorbar(t, fwhm_A, yerr=fwhm_err, fmt='o', ms=4, color='tab:blue',
            ecolor='k', elinewidth=0.8, capsize=2, label='Star A (unsaturated)')
ax.errorbar(t, fwhm_B, yerr=fwhm_err, fmt='s', ms=4, color='tab:orange',
            ecolor='k', elinewidth=0.8, capsize=2, label='Star B (saturates)')
ax.axvline(t_sat, color='gray', ls=':', lw=0.8)
ax.set_xlabel('Exposure time (s)')
ax.set_ylabel('FWHM (pixels)')
ax.set_title('FWHM vs. Exposure Time')
ax.legend(fontsize=8, loc='upper right')
plt.tight_layout()
plt.savefig('fwhm_vs_exptime.pdf'); plt.savefig('fwhm_vs_exptime.png', dpi=300)
plt.show()

# ================= PART 6: count rate vs exposure time =================
fig, ax = plt.subplots(figsize=(5, 3.8))
ax.errorbar(t, rate_A, yerr=err_A, fmt='o', ms=4, color='tab:blue',
            ecolor='k', elinewidth=0.8, capsize=2, label='Star A (unsaturated)')
ax.errorbar(t, rate_B, yerr=err_B, fmt='s', ms=4, color='tab:orange',
            ecolor='k', elinewidth=0.8, capsize=2, label='Star B (saturates)')
ax.axvline(t_sat, color='gray', ls=':', lw=0.8)
ax.set_xlabel('Exposure time (s)')
ax.set_ylabel('Count rate (DN/s)')
ax.set_title('Count Rate vs. Exposure Time')
ax.legend(fontsize=8, loc='center right')
plt.tight_layout()
plt.savefig('countrate_vs_exptime.pdf'); plt.savefig('countrate_vs_exptime.png', dpi=300)
plt.show()
