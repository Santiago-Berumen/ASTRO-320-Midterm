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

scale = 1.239                                   # arcsec/pixel (2x2 binning)
focus = np.array([17000, 16500, 16000, 15500, 15000, 14500, 14000])   # microns
fwhm  = np.array([12.53, 10.18, 8.50, 8.49, 7.86, 13.96, 14.66])      # pixels
err   = 0.27                                    # px, scatter of repeated frames

order = np.argsort(focus)
fig, ax = plt.subplots(figsize=(5, 3.8))
ax.plot(focus[order], fwhm[order], '--', color='k', lw=0.8)
ax.errorbar(focus, fwhm, yerr=err, fmt='o', ms=4, color='tab:blue',
            ecolor='k', elinewidth=0.8, capsize=2)
ax.set_xlabel(r'Focus position ($\mu$m)')
ax.set_ylabel('FWHM (pixels)')
ax.set_title(r'FWHM vs. Focus Position ($o$ Aql)')

sec = ax.secondary_yaxis('right', functions=(lambda p: p * scale, lambda a: a / scale))
sec.set_ylabel('FWHM (arcsec)')

plt.tight_layout()
plt.savefig('fwhm_vs_focus.pdf'); plt.savefig('fwhm_vs_focus.png', dpi=300)
plt.show()
