import glob
import re
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


dark_folder = '/Users/santiago/Documents/Astronomy/Astro 320/Exercise A/dark_frames/'
mbias_path  = '/Users/santiago/Documents/Astronomy/Astro 320/Exercise A/Bias_Subtracted/masterbias.fits'
edge = 50               
g, g_err = 1.201, 0.003  
k = 8.617e-5          
T_cold = -15    

def num(f):
    m = re.search(r'(\d+)\.fits?$', f)
    return int(m.group(1)) if m else 0


mbias = fits.getdata(mbias_path).astype(float)
print(f"Master bias: {mbias.shape}, median = {np.median(mbias):.1f} DN\n")


files = sorted(glob.glob(dark_folder + '*.fit*'), key=num)
print(f"{'file':<22}{'EXPTIME':>9}{'CCD-TEMP':>10}{'CALSTAT':>9}{'XBIN':>6}"
      f"{'raw med':>10}{'dark med':>10}{'dark mean':>11}")

rows = []
for f in files:
    with fits.open(f) as h:
        data = h[0].data.astype(float)
        hdr = h[0].header
    t    = hdr.get('EXPTIME', np.nan)
    temp = hdr.get('CCD-TEMP', hdr.get('CCDTEMP', np.nan))
    cal  = str(hdr.get('CALSTAT', '')).strip()
    xbin = hdr.get('XBINNING', '?')

    if data.shape != mbias.shape:
        print(f"{f.split('/')[-1]:<22}  SKIPPED: size {data.shape} doesn't match master bias")
        continue

    dark = (data - mbias)[edge:-edge, edge:-edge]
    raw_med = np.median(data[edge:-edge, edge:-edge])
    med, mean = np.median(dark), np.mean(dark)
    rows.append((f.split('/')[-1], t, temp, med, mean))
    print(f"{f.split('/')[-1]:<22}{t:>9.2f}{temp:>10.2f}{cal:>9}{str(xbin):>6}"
          f"{raw_med:>10.1f}{med:>10.2f}{mean:>11.2f}")

names    = [r[0] for r in rows]
t_all    = np.array([r[1] for r in rows])
temp_all = np.array([r[2] for r in rows])
med_all  = np.array([r[3] for r in rows])

np.save('dark_measurements.npy', np.array([r[1:] for r in rows], dtype=float))
print("\nSaved dark_measurements.npy (EXPTIME, CCD-TEMP, median, mean)")


ts = t_all <= 60

(slope, intercept), cov = np.polyfit(t_all[ts], med_all[ts], 1, cov=True)
slope_err = np.sqrt(cov[0, 0])

I = slope * g
I_err = np.sqrt((slope_err * g) ** 2 + (slope * g_err) ** 2)

print(f"\nslope = {slope:.4f} +/- {slope_err:.4f} DN/s, intercept = {intercept:.2f} DN")
print(f"dark current at ~7 C = {I:.4f} +/- {I_err:.4f} e-/pixel/s")


tt = t_all == 300             
T_C = temp_all[tt]
T_K = T_C + 273.15
med_T = med_all[tt]

# The master bias was taken at a different temperature than the darks, so each
# bias-subtracted dark has a leftover offset. We know that offset in two places:
#   - the 7 C time-series intercept (offset with zero exposure time)
#   - the cold 300 s frames, whose dark signal is negligible, so their median IS the offset
# Interpolate between those to get the offset at every temperature.
cold  = T_C < T_cold
T_off = np.r_[T_C[cold], 6.97]
N_off = np.r_[med_T[cold], intercept]
o = np.argsort(T_off)
offset = np.interp(T_C, T_off[o], N_off[o])

signal = med_T - offset              
I_T = g * signal / 300               
use = (~cold) & (signal > 3)        

print("\n   T (C)   offset (DN)   signal (DN)   I (e-/px/s)")
for Tc, off, s, i, c in zip(T_C, offset, signal, I_T, cold):
    tag = '  (offset reference)' if c else ''
    print(f"{Tc:8.2f}{off:13.2f}{s:14.2f}{i:14.4f}{tag}")

def arrhenius(TK, I):
    (m, b), c = np.polyfit(1 / TK, np.log(I), 1, cov=True)
    return m, b, np.sqrt(c[0, 0])

m, b, m_err = arrhenius(T_K[use], I_T[use])
A, A_err = -k * m, k * m_err
print(f"\nslope = {m:.0f} +/- {m_err:.0f} K")
print(f"activation energy A = {A:.2f} +/- {A_err:.2f} eV (stat)  (N = {use.sum()} points)")


off_old = np.interp(T_C, [T_C.min(), 6.97], [med_T[np.argmin(T_C)], intercept])
sig_old = med_T - off_old
use_old = sig_old > 3
m_old, _, _ = arrhenius(T_K[use_old], g * sig_old[use_old] / 300)
A_sys = abs(-k * m_old - A)
print(f"old offset model gives A = {-k * m_old:.2f} eV -> systematic ~ {A_sys:.2f} eV")
print(f"A = {A:.2f} +/- {A_err:.2f} (stat) +/- {A_sys:.2f} (sys) eV")


plt.rcParams.update({
    'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 10,
    'axes.linewidth': 0.8,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'legend.frameon': False,
})


charge_t = g * (med_all[ts] - intercept)         
tt_fit = np.linspace(0, 65, 100)
fig, ax = plt.subplots(figsize=(5, 3.8))
ax.plot(tt_fit, g * slope * tt_fit, 'k', lw=1,
        label=rf'fit, $I = {I:.3f} \pm {I_err:.3f}$ e$^-$/pixel/s')
ax.plot(t_all[ts], charge_t, 'o', ms=4, color='tab:blue', label=r'darks at $\approx 7\,^\circ$C')
ax.set_xlabel('Exposure time (s)')
ax.set_ylabel(r'Dark charge (e$^-$/pixel)')
ax.set_title('Dark Charge vs. Exposure Time')
ax.legend(fontsize=8, loc='upper left')
plt.tight_layout()
plt.savefig('dark_vs_time.pdf'); plt.savefig('dark_vs_time.png', dpi=300)
plt.show()


charge_T = g * signal                           
TT = np.linspace(-36, 9, 200)
fig, ax = plt.subplots(figsize=(5, 3.8))
ax.plot(TT, 300 * np.exp(b + m / (TT + 273.15)), 'k', lw=1, label=rf'model, $A = {A:.2f}$ eV')
ax.plot(T_C[use], charge_T[use], 'o', ms=4, color='tab:blue', label='used in fit')
ax.plot(T_C[~use], charge_T[~use], 'o', ms=4, mfc='none', color='tab:blue',
        label='offset reference')
ax.set_xlabel(r'CCD temperature ($^\circ$C)')
ax.set_ylabel(r'Dark charge in 300 s (e$^-$/pixel)')
ax.set_title('Dark Charge vs. Temperature')
ax.legend(fontsize=8, loc='upper left')
plt.tight_layout()
plt.savefig('dark_vs_temp.pdf'); plt.savefig('dark_vs_temp.png', dpi=300)
plt.show()


x = 1 / T_K[use]
xx = np.linspace(x.min() * 0.995, x.max() * 1.005, 100)
fig, ax = plt.subplots(figsize=(5, 3.8))
ax.plot(xx, m * xx + b, 'k', lw=1, label=rf'fit, $A = {A:.2f} \pm {A_err:.2f}$ eV')
ax.plot(x, np.log(I_T[use]), 'o', ms=4, color='tab:blue', label='300 s darks')
ax.set_xlabel(r'$1/T$ (K$^{-1}$)')
ax.set_ylabel(r'$\ln(I)$ [$I$ in e$^-$/pixel/s]')
ax.set_title('Linearized Dark Current vs. Temperature')
ax.legend(fontsize=8, loc='upper right')
plt.tight_layout()
plt.savefig('dark_arrhenius.pdf'); plt.savefig('dark_arrhenius.png', dpi=300)
plt.show()
