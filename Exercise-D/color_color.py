#!/usr/bin/env python3
"""
Astro 320 - Exercise D, step 4: color-color diagram of NGC 6811
================================================================
Reads the three AstroImageJ Multi-Aperture "Measurements" tables (one per filter, measured
on the summed image, same stars clicked in the same order), builds a photometry table and
plots V-R vs R-I.

- If you entered reference-star magnitudes in AstroImageJ, the calibrated magnitudes
  (Source_AMag_Tn) are used.
- Otherwise instrumental magnitudes m = -2.5 log10(Source-Sky / t) are used. Colors don't
  need the absolute calibration (the zero-point offsets just shift the whole diagram).

pip install numpy pandas matplotlib
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ================= SETTINGS =================
BASE = Path("/Users/santiago/Documents/Astronomy/Astro 320/Exercise D")

TABLES = {                                   # AstroImageJ Measurements tables (.csv or .xls)
    "V": BASE / "MagVmap.csv",
    "R": BASE / "MagRmap.csv",
    "I": BASE / "MagImap.csv",
}
EXPTIME = {"V": 48.0, "R": 48.0, "I": 48.0}  # total exposure of each summed image [s] (e.g. 8 x 6 s)
USE_CALIBRATED = False                       # False: instrumental mags from Source-Sky for every star.
# (AstroImageJ's Source_AMag for comparison stars C2, C3, ... just repeats the reference magnitude you
#  typed in, so it is NOT a measurement. Instrumental mags are measured for all stars; colors are fine.)

# Optional: names for the star numbers (T1, T2, ...). Leave empty to label 1, 2, 3, ...
NAMES = {                                    # matched from the reference mags typed into AstroImageJ
    1: "KIC 9716045",   # target T1 (its measured I = 8.92 matches the expected 8.91)
    2: "KIC 9716522",
    3: "KIC 9654907",
    4: "V2502 Cyg",
    5: "KIC 9655101?",  # reference I = 10.01 doesn't match the table: check which star this is
    6: "KIC 9655167",
    7: "V2503 Cyg",
    8: "KIC 9594437",
    9: "KIC 9716431",
    10: "KIC 9655187",
}
SKIP = []                                    # star numbers to leave out of the plot (e.g. variables)

# ================= READ =================
def read_aij(path):
    """AstroImageJ saves tab-separated text even with a .xls extension; this handles both."""
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = [c.strip() for c in df.columns]
    return df


def star_numbers(df):
    """Aperture numbers, whether target (T1) or comparison (C2, C3, ...)."""
    return sorted({int(m.group(2)) for c in df.columns if (m := re.fullmatch(r"Source-Sky_([TC])(\d+)", c))})


def col(df, stem, n):
    """Column name for aperture n, e.g. Source-Sky_T1 or Source-Sky_C7."""
    for p in ("T", "C"):
        if f"{stem}_{p}{n}" in df.columns:
            return f"{stem}_{p}{n}"
    return None


def magnitudes(df, band):
    """One row per star: magnitude and error (averaged over rows if the table has several)."""
    stars = star_numbers(df)
    calibrated = (USE_CALIBRATED is True) or (USE_CALIBRATED == "auto" and col(df, "Source_AMag", stars[0]) is not None)
    out = {}
    if len(stars) == 1 and len(df) > 1:  # (not your files: they have T1 + C2..C10 in one row)
        # one aperture column (T1) but one ROW per star: each row is a different star
        n0 = stars[0]
        for i, (_, row) in enumerate(df.iterrows(), start=1):
            if calibrated:
                m = float(row[f"Source_AMag_T{n0}"])
                e = float(row[f"Source_AMag_Err_T{n0}"]) if f"Source_AMag_Err_T{n0}" in df.columns else 0.0
            else:
                f, ferr = float(row[f"Source-Sky_T{n0}"]), float(row[f"Source_Error_T{n0}"])
                m, e = -2.5 * np.log10(f / EXPTIME[band]), 2.5 / np.log(10) * ferr / f
            out[i] = (m, e)
        return out, calibrated
    for n in stars:
        if calibrated:
            m = df[col(df, "Source_AMag", n)].astype(float)
            ce = col(df, "Source_AMag_Err", n)
            e = df[ce].astype(float) if ce else m * 0
        else:
            f = df[col(df, "Source-Sky", n)].astype(float)
            ferr = df[col(df, "Source_Error", n)].astype(float)
            m = -2.5 * np.log10(f / EXPTIME[band])
            e = 2.5 / np.log(10) * ferr / f
        out[n] = (m.mean(), np.sqrt((e ** 2).sum()) / len(e))
    return out, calibrated


mags, kinds = {}, {}
for band, path in TABLES.items():
    mags[band], kinds[band] = magnitudes(read_aij(path), band)
    print(f"{band}: {len(mags[band])} stars from {path.name} "
          f"({'calibrated Source_AMag' if kinds[band] else 'instrumental'} magnitudes)")

if len(set(kinds.values())) > 1:
    print("(!) some filters are calibrated and some instrumental: colors will have an arbitrary offset")

common = sorted(set(mags["V"]) & set(mags["R"]) & set(mags["I"]))
if len({len(m) for m in mags.values()}) > 1:
    print(f"(!) the tables have different numbers of stars; using the {len(common)} present in all three")

# ================= TABLE =================
rows = []
for n in common:
    (V, eV), (R, eR), (I, eI) = mags["V"][n], mags["R"][n], mags["I"][n]
    rows.append(dict(star=n, name=NAMES.get(n, ""), V=V, V_err=eV, R=R, R_err=eR, I=I, I_err=eI,
                     V_R=V - R, V_R_err=np.hypot(eV, eR), R_I=R - I, R_I_err=np.hypot(eR, eI),
                     V_I=V - I, V_I_err=np.hypot(eV, eI)))
tab = pd.DataFrame(rows)
tab.to_csv(BASE / "ngc6811_photometry.csv", index=False, float_format="%.3f")
print("\n" + tab.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# LaTeX table for the report
with open(BASE / "ngc6811_photometry.tex", "w") as f:
    f.write("\\begin{tabular}{clccccc}\n\\hline\n"
            "\\# & Star & $V$ & $R$ & $I$ & $V-R$ & $R-I$ \\\\\n\\hline\n")
    for r in rows:
        f.write(f"{r['star']} & {r['name']} & ${r['V']:.3f} \\pm {r['V_err']:.3f}$ & "
                f"${r['R']:.3f} \\pm {r['R_err']:.3f}$ & ${r['I']:.3f} \\pm {r['I_err']:.3f}$ & "
                f"${r['V_R']:.3f} \\pm {r['V_R_err']:.3f}$ & ${r['R_I']:.3f} \\pm {r['R_I_err']:.3f}$ \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")
print(f"\nSaved {BASE / 'ngc6811_photometry.csv'} and ngc6811_photometry.tex")

# ================= PLOT =================
plt.rcParams.update({
    'font.family': 'serif', 'mathtext.fontset': 'cm', 'font.size': 10,
    'axes.linewidth': 0.8,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'legend.frameon': False,
    'axes.grid': True, 'axes.axisbelow': True,          # grid behind the data
    'grid.color': '0.85', 'grid.linestyle': '--', 'grid.linewidth': 0.6,
})

plot = tab[~tab["star"].isin(SKIP)]
fig, ax = plt.subplots(figsize=(5, 4.2))
ax.errorbar(plot["R_I"], plot["V_R"], xerr=plot["R_I_err"], yerr=plot["V_R_err"], fmt='none',
            ecolor='k', elinewidth=0.8, capsize=2, zorder=1)
# color each star by its V-R, blue (hot) -> red (cool), same scale as the color-magnitude diagram
sc0 = ax.scatter(plot["R_I"], plot["V_R"], c=plot["V_R"], cmap="RdYlBu_r", s=45,
                 edgecolors='k', linewidths=0.5, zorder=2)
cb0 = fig.colorbar(sc0, ax=ax, pad=0.02)
cb0.set_label(r"$V - R$ (instrumental):  hot $\rightarrow$ cool")
for _, r in plot.iterrows():
    ax.annotate(str(int(r["star"])), (r["R_I"], r["V_R"]), xytext=(4, 3), textcoords="offset points", fontsize=8)
label = "" if all(kinds.values()) else r" (instrumental)"
ax.set_xlabel(r"$R - I$" + label)
ax.set_ylabel(r"$V - R$" + label)
ax.set_title("Color-Color Diagram of NGC 6811")
plt.tight_layout()
plt.savefig(BASE / "color_color.pdf"); plt.savefig(BASE / "color_color.png", dpi=300)
plt.show()

# ================= COLOR-MAGNITUDE DIAGRAM: V vs V-R =================
# Put V on a real magnitude scale with one zero point: the median offset between the catalog V you
# typed into AstroImageJ (Source_AMag of the comparison stars) and their instrumental V. The same
# offset is applied to every star, so the relative brightnesses are still the measured ones.
VAR_STARS = [4, 7]                           # variables: left out of the zero point
dfV = read_aij(TABLES["V"])
ref = []
for n in common:
    c = col(dfV, "Source_AMag", n)
    if c and c.endswith(f"_C{n}") and n not in VAR_STARS:
        ref.append(float(dfV[c].iloc[0]) - mags["V"][n][0])
zp = float(np.median(ref)) if ref else 0.0
zp_err = 1.2533 * np.std(ref, ddof=1) / np.sqrt(len(ref)) if len(ref) > 1 else 0.0
print(f"\nV zero point = {zp:.3f} +/- {zp_err:.3f} mag (from {len(ref)} comparison stars)")
tab["V_cal"] = tab["V"] + zp
tab.to_csv(BASE / "ngc6811_photometry.csv", index=False, float_format="%.3f")

plot = tab[~tab["star"].isin(SKIP)]
fig, ax = plt.subplots(figsize=(5, 4.2))
ax.errorbar(plot["V_R"], plot["V_cal"], xerr=plot["V_R_err"], yerr=np.hypot(plot["V_err"], zp_err),
            fmt='none', ecolor='k', elinewidth=0.8, capsize=2, zorder=1)
# color each star by its V-R, blue (hot) -> red (cool), like an HR diagram
sc = ax.scatter(plot["V_R"], plot["V_cal"], c=plot["V_R"], cmap="RdYlBu_r", s=45,
                edgecolors='k', linewidths=0.5, zorder=2)
cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label(r"$V - R$ (instrumental):  hot $\rightarrow$ cool")
for _, r in plot.iterrows():
    ax.annotate(str(int(r["star"])), (r["V_R"], r["V_cal"]), xytext=(4, 3), textcoords="offset points", fontsize=8)
ax.invert_yaxis()                            # brighter stars at the top
ax.set_xlabel(r"$V - R$ (instrumental)")
ax.set_ylabel(r"$V$ (mag)")
ax.set_title("Color-Magnitude Diagram of NGC 6811")
plt.tight_layout()
plt.savefig(BASE / "cmd_V_VR.pdf"); plt.savefig(BASE / "cmd_V_VR.png", dpi=300)
plt.show()
