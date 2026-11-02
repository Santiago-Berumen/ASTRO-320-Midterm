#!/usr/bin/env python3
"""
Astro 320 - calibrate star images for photometry
calibrated = (raw - pedestal - bias - dark), hot pixels/trails repaired, / normalised flat, + offset

Every FITS file in the folder is saved as <name>_calibrated.fit in <folder>/<folder>_calibrated,
together with calibration_log.csv, hot_pixels.csv and calibration_summary.txt.
"""
import csv
import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy import ndimage


BASE = Path("/Users/santiago/Documents/Astronomy/Astro 320")


INPUT_FOLDER = BASE / "Exercise D" / "partD(NGC6811)" / "I"   
MASTER_FLAT  = BASE / "Exercise D" / "masterflatI"          
MASTER_BIAS  = BASE / "Exercise D" / "masterbias"             
DARK_FRAMES_FOLDER = None  

GAIN, GAIN_ERR = 1.201, 0.003    
RN_DN, RN_ERR = 9.637, 0.014    


TB, TB_ERR = 11778.0, 212.0
I0 = 6.596e17
I0_ERR = 0.006 * I0         
DARK_7C, DARK_7C_ERR = 0.4358, 0.0025  
DARK_MODE = "law"                

SUBTRACT_PEDESTAL = True
BIAS_SMOOTH = 31                  
REF_BIAS_LEVEL_RAW = 1026.9      
FLAT_MIN = 0.1                  


FLAT_ORIENT = "auto"            
BIAS_ORIENT = "auto"

# hot pixels / trails
FIX_HOT_PIXELS = True
HOT_SIGMA = 8.0
SHARP_FRACTION = 0.3
HOT_FRACTION = 0.5
TRAIL_SIGMA = 1.5
TRAIL_AXIS = "auto"               

ADD_OFFSET = "auto"              

TEMP_KEYS = ["CCD-TEMP", "CCDTEMP", "TEMPERAT"]
EXPTIME_KEYS = ["EXPTIME", "EXPOSURE"]
ADC_DN, FULL_WELL_DN = 65535.0, 44000.0   



def natural_key(p):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", p.name)]


def find_fits(folder):
    files = [p for p in Path(folder).iterdir()
             if p.is_file() and p.suffix.lower() in (".fit", ".fits", ".fts") and "_calibrated" not in p.stem]
    if not files:
        raise FileNotFoundError(f"No FITS files in {folder}")
    return sorted(files, key=natural_key)


def read_fits(path):
    with fits.open(path) as h:
        return h[0].data.astype(np.float64), h[0].header.copy()


def get(header, keys, default=None):
    for k in keys:
        if k in header:
            return float(header[k])
    return default


def rebin_mean(a, fy, fx):
    """Block-average by integer factors (a bias level is an offset, so it is averaged, not summed)."""
    ny, nx = a.shape[0] // fy, a.shape[1] // fx
    return a[:ny * fy, :nx * fx].reshape(ny, fy, nx, fx).mean(axis=(1, 3))


def resolve_fits(path):
    """Accept a FITS file, or a folder that holds exactly one FITS file."""
    p = Path(path)
    if p.is_dir():
        cands = [f for f in p.iterdir() if f.suffix.lower() in (".fit", ".fits", ".fts")]
        if len(cands) != 1:
            raise FileNotFoundError(f"Expected exactly one FITS file in {p}, found {len(cands)}: "
                                    f"{[c.name for c in cands]}")
        p = cands[0]
    return p


def load_master_flat():
    p = resolve_fits(MASTER_FLAT)
    flat, hdr = read_fits(p)
    good = np.isfinite(flat) & (flat > 0)
    return flat / flat[good].mean(), p, str(hdr.get("FILTER", "?")).strip()


def dark_rate_dn(T_celsius):
    """Dark rate [DN / unbinned pixel / s] and its relative uncertainty for the chosen DARK_MODE."""
    if DARK_MODE == "none":
        return 0.0, 0.0
    if DARK_MODE == "static":
        return DARK_7C / GAIN, np.hypot(DARK_7C_ERR / DARK_7C, GAIN_ERR / GAIN)
    T = T_celsius + 273.15
    return I0 * np.exp(-TB / T), np.hypot(I0_ERR / I0, TB_ERR / T)


def robust_sigma(a):
    a = a[np.isfinite(a)]
    return 1.4826 * np.median(np.abs(a - np.median(a)))


def load_master_darks(folder):
    if folder is None:
        return {}
    groups = {}
    for p in find_fits(folder):
        d, h = read_fits(p)
        groups.setdefault(round(get(h, EXPTIME_KEYS, -1.0), 3), []).append(d)
    return {t: np.median(np.stack(v), axis=0) for t, v in groups.items()}


ORIENTATIONS = {"as is": lambda a: a, "flip up-down": np.flipud,
                "flip left-right": np.fliplr, "rotate 180": lambda a: a[::-1, ::-1]}


def pearson(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m] - a[m].mean(), b[m] - b[m].mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else 0.0


def orient(name, frame, reference, choice, margin=0.02, min_corr=0.0):
    if choice != "auto":
        print(f"{name:12s} orientation: {choice} (set by hand)")
        return ORIENTATIONS[choice](frame)
    fill = lambda a: np.where(np.isfinite(a), a, np.nanmedian(a))
    smooth = lambda a: ndimage.median_filter(fill(a), size=9, mode="nearest")
    ref, frm = smooth(reference), smooth(frame)
    scores = {k: pearson(f(frm), ref) for k, f in ORIENTATIONS.items()}
    best = max(scores, key=scores.get)
    if scores[best] - scores["as is"] <= margin or scores[best] < min_corr:
        best = "as is"
    print(f"{name:12s} orientation: {best:16s} (correlation with star images: "
          + ", ".join(f"{k} {v:+.3f}" for k, v in scores.items()) + ")")
    return ORIENTATIONS[best](frame)


def hot_candidates(pre):
    resid = pre - ndimage.median_filter(pre, size=7, mode="nearest")
    sig = robust_sigma(resid)
    nb_x = np.maximum(np.roll(resid, 1, axis=1), np.roll(resid, -1, axis=1))
    nb_y = np.maximum(np.roll(resid, 1, axis=0), np.roll(resid, -1, axis=0))
    return (resid > HOT_SIGMA * sig) & (np.minimum(nb_x, nb_y) < SHARP_FRACTION * resid), sig


def trail_axis(img, hot):
    resid = img - ndimage.median_filter(img, size=7, mode="nearest")
    ny, nx = img.shape
    v = h = 0.0
    for y, x in zip(*np.nonzero(hot)):
        v += sum(resid[y + d, x] for d in range(-8, 9) if d and 0 <= y + d < ny)
        h += sum(resid[y, x + d] for d in range(-8, 9) if d and 0 <= x + d < nx)
    return 0 if v >= h else 1


def repair(pre, hot, sig, axis=0):
    if axis == 1:
        return repair(pre.T, hot.T, sig, 0).T
    out = pre.copy()
    ny, nx = pre.shape
    ys, xs = np.nonzero(hot)
    for y, x in zip(ys, xs):
        nb = [c for c in (x - 3, x - 2, x + 2, x + 3) if 0 <= c < nx and not hot[y, c]]
        out[y, x] = np.median(pre[y, nb]) if nb else np.median(pre[max(y - 3, 0):y + 4, max(x - 3, 0):x + 4])
    block = {c for x in np.unique(xs) for c in (x - 1, x, x + 1) if 0 <= c < nx}
    for xc in sorted(block):
        refs = [c for c in (xc - 4, xc - 3, xc + 3, xc + 4) if 0 <= c < nx and c not in block]
        if len(refs) < 2:
            continue
        excess = out[:, xc] - np.median(out[:, refs], axis=1)
        trail = ndimage.median_filter(excess, size=9, mode="nearest")
        out[:, xc] -= np.where(trail > TRAIL_SIGMA * sig, trail, 0.0)
    return out


def repair_flat(flat, hot, axis):
    f, h = (flat, hot) if axis == 0 else (flat.T, hot.T)
    f = f.copy()
    nx = f.shape[1]
    ys, xs = np.nonzero(h)
    before = f[ys, xs].copy()
    block = {c for x in np.unique(xs) for c in (x - 1, x, x + 1) if 0 <= c < nx}
    for xc in sorted(block):
        refs = [c for c in (xc - 4, xc - 3, xc + 3, xc + 4) if 0 <= c < nx and c not in block]
        if len(refs) >= 2:
            f[:, xc] = np.nanmedian(f[:, refs], axis=1)
    ratio = before / f[ys, xs]
    print(f"master flat at the hot pixels / neighbours: median {np.nanmedian(ratio):.3f}, "
          f"lowest {np.nanmin(ratio):.3f}, highest {np.nanmax(ratio):.3f}")
    return f if axis == 0 else f.T


def calibrate_folder(folder):
    folder = Path(folder)
    files = find_fits(folder)
    out_dir = folder / f"{folder.name}_calibrated"
    out_dir.mkdir(exist_ok=True)
    flat, flat_path, flat_filter = load_master_flat()
    bias_path = resolve_fits(MASTER_BIAS)
    bias_full, _ = read_fits(bias_path)
    master_darks = load_master_darks(DARK_FRAMES_FOLDER)
    print(f"{len(files)} images in {folder}\n"
          f"master flat : {flat_path.name} {flat.shape}\n"
          f"master bias : {bias_path.name} {bias_full.shape}, median {np.median(bias_full):.1f} DN\n"
          f"master darks: {sorted(master_darks) if master_darks else 'none (' + DARK_MODE + ' dark model)'}")

    # filter check: every science frame should match the flat's filter
    filters = {str(read_fits(p)[1].get("FILTER", "?")).strip() for p in files}
    print(f"filters     : science {sorted(filters)}, flat '{flat_filter}'")
    if flat_filter != "?" and any(f != flat_filter for f in filters):
        print("  (!) some science frames are in a different filter than the master flat")

    # ---------------- orientation of flat and bias ----------------
    sample = files[::max(1, len(files) // 25)]
    raw_med = np.median(np.stack([read_fits(p)[0].astype(np.float32) for p in sample]), axis=0).astype(np.float64)
    if raw_med.shape != flat.shape:
        raise ValueError(f"star images are {raw_med.shape} but the master flat is {flat.shape}: same binning needed")
    fy, fx = bias_full.shape[0] // flat.shape[0], bias_full.shape[1] // flat.shape[1]
    if fy < 1 or fx < 1 or (bias_full.shape[0] // fy, bias_full.shape[1] // fx) != flat.shape:
        raise ValueError(f"Cannot match the master bias {bias_full.shape} to the images {flat.shape}")
    bias_binned = rebin_mean(bias_full, fy, fx) if (fy, fx) != (1, 1) else bias_full
    flat_raw = flat
    flat = orient("master flat", flat_raw, raw_med - bias_binned, FLAT_ORIENT)
    g = np.isfinite(flat) & (flat >= FLAT_MIN)
    sky_scale = np.linalg.lstsq(np.c_[np.ones(g.sum()), flat[g]], raw_med[g], rcond=None)[0][1]
    bias_ref = raw_med - sky_scale * np.nan_to_num(flat, nan=1.0)
    bias_oriented = orient("master bias", bias_binned, bias_ref, BIAS_ORIENT, margin=0.1, min_corr=0.5)
    if bias_oriented is not bias_binned and FLAT_ORIENT == "auto":
        flat = orient("master flat", flat_raw, raw_med - bias_oriented, FLAT_ORIENT)


    bias_cache, meta, pres = {}, [], []
    hot_count = np.zeros(flat.shape, int)
    for path in files:
        raw, hdr = read_fits(path)
        if raw.shape != flat.shape:
            raise ValueError(f"{path.name} is {raw.shape} but the master flat is {flat.shape}: same binning needed")
        xb, yb = int(get(hdr, ["XBINNING"], 1)), int(get(hdr, ["YBINNING"], 1))
        t, T = get(hdr, EXPTIME_KEYS), get(hdr, TEMP_KEYS)
        ped = abs(float(hdr.get("PEDESTAL", 0.0))) if SUBTRACT_PEDESTAL else 0.0
        t_key = round(t, 3) if t is not None else None
        if t_key in master_darks:
            pre = raw - master_darks[t_key]
            rate, rel, dark_dn, how = np.nan, 0.0, np.nan, "dark frame"
        else:
            if (xb, yb) not in bias_cache:
                b = bias_oriented
                bias_cache[(xb, yb)] = ndimage.median_filter(b, size=BIAS_SMOOTH, mode="nearest") if BIAS_SMOOTH > 1 else b
            bias = bias_cache[(xb, yb)]
            bm = np.median(bias)
            bias_ped = ped if abs(bm - REF_BIAS_LEVEL_RAW) < abs(bm - (REF_BIAS_LEVEL_RAW - ped)) else 0.0
            if DARK_MODE != "none" and (t is None or (T is None and DARK_MODE == "law")):
                print(f"  (!) {path.name}: exposure time or CCD temperature missing, no dark subtracted")
                rate, rel = 0.0, 0.0
            else:
                rate, rel = dark_rate_dn(T if T is not None else 0.0)
            dark_dn = rate * (t or 0.0) * xb * yb
            pre = (raw - ped) - (bias - bias_ped) - dark_dn
            rate = rate * GAIN * xb * yb
            how = f"bias + {DARK_MODE} dark"
        cand, sig = (hot_candidates(pre) if FIX_HOT_PIXELS else (np.zeros(raw.shape, bool), robust_sigma(pre)))
        hot_count += cand
        pres.append(pre.astype(np.float32))
        meta.append(dict(path=path, hdr=hdr, t=t, T=T, xb=xb, yb=yb, ped=ped, rate=rate, rel=rel, dark_dn=dark_dn,
                         how=how, raw_max=raw.max(), raw_med=np.median(raw), sig=sig))
        print(f"  pre-calibrated {path.name:26s} t={t} T={T}  {how}")


    n = len(files)
    hot = hot_count >= max(1, int(np.ceil(HOT_FRACTION * n))) if FIX_HOT_PIXELS else np.zeros(flat.shape, bool)
    hys, hxs = np.nonzero(hot)
    axis = 0
    if hot.any():
        axis = TRAIL_AXIS if TRAIL_AXIS != "auto" else \
            trail_axis(np.median(np.stack(pres[::max(1, n // 15)]), axis=0).astype(np.float64), hot)
        print(f"\nhot pixels fixed on the detector: {hot.sum()} at (x, y) 1-based "
              f"{list(zip((hxs + 1).tolist(), (hys + 1).tolist()))[:10]}"
              f"\ntrails run along numpy axis {axis} ({'up/down a column' if axis == 0 else 'along a row'})")
        flat = repair_flat(flat, hot, axis)
    else:
        print("\nno persistent hot pixels found")
    with open(out_dir / "hot_pixels.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_numpy_col", "y_numpy_row", "x_1based", "y_1based", "frames_flagged"])
        for y, x in zip(hys, hxs):
            w.writerow([x, y, x + 1, y + 1, hot_count[y, x]])

    
    cals = []
    for pre, m in zip(pres, meta):
        fixed = repair(pre.astype(np.float64), hot, m["sig"], axis) if hot.any() else pre.astype(np.float64)
        cals.append(np.where(flat >= FLAT_MIN, fixed / flat, np.nan).astype(np.float32))
    p_low = min(np.nanmin(c) for c in cals)
    offset = float(ADD_OFFSET) if ADD_OFFSET != "auto" else (float(np.ceil(-p_low / 10) * 10) if p_low < 0 else 0.0)
    print(f"lowest calibrated pixel {p_low:.1f} DN -> offset {offset:.0f} DN")

    rows = []
    for cal, m in zip(cals, meta):
        sky, sigma = np.nanmedian(cal), robust_sigma(cal)
        frac_neg = float(np.mean(cal[np.isfinite(cal)] < 0))
        out = (cal + offset).astype(np.float32)
        hdr = m["hdr"]
        for k in ("BZERO", "BSCALE", "BLANK"):
            hdr.remove(k, ignore_missing=True)
        hdr["EGAIN"] = (GAIN, "gain [e-/DN] (Exercise B)")
        hdr["RDNOISE"] = (round(RN_DN * GAIN, 3), "read noise [e-] (Exercise A)")
        hdr["CALOFFS"] = (offset, "constant added after calibration [DN]")
        if np.isfinite(m["rate"]):
            hdr["DARKRATE"] = (round(m["rate"], 6), "dark rate at CCD temp [e-/s per pixel]")
        hdr.add_history("Calibrated: (raw - pedestal - bias - dark), hot pixels repaired, / flat, + offset")
        hdr.add_history(f"{m['how']}; flat {flat_path.name}; bias {bias_path.name} (smoothed {BIAS_SMOOTH} px)")
        fits.PrimaryHDU(out, header=hdr).writeto(out_dir / f"{m['path'].stem}_calibrated.fit", overwrite=True)
        flag = "ADC SATURATED" if m["raw_max"] >= ADC_DN else ("near/over full well" if m["raw_max"] >= FULL_WELL_DN else "")
        rows.append(dict(file=m["path"].name, exptime_s=m["t"], ccd_temp_C=m["T"], xbin=m["xb"], ybin=m["yb"],
                         method=m["how"], dark_rate_e_per_s_pix=m["rate"], dark_sub_DN=m["dark_dn"],
                         sky_median_DN=sky, sky_sigma_DN=sigma, fraction_negative_before_offset=frac_neg,
                         min_before_offset_DN=np.nanmin(cal), max_DN=np.nanmax(cal), offset_added_DN=offset,
                         raw_max_DN=m["raw_max"], flag=flag))
        print(f"  {m['path'].name:26s} sky {sky:8.1f} +/- {sigma:5.1f} DN  max {np.nanmax(cal):8.0f}  "
              f"negative before offset: {100 * frac_neg:5.1f}%  {flag}")

    with open(out_dir / "calibration_log.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    rn_e = RN_DN * GAIN
    d_rn_e = rn_e * np.hypot(RN_ERR / RN_DN, GAIN_ERR / GAIN)
    dark_e = np.array([r["dark_rate_e_per_s_pix"] for r in rows], float)
    dark_txt = (f"{np.nanmean(dark_e):.4f} (range {np.nanmin(dark_e):.4f} to {np.nanmax(dark_e):.4f}; per binned pixel)"
                if np.isfinite(dark_e).any() else "0 (matching dark frames were subtracted)")
    summary = [
        f"Calibrated {len(rows)} images from {folder}", f"Output folder: {out_dir}", "",
        "Calibration: (raw - pedestal - bias - dark) -> hot pixels/trails repaired -> / flat -> + constant offset",
        f"  master flat: {flat_path}", f"  master bias: {bias_path} (large-scale structure only)",
        f"  dark: {'frame(s) from ' + str(DARK_FRAMES_FOLDER) if master_darks else DARK_MODE}",
        f"  hot pixels repaired: {int(hot.sum())} (list in hot_pixels.csv); offset added: {offset:.0f} DN (CALOFFS)",
        "", "Numbers for the AstroImageJ Multi-Aperture window:",
        f"  CCD gain (e-/ADU)          : {GAIN} +/- {GAIN_ERR}",
        f"  CCD readout noise (e-)     : {rn_e:.2f} +/- {d_rn_e:.2f}",
        f"  CCD dark current (e-/s/pix): {dark_txt}",
    ]
    (out_dir / "calibration_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n" + "\n".join(summary[9:]) + f"\nAll files written to {out_dir}")


if __name__ == "__main__":
    calibrate_folder(sys.argv[1] if len(sys.argv) > 1 else INPUT_FOLDER)
