#!/usr/bin/env python3
"""
Phase 35b — wind stress curl 鞍點：格點 mask corr + null + 相位骨架公平對照
============================================================================
補齊 Phase 32 方法學三教訓：binned Spearman 可能係 K=det(H) 同 saddle mask
數學耦合 → 需要格點 mask corr + null 對照；相位骨架需要公平對照（30v30 月 vs
隨機月 null）。
"""
import numpy as np
import xarray as xr
import json, os
from scipy import stats

OUT_DIR = '/tmp/enso_saddle_3d'
RESULTS_DIR = '/app/working/workspaces/tygtDc/projects/atmosphere/vorticity_dfold'
YEARS = {'2015': 'El Niño (strong)', '2016': 'El Niño→Neutral', '2020': 'La Niña'}
RHO_A = 1.225
CD = 1.3e-3


def load_u10(year):
    ds = xr.open_dataset(os.path.join(OUT_DIR, f'era5_u10_{year}.nc'))
    lat = ds['latitude'].values
    lon = ds['longitude'].values
    u = ds['u10'].values.squeeze()
    ds.close()
    return u, lat, lon


def load_v10(year):
    ds = xr.open_dataset(os.path.join(OUT_DIR, f'era5_v10_{year}.nc'))
    v = ds['v10'].values.squeeze()
    ds.close()
    return v


def wind_stress(u, v):
    speed = np.hypot(u, v)
    return RHO_A * CD * speed * u, RHO_A * CD * speed * v


def curl_sphere(tau_x, tau_y, lat, lon):
    R = 6371e3
    dphi = np.radians(np.abs(lat[1] - lat[0]))
    dlam = np.radians(np.abs(lon[1] - lon[0]))
    phi = np.radians(lat)
    dtauy_dx = np.gradient(tau_y, dlam, axis=1) / (R * np.cos(phi[:, None]))
    dtautx_dy = np.gradient(tau_x, dphi, axis=0) / R
    return dtauy_dx - dtautx_dy


def hessian2d(f, spacing):
    dy, dx = spacing
    fy = np.gradient(f, dy, axis=0)
    fx = np.gradient(f, dx, axis=1)
    fyy = np.gradient(fy, dy, axis=0)
    fxx = np.gradient(fx, dx, axis=1)
    fxy = np.gradient(fx, dy, axis=0)
    return fxx, fyy, fxy


def morse_classify_2d(field, spacing, parabolic_frac=0.10):
    valid = ~np.isnan(field)
    fxx, fyy, fxy = hessian2d(field, spacing)
    det = fxx * fyy - fxy**2
    trace = fxx + fyy
    adet = np.abs(det)
    thr = np.nanpercentile(adet[valid], parabolic_frac * 100) if valid.any() else 0.0
    parabolic = (~valid) | (adet < thr)
    saddle = (~parabolic) & (det < 0)
    return {'saddle': saddle, 'valid': valid}


def saddle_density_map(saddle, valid, box=9):
    from scipy.ndimage import uniform_filter
    d = uniform_filter(saddle.astype(float), size=box, mode='constant')
    v = uniform_filter(valid.astype(float), size=box, mode='constant')
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(v > 0, d / np.maximum(v, 1e-9), np.nan)


def main():
    print("=" * 72)
    print("Phase 35b — curl τ 鞍點：格點 mask corr + null + 相位骨架公平對照")
    print("=" * 72)
    spacing = (0.5, 0.5)

    summary = {}
    fields = {}
    for year in YEARS:
        u, lat, lon = load_u10(year)
        v = load_v10(year)
        tx, ty = wind_stress(u, v)
        ct = curl_sphere(tx, ty, lat, lon)
        ct = (ct - np.nanmean(ct)) / np.nanstd(ct)
        M = morse_classify_2d(ct, spacing)
        fields[year] = ct
        summary[year] = {'M': M}

    # ① 格點 mask corr（鞍點密度 × |curl τ|，逐格點）
    print("\n--- ① 格點 mask corr（saddle density × |curl τ|）---")
    for year in YEARS:
        ct = fields[year]
        M = summary[year]['M']
        rho = saddle_density_map(M['saddle'], M['valid'])
        mag = np.abs(ct)
        valid2 = M['valid'] & ~np.isnan(rho)
        r, p = stats.spearmanr(mag[valid2], rho[valid2])
        print(f"  {year}: grid Spearman r={r:.4f} (p={p:.2g}, n={valid2.sum():,})")
        summary[year]['grid_spearman'] = float(r)
        summary[year]['grid_spearman_p'] = float(p)

    # ② null 對照（置換場 mask corr）
    print("\n--- ② null 對照（置換場）---")
    for year in YEARS:
        ct = fields[year]
        M = summary[year]['M']
        valid = M['valid']
        rng = np.random.default_rng(int(year) + 1000)
        null_corrs = []
        for s in range(20):
            fs = ct.copy()
            fs[valid] = rng.permutation(ct[valid])
            Ms = morse_classify_2d(fs, spacing)
            rho = saddle_density_map(Ms['saddle'], Ms['valid'])
            valid2 = Ms['valid'] & ~np.isnan(rho)
            r, _ = stats.spearmanr(np.abs(fs)[valid2], rho[valid2])
            null_corrs.append(r)
        obs = summary[year]['grid_spearman']
        z = (obs - np.mean(null_corrs)) / np.std(null_corrs)
        print(f"  {year}: null {np.mean(null_corrs):.4f}±{np.std(null_corrs):.4f} → z={z:+.1f}")
        summary[year]['null_corr_mean'] = float(np.mean(null_corrs))
        summary[year]['null_corr_std'] = float(np.std(null_corrs))
        summary[year]['corr_z'] = float(z)

    # ③ 相位骨架公平對照：30 日骨架 corr vs 隨機日 null
    # （此處用年度快照做框架示範 — 真正需要逐日數據，但先跑年度對照）
    print("\n--- ③ 相位骨架對照（年度場 corr）---")
    s15, s16, s20 = fields['2015'], fields['2016'], fields['2020']
    M15, M16, M20 = summary['2015']['M'], summary['2016']['M'], summary['2020']['M']
    # 骨架 = saddle 密度 map
    rho15 = saddle_density_map(M15['saddle'], M15['valid'])
    rho16 = saddle_density_map(M16['saddle'], M16['valid'])
    rho20 = saddle_density_map(M20['saddle'], M20['valid'])
    common = M15['valid'] & M16['valid'] & M20['valid'] & ~np.isnan(rho15) & ~np.isnan(rho16) & ~np.isnan(rho20)
    c1516, _ = stats.spearmanr(rho15[common], rho16[common])
    c1520, _ = stats.spearmanr(rho15[common], rho20[common])
    c1620, _ = stats.spearmanr(rho16[common], rho20[common])
    print(f"  骨架 corr: 2015v2016 {c1516:.3f} | 2015v2020 {c1520:.3f} | 2016v2020 {c1620:.3f}")
    # null：隨機日骨架（用置換場）
    rng = np.random.default_rng(42)
    null_corrs = []
    for s in range(30):
        fs = s15.copy(); fs[M15['valid']] = rng.permutation(s15[M15['valid']])
        Ms = morse_classify_2d(fs, spacing)
        rs = saddle_density_map(Ms['saddle'], Ms['valid'])
        fs2 = s16.copy(); fs2[M16['valid']] = rng.permutation(s16[M16['valid']])
        Ms2 = morse_classify_2d(fs2, spacing)
        rs2 = saddle_density_map(Ms2['saddle'], Ms2['valid'])
        m = M15['valid'] & M16['valid'] & ~np.isnan(rs) & ~np.isnan(rs2)
        r, _ = stats.spearmanr(rs[m], rs2[m])
        null_corrs.append(r)
    z15 = (c1516 - np.mean(null_corrs)) / np.std(null_corrs)
    print(f"  2015v2016 null: {np.mean(null_corrs):.3f}±{np.std(null_corrs):.3f} → z={z15:+.1f}")
    summary['skeleton'] = {'c1516': float(c1516), 'c1520': float(c1520), 'c1620': float(c1620),
                           'null_mean': float(np.mean(null_corrs)), 'null_std': float(np.std(null_corrs)),
                           'z_1516': float(z15)}

    # 存
    out = {}
    for y in YEARS:
        out[y] = {k: v for k, v in summary[y].items() if k != 'M'}
    out['skeleton'] = summary['skeleton']
    out_path = os.path.join(RESULTS_DIR, 'phase35b_enso_curl_corr_null.json')
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n  ✅ Saved: {out_path}")


if __name__ == '__main__':
    main()
