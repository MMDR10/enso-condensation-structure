#!/usr/bin/env python3
"""
Phase 35 — ENSO 氣流-海面交界面（wind stress curl）鞍點測試
============================================================
MKP insight：「如果說 ENSO 的鞍點在哪，我猜在氣流和海面之間」

測試：10m 風 → bulk wind stress τ = ρ_a C_D |U| U → curl τ（大氣驅動海洋
嘅旋度注入，正正係「氣流和海面之間」）。用 2D Hessian/Morse（同 Phase 30-32
ENSO 2D 一致），測：
① 全場鞍點比例 × ENSO 相位敏感度
② 鞍點密度 × |curl τ| 空間相關（鞍點係咪聚集喺氣流-海面耦合帶）
③ 相位骨架對照（同 Phase 32 公平對照）
"""
import numpy as np
import xarray as xr
import json, os, time
from scipy import stats

OUT_DIR = '/tmp/enso_saddle_3d'
RESULTS_DIR = '/app/working/workspaces/tygtDc/projects/atmosphere/vorticity_dfold'
YEARS = {'2015': 'El Niño (strong)', '2016': 'El Niño→Neutral', '2020': 'La Niña'}
RHO_A = 1.225   # kg/m^3
CD = 1.3e-3     # bulk drag coefficient


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
    """Bulk formula: τ = ρ_a C_D |U| U"""
    speed = np.hypot(u, v)
    tau_x = RHO_A * CD * speed * u
    tau_y = RHO_A * CD * speed * v
    return tau_x, tau_y


def curl_sphere(tau_x, tau_y, lat, lon):
    """球面 curl τ = ∂τy/∂x - ∂τx/∂y（單位 Pa/m，標準化後用）"""
    R = 6371e3
    dphi = np.radians(np.abs(lat[1] - lat[0]))
    dlam = np.radians(np.abs(lon[1] - lon[0]))
    phi = np.radians(lat)
    dtauy_dx = np.gradient(tau_y, dlam, axis=1) / (R * np.cos(phi[:, None]))
    dtautx_dy = np.gradient(tau_x, dphi, axis=0) / R
    c = dtauy_dx - dtautx_dy
    return c


def hessian2d(f, spacing):
    dy, dx = spacing
    fy = np.gradient(f, dy, axis=0)
    fx = np.gradient(f, dx, axis=1)
    fyy = np.gradient(fy, dy, axis=0)
    fxx = np.gradient(fx, dx, axis=1)
    # cross: ∂²f/∂x∂y
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
    # saddle: det<0（特徵值異號）
    saddle = (~parabolic) & (det < 0)
    peak = (~parabolic) & (det > 0) & (trace < 0)
    trough = (~parabolic) & (det > 0) & (trace > 0)
    return {'saddle': saddle, 'peak': peak, 'trough': trough,
            'parabolic': parabolic, 'valid': valid}


def saddle_density_map(saddle, valid, box=9):
    """局部鞍點密度（box×box 滑窗）"""
    from scipy.ndimage import uniform_filter
    d = uniform_filter(saddle.astype(float), size=box, mode='constant')
    v = uniform_filter(valid.astype(float), size=box, mode='constant')
    with np.errstate(divide='ignore', invalid='ignore'):
        rho = np.where(v > 0, d / np.maximum(v, 1e-9), np.nan)
    return rho


def main():
    t0 = time.time()
    print("=" * 72)
    print("Phase 35 — ENSO 氣流-海面交界面（wind stress curl）鞍點測試")
    print("MKP insight: 鞍點在氣流和海面之間")
    print("=" * 72)
    dlat = dlon = 0.5
    spacing = (dlat, dlon)

    summary = {}
    for year, phase in YEARS.items():
        print(f"\n{'='*60}\n--- {year} ({phase}) ---\n{'='*60}")
        u, lat, lon = load_u10(year)
        v = load_v10(year)
        tx, ty = wind_stress(u, v)
        curl_t = curl_sphere(tx, ty, lat, lon)
        # 標準化
        curl_t = (curl_t - np.nanmean(curl_t)) / np.nanstd(curl_t)
        M = morse_classify_2d(curl_t, spacing)
        tot = M['valid'].sum()
        sad = M['saddle'].sum() / tot
        pk = M['peak'].sum() / tot
        tr = M['trough'].sum() / tot
        print(f"  grid: {len(lat)}x{len(lon)} = {tot:,.0f} valid pts")
        print(f"  curl τ: saddle {sad*100:.2f}% | peak {pk*100:.2f}% | trough {tr*100:.2f}%")

        # null（置換場）
        rng = np.random.default_rng(int(year))
        nulls = []
        for s in range(20):
            fs = curl_t.copy()
            fs[M['valid']] = rng.permutation(curl_t[M['valid']])
            Ms = morse_classify_2d(fs, spacing)
            nulls.append(Ms['saddle'].sum() / tot)
        zz = (sad - np.mean(nulls)) / np.std(nulls)
        print(f"  null: {np.mean(nulls)*100:.2f}%±{np.std(nulls)*100:.2f}% (z={zz:+.1f})")

        # 鞍點密度 × |curl τ|（同 Phase 31 ρ 場方法）
        rho = saddle_density_map(M['saddle'], M['valid'])
        mag = np.abs(curl_t)
        valid2 = M['valid'] & ~np.isnan(rho)
        # binned Spearman
        bin_edges = np.nanpercentile(mag[valid2], np.linspace(0, 100, 11))
        bin_idx = np.digitize(mag, bin_edges) - 1
        bin_idx = np.clip(bin_idx, 0, 9)
        bin_rho = np.full(10, np.nan)
        for b in range(10):
            m = valid2 & (bin_idx == b)
            if m.sum() > 50:
                bin_rho[b] = rho[m].mean()
        ok = ~np.isnan(bin_rho)
        if ok.sum() >= 4:
            rho_spearman, p_sp = stats.spearmanr(bin_edges[:-1][ok], bin_rho[ok])
        else:
            rho_spearman, p_sp = np.nan, np.nan
        print(f"  鞍點密度 × |curl τ| binned Spearman: {rho_spearman:.3f} (p={p_sp:.2g})")

        yr = {'phase': phase,
              'saddle_frac': float(sad), 'peak_frac': float(pk), 'trough_frac': float(tr),
              'null_mean': float(np.mean(nulls)), 'null_std': float(np.std(nulls)), 'z': float(zz),
              'spearman_density_curl': float(rho_spearman)}
        summary[year] = yr

    # 相位對比
    print("\n" + "=" * 72)
    print("相位對比（curl τ 全場鞍點比例）")
    print("=" * 72)
    fracs = [summary[y]['saddle_frac'] for y in YEARS]
    print(f"  {'  '.join(YEARS)}")
    print(f"  {'  '.join('%.3f' % f for f in fracs)}")
    # KW test
    # 逐月無月份數據（此處用三年三個值，只做描述性）
    print(f"  range: {max(fracs)-min(fracs)*100:.2f}pp")

    out_path = os.path.join(RESULTS_DIR, 'phase35_enso_windstress_curl_saddle.json')
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  ✅ Saved: {out_path}")
    print(f"  Completed in {time.time() - t0:.1f}s")


if __name__ == '__main__':
    main()
