#!/usr/bin/env python3
"""
ENSO 3D Saddle Ring — Hessian/Morse 3D 跨域測試（第二域 3D 版）
=============================================================
MKP「之前enso好像也用了2d測，試試再3d測返」→ 同 Greenland 3D 平行，
用熱帶太平洋 ERA5 8 層壓力層 geopotential + 渦度，3D Hessian/Morse，
測 ENSO 相位（2015 El Niño / 2016 轉中性 / 2020 La Niña）鞍點環結構。

對照：
- ENSO 2D（Phase 30-32）：場級鞍點比例 trivial（≈0.564，相位唔敏感）
- Greenland 3D（Phase 33）：冇鞍點環（逐月 sign test 拆穿季節 artifacts）
- 颱風 2D/3D：眼牆鞍點環普適（saddle>55% + Rayleigh<0.08 + 時間穩定）
"""
import numpy as np
import xarray as xr
import json, os, time
from scipy import stats

OUT_DIR = '/tmp/enso_saddle_3d'
RESULTS_DIR = '/app/working/workspaces/tygtDc/projects/atmosphere/vorticity_dfold'
YEARS = {'2015': 'El Niño (strong)', '2016': 'El Niño→Neutral', '2020': 'La Niña'}
# ENSO 候選中心：Nino3.4 區域 (-5..5N, 190..240E)、東端 (0, 280E)、暖池 (0, 160E)
CENTERS = {'nino34': (0.0, 210.0), 'east_end': (0.0, 280.0), 'warm_pool': (0.0, 160.0)}


def hessian_3d(f, spacing):
    dz, dy, dx = spacing
    gz = np.gradient(f, dz, axis=0)
    gy = np.gradient(f, dy, axis=1)
    gx = np.gradient(f, dx, axis=2)
    fzz = np.gradient(gz, dz, axis=0)
    fyy = np.gradient(gy, dy, axis=1)
    fxx = np.gradient(gx, dx, axis=2)
    fxz = (np.gradient(gx, dz, axis=0) + np.gradient(gz, dx, axis=2)) / 2.0
    fyz = (np.gradient(gy, dz, axis=0) + np.gradient(gz, dy, axis=1)) / 2.0
    fxy = (np.gradient(gx, dy, axis=1) + np.gradient(gy, dx, axis=2)) / 2.0
    return fxx, fyy, fzz, fxy, fxz, fyz


def eigvals_sym3(a, b, c, d, e, f):
    p1 = a + b + c
    p2 = a*b + b*c + a*c - (d**2 + e**2 + f**2)
    p3 = a*b*c + 2*d*e*f - a*f**2 - b*e**2 - c*d**2
    Q = np.maximum((p1**2 - 3*p2) / 9.0, 0)
    R = (2*p1**3 - 9*p1*p2 + 27*p3) / 54.0
    theta = np.arccos(np.clip(R / np.sqrt(Q**3 + 1e-300), -1, 1))
    sqrtQ = np.sqrt(Q)
    l1 = 2*sqrtQ*np.cos(theta/3) + p1/3.0
    l2 = 2*sqrtQ*np.cos((theta - 2*np.pi)/3) + p1/3.0
    l3 = 2*sqrtQ*np.cos((theta + 2*np.pi)/3) + p1/3.0
    l = np.stack([l1, l2, l3], axis=-1)
    l.sort(axis=-1)
    return l[..., 0], l[..., 1], l[..., 2]


def morse_classify_3d(field, spacing, parabolic_frac=0.10):
    valid = ~np.isnan(field)
    fxx, fyy, fzz, fxy, fxz, fyz = hessian_3d(field, spacing)
    l1, l2, l3 = eigvals_sym3(fxx, fyy, fzz, fxy, fxz, fyz)
    det = l1 * l2 * l3
    adet = np.abs(det)
    thr = np.nanpercentile(adet[valid], parabolic_frac * 100) if valid.any() else 0.0
    parabolic = (~valid) | (adet < thr)
    n_neg = (l1 < 0).astype(int) + (l2 < 0).astype(int) + (l3 < 0).astype(int)
    return {'idx0_peak': (~parabolic) & (n_neg == 3),
            'idx1_saddle': (~parabolic) & (n_neg == 2),
            'idx2_saddle': (~parabolic) & (n_neg == 1),
            'idx3_trough': (~parabolic) & (n_neg == 0),
            'parabolic': parabolic, 'valid': valid}


def load_field(var, year, field_type='z'):
    """field_type='z' → geopotential (m); 'vort' → relative vorticity"""
    path = os.path.join(OUT_DIR, f'era5_{var}_{year}_8lev.nc')
    ds = xr.open_dataset(path)
    lat = ds['latitude'].values
    lon = ds['longitude'].values
    lev = ds['pressure_level'].values
    if field_type == 'z':
        f = ds['z'].values.mean(axis=0) / 9.80665
    else:
        f = ds[var].values.mean(axis=0)
    ds.close()
    return f, lat, lon, lev


def load_vorticity(year):
    """相對渦度 ζ = dv/dx - du/dy（球面有限差分）"""
    u, lat, lon, lev = load_field('u', year, 'u')
    v, _, _, _ = load_field('v', year, 'v')
    R = 6371e3
    dphi = np.radians(np.abs(lat[1] - lat[0]))
    dlam = np.radians(np.abs(lon[1] - lon[0]))
    phi = np.radians(lat)
    du = np.gradient(u, dphi, axis=1)
    dv = np.gradient(v, dlam, axis=2)
    zeta = (1.0 / (R * np.cos(phi[None, :, None]))) * dv - (1.0 / R) * du
    # 每層標準化（垂直可比）
    zeta = (zeta - zeta.mean(axis=(1, 2), keepdims=True)) / zeta.std(axis=(1, 2), keepdims=True)
    return zeta, lat, lon, lev


def band_stats(field, M, lat, lon, clat, clon, spacing):
    """3D 環帶統計：core<1.5° / ring 1.5-3° / inner 3-6° / outer 6-12°"""
    dist2d = np.sqrt(np.abs(lat - clat)[:, None]**2 +
                     np.abs(((lon - clon + 180) % 360) - 180)[None, :]**2)
    valid = M['valid']
    saddle = M['idx1_saddle'] | M['idx2_saddle']
    out = {}
    for name, r0, r1 in [('core', 0, 1.5), ('ring', 1.5, 3.0), ('inner', 3.0, 6.0), ('outer', 6.0, 12.0)]:
        mask3d = np.broadcast_to((dist2d >= r0) & (dist2d < r1), field.shape) & valid
        n = mask3d.sum()
        ns = (mask3d & saddle).sum()
        frac = ns / n if n else np.nan
        zs, ys, xs = np.where(mask3d & saddle)
        if len(ys) > 30:
            az = np.degrees(np.arctan2(lat[ys] - clat,
                                       ((lon[xs] - clon + 180) % 360) - 180)) % 360
            rv = np.abs(np.mean(np.exp(1j * np.radians(az))))
        else:
            rv = None
        out[name] = {'n': int(n), 'n_saddle': int(ns), 'saddle_frac': float(frac),
                     'rayleigh_r': rv}
    return out


def main():
    t0 = time.time()
    print("=" * 72)
    print("ENSO 3D Saddle Ring — Hessian/Morse 3D 跨域測試（Phase 34）")
    print("=" * 72)
    dlat = dlon = 0.5
    logp = np.log10(np.array([850, 700, 600, 500, 400, 300, 250, 200]))
    dlogp = np.abs(np.diff(logp).mean())
    spacing = (dlogp, dlat, dlon)

    summary = {}
    for year, phase in YEARS.items():
        print(f"\n{'='*60}\n--- {year} ({phase}) ---\n{'='*60}")
        # geopotential
        fz, lat, lon, lev = load_field('z', year, 'z')
        Mz = morse_classify_3d(fz, spacing)
        nz_, ny_, nx_ = fz.shape
        tot = Mz['valid'].sum()
        sad_z = (Mz['idx1_saddle'].sum() + Mz['idx2_saddle'].sum()) / tot
        # 渦度
        fv, _, _, _ = load_vorticity(year)
        Mv = morse_classify_3d(fv, spacing)
        sad_v = (Mv['idx1_saddle'].sum() + Mv['idx2_saddle'].sum()) / tot
        print(f"  3D volume: {nz_} levels x {ny_} lat x {nx_} lon = {fz.size:,} pts")
        print(f"  全場鞍點: geopotential {sad_z*100:.2f}% | vorticity {sad_v*100:.2f}%")
        # null
        rng = np.random.default_rng(int(year))
        nulls_z, nulls_v = [], []
        for s in range(20):
            fs = fz.copy(); fs[Mz['valid']] = rng.permutation(fz[Mz['valid']])
            Ms = morse_classify_3d(fs, spacing)
            nulls_z.append((Ms['idx1_saddle'].sum() + Ms['idx2_saddle'].sum()) / tot)
            fs = fv.copy(); fs[Mv['valid']] = rng.permutation(fv[Mv['valid']])
            Ms = morse_classify_3d(fs, spacing)
            nulls_v.append((Ms['idx1_saddle'].sum() + Ms['idx2_saddle'].sum()) / tot)
        zz = (sad_z - np.mean(nulls_z)) / np.std(nulls_z)
        zv = (sad_v - np.mean(nulls_v)) / np.std(nulls_v)
        print(f"  null: z {np.mean(nulls_z)*100:.2f}%±{np.std(nulls_z)*100:.2f}% (z={zz:+.1f}) | "
              f"vort {np.mean(nulls_v)*100:.2f}%±{np.std(nulls_v)*100:.2f}% (z={zv:+.1f})")

        yr = {'phase': phase, 'n': int(nz_ * ny_ * nx_),
              'saddle_geo_frac': float(sad_z), 'saddle_vort_frac': float(sad_v),
              'null_geo_mean': float(np.mean(nulls_z)), 'null_geo_std': float(np.std(nulls_z)),
              'null_vort_mean': float(np.mean(nulls_v)), 'null_vort_std': float(np.std(nulls_v)),
              'z_geo': float(zz), 'z_vort': float(zv)}
        # 環帶（渦度場 — 同颱風域可比）
        print("  --- 渦度場 3D 環帶（ENSO 候選中心）---")
        for cname, (clat, clon) in CENTERS.items():
            bs = band_stats(fv, Mv, lat, lon, clat, clon, spacing)
            ring = bs['ring']['saddle_frac']; core = bs['core']['saddle_frac']
            rr = bs['ring']['rayleigh_r']
            print(f"  {cname:<10s} ({clat:.0f},{clon:.0f}): core={core*100:5.1f}%  "
                  f"ring={ring*100:5.1f}% (Δ{100*(ring-core):+5.1f}pp)  R={rr if rr is None else f'{rr:.3f}'}")
            yr[cname] = bs
        summary[year] = yr

    # 相位對比：全場鞍點比例
    print("\n" + "=" * 72)
    print("相位對比（全場鞍點比例）")
    print("=" * 72)
    geo = [summary[y]['saddle_geo_frac'] for y in YEARS]
    vort = [summary[y]['saddle_vort_frac'] for y in YEARS]
    print(f"  geopotential: {['%.3f' % g for g in geo]}")
    print(f"  vorticity:    {['%.3f' % v for v in vort]}")

    out_path = os.path.join(RESULTS_DIR, 'phase34_enso_saddle_ring_3d_results.json')
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  ✅ Saved: {out_path}")
    print(f"  Completed in {time.time() - t0:.1f}s")


if __name__ == '__main__':
    main()
