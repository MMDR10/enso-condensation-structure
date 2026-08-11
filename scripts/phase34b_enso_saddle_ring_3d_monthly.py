#!/usr/bin/env python3
"""
ENSO 3D Saddle Ring — 逐月穩定性測試（nino34 中心）
=====================================================
Greenland 教訓：年平均 ring>core 可能係季節混合 artifacts，
必須拆逐月 sign test 驗證。2015/2016/2020 逐月 ring-core diff。
"""
import numpy as np
import xarray as xr
import json, os
from scipy import stats as ss

OUT_DIR = '/tmp/enso_saddle_3d'
RESULTS_DIR = '/app/working/workspaces/tygtDc/projects/atmosphere/vorticity_dfold'
YEARS = ['2015', '2016', '2020']
CLAT, CLON = 0.0, 210.0  # nino34


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
    return {'saddle': (~parabolic) & ((n_neg == 2) | (n_neg == 1)),
            'parabolic': parabolic, 'valid': valid}


def load_vorticity_month(year, month_idx):
    """單月 3D 渦度場（標準化）"""
    R = 6371e3
    fields = {}
    for var in ('u', 'v'):
        path = os.path.join(OUT_DIR, f'era5_{var}_{year}_8lev.nc')
        ds = xr.open_dataset(path)
        fields[var] = ds[var].values[month_idx]
        if var == 'u':
            lat = ds['latitude'].values
            lon = ds['longitude'].values
            lev = ds['pressure_level'].values
        ds.close()
    dphi = np.radians(np.abs(lat[1] - lat[0]))
    dlam = np.radians(np.abs(lon[1] - lon[0]))
    phi = np.radians(lat)
    du = np.gradient(fields['u'], dphi, axis=1)
    dv = np.gradient(fields['v'], dlam, axis=2)
    zeta = (1.0 / (R * np.cos(phi[None, :, None]))) * dv - (1.0 / R) * du
    zeta = (zeta - zeta.mean(axis=(1, 2), keepdims=True)) / zeta.std(axis=(1, 2), keepdims=True)
    return zeta, lat, lon, lev


def main():
    print("=" * 72)
    print("ENSO 3D nino34 環帶 — 逐月穩定性（Greenland 教訓）")
    print("=" * 72)
    dlat = dlon = 0.5
    logp = np.log10(np.array([850, 700, 600, 500, 400, 300, 250, 200]))
    dlogp = np.abs(np.diff(logp).mean())
    spacing = (dlogp, dlat, dlon)

    results = {}
    for year in YEARS:
        print(f"\n--- {year} 逐月 ring-core（nino34 中心）---")
        diffs = []
        for m in range(12):
            fv, lat, lon, lev = load_vorticity_month(year, m)
            M = morse_classify_3d(fv, spacing)
            dist2d = np.sqrt(np.abs(lat - CLAT)[:, None]**2 +
                             np.abs(((lon - CLON + 180) % 360) - 180)[None, :]**2)
            valid = M['valid']
            core = np.broadcast_to(dist2d < 1.5, fv.shape) & valid
            ring = np.broadcast_to((dist2d >= 1.5) & (dist2d < 3.0), fv.shape) & valid
            fc = (core & M['saddle']).sum() / core.sum()
            fr = (ring & M['saddle']).sum() / ring.sum()
            diffs.append(fr - fc)
            print(f"  month {m+1:2d}: ring-core={100*(fr-fc):+6.2f}pp  "
                  f"(ring {fr*100:.1f}% / core {fc*100:.1f}%)")
        d = np.array(diffs)
        npos = (d > 0).sum()
        p = ss.binomtest(npos, 12).pvalue
        print(f"  → 正值 {npos}/12, mean {d.mean()*100:+.2f}pp, sign test p={p:.4f}")
        results[year] = {'monthly_ring_minus_core_pp': [float(x*100) for x in diffs],
                         'n_pos': int(npos), 'mean_pp': float(d.mean()*100),
                         'sign_test_p': float(p)}

    # 三年合併 sign test
    all_d = np.concatenate([np.array(results[y]['monthly_ring_minus_core_pp']) for y in YEARS])
    npos_all = (all_d > 0).sum()
    p_all = ss.binomtest(npos_all, len(all_d)).pvalue
    print(f"\n  三年合併: 正值 {npos_all}/36, mean {all_d.mean():+.2f}pp, sign test p={p_all:.6f}")
    results['_combined'] = {'n_pos': int(npos_all), 'n': int(len(all_d)),
                            'mean_pp': float(all_d.mean()), 'sign_test_p': float(p_all)}

    out = os.path.join(RESULTS_DIR, 'phase34_enso_saddle_ring_3d_monthly.json')
    with open(out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  ✅ Saved: {out}")


if __name__ == '__main__':
    main()
