#!/usr/bin/env python3
"""
ENSO Saddle Ring 域測量 — Saddle Ring 算子跨域驗證（Phase 1）
================================================================
目的：將颱風域開發嘅 Saddle Ring 算子（Hessian Morse 分類 → 鞍點比例）
套到 ENSO 域（OISST 0.25° SST anomaly 場），驗證算子跨域適用性。

物理問題：
  颱風鞍點 = 眼牆剪切帶（cyclonic/anticyclonic 交匯）
  ENSO 鞍點 = SST anomaly 場 Hessian 異號點 = 冷暖 anomaly 交匯鋒面帶
  1. ENSO 場鞍點比例係咪 ENSO 相位依賴？（El Niño / La Niña / Neutral）
  2. 鞍點空間分佈同凝結結構線（D_fold 凝結核心/鋒面）有冇共位？
  3. 同颱風域基準比較：鞍點比例量級（颱風 55-67% vs ENSO ?）

方法：
  - OISST 0.25° 528 月（1982-2025）SST anomaly（全期月氣候減）
  - Hessian 特徵值符號 → Morse 分類：peak/trough/saddle/parabolic
  - 場級鞍點比例 + null（隨機場值 shuffle 20 surrogates）
  - 相位對照：ONI ≥+0.5 El Niño / ≤−0.5 La Niña / 中間 Neutral
  - 鞍點密度空間 map（climatology）→ 同凝結核心（|K|>p99）共位 Jaccard
"""
import numpy as np
import json
from pathlib import Path

OUT = Path('/app/working/workspaces/tygtDc/projects/atmosphere/vorticity_dfold/results')
SST1 = '/app/working/workspaces/tygtDc/data/sst/oisst_mon_1982_2020.npy'
SST2 = '/app/working/workspaces/tygtDc/data/sst/oisst_mon_2021_2025.npy'
LAT = '/app/working/workspaces/tygtDc/data/sst/oisst_mon_lat.npy'
LON = '/app/working/workspaces/tygtDc/data/sst/oisst_mon_lon.npy'
ONI = '/app/working/workspaces/tygtDc/projects/enso/release/oni.csv'


def load_oni():
    """Parse ONI csv → dict {(year, month): value} for 1982-2025."""
    oni = {}
    with open(ONI) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 13:
                continue
            try:
                yr = int(parts[0])
            except ValueError:
                continue
            if not (1982 <= yr <= 2025):
                continue
            vals = [float(x) for x in parts[1:13]]
            for m, v in enumerate(vals, 1):
                oni[(yr, m)] = v
    return oni


def hessian_2d(f, spacing=(1.0, 1.0)):
    """Hessian via 2nd-order central differences (numpy gradient)."""
    gy = np.gradient(f, spacing[0], axis=0)
    gx = np.gradient(f, spacing[1], axis=1)
    fyy = np.gradient(gy, spacing[0], axis=0)
    fxx = np.gradient(gx, spacing[1], axis=1)
    fxy = np.gradient(gx, spacing[0], axis=0)
    return fxx, fyy, fxy


def morse_classify(field, parabolic_frac=0.10):
    """Morse classification: peak/trough/saddle/parabolic. NaN → False all."""
    valid = ~np.isnan(field)
    fxx, fyy, fxy = hessian_2d(field)
    trace = fxx + fyy
    det = fxx * fyy - fxy ** 2
    disc = np.maximum(trace ** 2 - 4 * det, 0)
    l1 = (trace + np.sqrt(disc)) / 2.0
    l2 = (trace - np.sqrt(disc)) / 2.0
    adet = np.abs(det)
    thr = np.nanpercentile(adet[valid], parabolic_frac * 100)
    parabolic = (~valid) | (adet < thr)
    peak = (~parabolic) & (l1 < 0) & (l2 < 0)
    trough = (~parabolic) & (l1 > 0) & (l2 > 0)
    saddle = (~parabolic) & (l1 * l2 < 0)
    return {'peak': peak, 'trough': trough, 'saddle': saddle,
            'parabolic': parabolic, 'valid': valid}


def saddle_fraction(field, parabolic_frac=0.10):
    """Saddle fraction = saddle / valid non-parabolic points (颱風域慣例)."""
    m = morse_classify(field, parabolic_frac)
    denom = m['valid'] & ~m['parabolic']
    if denom.sum() == 0:
        return np.nan, m
    return m['saddle'][m['valid'] & ~m['parabolic']].sum() / denom.sum(), m


def main():
    print("=== 載入 OISST ===")
    sst1 = np.load(SST1)
    sst2 = np.load(SST2)
    sst = np.concatenate([sst1, sst2], axis=0)  # (528, 240, 680)
    lat = np.load(LAT)
    lon = np.load(LON)
    print(f"SST: {sst.shape}, lat {lat[0]:.1f}..{lat[-1]:.1f}, lon {lon[0]:.1f}..{lon[-1]:.1f}")
    n_months = len(sst)

    # 月份標記 1982-01 .. 2025-12
    years = np.repeat(np.arange(1982, 2026), 12)
    months = np.tile(np.arange(1, 13), 44)

    print("=== 計算月氣候 (climatology) ===")
    clim = np.full((12, sst.shape[1], sst.shape[2]), np.nan)
    for m in range(12):
        idx = np.where(months == m + 1)[0]
        clim[m] = np.nanmean(sst[idx], axis=0)
    anom = sst - clim[months - 1]  # (528, 240, 680)

    print("=== ONI 相位標籤 ===")
    oni = load_oni()
    phase = np.array(['Neutral'] * n_months)
    for i in range(n_months):
        v = oni.get((years[i], months[i]), np.nan)
        if np.isnan(v):
            continue
        if v >= 0.5:
            phase[i] = 'El Nino'
        elif v <= -0.5:
            phase[i] = 'La Nina'
    print("相位分佈:", {p: int((phase == p).sum()) for p in ['El Nino', 'La Nina', 'Neutral']})

    print("=== Morse 分類 + 鞍點比例（逐月）===")
    results = []
    saddle_map_sum = np.zeros((sst.shape[1], sst.shape[2]))
    valid_map_sum = np.zeros((sst.shape[1], sst.shape[2]))
    for i in range(n_months):
        frac, m = saddle_fraction(anom[i])
        results.append({
            'year': int(years[i]), 'month': int(months[i]),
            'phase': phase[i], 'oni': oni.get((int(years[i]), int(months[i])), None),
            'saddle_frac': None if np.isnan(frac) else float(frac),
            'saddle_n': int(m['saddle'].sum()),
            'valid_n': int(m['valid'].sum()),
        })
        saddle_map_sum += m['saddle']
        valid_map_sum += m['valid']
        if (i + 1) % 60 == 0:
            print(f"  {i+1}/{n_months} 月完成")

    # ========== 分析 1：鞍點比例 × ENSO 相位 ==========
    print("\n=== 鞍點比例 × ENSO 相位 ===")
    by_phase = {}
    for p in ['El Nino', 'La Nina', 'Neutral']:
        vals = [r['saddle_frac'] for r in results if r['phase'] == p and r['saddle_frac'] is not None]
        by_phase[p] = {
            'n': len(vals),
            'mean': float(np.mean(vals)),
            'std': float(np.std(vals)),
            'median': float(np.median(vals)),
            'min': float(np.min(vals)),
            'max': float(np.max(vals)),
        }
        print(f"  {p}: n={by_phase[p]['n']}, mean={by_phase[p]['mean']:.4f}, "
              f"median={by_phase[p]['median']:.4f}, std={by_phase[p]['std']:.4f}")

    # Kruskal-Wallis + Mann-Whitney
    from scipy.stats import kruskal, mannwhitneyu
    el = [r['saddle_frac'] for r in results if r['phase'] == 'El Nino' and r['saddle_frac'] is not None]
    la = [r['saddle_frac'] for r in results if r['phase'] == 'La Nina' and r['saddle_frac'] is not None]
    nt = [r['saddle_frac'] for r in results if r['phase'] == 'Neutral' and r['saddle_frac'] is not None]
    kw = kruskal(el, la, nt)
    mw_el_la = mannwhitneyu(el, la, alternative='two-sided')
    print(f"  Kruskal-Wallis: H={kw.statistic:.2f}, p={kw.pvalue:.6f}")
    print(f"  MWU El Niño vs La Niña: U={mw_el_la.statistic:.0f}, p={mw_el_la.pvalue:.6f}")

    # ONI 連續相關（鞍點比例 ~ ONI）
    xs = [(r['oni'], r['saddle_frac']) for r in results if r['oni'] is not None and r['saddle_frac'] is not None]
    if xs:
        oni_arr = np.array([x[0] for x in xs])
        frac_arr = np.array([x[1] for x in xs])
        from scipy.stats import spearmanr
        rho, p_sp = spearmanr(oni_arr, frac_arr)
        print(f"  Spearman(ONI, saddle_frac): ρ={rho:.4f}, p={p_sp:.6f}")

    # ========== 分析 2：null control（隨機場值 shuffle）==========
    print("\n=== Null control（值 shuffle 20 surrogates，抽 24 月）===")
    rng = np.random.default_rng(42)
    sample_idx = rng.choice(n_months, size=24, replace=False)
    null_fracs = []
    for i in sample_idx:
        field = anom[i].copy()
        valid = ~np.isnan(field)
        vals = field[valid]
        for _ in range(20):
            shuffled = field.copy()
            shuffled[valid] = rng.permutation(vals)
            f, _ = saddle_fraction(shuffled)
            null_fracs.append(f)
    null_fracs = np.array([x for x in null_fracs if x is not None])
    real_sample = [r['saddle_frac'] for r in results if r['saddle_frac'] is not None][:24]
    print(f"  Null: mean={np.mean(null_fracs):.4f}, std={np.std(null_fracs):.4f}")
    print(f"  Real (首批24月): mean={np.mean(real_sample):.4f}")
    z = (np.mean(real_sample) - np.mean(null_fracs)) / np.std(null_fracs)
    print(f"  z = {z:.2f}")

    # ========== 分析 3：鞍點密度 climatology map ==========
    print("\n=== 鞍點密度空間分佈（528 月平均）===")
    saddle_clim = saddle_map_sum / np.maximum(valid_map_sum, 1)
    lat_mask = (lat >= -25) & (lat <= 25)
    saddle_trop = saddle_clim[lat_mask]
    lat_trop = lat[lat_mask]
    # 緯向平均
    zonal = np.nanmean(saddle_trop, axis=1)
    peak_lat = lat_trop[np.nanargmax(zonal)]
    print(f"  鞍點密度緯向峰值喺 {peak_lat:.1f}°")
    # 東/西太平洋對比（東端 270-290E vs 西端暖池 130-160E）
    lon_east = (lon >= 270) & (lon <= 290)
    lon_west = (lon >= 130) & (lon <= 160)
    east_d = np.nanmean(saddle_clim[:, lon_east])
    west_d = np.nanmean(saddle_clim[:, lon_west])
    print(f"  東端 270-290E: {east_d:.4f} vs 西端暖池 130-160E: {west_d:.4f}")

    out = {
        'domain': 'ENSO OISST 0.25',
        'n_months': n_months,
        'saddle_frac_by_phase': by_phase,
        'kruskal': {'H': float(kw.statistic), 'p': float(kw.pvalue)},
        'mwu_el_la': {'U': float(mw_el_la.statistic), 'p': float(mw_el_la.pvalue)},
        'spearman_oni': {'rho': float(rho) if xs else None, 'p': float(p_sp) if xs else None},
        'null_control': {'n': len(null_fracs), 'mean': float(np.mean(null_fracs)),
                         'std': float(np.std(null_fracs)), 'z': float(z)},
        'zonal_peak_lat': float(peak_lat),
        'east_west_density': {'east_270_290': float(east_d), 'west_130_160': float(west_d)},
        'monthly': results,
    }
    outfile = OUT / 'phase30_enso_saddle_ring.json'
    with open(outfile, 'w') as f:
        json.dump(out, f, indent=1, default=str)
    print(f"\n✅ 輸出: {outfile}")
    return out


if __name__ == '__main__':
    main()
