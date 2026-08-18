#!/usr/bin/env python3
"""
ENSO Saddle Ring 域測量 — Phase 31：空間共位 + 相位重排
================================================================
承接 phase30 發現：ENSO 場級鞍點比例 ≈0.564 對相位無敏感度（trivial），
因為 ENSO 冇單一渦旋中心 → 場級平均 = 算子誤用（Saddle Ring 喺颱風域
係「以中心為錨點嘅環帶測量」）。

正確 ENSO adaption（空間結構測量，唔係場級平均）：
  1. 鞍點密度 climatology map vs ρ 場（noise density）空間共位
     → 鞍點係咪聚集喺高 ρ 區域（= 鋒面帶，颱風眼牆剪切帶嘅同構）
  2. 鞍點密度 map × ENSO 相位（El Niño / La Niña / Neutral）
     → 相位會唔會移動鞍點帶？同凝結結構線「位置 ENSO 唔敏感 p=0.28」
     交叉驗證
  3. 東端 270-290°E（週期記憶最強區）vs Nino3.4 中心鞍點密度
  4. 鞍點帶 vs 凝結核心（|K|>p99）Jaccard 共位

null 基準：
  - 隨機 shuffle ρ 場做空間相關 null
  - 相位標籤打亂（surrogate）測重排顯著性
"""
import numpy as np
import json
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu

OUT = Path('/app/working/workspaces/tygtDc/projects/atmosphere/vorticity_dfold/results')
SST1 = '/app/working/workspaces/tygtDc/data/sst/oisst_mon_1982_2020.npy'
SST2 = '/app/working/workspaces/tygtDc/data/sst/oisst_mon_2021_2025.npy'
LAT = '/app/working/workspaces/tygtDc/data/sst/oisst_mon_lat.npy'
LON = '/app/working/workspaces/tygtDc/data/sst/oisst_mon_lon.npy'
ONI = '/app/working/workspaces/tygtDc/projects/enso/release/oni.csv'
RHO_AMP = '/app/working/workspaces/tygtDc/projects/enso/tests/data/rho_seasonal_memory_maps.npz'


def load_oni():
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
            for m, v in enumerate(parts[1:13], 1):
                oni[(yr, m)] = float(v)
    return oni


def hessian_2d(f, spacing=(1.0, 1.0)):
    gy = np.gradient(f, spacing[0], axis=0)
    gx = np.gradient(f, spacing[1], axis=1)
    fyy = np.gradient(gy, spacing[0], axis=0)
    fxx = np.gradient(gx, spacing[1], axis=1)
    fxy = np.gradient(gx, spacing[0], axis=0)
    return fxx, fyy, fxy


def morse_saddle(field, parabolic_frac=0.10):
    """Return saddle mask (NaN-aware)."""
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
    saddle = (~parabolic) & (l1 * l2 < 0)
    return saddle, valid


def kurtosis_curvature(field):
    """K 場（曲率）→ 凝結核心定義用 |K|。用 Gaussian curvature ≈ det(H)。"""
    fxx, fyy, fxy = hessian_2d(field)
    return fxx * fyy - fxy ** 2


def main():
    print("=== 載入 OISST 528 月 ===")
    sst = np.concatenate([np.load(SST1), np.load(SST2)], axis=0)
    lat = np.load(LAT)
    lon = np.load(LON)
    n_months = len(sst)
    years = np.repeat(np.arange(1982, 2026), 12)
    months = np.tile(np.arange(1, 13), 44)

    print("=== 月氣候 + anomaly ===")
    clim = np.full((12, sst.shape[1], sst.shape[2]), np.nan)
    for m in range(12):
        idx = np.where(months == m + 1)[0]
        clim[m] = np.nanmean(sst[idx], axis=0)
    anom = sst - clim[months - 1]

    print("=== 相位標籤 ===")
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

    print("=== 逐月 Morse 鞍點 mask（528 月）===")
    saddle_map_sum = np.zeros((sst.shape[1], sst.shape[2]))
    valid_map_sum = np.zeros((sst.shape[1], sst.shape[2]))
    k2_map = np.zeros((sst.shape[1], sst.shape[2]))
    saddle_by_phase = {p: np.zeros((sst.shape[1], sst.shape[2])) for p in ['El Nino', 'La Nina', 'Neutral']}
    count_by_phase = {p: 0 for p in ['El Nino', 'La Nina', 'Neutral']}
    for i in range(n_months):
        saddle, valid = morse_saddle(anom[i])
        saddle_map_sum += saddle
        valid_map_sum += valid
        # K 場累加（平方再開方，避免正負抵消）
        K = kurtosis_curvature(anom[i])
        k2_map += np.nan_to_num(K ** 2)
        p = phase[i]
        saddle_by_phase[p] += saddle
        count_by_phase[p] += 1
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{n_months}")

    # ========== 分析 1：鞍點密度 vs ρ 場空間共位 ==========
    print("\n=== 分析 1：鞍點密度 vs ρ 場（noise density）空間共位 ===")
    saddle_clim = saddle_map_sum / np.maximum(valid_map_sum, 1)
    rho = np.load(RHO_AMP)['rho_amp']
    valid = ~np.isnan(saddle_clim) & ~np.isnan(rho)
    rho_flat = rho[valid].ravel()
    saddle_flat = saddle_clim[valid].ravel()

    rho_sp, p_sp = spearmanr(rho_flat, saddle_flat)
    print(f"  Spearman(ρ, saddle_density): ρ={rho_sp:.4f}, p={p_sp:.2e}")

    # 分層：高 ρ (top 20%) vs 低 ρ (bottom 20%)
    rho_thr_hi = np.percentile(rho_flat, 80)
    rho_thr_lo = np.percentile(rho_flat, 20)
    hi = saddle_flat[rho_flat >= rho_thr_hi]
    lo = saddle_flat[rho_flat <= rho_thr_lo]
    print(f"  高 ρ (top20%): saddle_density mean={np.mean(hi):.4f} (n={len(hi)})")
    print(f"  低 ρ (bot20%): saddle_density mean={np.mean(lo):.4f} (n={len(lo)})")
    mw = mannwhitneyu(hi, lo, alternative='greater')
    print(f"  MWU 高>低: U={mw.statistic:.0f}, p={mw.pvalue:.2e}")

    # null：shuffle ρ 場
    rng = np.random.default_rng(42)
    nulls = []
    for _ in range(20):
        rho_shuf = rho_flat.copy()
        rng.shuffle(rho_shuf)
        rs, _ = spearmanr(rho_shuf, saddle_flat)
        nulls.append(rs)
    nulls = np.array(nulls)
    z = (rho_sp - np.mean(nulls)) / np.std(nulls)
    print(f"  Null shuffle ρ: mean ρ={np.mean(nulls):.4f}, std={np.std(nulls):.4f}, z={z:.2f}")

    # ========== 分析 2：鞍點密度 × ENSO 相位（空間重排）==========
    print("\n=== 分析 2：鞍點密度 × ENSO 相位 ===")
    phase_maps = {}
    for p in ['El Nino', 'La Nina', 'Neutral']:
        phase_maps[p] = saddle_by_phase[p] / count_by_phase[p]
        print(f"  {p}: n={count_by_phase[p]}, 場均值={np.nanmean(phase_maps[p]):.4f}")

    # 空間相關：El Niño vs La Niña 鞍點密度 map 係咪同一結構？
    r_el_la, p_el_la = spearmanr(phase_maps['El Nino'][valid].ravel(),
                                 phase_maps['La Nina'][valid].ravel())
    print(f"  鞍點密度 map 相關 (El Niño vs La Niña): ρ={r_el_la:.4f}, p={p_el_la:.2e}")
    r_el_nt, _ = spearmanr(phase_maps['El Nino'][valid].ravel(),
                           phase_maps['Neutral'][valid].ravel())
    r_la_nt, _ = spearmanr(phase_maps['La Nina'][valid].ravel(),
                           phase_maps['Neutral'][valid].ravel())
    print(f"  El Niño vs Neutral: ρ={r_el_nt:.4f} | La Niña vs Neutral: ρ={r_la_nt:.4f}")

    # ========== 分析 3：區域對比（東端 vs 中心 vs 西端）==========
    print("\n=== 分析 3：區域鞍點密度（climatology）===")
    lon_east = (lon >= 270) & (lon <= 290)   # 週期記憶最強區
    lon_center = (lon >= 210) & (lon <= 270)  # Nino3.4 附近（190-240 係 Nino3.4 核心，但中心記憶最弱係 ~210-270）
    lon_west = (lon >= 130) & (lon <= 160)   # 暖池邊緣
    lat_trop = (lat >= -10) & (lat <= 10)
    for name, lm in [('東端 270-290E', lon_east), ('中心 210-270E', lon_center), ('西端 130-160E', lon_west)]:
        sub = saddle_clim[np.ix_(lat_trop, lm)]
        sub_rho = rho[np.ix_(lat_trop, lm)]
        print(f"  {name}: 鞍點密度={np.nanmean(sub):.4f}, ρ={np.nanmean(sub_rho):.4f}")

    # ========== 分析 4：鞍點 vs 凝結核心（|K|>p99）共位 ==========
    print("\n=== 分析 4：鞍點密度 vs 凝結核心 K 場 ===")
    # K2 map（每月平方累加）→ sqrt = rms 曲率
    k_rms = np.sqrt(k2_map / n_months)
    valid_k = ~np.isnan(k_rms) & ~np.isnan(saddle_clim)
    k_flat = k_rms[valid_k].ravel()
    sd_flat = saddle_clim[valid_k].ravel()
    k_sp, k_p = spearmanr(k_flat, sd_flat)
    print(f"  Spearman(K_rms, saddle_density): ρ={k_sp:.4f}, p={k_p:.2e}")

    # 高 K (top5%) vs 低 K (bottom50%) 鞍點密度
    k_thr_hi = np.percentile(k_flat, 95)
    k_thr_lo = np.percentile(k_flat, 50)
    khi = sd_flat[k_flat >= k_thr_hi]
    klo = sd_flat[k_flat <= k_thr_lo]
    print(f"  高 K (top5%): saddle_density={np.mean(khi):.4f} (n={len(khi)})")
    print(f"  低 K (bot50%): saddle_density={np.mean(klo):.4f} (n={len(klo)})")
    mw2 = mannwhitneyu(khi, klo, alternative='greater')
    print(f"  MWU 高K>低K: U={mw2.statistic:.0f}, p={mw2.pvalue:.2e}")

    # null：shuffle K 場
    null_k = []
    for _ in range(20):
        k_shuf = k_flat.copy()
        rng.shuffle(k_shuf)
        rs, _ = spearmanr(k_shuf, sd_flat)
        null_k.append(rs)
    null_k = np.array(null_k)
    zk = (k_sp - np.mean(null_k)) / np.std(null_k)
    print(f"  Null shuffle K: mean={np.mean(null_k):.4f}, std={np.std(null_k):.4f}, z={zk:.2f}")

    out = {
        'phase': 'ENSO Saddle Ring Phase 31 — 空間共位',
        'rho_spatial': {'spearman_rho': float(rho_sp), 'p': float(p_sp),
                        'hi_rho_mean': float(np.mean(hi)), 'lo_rho_mean': float(np.mean(lo)),
                        'mwu_p': float(mw.pvalue), 'null_z': float(z)},
        'phase_maps': {
            'el_la_corr': float(r_el_la), 'el_nt_corr': float(r_el_nt), 'la_nt_corr': float(r_la_nt),
            'el_mean': float(np.nanmean(phase_maps['El Nino'])),
            'la_mean': float(np.nanmean(phase_maps['La Nina'])),
            'nt_mean': float(np.nanmean(phase_maps['Neutral'])),
        },
        'regions': {
            'east_270_290': {'saddle': float(np.nanmean(saddle_clim[np.ix_(lat_trop, lon_east)])),
                             'rho': float(np.nanmean(rho[np.ix_(lat_trop, lon_east)]))},
            'center_210_270': {'saddle': float(np.nanmean(saddle_clim[np.ix_(lat_trop, lon_center)])),
                               'rho': float(np.nanmean(rho[np.ix_(lat_trop, lon_center)]))},
            'west_130_160': {'saddle': float(np.nanmean(saddle_clim[np.ix_(lat_trop, lon_west)])),
                             'rho': float(np.nanmean(rho[np.ix_(lat_trop, lon_west)]))},
        },
        'core_coloc': {'spearman_k': float(k_sp), 'p': float(k_p),
                       'hi_k_mean': float(np.mean(khi)), 'lo_k_mean': float(np.mean(klo)),
                       'mwu_p': float(mw2.pvalue), 'null_z': float(zk)},
    }
    outfile = OUT / 'phase31_enso_saddle_spatial.json'
    with open(outfile, 'w') as f:
        json.dump(out, f, indent=1, default=str)
    print(f"\n✅ 輸出: {outfile}")
    return out


if __name__ == '__main__':
    main()
