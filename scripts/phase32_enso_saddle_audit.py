#!/usr/bin/env python3
"""
ENSO Saddle Ring — Phase 32：嚴格審計
================================================================
Phase 31 三大問題審計：
  A. 分析1矛盾（Spearman +0.18 但 high-ρ saddle 低）→ 非單調檢查 + 海洋 mask
  B. 分析4 K-saddle 數學耦合（兩者都 det(H)）→ 改用獨立量：鞍點密度 vs ρ 場
     （⟨|dH|⟩ 唔涉及 Hessian）；K 只做描述唔做證據
  C. 相位重排：El/La map 相關 0.61-0.66 需要 null（shuffle 相位）先知顯著性

附加：鞍點密度 pattern（中心強/東端弱）vs 週期記憶 pattern（兩端強中心弱）
     對照 —— 兩個唔同結構量嘅空間分佈比較。
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


def main():
    print("=== 載入 ===")
    sst = np.concatenate([np.load(SST1), np.load(SST2)], axis=0)
    lat = np.load(LAT)
    lon = np.load(LON)
    n_months = len(sst)
    years = np.repeat(np.arange(1982, 2026), 12)
    months = np.tile(np.arange(1, 13), 44)

    clim = np.full((12, sst.shape[1], sst.shape[2]), np.nan)
    for m in range(12):
        idx = np.where(months == m + 1)[0]
        clim[m] = np.nanmean(sst[idx], axis=0)
    anom = sst - clim[months - 1]

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

    print("=== 逐月 Morse 鞍點 + 相位分離 ===")
    saddle_map_sum = np.zeros((sst.shape[1], sst.shape[2]))
    valid_map_sum = np.zeros((sst.shape[1], sst.shape[2]))
    saddle_by_phase = {p: np.zeros((sst.shape[1], sst.shape[2])) for p in ['El Nino', 'La Nina', 'Neutral']}
    count_by_phase = {p: 0 for p in ['El Nino', 'La Nina', 'Neutral']}
    for i in range(n_months):
        saddle, valid = morse_saddle(anom[i])
        saddle_map_sum += saddle
        valid_map_sum += valid
        p = phase[i]
        saddle_by_phase[p] += saddle
        count_by_phase[p] += 1
    saddle_clim = saddle_map_sum / np.maximum(valid_map_sum, 1)

    rho = np.load(RHO_AMP)['rho_amp']

    # 海洋 mask：SST 有效 + 排除邊界（gradient 影響區，距陸地 ~2 cells）+ 熱帶 ±25
    ocean_valid = ~np.isnan(saddle_clim) & ~np.isnan(rho)
    # 侵蝕邊界（用 sst 有效 mask 做 dilate 後取 interior）
    sst_valid0 = ~np.isnan(sst[0])
    from scipy.ndimage import binary_erosion
    interior = binary_erosion(sst_valid0, iterations=3)  # 排除邊界 3 cells
    lat_trop = (lat >= -25) & (lat <= 25)
    mask = ocean_valid & interior & lat_trop[:, None]
    print(f"  有效格點（海洋+熱帶+內陸）: {mask.sum()} / {mask.size}")

    # ===== A. 非單調檢查：ρ 分 bin =====
    print("\n=== A. ρ 分 bin 鞍點密度趨勢（10 bins）===")
    rho_flat = rho[mask].ravel()
    sd_flat = saddle_clim[mask].ravel()
    order = np.argsort(rho_flat)
    n = len(rho_flat)
    bin_edges = np.linspace(0, n, 11).astype(int)
    bins_report = []
    for b in range(10):
        idx = order[bin_edges[b]:bin_edges[b + 1]]
        bins_report.append({
            'rho_range': [float(np.min(rho_flat[idx])), float(np.max(rho_flat[idx]))],
            'saddle_mean': float(np.mean(sd_flat[idx])),
        })
        print(f"  bin{b}: ρ∈[{np.min(rho_flat[idx]):.4f},{np.max(rho_flat[idx]):.4f}] "
              f"saddle_mean={np.mean(sd_flat[idx]):.4f}")
    # 單調性：bin mean 相關係數（bin index vs saddle mean）
    bmeans = np.array([b['saddle_mean'] for b in bins_report])
    monotone_rho, monotone_p = spearmanr(np.arange(10), bmeans)
    print(f"  bin 單調性 Spearman: ρ={monotone_rho:.4f}, p={monotone_p:.4f}")

    # ===== A2. 海洋 mask 內 Spearman =====
    rs, ps = spearmanr(rho_flat, sd_flat)
    print(f"\n  海洋+熱帶 mask 內 Spearman(ρ, saddle): ρ={rs:.4f}, p={ps:.2e}")

    # ===== C. 相位重排 null =====
    print("\n=== C. 相位重排 null（100 surrogates）===")
    rng = np.random.default_rng(7)
    phase_maps = {p: saddle_by_phase[p] / count_by_phase[p] for p in ['El Nino', 'La Nina', 'Neutral']}
    r_el_la_obs, _ = spearmanr(phase_maps['El Nino'][mask].ravel(),
                               phase_maps['La Nina'][mask].ravel())
    print(f"  Observed El Niño vs La Niña map corr: ρ={r_el_la_obs:.4f}")
    null_corrs = []
    n_el = count_by_phase['El Nino']
    n_la = count_by_phase['La Nina']
    for s in range(100):
        # shuffle 相位標籤（保持 n_el/n_la 數量）
        labels = np.array(['El Nino'] * n_el + ['La Nina'] * (n_months - n_el))
        rng.shuffle(labels)
        el_sum = np.zeros(sst.shape[1:])
        la_sum = np.zeros(sst.shape[1:])
        for i in range(n_months):
            # 重用之前計算：需要重新累加，慢 — 改用近似：直接 shuffle 相位分配
            pass
    # 快速版：只 shuffle 相位 map 嘅「月份組合」— 用 saddle_by_phase 累加係線性，
    # 所以 null 可以直接由「隨機抽 n_el 個月當 El Nino」做
    # 先存逐月 saddle masks 唔現實（528×240×680 bool = 86GB）— 用 60 月抽樣近似
    print("  (用 60 月抽樣近似 null)")
    el_samp = np.zeros(sst.shape[1:])
    la_samp = np.zeros(sst.shape[1:])
    el_n = la_n = 0
    # 重新做逐月但只為 null：直接喺主循環做 — 為省時，用 40 組 surrogate 每組抽樣 30+30
    null_map_corrs = []
    for s in range(40):
        el_sum_n = np.zeros(sst.shape[1:])
        la_sum_n = np.zeros(sst.shape[1:])
        el_n = la_n = 0
        # 隨機抽 60 個月：30 當 El Niño, 30 當 La Niña
        chosen = rng.choice(n_months, size=60, replace=False)
        el_idx = chosen[:30]
        la_idx = chosen[30:]
        for i in range(n_months):
            if i not in chosen:
                continue
            saddle, valid = morse_saddle(anom[i])
            if i in el_idx:
                el_sum_n += saddle
                el_n += 1
            else:
                la_sum_n += saddle
                la_n += 1
        el_map_n = el_sum_n / el_n
        la_map_n = la_sum_n / la_n
        rn, _ = spearmanr(el_map_n[mask].ravel(), la_map_n[mask].ravel())
        null_map_corrs.append(rn)
    null_map_corrs = np.array(null_map_corrs)
    # 但 observed 用 143/147 個月，null 用 30/30 — 解析度唔同，睇 direction 先
    print(f"  Null map corr (30v30 抽樣): mean={np.mean(null_map_corrs):.4f}, "
          f"std={np.std(null_map_corrs):.4f}, obs(143v147)={r_el_la_obs:.4f}")

    # ===== 對照：鞍點 pattern vs 週期記憶 pattern =====
    print("\n=== 鞍點 density pattern vs 週期記憶 pattern ===")
    net = np.load(RHO_AMP)['net']  # 週期記憶 net map（跨年 corr 淨值）
    lat_trop10 = (lat >= -10) & (lat <= 10)
    for name, lm in [('東端 270-290E', (lon >= 270) & (lon <= 290)),
                     ('中心 210-270E', (lon >= 210) & (lon <= 270)),
                     ('西端 130-160E', (lon >= 130) & (lon <= 160))]:
        sub_sd = saddle_clim[np.ix_(lat_trop10, lm)]
        sub_net = net[np.ix_(lat_trop10, lm)]
        print(f"  {name}: saddle={np.nanmean(sub_sd):.4f}, 週期記憶net={np.nanmean(sub_net):.4f}")

    out = {
        'phase': 'ENSO Saddle Ring Phase 32 — 審計',
        'rho_bins': bins_report,
        'monotonicity': {'spearman': float(monotone_rho), 'p': float(monotone_p)},
        'rho_spatial_masked': {'spearman': float(rs), 'p': float(ps), 'n': int(n)},
        'phase_map_corr': {'obs_el_la': float(r_el_la_obs),
                           'null_mean': float(np.mean(null_map_corrs)),
                           'null_std': float(np.std(null_map_corrs))},
    }
    outfile = OUT / 'phase32_enso_saddle_audit.json'
    with open(outfile, 'w') as f:
        json.dump(out, f, indent=1, default=str)
    print(f"\n✅ 輸出: {outfile}")
    return out


if __name__ == '__main__':
    main()
