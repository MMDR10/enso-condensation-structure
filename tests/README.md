# ENSO 凝結結構 — 測試檔（2026-08-04 更新）

用嚟重跑 ENSO 凝結結構論文（`docs/ENSO_CONDENSATION_PAPER.md`）嘅決定性測試。
**MKP 可以直接喺自己電腦跑，唔使再下載原始數據。**

## 內容

```
tests/
├── scripts/
│   ├── enso_etopo_anchor.py      # 地形錨定決定性測試（ETOPO × 記憶場）
│   └── prep_etopo_enso_domain.py # 抽 ETOPO 子區域 + coarsen（可選，數據已預備）
└── data/
    ├── rho_seasonal_memory_maps.npz  # 季節記憶場（240×680, 0.25°, 39年 OISST）
    ├── etopo_0p25_mean.npy           # ETOPO1 平均海拔（0.25°）
    ├── etopo_0p25_std.npy            # 粗糙度（0.25° std）
    ├── etopo_0p25_grad.npy           # 梯度 magnitude
    ├── etopo_0p25_lat.npy            # 緯度座標
    └── etopo_0p25_lon.npy            # 經度座標
```

## 點跑

```bash
cd tests
python scripts/enso_etopo_anchor.py
```

**冇依賴原始大數據**（ETOPO1 2GB grd、OISST 583MB 已處理成上面嘅 npy/npz）。
需要：numpy, scipy。

## 預期輸出（同論文 §3.6 一致）

- 東端 270–290°E：net memory 0.109、海拔 −2636m、粗糙度 183、梯度 1148
- Nino3.4：net 0.033、海拔 −4629m、梯度 467
- memory~grad r=+0.038、shuffle null z≈+13.9
- 陸架排除（>1000m）：r=+0.040, p=5.6e-49

## 其他 ENSO 測試（完整 scripts 喺 workspace `scripts/`）

| Script | 測咩 | 數據 |
|--------|------|------|
| `enso_phase_condensation.py` | La Niña 凝結最強（p=0.0011） | OISST 39 年 |
| `oisst_0p25_60mo_rho.py` | ρ 場凝結 D=1.334, z=−46.4（60 月） | OISST 0.25° |
| `rho_seasonal_memory_39yr.py` | 週期記憶「兩端強、中心弱」 | OISST 39 年 |
| `rho_seasonal_memory_localization.py` | 東端 270–290°E 記憶最強（0.181 vs 0.142） | OISST 39 年 |
| `report_enso_etopo_anchor.py` | 地形錨定（本包） | 見上 |

## 論文結論一覽（v1.0, 2026-08-04）

1. ρ 場凝結常態：D=1.334±0.016, 60/60 月, z=−46.4
2. 時間記憶微弱：相鄰月/跨年 = 1.12×, p=0.004（86% 重排）
3. 週期記憶「兩端強、中心弱」：東端 net 0.096–0.117 vs 中心最弱
4. 記憶位置對 ENSO 唔敏感：p=0.28 → 由非 ENSO 因素錨定
5. 火山活動排除：EPR 三線 NULL（時間 r<0.15 / 高低月 p=0.78,0.87 / 空間 r=−0.29）
6. 地形錨定支持：ETOPO z=+13.9、陸架排除後保持（p=5.6e-49）
7. La Niña 凝結最強：p=0.0011

**待驗證：** EUC 上湧→鋒面錨定機制（相關支持、機制未直接測）。

---
*tygtDc, Deep Research · 2026-08-04 · CC BY 4.0*
