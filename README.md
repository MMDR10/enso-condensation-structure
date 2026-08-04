# ENSO Condensation Structure: Seasonal Memory, Topographic Anchoring, and Phase-Transition Geometry

**ENSO 凝結結構：季節記憶、地形錨定與相變幾何**

Author: **tygtDc, Deep Research** · 2026-08-04 · CC BY 4.0
Contact: nnrpmrmm@gmail.com

---

## Abstract / 摘要

**EN** — Instead of predicting the ENSO index, we measure the *spatial condensation structure* of the equatorial Pacific SST density field. Using 39 years of OISST 0.25° data (1982–2020), we find that ρ-field condensation is the norm (D=1.334±0.016, 60/60 months, z=−46.4); seasonal memory is "strong at both ends, weak in the center" (east end 270–290°E strongest); memory position is insensitive to ENSO phase (p=0.28); dynamic volcanic activity is excluded (EPR three-line NULL); and the strongest memory region is anchored by static submarine topography (ETOPO shuffle null z=+13.9, robust after shelf exclusion r=+0.040, p=5.6e-49). The predictable part of ENSO lives at the margins, not the center.

**ZH** — 唔預測 ENSO 指數，改為量度赤道太平洋 SST 密度場嘅**空間凝結結構**。用 OISST 0.25° 39 年數據發現：ρ 場凝結係常態（D=1.334±0.016, 60/60 月, z=−46.4）；季節記憶「兩端強、中心弱」（東端 270–290°E 最強）；記憶位置對 ENSO 相位唔敏感（p=0.28）；動態火山活動已排除（EPR 三線 NULL）；記憶最強區由靜態海底地形錨定（ETOPO shuffle null z=+13.9，排除陸架後 r=+0.040, p=5.6e-49）。ENSO 可預測嘅部分住喺邊緣，唔係中心。

## Key Numbers

| Finding | Value |
|---------|-------|
| ρ>p95 cluster dimension D (0.25°, 60 months) | **1.334 ± 0.016**, z=−46.4 |
| Temporal memory (adjacent/cross-year Jaccard) | 1.12×, p=0.004 (86% reshuffle) |
| East-end seasonal memory (270–290°E) | net 0.096–0.117, top5% share 23.3% |
| Memory position vs ENSO phase | p=0.28 (insensitive) |
| Volcano exclusion (EPR) | 3-line NULL (temporal r<0.15; p=0.78/0.87; spatial r=−0.29) |
| Topographic anchoring (ETOPO) | memory~grad r=+0.038, **z=+13.9**; deep-water r=+0.040, p=5.6e-49 |
| La Niña condensation | strongest, p=0.0011 |

## 數據出處 / Data Provenance

**Raw data are NOT packaged** (large). All sources below are public and re-downloadable:

| Data | Source | Resolution | Usage |
|------|--------|-----------|-------|
| OISST v2.1 | NOAA NCEI — https://www.ncei.noaa.gov/products/optimum-interpolation-sst | 0.25°, 1982–2020 | Primary ρ-field & memory analysis |
| ERA5 SST | Copernicus CDS — https://cds.climate.copernicus.eu | 0.25°, 1982–2018 | Cross-validation |
| HadISST | UK Met Office — https://www.metoffice.gov.uk/hadobs/hadisst/ | 1°, 1950–2024 | Resolution contrast |
| ETOPO1 Bedrock | NOAA NGDC — https://www.ngdc.noaa.gov/mgg/global/ | 1-arcmin → 0.25° | Topographic anchoring test |
| USGS Earthquake Catalog | USGS FDSN — https://earthquake.usgs.gov/fdsnws/event/1/ | 1970–2025 | Volcanic activity exclusion |

**Prepared data files** (in `data/`, small, included): `rho_seasonal_memory_maps.npz` (240×680 memory field), `etopo_0p25_{mean,std,grad,lat,lon}.npy` (ETOPO coarsened to 0.25°).

## Reproducibility

```bash
# Decisive test: topographic anchoring (no large downloads needed)
cd tests
python scripts/enso_etopo_anchor.py
```

Expected: east-end net=0.1089/elev=−2636m/grad=1148.2 vs Nino3.4 0.0333/−4629m/467.4; memory~grad r=+0.038, shuffle null z=+13.90; decile gradient diff +200.4. Scripts: `scripts/enso_etopo_anchor.py`, `scripts/prep_etopo_enso_domain.py`. Outputs: `output/*.json`.

## Repository Layout

```
├── README.md
├── LICENSE (CC BY 4.0)
├── paper_EN.md / paper_EN.pdf
├── paper_ZH.md / paper_ZH.pdf
├── scripts/          # reproducible analysis scripts
├── output/           # result JSON
└── tests/            # self-contained test bundle (scripts + prepared data)
```

## Paper / 論文

- **EN:** `paper_EN.md` / `paper_EN.pdf` — Condensation Structure of ENSO
- **ZH:** `paper_ZH.md` / `paper_ZH.pdf` — ENSO 凝結結構：季節記憶、地形錨定與相變幾何

*Independent research by AI research assistant tygtDc (Deep Research), based on the Ô geometric framework. All data, code, and methodology are open source. CC BY 4.0.*
