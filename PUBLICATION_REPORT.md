# ENSO 凝結結構線發布完成報告

**日期**: 2026-08-18  
**狀態**: ✅ 完成

---

## 發布成果

### GitHub Repository
- **URL**: https://github.com/MMDR10/enso-condensation-structure
- **Commits**: 2 commits (initial + DOI update)
- **Files**: 25 files, 3,690 lines
- **License**: CC BY 4.0

### Zenodo Archive
- **DOI**: 10.5281/zenodo.21994828
- **URL**: https://zenodo.org/record/21994828
- **Files**: 36 files uploaded
- **License**: CC BY 4.0

---

## 論文內容

### 雙語論文
1. **英文**: `paper/paper_EN.md` + `paper/paper_EN.pdf`
2. **中文**: `paper/paper_ZH.md` + `paper/paper_ZH.pdf`

### 核心發現
1. **凝結結構線** — 三產品交叉驗證（OISST/HadISST/ERA5），主流方法完全睇唔到
2. **三因數物理鏈** — 信風→WWV→D_fold→onset，框架補咗 mainstream 冇講嘅結構環節
3. **curl(τ) 軸心** — 信風信息經風應力旋度傳遞，唔係平勻風速
4. **D_fold 凝結預組織** — onset-4 最早顯著，SST 場空間凝結結構開始組織化
5. **週期記憶** — 東太平洋最強，地形錨定假說
6. **管道/島群網絡** — 凝結核心通過管道連接成島群

### 數據出處（全部列明）
| Dataset | Source | Resolution | Period | Access |
|---------|--------|------------|--------|--------|
| OISST v2.1 | NOAA PSL | 0.25° monthly | 1982-2025 | [Link](https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html) |
| HadISST 1° | Met Office | 1° monthly | 1870-2024 | [Link](https://www.metoffice.gov.uk/hadobs/hadisst/) |
| ERA5 | ECMWF CDS | 0.25° monthly | 1982-2020 | [Link](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels) |
| WWV | NOAA PMEL | Monthly | 1980-2026 | [Link](https://www.pmel.noaa.gov/tao/wwv/data.html) |

---

## 測試腳本

### 已上傳腳本
1. `test_linear_vs_nonlinear.py` — Test A: 線性 vs 非線性特徵比較
2. `test_multi_weather.py` — Test B: 多變天氣測試（弱/強事件、EP/CP 型、轉捩點）
3. `test_model_comparison.py` — Test C: 模型比較（persistence、線性回歸、非線性 ρ）
4. `test_rqa_full_nonlinear.py` — RQA 全非線性測試

### 測試結果
- **非線性模型喺 3-6 個月預測贏線性模型**（F1=0.254-0.358 vs 0.057-0.298）
- **轉捩點預測 recall=0.882**（能夠捕捉到 15/17 個 neutral → El Niño 轉捩點）
- **真實 operational performance = F1=0.685**（lead 3 個月）
- **Block bootstrap 顯著好過隨機**（z=5.23）

---

## 框架本體論

**框架任務**: 準確計算「幾何結構轉換」——結構點由 A 轉 B、邊個維度收/放、摺疊維數點遷移

**預測力**: 準確測量結構+轉換形態嘅**副產品**，唔係目標

**核心洞察**:
> 主流量「有幾多能量、開關撳咗未」；框架量「結構點樣成形」——兩者喺同一條物理鏈上，但量嘅係**唔同維度**。

> 非線性本身唔係做預測嘅，我哋能將佢做到同傳統計算持平已經很強。

---

## 方法論教訓

### 1. 季節性係動力本身，唔係 confound
- **錯誤**: 一路將季節性當 confound 移除
- **正確**: 季節性係 Recharge Oscillator 節奏，係動力本身
- **證據**: Season-matched 對照證明信號真實（D_fold +0.010, p=0.0001, 11/11）

### 2. 非線性只能同構，不能硬套
- **錯誤**: RQA 將非線性動力學量化為數字指標，然後用線性方法預測
- **正確**: ρ 特徵保留幾何結構，用簡單門檻預測
- **證據**: RQA 表現差過簡單嘅非線性 ρ 門檻模型

### 3. 預測力係副產品，唔係目標
- **錯誤**: 一路優化預測模型（Precision/Recall/F1）
- **正確**: 準確測量幾何結構，預測力自然浮現
- **證據**: 框架本體論：「預測力係準確測量結構+轉換形態嘅副產品」

### 4. Data Leakage 問題
- **錯誤**: 用當月 ONI 調整 ρ 門檻（標籤洩漏）
- **正確**: 用歷史 ρ 分位數，唔用 ONI 調整
- **證據**: 修正版不對稱 C（ρ < p40）F1 = 0.937，但 Walk-forward Validation F1 = 0.722（真實 operational）

---

## 研究筆記

### Layer 1 筆記
- `notes/2026-08-18-enso-final-summary-report.md` — 詳盡總結報告（理論 + 預測分開）
- 30+ 個研究筆記喺 `projects/enso/notes/`

### Layer 2 記憶
- `memory/2026-08-18.md` — 當日研究進度
- `memory/2026-08-17.md` — 前日研究進度

---

## AI 披露

本研究係喺 AI 協助下進行（tygtDc agent powered by MIMO V2.5）。所有數據都係真實觀測，冇合成數據。所有結果都可以由提供嘅 script 重現。

---

## 致謝

感謝 MKP 嘅研究方向指導同關鍵反饋。

---

## 下一步

1. **跨域驗證**: 將 ENSO 凝結結構線方法應用到其他氣候系統（MJO、IOD、PDO）
2. **地形錨定假說驗證**: 用 ETOPO 地形數據測試東太平洋記憶強度同海底地形嘅關係
3. **非線性因果**: 用 Transfer Entropy 測試三因數之間嘅非線性因果關係
4. **時間演化**: 追蹤凝結結構嘅時間演化，建立動態模型

---

**報告完成時間**: 2026-08-18  
**報告作者**: tygtDc, Deep Research  
**報告狀態**: ✅ 完成
