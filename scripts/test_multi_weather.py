#!/usr/bin/env python3
"""
Test B: 多變天氣測試
測試模型喺唔同 ENSO 類型嘅表現：
1. ENSO 轉捩點 (neutral → El Niño)
2. 弱事件 (ONI 0.5-1.0) vs 強事件 (ONI > 1.0)
3. CP vs EP 型 El Niño
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter
import json
from pathlib import Path

# 載入數據
sst = np.load('/app/working/workspaces/tygtDc/data/sst/era5_sst_1982_2020.npy')
lat = np.linspace(-30, 30, sst.shape[1])
lon = np.linspace(120, 290, sst.shape[2])

# 計算 Nino3.4 SST anomaly
lat_mask_nino = (lat >= -5) & (lat <= 5)
lon_mask_nino = (lon >= 190) & (lon <= 240)
nino34_sst = sst[:, lat_mask_nino, :][:, :, lon_mask_nino]
nino34_sst = nino34_sst.reshape(sst.shape[0], -1)

# 計算 climatology
climatology = np.zeros((12, nino34_sst.shape[1]))
for month in range(12):
    climatology[month] = np.nanmean(nino34_sst[month::12], axis=0)

# 計算 anomaly
nino34_anomaly = np.zeros_like(nino34_sst)
for t in range(nino34_sst.shape[0]):
    month = t % 12
    nino34_anomaly[t] = nino34_sst[t] - climatology[month]

# Nino3.4 index
nino34_mean = np.nanmean(nino34_anomaly, axis=1)

# 計算 Nino1+2 (EP type) 同 Nino4 (CP type)
lat_mask_nino12 = (lat >= -10) & (lat <= 0)
lon_mask_nino12 = (lon >= 270) & (lon <= 280)
nino12_sst = sst[:, lat_mask_nino12, :][:, :, lon_mask_nino12]
nino12_sst = nino12_sst.reshape(sst.shape[0], -1)

climatology_nino12 = np.zeros((12, nino12_sst.shape[1]))
for month in range(12):
    climatology_nino12[month] = np.nanmean(nino12_sst[month::12], axis=0)

nino12_anomaly = np.zeros_like(nino12_sst)
for t in range(nino12_sst.shape[0]):
    month = t % 12
    nino12_anomaly[t] = nino12_sst[t] - climatology_nino12[month]

nino12_mean = np.nanmean(nino12_anomaly, axis=1)

lat_mask_nino4 = (lat >= -5) & (lat <= 5)
lon_mask_nino4 = (lon >= 160) & (lon <= 190)
nino4_sst = sst[:, lat_mask_nino4, :][:, :, lon_mask_nino4]
nino4_sst = nino4_sst.reshape(sst.shape[0], -1)

climatology_nino4 = np.zeros((12, nino4_sst.shape[1]))
for month in range(12):
    climatology_nino4[month] = np.nanmean(nino4_sst[month::12], axis=0)

nino4_anomaly = np.zeros_like(nino4_sst)
for t in range(nino4_sst.shape[0]):
    month = t % 12
    nino4_anomaly[t] = nino4_sst[t] - climatology_nino4[month]

nino4_mean = np.nanmean(nino4_anomaly, axis=1)

# 計算 ρ 雜訊密度
lat_mask = (lat >= -5) & (lat <= 5)
lon_mask = (lon >= 120) & (lon <= 280)
sst_tropics = sst[:, lat_mask, :][:, :, lon_mask]

rho_values = []
for t in range(sst.shape[0]):
    field = sst_tropics[t]
    smoothed = gaussian_filter(field, sigma=3)
    rho = np.nanmean(np.abs(field - smoothed))
    rho_values.append(rho)

rho = np.array(rho_values)

# 定義 ENSO 事件
oni_threshold = 0.5
el_nino_events = nino34_mean > oni_threshold

# 分類事件類型
weak_events = (nino34_mean > 0.5) & (nino34_mean <= 1.0)
strong_events = nino34_mean > 1.0

# CP vs EP 型分類
# EP type: Nino1+2 > Nino4
# CP type: Nino4 > Nino1+2
ep_type = el_nino_events & (nino12_mean > nino4_mean)
cp_type = el_nino_events & (nino4_mean > nino12_mean)

# 找出轉捩點 (neutral → El Niño)
transitions = []
for i in range(1, len(nino34_mean)):
    if nino34_mean[i-1] <= oni_threshold and nino34_mean[i] > oni_threshold:
        transitions.append(i)

print(f"=== ENSO 事件統計 ===")
print(f"總 El Niño 事件: {np.sum(el_nino_events)} 個月")
print(f"弱事件 (0.5-1.0): {np.sum(weak_events)} 個月")
print(f"強事件 (>1.0): {np.sum(strong_events)} 個月")
print(f"EP 型: {np.sum(ep_type)} 個月")
print(f"CP 型: {np.sum(cp_type)} 個月")
print(f"轉捩點: {len(transitions)} 個")

# 測試模型喺唔同類型事件嘅表現
results = {
    'weak_events': {},
    'strong_events': {},
    'ep_type': {},
    'cp_type': {},
    'transitions': {}
}

# 用 ρ 門檻做預測
rho_threshold = np.percentile(rho, 40)
pred = (rho < rho_threshold).astype(int)

# 測試 1: 弱事件 vs 強事件
for event_type, events, name in [
    (weak_events, weak_events, 'weak_events'),
    (strong_events, strong_events, 'strong_events')
]:
    tp = np.sum((pred == 1) & (events == 1))
    fp = np.sum((pred == 1) & (events == 0))
    fn = np.sum((pred == 0) & (events == 1))
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    results[name] = {
        'tp': int(tp),
        'fp': int(fp),
        'fn': int(fn),
        'precision': float(precision),
        'recall': float(recall),
        'f1': float(f1),
        'n_events': int(np.sum(events))
    }
    
    print(f"\n{name}: F1 = {f1:.3f} (P={precision:.3f}, R={recall:.3f})")

# 測試 2: CP vs EP 型
for event_type, events, name in [
    (ep_type, ep_type, 'ep_type'),
    (cp_type, cp_type, 'cp_type')
]:
    tp = np.sum((pred == 1) & (events == 1))
    fp = np.sum((pred == 1) & (events == 0))
    fn = np.sum((pred == 0) & (events == 1))
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    results[name] = {
        'tp': int(tp),
        'fp': int(fp),
        'fn': int(fn),
        'precision': float(precision),
        'recall': float(recall),
        'f1': float(f1),
        'n_events': int(np.sum(events))
    }
    
    print(f"\n{name}: F1 = {f1:.3f} (P={precision:.3f}, R={recall:.3f})")

# 測試 3: 轉捩點預測
# 檢查轉捩點前 6 個月有冇 ρ 信號
transition_hits = 0
for trans_idx in transitions:
    # 檢查轉捩點前 6 個月
    window_start = max(0, trans_idx - 6)
    window_end = trans_idx
    if np.any(pred[window_start:window_end] == 1):
        transition_hits += 1

transition_recall = transition_hits / len(transitions) if len(transitions) > 0 else 0
results['transitions'] = {
    'n_transitions': len(transitions),
    'hits': transition_hits,
    'recall': float(transition_recall)
}

print(f"\n轉捩點預測: {transition_hits}/{len(transitions)} = {transition_recall:.3f}")

# 保存結果
output_file = Path('/app/working/workspaces/tygtDc/projects/enso/notes/2026-08-18-multi-weather-test.json')
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n結果已保存到: {output_file}")
