#!/usr/bin/env python3
"""
非線性 vs 線性模型比較測試
Test A: 線性特徵 (Nino3.4 SST) vs 非線性特徵 (ρ 雜訊密度)
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

# 計算 Nino3.4 SST anomaly (線性特徵)
lat_mask_nino = (lat >= -5) & (lat <= 5)
lon_mask_nino = (lon >= 190) & (lon <= 240)
nino34_sst = sst[:, lat_mask_nino, :][:, :, lon_mask_nino]
nino34_sst = nino34_sst.reshape(sst.shape[0], -1)

# 計算 climatology (每月嘅平均值)
climatology = np.zeros((12, nino34_sst.shape[1]))
for month in range(12):
    climatology[month] = np.nanmean(nino34_sst[month::12], axis=0)

# 計算 anomaly
nino34_anomaly = np.zeros_like(nino34_sst)
for t in range(nino34_sst.shape[0]):
    month = t % 12
    nino34_anomaly[t] = nino34_sst[t] - climatology[month]

# 計算 Nino3.4 index (區域平均 anomaly)
nino34_mean = np.nanmean(nino34_anomaly, axis=1)

# 計算 ρ 雜訊密度 (非線性特徵)
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

# 定義 El Niño 事件 (ONI > 0.5)
oni_threshold = 0.5
el_nino_events = nino34_mean > oni_threshold

# Walk-forward validation (滾動窗口)
window_size = 36  # 3 年訓練窗口
lead_times = [0, 3, 6]  # 預測 0, 3, 6 個月後

results = {
    'linear': {},
    'nonlinear': {},
    'combined': {}
}

for lead_time in lead_times:
    print(f"\n=== Lead Time: {lead_time} 個月 ===")
    
    # 準備目標變量
    target = np.zeros(len(nino34_mean) - lead_time)
    for i in range(len(target)):
        target[i] = 1 if nino34_mean[i + lead_time] > oni_threshold else 0
    
    print(f"Target 分佈: {np.sum(target)} El Niño / {len(target) - np.sum(target)} Non-El Niño")
    
    # 準備特徵 (需要延遲 lead_time 個月)
    n_samples = len(target)
    
    # 線性模型: Nino3.4 SST
    linear_feature = nino34_mean[:n_samples]
    
    # 非線性模型: ρ 雜訊密度
    nonlinear_feature = rho[:n_samples]
    
    print(f"特徵數量: {n_samples}")
    
    # Walk-forward validation
    linear_scores = []
    nonlinear_scores = []
    combined_scores = []
    
    test_count = 0
    for start_idx in range(window_size, n_samples - 12, 12):  # 每年預測一次
        end_idx = start_idx + 12
        test_count += 1
        
        # 訓練數據
        train_linear = linear_feature[:start_idx]
        train_nonlinear = nonlinear_feature[:start_idx]
        train_target = target[:start_idx]
        
        # 測試數據
        test_linear = linear_feature[start_idx:end_idx]
        test_nonlinear = nonlinear_feature[start_idx:end_idx]
        test_target = target[start_idx:end_idx]
        
        # 線性模型: 簡單門檻 (基於訓練數據的百分位數)
        threshold_linear = np.percentile(train_linear, 60)
        pred_linear = (test_linear > threshold_linear).astype(int)
        
        # 非線性模型: ρ 門檻 (基於訓練數據的百分位數)
        threshold_nonlinear = np.percentile(train_nonlinear, 40)
        pred_nonlinear = (test_nonlinear < threshold_nonlinear).astype(int)
        
        # 組合模型: 兩者都用
        pred_combined = ((test_linear > threshold_linear) & (test_nonlinear < threshold_nonlinear)).astype(int)
        
        # 計算 F1 score
        def calc_f1(pred, true):
            tp = np.sum((pred == 1) & (true == 1))
            fp = np.sum((pred == 1) & (true == 0))
            fn = np.sum((pred == 0) & (true == 1))
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            return {
                'precision': precision,
                'recall': recall,
                'f1': f1,
                'tp': int(tp),
                'fp': int(fp),
                'fn': int(fn)
            }
        
        linear_scores.append(calc_f1(pred_linear, test_target))
        nonlinear_scores.append(calc_f1(pred_nonlinear, test_target))
        combined_scores.append(calc_f1(pred_combined, test_target))
    
    print(f"Walk-forward 測試次數: {test_count}")
    print(f"線性模型 F1 scores: {[s['f1'] for s in linear_scores[:5]]} ...")
    print(f"非線性模型 F1 scores: {[s['f1'] for s in nonlinear_scores[:5]]} ...")
    
    # 平均結果
    results['linear'][f'lead_{lead_time}m'] = {
        'avg_precision': np.mean([s['precision'] for s in linear_scores]),
        'avg_recall': np.mean([s['recall'] for s in linear_scores]),
        'avg_f1': np.mean([s['f1'] for s in linear_scores]),
        'std_f1': np.std([s['f1'] for s in linear_scores]),
        'n_tests': len(linear_scores)
    }
    
    results['nonlinear'][f'lead_{lead_time}m'] = {
        'avg_precision': np.mean([s['precision'] for s in nonlinear_scores]),
        'avg_recall': np.mean([s['recall'] for s in nonlinear_scores]),
        'avg_f1': np.mean([s['f1'] for s in nonlinear_scores]),
        'std_f1': np.std([s['f1'] for s in nonlinear_scores]),
        'n_tests': len(nonlinear_scores)
    }
    
    results['combined'][f'lead_{lead_time}m'] = {
        'avg_precision': np.mean([s['precision'] for s in combined_scores]),
        'avg_recall': np.mean([s['recall'] for s in combined_scores]),
        'avg_f1': np.mean([s['f1'] for s in combined_scores]),
        'std_f1': np.std([s['f1'] for s in combined_scores]),
        'n_tests': len(combined_scores)
    }
    
    print(f"線性模型 (Nino3.4 SST): F1 = {results['linear'][f'lead_{lead_time}m']['avg_f1']:.3f} ± {results['linear'][f'lead_{lead_time}m']['std_f1']:.3f}")
    print(f"非線性模型 (ρ 雜訊): F1 = {results['nonlinear'][f'lead_{lead_time}m']['avg_f1']:.3f} ± {results['nonlinear'][f'lead_{lead_time}m']['std_f1']:.3f}")
    print(f"組合模型: F1 = {results['combined'][f'lead_{lead_time}m']['avg_f1']:.3f} ± {results['combined'][f'lead_{lead_time}m']['std_f1']:.3f}")

# 保存結果
output_file = Path('/app/working/workspaces/tygtDc/projects/enso/notes/2026-08-18-linear-vs-nonlinear-test.json')
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n結果已保存到: {output_file}")
