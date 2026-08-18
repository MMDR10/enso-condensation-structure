#!/usr/bin/env python3
"""
Test C: 同傳統 ENSO 預測模型比較
比較對象：
1. Persistence model (假設下個月同今個月一樣)
2. Climatology model (用歷史平均)
3. 簡單線性回歸 (Nino3.4 自回歸)
4. 我哋嘅非線性模型 (ρ 雜訊)
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter
from sklearn.linear_model import LinearRegression
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

# Nino3.4 index (ONI)
oni = np.nanmean(nino34_anomaly, axis=1)

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

# 定義 El Niño 事件
oni_threshold = 0.5
el_nino = (oni > oni_threshold).astype(int)

print(f"=== 模型比較 ===")
print(f"數據長度: {len(oni)} 個月")
print(f"El Niño 事件: {np.sum(el_nino)} 個月 ({np.sum(el_nino)/len(el_nino)*100:.1f}%)")

# Walk-forward validation
window_size = 36  # 3 年訓練窗口
lead_times = [0, 3, 6]

results = {
    'persistence': {},
    'climatology': {},
    'linear_regression': {},
    'nonlinear_rho': {}
}

for lead_time in lead_times:
    print(f"\n=== Lead Time: {lead_time} 個月 ===")
    
    # 準備目標
    target = np.zeros(len(oni) - lead_time)
    for i in range(len(target)):
        target[i] = 1 if oni[i + lead_time] > oni_threshold else 0
    
    n_samples = len(target)
    
    # 準備特徵
    oni_feature = oni[:n_samples]
    rho_feature = rho[:n_samples]
    
    # Walk-forward validation
    model_scores = {
        'persistence': [],
        'climatology': [],
        'linear_regression': [],
        'nonlinear_rho': []
    }
    
    for start_idx in range(window_size, n_samples - 12, 12):
        end_idx = start_idx + 12
        
        # 訓練數據
        train_oni = oni_feature[:start_idx]
        train_rho = rho_feature[:start_idx]
        train_target = target[:start_idx]
        
        # 測試數據
        test_oni = oni_feature[start_idx:end_idx]
        test_rho = rho_feature[start_idx:end_idx]
        test_target = target[start_idx:end_idx]
        
        # 模型 1: Persistence (假設下個月同今個月一樣)
        if lead_time == 0:
            pred_persistence = (test_oni > oni_threshold).astype(int)
        else:
            # 對於 lead time > 0，persistence 假設未來同現在一樣
            pred_persistence = (test_oni > oni_threshold).astype(int)
        
        # 模型 2: Climatology (用歷史平均)
        # 計算訓練數據中每個月的 El Niño 概率
        clim_prob = np.zeros(12)
        for month in range(12):
            month_indices = [i for i in range(len(train_target)) if i % 12 == month]
            if len(month_indices) > 0:
                clim_prob[month] = np.mean(train_target[month_indices])
        
        # 用 climatology 預測
        pred_climatology = np.zeros(len(test_target))
        for i in range(len(test_target)):
            month = (start_idx + i) % 12
            pred_climatology[i] = 1 if clim_prob[month] > 0.5 else 0
        
        # 模型 3: 線性回歸 (ONI 自回歸)
        if lead_time == 0:
            # Nowcasting: 用當前 ONI 預測
            pred_linear = (test_oni > oni_threshold).astype(int)
        else:
            # Forecasting: 用 AR(1) 模型
            # 訓練 AR(1) 模型: ONI(t+lead) = a * ONI(t) + b
            X_train = train_oni[:-lead_time].reshape(-1, 1)
            y_train = train_oni[lead_time:]
            
            if len(X_train) > 0:
                model = LinearRegression()
                model.fit(X_train, y_train)
                
                # 預測
                X_test = test_oni.reshape(-1, 1)
                oni_pred = model.predict(X_test)
                pred_linear = (oni_pred > oni_threshold).astype(int)
            else:
                pred_linear = np.zeros(len(test_target))
        
        # 模型 4: 非線性 ρ 模型
        rho_threshold = np.percentile(train_rho, 40)
        pred_nonlinear = (test_rho < rho_threshold).astype(int)
        
        # 計算 F1 score
        def calc_f1(pred, true):
            tp = np.sum((pred == 1) & (true == 1))
            fp = np.sum((pred == 1) & (true == 0))
            fn = np.sum((pred == 0) & (true == 1))
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            return {
                'precision': float(precision),
                'recall': float(recall),
                'f1': float(f1),
                'tp': int(tp),
                'fp': int(fp),
                'fn': int(fn)
            }
        
        model_scores['persistence'].append(calc_f1(pred_persistence, test_target))
        model_scores['climatology'].append(calc_f1(pred_climatology, test_target))
        model_scores['linear_regression'].append(calc_f1(pred_linear, test_target))
        model_scores['nonlinear_rho'].append(calc_f1(pred_nonlinear, test_target))
    
    # 平均結果
    for model_name in model_scores:
        results[model_name][f'lead_{lead_time}m'] = {
            'avg_precision': np.mean([s['precision'] for s in model_scores[model_name]]),
            'avg_recall': np.mean([s['recall'] for s in model_scores[model_name]]),
            'avg_f1': np.mean([s['f1'] for s in model_scores[model_name]]),
            'std_f1': np.std([s['f1'] for s in model_scores[model_name]]),
            'n_tests': len(model_scores[model_name])
        }
        
        print(f"{model_name}: F1 = {results[model_name][f'lead_{lead_time}m']['avg_f1']:.3f} ± {results[model_name][f'lead_{lead_time}m']['std_f1']:.3f}")

# 保存結果
output_file = Path('/app/working/workspaces/tygtDc/projects/enso/notes/2026-08-18-model-comparison.json')
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n結果已保存到: {output_file}")

# 總結
print("\n=== 總結 ===")
print("Lead 0 個月 (Nowcasting):")
print(f"  Persistence: {results['persistence']['lead_0m']['avg_f1']:.3f}")
print(f"  Climatology: {results['climatology']['lead_0m']['avg_f1']:.3f}")
print(f"  Linear Regression: {results['linear_regression']['lead_0m']['avg_f1']:.3f}")
print(f"  Nonlinear ρ: {results['nonlinear_rho']['lead_0m']['avg_f1']:.3f}")

print("\nLead 3 個月:")
print(f"  Persistence: {results['persistence']['lead_3m']['avg_f1']:.3f}")
print(f"  Climatology: {results['climatology']['lead_3m']['avg_f1']:.3f}")
print(f"  Linear Regression: {results['linear_regression']['lead_3m']['avg_f1']:.3f}")
print(f"  Nonlinear ρ: {results['nonlinear_rho']['lead_3m']['avg_f1']:.3f}")

print("\nLead 6 個月:")
print(f"  Persistence: {results['persistence']['lead_6m']['avg_f1']:.3f}")
print(f"  Climatology: {results['climatology']['lead_6m']['avg_f1']:.3f}")
print(f"  Linear Regression: {results['linear_regression']['lead_6m']['avg_f1']:.3f}")
print(f"  Nonlinear ρ: {results['nonlinear_rho']['lead_6m']['avg_f1']:.3f}")
