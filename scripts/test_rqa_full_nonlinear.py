#!/usr/bin/env python3
"""
Recurrence Quantification Analysis (RQA) for ENSO
完全非線性方法，唔降維，保留完整動力學結構

RQA 指標：
1. Recurrence Rate (RR): 系統返回相近狀態嘅概率
2. Determinism (DET): 可預測性（對角線結構比例）
3. Laminarity (LAM): 層流結構（垂直線比例）
4. Entropy: 複雜度
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter
import json
from pathlib import Path
import time

print("=== 載入數據 ===")
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

print(f"數據長度: {len(oni)} 個月")
print(f"El Niño 事件: {np.sum(el_nino)} 個月")

# === RQA 計算 ===
print("\n=== Recurrence Quantification Analysis ===")

def compute_recurrence_matrix(x, epsilon):
    """計算 recurrence matrix"""
    N = len(x)
    R = np.zeros((N, N), dtype=bool)
    
    for i in range(N):
        for j in range(i, N):
            if abs(x[i] - x[j]) < epsilon:
                R[i, j] = True
                R[j, i] = True
    
    return R

def compute_rqa_metrics(R):
    """計算 RQA 指標"""
    N = R.shape[0]
    
    # 1. Recurrence Rate (RR)
    rr = np.sum(R) / (N * N)
    
    # 2. Determinism (DET) - 對角線結構比例
    # 計算每條對角線嘅長度
    diag_lengths = []
    for k in range(-N+1, N):
        diag = np.diagonal(R, offset=k)
        if len(diag) > 0:
            # 計算連續 1 嘅長度
            consecutive = 0
            for val in diag:
                if val:
                    consecutive += 1
                else:
                    if consecutive >= 2:  # 至少 2 個連續點
                        diag_lengths.append(consecutive)
                    consecutive = 0
            if consecutive >= 2:
                diag_lengths.append(consecutive)
    
    if len(diag_lengths) > 0:
        det = np.sum(diag_lengths) / np.sum(R)
    else:
        det = 0
    
    # 3. Laminarity (LAM) - 垂直線比例
    vert_lengths = []
    for j in range(N):
        consecutive = 0
        for i in range(N):
            if R[i, j]:
                consecutive += 1
            else:
                if consecutive >= 2:
                    vert_lengths.append(consecutive)
                consecutive = 0
        if consecutive >= 2:
            vert_lengths.append(consecutive)
    
    if len(vert_lengths) > 0:
        lam = np.sum(vert_lengths) / np.sum(R)
    else:
        lam = 0
    
    # 4. Entropy - 對角線長度分佈嘅 entropy
    if len(diag_lengths) > 0:
        hist, _ = np.histogram(diag_lengths, bins=20)
        hist = hist[hist > 0]
        prob = hist / np.sum(hist)
        entropy = -np.sum(prob * np.log(prob))
    else:
        entropy = 0
    
    return {
        'RR': float(rr),
        'DET': float(det),
        'LAM': float(lam),
        'Entropy': float(entropy)
    }

# 用 ONI 做 RQA
print("\n--- RQA on ONI ---")
epsilon_oni = 0.5  # 閾值
print(f"計算 recurrence matrix (epsilon={epsilon_oni})...")
t0 = time.time()
R_oni = compute_recurrence_matrix(oni, epsilon_oni)
print(f"Recurrence matrix 計算完成: {time.time()-t0:.2f} 秒")

print("計算 RQA 指標...")
rqa_oni = compute_rqa_metrics(R_oni)
print(f"ONI RQA: RR={rqa_oni['RR']:.4f}, DET={rqa_oni['DET']:.4f}, LAM={rqa_oni['LAM']:.4f}, Entropy={rqa_oni['Entropy']:.4f}")

# 用 ρ 做 RQA
print("\n--- RQA on ρ ---")
epsilon_rho = np.std(rho) * 0.5  # 用標準差嘅 50% 做閾值
print(f"計算 recurrence matrix (epsilon={epsilon_rho:.4f})...")
t0 = time.time()
R_rho = compute_recurrence_matrix(rho, epsilon_rho)
print(f"Recurrence matrix 計算完成: {time.time()-t0:.2f} 秒")

print("計算 RQA 指標...")
rqa_rho = compute_rqa_metrics(R_rho)
print(f"ρ RQA: RR={rqa_rho['RR']:.4f}, DET={rqa_rho['DET']:.4f}, LAM={rqa_rho['LAM']:.4f}, Entropy={rqa_rho['Entropy']:.4f}")

# === RQA 特徵作為預測因子 ===
print("\n=== RQA 特徵作為預測因子 ===")

# 用 rolling window 計算 RQA 特徵
window_size = 12  # 12 個月窗口

rqa_features_oni = []
rqa_features_rho = []

for i in range(window_size, len(oni)):
    # ONI RQA
    oni_window = oni[i-window_size:i]
    R_window_oni = compute_recurrence_matrix(oni_window, epsilon_oni)
    rqa_oni_window = compute_rqa_metrics(R_window_oni)
    rqa_features_oni.append(rqa_oni_window)
    
    # ρ RQA
    rho_window = rho[i-window_size:i]
    R_window_rho = compute_recurrence_matrix(rho_window, epsilon_rho)
    rqa_rho_window = compute_rqa_metrics(R_window_rho)
    rqa_features_rho.append(rqa_rho_window)

rqa_features_oni = pd.DataFrame(rqa_features_oni)
rqa_features_rho = pd.DataFrame(rqa_features_rho)

print(f"RQA 特徵數量: {len(rqa_features_oni)}")

# 用 RQA 特徵做預測
# 目標：預測下個月嘅 El Niño
target = el_nino[window_size:]

# Walk-forward validation
window_train = 36  # 3 年訓練窗口
lead_times = [0, 3, 6]

results = {
    'rqa_oni': {},
    'rqa_rho': {},
    'rqa_combined': {}
}

for lead_time in lead_times:
    print(f"\n=== Lead Time: {lead_time} 個月 ===")
    
    # 準備目標
    target_shifted = el_nino[window_size + lead_time:]
    n_samples = len(target_shifted)
    
    # 準備特徵
    features_oni = rqa_features_oni[:n_samples]
    features_rho = rqa_features_rho[:n_samples]
    
    # Walk-forward validation
    scores_oni = []
    scores_rho = []
    scores_combined = []
    
    for start_idx in range(window_train, n_samples - 12, 12):
        end_idx = start_idx + 12
        
        # 訓練數據
        train_oni = features_oni[:start_idx]
        train_rho = features_rho[:start_idx]
        train_target = target_shifted[:start_idx]
        
        # 測試數據
        test_oni = features_oni[start_idx:end_idx]
        test_rho = features_rho[start_idx:end_idx]
        test_target = target_shifted[start_idx:end_idx]
        
        # 簡單門檻模型
        # ONI RQA: 用 RR 做預測（高 RR = 系統穩定 = 可能 El Niño）
        rr_threshold_oni = np.percentile(train_oni['RR'], 60)
        pred_oni = (test_oni['RR'] > rr_threshold_oni).astype(int)
        
        # ρ RQA: 用 DET 做預測（高 DET = 高可預測性 = 可能 El Niño）
        det_threshold_rho = np.percentile(train_rho['DET'], 60)
        pred_rho = (test_rho['DET'] > det_threshold_rho).astype(int)
        
        # Combined: 兩者都用
        pred_combined = ((test_oni['RR'] > rr_threshold_oni) & (test_rho['DET'] > det_threshold_rho)).astype(int)
        
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
                'f1': float(f1)
            }
        
        scores_oni.append(calc_f1(pred_oni, test_target))
        scores_rho.append(calc_f1(pred_rho, test_target))
        scores_combined.append(calc_f1(pred_combined, test_target))
    
    results['rqa_oni'][f'lead_{lead_time}m'] = {
        'avg_precision': np.mean([s['precision'] for s in scores_oni]),
        'avg_recall': np.mean([s['recall'] for s in scores_oni]),
        'avg_f1': np.mean([s['f1'] for s in scores_oni]),
        'std_f1': np.std([s['f1'] for s in scores_oni])
    }
    
    results['rqa_rho'][f'lead_{lead_time}m'] = {
        'avg_precision': np.mean([s['precision'] for s in scores_rho]),
        'avg_recall': np.mean([s['recall'] for s in scores_rho]),
        'avg_f1': np.mean([s['f1'] for s in scores_rho]),
        'std_f1': np.std([s['f1'] for s in scores_rho])
    }
    
    results['rqa_combined'][f'lead_{lead_time}m'] = {
        'avg_precision': np.mean([s['precision'] for s in scores_combined]),
        'avg_recall': np.mean([s['recall'] for s in scores_combined]),
        'avg_f1': np.mean([s['f1'] for s in scores_combined]),
        'std_f1': np.std([s['f1'] for s in scores_combined])
    }
    
    print(f"RQA ONI: F1 = {results['rqa_oni'][f'lead_{lead_time}m']['avg_f1']:.3f} ± {results['rqa_oni'][f'lead_{lead_time}m']['std_f1']:.3f}")
    print(f"RQA ρ: F1 = {results['rqa_rho'][f'lead_{lead_time}m']['avg_f1']:.3f} ± {results['rqa_rho'][f'lead_{lead_time}m']['std_f1']:.3f}")
    print(f"RQA Combined: F1 = {results['rqa_combined'][f'lead_{lead_time}m']['avg_f1']:.3f} ± {results['rqa_combined'][f'lead_{lead_time}m']['std_f1']:.3f}")

# 保存結果
output = {
    'rqa_metrics_oni': rqa_oni,
    'rqa_metrics_rho': rqa_rho,
    'prediction_results': results
}

output_file = Path('/app/working/workspaces/tygtDc/projects/enso/notes/2026-08-18-rqa-analysis.json')
with open(output_file, 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n結果已保存到: {output_file}")
