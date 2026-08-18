"""ENSO memory map x ETOPO topographic anchoring test.
Question: is the seasonal memory field (rho net memory) anchored to static
bathymetric features (EPR ridge / Galapagos platform / continental shelves)?
Prediction: memory strong where topography is anomalous (shallow ridges,
steep gradients), weak over deep abyssal plains.
"""
import os
import numpy as np
from scipy.stats import spearmanr

D = os.path.join(os.path.dirname(__file__), '..', 'data') + '/'
mem = np.load(os.path.join(os.path.dirname(__file__), '..', 'data', 'rho_seasonal_memory_maps.npz'))
net = mem['net']            # net memory (12-month recurrence - random baseline?)
corr_yr = mem['corr_yr']
corr_adj = mem['corr_adj']
corr_rnd = mem['corr_rnd']
rho = mem['rho_amp']

zmean = np.load(D + 'etopo_0p25_mean.npy')
zstd = np.load(D + 'etopo_0p25_std.npy')
grad = np.load(D + 'etopo_0p25_grad.npy')
lats = np.load(D + 'etopo_0p25_lat.npy')
lons = np.load(D + 'etopo_0p25_lon.npy')

# ocean mask: memory map valid (non-nan) & deep ocean (ETOPO mean < -50 m)
ocean = (~np.isnan(net)) & (zmean < -50)
print(f"ocean cells: {ocean.sum()}/{net.size} ({ocean.sum()/net.size*100:.1f}%)")

def report(name, target):
    t = target[ocean]
    # 1) simple Spearman with topography features
    for feat_name, feat in [('elev', zmean), ('rough', zstd), ('grad', grad)]:
        f = feat[ocean]
        r, p = spearmanr(t, f)
        print(f"  [{name}] spearman memory~{feat_name}: r={r:+.3f} (p={p:.2e})")

print("== simple spearman (ocean only) ==")
report('net', net)
report('corr_yr', corr_yr)
report('corr_adj', corr_adj)
report('corr_rnd', corr_rnd)
report('rho_amp', rho)

# 2) East-end (270-290E) vs Nino3.4 (190-240E) topographic comparison
print("\n== regional topography ==")
def region(lon_a, lon_b, lat_a=-5, lat_b=5):
    lonm, latm = np.meshgrid(lons, lats)
    m = (lonm >= lon_a) & (lonm <= lon_b) & (latm >= lat_a) & (latm <= lat_b)
    return m
east = region(270, 290) & ocean
nino = region(190, 240) & ocean
print(f"east-end cells: {east.sum()}, nino3.4 cells: {nino.sum()}")
print(f"  east-end: net={np.nanmean(net[east]):.4f}, elev={np.nanmean(zmean[east]):.0f}m, rough={np.nanmean(zstd[east]):.0f}, grad={np.nanmean(grad[east]):.1f}")
print(f"  nino3.4: net={np.nanmean(net[nino]):.4f}, elev={np.nanmean(zmean[nino]):.0f}m, rough={np.nanmean(zstd[nino]):.0f}, grad={np.nanmean(grad[nino]):.1f}")

# 3) density-weighted: top-10% memory cells vs bottom-10% memory cells topography
print("\n== top vs bottom memory decile topography ==")
t = net[ocean]
thr_hi = np.nanpercentile(t, 90)
thr_lo = np.nanpercentile(t, 10)
hi = ocean & (net >= thr_hi)
lo = ocean & (net <= thr_lo)
print(f"top decile: n={hi.sum()}, bottom decile: n={lo.sum()}")
for feat_name, feat in [('elev', zmean), ('rough', zstd), ('grad', grad)]:
    fh = feat[hi]; fl = feat[lo]
    print(f"  {feat_name}: top={np.nanmean(fh):.1f} vs bottom={np.nanmean(fl):.1f} -> diff={np.nanmean(fh)-np.nanmean(fl):+.1f}")

# 4) null: shuffle memory field spatially, recompute spearman vs elev (200 reps)
print("\n== null (shuffle memory, 200 reps) ==")
rng = np.random.default_rng(0)
idx_ocean = np.argwhere(ocean)
obs_r, _ = spearmanr(net[ocean], zmean[ocean])
obs_rg, _ = spearmanr(net[ocean], grad[ocean])
null_r, null_rg = [], []
for _ in range(200):
    perm = rng.permutation(len(idx_ocean))
    t_shuf = net[ocean][perm]
    null_r.append(spearmanr(t_shuf, zmean[ocean]).statistic)
    null_rg.append(spearmanr(t_shuf, grad[ocean]).statistic)
null_r = np.array(null_r); null_rg = np.array(null_rg)
print(f"  memory~elev: obs={obs_r:+.3f}, null={null_r.mean():+.3f}+-{null_r.std():.3f}, z={(obs_r-null_r.mean())/null_r.std():.2f}")
print(f"  memory~grad: obs={obs_rg:+.3f}, null={null_rg.mean():+.3f}+-{null_rg.std():.3f}, z={(obs_rg-null_rg.mean())/null_rg.std():.2f}")
