"""Pre-extract ETOPO1 bathymetry/topography over the ENSO memory-map domain
(lat -30..30, lon 120..290) at 0.25deg, matching rho_seasonal_memory_maps.npz grid.
Outputs small npy files for downstream correlation (no need for full 2GB ETOPO).
"""
import numpy as np, xarray as xr, os

DS_PATH = 'output/historical/data/etopo1_bed.grd'
OUT_DIR = 'output/historical/data/'
os.makedirs(OUT_DIR, exist_ok=True)

# memory map grid
LAT0, LAT1, LON0, LON1, STEP = -30.0, 30.0, 120.0, 290.0, 0.25
nlat = int((LAT1 - LAT0) / STEP)  # 240
nlon = int((LON1 - LON0) / STEP)  # 680

# ETOPO1: 1 arc-min. y: -90..90 (10801), x: -180..180 (21601)
# lat index: i = (lat+90)*60 ; lon index: j = (lon+180)*60
# We need lat -30..30 -> y 3600..5401 ; lon 120..290 wraps: 120..180 -> x 18000..21601,
# 180..290 -> x 0..6601
y0, y1 = int((LAT0 + 90) * 60), int((LAT1 + 90) * 60) + 1   # 3600, 5401
xa0, xa1 = int((120 + 180) * 60), int((180 + 180) * 60) + 1  # 18000, 21601
xb0, xb1 = int((-180 + 180) * 60), int((-70 + 180) * 60) + 1  # 0, 6601

print(f"ETOPO slice: y[{y0}:{y1}] ({y1-y0} pts), x[{xa0}:{xa1}] + x[{xb0}:{xb1}]")

ds = xr.open_dataset(DS_PATH, chunks={'x': 21601, 'y': 300})
z = ds['z']

# Build wrapped lon coordinate array for the whole needed band
za = z.isel(y=slice(y0, y1), x=slice(xa0, xa1))   # lon 120..180 (east)
zb = z.isel(y=slice(y0, y1), x=slice(xb0, xb1))   # lon -180..-70 (=180..290E)
zband = xr.concat([za, zb], dim='x')
print("zband shape:", zband.shape)

# coarsen factor: 0.25deg / (1/60 deg) = 15
F = 15
zmean = zband.coarsen(x=F, y=F, boundary='trim').mean().compute()
zstd = zband.coarsen(x=F, y=F, boundary='trim').std().compute()
print("coarsened:", zmean.shape)

# lat/lon arrays
lats = np.linspace(LAT0 + STEP/2, LAT1 - STEP/2, nlat)
lons = np.linspace(LON0 + STEP/2, LON1 - STEP/2, nlon)

# gradient magnitude on the 0.25deg grid
gx, gy = np.gradient(zmean.values, STEP, STEP)
grad = np.hypot(gx, gy)

np.save(os.path.join(OUT_DIR, 'etopo_0p25_mean.npy'), zmean.values)
np.save(os.path.join(OUT_DIR, 'etopo_0p25_std.npy'), zstd.values)
np.save(os.path.join(OUT_DIR, 'etopo_0p25_grad.npy'), grad)
np.save(os.path.join(OUT_DIR, 'etopo_0p25_lat.npy'), lats)
np.save(os.path.join(OUT_DIR, 'etopo_0p25_lon.npy'), lons)

print("Saved etopo_0p25_mean/std/grad/lat/lon")
print("mean elev range:", np.nanmin(zmean.values), np.nanmax(zmean.values))
print("land cells (mean>-50):", int((zmean.values > -50).sum()), "/", zmean.values.size)
