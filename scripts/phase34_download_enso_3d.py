#!/usr/bin/env python3
"""下載熱帶太平洋 ERA5 8 層 geopotential + u/v（同 Greenland 平行）"""
import cdsapi, os

OUT_DIR = '/tmp/enso_saddle_3d'
os.makedirs(OUT_DIR, exist_ok=True)

# 2020（La Niña 年尾）+ 2015（強 El Niño）+ 2016（El Niño→中性）
for var, name in [('geopotential', 'z'), ('u_component_of_wind', 'u'), ('v_component_of_wind', 'v')]:
    for year in ['2015', '2016', '2020']:
        out = os.path.join(OUT_DIR, f'era5_{name}_{year}_8lev.nc')
        if os.path.exists(out):
            print("已存在:", out)
            continue
        c = cdsapi.Client()
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'variable': [var],
                'pressure_level': ['850', '700', '600', '500', '400', '300', '250', '200'],
                'year': [year],
                'month': ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'],
                'day': ['01'],
                'time': '00:00',
                'area': [30.0, 120.0, -30.0, 290.0],
                'grid': [0.5, 0.5],
                'format': 'netcdf',
            },
            out)
        print("✅ 下載完成:", out)
