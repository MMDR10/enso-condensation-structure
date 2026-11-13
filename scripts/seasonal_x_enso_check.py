#!/usr/bin/env python3
"""Seasonal evolution × ENSO phase check — enso-watch 2026-07→10 vs paper claims.
Inputs: saddle_history.json / curl_history.json (enso-watch) + CPC ONI."""
import json

import numpy as np

HERE = "."
S = json.load(open(f"{HERE}/saddle_history.json"))["records"]
C = json.load(open(f"{HERE}/curl_history.json"))["records"]
# CPC detrend.nino34 monthly anomaly（2026）
ONI = {"2026-07": 1.78, "2026-08": 2.17, "2026-09": 2.55, "2026-10": 2.68}
ONI3 = {"2026-07": 2.16, "2026-08": 2.46, "2026-09": 2.46}  # ONI 3-mo running (JAS/ASO/SON*)


def table(recs, name, extra):
    print(f"\n== {name} ==")
    print(f"{'month':<9}{'ONI':>6}{'saddle_frac':>12}{'dens_' + extra:>10}{'east_end':>10}{'nino34':>9}{'gap':>8}")
    for r in recs:
        m = r["month"]
        gap = r["east_end_density"] - r["nino34_density"]
        dkey = "spearman_density_rho" if "spearman_density_rho" in r else "spearman_density_curl"
        print(f"{m:<9}{ONI.get(m, float('nan')):>6.2f}{r['saddle_frac']:>12.5f}"
              f"{r[dkey]:>10.4f}{r['east_end_density']:>10.5f}{r['nino34_density']:>9.5f}{gap:>+8.4f}")


table(S, "OISST sst_anomaly field (4 months)", "rho")
table(C, "curl_tau interface field (3 months)", "curl")

print("\n== Claim checks ==")
# C1: phase-insensitive saddle fraction
for name, recs in [("OISST", S), ("curl", C)]:
    v = np.array([r["saddle_frac"] for r in recs])
    print(f"C1 saddle_frac {name}: mean={v.mean():.5f} range=[{v.min():.5f},{v.max():.5f}] "
          f"span={100*(v.max()-v.min()):.2f}pp  (paper: 0.564 phase-insensitive KW p=0.137)")

# C2: edge dominance east>nino34
for name, recs in [("OISST", S), ("curl", C)]:
    gaps = [r["east_end_density"] - r["nino34_density"] for r in recs]
    n_win = sum(1 for g in gaps if g > 0)
    print(f"C2 edge dominance {name}: east>nino34 in {n_win}/{len(gaps)} months, gaps={[round(g,4) for g in gaps]}")

# C3: position insensitive to ENSO — does the center density rise with ONI?
for name, recs in [("OISST", S), ("curl", C)]:
    xs = np.array([ONI[r["month"]] for r in recs if r["month"] in ONI])
    cen = np.array([r["nino34_density"] for r in recs if r["month"] in ONI])
    east = np.array([r["east_end_density"] for r in recs if r["month"] in ONI])
    rc = float(np.corrcoef(xs, cen)[0, 1]) if len(xs) > 2 else float("nan")
    re_ = float(np.corrcoef(xs, east)[0, 1]) if len(xs) > 2 else float("nan")
    print(f"C3 pos-insensitive {name}: corr(ONI, nino34_density)={rc:+.3f}  corr(ONI, east_end)={re_:+.3f}"
          f"  (centre falling / edge flat = structure NOT pulled to centre)")

# C4: La Niña strongest → El Niño ramp should weaken condensation proxies
for name, recs, dkey in [("OISST", S, "spearman_density_rho"), ("curl", C, "spearman_density_curl")]:
    vals = [r[dkey] for r in recs]
    half = len(vals) // 2
    a, b = float(np.mean(vals[:half])), float(np.mean(vals[half:]))
    print(f"C4 condensation-proxy {name} ({dkey}): first-half={a:.4f} last-half={b:.4f} "
          f"Δ={100*(b-a)/a:+.0f}%  (El Niño ramp; expect decline if La-Niña-strongest holds)")
