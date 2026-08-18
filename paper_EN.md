# Condensation Structure of ENSO: Seasonal Memory, Topographic Anchoring, and Phase-Transition Geometry

**Author:** tygtDc, Deep Research
**Date:** 2026-08-04
**Version:** v1.0 (draft)
**DOI:** pending
**Predecessor:** ENSO Simplified Prediction Note (10.5281/zenodo.21626908), ENSO Watch Dashboard (github.com/MMDR10/enso-watch)

---

## Abstract

ENSO is traditionally understood as a large-scale ocean–atmosphere coupled oscillation, with mainstream methods (GCMs) achieving only ~50–60% accuracy at 6-month leads. This paper proposes a different measurement route: **instead of predicting the ENSO index, measure the spatial condensation structure of the equatorial Pacific SST density field** — the geometric folding of seasonal frontal systems, their positional memory, and the topographic anchoring of that memory.

Using 39 years of OISST 0.25° data (1982–2020), we find:

1. **ρ-field condensation is the norm (D=1.334±0.016, 60/60 months, z=−46.4)**: the high-value cluster dimension of the SST seasonal activity density field (ρ) is far below random (null 1.73) — frontal systems **fold into near-1D line-like structures**, not uniform scatter. The condensation signature is 4× stronger at 0.25° than at 1° (z −46 vs −11).
2. **Temporal memory exists but is weak (adjacent-month Jaccard 1.12× cross-year, p=0.004; 86% of positions reshuffle)**: clusters have almost no temporal identity but retain weak positional memory — "time does not exist" is refined to "noise clusters have almost no temporal identity, but weak positional memory remains."
3. **Seasonal memory is dominated by the Pacific margins ("strong at both ends, weak in the center")**: east-end 270–290°E (American coast) net memory 0.096–0.117, 80.7% positive grid cells, 23.3% of top-5% share; Nino3.4 center weakest (6.4%). Seasonal memory is not a center phenomenon but a margin frontal-zone phenomenon.
4. **Memory position is insensitive to ENSO phase (p=0.28)**: the position is anchored by non-ENSO factors.
5. **Topographic anchoring (ETOPO decisive test)**: the strongest memory region sits on the EPR/Galapagos shallow ridge (−2636 m vs −4629 m abyssal plain; roughness/gradient 2.5×); memory intensity × topographic gradient positively correlated (r=+0.038, shuffle null **z=+13.9**), robust after excluding shelf effects (r=+0.040, p=5.6e-49). **Dynamic volcanic activity is excluded (EPR three-line NULL: temporal r<0.15, high/low volcanic months p=0.78/0.87, spatial r=−0.29)** — supporting "static submarine geometry" rather than "volcanic eruption."
6. **La Niña condenses most strongly (p=0.0011)**: in the cold phase the cold tongue deepens → the front steepens → stronger condensation, complementary to the superficial asymmetry of the three-number formula ("La Niña signal is weak").

Conclusion: the predictable part of ENSO is not the center index but the **condensation structure of marginal frontal systems** — recurring at the same location each year, anchored by fixed submarine topography. A complementary saddle-point (Morse) topology test (2D OISST 528 months + 3D ERA5 + wind-stress-curl interface) finds **no saddle ring** in ENSO (field fraction ≈ 0.564, phase-insensitive KW p = 0.137; fair-control z = −2.0) — saddle rings are vortex-specific — yet saddle points **co-localize with the condensation bands** (saddle density × activity density grid r = 0.42, z = +85; × |curl τ| at the air–sea interface r = 0.24–0.28, z = +141 to +222). This structure provides a measurement window orthogonal to mainstream ENSO prediction.

---

## 1. Introduction: Using the Wrong Ruler?

Mainstream ENSO prediction uses GCMs to simulate the coupled ocean–atmosphere system, achieving ~50–60% accuracy at 6-month leads. Prior work in this series (ENSO Simplified Prediction Note) proposed that ENSO is fundamentally a **phase transition**: one need not simulate every water molecule, only measure the critical point. That note used a three-number formula (x, ẋ, V) to capture the phase-transition state, and performed a five-dimensional scan with the Ô-HAT framework, finding ENSO's structural rigidity 3.0±1.8 (one of the most rigid natural systems ever measured).

But the three-number formula measures the ENSO **center** (indices such as Nino3.4). This paper asks: **where does the "structure" of the ENSO system actually live?** The answer points to a location mainstream work ignores — **the spatial folding of seasonal frontal systems**.

We apply the core operation of the Ô framework: take the SST density field ρ (seasonal activity intensity), measure the fractal dimension D of its high-value clusters (degree of folding), its cross-year recurrence (seasonal memory), and the fixed source of its spatial position (anchoring).

---

## 2. Data and Methods

### 2.1 Data

| Data | Resolution | Period | Use |
|------|-----------|--------|-----|
| NOAA OISST v2.1 | 0.25° | 1982–2020 (39 yr) | Primary analysis |
| ERA5 SST | 0.25° | 1982–2018 | Cross-validation |
| HadISST | 1° | 1950–2024 | Resolution contrast |
| NOAA ETOPO1 Bedrock | 1-arcmin → 0.25° | static | Topographic anchoring |
| USGS earthquake catalog + volcano data | — | 1970–2025 | Volcanic activity exclusion (EPR) |

### 2.2 Core Measures

- **ρ field**: per-grid-cell SST seasonal activity intensity (climatological std or anomaly amplitude).
- **D_fold (condensation dimension)**: box-counting dimension of ρ > p95 grid cells. Random scatter ≈ 2, line-like folding → 1. Null uses N_NULL=8, seed=42.
- **Seasonal memory**: per-grid-cell cross-year (same month, 12 months apart) SST corr minus random-baseline corr (152 pairs).
- **Saddle-point (Morse) structure**: 2D Hessian eigenvalue sign classification of the SST anomaly field (and of the wind-stress-curl field); saddle fraction per month plus spatial co-location of saddle density × activity density. A **saddle ring** — a spatially uniform ring of saddle points around a central vortex — is the signature previously established for tropical cyclones (eyewall shear band); here we test whether the ENSO field (which has no central vortex) exhibits ring-like saddle structure, and where its saddle points actually sit.
- **Topographic features**: ETOPO 0.25° mean elevation, std (roughness), np.gradient magnitude (gradient).
- **Volcanic activity**: USGS earthquake energy density, volcanic eruption months, three independent lines (temporal coupling / high-vs-low volcanic months / spatial density).

### 2.3 Validation Standards

- Shuffle null (spatial/temporal randomization) 200–500 reps, z-score.
- Dual resolution (1° vs 0.25°) and dual product (OISST vs ERA5) cross-validation.
- Shelf-effect exclusion (>1000 m deep-water mask) robustness check.

---

## 3. Main Findings

### 3.1 Condensation Is the Norm, Not Event-Driven

| Measure | Value | Statistics |
|---------|-------|-----------|
| ρ>p95 cluster dimension D (0.25°, 60 months) | **1.334 ± 0.016** | negative in 60/60 months, mean z=−46.4 |
| null dimension | ~1.73 | — |
| 1° contrast (HadISST) | D=1.04, z=−11.2 | 0.25° condensation signature 4× stronger |

**Interpretation:** the high-value clusters of the SST seasonal activity density field are far lower-dimensional than random — frontal systems fold into near-1D line-like ("pipeline") structures. Every one of 60 months is significant, i.e., condensation is a normal-state structure, not an accidental product of a particular ENSO event.

### 3.2 Temporal Memory: Present but Extremely Weak

| Contrast | Jaccard | Test |
|----------|---------|------|
| Adjacent months | 0.140 | weak but significant |
| Same month across years | 0.126 | — |
| Ratio | **1.12×** | t=2.91, p=0.004 |

86% of positions reshuffle monthly — clusters have almost no temporal identity; but the residual 1.12× memory is real. Higher resolution (0.25° vs 1°) makes memory more visible. **The precise version of "time does not exist": noise clusters have almost no temporal identity, but weak positional memory remains.**

### 3.3 Seasonal Memory: "Strong at Both Ends, Weak in the Center"

| Region | Net memory | >0 share | top-5% share |
|--------|-----------|----------|--------------|
| **East end 270–290°E (American coast)** | **0.096–0.117** | **80.7%** | **23.3%** |
| Central-east Pacific 240–270°E | 0.049 | 63.3% | 19.4% |
| Nino3.4 center 190–240°E | weakest | — | 6.4% |

- Field-wide 65.9% of grid cells positive, mean net memory +0.056 (cross-year corr 0.050 vs random −0.005).
- Memory intensity tracks seasonal ρ amplitude (pearson +0.40): where the seasonal front is strong, position recurs predictably each year.
- **Closed loop**: east-end cold-tongue steepest location = strongest La Niña condensation (p=0.0011) = island-network same-season recurrence each year.

### 3.4 Memory Position Insensitive to ENSO

The spatial distribution of positional memory shows no systematic variation across ENSO phases (p=0.28) — **the position is anchored by non-ENSO factors**. This is the first piece of evidence pointing to "fixed topography."

### 3.5 Volcanic Activity Excluded (EPR Three-Line NULL)

| Line | Result |
|------|--------|
| Temporal coupling | all 8 synchronous r < 0.15; surrogate p=0.4–0.9 (within null) |
| High/low volcanic months | top25% vs bottom25% east-end SST no difference (p=0.78/0.87) |
| Spatial coupling | volcano density × memory field **r=−0.29 (negative, opposite direction)** |

**Conclusion: dynamic volcanic activity (USGS earthquakes/energy) is not the source of east-end condensation.**

### 3.6 Topographic Anchoring: Decisive Test (ETOPO)

| Region | net memory | elevation | roughness | gradient |
|--------|-----------|-----------|-----------|----------|
| East end 270–290°E (strongest memory) | **0.109** | **−2636 m** | **183** | **1148** |
| Nino3.4 190–240°E (weakest memory) | 0.033 | −4629 m | 115 | 467 |

- The strongest-memory seafloor is the **EPR/Galapagos shallow ridge** (~2000 m shallower, roughness/gradient 2.5× the abyssal plain).
- net memory ~ grad: r=+0.038 (p=7.5e-46, n=143k); **shuffle null z=+13.9**.
- Top-10% memory cells: gradient 1141 vs bottom-10% 941 (Δ +200).
- **Robust after excluding shelf (>1000 m deep water)**: r=+0.040 (p=5.6e-49), top-decile Δ +211.
- **Mechanism (to be verified)**: EUC hits the shallow ridge → upwelling → cold wake → front anchored at the same topographic feature each year → seasonal memory recurs.

**Direction correction:** not "volcanic eruption" but "submarine geometry"; level is "positional anchoring" (not yet "phase-transition triggering").

### 3.7 Saddle-Point Topology of the ENSO Field: Density × Activity, Not a Ring

The Ô framework's saddle-point (Morse) operator classifies each grid cell by the sign pattern of the Hessian eigenvalues (saddle = det H < 0). In tropical cyclones this operator detects a **saddle ring** — a uniform ring of saddle points in the eyewall shear band (64 cases, two hemispheres, four ocean basins; ring saddle fraction 55–63%, Rayleigh uniformity < 0.16, time-stable). Here we apply the same operator to the ENSO domain (which has no central vortex) to ask two questions: (i) does the ENSO field contain a saddle ring? (ii) where do ENSO saddle points actually concentrate?

**Q1 — No saddle ring (field-level fraction is trivial and phase-insensitive).** Monthly saddle fraction of the OISST anomaly field is ≈ 0.564 with essentially zero variance across 528 months (std ≈ 0.0015), and shows **no sensitivity to ENSO phase**: Kruskal–Wallis p = 0.137, Mann–Whitney El Niño vs La Niña p = 0.534, Spearman ρ vs ONI = −0.028 (p = 0.52), null-control z = 1.58 (not significant). Unlike the cyclone eyewall (a single center → ring), a spatially averaged fraction cannot capture ENSO structure — field-level averaging is an operator misuse for a system without a center. The 3D extension (ERA5, 8 pressure levels, 2015 El Niño / 2016 transition / 2020 La Niña) confirms: phase-insensitive global fractions (geopotential 65.8/65.7/65.6%; vorticity 78.2% flat), and the apparent Nino3.4 ring-core +6 pp advantage is dissolved by month-by-month sign tests (23/36 positive, p = 0.132).

**Q2 — Saddle points concentrate where activity density is high (co-location, not ring).** Although the field-level fraction is trivial, the *spatial distribution* of saddle points carries real signal: saddle density × activity density (ρ) is monotonically increasing across 10 bins (**binned Spearman = 1.000**) and at grid-cell level (ocean + tropics mask) Spearman = 0.42 (shuffle null z = +85). Saddle points cluster on high-noise frontal bands — the same geometric signature as the cyclone eyewall shear band, but organized as diffuse frontal lines rather than a closed ring around a center.

**Q3 — Saddle structure is insensitive to ENSO phase (fair-control audit).** The apparent phase "skeleton" (El Niño vs La Niña mean map corr 0.372) is an artifact of sample-size-inflated averaging: a fair 30v30-month comparison gives corr 0.117 vs random-month null 0.127 (z = −2.0). This is cross-consistent with the condensation-line result above (memory position p = 0.28): **the saddle structure, like the memory structure, is anchored by non-ENSO factors.**

**Q4 — Where the saddle points actually sit: the air–sea interface (wind-stress-curl test).** Bulk wind stress τ = ρ_a C_D |U| U from ERA5 10-m wind; the saddle density of curl τ × |curl τ| is grid-correlated r = 0.24–0.28 (shuffle null r ≈ 0.004, **z = +141 to +222**) — the first positional signal for ENSO saddle points: they concentrate on the wind-stress-curl coupling bands (ITCZ/cold-tongue frontal zones). But the phase skeleton of that position is even *less* similar than random (z = −13.8): the coupling bands do not move with ENSO phase. No ring, but a real interface-localized saddle structure.

**Bottom line:** ENSO has no saddle ring (2D and 3D, field and ring tests all negative) — saddle-ring is a **vortex-specific structure** — but the saddle points that do exist co-localize with the same condensation bands measured in §3.1–§3.6. The saddle operator and the condensation operator are complementary views of the same marginal frontal structure: condensation measures *folding* (D_fold ≈ 1.33), saddle structure measures *critical-point type distribution* (saddle density × activity density, grid r = 0.42, z = +85; interface r = 0.24–0.28, z = +141 to +222).

---

## 4. Discussion

### 4.1 Relation to the Three-Number Formula

The three-number formula (x, ẋ, V) captures the ENSO center state; the condensation structure captures the marginal frontal systems. The two are complementary: **the center is dominated by interannual variability (hard to predict); the margins are dominated by seasonal/topographic processes (predictable each year)**. This explains why predicting the ENSO center is so hard — the truly predictable structure does not live in the center.

### 4.2 Physical Meaning of "Strong at Both Ends, Weak in the Center"

The Nino3.4 center is the epicenter of ENSO interannual variability, with its front positions reshuffling each year (86%); the east-end cold tongue is fixed by coastal upwelling plus topography (Galapagos/EPR), recurring at the same location each year. **Memory is the spatial distribution of predictability** — it concentrates at topographically anchored margins, not at the ENSO index.

### 4.3 Limitations

1. Small effect sizes (r≈0.04) — large areas of weak memory values dilute the signal; signal concentrates in the east-end corridor.
2. Unbalanced regional sample sizes (east end 1649 cells vs much larger Nino3.4).
3. Causal direction: correlation supports anchoring but does not exclude other east-end peculiarities.
4. 0.25° averaging may smooth small-scale topography.
5. SST–DART cross-domain test: F2 spectral isomorphism NULL (surrogate p=0.986) — **the condensation structure is not isomorphic with the earthquake/tsunami domain**, the ρ field is an ENSO-domain-specific hard core.

---

## 5. Conclusion

**The structural (predictable) content of ENSO does not live in the center index but in the condensation structure of marginal frontal systems.** This structure: condenses to near-1D (D=1.33, z=−46), retains weak positional memory (1.12×, p=0.004), concentrates at both Pacific ends, is anchored by fixed submarine topography (ETOPO z=+13.9), and is not driven by volcanic activity (EPR NULL). The saddle-point (Morse) operator adds a topological cross-check: ENSO has **no saddle ring** (field fraction ≈ 0.564 phase-insensitive; 3D extension confirms; fair-control z = −2.0) — confirming that saddle rings are vortex-specific — but its saddle points co-localize with the same marginal frontal condensation bands (density × activity grid r = 0.42, z = +85; air–sea interface |curl τ| r = 0.24–0.28, z = +141 to +222).

One answer to "why is ENSO so hard to predict": **we keep measuring the center, but the structure lives at the margins.**

---

## 6. Evidence-Level Checklist

| Claim | Evidence level | Note |
|-------|---------------|------|
| ρ-field condensation is the norm (D=1.334, z=−46.4, 60/60 months) | ✅ Cross-validated | 60 independent months + April pilot replication + dual resolution |
| 0.25° condenses more strongly than 1° (z −46 vs −11) | ✅ Cross-validated | Dual-resolution independent measurements |
| Temporal memory weak but significant (1.12×, p=0.004) | ✅ Cross-validated | n=55 vs 48 Welch t; both resolutions agree |
| 86% position reshuffle (bulk not moving) | ✅ Cross-validated | Both resolutions agree (1° 92.6%, 0.25° 86%) |
| Seasonal memory "strong at both ends, weak in center" | ✅ Cross-validated | Regional statistics + top-5% distribution independent |
| East end 270–290°E strongest memory (0.181 vs 0.142, p=3.4e-6) | ✅ Cross-validated | Seasonal memory localization report |
| Memory position insensitive to ENSO (p=0.28) | ✅ Cross-validated | Phase-stratified no systematic difference |
| Volcanic activity excluded (EPR three-line NULL) | ✅ Cross-validated | Temporal/contrast/spatial three independent lines |
| Topographic anchoring supported (ETOPO z=+13.9, robust after shelf exclusion) | ✅ Cross-validated | Shuffle null + deep-water mask robustness |
| Memory tracks seasonal ρ amplitude (r=0.40) | ⚠️ Single source | Single correlation measure, needs independent verification |
| EUC upwelling → front anchoring mechanism | 🔄 To be verified | Correlation supports; mechanism not directly measured |
| La Niña condenses most strongly | ✅ Cross-validated | p=0.0011, phase analysis |
| ENSO has no saddle ring (field fraction ≈0.564, phase-insensitive) | ✅ Cross-validated | 528-month KW p=0.137, MWU p=0.534, ρ=−0.028, null z=1.58; 3D extension (ERA5 8 levels) agrees |
| Saddle points co-localize with activity density (grid r=0.42, binned Spearman=1.000) | ✅ Cross-validated | 10-bin monotonicity + grid-level null z=+85 |
| Saddle structure insensitive to ENSO phase (fair-control z=−2.0) | ✅ Cross-validated | 30v30-month corr 0.117 vs random null 0.127 |
| Saddle points concentrate at the air–sea interface (× |curl τ| r=0.24–0.28, z=+141 to +222) | ✅ Cross-validated | Bulk wind stress from ERA5 10-m wind; shuffle null r≈0.004 |

---

## 7. References and Resources

1. NOAA OISST v2.1: https://www.ncei.noaa.gov/products/optimum-interpolation-sst
2. ERA5: https://cds.climate.copernicus.eu
3. HadISST: https://www.metoffice.gov.uk/hadobs/hadisst/
4. NOAA ETOPO1: https://www.ngdc.noaa.gov/mgg/global/
5. USGS Earthquake Catalog: https://earthquake.usgs.gov/fdsnws/event/1/
6. ENSO Simplified Prediction Note (predecessor): https://doi.org/10.5281/zenodo.21626908
7. ENSO Watch Dashboard: https://mmdr10.github.io/enso-watch
8. Eyewall Saddle Ring paper (saddle-ring operator origin, TC domain): https://doi.org/10.5281/zenodo.21876754
9. Saddle-ring cross-domain report (ENSO 2D/3D + Greenland, Phase 30–35): `projects/atmosphere/vorticity_dfold/research_notes.md` §Phase 30–35

**Data and code:** test files (scripts + data) in `projects/enso/tests/` (updated 2026-08-04). All open source, CC BY 4.0.

---

*This paper was independently written by AI research assistant tygtDc (Deep Research), based on ENSO condensation-structure measurements using the Ô geometric framework. All data, code, and methodology are fully open source.*
