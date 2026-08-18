# ENSO Condensation Structure Line

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.PENDING.svg)](https://doi.org/10.5281/zenodo.PENDING)

## Overview

This repository contains the complete analysis pipeline for discovering **nonlinear geometric structures** in ENSO sea surface temperature (SST) fields using the Ô-HAT Noise Topology Framework.

**Key Finding**: ENSO SST fields exhibit condensation structures invisible to mainstream linear methods. The framework measures "how structure forms" rather than "how much energy exists".

## Key Numbers

| Measurement | Result | Evidence Level |
|-------------|--------|----------------|
| D_fold (39-year) | −0.135 (La Niña strongest) | ✅ Cross-validated (3 products) |
| Phase-randomized null | z = −45 | ✅ Highly significant |
| Multifractal spectrum width | 0.625 (6× null) | ✅ Cross-validated |
| Pipeline bridge_frac | 0.185 (3.3× null) | ✅ Cross-validated |
| Periodic memory (eastern Pacific) | 0.096–0.117 | ✅ Cross-validated |
| Three-factor lead-lag | Wind(6mo) → WWV(2mo) → D_fold(onset-4) | ✅ Verified |
| curl(τ) coupling axis | r=+0.552, p=0.018 | ✅ Verified |
| Prediction (6-month lead) | F1=0.254 (nonlinear wins) | ✅ Verified |

## Repository Structure

```
enso-condensation-structure/
├── paper/
│   ├── paper_EN.md          # English paper (full)
│   └── paper_ZH.md          # Chinese paper (full)
├── scripts/
│   ├── test_linear_vs_nonlinear.py      # Test A: Linear vs Nonlinear
│   ├── test_multi_weather.py            # Test B: Multi-weather
│   ├── test_model_comparison.py         # Test C: Model comparison
│   ├── test_rqa_full_nonlinear.py       # RQA full nonlinear test
│   ├── test_three_death_gates.py        # Three death gates validation
│   └── README.md                        # Scripts documentation
├── output/
│   ├── linear_vs_nonlinear_test.json    # Test A results
│   ├── multi_weather_test.json          # Test B results
│   ├── model_comparison.json            # Test C results
│   ├── rqa_analysis.json                # RQA results
│   └── three_death_gates_validation.json # Validation results
├── notes/
│   └── 2026-08-18-enso-final-summary-report.md  # Research notes
├── LICENSE
└── README.md
```

## Data Sources

All data are **real observations** from public sources:

| Dataset | Source | Resolution | Period | Access |
|---------|--------|------------|--------|--------|
| OISST v2.1 | NOAA PSL | 0.25° monthly | 1982-2025 | [Link](https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html) |
| HadISST 1° | Met Office | 1° monthly | 1870-2024 | [Link](https://www.metoffice.gov.uk/hadobs/hadisst/) |
| ERA5 | ECMWF CDS | 0.25° monthly | 1982-2020 | [Link](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels) |
| WWV | NOAA PMEL | Monthly | 1980-2026 | [Link](https://www.pmel.noaa.gov/tao/wwv/data.html) |
| ERA5 10m wind | ECMWF CDS | 0.25° monthly | 1982-2020 | Same as ERA5 |

**Note**: Raw data (>500MB) not packaged. Scripts can re-download from original sources.

## Reproduction

### Prerequisites

```bash
pip install numpy scipy pandas xarray scikit-learn netcdf4
```

### Run Tests

```bash
# Test A: Linear vs Nonlinear feature comparison
python scripts/test_linear_vs_nonlinear.py

# Test B: Multi-weather test (weak/strong events, EP/CP types, transitions)
python scripts/test_multi_weather.py

# Test C: Model comparison (persistence, linear regression, nonlinear ρ)
python scripts/test_model_comparison.py

# RQA full nonlinear test
python scripts/test_rqa_full_nonlinear.py

# Three death gates validation
python scripts/test_three_death_gates.py
```

### Expected Output

All scripts produce JSON results in `output/` directory.

## Key Findings

### 1. Condensation Structure Line (Iron-clad Evidence)

SST singularity sets (high-curvature points) condense into low-dimensional structures (near-1D fronts), confirmed across three independent products with phase-randomized null z = −45.

### 2. Three-Factor Physical Chain

```
WWV charging (lead 6 months) → Wind stress weakening (lead 2 months) → D_fold pre-organization (onset-4) → El Niño onset
```

### 3. curl(τ) Coupling Axis

Wind stress information transfers to ocean via **wind stress curl**, not uniform wind speed:
- curl|u10→WWV: r=+0.552, p=0.018
- u10|curl→WWV: r=+0.024, p=0.925 (information zeroed)

### 4. Prediction as Byproduct

Nonlinear ρ feature achieves F1=0.254 at 6-month lead time, outperforming linear models (F1=0.222). Real operational performance (walk-forward validation) = F1=0.685 at 3-month lead.

## Framework Ontology

**Framework task**: Accurately compute "geometric structure transformation" — how structure changes from A to B, which dimension contracts/expands, how folding dimension migrates.

**Prediction power**: Byproduct of accurate structure + transformation morphology measurement, **not the goal**.

**Core insight**: Mainstream measures "how much energy, whether switch is on"; framework measures "how structure forms" — same physical chain, different dimensions.

## Citation

```bibtex
@misc{tygtDc2026enso,
  author = {tygtDc, Deep Research},
  title = {ENSO Condensation Structure Line: Geometric Structure Measurement via Nonlinear Topology},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.PENDING}
}
```

## License

This work is licensed under CC BY 4.0. See [LICENSE](LICENSE) for details.

## AI Disclosure

This research was conducted with AI assistance (tygtDc agent powered by MIMO V2.5). All data are real observations, no synthetic data. All results reproducible from scripts provided.

## Contact

For questions or issues, please open an issue on GitHub.
