#!/usr/bin/env python3
"""
Zenodo upload script for ENSO condensation structure line
"""
import requests
import os
import json
from pathlib import Path

ZENODO_TOKEN = "ukDixXoEjwQj28lIMKXcwUZmf3Wbl4oh8kUbu24eJ0pN13AJKIPwhudbMvDD"
GITHUB_REPO = "https://github.com/MMDR10/enso-condensation-structure"

# Headers
headers = {
    "Authorization": f"Bearer {ZENODO_TOKEN}",
    "Content-Type": "application/json"
}

# Create new deposition
print("Creating Zenodo deposition...")
resp = requests.post(
    "https://zenodo.org/api/deposit/depositions",
    headers=headers,
    json={}
)
deposition = resp.json()
deposition_id = deposition['id']
bucket_url = deposition['links']['bucket']
print(f"Deposition ID: {deposition_id}")
print(f"Bucket URL: {bucket_url}")

# Collect files with flat naming
repo_dir = Path("/app/working/workspaces/tygtDc/projects/enso/release-repo")
files_to_upload = []

# Paper files
for f in repo_dir.glob("paper/*"):
    if f.is_file():
        flat_name = f"paper_{f.name}"
        files_to_upload.append((f, flat_name))

# Scripts
for f in repo_dir.glob("scripts/*.py"):
    if f.is_file():
        flat_name = f"scripts_{f.name}"
        files_to_upload.append((f, flat_name))

# Output
for f in repo_dir.glob("output/*.json"):
    if f.is_file():
        flat_name = f"output_{f.name}"
        files_to_upload.append((f, flat_name))

# Notes
for f in repo_dir.glob("notes/*.md"):
    if f.is_file():
        flat_name = f"notes_{f.name}"
        files_to_upload.append((f, flat_name))

# Root files
for f in ["README.md", "LICENSE", ".gitignore"]:
    fp = repo_dir / f
    if fp.exists():
        files_to_upload.append((fp, f))

print(f"\nFound {len(files_to_upload)} files to upload")

# Upload files
print("\nUploading files...")
for filepath, flat_name in files_to_upload:
    print(f"  Uploading {flat_name} ({filepath.stat().st_size} bytes)...")
    with open(filepath, 'rb') as f:
        resp = requests.put(
            f"{bucket_url}/{flat_name}",
            headers={
                "Authorization": f"Bearer {ZENODO_TOKEN}",
                "Content-Type": "application/octet-stream"
            },
            data=f
        )
        if resp.status_code not in [200, 201]:
            print(f"    ERROR: {resp.status_code} - {resp.text}")
        else:
            print(f"    OK")

# Update metadata
print("\nUpdating metadata...")
metadata = {
    "metadata": {
        "title": "ENSO Condensation Structure Line: Geometric Structure Measurement via Nonlinear Topology",
        "upload_type": "publication",
        "publication_type": "article",
        "description": """<p>Discovery of <strong>nonlinear geometric structures</strong> in ENSO sea surface temperature (SST) fields invisible to mainstream linear methods.</p>

<h3>Key Findings</h3>
<ul>
<li><strong>Condensation structure line</strong>: SST singularity sets condense into low-dimensional structures (near-1D fronts), confirmed across three independent products (OISST/HadISST/ERA5) with phase-randomized null z = −45</li>
<li><strong>Three-factor physical chain</strong>: Wind stress weakening (lead 6 months) → Warm Water Volume charging (lead 2 months) → D_fold condensation pre-organization (onset-4) → El Niño onset</li>
<li><strong>curl(τ) coupling axis</strong>: Wind stress information transfers via wind stress curl, not uniform wind speed (mediation test: curl|u10→WWV r=+0.552, p=0.018)</li>
<li><strong>Periodic memory</strong>: Strongest at eastern Pacific (270-290°E), suggesting topographic anchoring</li>
<li><strong>Island network topology</strong>: Condensation cores connect through pipelines into fully-connected island clusters (K₉-K₁₁), with 53% of cores entering the network</li>
</ul>

<h3>Prediction as Byproduct</h3>
<p>Nonlinear ρ feature achieves F1=0.254 at 6-month lead time, outperforming linear models (F1=0.222). Real operational performance (walk-forward validation) = F1=0.685 at 3-month lead.</p>

<h3>Core Insight</h3>
<p>Mainstream measures "how much energy, whether switch is on"; framework measures "how structure forms" — same physical chain, different dimensions.</p>

<h3>Data Sources</h3>
<ul>
<li>OISST v2.1 (NOAA PSL): 0.25° monthly 1982-2025</li>
<li>HadISST 1° (Met Office): 1° monthly 1870-2024</li>
<li>ERA5 (ECMWF CDS): 0.25° monthly 1982-2020</li>
<li>WWV (NOAA PMEL): Monthly 1980-2026</li>
</ul>

<h3>AI Disclosure</h3>
<p>This research was conducted with AI assistance (tygtDc agent powered by MIMO V2.5). All data are real observations, no synthetic data. All results reproducible from scripts provided.</p>""",
        "creators": [
            {
                "name": "tygtDc, Deep Research",
                "affiliation": "Independent Researcher"
            }
        ],
        "license": "CC-BY-4.0",
        "keywords": [
            "ENSO",
            "El Niño",
            "nonlinear dynamics",
            "topology",
            "condensation structure",
            "geometric structure",
            "climate",
            "SST",
            "sea surface temperature",
            "Ô-HAT framework"
        ],
        "related_identifiers": [
            {
                "identifier": GITHUB_REPO,
                "relation": "isSupplementTo",
                "scheme": "url"
            }
        ],
        "communities": [
            {"identifier": "ohat-framework"}
        ]
    }
}

resp = requests.put(
    f"https://zenodo.org/api/deposit/depositions/{deposition_id}",
    headers=headers,
    json=metadata
)
print(f"Metadata update status: {resp.status_code}")
if resp.status_code != 200:
    print(f"Error: {resp.text}")

# Publish
print("\nPublishing...")
resp = requests.post(
    f"https://zenodo.org/api/deposit/depositions/{deposition_id}/actions/publish",
    headers=headers
)
print(f"Publish status: {resp.status_code}")
if resp.status_code == 202:
    result = resp.json()
    doi = result.get('doi', 'N/A')
    doi_url = result.get('doi_url', 'N/A')
    print(f"\n✅ Published successfully!")
    print(f"DOI: {doi}")
    print(f"DOI URL: {doi_url}")
    print(f"Record URL: https://zenodo.org/record/{deposition_id}")
else:
    print(f"Error: {resp.status_code}")
    print(f"Response: {resp.text}")
