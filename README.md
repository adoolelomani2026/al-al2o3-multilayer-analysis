# Al/Al2O3 Multilayer Analysis

Python analysis scripts and selected output figures for Al/Al2O3 multilayer characterization data from STEM-EELS, UV-Vis, cathodoluminescence, and photoluminescence measurements.

Repository target:

```powershell
https://github.com/adoolelomani2026/al-al2o3-multilayer-analysis
```

## Folder Structure

```text
.
├── data/
│   └── raw/
│       ├── cl/
│       ├── eels/
│       ├── pl/
│       └── uvvis/
├── docs/
│   ├── extracted_text/
│   └── report/
├── outputs/
│   ├── cl/
│   ├── eels/
│   ├── pl/
│   └── uvvis/
└── scripts/
```

## Outputs

The scripts are configured to regenerate only these selected PNG figures:

- `outputs/cl/cl_spectrum_ev_normalized.png`
- `outputs/eels/EELS_summary_combined.png`
- `outputs/pl/pl_emission_publication.png`
- `outputs/uvvis/uvvis_reflectance_corrected_publication_v2.png`

## Setup

Create and activate a virtual environment:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## Run

Run all analyses:

```powershell
python scripts\run_all.py
```

Or run one analysis at a time:

```powershell
python scripts\analyze_stem_eels.py
python scripts\analyze_uvvis.py
python scripts\analyze_cl.py
python scripts\analyze_pl.py
```

## Upload To GitHub

From this folder:

```powershell
git init
git branch -M main
git add .
git commit -m "Initial analysis repository"
git remote add origin https://github.com/adoolelomani2026/al-al2o3-multilayer-analysis.git
git push -u origin main
```
