# Al/Al₂O₃ Multilayer Analysis

Reproducible analysis workflows and selected processed outputs supporting
*Ultraviolet Plasmonic Response and Model-Predicted Type-II Hyperbolicity in
Al/Al₂O₃ Multilayers*.

This repository contains Python workflows for:

- STEM-EELS dielectric reconstruction
- UV–Vis reflectance
- Cathodoluminescence spectroscopy
- Photoluminescence emission
- VASE-workbook constituent optical functions and local effective-medium-theory
  (EMT) tensor estimates

## Reproducibility scope

The repository includes the processed numerical data underlying the manuscript's
EELS dielectric-response figure and EMT tensor figure. The EMT calculation uses
the fitted Al B-spline optical function and fixed handbook Al₂O₃ optical
constants retained in `Aneesa.xlsx`.

The primary EMT result uses the nominal design filling fraction
`f_Al = 5 / (5 + 7) = 0.417`. A secondary `f_Al = 0.497` case is retained only
as a geometry-sensitivity calculation. Because the exact workbook-to-model
export lineage could not be uniquely established, the secondary case is not
described as a workbook-matched fitted geometry. Both cases preserve the
model-supported type-II sign condition over 245–600 nm.

The VASE fit was isotropic; it did not directly retrieve the uniaxial tensor
components.

## Repository structure

```text
data/
  processed/
    emt/      Wavelength-resolved EMT components and sign-test summary
  raw/
    cl/       Cathodoluminescence raw spectrum
    eels/     STEM-EELS dielectric and energy-loss CSV files
    ellipsometry/  Optical-function workbook used by the EMT script
    pl/       Photoluminescence emission scans
    uvvis/    UV–Vis raw scan exports

outputs/
  cl/       Selected CL figure
  eels/     Manuscript-aligned EELS figures
  emt/      EMT tensor figure in PNG and PDF formats
  pl/       Selected PL figure
  uvvis/    Selected UV–Vis figure

scripts/
  analyze_cl.py
  analyze_emt.py
  analyze_pl.py
  analyze_stem_eels.py
  analyze_uvvis.py
  run_all.py
```

## Run the workflows

```bash
python -m pip install -r requirements.txt
python scripts/run_all.py
```

Individual workflows can also be run directly. Generated figures are written
to `outputs/`; the EMT script additionally refreshes the processed CSV files in
`data/processed/emt/`.

## Retained PL metadata

The archived emission record identifies 565 nm xenon-lamp excitation,
300–800 nm acquisition in 1 nm steps, nominal 5 nm excitation/emission
bandpasses, a 0.10 s dwell per point, one repeat, and a visible PMT-900
detector. Excitation irradiance and illuminated area were not retained, so no
absolute quantum yield or calibrated radiance is claimed.
