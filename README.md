# Al/Al2O3 Multilayer Analysis

Analysis scripts and selected output figures for Al/Al2O3 multilayer thin-film characterization data.

This repository contains Python scripts used to process and visualize:

- STEM-EELS dielectric reconstruction data
- UV-Vis reflectance data
- Cathodoluminescence spectrum data
- Photoluminescence emission data

## Repository Structure

```text
data/
  raw/
    cl/       Cathodoluminescence raw spectrum file
    eels/     STEM-EELS dielectric and energy-loss CSV files
    pl/       Photoluminescence emission scan CSV files
    uvvis/    UV-Vis raw scan exports

outputs/
  cl/       Selected CL figure
  eels/     Selected EELS figure
  pl/       Selected PL figure
  uvvis/    Selected UV-Vis figure

scripts/
  analyze_cl.py
  analyze_pl.py
  analyze_stem_eels.py
  analyze_uvvis.py
  run_all.py
