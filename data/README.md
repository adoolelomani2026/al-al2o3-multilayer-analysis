# Data

Raw inputs and processed numerical data used by the analysis scripts.

- `raw/eels/`: dielectric function and energy-loss CSV files.
- `raw/ellipsometry/`: optical-function workbook used by the EMT calculation.
- `raw/uvvis/`: Lambda 1050 UV-Vis scan exports.
- `raw/cl/`: cathodoluminescence spectrum export.
- `raw/pl/`: photoluminescence emission scan exports.
- `processed/emt/`: wavelength-resolved complex EMT tensor components, real-part
  sign products, and the geometry/classification summary.

The nominal `f_Al = 0.417` calculation is the primary EMT analysis. The
`f_Al = 0.497` result is a secondary geometry-sensitivity case and is not
presented as a uniquely workbook-matched fitted geometry.
