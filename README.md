# Snapshot Brasil — Camera-Trap Analysis Pipeline

Teaching modules, data, and analysis code for the Snapshot Brasil camera-trap
network, covering a merged dataset spanning three camera-trap networks
(Snapshot Brasil, the Atlantic Forest camera-trap compilation, and public
Wildlife Insights projects).

**Site:** once GitHub Pages is enabled for this repo (Settings → Pages →
Source: `main` branch, `/ (root)`), the landing page is served from
`index.html` at the repo's Pages URL. English/Portuguese are both on the same
page via a language toggle.

## Contents

- `index.html` — bilingual landing page linking every module and dataset.
- `modules/module1/` — Foundation & Detection Histories (`en/`, `pt/`: HTML
  report + annotated R script each).
- `modules/module2/` — Covariates, Sampling & Protected Areas (`en/`, `pt/`:
  HTML report + annotated R script each).
- `data/` — final combined dataset (camera deployments, mammal detections,
  remote-sensing covariates) as a zip, plus a data dictionary README.

Modules 3–7 (single-species occupancy, community/joint models, prediction &
mapping, species interactions, temporal analysis) are still under revision
and will be added once finalized.

## Running a module script

Each module ships as a rendered HTML report (open directly, no setup needed)
and a standalone, heavily-annotated R script (`Rscript module<N>_<name>.R`
from inside that module's folder, with `data.table` and `ggplot2` installed).
