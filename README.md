# Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot

Reproducible R workflow for the spatially explicit, individual-based population
simulations (RangeShiftR v3.0.0) and figures found in our study named above. The model
projects range dynamics of *Calophyllum paniculatum* across the CAZ+ corridor,
Madagascar, under combined climate change, land cover change (deforestation), defaunation (loss of lemur-mediated
seed dispersal) and fungal-wilt pathogen pressure.

## Authors


Emma L Underwood^1,*^, Kerry A Brown^1^, Rebekka Allgayer^2^, Mark Mulligan^3^, Nigel Walford^1^, Jette Wolff^4^.


1. Department of Geography, Geology, and the Environment, Centre for Engineering, Environment and Society Research (CEESR), Kingston University, Penrhyn Rd, Kingston upon Thames, KT1 2EE, UK

2. School of Biological Sciences, University of Aberdeen, King's College, Aberdeen, AB24 3FX, UK

3. Department of Geography, King’s College London, 30 Aldwych, London, WC2B 4BG, UK 

4. Institute of Biochemistry and Biology, University of Potsdam, Maulbeerallee 3, 14469, Potsdam, Germany


*Corresponding author: Emma L Underwood
Email: [elmah2707@gmail.com](mailto:elmah2707@gmail.com)


Code was implemented by **[ELU]** with contributions from **[JW]** and **[RA]**.


## Funding


- Royal Geographical Society (with IBG), Researcher Award Grant (including Albert Reckitt Award and Dudley Stamp 
Memorial Award), PRA47.22, Emma L Underwood.
- Kingston University London, Research Degree Studentship in the Faculty of Science, Engineering and Computing, School of 
Engineering and Environment, Emma L Underwood (formerly Hall).

## Related publication

This repository contains the R scripts needed to reproduce the analyses and
figures in manuscript *Climate change, defaunation and disease interact to accelerate biodiversity loss in a tropical hotspot* by Underwood, E.L., Brown, K.A., Allgayer, R., Mulligan, M., Walford, N., Wolff, J., submitted in 2026.
DOI: * [TBC pending pre-print/peer review]

## Abstract

**Context**
Tropical endemic forests face multiple stressors of climate change, deforestation, defaunation and emerging pathogens and diseases, yet most biodiversity forecasts model these pressures in isolation or with non-spatially explicit methods, limiting their ability to support comprehensive conservation efforts.

**Objectives**
We quantify how multiple global-change stressors acting through habitat change, dispersal limitation and emerging disease jointly shape population persistence and test whether their combined effects generate greater declines than would be expected from individual drivers acting alone.

**Methods**
Using the individual-based modelling platform RangeShifter, we simulated population and occupancy dynamics across 30 scenarios combining future changes to habitat suitability based on combined climate and deforestation predictions, three levels of lemur defaunation, and spatially explicit wilt pathogen spread.

**Results**
Climate change and deforestation alone drove occupancy declines exceeding 90% under high-emissions scenarios, with rapid contractions tracking landscape step-changes in habitat suitability. Lemur defaunation independently reduced starting occupancy by 78–93% relative to baseline, irrespective of climate trajectory, and accelerated spatial fragmentation of reproductive populations. Where wilt pathogen dynamics were incorporated, local extirpation occurred within 13–34 years of introduction across all scenario combinations, with populations subject to prior climate and defaunation stress collapsing up to seven years faster.

**Conclusion**
Our results support a shift from single-threat forecasting to multi-threat ecological assessment, in which persistence is determined by the outcome of interactive abiotic and biotic processes across space and time. Long-term persistence in tropical biodiversity hotspots will depend on retaining climatically suitable refugia, as well as maintaining species interactions and demographic processes that allow populations to function within them, emphasising the need for multi-threat frameworks in future conservation risk assessments.


**Keywords**: climate change, tropics, endemism, biodiversity hotspot, mechanistic model, scenario


## Workflow requirements

The workflow was developed and run in R. To reproduce the computational
environment:

**1. Install RangeShiftR v3.0.0** (this is currently not possible on CRAN; use pak package):

```r
# install.packages("pak")
pak::pak("RangeShifter/RangeShiftR-pkg/RangeShiftR@v3.0.0")
```

**2. Open the project via `RProj_Underwood-global-change-drivers-accelerate-biodiversity-loss-2026.Rproj`** so that all paths resolve through
`here::here()` from the repository root folder.

**Programming environment used**

- Platform: x86_64-pc-linux-gnu
- R version: *R version 4.5.2 (2025-10-31)*
- Attached R packages: RangeShiftR (3.0.0), terra (1.9-11), data.table (1.18.2.1), dplyr (1.2.0), tidyr (1.3.2), ggplot2 (4.0.2),
  patchwork (1.3.2), cowplot (1.2.0), scales (1.4.0), viridis (0.6.5), foreach (1.5.2), doParallel (1.0.17).


**Additional (non-R) software used**

- **QGIS** *[3.34.6]* — final spatial figures (manuscript Figure 1, Figure 4,
  Supplementary Figure S4) were composed in QGIS from the rasters and point
  layers written by `03__occupancy_probability.R`.
- **Canva** — final layout and annotation of several figures (see *Details*).

**Repository version / release**

The exact repository state accompanying the manuscript is tagged in Git as a release.

## Data availability

**Raw input data (not in this repository)**

- Ensemble species-distribution-model rasters (habitat suitability for tree,
  and pathogen suitability), current and future SSP projections:
  **Zenodo DOI [10.5281/zenodo.21160109](https://doi.org/10.5281/zenodo.21160109)**.
  Download and place under `data/SDM_ensemble/` before running
  `00__data_preparation.R` (this path can be changed in `config.R`).

**Additional input data (small, included in this repository under `data/`)**

- `data/StudyArea/CAZ_AOI.shp` — CAZ corridor area of interest (protected areas
  joined via convex hull).
- `data/StudyArea/calo_current_CAZ_100m.tif`, `CAZ_100m.tif`, `MDG_template_1km.tif`
  — study-area templates (extent, resolution, CRS ESRI:102022) used to align all
  layers and to place model output in real-world coordinates.
- `data/StudyArea/WDPA_clipCAZ.shp` — protected-area boundaries (spatial context) clipped to study exent.
    https://www.protectedplanet.net/. ©ProtectedPlanet 2014-2026. All rights reserved
- `data/StudyArea/calo_ecoregsA_intersect.shp` — ecoregion polygons (context for map background).
    https://rcmrd.africageoportal.com/datasets/rcmrd::africa-ecoregions

**Generated data (produced by the scripts; not tracked in Git)**

- `outputs/simulations/…` — raw RangeShiftR `_Pop.txt` / `_Range.txt` files.
- `outputs/occupancy/OccRasters/`, `outputs/occupancy/OccPoints/` — occupancy
  probability GeoTiffs and GeoPackages.
- `outputs/analysis/…`, `outputs/figures/…` — summary tables and figure files.
  *DOI for outputs TBC

## Repository structure

```bash
├── README.md
├── LICENSE
├── RProj_Underwood-global-change-drivers-2026.Rproj # RStudio project (sets working directory)
├── .gitignore                    # excludes large generated outputs
├── config.R                      # project root (here::here) + package loading
├── R/                            # all reusable functions, grouped by script / stage
│   ├── functions_utils.R         #   generic helpers (load, cleanup, timing)
│   ├── functions_data_prep.R     #   landscape / pathogen layer building
│   ├── functions_simulation.R    #   simulation setup + memory management for HPC
│   ├── functions_results.R       #   aggregation, refugial + range-shift metrics
│   ├── functions_occupancy.R     #   occupancy rasters + points (stage 03)
│   └── functions_figures.R       #   publication figures
├── 00__data_preparation.R        # crop/align SDMs -> RangeShiftR input matrices
├── 01__simulations.R             # run the 30 scenario simulations
├── 01b__simulation_validation.R  # summarise runs -> simulation_validation.csv
├── 02__process_results.R         # range aggregation, metrics, Supp. Figs
├── 03__occupancy_probability.R   # Pop files -> occupancy rasters + points
├── 04__publication_figures.R     # main manuscript figures (Fig 3/4/5)
├── data/                         # small inputs
│   ├── StudyArea/                #   templates, AOI, PA + ecoregion layers
│   └── SDM_ensemble/             #   (populate from Zenodo before running 00)
└── outputs/                      # generated; git-ignored except final figures
    ├── simulations/
    ├── occupancy/{OccRasters,OccPoints}/
    ├── analysis/
    └── figures/
```

## Simulation details

**Experimental design (30 simulations).** Three climate scenarios
(`current`, SSP1-2.6 = "Low", SSP5-8.5 = "High"), three defaunation levels set
via emigration probability (0.20 = Low, 0.10 = Medium, 0.05 = High), and
pathogen-introduction timing across four batches:

| Batch | Climate | Defaunation | Pathogen introduced | n |
|-------|---------|-------------|---------------------|---|
| 1 | Current / SSP1 / SSP5 | Low / Medium / High | none (baseline) | 9 |
| 2 | Current / SSP1 / SSP5 | Low / Medium / High | year 100 (current) / 101 (SSP) | 9 |
| 3 | SSP1 / SSP5 | Low / Medium / High | year 131 (31-yr pre-conditioning) | 6 |
| 4 | SSP1 / SSP5 | Low / Medium / High | year 161 (61-yr pre-conditioning) | 6 |

SSP3 was evaluated during development and excluded from our final study.
The scenario-to-`sim_id` map is written to
`outputs/simulations/Master_simulation_mapping.csv` and links every output file
to its scenario.

**Workflow order:**

```r
source("00__data_preparation.R")      # build model inputs (needs Zenodo data first)
source("01__simulations.R")           # run simulations (requires HPC access)
source("01b__simulation_validation.R")# summarise / validate runs
source("02__process_results.R")       # aggregate results + Supp. Figs
source("03__occupancy_probability.R") # occupancy rasters + points for mapping
source("04__publication_figures.R")   # main manuscript figures
```

`01__simulations.R` is memory-intensive; it is designed to run in sim-index
chunks (e.g. `1:9`, `10:18`, `19:30`) as separate HPC jobs. 
Set any cluster-specific library paths for your HPC.

**Figures (transparency for review).** Final manuscript figures were
assembled from a mix of R, QGIS and Canva; R outputs are the analytical
originals and differ only in fonts, colours, annotation, and cropping.

| Manuscript figure | Produced with | Source in this repository |
|-------------------|---------------|---------------------------|
| Figure 1 | QGIS only | *(no R code)* — study-area / SDM map |
| Figure 2 | Canva only | *(no R code)* — conceptual diagram |
| Figure 3 | R → Canva | `04__publication_figures.R` (occupancy trajectories + decline-rate heatmap) |
| Figure 4 | R rasters → QGIS | `03__occupancy_probability.R` (`OccRasters` / `OccPoints`) |
| Figure 5 | R → Canva (Medium climate removed) | `04__publication_figures.R` (years-to-extinction heatmap, computed from Range files) |
| Figure S1 | Canva only | *(no R code)* |
| Figure S2 | R | baseline population dynamics ± SD, Sim 1 (`02__process_results.R`) |
| Figure S3 | R | stage composition, Sims 1,2,4,5,6,8,9,10,12 (`02__process_results.R`) |
| Figure S4 | R rasters → QGIS | starting occupancy, Sims 1,5,9 (`OccRasters`) |
| Figure S5 | R | baseline population dynamics ± SD, Sim 1 with High Climate & Deforestation scenarios (`02__process_results.R`) |

**Coordinates** RangeShiftR uses a false origin (0,0); `x`/`y` in the
Pop files are 0-based cell indices. `03__occupancy_probability.R` converts these
to projected cell centres (ESRI:102022, Africa Albers Equal Area Conic) against
the 100 m study area template before writing rasters.

## License

> This project is released under the *GNU General Public License v3.0* license. See `LICENSE`.
