# casper_rt_paper

Repository for paper estimating Rt and rt from CASPER wastewater metagenomic
sequencing data (compared against clinical and targeted wastewater
surveillance data).

<!-- add citation / preprint link once available -->

## Repository structure

```
code/                       Analysis scripts to be run in numbered order (see pipeline below)
  functions/                Shared helper functions sourced by the analysis scripts
data/
  original/                 CASPER exports and other source data (read-only inputs)
  clean/                    Processed data produced by code/01-03 (git-ignored -- see Data & outputs below)
documents/                  Supporting reference files (run settings, digitized shedding data)
output/                     All generated results: figures, time-varying Rt/rt estimates, concordance tables 
                            (git-ignored -- see Data & outputs below)
  figures/                  All main and supplementary figures
  input_time_series/        All time series plots of processed data used for analysis
  output_concordance/       
    Rt_ERN/                 All concordance plots comparing derived Rt estimates and a summary table
    rt_growth_rate/         All concordance plots comparing derived rt estimates and a summary table
  output_time_series/       
    Rt_ERN/                 All input time series plots alongside Rt estimates
    rt_growth_rate/         All input time series plots alongside rt estimates
  time_varying_estimates/   All derived Rt/rt estimates
renv/, renv.lock            Package environment (see Setup below)
```

## Setup

This project uses [`renv`](https://rstudio.github.io/renv/) to pin R package
versions.

1. **Install/run the matching R version.** `renv.lock` records the R version the
   project was developed with (`R.Version` field at the top of the file --
   currently R 4.6.0). **renv does not install or switch R**; it
   only manages packages.
2. Open `casper_rt_paper.Rproj` in RStudio (or otherwise set the working
   directory to the repo root).
3. Restore the package library (also done in `code/00_config.R`):
   ```r
   renv::restore()
   ```
4. Every script sources `code/00_config.R` first, which loads packages,
   shared helper functions (`code/functions/`), and global plotting/analysis
   settings.

### Reproducibility notes

- `renv.lock` guarantees package *versions* match, but not the R
  *interpreter*. Ensure the R version matches the R version the project was 
  developed with (see Setup above).
- Processed data under `data/clean/` are git-ignored (see below) and are 
  reproduced by code/01-03.
- Large derived outputs under `output/` are git-ignored (see below) and are
  reproduced by code/03-06.
- Ensure the repository follows the folder structure outlined above so processed data and all outputs may be reproduced.

## Data & outputs

- `data/original/` contains read-only inputs required to run the pipeline from scratch. 
  - `casper_export_20260604/`: Untargeted wastewater monitoring data from Justen et al. 2026 (https://doi.org/10.64898/2026.07.14.26356442). Refer to the corresponding data descriptor for data and methodological details (https://doi.org/10.64898/2026.03.05.26345726). 
  - `all_epic.xlsx`: Clinical surveillance data obtained from Epic Cosmos. Disease diagnoses (COVID-19, norovirus, rotavirus, RSV) and lab results (influenza A and B) at the state and county level. git-ignored as Epic Cosmos data cannot be made publicly available and can only be accessed through approved users at participating organizations. 
  - `WWSCAN_AllPlants_SDR_2025.csv`: Targeted wastewater monitoring data from Boehm et al. 2026 (https://doi.org/10.1016/j.dib.2026.112503). Refer to the same publication for data and methodological details.
- `data/clean/` is git-ignored. It contains processed inputs produced/used throughout the pipeline.
- `output/` is git-ignored. It holds everything the pipeline generates including 
  large, computationally expensive results such as the time-varying Rt/rt 
  estimates. Look for commented out `saveRDS()` lines and `TRUE/FALSE` run flags in 
  `code/04-05`. Once large outputs are generated, later runs use `readRDS()`.
  <!-- TODO: note where a shared/backup copy of output/ lives, if any,
       so collaborators don't have to recompute everything from scratch -->

## Pipeline

Run the numbered scripts in `code/` in order:

| Script | Purpose |
|---|---|
| `00_config.R` | Loads packages, helper functions, and global settings. Sourced by every other script and not run directly. |
| `01_clean_data.R` | Cleans read-only data (clinical, targeted wastewater monitoring, and untargeted wastewater monitoring) from `data/original/` for use in the analysis. Saves output to `data/clean/`. |
| `02_subset_data.R` | Subsets cleaned data to be used in the analysis. Saves output to `data/clean/subsetted/`.|
| `03_group_data.R` | Groups data by virus and builds per-site data frames; saves output to `data/clean/subsetted/grouped/`. Also builds and saves input time series plots to `output/input_time_series/`. |
| `04_main_analyses.R` | Estimates time-varying Rt/rt and calculates concordance metrics. Saves results, output time series plots, and concordance figures/tables to `output/time_varying_estimates/`, `output/output/time_series/`, and `output/output_concordance/`. |
| `05_supplementary_analyses.R` | Supplementary analyses. Associated figures saved to `output/figures/`. |
| `06_final_figs.R` | Builds the final manuscript and supporting figures from the saved results in `output/` and saves to `output/figures/`. |

Other standalone scripts in `code/` are described below:

- `concordance_schematic.R`: Schematic depicting permutation used for concordance assessment. Can be run once the pipeline 
   above has produced its input.
- `format_big_tables.R`: Formats large concordance tables to be inserted into Supplementary Information. Can be run once the    pipeline above has produced its input.
- `map_figure.R`: Produces the map of study sites for Fig. 1.
- `shedding_pmfs.R`: Derives viral shedding probability mass functions (pmfs) used as EpiNow2 priors.

<!-- TODO: fill in author/contact and license info -->
