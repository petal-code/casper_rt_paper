# Project: Estimate rt and Rt using WW-MGS data from CASPER (+ clinical and targeted WW data) and assess concordance
# Author: Elana Chan
# Configuration file for project setup


# Libraries
renv::restore()

library(EpiNow2)
library(fs)
library(furrr)
library(future.apply)
library(ggh4x)
library(gratia)
library(here)
library(mgcv)
library(parallel)
library(patchwork)
library(readxl)
library(slider)
library(tidyverse)
library(truncdist)
library(zoo)

i_am("code/00_config.R")

# Functions
source(here("code/functions/compare_estimates.R"))
source(here("code/functions/estimation_functions.R"))
source(here("code/functions/plotting_helpers.R"))

# Global defaults
theme_set(theme_classic(base_size = 8))

options(scipen = 999)

clinical_col <- "#fc8d62"
targeted_ww_col <- "#8da0cb"
untargeted_ww_col <- "#66c2a5"
mae_col <- "#e78ac3"
aob_col <- "#a6d854"
cor_col <- "#ffd92f"

# Non-detect imputation: replace 0s with a small value before PMMoV normalization
impute_nondetects <- TRUE  # TRUE: WWSCAN 0s -> 500 gc/g; CASPER 0s -> 0.5 reads
nd_suffix <- if (impute_nondetects) "_NDs-imputed" else ""

# Reproducibility
seed <- 12345
set.seed(seed)

# EpiNow2 settings
settings_path <- here("documents/EpiNow2_run_settings.xlsx")
sheets <- excel_sheets(settings_path)
all_virus_settings <- lapply(sheets, function(x) read_excel(as.character(settings_path), sheet = x))
names(all_virus_settings) <- sheets
