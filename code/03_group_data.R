# Group data sources by virus/site for analysis

source(here::here("code/00_config.R"))


# ----- Load data -----
wwscan_wwtp_dfs <- read_rds(here(paste0("data/clean/subsetted/wwscan_wwtps_viruses", nd_suffix, ".rds")))
casper_wwtp_dfs <- read_rds(here(paste0("data/clean/subsetted/casper_wwtps_viruses", nd_suffix, ".rds")))
epic_state_dfs <- read_rds(here("data/clean/subsetted/epic_state_viruses.rds"))
epic_county_dfs <- read_rds(here("data/clean/subsetted/epic_county_viruses.rds"))


# ---- Group and plot time series ----
# Fixed start/end date for each site's analysis window
site_date_ranges <- list(
  Boston    = list(min = "2024-09-01", max = "2026-03-31"),
  Boise     = list(min = "2024-11-01", max = "2026-03-31"),
  Ottumwa   = list(min = "2025-02-01", max = "2026-03-31"),
  Riverside = list(min = "2024-12-01", max = "2026-03-31")
)

#' Assemble a site's data sources into one named list, each trimmed to that site's fixed date range
#' @param ... Named data sources for the site (e.g. clinical = ..., targeted_ww = ..., untargeted_ww = ...); any subset of sources may be supplied
#' @param site Name of the site
#'
#' @return A named list matching the `...` names, with each df filtered to the site's date range
make_site_dfs <- function(..., site) {
  dfs_list <- list(...)

  range <- site_date_ranges[[site]]
  min_date <- as.Date(range$min)
  max_date <- as.Date(range$max)

  lapply(dfs_list, function(df) {
    df[df$date >= min_date & df$date <= max_date, ]
  })
}

#' Build a site's input time series plot and save 
#' @param site_dfs A list with elements clinical, targeted_ww, and/or untargeted_ww (as
#'   returned by make_site_dfs()); missing sources should be NULL/absent
#' @param filepath Output PNG path, relative to the project root; nd_suffix is inserted before the .png extension
#'   
#' @return NULL; saves the plot
plot_and_save <- function(site_dfs, filepath) {
  clinical_plot <- if (!is.null(site_dfs$clinical))
    plot_ts_panel(site_dfs$clinical, clinical_col, "Clinical", FALSE, "Virus signal")
  
  targeted_ww_plot <- if (!is.null(site_dfs$targeted_ww))
    plot_ts_panel(site_dfs$targeted_ww, targeted_ww_col, "Targeted ww", TRUE, "Virus signal")
  
  untargeted_ww_plot <- if (!is.null(site_dfs$untargeted_ww))
    plot_ts_panel(site_dfs$untargeted_ww, untargeted_ww_col, "Untargeted ww", TRUE, "Virus signal")
  
  all_dates <- do.call(c, lapply(
    Filter(Negate(is.null), list(site_dfs$clinical, site_dfs$targeted_ww, site_dfs$untargeted_ww)),
    `[[`, "date"
  ))
  x_limits <- range(all_dates)
  
  ts_plots <- list(clinical_plot, targeted_ww_plot, untargeted_ww_plot)
  ts_plots <- Filter(Negate(is.null), ts_plots)
  time_series_all <- Reduce(`/`, ts_plots) +
    plot_layout(guides = "collect",
                axes = "collect",
                axis_titles = "collect") &
    scale_x_date(date_labels = "%b %Y", limits = x_limits) &
    theme(legend.position = "bottom")
  
  filepath_nd <- sub("\\.png$", paste0(nd_suffix, ".png"), filepath)
  ggsave(
    here(filepath_nd),
    dpi    = 300,
    width  = 5,
    height = 3,
    units  = "in"
  )
}

## SARS-CoV-2 ----
### Boston, MA----
BOS_SC2_dfs <- make_site_dfs(
  clinical      = epic_county_dfs[["Suffolk_MA"]][["SC2"]],
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["SC2_N"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["SARS_CoV_2"]],
  site = "Boston"
)

plot_and_save(BOS_SC2_dfs, "output/input_time_series/SC2_Boston_MA.png")

### Boise, ID ----
boise_SC2_dfs <- make_site_dfs(
  clinical = epic_county_dfs[["Ada_ID"]][["SC2"]],
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["SC2_N"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["SARS_CoV_2"]],
  site = "Boise"
)

plot_and_save(boise_SC2_dfs, "output/input_time_series/SC2_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_SC2_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["IA"]][["SC2"]],
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["SC2_N"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["SARS_CoV_2"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_SC2_dfs, "output/input_time_series/SC2_Ottumwa_IA.png")

### Riverside, CA ----
riverside_SC2_dfs <- make_site_dfs(
  clinical = epic_county_dfs[["Riverside_CA"]][["SC2"]],
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["SC2_N"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["SARS_CoV_2"]],
  site = "Riverside"
)

plot_and_save(riverside_SC2_dfs, "output/input_time_series/SC2_Riverside_CA.png")

## IAV ----
### Boston, MA----
BOS_IAV_dfs <- make_site_dfs(
  clinical      = epic_county_dfs[["Suffolk_MA"]][["IAV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["Influenza_A"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Influenza_A_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_IAV_dfs, "output/input_time_series/IAV_Boston_MA.png")

### Boise, ID ----
boise_IAV_dfs <- make_site_dfs(
  clinical = epic_county_dfs[["Ada_ID"]][["IAV"]],
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["Influenza_A"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Influenza_A_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_IAV_dfs, "output/input_time_series/IAV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_IAV_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["IA"]][["IAV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["Influenza_A"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Influenza_A_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_IAV_dfs, "output/input_time_series/IAV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_IAV_dfs <- make_site_dfs(
  clinical = epic_county_dfs[["Riverside_CA"]][["IAV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["Influenza_A"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Influenza_A_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_IAV_dfs, "output/input_time_series/IAV_Riverside_CA.png")

## IBV ----
### Boston, MA----
BOS_IBV_dfs <- make_site_dfs(
  clinical      = epic_county_dfs[["Suffolk_MA"]][["IBV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["Influenza_B"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Influenza_B"]],
  site = "Boston"
)

plot_and_save(BOS_IBV_dfs, "output/input_time_series/IBV_Boston_MA.png")

### Boise, ID ----
boise_IBV_dfs <- make_site_dfs(
  clinical = epic_county_dfs[["Ada_ID"]][["IBV"]],
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["Influenza_B"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Influenza_B"]],
  site = "Boise"
)

plot_and_save(boise_IBV_dfs, "output/input_time_series/IBV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_IBV_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["IA"]][["IBV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["Influenza_B"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Influenza_B"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_IBV_dfs, "output/input_time_series/IBV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_IBV_dfs <- make_site_dfs(
  clinical = epic_county_dfs[["Riverside_CA"]][["IBV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["Influenza_B"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Influenza_B"]],
  site = "Riverside"
)

plot_and_save(riverside_IBV_dfs, "output/input_time_series/IBV_Riverside_CA.png")

## RSV ----
### Boston, MA----
BOS_RSV_dfs <- make_site_dfs(
  clinical      = epic_county_dfs[["Suffolk_MA"]][["RSV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["RSV"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["RSV_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_RSV_dfs, "output/input_time_series/RSV_Boston_MA.png")

### Boise, ID ----
boise_RSV_dfs <- make_site_dfs(
  clinical = epic_county_dfs[["Ada_ID"]][["RSV"]],
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["RSV"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["RSV_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_RSV_dfs, "output/input_time_series/RSV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_RSV_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["IA"]][["RSV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["RSV"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["RSV_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_RSV_dfs, "output/input_time_series/RSV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_RSV_dfs <- make_site_dfs(
  clinical = epic_county_dfs[["Riverside_CA"]][["RSV"]],
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["RSV"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["RSV_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_RSV_dfs, "output/input_time_series/RSV_Riverside_CA.png")

## Noro G2 ----
# Note: group NoV Epic data with NoV GII (GII is predominant genogroup)

### Boston, MA----
BOS_NoroG2_dfs <- make_site_dfs(
  clinical      = epic_state_dfs[["MA"]][["noro"]],
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["Noro_G2"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Norovirus_genogroup_II"]],
  site = "Boston"
)

plot_and_save(BOS_NoroG2_dfs, "output/input_time_series/NoroG2_Boston_MA.png")

### Boise, ID ----
boise_NoroG2_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["ID"]][["noro"]],
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["Noro_G2"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Norovirus_genogroup_II"]],
  site = "Boise"
)

plot_and_save(boise_NoroG2_dfs, "output/input_time_series/NoroG2_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_NoroG2_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["IA"]][["noro"]],
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["Noro_G2"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Norovirus_genogroup_II"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_NoroG2_dfs, "output/input_time_series/NoroG2_Ottumwa_IA.png")

### Riverside, CA ----
riverside_NoroG2_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["CA"]][["noro"]],
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["Noro_G2"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Norovirus_genogroup_II"]],
  site = "Riverside"
)

plot_and_save(riverside_NoroG2_dfs, "output/input_time_series/NoroG2_Riverside_CA.png")

## Norwalk ----
# Note: choose to group NoV Epic data with Noro GII (above)

### Boston, MA----
BOS_Norwalk_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Norwalk"]],
  site = "Boston"
)

plot_and_save(BOS_NoroG2_dfs, "output/input_time_series/Norwalk_Boston_MA.png")

### Boise, ID ----
boise_Norwalk_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Norwalk"]],
  site = "Boise"
)

plot_and_save(boise_NoroG2_dfs, "output/input_time_series/Norwalk_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Norwalk_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Norwalk"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_NoroG2_dfs, "output/input_time_series/Norwalk_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Norwalk_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Norwalk"]],
  site = "Riverside"
)

plot_and_save(riverside_NoroG2_dfs, "output/input_time_series/Norwalk_Riverside_CA.png")

## Rota ----
### Boston, MA----
BOS_Rota_dfs <- make_site_dfs(
  clinical      = epic_state_dfs[["MA"]][["rota"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Rotavirus_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_Rota_dfs, "output/input_time_series/Rota_Boston_MA.png")

### Boise, ID ----
boise_Rota_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["ID"]][["rota"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Rotavirus_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_Rota_dfs, "output/input_time_series/Rota_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Rota_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["IA"]][["rota"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Rotavirus_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Rota_dfs, "output/input_time_series/Rota_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Rota_dfs <- make_site_dfs(
  clinical = epic_state_dfs[["CA"]][["rota"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Rotavirus_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_Rota_dfs, "output/input_time_series/Rota_Riverside_CA.png")

## HMPV ----
### Boston, MA----
BOS_HMPV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["HMPV_4"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Metapneumovirus"]],
  site = "Boston"
)

plot_and_save(BOS_HMPV_dfs, "output/input_time_series/HMPV_Boston_MA.png")

### Boise, ID ----
boise_HMPV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["HMPV_4"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Metapneumovirus"]],
  site = "Boise"
)

plot_and_save(boise_HMPV_dfs, "output/input_time_series/HMPV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HMPV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["HMPV_4"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Metapneumovirus"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HMPV_dfs, "output/input_time_series/HMPV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HMPV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["HMPV_4"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Metapneumovirus"]],
  site = "Riverside"
)

plot_and_save(riverside_HMPV_dfs, "output/input_time_series/HMPV_Riverside_CA.png")

## EVD68 ----
### Boston, MA----
BOS_EVD68_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["EVD68"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Enterovirus_D68"]],
  site = "Boston"
)

plot_and_save(BOS_EVD68_dfs, "output/input_time_series/EVD68_Boston_MA.png")

### Boise, ID ----
boise_EVD68_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["EVD68"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Enterovirus_D68"]],
  site = "Boise"
)

plot_and_save(boise_EVD68_dfs, "output/input_time_series/EVD68_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_EVD68_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["EVD68"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Enterovirus_D68"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_EVD68_dfs, "output/input_time_series/EVD68_Ottumwa_IA.png")

### Riverside, CA ----
riverside_EVD68_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["EVD68"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Enterovirus_D68"]],
  site = "Riverside"
)

plot_and_save(riverside_EVD68_dfs, "output/input_time_series/EVD68_Riverside_CA.png")

## Parvo ----
### Boston, MA----
BOS_Parvo_B19_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Parvovirus_B19"]],
  site = "Boston"
)

plot_and_save(BOS_Parvo_B19_dfs, "output/input_time_series/Parvo_B19_Boston_MA.png")

### Boise, ID ----
boise_Parvo_B19_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Parvovirus_B19"]],
  site = "Boise"
)

plot_and_save(boise_Parvo_B19_dfs, "output/input_time_series/Parvo_B19_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Parvo_B19_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Parvovirus_B19"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Parvo_B19_dfs, "output/input_time_series/Parvo_B19_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Parvo_B19_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Parvovirus_B19"]],
  site = "Riverside"
)

plot_and_save(riverside_Parvo_B19_dfs, "output/input_time_series/Parvo_B19_Riverside_CA.png")

## HAdVF ----
### Boston, MA----
BOS_HAdVF_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_mastadenovirus_F"]],
  site = "Boston"
)

plot_and_save(BOS_HAdVF_dfs, "output/input_time_series/HAdVF_Boston_MA.png")

### Boise, ID ----
boise_HAdVF_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_mastadenovirus_F"]],
  site = "Boise"
)

plot_and_save(boise_HAdVF_dfs, "output/input_time_series/HAdVF_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HAdVF_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_mastadenovirus_F"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HAdVF_dfs, "output/input_time_series/HAdVF_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HAdVF_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_mastadenovirus_F"]],
  site = "Riverside"
)

plot_and_save(riverside_HAdVF_dfs, "output/input_time_series/HAdVF_Riverside_CA.png")

## HAV ----
### Boston, MA----
BOS_HAV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["HAV"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Hepatitis_A"]],
  site = "Boston"
)

plot_and_save(BOS_HAV_dfs, "output/input_time_series/HAV_Boston_MA.png")

### Boise, ID ----
boise_HAV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["HAV"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Hepatitis_A"]],
  site = "Boise"
)

plot_and_save(boise_HAV_dfs, "output/input_time_series/HAV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HAV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["HAV"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Hepatitis_A"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HAV_dfs, "output/input_time_series/HAV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HAV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["HAV"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Hepatitis_A"]],
  site = "Riverside"
)

plot_and_save(riverside_HAV_dfs, "output/input_time_series/HAV_Riverside_CA.png")

## HBV ----
### Boston, MA----
BOS_HBV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["HBV"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Hepatitis_B"]],
  site = "Boston"
)

plot_and_save(BOS_HBV_dfs, "output/input_time_series/HBV_Boston_MA.png")

### Boise, ID ----
boise_HBV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["HBV"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Hepatitis_B"]],
  site = "Boise"
)

plot_and_save(boise_HBV_dfs, "output/input_time_series/HBV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HBV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["HBV"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Hepatitis_B"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HBV_dfs, "output/input_time_series/HBV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HBV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["HBV"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Hepatitis_B"]],
  site = "Riverside"
)

plot_and_save(riverside_HBV_dfs, "output/input_time_series/HBV_Riverside_CA.png")

## HCV ----
### Boston, MA----
BOS_HCV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Boston_MA"]][["HCV"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Hepatitis_C"]],
  site = "Boston"
)

plot_and_save(BOS_HCV_dfs, "output/input_time_series/HCV_Boston_MA.png")

### Boise, ID ----
boise_HCV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["All_Boise"]][["HCV"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Hepatitis_C"]],
  site = "Boise"
)

plot_and_save(boise_HCV_dfs, "output/input_time_series/HCV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HCV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Ottumwa_IA"]][["HCV"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Hepatitis_C"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HCV_dfs, "output/input_time_series/HCV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HCV_dfs <- make_site_dfs(
  targeted_ww   = wwscan_wwtp_dfs[["Riverside_CA"]][["HCV"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Hepatitis_C"]],
  site = "Riverside"
)

plot_and_save(riverside_HCV_dfs, "output/input_time_series/HCV_Riverside_CA.png")

## Aichivirus ----
### Boston, MA----
BOS_Aichivirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Aichivirus"]],
  site = "Boston"
)

plot_and_save(BOS_Aichivirus_dfs, "output/input_time_series/Aichivirus_Boston_MA.png")

### Boise, ID ----
boise_Aichivirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Aichivirus"]],
  site = "Boise"
)

plot_and_save(boise_Aichivirus_dfs, "output/input_time_series/Aichivirus_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Aichivirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Aichivirus"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Aichivirus_dfs, "output/input_time_series/Aichivirus_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Aichivirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Aichivirus"]],
  site = "Riverside"
)

plot_and_save(riverside_Aichivirus_dfs, "output/input_time_series/Aichivirus_Riverside_CA.png")

## Astrovirus ----
### Boston, MA----
BOS_Astrovirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Astrovirus_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_Astrovirus_dfs, "output/input_time_series/Astrovirus_Boston_MA.png")

### Boise, ID ----
boise_Astrovirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Astrovirus_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_Astrovirus_dfs, "output/input_time_series/Astrovirus_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Astrovirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Astrovirus_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Astrovirus_dfs, "output/input_time_series/Astrovirus_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Astrovirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Astrovirus_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_Astrovirus_dfs, "output/input_time_series/Astrovirus_Riverside_CA.png")

## ---- Sapovirus ----
### Boston, MA----
BOS_Sapovirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Sapovirus"]],
  site = "Boston"
)

plot_and_save(BOS_Sapovirus_dfs, "output/input_time_series/Sapovirus_Boston_MA.png")

### Boise, ID ----
boise_Sapovirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Sapovirus"]],
  site = "Boise"
)

plot_and_save(boise_Sapovirus_dfs, "output/input_time_series/Sapovirus_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Sapovirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Sapovirus"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Sapovirus_dfs, "output/input_time_series/Sapovirus_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Sapovirus_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Sapovirus"]],
  site = "Riverside"
)

plot_and_save(riverside_Sapovirus_dfs, "output/input_time_series/Sapovirus_Riverside_CA.png")

## ---- Bocavirus ----
### Boston, MA----
BOS_HBoV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_bocavirus"]],
  site = "Boston"
)

plot_and_save(BOS_HBoV_dfs, "output/input_time_series/HBoV_Boston_MA.png")

### Boise, ID ----
boise_HBoV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_bocavirus"]],
  site = "Boise"
)

plot_and_save(boise_HBoV_dfs, "output/input_time_series/HBoV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HBoV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_bocavirus"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HBoV_dfs, "output/input_time_series/HBoV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HBoV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_bocavirus"]],
  site = "Riverside"
)

plot_and_save(riverside_HBoV_dfs, "output/input_time_series/HBoV_Riverside_CA.png")

## ---- IAV H1 ----
### Boston, MA----
BOS_IAVH1_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Boston_MA"]][["InfA_H1"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Influenza_A_H1N1"]],
  site = "Boston"
)

plot_and_save(BOS_IAVH1_dfs, "output/input_time_series/IAVH1_Boston_MA.png")

### Boise, ID ----
boise_IAVH1_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["All_Boise"]][["InfA_H1"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Influenza_A_H1N1"]],
  site = "Boise"
)

plot_and_save(boise_IAVH1_dfs, "output/input_time_series/IAVH1_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_IAVH1_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Ottumwa_IA"]][["InfA_H1"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Influenza_A_H1N1"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_IAVH1_dfs, "output/input_time_series/IAVH1_Ottumwa_IA.png")

### Riverside, CA ----
riverside_IAVH1_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Riverside_CA"]][["InfA_H1"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Influenza_A_H1N1"]],
  site = "Riverside"
)

plot_and_save(riverside_IAVH1_dfs, "output/input_time_series/IAVH1_Riverside_CA.png")

## ---- IAV H3 ----
### Boston, MA----
BOS_IAVH3_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Boston_MA"]][["InfA_H3"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Influenza_A_H3N1"]],
  site = "Boston"
)

plot_and_save(BOS_IAVH3_dfs, "output/input_time_series/IAVH3_Boston_MA.png")

### Boise, ID ----
boise_IAVH3_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["All_Boise"]][["InfA_H3"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Influenza_A_H3N1"]],
  site = "Boise"
)

plot_and_save(boise_IAVH3_dfs, "output/input_time_series/IAVH3_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_IAVH3_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Ottumwa_IA"]][["InfA_H3"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Influenza_A_H3N1"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_IAVH3_dfs, "output/input_time_series/IAVH3_Ottumwa_IA.png")

### Riverside, CA ----
riverside_IAVH3_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Riverside_CA"]][["InfA_H3"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Influenza_A_H3N1"]],
  site = "Riverside"
)

plot_and_save(riverside_IAVH3_dfs, "output/input_time_series/IAVH3_Riverside_CA.png")

## ---- IAV H5 ----
### Boston, MA----
BOS_IAVH5_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Boston_MA"]][["InfA_H5"]],
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Influenza_A_H5N1"]],
  site = "Boston"
)

plot_and_save(BOS_IAVH5_dfs, "output/input_time_series/IAVH5_Boston_MA.png")

### Boise, ID ----
boise_IAVH5_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["All_Boise"]][["InfA_H5"]],
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Influenza_A_H5N1"]],
  site = "Boise"
)

plot_and_save(boise_IAVH5_dfs, "output/input_time_series/IAVH5_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_IAVH5_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Ottumwa_IA"]][["InfA_H5"]],
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Influenza_A_H5N1"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_IAVH5_dfs, "output/input_time_series/IAVH5_Ottumwa_IA.png")

### Riverside, CA ----
riverside_IAVH5_dfs <- make_site_dfs(
  targeted_ww = wwscan_wwtp_dfs[["Riverside_CA"]][["InfA_H5"]],
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Influenza_A_H5N1"]],
  site = "Riverside"
)

plot_and_save(riverside_IAVH5_dfs, "output/input_time_series/IAVH5_Riverside_CA.png")

## ---- RSV A ----
### Boston, MA----
BOS_RSV_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["RSV_A"]],
  site = "Boston"
)

plot_and_save(BOS_RSV_A_dfs, "output/input_time_series/RSV_A_Boston_MA.png")

### Boise, ID ----
boise_RSV_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["RSV_A"]],
  site = "Boise"
)

plot_and_save(boise_RSV_A_dfs, "output/input_time_series/RSV_A_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_RSV_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["RSV_A"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_RSV_A_dfs, "output/input_time_series/RSV_A_Ottumwa_IA.png")

### Riverside, CA ----
riverside_RSV_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["RSV_A"]],
  site = "Riverside"
)

plot_and_save(riverside_RSV_A_dfs, "output/input_time_series/RSV_A_Riverside_CA.png")

## ---- RSV B ----
### Boston, MA----
BOS_RSV_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["RSV_B"]],
  site = "Boston"
)

plot_and_save(BOS_RSV_B_dfs, "output/input_time_series/RSV_B_Boston_MA.png")

### Boise, ID ----
boise_RSV_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["RSV_B"]],
  site = "Boise"
)

plot_and_save(boise_RSV_B_dfs, "output/input_time_series/RSV_B_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_RSV_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["RSV_B"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_RSV_B_dfs, "output/input_time_series/RSV_B_Ottumwa_IA.png")

### Riverside, CA ----
riverside_RSV_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["RSV_B"]],
  site = "Riverside"
)

plot_and_save(riverside_RSV_B_dfs, "output/input_time_series/RSV_B_Riverside_CA.png")

## ---- Rota A ----
### Boston, MA----
BOS_RotaA_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Rotavirus_A"]],
  site = "Boston"
)

plot_and_save(BOS_RotaA_dfs, "output/input_time_series/RotaA_Boston_MA.png")

### Boise, ID ----
boise_RotaA_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Rotavirus_A"]],
  site = "Boise"
)

plot_and_save(boise_RotaA_dfs, "output/input_time_series/RotaA_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_RotaA_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Rotavirus_A"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_RotaA_dfs, "output/input_time_series/RotaA_Ottumwa_IA.png")

### Riverside, CA ----
riverside_RotaA_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Rotavirus_A"]],
  site = "Riverside"
)

plot_and_save(riverside_RotaA_dfs, "output/input_time_series/RotaA_Riverside_CA.png")

## ---- Rota B ----
### Boston, MA----
BOS_RotaB_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Rotavirus_B"]],
  site = "Boston"
)

plot_and_save(BOS_RotaB_dfs, "output/input_time_series/RotaB_Boston_MA.png")

### Boise, ID ----
boise_RotaB_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Rotavirus_B"]],
  site = "Boise"
)

plot_and_save(boise_RotaB_dfs, "output/input_time_series/RotaB_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_RotaB_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Rotavirus_B"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_RotaB_dfs, "output/input_time_series/RotaB_Ottumwa_IA.png")

### Riverside, CA ----
riverside_RotaB_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Rotavirus_B"]],
  site = "Riverside"
)

plot_and_save(riverside_RotaB_dfs, "output/input_time_series/RotaB_Riverside_CA.png")

## ---- Rota C ----
### Boston, MA----
BOS_RotaC_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Rotavirus_C_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_RotaC_dfs, "output/input_time_series/RotaC_Boston_MA.png")

### Boise, ID ----
boise_RotaC_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Rotavirus_C_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_RotaC_dfs, "output/input_time_series/RotaC_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_RotaC_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Rotavirus_C_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_RotaC_dfs, "output/input_time_series/RotaC_Ottumwa_IA.png")

### Riverside, CA ----
riverside_RotaC_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Rotavirus_C_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_RotaC_dfs, "output/input_time_series/RotaC_Riverside_CA.png")

## ---- HAdV F40 ----
### Boston, MA----
BOS_HAdVF40_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_adenovirus_F40"]],
  site = "Boston"
)

plot_and_save(BOS_HAdVF40_dfs, "output/input_time_series/HAdVF40_Boston_MA.png")

### Boise, ID ----
boise_HAdVF40_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_adenovirus_F40"]],
  site = "Boise"
)

plot_and_save(boise_HAdVF40_dfs, "output/input_time_series/HAdVF40_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HAdVF40_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_adenovirus_F40"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HAdVF40_dfs, "output/input_time_series/HAdVF40_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HAdVF40_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_adenovirus_F40"]],
  site = "Riverside"
)

plot_and_save(riverside_HAdVF40_dfs, "output/input_time_series/HAdVF40_Riverside_CA.png")

## ---- HAdV F41 ----
### Boston, MA----
BOS_HAdVF41_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_adenovirus_F41"]],
  site = "Boston"
)

plot_and_save(BOS_HAdVF41_dfs, "output/input_time_series/HAdVF41_Boston_MA.png")

### Boise, ID ----
boise_HAdVF41_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_adenovirus_F41"]],
  site = "Boise"
)

plot_and_save(boise_HAdVF41_dfs, "output/input_time_series/HAdVF41_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HAdVF41_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_adenovirus_F41"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HAdVF41_dfs, "output/input_time_series/HAdVF41_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HAdVF41_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_adenovirus_F41"]],
  site = "Riverside"
)

plot_and_save(riverside_HAdVF41_dfs, "output/input_time_series/HAdVF41_Riverside_CA.png")

## ---- ICV ----
### Boston, MA----
BOS_ICV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Influenza_C"]],
  site = "Boston"
)

plot_and_save(BOS_ICV_dfs, "output/input_time_series/ICV_Boston_MA.png")

### Boise, ID ----
boise_ICV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Influenza_C"]],
  site = "Boise"
)

plot_and_save(boise_ICV_dfs, "output/input_time_series/ICV_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_ICV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Influenza_C"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_ICV_dfs, "output/input_time_series/ICV_Ottumwa_IA.png")

### Riverside, CA ----
riverside_ICV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Influenza_C"]],
  site = "Riverside"
)

plot_and_save(riverside_ICV_dfs, "output/input_time_series/ICV_Riverside_CA.png")

## ---- HPIV 1 ----
### Boston, MA----
BOS_HPIV1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Parainfluenza_virus_1"]],
  site = "Boston"
)

plot_and_save(BOS_HPIV1_dfs, "output/input_time_series/HPIV1_Boston_MA.png")

### Boise, ID ----
boise_HPIV1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Parainfluenza_virus_1"]],
  site = "Boise"
)

plot_and_save(boise_HPIV1_dfs, "output/input_time_series/HPIV1_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HPIV1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Parainfluenza_virus_1"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HPIV1_dfs, "output/input_time_series/HPIV1_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HPIV1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Parainfluenza_virus_1"]],
  site = "Riverside"
)

plot_and_save(riverside_HPIV1_dfs, "output/input_time_series/HPIV1_Riverside_CA.png")

## ---- HPIV 2 ----
### Boston, MA----
BOS_HPIV2_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Parainfluenza_virus_2"]],
  site = "Boston"
)

plot_and_save(BOS_HPIV2_dfs, "output/input_time_series/HPIV2_Boston_MA.png")

### Boise, ID ----
boise_HPIV2_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Parainfluenza_virus_2"]],
  site = "Boise"
)

plot_and_save(boise_HPIV2_dfs, "output/input_time_series/HPIV2_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HPIV2_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Parainfluenza_virus_2"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HPIV2_dfs, "output/input_time_series/HPIV2_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HPIV2_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Parainfluenza_virus_2"]],
  site = "Riverside"
)

plot_and_save(riverside_HPIV2_dfs, "output/input_time_series/HPIV2_Riverside_CA.png")

## ---- HPIV 3 ----
### Boston, MA----
BOS_HPIV3_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Parainfluenza_virus_3"]],
  site = "Boston"
)

plot_and_save(BOS_HPIV3_dfs, "output/input_time_series/HPIV3_Boston_MA.png")

### Boise, ID ----
boise_HPIV3_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Parainfluenza_virus_3"]],
  site = "Boise"
)

plot_and_save(boise_HPIV3_dfs, "output/input_time_series/HPIV3_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HPIV3_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Parainfluenza_virus_3"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HPIV3_dfs, "output/input_time_series/HPIV3_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HPIV3_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Parainfluenza_virus_3"]],
  site = "Riverside"
)

plot_and_save(riverside_HPIV3_dfs, "output/input_time_series/HPIV3_Riverside_CA.png")

## ---- HPIV 4 ----
### Boston, MA----
BOS_HPIV4_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Parainfluenza_virus_4"]],
  site = "Boston"
)

plot_and_save(BOS_HPIV4_dfs, "output/input_time_series/HPIV4_Boston_MA.png")

### Boise, ID ----
boise_HPIV4_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Parainfluenza_virus_4"]],
  site = "Boise"
)

plot_and_save(boise_HPIV4_dfs, "output/input_time_series/HPIV4_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HPIV4_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Parainfluenza_virus_4"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HPIV4_dfs, "output/input_time_series/HPIV4_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HPIV4_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Parainfluenza_virus_4"]],
  site = "Riverside"
)

plot_and_save(riverside_HPIV4_dfs, "output/input_time_series/HPIV4_Riverside_CA.png")

## ---- Noro GI ----
### Boston, MA----
BOS_NoroG1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Norovirus_genogroup_I_"]],
  site = "Boston"
)

plot_and_save(BOS_NoroG1_dfs, "output/input_time_series/NoroG1__Boston_MA.png")

### Boise, ID ----
boise_NoroG1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Norovirus_genogroup_I_"]],
  site = "Boise"
)

plot_and_save(boise_NoroG1_dfs, "output/input_time_series/NoroG1__Boise_ID.png")

### Ottumwa, IA ----
ottumwa_NoroG1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Norovirus_genogroup_I_"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_NoroG1_dfs, "output/input_time_series/NoroG1__Ottumwa_IA.png")

### Riverside, CA ----
riverside_NoroG1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Norovirus_genogroup_I_"]],
  site = "Riverside"
)

plot_and_save(riverside_NoroG1_dfs, "output/input_time_series/NoroG1__Riverside_CA.png")

## ---- Noro GIV ----
### Boston, MA----
BOS_NoroG4_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Norovirus_genogroup_IV"]],
  site = "Boston"
)

plot_and_save(BOS_NoroG4_dfs, "output/input_time_series/NoroG4_Boston_MA.png")

### Boise, ID ----
boise_NoroG4_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Norovirus_genogroup_IV"]],
  site = "Boise"
)

plot_and_save(boise_NoroG4_dfs, "output/input_time_series/NoroG4_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_NoroG4_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Norovirus_genogroup_IV"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_NoroG4_dfs, "output/input_time_series/NoroG4_Ottumwa_IA.png")

### Riverside, CA ----
riverside_NoroG4_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Norovirus_genogroup_IV"]],
  site = "Riverside"
)

plot_and_save(riverside_NoroG4_dfs, "output/input_time_series/NoroG4_Riverside_CA.png")

## ---- Rhino A ----
### Boston, MA----
BOS_Rhinovirus_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Rhinovirus_A"]],
  site = "Boston"
)

plot_and_save(BOS_Rhinovirus_A_dfs, "output/input_time_series/Rhinovirus_A_Boston_MA.png")

### Boise, ID ----
boise_Rhinovirus_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Rhinovirus_A"]],
  site = "Boise"
)

plot_and_save(boise_Rhinovirus_A_dfs, "output/input_time_series/Rhinovirus_A_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Rhinovirus_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Rhinovirus_A"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Rhinovirus_A_dfs, "output/input_time_series/Rhinovirus_A_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Rhinovirus_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Rhinovirus_A"]],
  site = "Riverside"
)

plot_and_save(riverside_Rhinovirus_A_dfs, "output/input_time_series/Rhinovirus_A_Riverside_CA.png")

## ---- Rhino B ----
### Boston, MA----
BOS_Rhinovirus_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Rhinovirus_B"]],
  site = "Boston"
)

plot_and_save(BOS_Rhinovirus_B_dfs, "output/input_time_series/Rhinovirus_B_Boston_MA.png")

### Boise, ID ----
boise_Rhinovirus_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Rhinovirus_B"]],
  site = "Boise"
)

plot_and_save(boise_Rhinovirus_B_dfs, "output/input_time_series/Rhinovirus_B_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Rhinovirus_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Rhinovirus_B"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Rhinovirus_B_dfs, "output/input_time_series/Rhinovirus_B_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Rhinovirus_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Rhinovirus_B"]],
  site = "Riverside"
)

plot_and_save(riverside_Rhinovirus_B_dfs, "output/input_time_series/Rhinovirus_B_Riverside_CA.png")


## ---- Rhino C ----
### Boston, MA----
BOS_Rhinovirus_C_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Rhinovirus_C_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_Rhinovirus_C_dfs, "output/input_time_series/Rhinovirus_C_Boston_MA.png")

### Boise, ID ----
boise_Rhinovirus_C_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Rhinovirus_C_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_Rhinovirus_C_dfs, "output/input_time_series/Rhinovirus_C_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Rhinovirus_C_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Rhinovirus_C_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Rhinovirus_C_dfs, "output/input_time_series/Rhinovirus_C_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Rhinovirus_C_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Rhinovirus_C_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_Rhinovirus_C_dfs, "output/input_time_series/Rhinovirus_C_Riverside_CA.png")

## ---- Corona 229E ----
### Boston, MA----
BOS_HCoV229E_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_coronavirus_229E"]],
  site = "Boston"
)

plot_and_save(BOS_HCoV229E_dfs, "output/input_time_series/HCoV229E_Boston_MA.png")

### Boise, ID ----
boise_HCoV229E_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_coronavirus_229E"]],
  site = "Boise"
)

plot_and_save(boise_HCoV229E_dfs, "output/input_time_series/HCoV229E_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HCoV229E_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_coronavirus_229E"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HCoV229E_dfs, "output/input_time_series/HCoV229E_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HCoV229E_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_coronavirus_229E"]],
  site = "Riverside"
)

plot_and_save(riverside_HCoV229E_dfs, "output/input_time_series/HCoV229E_Riverside_CA.png")

## ---- Corona HKU1 ----
### Boston, MA----
BOS_HCoVHKU1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_coronavirus_HKU1"]],
  site = "Boston"
)

plot_and_save(BOS_HCoVHKU1_dfs, "output/input_time_series/HCoVHKU1_Boston_MA.png")

### Boise, ID ----
boise_HCoVHKU1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_coronavirus_HKU1"]],
  site = "Boise"
)

plot_and_save(boise_HCoVHKU1_dfs, "output/input_time_series/HCoVHKU1_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HCoVHKU1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_coronavirus_HKU1"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HCoVHKU1_dfs, "output/input_time_series/HCoVHKU1_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HCoVHKU1_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_coronavirus_HKU1"]],
  site = "Riverside"
)

plot_and_save(riverside_HCoVHKU1_dfs, "output/input_time_series/HCoVHKU1_Riverside_CA.png")

## ---- Corona NL63 ----
### Boston, MA----
BOS_HCoVNL63_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_coronavirus_NL63"]],
  site = "Boston"
)

plot_and_save(BOS_HCoVNL63_dfs, "output/input_time_series/HCoVNL63_Boston_MA.png")

### Boise, ID ----
boise_HCoVNL63_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_coronavirus_NL63"]],
  site = "Boise"
)

plot_and_save(boise_HCoVNL63_dfs, "output/input_time_series/HCoVNL63_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HCoVNL63_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_coronavirus_NL63"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HCoVNL63_dfs, "output/input_time_series/HCoVNL63_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HCoVNL63_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_coronavirus_NL63"]],
  site = "Riverside"
)

plot_and_save(riverside_HCoVNL63_dfs, "output/input_time_series/HCoVNL63_Riverside_CA.png")

## ---- Corona OC43 ----
### Boston, MA----
BOS_HCoVOC43_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_coronavirus_OC43"]],
  site = "Boston"
)

plot_and_save(BOS_HCoVOC43_dfs, "output/input_time_series/HCoVOC43_Boston_MA.png")

### Boise, ID ----
boise_HCoVOC43_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_coronavirus_OC43"]],
  site = "Boise"
)

plot_and_save(boise_HCoVOC43_dfs, "output/input_time_series/HCoVOC43_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HCoVOC43_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_coronavirus_OC43"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HCoVOC43_dfs, "output/input_time_series/HCoVOC43_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HCoVOC43_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_coronavirus_OC43"]],
  site = "Riverside"
)

plot_and_save(riverside_HCoVOC43_dfs, "output/input_time_series/HCoVOC43_Riverside_CA.png")

## Enterovirus ----
### Boston, MA----
BOS_EV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Enterovirus_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_NoroG2_dfs, "output/input_time_series/Enterovirus_Boston_MA.png")

### Boise, ID ----
boise_EV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Enterovirus_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_NoroG2_dfs, "output/input_time_series/Enterovirus_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_EV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Enterovirus_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_NoroG2_dfs, "output/input_time_series/Enterovirus_Ottumwa_IA.png")

### Riverside, CA ----
riverside_EV_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Enterovirus_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_NoroG2_dfs, "output/input_time_series/Enterovirus_Riverside_CA.png")

## ---- Enterovirus A ----
### Boston, MA----
BOS_Enterovirus_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Enterovirus_A"]],
  site = "Boston"
)

plot_and_save(BOS_Enterovirus_A_dfs, "output/input_time_series/Enterovirus_A_Boston_MA.png")

### Boise, ID ----
boise_Enterovirus_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Enterovirus_A"]],
  site = "Boise"
)

plot_and_save(boise_Enterovirus_A_dfs, "output/input_time_series/Enterovirus_A_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Enterovirus_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Enterovirus_A"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Enterovirus_A_dfs, "output/input_time_series/Enterovirus_A_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Enterovirus_A_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Enterovirus_A"]],
  site = "Riverside"
)

plot_and_save(riverside_Enterovirus_A_dfs, "output/input_time_series/Enterovirus_A_Riverside_CA.png")

## ---- Enterovirus B ----
### Boston, MA----
BOS_Enterovirus_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Enterovirus_B"]],
  site = "Boston"
)

plot_and_save(BOS_Enterovirus_B_dfs, "output/input_time_series/Enterovirus_B_Boston_MA.png")

### Boise, ID ----
boise_Enterovirus_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Enterovirus_B"]],
  site = "Boise"
)

plot_and_save(boise_Enterovirus_B_dfs, "output/input_time_series/Enterovirus_B_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Enterovirus_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Enterovirus_B"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Enterovirus_B_dfs, "output/input_time_series/Enterovirus_B_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Enterovirus_B_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Enterovirus_B"]],
  site = "Riverside"
)

plot_and_save(riverside_Enterovirus_B_dfs, "output/input_time_series/Enterovirus_B_Riverside_CA.png")

## ---- Enterovirus C ----
### Boston, MA----
BOS_Enterovirus_C_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Enterovirus_C_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_Enterovirus_C_dfs, "output/input_time_series/Enterovirus_C_Boston_MA.png")

### Boise, ID ----
boise_Enterovirus_C_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Enterovirus_C_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_Enterovirus_C_dfs, "output/input_time_series/Enterovirus_C_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_Enterovirus_C_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Enterovirus_C_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_Enterovirus_C_dfs, "output/input_time_series/Enterovirus_C_Ottumwa_IA.png")

### Riverside, CA ----
riverside_Enterovirus_C_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Enterovirus_C_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_Enterovirus_C_dfs, "output/input_time_series/Enterovirus_C_Riverside_CA.png")

## ---- HAdV A ----
### Boston, MA----
BOS_HAdVA_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_mastadenovirus_A"]],
  site = "Boston"
)

plot_and_save(BOS_HAdVA_dfs, "output/input_time_series/HAdVA_Boston_MA.png")

### Boise, ID ----
boise_HAdVA_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_mastadenovirus_A"]],
  site = "Boise"
)

plot_and_save(boise_HAdVA_dfs, "output/input_time_series/HAdVA_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HAdVA_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_mastadenovirus_A"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HAdVA_dfs, "output/input_time_series/HAdVA_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HAdVA_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_mastadenovirus_A"]],
  site = "Riverside"
)

plot_and_save(riverside_HAdVA_dfs, "output/input_time_series/HAdVA_Riverside_CA.png")

## ---- HAdV B ----
### Boston, MA----
BOS_HAdVB_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_mastadenovirus_B"]],
  site = "Boston"
)

plot_and_save(BOS_HAdVB_dfs, "output/input_time_series/HAdVB_Boston_MA.png")

### Boise, ID ----
boise_HAdVB_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_mastadenovirus_B"]],
  site = "Boise"
)

plot_and_save(boise_HAdVB_dfs, "output/input_time_series/HAdVB_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HAdVB_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_mastadenovirus_B"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HAdVB_dfs, "output/input_time_series/HAdVB_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HAdVB_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_mastadenovirus_B"]],
  site = "Riverside"
)

plot_and_save(riverside_HAdVB_dfs, "output/input_time_series/HAdVB_Riverside_CA.png")

## ---- HAdV C ----
### Boston, MA----
BOS_HAdVC_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_mastadenovirus_C_clade_counts"]],
  site = "Boston"
)

plot_and_save(BOS_HAdVC_dfs, "output/input_time_series/HAdVC_Boston_MA.png")

### Boise, ID ----
boise_HAdVC_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_mastadenovirus_C_clade_counts"]],
  site = "Boise"
)

plot_and_save(boise_HAdVC_dfs, "output/input_time_series/HAdVC_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HAdVC_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_mastadenovirus_C_clade_counts"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HAdVC_dfs, "output/input_time_series/HAdVC_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HAdVC_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_mastadenovirus_C_clade_counts"]],
  site = "Riverside"
)

plot_and_save(riverside_HAdVC_dfs, "output/input_time_series/HAdVC_Riverside_CA.png")

## ---- HAdV D ----
### Boston, MA----
BOS_HAdVD_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boston_DITP_ALL_MA__MU_"]][["Human_mastadenovirus_D"]],
  site = "Boston"
)

plot_and_save(BOS_HAdVD_dfs, "output/input_time_series/HAdVD_Boston_MA.png")

### Boise, ID ----
boise_HAdVD_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Boise_RWF_ID"]][["Human_mastadenovirus_D"]],
  site = "Boise"
)

plot_and_save(boise_HAdVD_dfs, "output/input_time_series/HAdVD_Boise_ID.png")

### Ottumwa, IA ----
ottumwa_HAdVD_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Ottumwa_WPCF_IA"]][["Human_mastadenovirus_D"]],
  site = "Ottumwa"
)

plot_and_save(ottumwa_HAdVD_dfs, "output/input_time_series/HAdVD_Ottumwa_IA.png")

### Riverside, CA ----
riverside_HAdVD_dfs <- make_site_dfs(
  untargeted_ww = casper_wwtp_dfs[["Riverside_WQCP_CA"]][["Human_mastadenovirus_D"]],
  site = "Riverside"
)

plot_and_save(riverside_HAdVD_dfs, "output/input_time_series/HAdVD_Riverside_CA.png")

# ---- Export grouped lists ----
SC2_dfs <- list(
  Boston_MA = BOS_SC2_dfs,
  Boise_ID = boise_SC2_dfs,
  Ottumwa_IA = ottumwa_SC2_dfs,
  Riverside_CA = riverside_SC2_dfs
)

IAV_dfs <- list(
  Boston_MA = BOS_IAV_dfs,
  Boise_ID = boise_IAV_dfs,
  Ottumwa_IA = ottumwa_IAV_dfs,
  Riverside_CA = riverside_IAV_dfs
)

IBV_dfs <- list(
  Boston_MA = BOS_IBV_dfs,
  Boise_ID = boise_IBV_dfs,
  Ottumwa_IA = ottumwa_IBV_dfs,
  Riverside_CA = riverside_IBV_dfs
)

RSV_dfs <- list(
  Boston_MA = BOS_RSV_dfs,
  Boise_ID = boise_RSV_dfs,
  Ottumwa_IA = ottumwa_RSV_dfs,
  Riverside_CA = riverside_RSV_dfs
)

NoroG2_dfs <- list(
  Boston_MA = BOS_NoroG2_dfs,
  Boise_ID = boise_NoroG2_dfs,
  Ottumwa_IA = ottumwa_NoroG2_dfs,
  Riverside_CA = riverside_NoroG2_dfs
)

Norwalk_dfs <- list(
  Boston_MA = BOS_Norwalk_dfs,
  Boise_ID = boise_Norwalk_dfs,
  Ottumwa_IA = ottumwa_Norwalk_dfs,
  Riverside_CA = riverside_Norwalk_dfs
)

Rota_dfs <- list(
  Boston_MA = BOS_Rota_dfs,
  Boise_ID = boise_Rota_dfs,
  Ottumwa_IA = ottumwa_Rota_dfs,
  Riverside_CA = riverside_Rota_dfs
)

HMPV_dfs <- list(
  Boston_MA = BOS_HMPV_dfs,
  Boise_ID = boise_HMPV_dfs,
  Ottumwa_IA = ottumwa_HMPV_dfs,
  Riverside_CA = riverside_HMPV_dfs
)

EVD68_dfs <- list(
  Boston_MA = BOS_EVD68_dfs,
  Boise_ID = boise_EVD68_dfs,
  Ottumwa_IA = ottumwa_EVD68_dfs,
  Riverside_CA = riverside_EVD68_dfs
)

Parvo_B19_dfs <- list(
  Boston_MA = BOS_Parvo_B19_dfs,
  Boise_ID = boise_Parvo_B19_dfs,
  Ottumwa_IA = ottumwa_Parvo_B19_dfs,
  Riverside_CA = riverside_Parvo_B19_dfs
)

HAdVF_dfs <- list(
  Boston_MA = BOS_HAdVF_dfs,
  Boise_ID = boise_HAdVF_dfs,
  Ottumwa_IA = ottumwa_HAdVF_dfs,
  Riverside_CA = riverside_HAdVF_dfs
)

HAV_dfs <- list(
  Boston_MA = BOS_HAV_dfs,
  Boise_ID = boise_HAV_dfs,
  Ottumwa_IA = ottumwa_HAV_dfs,
  Riverside_CA = riverside_HAV_dfs
)

HBV_dfs <- list(
  Boston_MA = BOS_HBV_dfs,
  Boise_ID = boise_HBV_dfs,
  Ottumwa_IA = ottumwa_HBV_dfs,
  Riverside_CA = riverside_HBV_dfs
)

HCV_dfs <- list(
  Boston_MA = BOS_HCV_dfs,
  Boise_ID = boise_HCV_dfs,
  Ottumwa_IA = ottumwa_HCV_dfs,
  Riverside_CA = riverside_HCV_dfs
)

Aichivirus_dfs <- list(
  Boston_MA = BOS_Aichivirus_dfs,
  Boise_ID = boise_Aichivirus_dfs,
  Ottumwa_IA = ottumwa_Aichivirus_dfs,
  Riverside_CA = riverside_Aichivirus_dfs
)

Astrovirus_dfs <- list(
  Boston_MA = BOS_Astrovirus_dfs,
  Boise_ID = boise_Astrovirus_dfs,
  Ottumwa_IA = ottumwa_Astrovirus_dfs,
  Riverside_CA = riverside_Astrovirus_dfs
)

Sapovirus_dfs <- list(
  Boston_MA = BOS_Sapovirus_dfs,
  Boise_ID = boise_Sapovirus_dfs,
  Ottumwa_IA = ottumwa_Sapovirus_dfs,
  Riverside_CA = riverside_Sapovirus_dfs
)

HBoV_dfs <- list(
  Boston_MA = BOS_HBoV_dfs,
  Boise_ID = boise_HBoV_dfs,
  Ottumwa_IA = ottumwa_HBoV_dfs,
  Riverside_CA = riverside_HBoV_dfs
)

IAVH1_dfs <- list(
  Boston_MA = BOS_IAVH1_dfs,
  Boise_ID = boise_IAVH1_dfs,
  Ottumwa_IA = ottumwa_IAVH1_dfs,
  Riverside_CA = riverside_IAVH1_dfs
)

IAVH3_dfs <- list(
  Boston_MA = BOS_IAVH3_dfs,
  Boise_ID = boise_IAVH3_dfs,
  Ottumwa_IA = ottumwa_IAVH3_dfs,
  Riverside_CA = riverside_IAVH3_dfs
)

IAVH5_dfs <- list(
  Boston_MA = BOS_IAVH5_dfs,
  Boise_ID = boise_IAVH5_dfs,
  Ottumwa_IA = ottumwa_IAVH5_dfs,
  Riverside_CA = riverside_IAVH5_dfs
)

RSV_A_dfs <- list(
  Boston_MA = BOS_RSV_A_dfs,
  Boise_ID = boise_RSV_A_dfs,
  Ottumwa_IA = ottumwa_RSV_A_dfs,
  Riverside_CA = riverside_RSV_A_dfs
)

RSV_B_dfs <- list(
  Boston_MA = BOS_RSV_B_dfs,
  Boise_ID = boise_RSV_B_dfs,
  Ottumwa_IA = ottumwa_RSV_B_dfs,
  Riverside_CA = riverside_RSV_B_dfs
)

RotaA_dfs <- list(
  Boston_MA = BOS_RotaA_dfs,
  Boise_ID = boise_RotaA_dfs,
  Ottumwa_IA = ottumwa_RotaA_dfs,
  Riverside_CA = riverside_RotaA_dfs
)

RotaB_dfs <- list(
  Boston_MA = BOS_RotaB_dfs,
  Boise_ID = boise_RotaB_dfs,
  Ottumwa_IA = ottumwa_RotaB_dfs,
  Riverside_CA = riverside_RotaB_dfs
)

RotaC_dfs <- list(
  Boston_MA = BOS_RotaC_dfs,
  Boise_ID = boise_RotaC_dfs,
  Ottumwa_IA = ottumwa_RotaC_dfs,
  Riverside_CA = riverside_RotaC_dfs
)

HAdVF40_dfs <- list(
  Boston_MA = BOS_HAdVF40_dfs,
  Boise_ID = boise_HAdVF40_dfs,
  Ottumwa_IA = ottumwa_HAdVF40_dfs,
  Riverside_CA = riverside_HAdVF40_dfs
)

HAdVF41_dfs <- list(
  Boston_MA = BOS_HAdVF41_dfs,
  Boise_ID = boise_HAdVF41_dfs,
  Ottumwa_IA = ottumwa_HAdVF41_dfs,
  Riverside_CA = riverside_HAdVF41_dfs
)

ICV_dfs <- list(
  Boston_MA = BOS_ICV_dfs,
  Boise_ID = boise_ICV_dfs,
  Ottumwa_IA = ottumwa_ICV_dfs,
  Riverside_CA = riverside_ICV_dfs
)

HPIV1_dfs <- list(
  Boston_MA = BOS_HPIV1_dfs,
  Boise_ID = boise_HPIV1_dfs,
  Ottumwa_IA = ottumwa_HPIV1_dfs,
  Riverside_CA = riverside_HPIV1_dfs
)

HPIV2_dfs <- list(
  Boston_MA = BOS_HPIV2_dfs,
  Boise_ID = boise_HPIV2_dfs,
  Ottumwa_IA = ottumwa_HPIV2_dfs,
  Riverside_CA = riverside_HPIV2_dfs
)

HPIV3_dfs <- list(
  Boston_MA = BOS_HPIV3_dfs,
  Boise_ID = boise_HPIV3_dfs,
  Ottumwa_IA = ottumwa_HPIV3_dfs,
  Riverside_CA = riverside_HPIV3_dfs
)

HPIV4_dfs <- list(
  Boston_MA = BOS_HPIV4_dfs,
  Boise_ID = boise_HPIV4_dfs,
  Ottumwa_IA = ottumwa_HPIV4_dfs,
  Riverside_CA = riverside_HPIV4_dfs
)

NoroG1_dfs <- list(
  Boston_MA = BOS_NoroG1_dfs,
  Boise_ID = boise_NoroG1_dfs,
  Ottumwa_IA = ottumwa_NoroG1_dfs,
  Riverside_CA = riverside_NoroG1_dfs
)

NoroG4_dfs <- list(
  Boston_MA = BOS_NoroG4_dfs,
  Boise_ID = boise_NoroG4_dfs,
  Ottumwa_IA = ottumwa_NoroG4_dfs,
  Riverside_CA = riverside_NoroG4_dfs
)

Rhinovirus_A_dfs <- list(
  Boston_MA = BOS_Rhinovirus_A_dfs,
  Boise_ID = boise_Rhinovirus_A_dfs,
  Ottumwa_IA = ottumwa_Rhinovirus_A_dfs,
  Riverside_CA = riverside_Rhinovirus_A_dfs
)

Rhinovirus_B_dfs <- list(
  Boston_MA = BOS_Rhinovirus_B_dfs,
  Boise_ID = boise_Rhinovirus_B_dfs,
  Ottumwa_IA = ottumwa_Rhinovirus_B_dfs,
  Riverside_CA = riverside_Rhinovirus_B_dfs
)

Rhinovirus_C_dfs <- list(
  Boston_MA = BOS_Rhinovirus_C_dfs,
  Boise_ID = boise_Rhinovirus_C_dfs,
  Ottumwa_IA = ottumwa_Rhinovirus_C_dfs,
  Riverside_CA = riverside_Rhinovirus_C_dfs
)

HCoV229E_dfs <- list(
  Boston_MA = BOS_HCoV229E_dfs,
  Boise_ID = boise_HCoV229E_dfs,
  Ottumwa_IA = ottumwa_HCoV229E_dfs,
  Riverside_CA = riverside_HCoV229E_dfs
)

HCoVHKU1_dfs <- list(
  Boston_MA = BOS_HCoVHKU1_dfs,
  Boise_ID = boise_HCoVHKU1_dfs,
  Ottumwa_IA = ottumwa_HCoVHKU1_dfs,
  Riverside_CA = riverside_HCoVHKU1_dfs
)

HCoVNL63_dfs <- list(
  Boston_MA = BOS_HCoVNL63_dfs,
  Boise_ID = boise_HCoVNL63_dfs,
  Ottumwa_IA = ottumwa_HCoVNL63_dfs,
  Riverside_CA = riverside_HCoVNL63_dfs
)

HCoVOC43_dfs <- list(
  Boston_MA = BOS_HCoVOC43_dfs,
  Boise_ID = boise_HCoVOC43_dfs,
  Ottumwa_IA = ottumwa_HCoVOC43_dfs,
  Riverside_CA = riverside_HCoVOC43_dfs
)

Enterovirus_dfs <- list(
  Boston_MA = BOS_EV_dfs,
  Boise_ID = boise_EV_dfs,
  Ottumwa_IA = ottumwa_EV_dfs,
  Riverside_CA = riverside_EV_dfs
)

Enterovirus_A_dfs <- list(
  Boston_MA = BOS_Enterovirus_A_dfs,
  Boise_ID = boise_Enterovirus_A_dfs,
  Ottumwa_IA = ottumwa_Enterovirus_A_dfs,
  Riverside_CA = riverside_Enterovirus_A_dfs
)

Enterovirus_B_dfs <- list(
  Boston_MA = BOS_Enterovirus_B_dfs,
  Boise_ID = boise_Enterovirus_B_dfs,
  Ottumwa_IA = ottumwa_Enterovirus_B_dfs,
  Riverside_CA = riverside_Enterovirus_B_dfs
)

Enterovirus_C_dfs <- list(
  Boston_MA = BOS_Enterovirus_C_dfs,
  Boise_ID = boise_Enterovirus_C_dfs,
  Ottumwa_IA = ottumwa_Enterovirus_C_dfs,
  Riverside_CA = riverside_Enterovirus_C_dfs
)

HAdVA_dfs <- list(
  Boston_MA = BOS_HAdVA_dfs,
  Boise_ID = boise_HAdVA_dfs,
  Ottumwa_IA = ottumwa_HAdVA_dfs,
  Riverside_CA = riverside_HAdVA_dfs
)

HAdVB_dfs <- list(
  Boston_MA = BOS_HAdVB_dfs,
  Boise_ID = boise_HAdVB_dfs,
  Ottumwa_IA = ottumwa_HAdVB_dfs,
  Riverside_CA = riverside_HAdVB_dfs
)

HAdVC_dfs <- list(
  Boston_MA = BOS_HAdVC_dfs,
  Boise_ID = boise_HAdVC_dfs,
  Ottumwa_IA = ottumwa_HAdVC_dfs,
  Riverside_CA = riverside_HAdVC_dfs
)

HAdVD_dfs <- list(
  Boston_MA = BOS_HAdVD_dfs,
  Boise_ID = boise_HAdVD_dfs,
  Ottumwa_IA = ottumwa_HAdVD_dfs,
  Riverside_CA = riverside_HAdVD_dfs
)

all_virus_dfs <- list(
  SC2  = SC2_dfs,
  IAV  = IAV_dfs,
  IBV  = IBV_dfs,
  RSV  = RSV_dfs,
  NoroG2 = NoroG2_dfs,
  Norwalk = Norwalk_dfs,
  Rota = Rota_dfs,
  HMPV = HMPV_dfs,
  HAdVF = HAdVF_dfs,
  HAdVF40 = HAdVF40_dfs,
  HAdVF41 = HAdVF41_dfs,
  Aichivirus = Aichivirus_dfs,
  Astrovirus = Astrovirus_dfs,
  Sapovirus = Sapovirus_dfs,
  HBoV = HBoV_dfs,
  IAVH1 = IAVH1_dfs,
  IAVH3 = IAVH3_dfs,
  IAVH5 = IAVH5_dfs,
  RSV_A = RSV_A_dfs,
  RSV_B = RSV_B_dfs,
  RotaA = RotaA_dfs,
  RotaB = RotaB_dfs,
  RotaC = RotaC_dfs,
  ICV = ICV_dfs,
  HPIV1 = HPIV1_dfs,
  HPIV2 = HPIV2_dfs,
  HPIV3 = HPIV3_dfs,
  HPIV4 = HPIV4_dfs,
  NoroG1 = NoroG1_dfs,
  NoroG4 = NoroG4_dfs,
  Rhinovirus_A = Rhinovirus_A_dfs,
  Rhinovirus_B = Rhinovirus_B_dfs,
  Rhinovirus_C = Rhinovirus_C_dfs,
  HCoV229E = HCoV229E_dfs,
  HCoVHKU1 = HCoVHKU1_dfs,
  HCoVNL63 = HCoVNL63_dfs,
  HCoVOC43 = HCoVOC43_dfs,
  Enterovirus = Enterovirus_dfs,
  Enterovirus_A = Enterovirus_A_dfs,
  Enterovirus_B = Enterovirus_B_dfs,
  Enterovirus_C = Enterovirus_C_dfs,
  EVD68 = EVD68_dfs,
  HAdVA = HAdVA_dfs,
  HAdVB = HAdVB_dfs,
  HAdVC = HAdVC_dfs,
  HAdVD = HAdVD_dfs,
  Parvo_B19 = Parvo_B19_dfs,
  HAV = HAV_dfs,
  HBV = HBV_dfs,
  HCV = HCV_dfs
)
saveRDS(all_virus_dfs, file = here(paste0("data/clean/subsetted/grouped/all_virus_dfs", nd_suffix, ".rds")))
