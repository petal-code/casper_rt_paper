# Format big results tables

source(here::here("code/00_config.R"))

Rt_concordance <- read_csv(here(
  "output/output_concordance/Rt_ERN/Rt_concordance_table_NDs-imputed.csv"
))

growth_rate_concordance <- read_csv(here(
  "output/output_concordance/rt_growth_rate/rt_concordance_table_NDs-imputed.csv"
))

n_digits <- 4

Rt_table <- Rt_concordance %>%
  mutate(
    MAE_formatted = paste(
      format(round(obs_mae_median, n_digits), nsmall = n_digits),
      " (",
      format(round(obs_mae_lo95, n_digits), nsmall = n_digits),
      ", ",
      format(round(obs_mae_hi95, n_digits), nsmall = n_digits),
      ")",
      sep = ""
    ),
    MAE_p = format(round(mae_pval, n_digits), nsmall = n_digits),
    AOB_formatted = paste(
      format(round(obs_aob_median, n_digits), nsmall = n_digits),
      " (",
      format(round(obs_aob_lo95, n_digits), nsmall = n_digits),
      ", ",
      format(round(obs_aob_hi95, n_digits), nsmall = n_digits),
      ")",
      sep = ""
    ),
    AOB_p = format(round(aob_pval, n_digits), nsmall = n_digits),
    COR_formatted = paste(
      trimws(format(round(obs_cor_median, n_digits), nsmall = n_digits)),
      " (",
      trimws(format(round(obs_cor_lo95, n_digits), nsmall = n_digits)),
      ", ",
      trimws(format(round(obs_cor_hi95, n_digits), nsmall = n_digits)),
      ")",
      sep = ""
    ),
    COR_p = format(round(cor_pval, n_digits), nsmall = n_digits),
    comparison = case_when(
      comparison == "clinical vs targeted_ww" ~ "Clinical vs targeted wastewater",
      comparison == "clinical vs untargeted_ww" ~ "Clinical vs untargeted wastewater",
      comparison == "targeted_ww vs untargeted_ww" ~ "Targeted vs untargeted wastewater",
      T  ~ NA
    )
  ) %>% 
  select(virus,run,site,comparison,MAE_formatted,MAE_p,AOB_formatted,AOB_p,COR_formatted,COR_p)

clipr::write_clip(Rt_table) 


rt_table <- growth_rate_concordance %>%
  mutate(
    MAE_formatted = paste(
      format(round(obs_mae_median, n_digits), nsmall = n_digits),
      " (",
      format(round(obs_mae_lo95, n_digits), nsmall = n_digits),
      ", ",
      format(round(obs_mae_hi95, n_digits), nsmall = n_digits),
      ")",
      sep = ""
    ),
    MAE_p = format(round(mae_pval, n_digits), nsmall = n_digits),
    AOB_formatted = paste(
      format(round(obs_aob_median, n_digits), nsmall = n_digits),
      " (",
      format(round(obs_aob_lo95, n_digits), nsmall = n_digits),
      ", ",
      format(round(obs_aob_hi95, n_digits), nsmall = n_digits),
      ")",
      sep = ""
    ),
    AOB_p = format(round(aob_pval, n_digits), nsmall = n_digits),
    COR_formatted = paste(
      trimws(format(round(obs_cor_median, n_digits), nsmall = n_digits)),
      " (",
      trimws(format(round(obs_cor_lo95, n_digits), nsmall = n_digits)),
      ", ",
      trimws(format(round(obs_cor_hi95, n_digits), nsmall = n_digits)),
      ")",
      sep = ""
    ),
    COR_p = format(round(cor_pval, n_digits), nsmall = n_digits),
    comparison = case_when(
      comparison == "clinical vs targeted_ww" ~ "Clinical vs targeted wastewater",
      comparison == "clinical vs untargeted_ww" ~ "Clinical vs untargeted wastewater",
      comparison == "targeted_ww vs untargeted_ww" ~ "Targeted vs untargeted wastewater",
      T  ~ NA
    )
  ) %>% 
  select(virus,site,comparison,MAE_formatted,MAE_p,AOB_formatted,AOB_p,COR_formatted,COR_p)

clipr::write_clip(rt_table) 
