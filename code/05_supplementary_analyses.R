# Supplementary analyses

source(here::here("code/00_config.R")) 
impute_nondetects <- TRUE # ensure set to T bc supplementary analyses are based on this
nd_suffix <- if (impute_nondetects) "_NDs-imputed" else ""

# Load cleaned dfs
all_virus_dfs <- read_rds(here(
  paste0(
    "data/clean/subsetted/grouped/all_virus_dfs",
    nd_suffix,
    ".rds"
  )
))


# ---- No ND substitution ----
# # Example: IBV, Boise, untargeted
# all_virus_dfs_no_impute <- read_rds(here("data/clean/subsetted/grouped/all_virus_dfs.rds"))
# Rt_no_impute <- read_rds(here("output/time_varying_estimates/IBV_run1_Rt.rds"))
# Rt_impute <- read_rds(here("output/time_varying_estimates/IBV_run1_NDs-imputed_Rt.rds"))
# 
# # Time series with/without substitution
# impute_rt <- fit_smooths(all_virus_dfs$IBV$Boise_ID$untargeted_ww, family = tw())
# impute_ts <- plot_ts_panel(all_virus_dfs$IBV$Boise_ID$untargeted_ww,
#                                       untargeted_ww_col,"With ND imputation",T,"IBV clade counts, normalized")
# impute_smooth <- plot_smooth_ts_panel(all_virus_dfs$IBV$Boise_ID$untargeted_ww,
#                                       impute_rt$gam_pred,
#                                       untargeted_ww_col,"With ND imputation",T,"IBV clade counts, normalized")
# 
# no_impute_rt <- fit_smooths(all_virus_dfs_no_impute$IBV$Boise_ID$untargeted_ww, family = tw())
# no_impute_ts <- plot_ts_panel(all_virus_dfs_no_impute$IBV$Boise_ID$untargeted_ww,
#                                          "grey40","Without ND imputation",T,"IBV clade counts, normalized")
# no_impute_smooth <- plot_smooth_ts_panel(all_virus_dfs_no_impute$IBV$Boise_ID$untargeted_ww,
#                                          no_impute_rt$gam_pred,
#                                               "grey40","Without ND imputation",T,"IBV clade counts, normalized")
# 
# # rt with/without substitution
# impute_rt_plot_data <- impute_rt$rt_pred %>% mutate(ND="With ND imputation")
# no_impute_rt_plot_data <- no_impute_rt$rt_pred %>% mutate(ND="Without ND imputation")
# rt_plot_data_NDs <- rbind(impute_rt_plot_data,no_impute_rt_plot_data)
# 
# rt_subplot_NDs <- ggplot(rt_plot_data_NDs,aes(x,y)) +
#   geom_hline(yintercept = 0, linetype="dashed",linewidth=0.5) +
#   geom_ribbon(aes(ymin = ymin,ymax = ymax,fill=ND),alpha=0.2) +
#   geom_line(aes(col=ND),linewidth=0.5) +
#   labs(x="",y="rt") +
#   scale_color_manual(values=c(untargeted_ww_col,"grey40")) +
#   scale_fill_manual(values=c(untargeted_ww_col,"grey40"))
# 
# # Rt with/without substitution
# impute_Rt_plot_data <- summary(Rt_impute$Boise_ID$untargeted_ww, type = "parameters", params = "R") %>% mutate(ND="With ND imputation")
# no_impute_Rt_plot_data <- summary(Rt_no_impute$Boise_ID$untargeted_ww, type = "parameters", params = "R") %>% mutate(ND="Without ND imputation")
# Rt_plot_data_NDs <- rbind(impute_Rt_plot_data,no_impute_Rt_plot_data)
# 
# Rt_subplot_NDs <- ggplot(Rt_plot_data_NDs,aes(date,median)) +
#   geom_hline(yintercept = 1,linetype="dashed",linewidth=0.5) +
#   geom_ribbon(aes(ymin = lower_90,ymax = upper_90,fill=ND),alpha=0.2)+
#   geom_line(aes(col=ND),linewidth=0.5) +
#   labs(x="",y="Rt") +
#   scale_color_manual(values=c(untargeted_ww_col,"grey40")) +
#   scale_fill_manual(values=c(untargeted_ww_col,"grey40"))
# 
# # assemble
# impute_smooth / no_impute_smooth / rt_subplot_NDs  +
#   plot_annotation(tag_levels = "a", tag_suffix = "  ") +
#   plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
#   labs(col="",fill="") &
#   theme(legend.position = "bottom")
# ggsave(
#   filename = here("output/figures/IBV_Boise_ID_growth-rate_NDs.png"),
#   width = 160, height = 150, units = "mm", dpi = 300)
# 
# impute_ts / no_impute_ts / Rt_subplot_NDs  +
#   plot_annotation(tag_levels = "a", tag_suffix = "  ") +
#   plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
#   labs(col="",fill="") &
#   theme(legend.position = "bottom")
# ggsave(
#   filename = here("output/figures/IBV_Boise_ID_ERN_NDs.png"),
#   width = 160, height = 150, units = "mm", dpi = 300)


# ---- GAM without outlier removal ----
# Example: SC2, Boston, targeted ww
sc2_bos_targeted_rt_main <- fit_smooths(all_virus_dfs$SC2$Boston_MA$targeted_ww, family = tw()) # as done in main analysis

# Modified function - no outlier removal
fit_smooths <- function(df, family) {
  
  family_str <- family$family
  
  is_ww <- family_str=="Tweedie"#"Gamma"
  
  if(is_ww && mean(df$is_nd==1)>0.75){
    list(gam_fit = NULL, gam_pred = NULL, rt_pred = NULL, rt_draws = NULL, k = NULL)
  }
  else{
    
    if(is_ww){
      # do nothing
    }
    
    # Date range
    date_seq <- seq(min(df$date), max(df$date), by = "day")
    
    n_samples <- length(df$date)
    k_max <- ceiling(n_samples/4)
    gam_fit <- gam(confirm ~ s(as.numeric(date), k=k_max),
                   method = "REML",
                   optimizer = c("outer", "bfgs"),
                   data = df,
                   family = family)
    newdata <- data.frame(date = date_seq)
    
    draws <- suppressWarnings(
      fitted_samples(
        gam_fit,
        data = newdata,
        n = 1000, # 1000 posterior draws
        seed = seed
      )
    )
    
    # Smoothed time series (median and 90% CrI for each date)
    gam_pred <- draws %>%
      group_by(.row) %>%  # date
      summarise(
        x    = date_seq[cur_group_id()],
        y    = median(.fitted),
        ymin = quantile(.fitted, 0.05),
        ymax = quantile(.fitted, 0.95)
      )
    
    # rt estimations (calculate within each draw, then tabulate median and 90% CrI across draws for each date)
    rt_draws <- draws %>%
      group_by(.draw) %>%
      arrange(.row, .by_group = TRUE) %>%
      mutate(
        x  = date_seq[.row],
        rt = c(NA, diff(log(.fitted))) 
      ) %>%
      ungroup() %>%
      filter(!is.na(rt))
    
    rt_pred <- rt_draws %>%
      group_by(.row) %>%
      summarise(
        x    = date_seq[.row[1]],
        y    = median(rt, na.rm = TRUE),
        ymin = quantile(rt, 0.05, na.rm = TRUE),
        ymax = quantile(rt, 0.95, na.rm = TRUE)
      )
    
    list(gam_fit = gam_fit, gam_pred = gam_pred, rt_pred = rt_pred, rt_draws = rt_draws, k = k_max)
  }
}

sc2_bos_targeted_rt_with_out <- fit_smooths(all_virus_dfs$SC2$Boston_MA$targeted_ww, family = tw()) # with outlier retained

sc2_bos_targeted_data <- plot_ts_panel(
  all_virus_dfs$SC2$Boston_MA$targeted_ww,
  targeted_ww_col,
  "Targeted ww",
  T,
  "SARS-CoV-2 RNA, normalized"
)

GAM_main <- plot_smooth_ts_panel(
  all_virus_dfs$SC2$Boston_MA$targeted_ww,
  sc2_bos_targeted_rt_main$gam_pred,
  targeted_ww_col,
  "Targeted ww",
  T,
  "SARS-CoV-2 RNA, normalized"
) +
  ylim(c(0, 0.002))

GAM_with_out <- plot_smooth_ts_panel(
  all_virus_dfs$SC2$Boston_MA$targeted_ww,
  sc2_bos_targeted_rt_with_out$gam_pred,
  targeted_ww_col,
  "Targeted ww",
  T,
  "SARS-CoV-2 RNA, normalized"
) +
  ylim(c(0, 0.002))

sc2_bos_targeted_data / GAM_main / GAM_with_out +
  plot_layout(guides = "collect",
              axes = "collect",
              axis_titles = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = "  ") &
  scale_x_date(breaks = "4 months", date_labels = "%b %Y")
ggsave(
  filename = here("output/figures/SC2_Boston_MA_GAM-outliers.png"),
  width = 160, height = 140, units = "mm", dpi = 300
)


# ---- Gamma instead of tweedie for ww GAM ----
# # Example: SC2, Boston, targeted ww
# 
# # Modified function - Gamma distribution instead of tw
# fit_smooths <- function(df, family) {
#   family_str <- family$family
#   
#   is_ww <- family_str == "Gamma"
#   
#   if (is_ww && mean(df$is_nd == 1) > 0.75) {
#     list(
#       gam_fit = NULL,
#       gam_pred = NULL,
#       rt_pred = NULL,
#       rt_draws = NULL,
#       k = NULL
#     )
#   }
#   else{
#     if (is_ww) {
#       df <- identify_ww_outlier(df)
#       df <- df %>% filter(is.na(is_outlier) | !is_outlier)
#     }
#     
#     # Date range
#     date_seq <- seq(min(df$date), max(df$date), by = "day")
#     
#     n_samples <- length(df$date)
#     k_max <- ceiling(n_samples / 4)
#     gam_fit <- gam(
#       confirm ~ s(as.numeric(date), k = k_max),
#       method = "REML",
#       optimizer = c("outer", "bfgs"),
#       data = df,
#       family = family
#     )
#     newdata <- data.frame(date = date_seq)
#     
#     draws <- suppressWarnings(fitted_samples(
#       gam_fit,
#       data = newdata,
#       n = 1000, # 1000 posterior draws
#       seed = seed
#     ))
#     
#     # Smoothed time series (median and 90% CrI for each date)
#     gam_pred <- draws %>%
#       group_by(.row) %>%  # date
#       summarise(
#         x    = date_seq[cur_group_id()],
#         y    = median(.fitted),
#         ymin = quantile(.fitted, 0.05),
#         ymax = quantile(.fitted, 0.95)
#       )
#     
#     # rt estimations (calculate within each draw, then tabulate median and 90% CrI across draws for each date)
#     rt_draws <- draws %>%
#       group_by(.draw) %>%
#       arrange(.row, .by_group = TRUE) %>%
#       mutate(
#         x  = date_seq[.row],
#         rt = c(NA, diff(log(.fitted))) 
#       ) %>%
#       ungroup() %>%
#       filter(!is.na(rt))
#     
#     rt_pred <- rt_draws %>%
#       group_by(.row) %>%
#       summarise(
#         x    = date_seq[.row[1]],
#         y    = median(rt, na.rm = TRUE),
#         ymin = quantile(rt, 0.05, na.rm = TRUE),
#         ymax = quantile(rt, 0.95, na.rm = TRUE)
#       )
#     
#     list(
#       gam_fit = gam_fit,
#       gam_pred = gam_pred,
#       rt_pred = rt_pred,
#       rt_draws = rt_draws,
#       k = k_max
#     )
#   }
# }
# 
# sc2_bos_targeted_rt_gamma <- fit_smooths(all_virus_dfs$SC2$Boston_MA$targeted_ww, family = stats::Gamma(link = "log"))
# 
# GAM_gamma <- plot_smooth_ts_panel(
#   all_virus_dfs$SC2$Boston_MA$targeted_ww,
#   sc2_bos_targeted_rt_gamma$gam_pred,
#   targeted_ww_col,
#   "Targeted ww",
#   T,
#   "SARS-CoV-2 RNA, normalized"
# ) +
#   ylim(c(0, 0.002))
# 
# GAM_main / GAM_gamma +
#   plot_layout(guides = "collect",
#               axes = "collect",
#               axis_titles = "collect") +
#   plot_annotation(tag_levels = "a", tag_suffix = "  ") &
#   scale_x_date(breaks = "4 months", date_labels = "%b %Y")
# ggsave(
#   filename = here("output/figures/SC2_Boston_MA_GAM-family.png"),
#   width = 160, height = 140, units = "mm", dpi = 300
# )


# ---- Rt with different multiplication factors ----
# Example: SC2, Boston, untargeted ww
sc2_bos_5x <- read_rds(here("output/time_varying_estimates/SC2_run1_NDs-imputed_Rt.rds"))
sc2_bos_5x <- sc2_bos_5x$Boston_MA$untargeted_ww
sc2_run1 <- parse_run_settings(all_virus_settings$`sars-cov-2`[1,])

# SC2_bos_untarget_6x<- fill_missing(
#   data = all_virus_dfs$SC2$Boston_MA$untargeted_ww %>%
#     select(date,confirm) %>%
#     mutate(confirm=round(confirm*10^6)),
#   missing_dates = "ignore"
# )
# sc2_bos_6x <- epinow(
#   data            = SC2_bos_untarget_6x,
#   generation_time = gt_opts(sc2_run1$generation_time),
#   delays          = delay_opts(dist = sc2_run1$shedding),
#   rt              = rt_opts(prior = sc2_run1$initial_Rt),
#   gp              = gp_opts(ls = sc2_run1$ls, alpha = sc2_run1$alpha),
#   obs             = obs_opts(
#     scale       = sc2_run1$obs_scale_ww,
#     dispersion  = sc2_run1$dispersion,
#     week_effect = sc2_run1$week_effect
#   ),
#   truncation = trunc_opts(sc2_run1$truncation),
#   stan       = stan_opts(cores = 1, chains = 4, seed = seed)
# )
# saveRDS(sc2_bos_6x, here("output/time_varying_estimates/SC2_run1_NDs-imputed_Bos_untargeted-6x_Rt.rds"))
sc2_bos_6x <- read_rds(here("output/time_varying_estimates/SC2_run1_NDs-imputed_Bos_untargeted-6x_Rt.rds"))
  
# SC2_bos_untarget_7x<- fill_missing(
#   data = all_virus_dfs$SC2$Boston_MA$untargeted_ww %>%
#     select(date,confirm) %>%
#     mutate(confirm=round(confirm*10^7)),
#   missing_dates = "ignore"
# )
# sc2_bos_7x <- epinow(
#   data            = SC2_bos_untarget_7x,
#   generation_time = gt_opts(sc2_run1$generation_time),
#   delays          = delay_opts(dist = sc2_run1$shedding),
#   rt              = rt_opts(prior = sc2_run1$initial_Rt),
#   gp              = gp_opts(ls = sc2_run1$ls, alpha = sc2_run1$alpha),
#   obs             = obs_opts(
#     scale       = sc2_run1$obs_scale_ww,
#     dispersion  = sc2_run1$dispersion,
#     week_effect = sc2_run1$week_effect
#   ),
#   truncation = trunc_opts(sc2_run1$truncation),
#   stan       = stan_opts(samples = 4000, warmup = 1000, # increase MCMC samples so CrIs smooth
#                          cores = 1, chains = 4, seed = seed)
# )
# saveRDS(sc2_bos_7x, here("output/time_varying_estimates/SC2_run1_NDs-imputed_Bos_untargeted-7x_Rt.rds"))
sc2_bos_7x <- read_rds(here("output/time_varying_estimates/SC2_run1_NDs-imputed_Bos_untargeted-7x_Rt.rds"))

# SC2_bos_untarget_8x<- fill_missing(
#   data = all_virus_dfs$SC2$Boston_MA$untargeted_ww %>%
#     select(date,confirm) %>%
#     mutate(confirm=round(confirm*10^8)),
#   missing_dates = "ignore"
# )
# sc2_bos_8x <- epinow(
#   data            = SC2_bos_untarget_8x,
#   generation_time = gt_opts(sc2_run1$generation_time),
#   delays          = delay_opts(dist = sc2_run1$shedding),
#   rt              = rt_opts(prior = sc2_run1$initial_Rt),
#   gp              = gp_opts(ls = sc2_run1$ls, alpha = sc2_run1$alpha),
#   obs             = obs_opts(
#     scale       = sc2_run1$obs_scale_ww,
#     dispersion  = sc2_run1$dispersion,
#     week_effect = sc2_run1$week_effect
#   ),
#   truncation = trunc_opts(sc2_run1$truncation),
#   stan       = stan_opts(samples = 4000, warmup = 1000, # increase MCMC samples so CrIs smooth
#                          cores = 1, chains = 4, seed = seed)
# )
# saveRDS(sc2_bos_8x, here("output/time_varying_estimates/SC2_run1_NDs-imputed_Bos_untargeted-8x_Rt.rds"))
sc2_bos_8x <- read_rds(here("output/time_varying_estimates/SC2_run1_NDs-imputed_Bos_untargeted-8x_Rt.rds"))

# Plot
sc2_bos_untarget_ts <- plot_ts_panel(all_virus_dfs$SC2$Boston_MA$untargeted_ww,
                                     untargeted_ww_col,
                                     "Untargeted ww",T,"SARS-CoV-2, normalized")
 
mult5x_Rt_plot_data <- summary(sc2_bos_5x, type = "parameters", params = "R") %>% 
  mutate(mult="10^5")
mult6x_Rt_plot_data <- summary(sc2_bos_6x, type = "parameters", params = "R") %>% 
  mutate(mult="10^6")
mult7x_Rt_plot_data <- summary(sc2_bos_7x, type = "parameters", params = "R") %>% 
  mutate(mult="10^7")
mult8x_Rt_plot_data <- summary(sc2_bos_8x, type = "parameters", params = "R") %>% 
  mutate(mult="10^8")
Rt_mult_plot_data <- rbind(mult5x_Rt_plot_data,mult6x_Rt_plot_data,
                           mult7x_Rt_plot_data,mult8x_Rt_plot_data)

Rt_mult_subplot <- ggplot(Rt_mult_plot_data,aes(date,median)) +
  geom_hline(yintercept = 1,linetype="dashed",linewidth=1) +
  geom_ribbon(aes(ymin = lower_90,ymax = upper_90,fill=mult),alpha=0.2)+
  geom_line(aes(col=mult),linewidth=1) +
  labs(x="",y="Rt") +
  scale_color_viridis_d(option = "H") +
  scale_fill_viridis_d(option = "H") +
  scale_y_log10()

sc2_bos_untarget_ts / Rt_mult_subplot +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect",axes = "collect",axis_titles = "collect") &
  theme(legend.position = "bottom") &
  labs(col="",fill="") &
  scale_x_date(breaks = "4 months", date_labels = "%b %Y")
ggsave(
  filename = here("output/figures/SC2_Boston_MA_Rt-mult.png"),
  width = 160, height = 140, units = "mm", dpi = 300
)


# ---- Downsampled targeted ww data ----
# Example: SC2, Boston, targeted ww
nondownsampled_Rt <- read_rds(here("output/time_varying_estimates/SC2_run1_NDs-imputed_Rt.rds"))
nondownsampled_Rt <- nondownsampled_Rt$Boston_MA$targeted_ww

downsampled <- all_virus_dfs$SC2$Boston_MA$targeted_ww %>%
  mutate(date = floor_date(date, unit = "week", week_start = "Sunday") + 6) %>%
  group_by(date) %>%
  slice_sample(n = 1) %>%
  ungroup()
# downsampled_fill <- fill_missing(
#   data = downsampled %>% select(date,confirm) %>% mutate(confirm = round(confirm * 10^6)),
#   missing_dates = "ignore"
# )
# downsampled_Rt <- epinow(
#   data            = downsampled_fill,
#   generation_time = gt_opts(sc2_run1$generation_time),
#   delays          = delay_opts(dist = sc2_run1$shedding),
#   rt              = rt_opts(prior = sc2_run1$initial_Rt),
#   gp              = gp_opts(ls = sc2_run1$ls, alpha = sc2_run1$alpha),
#   obs             = obs_opts(
#     scale       = sc2_run1$obs_scale_ww,
#     dispersion  = sc2_run1$dispersion,
#     week_effect = sc2_run1$week_effect
#   ),
#   truncation = trunc_opts(sc2_run1$truncation),
#   stan       = stan_opts(cores = 1, chains = 4, seed = seed)
# )
# saveRDS(downsampled_Rt, here("output/time_varying_estimates/SC2_run1_NDs-imputed_Bos_targeted-downsampled_Rt.rds"))
downsampled_Rt <- read_rds(here("output/time_varying_estimates/SC2_run1_NDs-imputed_Bos_targeted-downsampled_Rt.rds"))

# resource original rt functions to estimate rt
source(here("code/functions/estimation_functions.R"))

downsampled_growth_rate <- fit_smooths(downsampled, family = tw())
nondownsampled_growth_rate <- fit_smooths(all_virus_dfs$SC2$Boston_MA$targeted_ww, family = tw())

# Plot
nondownsampled_smooth <- plot_smooth_ts_panel(all_virus_dfs$SC2$Boston_MA$targeted_ww,
                                       nondownsampled_growth_rate$gam_pred,
                                       targeted_ww_col,"Targeted ww - 3x/week",T,"SARS-CoV-2 RNA, normalized", plot_out=F) 
nondownsampled_ts <- plot_ts_panel(all_virus_dfs$SC2$Boston_MA$targeted_ww,
                                              targeted_ww_col,"Targeted ww - 3x/week",T,"SARS-CoV-2 RNA, normalized", plot_out=F) 

downsampled_smooth <- plot_smooth_ts_panel(downsampled,
                                       downsampled_growth_rate$gam_pred,
                                       "grey40","Targeted ww - 1x/week",T,"SARS-CoV-2 RNA, normalized")
downsampled_ts <- plot_ts_panel(downsampled,
                                           "grey40","Targeted ww - 1x/week",T,"SARS-CoV-2 RNA, normalized")

downsampled_Rt_plot_data <- summary(downsampled_Rt, type = "parameters", params = "R") %>% 
  mutate(freq="1 sample/week")
nondownsampled_Rt_plot_data <- summary(nondownsampled_Rt, type = "parameters", params = "R") %>% mutate(freq="3 samples/week")
Rt_plot_data <- rbind(downsampled_Rt_plot_data,nondownsampled_Rt_plot_data)

Rt_subplot <- ggplot(Rt_plot_data,aes(date,median)) +
  geom_hline(yintercept = 1,linetype="dashed",linewidth=0.5) +
  geom_ribbon(aes(ymin = lower_90,ymax = upper_90,fill=freq),alpha=0.2)+
  geom_line(aes(col=freq),linewidth=0.5) +
  labs(x="",y="Rt") +
  scale_color_manual(values=c("grey40",targeted_ww_col)) +
  scale_fill_manual(values=c("grey40",targeted_ww_col))

downsampled_rt_plot_data <- downsampled_growth_rate$rt_pred %>% mutate(freq="1 sample/week")
nondownsampled_rt_plot_data <- nondownsampled_growth_rate$rt_pred %>% mutate(freq="3 samples/week")
rt_plot_data <- rbind(downsampled_rt_plot_data,nondownsampled_rt_plot_data)

rt_subplot <- ggplot(rt_plot_data,aes(x,y)) +
  geom_hline(yintercept = 0, linetype="dashed",linewidth=0.5) +
  geom_ribbon(aes(ymin = ymin,ymax = ymax,fill=freq),alpha=0.2) +
  geom_line(aes(col=freq),linewidth=0.5) +
  labs(x="",y="rt") +
  scale_color_manual(values=c("grey40",targeted_ww_col)) +
  scale_fill_manual(values=c("grey40",targeted_ww_col))

# rt
nondownsampled_smooth / downsampled_smooth / rt_subplot  +
  plot_annotation(tag_levels = "a", tag_suffix = "  ") +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(col="",fill="") &
  theme(legend.position = "bottom")
ggsave(
  filename = here("output/figures/SC2_Boston_MA_growth-rate_freq.png"),
  width = 160, height = 150, units = "mm", dpi = 300
)

# Rt
nondownsampled_ts / downsampled_ts / Rt_subplot +
  plot_annotation(tag_levels = "a", tag_suffix = "  ") +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(col="",fill="") &
  theme(legend.position = "bottom")
ggsave(
  filename = here("output/figures/SC2_Boston_MA_ERN_freq.png"),
  width = 160, height = 150, units = "mm", dpi = 300
)
