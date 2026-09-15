# Manuscript figures

source(here::here("code/00_config.R"))
impute_nondetects <- TRUE # ensure set to T bc main analyses are based on this
nd_suffix <- if (impute_nondetects) "_NDs-imputed" else ""


# ---- Load data and results ----
all_virus_dfs <- read_rds(here(paste0("data/clean/subsetted/grouped/all_virus_dfs", nd_suffix, ".rds")))
all_rt_estimates <- read_rds(here("output/time_varying_estimates", paste0("all_rt", nd_suffix, ".rds")))
Rt_runs <- map(set_names(c("SC2", "IAV", "IBV", "RSV", "NoroG2", "Rota")), function(virus) {
  read_rds(here(paste0("output/time_varying_estimates/", virus, "_run1", nd_suffix,"_Rt.rds")))
})
Rt_runs$IBV$Ottumwa_IA$untargeted_ww <- NULL # IBV/Boise/untargeted data are mostly NDs --> don't actually use Rt estimates
rt_df <- read_csv(here(
  paste0(
    "output/output_concordance/rt_growth_rate/rt_concordance_table",
    nd_suffix,
    ".csv"
  )
))
ern_df <- read_csv(here(
  paste0(
    "output/output_concordance/Rt_ERN/Rt_concordance_table",
    nd_suffix,
    ".csv"
  )
))


# ---- Shared helpers ----
# Y-axis label per virus, and a)-d) panel titles keyed by 
# virus (for the per-virus 2x2 grids) or by site (for the per-site 2x2 grids)
virus_ylabs <- c(SC2 = "SARS-CoV-2 signal", IAV = "IAV signal", IBV = "IBV signal",
                  RSV = "RSV signal", NoroG2 = "NoV signal", Rota = "Rotavirus signal",
                  IAVH1 = "IAV H1N1 signal", IAVH3 = "IAV H3N1 signal", IAVH5 = "IAV H5N1 signal",
                  HMPV = "HMPV signal", EVD68 = "EVD68 signal", HAV = "HAV signal")

virus_panel_titles <- c(
  SC2 = "SARS-CoV-2",
  IAV = "Influenza A virus", IBV = "Influenza B virus",
  RSV = "Respiratory syncytial virus", NoroG2 = "Norovirus", Rota = "Rotavirus")

site_panel_titles <- c(Boise_ID = "a) Boise, ID", Boston_MA = "b) Boston, MA",
                        Ottumwa_IA = "c) Ottumwa, IA", Riverside_CA = "d) Riverside, CA")

#' Build one time series panel (one virus at one site) for the 2x2 manuscript figure grids --
#' either the EpiNow2 Rt estimate or the GAM-derived rt estimate. 
#' @param virus Virus key into all_virus_dfs / all_rt_estimates / Rt_runs (e.g. "SC2")
#' @param site Site key into all_virus_dfs[[virus]] / all_rt_estimates[[virus]] / Rt_runs[[virus]]
#' @param title Panel title (e.g. "a) Boise, ID" or "a) Influenza A virus")
#' @param type "Rt" (default) for the EpiNow2 reproduction number (from Rt_runs); "rt" for the
#'   GAM-derived growth rate (from all_rt_estimates)
#'
#' @return A patchwork plot, as built by plot_time_varying_estimates()
build_manuscript_estimate_panel <- function(virus, site, title, type = c("Rt", "rt")) {
  type <- match.arg(type)
  site_smooths <- all_rt_estimates[[virus]][[site]]

  if (type == "Rt") {
    estimate_data <- extract_Rt_data(Rt_runs[[virus]][[site]])
    estimate_ylab <- expression(R[t])
    yintercept    <- 1
    
    plot_time_varying_estimates(
      dfs           = all_virus_dfs[[virus]][[site]],
      estimate_data = estimate_data,
      ts_ylab       = virus_ylabs[[virus]],
      estimate_ylab = estimate_ylab,
      yintercept    = yintercept,
      yrange        = NULL,
      plot_out      = FALSE,
      title         = title
    )
  } else {
    estimate_data <- extract_rt_data(site_smooths)
    estimate_ylab <- expression(r[t])
    yintercept    <- 0
    
    plot_time_varying_estimates(
      dfs           = all_virus_dfs[[virus]][[site]],
      smooths       = lapply(site_smooths, `[[`, "gam_pred"),
      estimate_data = estimate_data,
      ts_ylab       = virus_ylabs[[virus]],
      estimate_ylab = estimate_ylab,
      yintercept    = yintercept,
      yrange        = NULL,
      plot_out      = FALSE,
      title         = title
    )
  }
}

#' Combine 4 panels into the standard 2x2 manuscript grid (shared legend, bold titles) and save
#' @param panels List of 4 panels in a)-d) reading order (top-left, top-right, bottom-left,
#'   bottom-right)
#' @param filename Output PNG filename, saved under output/figures/
#'
#' @return NULL; saves the combined figure
save_grid_figure <- function(panels, filename) {
  grid <- ((panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]])) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  ggsave(here("output/figures", filename), plot = grid,
         width = 160, height = 170, units = "mm", dpi = 400)
}


# ---- Fig 2: SC2 Rt, all sites ----
SC2_Rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("SC2", .y, .x, type = "Rt"))

# Annotate summer 2025 peaks
annotate_SC2_peak_Rt <- function(site_SC2_run, site_SC2_panel){
  peak_df <- extract_Rt_data(site_SC2_run) %>% 
    filter(x >= as.Date("2025-06-01") & x <= as.Date("2025-09-30")) %>% 
    group_by(source) %>% 
    summarise(max_y_x = x[which.max(y)],
              max_y = max(y))
  
  site_SC2_panel[[4]] <- site_SC2_panel[[4]] + geom_point(data = peak_df,
                                                          aes(x = max_y_x,y=max_y),pch = 17,size=2,guide="none")
  
  print(peak_df)
  
  return(site_SC2_panel[[4]])
}

SC2_Rt_panels$Boise_ID[[4]] <- annotate_SC2_peak_Rt(Rt_runs$SC2$Boise_ID,SC2_Rt_panels$Boise_ID) 
SC2_Rt_panels$Boston_MA[[4]] <- annotate_SC2_peak_Rt(Rt_runs$SC2$Boston_MA,SC2_Rt_panels$Boston_MA) 
SC2_Rt_panels$Ottumwa_IA[[4]] <- annotate_SC2_peak_Rt(Rt_runs$SC2$Ottumwa_IA,SC2_Rt_panels$Ottumwa_IA) 
SC2_Rt_panels$Riverside_CA[[4]] <- annotate_SC2_peak_Rt(Rt_runs$SC2$Riverside_CA,SC2_Rt_panels$Riverside_CA) 

# Make manual adjustments - adjust date breaks for Ottumwa
SC2_Rt_panels$Ottumwa_IA <- SC2_Rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")

save_grid_figure(SC2_Rt_panels, "Fig2.png")

## SI: Other viruses ----
IAV_Rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("IAV", .y, .x, type = "Rt"))
IAV_Rt_panels$Ottumwa_IA <- IAV_Rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(IAV_Rt_panels, "Fig2_SI_IAV.png")

IBV_Rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("IBV", .y, .x, type = "Rt"))
IBV_Rt_panels$Ottumwa_IA <- IBV_Rt_panels$Ottumwa_IA + guides(col="none",fill="none")
IBV_Rt_panels$Ottumwa_IA <- IBV_Rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
IBV_Rt_panels$Ottumwa_IA[[3]] <- IBV_Rt_panels$Ottumwa_IA[[3]] + scale_y_discrete(labels = "+/-")
IBV_Rt_panels$Boston_MA[[4]] <- IBV_Rt_panels$Boston_MA[[4]] + ylim(c(0.5,1.75)) # yscale
save_grid_figure(IBV_Rt_panels, "Fig2_SI_IBV.png")

RSV_Rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("RSV", .y, .x, type = "Rt"))
RSV_Rt_panels$Ottumwa_IA <- RSV_Rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(RSV_Rt_panels, "Fig2_SI_RSV.png")

Rota_Rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("Rota", .y, .x, type = "Rt"))
Rota_Rt_panels$Ottumwa_IA <- Rota_Rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(Rota_Rt_panels, "Fig2_SI_Rota.png")

NoroG2_Rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("NoroG2", .y, .x, type = "Rt"))
NoroG2_Rt_panels$Ottumwa_IA <- NoroG2_Rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(NoroG2_Rt_panels, "Fig2_SI_NoroG2.png")


# ---- Fig 3: All Rt, Bos ----
Bos_Rt_panels <- imap(virus_panel_titles, ~ build_manuscript_estimate_panel(.y, "Boston_MA", .x, type = "Rt"))

# Calculate peak winter 2025-25
compute_winter_peak_df <- function(site_virus_run){
  extract_Rt_data(site_virus_run) %>%
    filter(x >= as.Date("2024-09-01") & x <= as.Date("2025-03-01")) %>%
    group_by(source) %>%
    slice_max(y, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(source, max_y_x = x, max_y = y, max_y_min = ymin, max_y_max = ymax)
}

SC2_peak_Rt_Bos    <- compute_winter_peak_df(Rt_runs$SC2$Boston_MA)
IAV_peak_Rt_Bos    <- compute_winter_peak_df(Rt_runs$IAV$Boston_MA)
IBV_peak_Rt_Bos    <- compute_winter_peak_df(Rt_runs$IBV$Boston_MA)
RSV_peak_Rt_Bos    <- compute_winter_peak_df(Rt_runs$RSV$Boston_MA)
NoroG2_peak_Rt_Bos <- compute_winter_peak_df(Rt_runs$NoroG2$Boston_MA)
Rota_peak_Rt_Bos   <- compute_winter_peak_df(Rt_runs$Rota$Boston_MA)

# Annotate peak winter 2024-25
annotate_winter_peak <- function(peak_df, site_virus_panel){
  site_virus_panel + geom_point(data = peak_df,
                                 aes(x = max_y_x,y=max_y),pch = 17,size=2)
}

Bos_Rt_panels$SC2[[4]] <- annotate_winter_peak(SC2_peak_Rt_Bos, Bos_Rt_panels$SC2[[4]])
Bos_Rt_panels$IAV[[4]] <- annotate_winter_peak(IAV_peak_Rt_Bos, Bos_Rt_panels$IAV[[4]])
Bos_Rt_panels$IBV[[4]] <- annotate_winter_peak(IBV_peak_Rt_Bos, Bos_Rt_panels$IBV[[4]])
Bos_Rt_panels$RSV[[4]] <- annotate_winter_peak(RSV_peak_Rt_Bos, Bos_Rt_panels$RSV[[4]])
Bos_Rt_panels$NoroG2[[4]] <- annotate_winter_peak(NoroG2_peak_Rt_Bos, Bos_Rt_panels$NoroG2[[4]])
Bos_Rt_panels$Rota[[3]] <- annotate_winter_peak(Rota_peak_Rt_Bos, Bos_Rt_panels$Rota[[3]])

# Manual adjustments
Bos_Rt_panels$SC2[[4]] <- Bos_Rt_panels$SC2[[4]] +   guides(fill=guide_legend(nrow = 2),col=guide_legend(nrow=2))
Bos_Rt_panels$IAV[[4]] <- Bos_Rt_panels$IAV[[4]] +   guides(fill=guide_legend(nrow = 2),col=guide_legend(nrow=2))
Bos_Rt_panels$IBV[[4]] <- Bos_Rt_panels$IBV[[4]] +   guides(fill=guide_legend(nrow = 2),col=guide_legend(nrow=2)) & ylim(c(0.5,1.75)) # yscale
Bos_Rt_panels$RSV[[4]] <- Bos_Rt_panels$RSV[[4]] +   guides(fill=guide_legend(nrow = 2),col=guide_legend(nrow=2))
Bos_Rt_panels$NoroG2[[4]] <- Bos_Rt_panels$NoroG2[[4]] +   guides(fill=guide_legend(nrow = 2),col=guide_legend(nrow=2))
Bos_Rt_panels$Rota[[3]] <- Bos_Rt_panels$Rota[[3]] + guides(col="none",fill="none",pch="none") # so legend appears ok

# Add labels
row_title_theme <- theme(plot.title = element_text(size = 7, face = "plain", hjust = 0))
Bos_Rt_panels$SC2[[4]]    <- Bos_Rt_panels$SC2[[4]]    + labs(title = "SARS-CoV-2") + row_title_theme
Bos_Rt_panels$IAV[[4]]    <- Bos_Rt_panels$IAV[[4]]    + labs(title = "IAV") + row_title_theme
Bos_Rt_panels$IBV[[4]]    <- Bos_Rt_panels$IBV[[4]]    + labs(title = "IBV") + row_title_theme
Bos_Rt_panels$RSV[[4]]    <- Bos_Rt_panels$RSV[[4]]    + labs(title = "RSV") + row_title_theme
Bos_Rt_panels$NoroG2[[4]] <- Bos_Rt_panels$NoroG2[[4]] + labs(title = "Norovirus") + row_title_theme
Bos_Rt_panels$Rota[[3]]   <- Bos_Rt_panels$Rota[[3]]   + labs(title = "Rotavirus") + row_title_theme

# Assemble panel a
Fig3_a <- Bos_Rt_panels$SC2[[4]] / Bos_Rt_panels$IAV[[4]] / Bos_Rt_panels$IBV[[4]] / Bos_Rt_panels$RSV[[4]] /
  Bos_Rt_panels$NoroG2[[4]] / Bos_Rt_panels$Rota[[3]] +
  plot_annotation(title = "a) Rₜ estimates: Boston, MA",
                  theme = theme(plot.title = element_text(face = "bold", size = 10))) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") &
  theme(legend.position = "bottom",
        )

# Panel b: magnitude of peak Rt
peak_rt_df <- bind_rows(
  SC2_peak_Rt_Bos    %>% mutate(virus = "SARS-CoV-2"),
  IAV_peak_Rt_Bos    %>% mutate(virus = "IAV"),
  IBV_peak_Rt_Bos    %>% mutate(virus = "IBV"),
  RSV_peak_Rt_Bos    %>% mutate(virus = "RSV"),
  NoroG2_peak_Rt_Bos %>% mutate(virus = "Norovirus"),
  Rota_peak_Rt_Bos   %>% mutate(virus = "Rotavirus")
) %>%
  mutate(virus = factor(virus, levels = c("SARS-CoV-2","IAV","IBV","RSV","Norovirus","Rotavirus")))

peak_rt_df <- peak_rt_df %>% 
  complete(virus,source) # necessary for plot dodging below

peak_rt_df_date <- peak_rt_df %>% 
  group_by(virus) %>% 
  summarise(avg_max_y_x = mean(max_y_x,na.rm=T))

peak_rt_df_mag <- peak_rt_df %>% 
  group_by(virus) %>% 
  summarise(avg_max_y = mean(max_y,na.rm=T))

peak_dodge <- position_dodge(width = 0.7, preserve = "single")

peak_rt_panel <- ggplot(peak_rt_df, aes(x = virus, y = max_y, col=source)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_pointrange(aes(ymin = max_y_min, ymax = max_y_max),
                position = peak_dodge,
                linewidth = 0.5,pch = 17) +
  geom_point(data = peak_rt_df_mag,aes(virus,avg_max_y),col="black",pch = 8) +
  scale_color_manual(values = setNames(c(clinical_col, targeted_ww_col, untargeted_ww_col), data_source_labels)) +
  labs(x = "", y = expression(paste("2024-25 peak season ", R[t], " (90% CrI)")), fill = "",
       title = "c) Magnitude of seasonal peak Rₜ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 10))

peak_date_rt_panel <- ggplot(peak_rt_df, aes(x = virus, y = max_y_x, col=source)) +
  geom_hline(yintercept = as.Date("2025-01-01"), linetype = "dashed", linewidth = 0.5) +
  geom_point(position = peak_dodge,pch = 17,size=2) +
  geom_point(data = peak_rt_df_date,aes(virus,avg_max_y_x),col="black",pch = 8) +
  scale_color_manual(values = setNames(c(clinical_col, targeted_ww_col, untargeted_ww_col), data_source_labels)) +
  labs(y = expression(paste("2024-25 peak season ", R[t])), x = "", 
       title = "b) Timing of seasonal peak Rₜ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 10))

# Assemble all
(wrap_elements(full = Fig3_a) | (peak_date_rt_panel / peak_rt_panel)) +
  plot_layout(widths = c(1.5, 1), guides = "collect", axes = "collect",axis_titles = "collect")

ggsave(here("output/figures/Fig3.png"),width = 7, height = 6.3, units = "in", dpi = 400)

## SI: rt ----
# all viruses with 2+ data sources (where most are above LOD)
### SC2 ----
SC2_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("SC2", .y, .x, type = "rt"))
SC2_rt_panels$Ottumwa_IA <- SC2_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(SC2_rt_panels, "Fig2_SI_rt_SC2.png")

compute_summer_peak_df <- function(site_virus_rt){
  extract_rt_data(site_virus_rt) %>%
    filter(x >= as.Date("2025-06-01") & x <= as.Date("2025-09-30")) %>% 
    group_by(source) %>%
    slice_max(y, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(source, max_y_x = x, max_y = y, max_y_min = ymin, max_y_max = ymax)
}
compute_summer_peak_df(all_rt_estimates$SC2$Boise_ID)
compute_summer_peak_df(all_rt_estimates$SC2$Boston_MA)
compute_summer_peak_df(all_rt_estimates$SC2$Ottumwa_IA)
compute_summer_peak_df(all_rt_estimates$SC2$Riverside_CA)

### IAV ----
compute_winter_peak_rt_df <- function(site_virus_rt){
  extract_rt_data(site_virus_rt) %>%
    filter(x >= as.Date("2024-09-01") & x <= as.Date("2025-03-01")) %>%
    group_by(source) %>%
    slice_max(y, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(source, max_y_x = x, max_y = y, max_y_min = ymin, max_y_max = ymax)
}

IAV_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("IAV", .y, .x, type = "rt"))
IAV_rt_panels$Ottumwa_IA <- IAV_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(IAV_rt_panels, "Fig2_SI_rt_IAV.png")

compute_winter_peak_rt_df(all_rt_estimates$IAV$Boise_ID)
compute_winter_peak_rt_df(all_rt_estimates$IAV$Boston_MA)
compute_winter_peak_rt_df(all_rt_estimates$IAV$Ottumwa_IA)
compute_winter_peak_rt_df(all_rt_estimates$IAV$Riverside_CA)

### IBV ----
IBV_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("IBV", .y, .x, type = "rt"))
IBV_rt_panels$Ottumwa_IA <- IBV_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
IBV_rt_panels$Ottumwa_IA <- IBV_rt_panels$Ottumwa_IA + guides(col = "none", fill = "none")
IBV_rt_panels$Ottumwa_IA[[3]] <- IBV_rt_panels$Ottumwa_IA[[3]] + scale_y_discrete(labels = "+/-")
save_grid_figure(IBV_rt_panels, "Fig2_SI_rt_IBV.png")

compute_winter_peak_rt_df(all_rt_estimates$IBV$Boise_ID)
compute_winter_peak_rt_df(all_rt_estimates$IBV$Boston_MA)
compute_winter_peak_rt_df(all_rt_estimates$IBV$Ottumwa_IA)
compute_winter_peak_rt_df(all_rt_estimates$IBV$Riverside_CA)

### RSV ----
RSV_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("RSV", .y, .x, type = "rt"))
RSV_rt_panels$Ottumwa_IA <- RSV_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(RSV_rt_panels, "Fig2_SI_rt_RSV.png")

compute_winter_peak_rt_df(all_rt_estimates$RSV$Boise_ID)
compute_winter_peak_rt_df(all_rt_estimates$RSV$Boston_MA)
compute_winter_peak_rt_df(all_rt_estimates$RSV$Ottumwa_IA)
compute_winter_peak_rt_df(all_rt_estimates$RSV$Riverside_CA)

### NoroG2 ----
NoroG2_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("NoroG2", .y, .x, type = "rt"))
NoroG2_rt_panels$Ottumwa_IA <- NoroG2_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(NoroG2_rt_panels, "Fig2_SI_rt_NoV.png")

compute_winter_peak_rt_df(all_rt_estimates$NoroG2$Boise_ID)
compute_winter_peak_rt_df(all_rt_estimates$NoroG2$Boston_MA)
compute_winter_peak_rt_df(all_rt_estimates$NoroG2$Ottumwa_IA)
compute_winter_peak_rt_df(all_rt_estimates$NoroG2$Riverside_CA)

### Rota ----
Rota_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("Rota", .y, .x, type = "rt"))
Rota_rt_panels$Ottumwa_IA <- Rota_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
save_grid_figure(Rota_rt_panels, "Fig2_SI_rt_Rota.png")

compute_winter_peak_rt_df(all_rt_estimates$Rota$Boise_ID)
compute_winter_peak_rt_df(all_rt_estimates$Rota$Boston_MA)
compute_winter_peak_rt_df(all_rt_estimates$Rota$Ottumwa_IA)
compute_winter_peak_rt_df(all_rt_estimates$Rota$Riverside_CA)

### IAV H1 ----
IAVH1_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("IAVH1", .y, .x, type = "rt"))
IAVH1_rt_panels$Ottumwa_IA <- IAVH1_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
IAVH1_rt_panels$Riverside_CA <- IAVH1_rt_panels$Riverside_CA + guides(col = "none", fill = "none")
IAVH1_rt_panels$Riverside_CA[[2]] <- IAVH1_rt_panels$Riverside_CA[[2]] + scale_y_discrete(labels = "+/-")
save_grid_figure(IAVH1_rt_panels, "Fig2_SI_rt_IAVH1.png")

### HMPV ----
HMPV_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("HMPV", .y, .x, type = "rt"))
HMPV_rt_panels$Ottumwa_IA <- HMPV_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
HMPV_rt_panels$Ottumwa_IA <- HMPV_rt_panels$Ottumwa_IA + guides(col = "none", fill = "none")
HMPV_rt_panels$Ottumwa_IA[[2]] <- HMPV_rt_panels$Ottumwa_IA[[2]] + scale_y_discrete(labels = "+/-")
HMPV_rt_panels$Riverside_CA <- HMPV_rt_panels$Riverside_CA + guides(col = "none", fill = "none")
HMPV_rt_panels$Riverside_CA[[1]] <- HMPV_rt_panels$Riverside_CA[[1]] + scale_y_discrete(labels = "+/-")
save_grid_figure(HMPV_rt_panels, "Fig2_SI_rt_HMPV.png")

### HAV ----
HAV_rt_panels <- imap(site_panel_titles, ~ build_manuscript_estimate_panel("HAV", .y, .x, type = "rt"))
HAV_rt_panels$Ottumwa_IA <- HAV_rt_panels$Ottumwa_IA &
  scale_x_date(date_breaks = "4 months",
               limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
               date_labels = "%b %Y")
HAV_rt_panels$Boise_ID <- HAV_rt_panels$Boise_ID + guides(col = "none", fill = "none")
HAV_rt_panels$Boise_ID[[1]] <- HAV_rt_panels$Boise_ID[[1]] + scale_y_discrete(labels = "+/-")
HAV_rt_panels$Boise_ID[[2]] <- HAV_rt_panels$Boise_ID[[2]] + scale_y_discrete(labels = "+/-")
HAV_rt_panels$Ottumwa_IA <- HAV_rt_panels$Ottumwa_IA + guides(col = "none", fill = "none")
HAV_rt_panels$Ottumwa_IA[[1]] <- HAV_rt_panels$Ottumwa_IA[[1]] + scale_y_discrete(labels = "+/-")
HAV_rt_panels$Riverside_CA <- HAV_rt_panels$Riverside_CA + guides(col = "none", fill = "none")
HAV_rt_panels$Riverside_CA[[1]] <- HAV_rt_panels$Riverside_CA[[1]] + scale_y_discrete(labels = "+/-")
save_grid_figure(HAV_rt_panels, "Fig2_SI_rt_HAV.png")


# ---- Fig 5: rt for virus groupings: Bos ----
extract_med_rt <- function(virus_list) {
  imap_dfr(virus_list, function(virus_data, virus_name) {
    rt <- virus_data[["Boston_MA"]][["untargeted_ww"]][["rt_pred"]]
    if (is.null(rt) || nrow(rt) == 0) return(NULL)
    
    # bucket each date into a season based on month, ignoring year
    rt <- rt %>%
      mutate(season = case_when(
        lubridate::month(x) %in% c(6, 7, 8)   ~ "jun_aug",
        lubridate::month(x) %in% c(9, 10, 11) ~ "sep_nov",
        lubridate::month(x) %in% c(12, 1, 2)  ~ "dec_feb",
        lubridate::month(x) %in% c(3, 4, 5)   ~ "mar_may"
      ))
    
    tibble(
      virus = virus_name,
      # med rt during Jun-Aug
      med_rt_jun_aug = median(rt$y[rt$season == "jun_aug"], na.rm = TRUE),
      # med rt during Sep-Nov
      med_rt_sep_nov = median(rt$y[rt$season == "sep_nov"], na.rm = TRUE),
      # med rt during Dec-Feb
      med_rt_dec_feb = median(rt$y[rt$season == "dec_feb"], na.rm = TRUE),
      # med rt during Mar-May
      med_rt_mar_may = median(rt$y[rt$season == "mar_may"], na.rm = TRUE),
      # q1 rt during Jun-Aug
      q1_rt_jun_aug = quantile(rt$y[rt$season == "jun_aug"], 0.25,na.rm = TRUE),
      # q1 rt during Sep-Nov
      q1_rt_sep_nov = quantile(rt$y[rt$season == "sep_nov"], 0.25,na.rm = TRUE),
      # q1 rt during Dec-Feb
      q1_rt_dec_feb = quantile(rt$y[rt$season == "dec_feb"], 0.25,na.rm = TRUE),
      # q1 rt during Mar-May
      q1_rt_mar_may = quantile(rt$y[rt$season == "mar_may"], 0.25,na.rm = TRUE),
      # q3 rt during Jun-Aug
      q3_rt_jun_aug = quantile(rt$y[rt$season == "jun_aug"], 0.75,na.rm = TRUE),
      # q3 rt during Sep-Nov
      q3_rt_sep_nov = quantile(rt$y[rt$season == "sep_nov"], 0.75,na.rm = TRUE),
      # q3 rt during Dec-Feb
      q3_rt_dec_feb = quantile(rt$y[rt$season == "dec_feb"], 0.75,na.rm = TRUE),
      # q3 rt during Mar-May
      q3_rt_mar_may = quantile(rt$y[rt$season == "mar_may"], 0.75,na.rm = TRUE)
    )
  })
}

# extract_peak_rt <- function(virus_list) {
#   imap_dfr(virus_list, function(virus_data, virus_name) {
#     rt <- virus_data[["Boston_MA"]][["untargeted_ww"]][["rt_pred"]]
#     if (is.null(rt) || nrow(rt) == 0) return(NULL)
#     
#     # just 2025
#     rt <- rt %>%
#       filter(year(x)==2025)
#     
#     tibble(
#       virus = virus_name,
#       max_rt = max(rt$y),
#       max_rt_min = max(rt$ymin),
#       max_rt_max = max(rt$ymax),
#       max_rt_date = rt$x[which.max(rt$y)]
#     )
#   })
# }

med_rt_untargeted <- extract_med_rt(all_rt_estimates) %>% 
  pivot_longer(
    cols = matches("^(med_rt|q1_rt|q3_rt)_"),
    names_to = c(".value", "season"),
    names_pattern = "(med_rt|q1_rt|q3_rt)_(.*)"
  )

med_rt_untargeted_group <- med_rt_untargeted %>% 
  filter(!(virus %in% c("Norwalk","Enterovirus","Rota","RSV","IAV","HAdVF"))) %>% 
  mutate(virus_group = case_when(
    str_starts(virus,"Enterovirus") ~ "Enterovirus*",
    str_starts(virus,"HPIV") ~ "Human parainfluenza virus",
    str_starts(virus,"HAdV") ~ "Human adenovirus*",
    str_starts(virus,"HCoV") | virus == "SC2" ~ "Human coronavirus*",
    str_starts(virus,"Rota") ~ "Rotavirus*",
    str_starts(virus,"Noro") ~ "Norovirus",
    str_starts(virus,"RSV") ~ "Respiratory syncytial virus*",
    virus %in% c("HAV","HBV","HCV") ~ "Hepatitis virus*",
    str_starts(virus,"Rhino") ~ "Rhinovirus",
    virus %in% c("IAV","IAVH1","IAVH3","IAVH5","IBV","ICV") ~ "Influenza virus*",
    virus %in% c("HBoV","HMPV","Parvo_B19") ~ "Other respiratory viruses*",
    virus %in% c("Aichivirus","Astrovirus","Sapovirus") ~ "Other enteric viruses",
    T ~ "other"
  ),
  virus_group = factor(virus_group,
                       levels = c("Enterovirus*","Hepatitis virus*","Human adenovirus*",
                                  "Human coronavirus*", "Human parainfluenza virus",
                                  "Influenza virus*","Norovirus","Respiratory syncytial virus*",
                                  "Rhinovirus","Rotavirus*","Other enteric viruses",
                                  "Other respiratory viruses*")))

ggplot(med_rt_untargeted_group,aes(season,med_rt,col=virus,group=virus)) +
  geom_hline(yintercept = 0,linetype="dashed") +
  facet_wrap(~virus_group) +
  geom_linerange(aes(ymin = q1_rt, ymax = q3_rt),show.legend = FALSE) +
  geom_line() +
  geom_point() +
  labs(x="",y=expression(paste("Median ", r[t], " (IQR)")),col="Untargeted wastewater\nmonitoring: Boston, MA") +
  scale_x_discrete(limits = c("jun_aug","sep_nov","dec_feb","mar_may"),
                   labels = c("Jun-Aug","Sep-Nov","Dec-Feb","Mar-May")) +
  scale_color_manual(limits = c("Enterovirus_A","Enterovirus_B","Enterovirus_C",
                                "HAV",
                                "HAdVA","HAdVC","HAdVD","HAdVF41",
                                "HCoV229E","HCoVHKU1","HCoVNL63","SC2",
                                "HPIV1","HPIV2","HPIV3","HPIV4",
                                "IAVH1","IAVH3","IBV","ICV",
                                "NoroG1","NoroG2","NoroG4",
                                "RSV_B",
                                "Rhinovirus_A","Rhinovirus_B","Rhinovirus_C",
                                "RotaA","RotaC",
                                "Aichivirus","Astrovirus","Sapovirus",
                                "HBoV","HMPV"),
                     labels = c("EV A","EV B","EV C",
                                "HAV",
                                "HAdV A","HAdV C","HAdV D","HAdV F41",
                                "HCoV 229E","HCoV HKU1","HCoV NL63","SARS-CoV-2",
                                "HPIV 1","HPIV 2","HPIV 3","HPIV 4",
                                "IAV H1N1","IAV H3N1","IBV","ICV",
                                "Noro GI","Noro GII","Noro GIV",
                                "RSV B",
                                "Rhino A","Rhino B","Rhino C",
                                "Rota A","Rota C",
                                "AiV","AstV","Sapo",
                                "HBoV","HMPV"),
                     values = c(color_gradient(3,"Blues"),
                                color_gradient(1,"GnBu"),
                                color_gradient(4,"YlOrBr"),
                                color_gradient(4,"BuPu"),
                                color_gradient(4,"PuRd"),
                                color_gradient(4,"YlGnBu"),
                                color_gradient(3,"YlGn"),
                                color_gradient(1,"Reds"),
                                color_gradient(3,"Greens"),
                                color_gradient(2,"Purples"),
                                color_gradient(3,"RdPu"),
                                color_gradient(2,"PuBuGn")
                                )) +
  theme(
    axis.text.x = element_text(angle = 90,vjust = 0.5),
        # legend.position = "none"
        )

ggsave(here("output/figures/Fig5.png"),width = 7, height = 6.3, units = "in", dpi = 400)

# peak_rt_untarget <- extract_peak_rt(all_rt_estimates) %>% 
#   filter(!(virus %in% c("Norwalk","Enterovirus","Rota","RSV","IAV","HAdVF"))) %>% 
#   mutate(virus = case_when(
#     virus == "Enterovirus_A" ~ "EV A",
#     virus == "Enterovirus_B" ~ "EV B",
#     virus == "Enterovirus_C" ~ "EV C",
#     virus == "HAV" ~ "HAV",
#     virus == "HAdVA" ~ "HAdV A",
#     virus == "HAdVC" ~ "HAdV C",
#     virus == "HAdVD" ~ "HAdV D",
#     virus == "HAdVF41" ~ "HAdV F41",
#     virus == "HCoV229E" ~ "HCoV 229E",
#     virus == "HCoVHKU1" ~ "HCoV HKU1",
#     virus == "HCoVNL63" ~ "HCoV NL63",
#     virus == "SC2" ~ "SARS-CoV-2",
#     virus == "HPIV1" ~ "HPIV 1",
#     virus == "HPIV2" ~ "HPIV 2",
#     virus == "HPIV3" ~ "HPIV 3",
#     virus == "HPIV4" ~ "HPIV 4",
#     virus == "IAVH1" ~ "IAV H1N1",
#     virus == "IAVH3" ~ "IAV H3N1",
#     virus == "IBV" ~ "IBV",
#     virus == "ICV" ~ "ICV",
#     virus == "NoroG1" ~ "Noro GI",
#     virus == "NoroG2" ~ "Noro GII",
#     virus == "NoroG4" ~ "Noro GIV",
#     virus == "RSV_B" ~ "RSV B",
#     virus == "Rhinovirus_A" ~ "Rhino A",
#     virus == "Rhinovirus_B" ~ "Rhino B",
#     virus == "Rhinovirus_C" ~ "Rhino C",
#     virus == "RotaA" ~ "Rota A",
#     virus == "RotaC" ~ "Rota B",
#     virus == "Aichivirus" ~ "AiV",
#     virus == "Astrovirus" ~ "AstV",
#     virus == "Sapovirus" ~ "Sapo",
#     virus == "HBoV" ~ "HBoV",
#     virus == "HMPV" ~ "HMPV",
#     T ~ ""
#   ))
#   
# ggplot(peak_rt_untarget,aes(max_rt_date,reorder(virus, max_rt_date))) +
#   geom_vline(xintercept = c(
#     as.Date("2025-03-01"),
#     as.Date("2025-06-01"),
#     as.Date("2025-09-01"),
#     as.Date("2025-12-01")
#   ), linetype="dashed") +
#   geom_point(size=2) +
#   scale_x_date(breaks = "1 month", date_labels = "%b",limits = c(
#     as.Date("2025-01-01"), as.Date("2025-12-31")
#   )) +
#   labs(x = "2025 peak rₜ timing",y="")
# 
# ggsave(here("output/figures/peak_rt_timing_bos.png"),width = 180, height = 170, units = "mm", dpi = 400) 


# ggplot(peak_rt_untarget,aes(max_rt,reorder(virus, max_rt))) +
#   geom_vline(xintercept = 0,linetype="dashed") +
#   geom_pointrange(aes(xmin=max_rt_min,xmax=max_rt_max),size=0.5,linewidth=0.5) +
#   labs(x = "2025 peak rₜ magnitude (95% CrI)",y="")
# 
# ggsave(here("output/figures/peak_rt_magnitude_bos.png"),width = 180, height = 170, units = "mm", dpi = 400) 

## SI: groups ----
bos_type <- plot_virus_types("Boston_MA")
boise_type <- plot_virus_types("Boise_ID")
ottumwa_type <- plot_virus_types("Ottumwa_IA")
riverside_type <- plot_virus_types("Riverside_CA")

ottumwa_type <- lapply(ottumwa_type, function(x) x &
                         scale_x_date(date_breaks = "4 months",
                                      limits = c(as.Date("2025-02-01"), as.Date("2026-03-31")),
                                      date_labels = "%b %Y"))

# HCoVs
boise_type$HCoV[[4]] <- boise_type$HCoV[[4]] + scale_y_discrete(label="HCoV OC43\n+/-")
bos_type$HCoV[[4]] <- bos_type$HCoV[[4]] + scale_y_discrete(label="HCoV OC43\n+/-")
riverside_type$HCoV[[4]] <- riverside_type$HCoV[[4]] + scale_y_discrete(label="HCoV OC43\n+/-")

((boise_type$HCoV & guides(col = "none", fill = "none")) |
    (bos_type$HCoV & guides(col = "none", fill = "none"))
) /
  (ottumwa_type$HCoV | riverside_type$HCoV & guides(col = "none", fill = "none")) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","","",
                                      "b) Boston, MA","","","","",
                                      "c) Ottumwa, IA","","","","",
                                      "d) Riverside, CA","","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/HCoVs.png"),width = 180, height = 170, units = "mm", dpi = 400)

# HPIVs
ottumwa_type$HPIV[[4]] <- ottumwa_type$HPIV[[4]] + scale_y_discrete(label="HPIV 4\n+/-")

((boise_type$HPIV & guides(col = "none", fill = "none")) |
    (bos_type$HPIV & guides(col = "none", fill = "none"))
) /
  (ottumwa_type$HPIV & guides(col = "none", fill = "none") | riverside_type$HPIV) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","","",
                                      "b) Boston, MA","","","","",
                                      "c) Ottumwa, IA","","","","",
                                      "d) Riverside, CA","","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/HPIVs.png"),width = 180, height = 170, units = "mm", dpi = 400)

# Rota
rota_boise_rt <- boise_type$Rota & guides(col = "none", fill = "none")
rota_boise_rt[[1]] <- rota_boise_rt[[1]] + labs(y="Normalized\nclade couts")
rota_boise_rt[[2]] <- rota_boise_rt[[2]] + labs(y="Rota B\n+/-")
rota_boise_rt[[3]] <- rota_boise_rt[[3]] + labs(y="Normalized\nclade couts")

rota_bos_rt <- bos_type$Rota & guides(col = "none", fill = "none")
rota_bos_rt[[1]] <- rota_bos_rt[[1]] + labs(y="Normalized\nclade couts")
rota_bos_rt[[2]] <- rota_bos_rt[[2]] + labs(y="Rota B\n+/-")
rota_bos_rt[[3]] <- rota_bos_rt[[3]] + labs(y="Normalized\nclade couts")

rota_ottumwa_rt <- ottumwa_type$Rota 

rota_river_rt <- riverside_type$Rota & guides(col = "none", fill = "none")
rota_river_rt[[1]] <- rota_river_rt[[1]] + labs(y="Normalized\nclade couts")
rota_river_rt[[2]] <- rota_river_rt[[2]] + labs(y="Rota B\n+/-")
rota_river_rt[[3]] <- rota_river_rt[[3]] + labs(y="Normalized\nclade couts")

(rota_boise_rt | rota_bos_rt) / (rota_ottumwa_rt | rota_river_rt) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","",
                                      "b) Boston, MA","","","",
                                      "c) Ottumwa, IA","","","",
                                      "d) Riverside, CA","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  ) 
ggsave(here("output/figures/Rota.png"),width = 180, height = 170, units = "mm", dpi = 400)

# EVs
boise_type$EV[[4]] <- boise_type$EV[[4]] + scale_y_discrete(label="EV D68\n+/-")
bos_type$EV[[4]] <- bos_type$EV[[4]] + scale_y_discrete(label="EV D68\n+/-")
ottumwa_type$EV[[4]] <- ottumwa_type$EV[[4]] + scale_y_discrete(label="EV D68\n+/-")

((boise_type$EV & guides(col = "none", fill = "none")) |
    (bos_type$EV & guides(col = "none", fill = "none"))
) /
  (ottumwa_type$EV & guides(col = "none", fill = "none") | riverside_type$EV) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","","",
                                      "b) Boston, MA","","","","",
                                      "c) Ottumwa, IA","","","","",
                                      "d) Riverside, CA","","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/EVs.png"),width = 180, height = 170, units = "mm", dpi = 400)

# RSVs
boise_type$RSV[[1]] <- boise_type$RSV[[1]] + scale_y_discrete(label="RSV A\n+/-")
bos_type$RSV[[1]] <- bos_type$RSV[[1]] + scale_y_discrete(label="RSV A\n+/-")
ottumwa_type$RSV[[1]] <- ottumwa_type$RSV[[1]] + scale_y_discrete(label="RSV A\n+/-")
ottumwa_type$RSV[[2]] <- ottumwa_type$RSV[[2]] + scale_y_discrete(label="RSV B\n+/-")
riverside_type$RSV[[1]] <- riverside_type$RSV[[1]] + scale_y_discrete(label="RSV A\n+/-")
riverside_type$RSV[[2]] <- riverside_type$RSV[[2]] + scale_y_discrete(label="RSV B\n+/-")

((boise_type$RSV)  |
    (bos_type$RSV & guides(col = "none", fill = "none"))
) /
  (ottumwa_type$RSV & guides(col = "none", fill = "none") | riverside_type$RSV) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","",
                                      "b) Boston, MA","","",
                                      "c) Ottumwa, IA","",
                                      "d) Riverside, CA","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/RSVs.png"),width = 180, height = 170, units = "mm", dpi = 400) # change color

# NoVs
((boise_type$NoV)  |
    (bos_type$NoV & guides(col = "none", fill = "none"))
) /
  (ottumwa_type$NoV & guides(col = "none", fill = "none") | riverside_type$NoV) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","",
                                      "b) Boston, MA","","","",
                                      "c) Ottumwa, IA","","","",
                                      "d) Riverside, CA","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/NoVs.png"),width = 180, height = 170, units = "mm", dpi = 400)

# Heps
hep_boise_rt <- boise_type$Hep & guides(col = "none", fill = "none")
hep_boise_rt[[1]] <- hep_boise_rt[[1]] + scale_y_discrete(label="Hep A\n+/-")
hep_boise_rt[[2]] <- hep_boise_rt[[2]] + scale_y_discrete(label="Hep B\n+/-")
hep_boise_rt[[3]] <- hep_boise_rt[[3]] + scale_y_discrete(label="Hep C\n+/-")

hep_bos_rt <- bos_type$Hep & guides(col = "none", fill = "none")
hep_bos_rt[[1]] <- hep_bos_rt[[1]] + labs(y="Normalized\nclade couts")
hep_bos_rt[[2]] <- hep_bos_rt[[2]] + scale_y_discrete(label="Hep B\n+/-")
hep_bos_rt[[3]] <- hep_bos_rt[[3]] + scale_y_discrete(label="Hep C\n+/-")

hep_ottumwa_rt <- ottumwa_type$Hep
hep_ottumwa_rt[[1]] <- hep_ottumwa_rt[[1]] + labs(y="Normalized\nclade couts")
hep_ottumwa_rt[[2]] <- hep_ottumwa_rt[[2]] + labs(y="Normalized\nclade couts")
hep_ottumwa_rt[[3]] <- hep_ottumwa_rt[[3]] + scale_y_discrete(label="Hep C\n+/-")

hep_river_rt <- riverside_type$Hep & guides(col = "none", fill = "none")
hep_river_rt[[1]] <- hep_river_rt[[1]] + labs(y="Normalized\nclade couts")
hep_river_rt[[2]] <- hep_river_rt[[2]] + scale_y_discrete(label="Hep B\n+/-")
hep_river_rt[[3]] <- hep_river_rt[[3]] + scale_y_discrete(label="Hep C\n+/-")

(hep_boise_rt | hep_bos_rt) / (hep_ottumwa_rt | hep_river_rt) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","",
                                      "b) Boston, MA","","","",
                                      "c) Ottumwa, IA","","","",
                                      "d) Riverside, CA","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/Heps.png"),width = 180, height = 170, units = "mm", dpi = 400)

# HAdVs
hadv_boise_rt <- boise_type$HAdV & guides(col = "none", fill = "none")
hadv_boise_rt[[1]] <- hadv_boise_rt[[1]] + labs(y="Normalized\nclade couts")
hadv_boise_rt[[2]] <- hadv_boise_rt[[2]] + scale_y_discrete(label="HAdV B\n+/-")
hadv_boise_rt[[3]] <- hadv_boise_rt[[3]] + labs(y="Normalized\nclade couts")
hadv_boise_rt[[4]] <- hadv_boise_rt[[4]] + labs(y="Normalized\nclade couts")
hadv_boise_rt[[5]] <- hadv_boise_rt[[5]] + scale_y_discrete(label="HAdV F40\n+/-")
hadv_boise_rt[[6]] <- hadv_boise_rt[[6]] + labs(y="Normalized\nclade couts")

hadv_bos_rt <- bos_type$HAdV & guides(col = "none", fill = "none")
hadv_bos_rt[[1]] <- hadv_bos_rt[[1]] + labs(y="Normalized\nclade couts")
hadv_bos_rt[[2]] <- hadv_bos_rt[[2]] + scale_y_discrete(label="HAdV B\n+/-")
hadv_bos_rt[[3]] <- hadv_bos_rt[[3]] + labs(y="Normalized\nclade couts")
hadv_bos_rt[[4]] <- hadv_bos_rt[[4]] + labs(y="Normalized\nclade couts")
hadv_bos_rt[[5]] <- hadv_bos_rt[[5]] + scale_y_discrete(label="HAdV F40\n+/-")
hadv_bos_rt[[6]] <- hadv_bos_rt[[6]] + labs(y="Normalized\nclade couts")

hadv_ottumwa_rt <- ottumwa_type$HAdV & guides(col = "none", fill = "none")
hadv_ottumwa_rt[[1]] <- hadv_ottumwa_rt[[1]] + labs(y="Normalized\nclade couts")
hadv_ottumwa_rt[[2]] <- hadv_ottumwa_rt[[2]] + scale_y_discrete(label="HAdV B\n+/-")
hadv_ottumwa_rt[[3]] <- hadv_ottumwa_rt[[3]] + labs(y="Normalized\nclade couts")
hadv_ottumwa_rt[[4]] <- hadv_ottumwa_rt[[4]] + scale_y_discrete(label="HAdV D\n+/-")
hadv_ottumwa_rt[[5]] <- hadv_ottumwa_rt[[5]] + scale_y_discrete(label="HAdV F40\n+/-")
hadv_ottumwa_rt[[6]] <- hadv_ottumwa_rt[[6]] + labs(y="Normalized\nclade couts")

hadv_river_rt <- riverside_type$HAdV
hadv_river_rt[[1]] <- hadv_river_rt[[1]] + labs(y="Normalized\nclade couts")
hadv_river_rt[[2]] <- hadv_river_rt[[2]] + scale_y_discrete(label="HAdV B\n+/-")
hadv_river_rt[[3]] <- hadv_river_rt[[3]] + labs(y="Normalized\nclade couts")
hadv_river_rt[[4]] <- hadv_river_rt[[4]] + labs(y="Normalized\nclade couts")
hadv_river_rt[[5]] <- hadv_river_rt[[5]] + labs(y="Normalized\nclade couts")
hadv_river_rt[[6]] <- hadv_river_rt[[6]] + labs(y="Normalized\nclade couts")

(hadv_boise_rt | hadv_bos_rt) / (hadv_ottumwa_rt | hadv_river_rt) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","","","","",
                                      "b) Boston, MA","","","","","","",
                                      "c) Ottumwa, IA","","","","","","",
                                      "d) Riverside, CA","","","","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    text = element_text(size=6),
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/HAdVs.png"),width = 180, height = 200, units = "mm", dpi = 400)

# Rhinos
((boise_type$Rhino)  & guides(col = "none", fill = "none") |
    (bos_type$Rhino & guides(col = "none", fill = "none"))
) /
  (ottumwa_type$Rhino  & guides(col = "none", fill = "none") | riverside_type$Rhino ) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","",
                                      "b) Boston, MA","","","",
                                      "c) Ottumwa, IA","","","",
                                      "d) Riverside, CA","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/Rhinos.png"),width = 180, height = 170, units = "mm", dpi = 400)

# Flus
flu_boise_rt <- boise_type$Flu
flu_boise_rt[[1]] <- flu_boise_rt[[1]] + labs(y="Normalized\nclade couts")
flu_boise_rt[[2]] <- flu_boise_rt[[2]] + labs(y="Normalized\nclade couts") 
flu_boise_rt[[3]] <- flu_boise_rt[[3]] + scale_y_discrete(label="IAV H5N1\n+/-")
flu_boise_rt[[4]] <- flu_boise_rt[[4]] + labs(y="Normalized\nclade couts")
flu_boise_rt[[5]] <- flu_boise_rt[[5]] + labs(y="Normalized\nclade couts")

flu_bos_rt <- bos_type$Flu & guides(col = "none", fill = "none")
flu_bos_rt[[1]] <- flu_bos_rt[[1]] + labs(y="Normalized\nclade couts")
flu_bos_rt[[2]] <- flu_bos_rt[[2]] + labs(y="Normalized\nclade couts")
flu_bos_rt[[3]] <- flu_bos_rt[[3]] + scale_y_discrete(label="IAV H5N1\n+/-")
flu_bos_rt[[4]] <- flu_bos_rt[[4]] + labs(y="Normalized\nclade couts")
flu_bos_rt[[5]] <- flu_bos_rt[[5]] + labs(y="Normalized\nclade couts")

flu_ottumwa_rt <- ottumwa_type$Flu & guides(col = "none", fill = "none")
flu_ottumwa_rt[[1]] <- flu_ottumwa_rt[[1]] + labs(y="Normalized\nclade couts")
flu_ottumwa_rt[[2]] <- flu_ottumwa_rt[[2]] + labs(y="Normalized\nclade couts") 
flu_ottumwa_rt[[3]] <- flu_ottumwa_rt[[3]] + scale_y_discrete(label="IAV H5N1\n+/-")
flu_ottumwa_rt[[4]] <- flu_ottumwa_rt[[4]] + scale_y_discrete(label="IBV\n+/-")
flu_ottumwa_rt[[5]] <- flu_ottumwa_rt[[5]] + labs(y="Normalized\nclade couts")

flu_river_rt <- riverside_type$Flu & guides(col = "none", fill = "none")
flu_river_rt[[1]] <- flu_river_rt[[1]] + scale_y_discrete(label="IAV H1N1\n+/-")
flu_river_rt[[2]] <- flu_river_rt[[2]] + labs(y="Normalized\nclade couts")
flu_river_rt[[3]] <- flu_river_rt[[3]] + scale_y_discrete(label="IAV H5N1\n+/-")
flu_river_rt[[4]] <- flu_river_rt[[4]] + labs(y="Normalized\nclade couts")
flu_river_rt[[5]] <- flu_river_rt[[5]] + labs(y="Normalized\nclade couts")

(flu_boise_rt | flu_bos_rt) / (flu_ottumwa_rt | flu_river_rt) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","","","",
                                      "b) Boston, MA","","","","","",
                                      "c) Ottumwa, IA","","","","","",
                                      "d) Riverside, CA","","","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/Flus.png"),width = 180, height = 200, units = "mm", dpi = 400) 

# Other respiratory 
boise_type$OtherResp[[3]] <-boise_type$OtherResp[[3]] + scale_y_discrete(label="Parvo B19\n+/-")
bos_type$OtherResp[[3]] <-bos_type$OtherResp[[3]] + scale_y_discrete(label="Parvo B19\n+/-")
ottumwa_type$OtherResp[[2]] <-ottumwa_type$OtherResp[[2]] + scale_y_discrete(label="HMPV\n+/-")
ottumwa_type$OtherResp[[3]] <-ottumwa_type$OtherResp[[3]] + scale_y_discrete(label="Parvo B19\n+/-")
riverside_type$OtherResp[[3]] <-riverside_type$OtherResp[[3]] + scale_y_discrete(label="Parvo B19\n+/-")

((boise_type$OtherResp)  |
    (bos_type$OtherResp)
) /
  (ottumwa_type$OtherResp  & guides(col = "none", fill = "none") | riverside_type$OtherResp ) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","",
                                      "b) Boston, MA","","","",
                                      "c) Ottumwa, IA","","","",
                                      "d) Riverside, CA","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/OtherResp.png"),width = 180, height = 170, units = "mm", dpi = 400)

# Other enteric
((boise_type$OtherEnteric)  |
    (bos_type$OtherEnteric)
) /
  (ottumwa_type$OtherEnteric  | riverside_type$OtherEnteric ) +
  plot_annotation(tag_levels = list(c("a) Boise, ID","","","",
                                      "b) Boston, MA","","","",
                                      "c) Ottumwa, IA","","","",
                                      "d) Riverside, CA","","",""))) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  labs(x = "") &
  theme(
    plot.tag.location = "panel",  
    plot.tag.position  = "topleft", 
    plot.tag = element_text(hjust = 0.4, vjust = 2, face = "bold"),
    legend.position = "bottom"
  )
ggsave(here("output/figures/OtherEnteric.png"),width = 180, height = 170, units = "mm", dpi = 400)


# ---- Fig 4: Concordance metrics ----
# Schematic
example_df <- extract_Rt_data(Rt_runs$SC2$Boston_MA) %>% 
  filter(source != data_source_labels[["targeted_ww"]],
         x >= as.Date("2024-10-01"))

p <- ggplot(example_df,aes(x,y,col = source)) +
  geom_hline(yintercept = 1,linetype="dashed",linewidth=0.5) +
  geom_line(linewidth=1,alpha=0.6) +
  labs(x="",y=expression(paste("Derived ",R[t])),col="",
       title = "a) Concordance metrics") +
  scale_color_manual(values = c(clinical_col,untargeted_ww_col),
                     labels = c("Data source A","Data source B")) +
  ylim(c(0.875,1.15)) +
  theme(axis.line = element_blank(),
        # axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(face="bold"))

# annotate
x_mae <- as.Date("2024-11-01") # a point where the two curves are visibly apart
y_mae_top <- example_df %>% filter(x==as.Date("2024-11-01") & source == data_source_labels[["clinical"]]) %>% pull(y)
y_mae_bot <- example_df %>% filter(x==as.Date("2024-11-01") & source == data_source_labels[["untargeted_ww"]]) %>% pull(y)

x_aob_start <- as.Date("2025-06-07")  # start of a region where both curves agree
x_aob_end   <-  as.Date("2025-08-14")

x_cor_mid <- as.Date("2025-12-01") # where along the curves to anchor the tangent arrows
cor_dx    <- 21 # how many days ahead each arrow points

get_tangent_arrow <- function(df, source_val, x0, dx) {
  df_s <- df %>% filter(source == source_val)
  i0 <- which.min(abs(as.numeric(df_s$x - x0)))
  i1 <- which.min(abs(as.numeric(df_s$x - (x0 + dx))))
  data.frame(x0 = df_s$x[i0], y0 = df_s$y[i0],
             x1 = df_s$x[i1], y1 = df_s$y[i1])
}

arrow_A   <- get_tangent_arrow(example_df, data_source_labels[["clinical"]], x_cor_mid, cor_dx)
arrow_B <- get_tangent_arrow(example_df, data_source_labels[["untargeted_ww"]], x_cor_mid, cor_dx)

p_annotated <- p +
  
  annotate("segment", x = x_mae, xend = x_mae, y = y_mae_bot, yend = y_mae_top,
           linewidth = 0.5, color = "black", arrow = arrow(ends = "both", length = unit(3, "mm"))) +
  annotate("label", x = x_mae +150, y = mean(c(y_mae_top, y_mae_bot)-0.06),
           label = "MAE: numerical discrepancy", fill = "lightyellow", fontface = "bold",
           size = 3, hjust = 1) +
  
  annotate("rect", xmin = x_aob_start, xmax = x_aob_end, ymin = 1, ymax = 1.09,
           fill = NA, alpha = 0.3, col="black", linewidth = 0.5) +
  annotate("label", x = x_aob_start-28, y = 1.125,
           label = "AOB: directional agreement",
           fill = "lightyellow", fontface = "bold", size = 3,
           hjust = 0) +
  
  geom_segment(data = arrow_A, aes(x = x0, xend = x1, y = y0, yend = y1),
               inherit.aes = FALSE, color = "black", linewidth = 0.5,
               arrow = arrow(ends = "last", length = unit(3, "mm"))) +
  geom_segment(data = arrow_B, aes(x = x0, xend = x1, y = y0, yend = y1),
               inherit.aes = FALSE, color = "black", linewidth = 0.5,
               arrow = arrow(ends = "last", length = unit(3, "mm"))) +
  annotate("label", x = x_cor_mid + 40, y = 0.94,
           label = "Tau: temporal association",
           fill = "lightyellow", fontface = "bold", size = 3) 

# p_annotated 


# Convert output dfs to long format
ern_df_long <- ern_df %>%
  filter(run == "run1" ) %>% # only get run 1
  rename(obs_null_mae_pval = mae_pval,
         obs_null_aob_pval = aob_pval,
         obs_null_cor_pval = cor_pval) %>% 
  pivot_longer(
    cols = matches("^(obs_null|obs|null)_(mae|aob|cor)_(median|lo95|hi95|pval)$"),
    names_pattern = "^(obs_null|obs|null)_(mae|aob|cor)_(median|lo95|hi95|pval)$",
    names_to = c("distribution", "metric", ".value")
  )

rt_df_long <- rt_df %>%
  rename(obs_null_mae_pval = mae_pval,
         obs_null_aob_pval = aob_pval,
         obs_null_cor_pval = cor_pval) %>% 
  pivot_longer(
    cols = matches("^(obs_null|obs|null)_(mae|aob|cor)_(median|lo95|hi95|pval)$"),
    names_pattern = "^(obs_null|obs|null)_(mae|aob|cor)_(median|lo95|hi95|pval)$",
    names_to = c("distribution", "metric", ".value")
  )

ern_df_long <- ern_df_long %>% 
  mutate(metric = case_when(
    metric=="mae" ~ "MAE",
    metric=="aob" ~ "AOB",
    metric=="cor" ~ "Tau"
  ),
  metric = factor(metric,levels=c("MAE","AOB","Tau")),
  virus = ifelse(virus=="NoroG2","NoV",virus),
  virus = factor(virus,levels = c("SC2","IAV","IBV","RSV","NoV","Rota"))
  ) 

rt_df_long <- rt_df_long %>% 
  mutate(metric = case_when(
    metric=="mae" ~ "MAE",
    metric=="aob" ~ "AOB",
    metric=="cor" ~ "Tau"
  ),
  metric = factor(metric,levels=c("MAE","AOB","Tau")),
  virus = ifelse(virus=="NoroG2","NoV",virus))

plt_1 <- ggplot(ern_df_long %>% filter(comparison == "clinical vs untargeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=site,alpha = pval < 0.05),position = position_dodge(width = 0.75),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "b) Clinical vs untargeted wastewater",col="",shape="") +
  scale_y_discrete(limits = rev, drop = FALSE) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Boise, ID","Boston, MA","Ottumwa, IA", "Riverside, CA")) +
  scale_alpha_manual(values = c(0.25,1)) +
  guides(shape = guide_legend(reverse = TRUE)) +
  theme(legend.position = "none",
        plot.title = element_text(face="bold"))

plt_2 <- ggplot(ern_df_long %>% filter(comparison == "targeted_ww vs untargeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=site,alpha = pval < 0.05),position = position_dodge(width = 0.75),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "c) Targeted vs untargeted wastewater",col="",shape="", alpha = "p < 0.05") +
  scale_y_discrete(limits = rev, drop = FALSE) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Boise, ID","Boston, MA","Ottumwa, IA", "Riverside, CA")) +
  scale_alpha_manual(values = c(0.25,1)) +
  guides(shape = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(face="bold"))
  
p_annotated / (plt_1 | plt_2) +
  plot_layout(heights = c(1,2.5),
              # guides = "collect",
              axes = "collect",
              axis_titles = "collect") &
  scale_x_continuous(n.breaks = 4)
ggsave(here("output/figures/Fig4.png"),width = 180, height = 170, units = "mm", dpi = 400)

## SI: clinical vs targeted ----
ggplot(ern_df_long %>% filter(comparison == "clinical vs targeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=site,alpha = pval < 0.05),position = position_dodge(width = 0.75),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "Clinical vs targeted wastewater",col="",shape="",alpha="p < 0.05") +
  scale_y_discrete(limits = rev) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Boise, ID","Boston, MA","Ottumwa, IA", "Riverside, CA")) +
  scale_alpha_manual(values = c(0.25,1)) +
  guides(shape = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(face="bold"))
ggsave(here("output/figures/Rt_clinical_targeted_concordance.png"),width = 150, height = 170, units = "mm", dpi = 400)

## SI: rt ----
plt_1_rt <- ggplot(rt_df_long %>% filter(comparison == "clinical vs untargeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=site,alpha = pval < 0.05),position = position_dodge(width = 0.5),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "a) Clinical vs untargeted wastewater",col="",shape="",alpha="p < 0.05") +
  scale_y_discrete(limits = rev) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Boise, ID","Boston, MA","Ottumwa, IA", "Riverside, CA")) +
  scale_alpha_manual(values = c(0.25,1)) +
  guides(shape = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(face="bold"))

plt_2_rt <- ggplot(rt_df_long %>% filter(comparison == "targeted_ww vs untargeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=site,alpha = pval < 0.05),position = position_dodge(width = 0.5),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "b) Targeted vs untargeted wastewater",col="",shape="", alpha = "p < 0.05") +
  scale_y_discrete(limits = rev) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Boise, ID","Boston, MA","Ottumwa, IA", "Riverside, CA")) +
  scale_alpha_manual(values = c(0.25,1), guide="none") +
  guides(shape = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(face="bold"))

plt_3_rt <- ggplot(rt_df_long %>% filter(comparison == "clinical vs targeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=site,alpha = pval < 0.05),position = position_dodge(width = 0.5),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "c) Clinical vs targeted wastewater",col="",shape="",alpha="p < 0.05") +
  scale_y_discrete(limits = rev) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Boise, ID","Boston, MA","Ottumwa, IA", "Riverside, CA")) +
  scale_alpha_manual(values = c(0.25,1)) +
  guides(shape = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(face="bold"))

(plt_1_rt | plt_2_rt) +
  plot_layout(guides = "collect") &
  scale_x_continuous(n.breaks = 4)
ggsave(here("output/figures/rt_concordance.png"),width = 180, height = 170, units = "mm", dpi = 400)

## SI: Shedding ----
ern_bos_df_long <- ern_df %>%
  filter(site == "Boston_MA" & virus %in% c("SC2","IAV","RSV")) %>% # only get boston, viruses with alt shedding
  rename(obs_null_mae_pval = mae_pval,
         obs_null_aob_pval = aob_pval,
         obs_null_cor_pval = cor_pval) %>% 
  pivot_longer(
    cols = matches("^(obs_null|obs|null)_(mae|aob|cor)_(median|lo95|hi95|pval)$"),
    names_pattern = "^(obs_null|obs|null)_(mae|aob|cor)_(median|lo95|hi95|pval)$",
    names_to = c("distribution", "metric", ".value")
  )

ern_bos_df_long <- ern_bos_df_long %>% 
  mutate(metric = case_when(
    metric=="mae" ~ "MAE",
    metric=="aob" ~ "AOB",
    metric=="cor" ~ "Tau"
  ),
  metric = factor(metric,levels=c("MAE","AOB","Tau")),
  virus = ifelse(virus=="NoroG2","NoV",virus),
  virus = factor(virus,levels = c("SC2","IAV","IBV","RSV","NoV","Rota"))
  )

plt_1_Rt_Bos <- ggplot(ern_bos_df_long %>% filter(comparison == "clinical vs untargeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=run,alpha = pval < 0.05),position = position_dodge(width = 0.5),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "a) Clinical vs untargeted wastewater",col="",shape="",alpha="p < 0.05") +
  scale_y_discrete(limits = rev) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Main","Alt A","Alt B", "Alt C")) +
  scale_alpha_manual(values = c(1), guide="none") +
  # guides(shape = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(face="bold"))

plt_2_Rt_Bos <- ggplot(ern_bos_df_long %>% filter(comparison == "targeted_ww vs untargeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=run,alpha = pval < 0.05),position = position_dodge(width = 0.5),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "b) Targeted vs untargeted wastewater",col="",shape="", alpha = "p < 0.05") +
  scale_y_discrete(limits = rev) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Main","Alt A","Alt B", "Alt C")) +
  scale_alpha_manual(values = c(0.25,1)) +
  # guides(shape = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(face="bold"))

plt_3_Rt_Bos <- ggplot(ern_bos_df_long %>% filter(comparison == "clinical vs targeted_ww", distribution=="obs_null")) +
  geom_vline(xintercept = 0,linetype="dashed",linewidth=0.25) +
  facet_wrap(~metric,scales="free_x") +
  geom_pointrange(aes(x=median,y=virus,xmin=lo95, xmax=hi95,col=metric,shape=run,alpha = pval < 0.05),position = position_dodge(width = 0.5),size=0.5,linewidth = 0.5) +
  labs(x="Median observed - null (95% CrI)",y="",title = "c) Clinical vs targeted wastewater",col="",shape="",alpha="p < 0.05") +
  scale_y_discrete(limits = rev) +
  scale_color_viridis_d(option = "E", begin = 0.2,end = 0.8, guide="none") +
  scale_shape_manual(values = c(15,16,17,18), labels = c("Main","Alt A","Alt B", "Alt C")) +
  scale_alpha_manual(values = c(1), guide="none") +
  # guides(shape = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(face="bold"))

(plt_1_Rt_Bos | plt_2_Rt_Bos) +
  plot_layout(guides = "collect") &
  scale_x_continuous(n.breaks = 4)
ggsave(here("output/figures/Rt_concordance_Bos_shedding.png"),width = 180, height = 170, units = "mm", dpi = 400)
