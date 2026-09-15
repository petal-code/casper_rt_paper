# Estimate Rt and rt and assess concordance

source(here::here("code/00_config.R"))

plan(multisession, workers = availableCores() - 1, # update as needed based on computer
     rscript_args = "--no-init-file",
     rscript_libs = .libPaths())


# ----- Load data -----
all_virus_dfs <- read_rds(here(paste0("data/clean/subsetted/grouped/all_virus_dfs", nd_suffix, ".rds")))


# ---- Temporal associations ----
#' Compute Kendall's tau + pval, plus each source's  Shapiro-Wilk normality pval
#' between the weekly-aggregated observations of every combination of data sources at a site.
#' @param site_dfs Named list of (date, confirm, and -- for wastewater sources -- is_nd) dfs per data source 
#' (as in all_virus_dfs[[virus]][[site]]); NULL/absent entries are skipped
#'
#' @return A tibble with one row per remaining source pair: comparison ("<source A> vs
#'   <source B>"), tau, p_value, n (number of shared weeks the pair was computed over), and
#'   shapiro_p_a/shapiro_p_b (each source's own Shapiro-Wilk p-value on its weekly series).
#'   tau/p_value/shapiro_p_* are NA where there are too few weeks (fewer than 3) to compute
#'   them. Zero rows if fewer than 2 sources remain after the LOD exclusion.
compute_site_tau <- function(site_dfs) {
  site_dfs <- Filter(Negate(is.null), site_dfs)

  # Drop wastewater sources that are mostly below the limit of detection
  site_dfs <- Filter(function(df) {
    !("is_nd" %in% names(df)) || mean(df$is_nd == 1, na.rm = TRUE) <= 0.75
  }, site_dfs)

  sources <- names(site_dfs)
  if (length(sources) < 2) return(tibble())

  # Aggregate every source to MMWR week (Sunday-Saturday) -- sampling days can differ
  # across sources within the same week, so match on week rather than raw date
  weekly_dfs <- lapply(site_dfs, function(df) {
    df %>%
      mutate(week = floor_date(date, unit = "week", week_start = "Sunday") + 6) %>%
      group_by(week) %>%
      summarise(confirm = mean(confirm, na.rm = TRUE), .groups = "drop")
  })

  shapiro_p <- function(df) {
    if (nrow(df) < 3 || nrow(df) > 5000) return(NA_real_)
    tryCatch(
      suppressWarnings(shapiro.test(df$confirm)$p.value),
      error = function(e) NA_real_
    )
  }
  shapiro_ps <- vapply(weekly_dfs, shapiro_p, numeric(1))

  pairs <- combn(sources, 2, simplify = FALSE)
  map_dfr(pairs, function(pair) {
    joined <- inner_join(
      weekly_dfs[[pair[1]]],
      weekly_dfs[[pair[2]]],
      by = "week", suffix = c("_a", "_b")
    )
    test <- if (nrow(joined) >= 3) {
      suppressWarnings(cor.test(joined$confirm_a, joined$confirm_b, method = "kendall"))
    } else {
      list(estimate = NA_real_, p.value = NA_real_)
    }
    tibble(
      comparison  = paste(pair[1], "vs", pair[2]),
      tau         = unname(test$estimate),
      p_value     = unname(test$p.value),
      n           = nrow(joined),
      shapiro_p_a = shapiro_ps[[pair[1]]],
      shapiro_p_b = shapiro_ps[[pair[2]]]
    )
  })
}

# For each virus/site, calculate tau between each pair of data sources and save tau and p val
temporal_associations <- imap_dfr(all_virus_dfs, function(virus_dfs, virus_name) {
  imap_dfr(virus_dfs, function(site_dfs, site_name) {
    compute_site_tau(site_dfs) %>% mutate(site = site_name)
  }) %>% mutate(virus = virus_name)
}) %>%
  select(virus, site, comparison, tau, p_value, n, shapiro_p_a, shapiro_p_b)

write_csv(temporal_associations, here(paste0("output/input_time_series/temporal_associations", nd_suffix, ".csv")))

# temporal_associations %>% 
#   group_by(comparison) %>%
#   summarise(med_tau = median(tau),
#             q1_tau = quantile(tau,0.25),
#             q3_tau = quantile(tau,0.75))


# ----- Estimations -----
## rt ----
run_rt <- FALSE  # TRUE = run estimation and save RDS; FALSE = load saved RDS

rt_rds_path <- here("output/time_varying_estimates/", paste0("all_rt", nd_suffix, ".rds"))

if (run_rt) {
  all_rt_estimates <- map(all_virus_dfs, function(virus_dfs) {
    future_lapply(virus_dfs, function(site_dfs) {
        list(
          clinical      = if (!is.null(site_dfs$clinical))      fit_smooths(site_dfs$clinical,      family = nb()),
          targeted_ww   = if (!is.null(site_dfs$targeted_ww))   fit_smooths(site_dfs$targeted_ww,   family = tw()),
          untargeted_ww = if (!is.null(site_dfs$untargeted_ww)) fit_smooths(site_dfs$untargeted_ww, family = tw())
        )
    }, future.seed = TRUE)
  }) 
  saveRDS(all_rt_estimates, rt_rds_path)
} else {
  all_rt_estimates <- read_rds(rt_rds_path)
}

## Rt ----
run_Rt <- FALSE  # TRUE = run estimation and save RDS; FALSE = load saved RDS

# Fill dfs and apply multiplication factors to ww data
all_dfs_fill <- lapply(all_virus_dfs[c("SC2","IAV","IBV","RSV","NoroG2","Rota")], 
                       function(virus_dfs) { 
                         lapply(virus_dfs, fill_time_series)
                         })

# Map short virus names to Excel settings sheet names
virus_settings_key <- c(SC2 = "sars-cov-2", IAV = "iav", IBV = "ibv", RSV = "rsv", NoroG2 = "norovirus", Rota = "rotavirus")

if(impute_nondetects==FALSE){
  all_virus_settings <- lapply(all_virus_settings,function(x) slice(x,1)) # only execute for run 1 settings
}

# Execute
all_Rt_runs <- imap(virus_settings_key, function(sheet, virus_name) {
  settings_df <- all_virus_settings[[sheet]]
  runs <- lapply(seq_len(nrow(settings_df)), function(run_i) {
    rds_path <- here("output/time_varying_estimates/", paste0(virus_name, "_run", run_i, nd_suffix, "_Rt.rds"))
    if (run_Rt) {
      priors    <- parse_run_settings(settings_df[run_i, ])
      estimates <- future_lapply(all_dfs_fill[[virus_name]], function(site_dfs) {
        estimate_Rt(
          clinical_df      = site_dfs$clinical,
          targeted_ww_df   = site_dfs$targeted_ww,
          untargeted_ww_df = site_dfs$untargeted_ww,
          virus_priors     = priors
        )
      }, future.seed = TRUE)
      saveRDS(estimates, rds_path)
      estimates
    } else {
      read_rds(rds_path)
    }
  })
  names(runs) <- paste0("run", seq_len(length(runs)))
  runs
})

all_Rt_runs$IBV$run1$Ottumwa_IA$untargeted_ww <- NULL # IBV/Boise/untargeted data is mostly NDs --> don't actually use Rt estimates below


# ----- Plot estimates -----
## rt ----
iwalk(all_rt_estimates, function(virus_estimates, virus_name) {
  iwalk(virus_estimates, function(site_smooths, site_name) {
    rt_plot <- plot_time_varying_estimates(
      dfs           = all_virus_dfs[[virus_name]][[site_name]],
      smooths       = lapply(site_smooths, function(s) if (!is.null(s)) s$gam_pred),
      estimate_data = extract_rt_data(site_smooths),
      estimate_ylab = "rt",
      yintercept    = 0
    )
    ggsave(
      filename = here("output/output_time_series/rt_growth_rate", paste0(virus_name, "_rt_", site_name, nd_suffix, ".png")),
      plot     = rt_plot,
      dpi = 300, width = 6, height = 4, units = "in"
    )
  })
})

## Rt ----
iwalk(all_Rt_runs, function(virus_runs, virus_name) {
  iwalk(virus_runs, function(run_estimates, run_name) {
    iwalk(run_estimates, function(site_estimates, site_name) {
      Rt_plot <- plot_time_varying_estimates(
        dfs           = all_virus_dfs[[virus_name]][[site_name]],
        estimate_data = extract_Rt_data(site_estimates),
        estimate_ylab = "Rt",
        yintercept    = 1,
        yrange        = NULL
      )
      ggsave(
        filename = here("output/output_time_series/Rt_ERN", paste0(virus_name, "_", run_name, "_ts_Rt_", site_name, nd_suffix, ".png")),
        plot     = Rt_plot,
        dpi = 300, width = 6, height = 4, units = "in"
      )
    })
  })
})


# ----- Assess concordance -----
#' Run compare_posteriors() on every pairwise combination of data sources at a site
#' @param site_draws Named list of tibble(date, sample, value) draws, one per data source
#'   (as produced by extract_Rt_draws() or format_rt_draws())
#' @param aob_threshold Threshold for AOB agreement (1 for Rt, 0 for rt); passed through
#'   to compare_posteriors()
#'
#' @return Named list of compare_posteriors() outputs, one per source pair, named
#'   "<source A>_v_<source B>"
compare_all_pairs <- function(site_draws, aob_threshold) {
  sources <- names(site_draws)
  if (length(sources) < 2) return(list())
  pairs   <- combn(sources, 2, simplify = FALSE)
  results <- lapply(pairs, function(pair) {
    compare_posteriors(
      site_draws[[pair[1]]],
      site_draws[[pair[2]]],
      label_a        = pair[1],
      label_b        = pair[2],
      aob_threshold = aob_threshold
    )
  })
  names(results) <- sapply(pairs, function(p) paste(p[1], "v", p[2], sep = "_"))
  results
}

#' Flatten compare_all_pairs() output across all sites for one virus into a single tidy
#' tibble of concordance summary statistics
#' @param comparisons Named list of compare_all_pairs() outputs, one per site
#' @param virus Virus name to record in the output (e.g. "SC2")
#'
#' @return A tibble with columns virus, site, comparison, and the MAE/AOB/correlation
#'   summary statistics from compare_posteriors()
extract_concordance_table <- function(comparisons, virus) {
  imap_dfr(comparisons, function(site_pairs, site_name) {
    map_dfr(site_pairs, function(pair) pair$summary_stats) %>%
      mutate(site = site_name)
  }) %>%
    mutate(virus = virus) %>%
    select(virus, site, everything())
}

#' Build and save both concordance figures (observed distributions, and observed - null
#' distributions) for every site in a compare_all_pairs() result
#' @param comparisons Named list of compare_all_pairs() outputs, one per site (each containing
#'   up to the three source-pair comparisons: clinical_v_targeted_ww, clinical_v_untargeted_ww,
#'   targeted_ww_v_untargeted_ww)
#' @param prefix Filename prefix (e.g. "SC2_rt")
#' @param folder Output folder, relative to the project root
#'
#' @return NULL; saves one obs plot and one obs-null plot per site
save_concordance_plots <- function(comparisons, prefix, folder) {
  comparisons <- Filter(function(site) length(site) > 0, comparisons)

  obs_plots <- lapply(comparisons, function(site) {
    do.call(plot_obs_concordance, list(
      clinical_v_targeted   = site[["clinical_v_targeted_ww"]],
      clinical_v_untargeted = site[["clinical_v_untargeted_ww"]],
      targeted_v_untargeted = site[["targeted_ww_v_untargeted_ww"]]
    ))
  })
  iwalk(obs_plots, function(site_plot, site_name) {
    ggsave(
      filename = here(folder, paste0(prefix, "_concordance_obs_", site_name, nd_suffix, ".png")),
      plot = site_plot, dpi = 300, width = 13, height = 9, units = "in"
    )
  })

  diff_plots <- lapply(comparisons, function(site) {
    do.call(plot_obs_null_concordance, list(
      clinical_v_targeted   = site[["clinical_v_targeted_ww"]],
      clinical_v_untargeted = site[["clinical_v_untargeted_ww"]],
      targeted_v_untargeted = site[["targeted_ww_v_untargeted_ww"]]
    ))
  })
  iwalk(diff_plots, function(site_plot, site_name) {
    ggsave(
      filename = here(folder, paste0(prefix, "_concordance_diff_", site_name, nd_suffix, ".png")),
      plot = site_plot, dpi = 300, width = 13, height = 9, units = "in"
    )
  })
}

## rt ----
rt_concordance_table <- imap_dfr(all_rt_estimates, function(virus_estimates, virus_name) {
  rt_draws <- lapply(virus_estimates, function(site) {
    site_valid <- Filter(function(stream) !is.null(stream) && !is.null(stream$rt_draws), site)
    lapply(site_valid, function(stream) {
      format_rt_draws(stream$rt_draws)
    })
  })
  comparisons <- lapply(rt_draws, compare_all_pairs, aob_threshold = 0)
  save_concordance_plots(comparisons, prefix = paste0(virus_name, "_rt"), folder = "output/output_concordance/rt_growth_rate")
  extract_concordance_table(comparisons, virus = virus_name)
})
write_csv(rt_concordance_table, here(paste0("output/output_concordance/rt_growth_rate/rt_concordance_table", nd_suffix, ".csv")))

## Rt ----
Rt_concordance_table <- imap_dfr(all_Rt_runs, function(virus_runs, virus_name) {
  imap_dfr(virus_runs, function(run_estimates, run_name) {
    post_samples <- lapply(run_estimates, function(site) {
      lapply(Filter(Negate(is.null), site), get_samples)
    })
    post_samples_Rt <- lapply(post_samples, function(site) {
      lapply(site, function(source) extract_Rt_draws(source, var = "R"))
    })
    comparisons <- lapply(post_samples_Rt, compare_all_pairs, aob_threshold = 1)
    save_concordance_plots(comparisons, prefix = paste0(virus_name, "_",run_name,"_Rt"), folder = "output/output_concordance/Rt_ERN")
    extract_concordance_table(comparisons, virus = virus_name) %>%
      mutate(run = run_name, .after = virus)
  })
})
write_csv(Rt_concordance_table, here(paste0("output/output_concordance/Rt_ERN/Rt_concordance_table", nd_suffix, ".csv")))

plan(sequential)
