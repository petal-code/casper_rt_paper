# All estimation functions: Rt (EpiNow2) and rt (GAM smooth), plus their
# shared data-prep helpers (outlier identification, time series fill/scale)


#' Flag outliers in a wastewater time series using a modified z-score: a point is an outlier
#' if (log10(`column`) - median) / MAD -- computed across all detects (`is_nd == 0`) in `data`
#' -- exceeds 3.5, the conventional modified z-score threshold
#' @param data A df containing at least `column`, `date_column`, and `is_nd`
#' @param column Name of the numeric column to check for outliers (default "confirm")
#' @param date_column Name of the date column (default "date"); used only to sort `data`
#'   before flagging
#'
#' @return `data` with an added logical column `is_outlier`
identify_ww_outlier <- function(data, column = "confirm", date_column = "date") {
  col_sym <- sym(column)
  
  data <- data %>% arrange(!!sym(date_column))  
  
  log_vals <- log10(data[[column]])
  dates <- data[[date_column]]
  
  data %>%
    mutate(
      med = median(log_vals[is_nd==0]),
      mad = mad(log_vals[is_nd==0]), # the scaling factor is accounted for with mad()
      mod_z = (log_vals - med) / mad,
      is_outlier = mod_z > 3.5
    ) %>% 
    select(-med, -mad, -mod_z)
}

#' Estimate growth rate (rt): fit a GAM smooth to a single data stream and derive the
#' instantaneous growth rate as the slope of the log-smooth (chain rule: f'(t)/f(t)).
#' @param df A df of data (date, confirm, and -- for wastewater streams -- is_nd)
#' @param family GAM family (e.g. nb() for clinical counts, tw() for ww concentrations); a
#'   Tweedie family triggers the wastewater-specific non-detect skip and outlier filtering
#'
#' @return A list with:
#'   \item{gam_fit}{The fitted GAM object, or NULL if no GAM was fit}
#'   \item{gam_pred}{Smoothed time series (median and 90% CrI for each date), or NULL}
#'   \item{rt_pred}{Growth rate rt (median and 90% CrI for each date), or NULL}
#'   \item{rt_draws}{All 1000 per-draw rt estimates before summarising to median and CrI, or NULL}
#'   \item{k}{Basis dimension used for the GAM smooth, or NULL}
fit_smooths <- function(df, family) {
  
  family_str <- family$family
  
  is_ww <- family_str=="Tweedie"
  
  if(is_ww && mean(df$is_nd==1)>0.75){
    list(gam_fit = NULL, gam_pred = NULL, rt_pred = NULL, rt_draws = NULL, k = NULL)
  }
  else{
    
    if(is_ww){
      df <- identify_ww_outlier(df)
      df <- df %>% filter(is.na(is_outlier) | !is_outlier)
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
        rt = c(NA, diff(log(.fitted))) # rt ~ f'(t)/f(t) -- chain rule!
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

#' Fill and scale clinical, targeted ww, and untargeted ww time series for EpiNow2. Each ww
#' source is scaled by the nearest power of 10 that brings its non-zero median in line with
#' clinical data's non-zero median (or, absent clinical data, untargeted ww is scaled to
#' targeted ww instead); if no reference source is available, or the result still leaves the
#' median below 1, the source is scaled up further until its own non-zero median clears 1. Any
#' columns besides date and confirm (e.g. is_nd) are dropped before scaling/filling
#' @param dfs_list A named list with elements clinical, targeted_ww, and/or untargeted_ww,
#'   each a df with columns (date, confirm); NULL if not available
#'
#' @return A list with filled dfs (clinical, targeted_ww, untargeted_ww), the multiplicative
#'   scaling factors applied to each ww source (mult_factors), and the proportion of ww
#'   observations rounded to zero after scaling (prop_zero)
fill_time_series <- function(dfs_list) {

  # Strip extra columns — EpiNow2 only accepts date and confirm
  dfs_list <- lapply(dfs_list, function(df) {
    if (is.data.frame(df)) select(df, date, confirm) else df
  })

  # Compute nonzero median of a df -- use nonzero median in case NDs kept as 0s
  nonzero_median <- function(df) {
    median(df$confirm[df$confirm > 0], na.rm = TRUE)
  }

  # Compute mult factor (nearest power of 10) to match source to reference
  compute_mult_factor <- function(reference_df, source_df) {
    if (is.null(reference_df) || is.null(source_df)) return(1)
    ref_median    <- nonzero_median(reference_df)
    source_median <- nonzero_median(source_df)
    if (is.na(source_median) || source_median == 0) return(1)
    if (is.na(ref_median)    || ref_median    <= 0) return(above_one_factor(source_df))
    raw_ratio <- ref_median / source_median
    10 ^ round(log10(raw_ratio))
  }

  # Compute the minimum power-of-10 factor needed to bring median above 1
  above_one_factor <- function(df) {
    if (is.null(df)) return(1)
    m <- nonzero_median(df)
    if (is.na(m) || m <= 0) {
      warning("Non-zero median is NA or 0 — cannot scale, returning factor 1")
      return(1)
    }
    if (m >= 1) return(1)
    10 ^ ceiling(-log10(m))
  }

  # Calculate proportion rounded to zero
  prop_rounded_to_zero <- function(df, f) {
    if (is.null(df)) return(NA)
    mean(round(df$confirm * f) == 0, na.rm = TRUE)
  }

  # Initialize so return() can always reference them
  mult_factors <- list(targeted_ww = 1, untargeted_ww = 1)
  zero_prop_targeted   <- NA
  zero_prop_untargeted <- NA

  # Fill / scale dfs
  if (!is.null(dfs_list$clinical)) {
    dfs_list$clinical <- fill_missing(
      data = dfs_list$clinical,
      missing_dates = "accumulate", # assuming weekly totals; no need to scale bc already integer counts
      initial_accumulate = 7
    )
  }

  if (!is.null(dfs_list$targeted_ww)) {
    f <- if (!is.null(dfs_list$clinical)) { # if clinical data present
      compute_mult_factor(dfs_list$clinical, dfs_list$targeted_ww) # scale to clinical data
    } else { # if no clinical data present
      above_one_factor(dfs_list$targeted_ww) # scale so median above 1 at least
    }
    mult_factors$targeted_ww <- f
    zero_prop_targeted <- prop_rounded_to_zero(dfs_list$targeted_ww, f)
    dfs_list$targeted_ww <- fill_missing(
      data = dfs_list$targeted_ww %>% mutate(confirm = round(confirm * f)),
      missing_dates = "ignore" # dates represent distinct samples
    )
  }

  if (!is.null(dfs_list$untargeted_ww)) {
    f <- if (!is.null(dfs_list$clinical)) { # if clinical data present
      compute_mult_factor(dfs_list$clinical, dfs_list$untargeted_ww) # scale to clinical data
    } else { # if no clinical data present
      compute_mult_factor(dfs_list$targeted_ww, dfs_list$untargeted_ww) # scale to targeted ww if available (compute_mult_factor returns 1 if reference is NULL)
    }
    # Ensure non-zero median is above 1 after scaling
    f <- f * above_one_factor(data.frame(confirm = dfs_list$untargeted_ww$confirm * f))
    mult_factors$untargeted_ww <- f
    zero_prop_untargeted <- prop_rounded_to_zero(dfs_list$untargeted_ww, f)
    dfs_list$untargeted_ww <- fill_missing(
      data = dfs_list$untargeted_ww %>% mutate(confirm = round(confirm * f)),
      missing_dates = "ignore" # dates represent distinct samples
    )
  }

  return(list(
    clinical          = dfs_list$clinical,
    targeted_ww       = dfs_list$targeted_ww,
    untargeted_ww     = dfs_list$untargeted_ww,
    mult_factors      = mult_factors,
    prop_zero         = list(
      targeted_ww   = zero_prop_targeted,
      untargeted_ww = zero_prop_untargeted
    )
  ))
}

#' Parse one row of an EpiNow2 run-settings sheet into the priors list expected by estimate_Rt()
#' @param settings_row A single row from one sheet of all_virus_settings (one sheet per virus, one row per run)
#'
#' @return A named list of parsed priors, passed to estimate_Rt()
parse_run_settings <- function(settings_row) {
  eval_str <- function(x) if (!is.na(x) && nzchar(x)) eval(parse(text = x)) else NULL
  list(
    generation_time   = eval_str(settings_row[["generation_time"]]),
    incubation_period = eval_str(settings_row[["incubation period"]]),
    reporting_delay   = eval_str(settings_row[["reporting delay"]]),
    shedding          = eval_str(settings_row[["shedding"]]),
    initial_Rt        = eval_str(settings_row[["initial_Rt"]]),
    obs_scale_case    = eval_str(settings_row[["obs_scale_cases"]]),
    obs_scale_ww      = eval_str(settings_row[["obs_scale_cases_ww"]]),
    week_effect       = eval_str(settings_row[["week_effect"]]),
    dispersion        = eval_str(settings_row[["dispersion"]]),
    ls                = eval_str(settings_row[["ls"]]),
    alpha             = eval_str(settings_row[["alpha"]]),
    truncation        = eval_str(settings_row[["truncation"]])
  )
}

#' Estimate Rt for clinical, targeted ww, and/or untargeted ww data via EpiNow2::epinow(). Each
#' available source is fit independently and with a different delay distribution: clinical
#' data uses an incubation period + reporting delay (run_epinow()), while both wastewater
#' sources use a shedding-profile delay instead (run_epinow_ww()); a source whose df argument
#' is NULL is skipped (returned as NULL rather than fit)
#' @param clinical_df A df of clinical data formatted for EpiNow2 (date, confirm), filled; NULL if not available
#' @param targeted_ww_df A df of targeted ww data formatted for EpiNow2 (date, confirm), filled with mult factor; NULL if not available
#' @param untargeted_ww_df A df of untargeted ww data formatted for EpiNow2 (date, confirm), filled with mult factor; NULL if not available
#' @param virus_priors A parsed settings list from parse_run_settings()
#'
#' @return A named list (clinical, targeted_ww, untargeted_ww) of epinow() output for each
#'   source; NULL for any source whose df argument was NULL
estimate_Rt <- function(clinical_df, targeted_ww_df, untargeted_ww_df, virus_priors) {

  run_epinow <- function(df, obs_scale) {
    epinow(
      data            = df,
      generation_time = gt_opts(virus_priors$generation_time),
      delays          = delay_opts(virus_priors$incubation_period + virus_priors$reporting_delay),
      rt              = rt_opts(prior = virus_priors$initial_Rt),
      gp              = gp_opts(ls = virus_priors$ls, alpha = virus_priors$alpha),
      obs             = obs_opts(
        scale       = obs_scale,
        dispersion  = virus_priors$dispersion,
        week_effect = virus_priors$week_effect
      ),
      truncation = trunc_opts(virus_priors$truncation),
      stan       = stan_opts(cores = 1, chains = 4, seed = seed)
    )
  }

  run_epinow_ww <- function(df, obs_scale) {
    epinow(
      data            = df,
      generation_time = gt_opts(virus_priors$generation_time),
      delays          = delay_opts(dist = virus_priors$shedding),
      rt              = rt_opts(prior = virus_priors$initial_Rt),
      gp              = gp_opts(ls = virus_priors$ls, alpha = virus_priors$alpha),
      obs             = obs_opts(
        scale       = obs_scale,
        dispersion  = virus_priors$dispersion,
        week_effect = virus_priors$week_effect
      ),
      truncation = trunc_opts(virus_priors$truncation),
      stan       = stan_opts(cores = 1, chains = 4, seed = seed)
    )
  }

  list(
    clinical      = if (!is.null(clinical_df))      run_epinow(clinical_df, virus_priors$obs_scale_case)       else NULL,
    targeted_ww   = if (!is.null(targeted_ww_df))   run_epinow_ww(targeted_ww_df, virus_priors$obs_scale_ww)   else NULL,
    untargeted_ww = if (!is.null(untargeted_ww_df)) run_epinow_ww(untargeted_ww_df, virus_priors$obs_scale_ww) else NULL
  )
}
