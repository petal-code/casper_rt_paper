# Functions to compare posterior estimates between data sources


#' Extract a variable's posterior draws from get_samples() output
#' @param posterior_samples Output of EpiNow2::get_samples()
#' @param var Variable name (default "R" for Rt)
#'
#' @return tibble with date, sample, value
extract_Rt_draws <- function(posterior_samples, var = "R") {
  posterior_samples %>%
    filter(variable == var) %>%
    select(date, sample, value)
}

#' Reformat rt_draws from fit_smooths() into (date, sample, value) format
#' for use with compare_posteriors()
#' @param rt_draws rt_draws element from fit_smooths() output
#'
#' @return tibble(date, sample, value) where value is the per-draw rt estimate
format_rt_draws <- function(rt_draws) {
  rt_draws %>%
    select(date = x, sample = .draw, value = rt) 
}


# ----- Per-draw metric calculations -----
# Each function returns a numeric vector of length J (one value per posterior draw),
# constituting the posterior distribution of that metric.

#' MAE posterior distribution: per-draw mean(|Rt_a - Rt_b|) over shared dates
#' Draws are matched by (date, sample index).
#' @param draws_a tibble(date, sample, value) from extract_Rt_draws() or format_rt_draws()
#' @param draws_b tibble(date, sample, value) from extract_Rt_draws() or format_rt_draws()
#'
#' @return Numeric vector of length J (one MAE value per posterior draw)
mae_posterior <- function(draws_a, draws_b) {
  inner_join(draws_a,
             draws_b,
             by = c("date", "sample"),
             suffix = c("_a", "_b")) %>%
    group_by(sample) %>%
    summarise(mae = mean(abs(value_a - value_b)), .groups = "drop") %>%
    pull(mae)
}

#' AOB posterior distribution: per-draw proportion of dates where both draws agree relative
#' to threshold -- both above it, both below it, or both exactly equal to it (returned as a
#' proportion in [0, 1])
#' @param draws_a tibble(date, sample, value) from extract_Rt_draws() or format_rt_draws()
#' @param draws_b tibble(date, sample, value) from extract_Rt_draws() or format_rt_draws()
#' @param threshold Threshold value for agreement (use 1 for Rt; use 0 for rt)
#'
#' @return Numeric vector of length J (one AOB proportion per posterior draw)
aob_posterior <- function(draws_a, draws_b, threshold) {
  inner_join(draws_a,
             draws_b,
             by = c("date", "sample"),
             suffix = c("_a", "_b")) %>%
    mutate(agree = (value_a > threshold & value_b > threshold) |
             (value_a < threshold & value_b < threshold) |
             (value_a == threshold & value_b == threshold)) %>%
    group_by(sample) %>%
    summarise(aob = mean(agree), .groups = "drop") %>%
    pull(aob)
}

#' Correlation posterior distribution: per-draw Kendall's tau between Rt_a and Rt_b
#' @param draws_a tibble(date, sample, value) from extract_Rt_draws() or format_rt_draws()
#' @param draws_b tibble(date, sample, value) from extract_Rt_draws() or format_rt_draws()
#'
#' @return Numeric vector of length J (one Kendall's tau value per posterior draw)
cor_posterior <- function(draws_a, draws_b) {
  inner_join(draws_a,
             draws_b,
             by = c("date", "sample"),
             suffix = c("_a", "_b")) %>%
    group_by(sample) %>%
    summarise(tau = cor(value_a, value_b, method = "kendall"),
              .groups = "drop") %>%
    pull(tau)
}


# ----- Permutation tests -----
# Within each draw, shuffle draws_b values across time points to break temporal
# dependence while preserving the marginal distribution.

#' Permute draws_b's values across dates within each draw, to build the null comparison for
#' compare_posteriors()'s permutation tests
#' @param draws_b tibble(date, sample, value) to permute
#'
#' @return draws_b with values shuffled across dates within each draw, preserving marginal distribution
permute_draws <- function(draws_b) {
  draws_b %>%
    group_by(sample) %>%
    mutate(value = sample(value, replace = F)) %>%
    ungroup()
}

#' Bayesian permutation p-value for MAE
#' p = proportion of draws where delta = observed - null >= 0
#' (i.e., null MAE <= observed MAE — observed did NOT improve over null)
#' Significant concordance: p < 0.05
#' @param obs_mae Numeric vector of observed MAE draws
#' @param null_mae Numeric vector of null (permuted) MAE draws
#'
#' @return Scalar p-value in [0, 1]
pval_mae <- function(obs_mae, null_mae) {
  delta <- obs_mae - null_mae
  mean(delta >= 0)
}

#' Bayesian permutation p-value for AOB
#' p = proportion of draws where delta = observed - null <= 0
#' (i.e., observed AOB <= null — observed did NOT improve over null)
#' Significant concordance: p < 0.05
#' @param obs_aob Numeric vector of observed AOB draws
#' @param null_aob Numeric vector of null (permuted) AOB draws
#'
#' @return Scalar p-value in [0, 1]
pval_aob <- function(obs_aob, null_aob) {
  delta <- obs_aob - null_aob
  mean(delta <= 0)
}

#' Bayesian permutation p-value for correlation
#' p = proportion of draws where delta = observed - null <= 0
#' (i.e., observed correlation <= null — observed did NOT improve over null)
#' Significant concordance: p < 0.05
#' @param obs_cor Numeric vector of observed Kendall's tau draws
#' @param null_cor Numeric vector of null (permuted) Kendall's tau draws
#'
#' @return Scalar p-value in [0, 1]
pval_cor <- function(obs_cor, null_cor) {
  delta <- obs_cor - null_cor
  mean(delta <= 0)
}


# ----- Summary wrapper -----
#' Summarize a posterior metric vector as median and 95% CrI
#' @param x Numeric vector of posterior draws for a metric
#'
#' @return Named list with median, lo95, and hi95
summarise_posterior <- function(x) {
  list(
    median = median(x, na.rm = TRUE),
    lo95   = quantile(x, 0.025, na.rm = TRUE),
    hi95   = quantile(x, 0.975, na.rm = TRUE)
  )
}

#' Compute MAE, AOB, and Kendall correlation between two posterior draw tibbles,
#' with permutation-based p-values. The null distribution for each metric comes from
#' permute_draws()-shuffling draws_b's values across dates within each draw, then recomputing
#' the metric against the unshuffled draws_a -- breaking any real date-to-date relationship
#' between the two models while preserving each model's own marginal distribution
#'
#' @param draws_a tibble(date, sample, value) from extract_Rt_draws() or format_rt_draws(), model A
#' @param draws_b tibble(date, sample, value) from extract_Rt_draws() or format_rt_draws(), model B
#' @param label_a Label for model A
#' @param label_b Label for model B
#' @param aob_threshold Threshold for AOB agreement -- required, no default (use 1 for Rt, 0 for rt)
#'
#' @return A named list with:
#'   \item{obs_mae_dist}{Numeric vector (length J): observed per-draw MAE, from mae_posterior()}
#'   \item{obs_aob_dist}{Numeric vector (length J): observed per-draw AOB, from aob_posterior()}
#'   \item{obs_cor_dist}{Numeric vector (length J): observed per-draw Kendall's tau, from cor_posterior()}
#'   \item{perm_draws_b}{draws_b with its values permuted within each draw, from permute_draws()}
#'   \item{null_mae_dist}{Numeric vector (length J): null per-draw MAE (draws_a vs perm_draws_b)}
#'   \item{null_aob_dist}{Numeric vector (length J): null per-draw AOB (draws_a vs perm_draws_b)}
#'   \item{null_cor_dist}{Numeric vector (length J): null per-draw Kendall's tau (draws_a vs perm_draws_b)}
#'   \item{summary_stats}{One-row tibble: comparison label, plus median/95% CrI for the
#'     observed, null, and observed-minus-null distributions above, and a permutation p-value,
#'     for each of MAE/AOB/correlation -- see summarise_posterior() / pval_mae() / pval_aob() /
#'     pval_cor()}
compare_posteriors <- function(draws_a,
                               draws_b,
                               label_a = "A",
                               label_b = "B",
                               aob_threshold) {
  obs_mae  <- mae_posterior(draws_a, draws_b)
  obs_aob <- aob_posterior(draws_a, draws_b, threshold = aob_threshold)
  obs_cor  <- cor_posterior(draws_a, draws_b)
  
  draws_b_perm <- permute_draws(draws_b)
  
  null_mae  <- mae_posterior(draws_a, draws_b_perm)
  null_aob <- aob_posterior(draws_a, draws_b_perm, threshold = aob_threshold)
  null_cor  <- cor_posterior(draws_a, draws_b_perm)
  
  obs_mae_s  <- summarise_posterior(obs_mae)
  obs_aob_s <- summarise_posterior(obs_aob)
  obs_cor_s  <- summarise_posterior(obs_cor)
  
  null_mae_s  <- summarise_posterior(null_mae)
  null_aob_s <- summarise_posterior(null_aob)
  null_cor_s  <- summarise_posterior(null_cor)
  
  obs_null_mae_s  <- summarise_posterior(obs_mae - null_mae)
  obs_null_aob_s <- summarise_posterior(obs_aob - null_aob)
  obs_null_cor_s  <- summarise_posterior(obs_cor - null_cor)
  
  summary_tbl <- tibble(
    comparison     = paste(label_a, "vs", label_b),
    obs_mae_median     = obs_mae_s$median,
    obs_mae_lo95       = obs_mae_s$lo95,
    obs_mae_hi95       = obs_mae_s$hi95,
    null_mae_median     = null_mae_s$median,
    null_mae_lo95       = null_mae_s$lo95,
    null_mae_hi95       = null_mae_s$hi95,
    obs_null_mae_median     = obs_null_mae_s$median,
    obs_null_mae_lo95       = obs_null_mae_s$lo95,
    obs_null_mae_hi95       = obs_null_mae_s$hi95,
    mae_pval       = pval_mae(obs_mae, null_mae),
    obs_aob_median     = obs_aob_s$median,
    obs_aob_lo95       = obs_aob_s$lo95,
    obs_aob_hi95       = obs_aob_s$hi95,
    null_aob_median     = null_aob_s$median,
    null_aob_lo95       = null_aob_s$lo95,
    null_aob_hi95       = null_aob_s$hi95,
    obs_null_aob_median     = obs_null_aob_s$median,
    obs_null_aob_lo95       = obs_null_aob_s$lo95,
    obs_null_aob_hi95       = obs_null_aob_s$hi95,
    aob_pval       = pval_aob(obs_aob, null_aob),
    obs_cor_median     = obs_cor_s$median,
    obs_cor_lo95       = obs_cor_s$lo95,
    obs_cor_hi95       = obs_cor_s$hi95,
    null_cor_median     = null_cor_s$median,
    null_cor_lo95       = null_cor_s$lo95,
    null_cor_hi95       = null_cor_s$hi95,
    obs_null_cor_median     = obs_null_cor_s$median,
    obs_null_cor_lo95       = obs_null_cor_s$lo95,
    obs_null_cor_hi95       = obs_null_cor_s$hi95,
    cor_pval       = pval_cor(obs_cor, null_cor),
  )
  
  return(list(obs_mae_dist = obs_mae,
              obs_aob_dist = obs_aob,
              obs_cor_dist = obs_cor,
              perm_draws_b = draws_b_perm,
              null_mae_dist = null_mae,
              null_aob_dist = null_aob,
              null_cor_dist = null_cor,
              summary_stats = summary_tbl))
  
}
