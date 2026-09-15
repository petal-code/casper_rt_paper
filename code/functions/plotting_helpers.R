# All plotting helper functions

# Display labels for the three data sources
data_source_labels <- c(clinical      = "Clinical surveillance",
                        targeted_ww   = "Targeted wastewater monitoring",
                        untargeted_ww = "Untargeted wastewater monitoring")

#' Plot a heatmap of detect/non-detect by date
#' @param detections A df with columns date and is_nd
#' @param col Fill color for detections
#' @param ylab Y-axis title
#' 
#' @return A ggplot object
plot_heatmap_panel <- function(detections, col, ylab){
  ggplot(detections) +
    geom_tile(aes(date, ylab, fill=is_nd), width=1) +
    labs(x = "", y="") +
    scale_fill_manual(values = c(`FALSE`=col,`TRUE`="gray90"), guide = "none")
}

#' Plot a single raw time series; falls back to plot_heatmap_panel() when mostly non-detect (for ww)
#' @param data_points A df with columns date, confirm, and (if is_ww) is_nd
#' @param col Point/line color
#' @param data_source Label used for the color legend
#' @param is_ww TRUE for wastewater sources; FALSE for clinical data
#' @param ylab Y-axis title
#'
#' @return A ggplot object
plot_ts_panel <- function(data_points, col, data_source, is_ww, ylab, plot_out = TRUE){
  if(is_ww && mean(data_points$is_nd) > 0.75){
    return(plot_heatmap_panel(data_points, col, data_source))
  }

  if(!is_ww){
    ggplot(data_points,aes(x = date, y = confirm, col = data_source)) +
      geom_point(size = 0.5) +
      geom_line(linewidth = 0.5) +
      labs(x = "", y=ylab, col="") +
      scale_y_continuous(n.breaks = 4) +
      scale_color_manual(values = c(col), guide = "none") +
      ylim(c(0,NA)) +
      theme(legend.position = "bottom")
  }
  else{
    data_points <- identify_ww_outlier(data_points)
    if(!plot_out) data_points <- data_points %>% filter(is.na(is_outlier) | !is_outlier)
    ggplot(data_points,aes(x = date, y = confirm, col = data_source)) +
      geom_point(aes(shape = is_nd), size = 0.5) +
      geom_line(linewidth = 0.5) +
      labs(x = "", y=ylab, col="") +
      scale_y_continuous(n.breaks = 4) +
      scale_color_manual(values = c(col), guide = "none") +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
      ylim(c(0,NA)) +
      theme(legend.position = "bottom")
  }
}

#' Plot a single GAM-smoothed time series; falls back to plot_heatmap_panel() when mostly non-detect (for ww)
#' @param data_points A df with columns date, confirm, and (if is_ww) is_nd
#' @param data_smooth gam_pred output from fit_smooths() (df with x, y, ymin, ymax)
#' @param col Point/line/ribbon color
#' @param data_source Label used for the color legend
#' @param is_ww TRUE for wastewater sources; FALSE for clinical data
#' @param ylab Y-axis title
#' @param plot_out TRUE (default) keeps ww outliers in the plot, shown transparent; FALSE filters
#'   them out before plotting so they aren't shown at all. No effect when !is_ww
#'
#' @return A ggplot object
plot_smooth_ts_panel <- function(data_points, data_smooth, col, data_source, is_ww, ylab, plot_out = TRUE){
  if(is_ww && mean(data_points$is_nd) > 0.75){
    return(plot_heatmap_panel(data_points, col, data_source))
  }

  smooth_layers <- if (!is.null(data_smooth)) {
    list(
      geom_ribbon(data = data_smooth, aes(x = x, y = y, ymin = ymin, ymax = ymax, fill = data_source), alpha = 0.2),
      geom_line(data = data_smooth, aes(x = x, y = y, col = data_source), size = 0.5)
    )
  }

  if(!is_ww){
    ggplot() +
      smooth_layers +
      geom_point(data = data_points, aes(x = date, y = confirm, col = data_source), size = 0.5) +
      labs(x = "", y=ylab, col="",fill="") +
      scale_y_continuous(n.breaks = 4) +
      scale_color_manual(values = c(col), guide = "none") +
      scale_fill_manual(values = c(col), guide = "none") +
      theme(legend.position = "bottom")
  }
  else{
    data_points <- identify_ww_outlier(data_points)
    if(!plot_out) data_points <- data_points %>% filter(is.na(is_outlier) | !is_outlier)
    ggplot() +
      smooth_layers +
      geom_point(data = data_points, aes(x = date, y = confirm, col = data_source, shape = is_nd, alpha = is_outlier), size = 0.5) +
      labs(x = "", y=ylab, col="",fill="") +
      scale_y_continuous(n.breaks = 4) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
      scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.25), na.value = 1, guide = "none") +
      scale_color_manual(values = c(col), guide = "none") +
      scale_fill_manual(values = c(col), guide = "none") +
      theme(legend.position = "bottom")
  }
}

#' Combine per-source fit_smooths() output into a long df for plot_time_varying_estimates
#' @param smooths Named list (clinical, targeted_ww, untargeted_ww) of fit_smooths() output;
#'   missing/NULL entries are skipped
#'
#' @return tibble(x, y, ymin, ymax, source)
extract_rt_data <- function(smooths) {
  bind_rows(lapply(names(data_source_labels), function(source_name) {
    rt_pred <- smooths[[source_name]]$rt_pred
    if (is.null(rt_pred)) return(NULL)
    rt_pred %>% mutate(source = data_source_labels[[source_name]])
  }))
}

#' Extract Rt summaries from estimate_Rt() output into a long df for plot_time_varying_estimates
#' @param estimates Named list (clinical, targeted_ww, untargeted_ww) of epinow() output;
#'   missing/NULL entries are skipped
#'
#' @return tibble(x, y, ymin, ymax, source)
extract_Rt_data <- function(estimates) {
  bind_rows(lapply(names(data_source_labels), function(source_name) {
    est <- estimates[[source_name]]
    if (is.null(est)) return(NULL)
    summary(est, type = "parameters", params = "R") %>%
      transmute(x = date, y = median, ymin = lower_90, ymax = upper_90, source = data_source_labels[[source_name]])
  }))
}

#' Build the ribbon + line estimate panel for plot_time_varying_estimates() (one line/ribbon per
#' group -- e.g. per data source, or per subtype within a virus group)
#' @param estimate_data Long tibble (x, y, ymin, ymax, `id_col`) of the time-varying estimate --
#'   see extract_rt_data() / extract_Rt_data(). No rows (e.g. no source produced an estimate)
#'   returns NULL rather than building a ggplot from an empty, columnless data frame
#' @param id_col Name of the column in estimate_data that distinguishes lines/ribbons (e.g.
#'   "source" when combining clinical/targeted/untargeted ww for one virus, or "type" when
#'   combining subtypes within a virus group)
#' @param cols Named vector of colors, names matching the values in estimate_data[[id_col]]
#' @param ylab Y-axis label (e.g. "rt" or "Rt")
#' @param yintercept Reference line (0 for rt, 1 for Rt)
#' @param yrange y-axis limits; NULL for automatic
#'
#' @return A ggplot object, or NULL if estimate_data has no rows
build_estimate_panel <- function(estimate_data, id_col, cols, ylab, yintercept, yrange) {
  if (nrow(estimate_data) == 0) return(NULL)

  p <- estimate_data %>%
    ggplot(aes(x = x, col = .data[[id_col]], fill = .data[[id_col]])) +
    geom_hline(yintercept = yintercept,
               linetype = "dashed",
               linewidth = 0.5) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax),
                alpha = 0.2, col = NA) +
    geom_line(aes(y = y), linewidth = 0.5) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(
      x="",
      y = ylab,
      col = "",
      fill = ""
    ) +
    theme(legend.position = "bottom")

  if (!is.null(yrange)) p <- p + ylim(yrange)
  p
}

#' Build one raw/smoothed time series panel per source, for plot_time_varying_estimates()
#' @param dfs Named list of raw (date, confirm) dfs, one per source/stream, names matching
#'   `labels`; NULL entries are skipped
#' @param smooths Named list of gam_pred dfs (x, y, ymin, ymax) to overlay, one per source/stream
#' @param smooth_ts TRUE overlays the GAM smooth on each raw ts (plot_smooth_ts_panel); FALSE
#'   plots the raw ts alone (plot_ts_panel)
#' @param labels Named character vector mapping `dfs`/`smooths` names to display labels
#' @param cols Named vector of colors, names matching the values in `labels`
#' @param is_ww Named logical vector (same names as `labels`) -- TRUE for wastewater streams
#'   (adds non-detect shape, may fall back to heatmap); FALSE for clinical data
#' @param ylab Y-axis title for the panels
#' @param plot_out TRUE keeps ww outliers in the panels, shown transparent; FALSE filters them
#'   out before plotting. Only applies when `smooth_ts` is TRUE; no effect on raw/non-ww panels
#'
#' @return A list of ggplot objects, one per source with non-NULL data in `dfs` (sources with no
#'   data are dropped, not returned as NULL)
build_ts_panels <- function(dfs, smooths, smooth_ts, labels, cols, is_ww, ylab, plot_out) {
  ts_plots <- lapply(names(labels), function(source_name) {
    if (is.null(dfs[[source_name]])) return(NULL)
    if (smooth_ts) {
      plot_smooth_ts_panel(dfs[[source_name]], smooths[[source_name]], cols[[labels[[source_name]]]], labels[[source_name]], is_ww[[source_name]], ylab, plot_out)
    } else {
      plot_ts_panel(dfs[[source_name]], cols[[labels[[source_name]]]], labels[[source_name]], is_ww[[source_name]], ylab, plot_out)
    }
  })
  Filter(Negate(is.null), ts_plots)
}

#' Build a generic time-varying-estimate panel (ribbon + line per group; rt, Rt, or similar),
#' optionally stacked under a time series panel per data source
#' @param dfs Named list of raw (date, confirm) dfs, one per source/stream -- see
#'   build_ts_panels(). Unused when `include_ts` is FALSE
#' @param smooths Named list of gam_pred dfs to overlay on the raw ts panels -- see
#'   build_ts_panels(). Supplying this argument at all (regardless of individual entries being
#'   NULL) switches the ts panel style from raw (plot_ts_panel) to GAM-smoothed
#'   (plot_smooth_ts_panel)
#' @param estimate_data Long tibble (x, y, ymin, ymax, `id_col`) of the time-varying estimate --
#'   see build_estimate_panel() / extract_rt_data() / extract_Rt_data()
#' @param labels Named character vector mapping `dfs`/`smooths` names to display labels (default
#'   `data_source_labels`)
#' @param cols Named vector of colors, names matching the values in `labels` -- see
#'   build_estimate_panel() / build_ts_panels()
#' @param is_ww Named logical vector (same names as `labels`) -- see build_ts_panels()
#' @param id_col Name of the column in estimate_data that distinguishes lines/ribbons -- see
#'   build_estimate_panel()
#' @param ts_ylab Y-axis title for the time series panels -- see build_ts_panels()
#' @param estimate_ylab Y-axis label for the estimate panel -- see build_estimate_panel()
#' @param yintercept Reference line for the estimate panel -- see build_estimate_panel()
#' @param yrange y-axis limits for the estimate panel -- see build_estimate_panel()
#' @param include_ts TRUE (default) to stack a time series panel per source above the estimate
#'   panel; FALSE to build and return just the estimate panel on its own
#' @param plot_out TRUE (default) keeps ww outliers in the ts panels, shown transparent; FALSE
#'   filters them out before plotting -- see build_ts_panels()
#' @param estimate_height_ratio How many times taller the estimate panel is than each ts panel
#'   (default 3); ignored when the estimate panel is omitted (no rows in estimate_data)
#' @param title Optional title, added via labs(title=) to the topmost panel (the first ts panel,
#'   or the estimate panel itself when include_ts is FALSE) rather than plot_annotation(). A
#'   plot_annotation() title is dropped once this composite gets nested inside another patchwork
#'   composition (e.g. combining several virus panels into a grid), so it has to live on a leaf
#'   ggplot instead -- same reasoning as the tag workaround on plot_virus_panel()
#'
#' @return If `include_ts`, a patchwork plot: one ts panel per available source, stacked above
#'   the estimate panel (the estimate panel is `estimate_height_ratio` times taller than each ts
#'   panel). Otherwise, just the estimate panel (a ggplot object), or NULL if estimate_data has no
#'   rows
plot_time_varying_estimates <- function(dfs = list(), smooths = list(), estimate_data,
                                         labels = data_source_labels,
                                         cols = setNames(c(clinical_col, targeted_ww_col, untargeted_ww_col), data_source_labels),
                                         is_ww = c(clinical = FALSE, targeted_ww = TRUE, untargeted_ww = TRUE),
                                         id_col = "source",
                                         ts_ylab = "Virus signal", estimate_ylab, yintercept, yrange = NULL,
                                         include_ts = TRUE, plot_out = TRUE, estimate_height_ratio = 2, title = NULL) {

  estimate_panel <- build_estimate_panel(estimate_data, id_col, cols, estimate_ylab, yintercept, yrange)

  if (!include_ts) {
    if (!is.null(title) && !is.null(estimate_panel)) estimate_panel <- estimate_panel + labs(title = title)
    return(estimate_panel)
  }

  ts_plots <- build_ts_panels(dfs, smooths, !missing(smooths), labels, cols, is_ww, ts_ylab, plot_out)

  if (!is.null(title) && length(ts_plots) > 0) ts_plots[[1]] <- ts_plots[[1]] + labs(title = title)

  x_limits <- range(do.call(c, lapply(Filter(Negate(is.null), dfs), `[[`, "date")))

  all_plots <- Filter(Negate(is.null), c(ts_plots, list(estimate_panel)))

  # Estimate panel (if present) is stacked last -- give it more vertical room than the
  # ts panels above it, which each get equal weight
  panel_heights <- if (!is.null(estimate_panel)) {
    c(rep(1, length(ts_plots)), estimate_height_ratio)
  } else {
    rep(1, length(ts_plots))
  }

  # guides is deliberately left unset (not "collect") here: only the estimate panel ever carries
  # a legend (ts panels set guide = "none"), so there's nothing to collect locally, and setting
  # "collect" at this level would resolve/seal the guide -- making it opaque to a guides =
  # "collect" on an *outer* plot_layout() if this composite later gets nested inside a bigger one
  # (e.g. a grid of several virus panels sharing one legend)
  Reduce(`/`, all_plots) +
    plot_layout(axes = "collect",
                axis_titles = "collect",
                heights = panel_heights) &
    labs(x = "") &
    scale_x_date(date_labels = "%b %Y", limits = x_limits) &
    theme(legend.position = "bottom")
}

#' Build a sequential palette of `n` colorblind-safe shades from a single RColorBrewer
#' sequential palette, for distinguishing levels within one group (e.g.
#' subtypes within a virus) while keeping them all drawn from the same, single-hue map. Always
#' pulls the palette's full 9-color version and drops the palest 2 shades (near-white, illegible
#' against a white background) before interpolating to `n`, so this works for any n >= 1
#' regardless of RColorBrewer::brewer.pal()'s own 3-color minimum
#' @param n Number of shades to generate (e.g. number of levels/subtypes present)
#' @param palette Name of a colorblind-safe RColorBrewer sequential palette (e.g. "Blues",
#'   "Oranges", "Greens" -- see subset(RColorBrewer::brewer.pal.info, category == "seq" &
#'   colorblind) for the full colorblind-safe list)
#'
#' @return Character vector of `n` hex colors, light to dark
color_gradient <- function(n, palette) {
  base <- RColorBrewer::brewer.pal(9, palette)[3:9]
  colorRampPalette(base)(n)
}

#' Build a stacked panel of smoothed untargeted-ww time series + a combined rt panel
#' for a group of related virus subtypes
#' @param site Site name -- key into all_virus_dfs[[virus_key]] / all_rt_estimates[[virus_key]]
#' @param virus_labels Named character vector: names are virus_key values (keys into
#'   all_virus_dfs / all_rt_estimates), values are display labels 
#' @param palette Name of a colorblind-safe RColorBrewer sequential palette (see color_gradient())
#'   this group's subtypes are shaded from -- gives each virus group its own distinct palette
#' @param ylab Y-axis title for the time series panels
#' @param rt_range y-axis limits for the rt panel
#' @param plot_out TRUE keeps ww outliers in the ts panels, shown transparent; FALSE (default)
#'   filters them out before plotting so they aren't shown at all
#' @param title Optional title, passed through to plot_time_varying_estimates() -- put on a
#'   leaf panel via labs(title=) rather than plot_annotation(), so it (and the panel's axes)
#'   survive being nested inside an outer patchwork composition without needing wrap_elements()
#'
#' @return A patchwork plot: one time series panel per virus in virus_labels, stacked above
#'   a combined rt panel (omitted if no virus in the group produced an rt estimate)
plot_virus_group <- function(site, virus_labels, palette, ylab = "Normalized clade counts", rt_range = NULL, plot_out = FALSE, title = NULL){

  colors <- setNames(color_gradient(length(virus_labels), palette), unname(virus_labels))
  is_ww  <- setNames(rep(TRUE, length(virus_labels)), names(virus_labels))

  # data is always present; only the GAM smooth and rt estimate can be missing
  # (e.g. a stream with >75% non-detects never gets a GAM fit)
  dfs     <- lapply(names(virus_labels), function(key) all_virus_dfs[[key]][[site]]$untargeted_ww)
  smooths <- lapply(names(virus_labels), function(key) all_rt_estimates[[key]][[site]]$untargeted_ww$gam_pred)
  names(dfs) <- names(smooths) <- names(virus_labels)

  rt_list <- lapply(names(virus_labels), function(key) {
    virus_rt <- all_rt_estimates[[key]][[site]]$untargeted_ww$rt_pred
    if (!is.null(virus_rt)) virus_rt %>% mutate(type = virus_labels[[key]])
  })
  rt_data <- bind_rows(rt_list)

  if (all(vapply(dfs, is.null, logical(1))) && nrow(rt_data) == 0) return(NULL)

  plot_time_varying_estimates(
    dfs, smooths, rt_data,
    labels = virus_labels, cols = colors, is_ww = is_ww, id_col = "type",
    ts_ylab = ylab, estimate_ylab = expression(r[t]), yintercept = 0, yrange = rt_range,
    plot_out = plot_out, title = title
  )
}

#' Build stacked smoothed time series + rt panels comparing subtypes within each virus group,
#' for all untargeted-ww virus groups at a given site. Wraps plot_virus_group() once per
#' predefined virus group (HCoV, EV, Rota, HAdV, Hep, RSV, Rhino, NoV, Flu, HPIV)
#' @param site Site name -- key into all_virus_dfs[[virus_key]] / all_rt_estimates[[virus_key]]
#' @param plot_out TRUE keeps ww outliers in the ts panels, shown transparent; FALSE (default)
#'   filters them out before plotting so they aren't shown at all
#' @param title Optional title, passed through to every plot_virus_group() call (e.g. a
#'   site label like "a) Boise, ID") -- see plot_virus_group()'s own `title` docs
#'
#' @return A named list of patchwork plots (or NULL), one per virus group
plot_virus_types <- function(site, plot_out = FALSE, title = NULL){

  # One colorblind-safe RColorBrewer sequential palette per virus group -- plot_virus_group()
  # shades each into a gradient sized to however many subtypes that particular group has
  group_palettes <- c(HCoV = "BuPu", EV = "Blues", Rota = "Purples", HAdV = "YlOrBr",
                       Hep = "GnBu", RSV = "Reds", Rhino = "Greens", NoV = "YlGn",
                       Flu = "YlGnBu", HPIV = "PuRd",OtherResp ="PuBuGn",OtherEnteric="RdPu")

  HCoV_labels  <- c(HCoV229E = "HCoV 229E", HCoVHKU1 = "HCoV HKU1", HCoVNL63 = "HCoV NL63", HCoVOC43 = "HCoV OC43")
  EV_labels    <- c(Enterovirus_A = "EV A", Enterovirus_B = "EV B", Enterovirus_C = "EV C", EVD68 = "EV D68")
  Rota_labels  <- c(RotaA = "Rotavirus A", RotaB = "Rotavirus B",RotaC = "Rotavirus C")
  HAdV_labels  <- c(HAdVA = "HAdV A", HAdVB = "HAdV B", HAdVC = "HAdV C", HAdVD = "HAdV D",
                     HAdVF40 = "HAdV F40", HAdVF41 = "HAdV F41")
  Hep_labels   <- c(HAV = "Hep A", HBV = "Hep B", HCV = "Hep C")
  RSV_labels   <- c(RSV_A = "RSV A", RSV_B = "RSV B")
  Rhino_labels <- c(Rhinovirus_A = "Rhino A", Rhinovirus_B = "Rhino B", Rhinovirus_C = "Rhino C")
  NoV_labels   <- c(NoroG1 = "NoV I", NoroG2 = "NoV II", NoroG4 = "NoV IV")
  Flu_labels   <- c(IAVH1 = "IAV H1N1", IAVH3 = "IAV H3N1", IAVH5 = "IAV H5N1", IBV = "IBV", ICV = "ICV")
  HPIV_labels  <- c(HPIV1 = "HPIV 1", HPIV2 = "HPIV 2", HPIV3 = "HPIV 3", HPIV4 = "HPIV 4")
  OtherEnteric_labels  <- c(Aichivirus="AiV",Astrovirus="AstV",Sapovirus="Sapo")
  OtherResp_labels  <- c(HBoV = "HBoV",HMPV="HMPV",Parvo_B19="Parvo B19")

  return(list(
    HCoV  = plot_virus_group(site, HCoV_labels, group_palettes[["HCoV"]], plot_out = plot_out, title = title),
    EV    = plot_virus_group(site, EV_labels, group_palettes[["EV"]], plot_out = plot_out, title = title),
    Rota  = plot_virus_group(site, Rota_labels, group_palettes[["Rota"]], plot_out = plot_out, title = title),
    HAdV  = plot_virus_group(site, HAdV_labels, group_palettes[["HAdV"]], plot_out = plot_out, title = title),
    Hep   = plot_virus_group(site, Hep_labels, group_palettes[["Hep"]], plot_out = plot_out, title = title),
    RSV   = plot_virus_group(site, RSV_labels, group_palettes[["RSV"]], plot_out = plot_out, title = title),
    Rhino = plot_virus_group(site, Rhino_labels, group_palettes[["Rhino"]], plot_out = plot_out, title = title),
    NoV   = plot_virus_group(site, NoV_labels, group_palettes[["NoV"]], plot_out = plot_out, title = title),
    Flu   = plot_virus_group(site, Flu_labels, group_palettes[["Flu"]], plot_out = plot_out, title = title),
    HPIV  = plot_virus_group(site, HPIV_labels, group_palettes[["HPIV"]], plot_out = plot_out, title = title),
    OtherEnteric  = plot_virus_group(site, OtherEnteric_labels, group_palettes[["OtherEnteric"]], plot_out = plot_out, title = title),
    OtherResp  = plot_virus_group(site, OtherResp_labels, group_palettes[["OtherResp"]], plot_out = plot_out, title = title)
  ))
}

#' Plot observed posterior distribution for each concordance metric
#' @param clinical_v_targeted Output from compare_posteriors() for clinical vs targeted ww data
#' @param clinical_v_untargeted Output from compare_posteriors() for clinical vs untargeted ww data
#' @param targeted_v_untargeted Output from compare_posteriors() for targeted vs untargeted ww data
#' 
#' @return A patchwork of all concordance metrics (observed distributions)
plot_obs_concordance <- function(clinical_v_targeted,
                                 clinical_v_untargeted,
                                 targeted_v_untargeted) {
  # ---- MAE ----
  ## Clinical vs targeted ----
  if (is.null(clinical_v_targeted)) {
    mae_clinical_targeted <- NULL
  }
  else{
    mae_clinical_targeted <- ggplot() +
      geom_histogram(aes(clinical_v_targeted$obs_mae_dist),
                     fill = mae_col,
                     col = "black") +
      labs(x = "MAE",
           title = "Clinical vs Targeted WW",
           fill = "")
  }
  
  ## Clinical vs untargeted ----
  if (is.null(clinical_v_untargeted)) {
    mae_clinical_untargeted <- NULL
  }
  else{
    mae_clinical_untargeted <- ggplot() +
      geom_histogram(
        aes(clinical_v_untargeted$obs_mae_dist),
        fill = mae_col,
        col = "black"
      ) +
      labs(x = "MAE",
           title = "Clinical vs Untargeted WW",
           fill = "")
    
  }
  
  ## Targeted vs untargeted ----
  if (is.null(targeted_v_untargeted)) {
    mae_targeted_untargeted <- NULL
  }
  else{
    mae_targeted_untargeted <- ggplot() +
      geom_histogram(
        aes(targeted_v_untargeted$obs_mae_dist),
        fill = mae_col,
        col = "black"
      ) +
      labs(x = "MAE",
           title = "Targeted vs Unargeted WW",
           fill = "")
    
  }
  
  ## Combine ----
  mae_plots <- list(mae_clinical_targeted,
                    mae_clinical_untargeted,
                    mae_targeted_untargeted)
  mae_plots <- Filter(Negate(is.null), mae_plots)
  mae_plot <- Reduce(`+`, mae_plots) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "bottom")
  
  # ---- AOB ----
  ## Clinical vs targeted ----
  if (is.null(clinical_v_targeted)) {
    aob1_clinical_targeted <- NULL
  }
  else{
    aob1_clinical_targeted <- ggplot() +
      geom_histogram(
        aes(clinical_v_targeted$obs_aob_dist),
        fill = aob_col,
        col = "black"
      ) +
      labs(x = "AOB",
           title = "Clinical vs Targeted WW",
           fill = "")
  }
  
  ## Clinical vs untargeted ----
  if (is.null(clinical_v_untargeted)) {
    aob1_clinical_untargeted <- NULL
  }
  else{
    aob1_clinical_untargeted <- ggplot() +
      geom_histogram(
        aes(clinical_v_untargeted$obs_aob_dist),
        fill = aob_col,
        col = "black"
      ) +
      labs(x = "AOB",
           title = "Clinical vs Untargeted WW",
           fill = "")
    
  }
  
  ## Targeted vs untargeted ----
  if (is.null(targeted_v_untargeted)) {
    aob1_targeted_untargeted <- NULL
  }
  else{
    aob1_targeted_untargeted <- ggplot() +
      geom_histogram(
        aes(targeted_v_untargeted$obs_aob_dist),
        fill = aob_col,
        col = "black"
      ) +
      labs(x = "AOB",
           title = "Targeted vs Unargeted WW",
           fill = "")
    
  }
  
  ## Combine ----
  aob1_plots <- list(aob1_clinical_targeted,
                     aob1_clinical_untargeted,
                     aob1_targeted_untargeted)
  aob1_plots <- Filter(Negate(is.null), aob1_plots)
  aob1_plot <- Reduce(`+`, aob1_plots) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "bottom")
  
  # ---- Cor ----
  ## Clinical vs targeted ----
  if (is.null(clinical_v_targeted)) {
    cor_clinical_targeted <- NULL
  }
  else{
    cor_clinical_targeted <- ggplot() +
      geom_histogram(aes(clinical_v_targeted$obs_cor_dist),
                     fill = cor_col,
                     col = "black") +
      labs(x = "Tau",
           title = "Clinical vs Targeted WW",
           fill = "")
  }
  
  ## Clinical vs untargeted ----
  if (is.null(clinical_v_untargeted)) {
    cor_clinical_untargeted <- NULL
  }
  else{
    cor_clinical_untargeted <- ggplot() +
      geom_histogram(
        aes(clinical_v_untargeted$obs_cor_dist),
        fill = cor_col,
        col = "black"
      ) +
      labs(x = "Tau",
           title = "Clinical vs Untargeted WW",
           fill = "")
    
  }
  
  ## Targeted vs untargeted ----
  if (is.null(targeted_v_untargeted)) {
    cor_targeted_untargeted <- NULL
  }
  else{
    cor_targeted_untargeted <- ggplot() +
      geom_histogram(
        aes(targeted_v_untargeted$obs_cor_dist),
        fill = cor_col,
        col = "black"
      ) +
      labs(x = "Tau",
           title = "Targeted vs Unargeted WW",
           fill = "")
    
  }
  
  ## Combine ----
  cor_plots <- list(cor_clinical_targeted,
                    cor_clinical_untargeted,
                    cor_targeted_untargeted)
  cor_plots <- Filter(Negate(is.null), cor_plots)
  cor_plot <- Reduce(`+`, cor_plots) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "bottom")
  
  # ---- ALL ----
  col_header <- function(label) {
    p <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = label,
        fontface = "bold",
        size = 5
      ) +
      theme_void()
    wrap_elements(full = p, ignore_tag = TRUE)
  }
  
  header_labels <- c(
    if (!is.null(clinical_v_targeted))
      "Clinical vs Targeted WW",
    if (!is.null(clinical_v_untargeted))
      "Clinical vs Untargeted WW",
    if (!is.null(targeted_v_untargeted))
      "Targeted WW vs Untargeted WW"
  )
  headers <- lapply(header_labels, col_header)
  
  all_plot <- Reduce(`|`, headers) /
    Reduce(`+`, mae_plots) / Reduce(`+`, aob1_plots) / Reduce(`+`, cor_plots) +
    plot_layout(heights = c(0.08, 1, 1, 1)) +
    plot_annotation(tag_levels = "a") &
    labs(title = "")
  
  return(all_plot)
}

#' Plot distribution of difference between obs and null posterior draws
#' @param clinical_v_targeted Output from compare_posteriors() for clinical vs targeted ww data
#' @param clinical_v_untargeted Output from compare_posteriors() for clinical vs untargeted ww data
#' @param targeted_v_untargeted Output from compare_posteriors() for targeted vs untargeted ww data
#'
#' @return A patchwork of all concordance metrics (observed - null distributions)
plot_obs_null_concordance <- function(clinical_v_targeted,
                                      clinical_v_untargeted,
                                      targeted_v_untargeted) {
  # ---- MAE ----
  ## Clinical vs targeted ----
  if (is.null(clinical_v_targeted)) {
    mae_clinical_targeted <- NULL
  }
  else{
    mae_clinical_targeted <- ggplot() +
      geom_histogram(
        aes(
          clinical_v_targeted$obs_mae_dist - clinical_v_targeted$null_mae_dist
        ),
        fill = mae_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed MAE - Null MAE",
           title = "Clinical vs Targeted WW",
           fill = "")
  }
  
  ## Clinical vs untargeted ----
  if (is.null(clinical_v_untargeted)) {
    mae_clinical_untargeted <- NULL
  }
  else{
    mae_clinical_untargeted <- ggplot() +
      geom_histogram(
        aes(
          clinical_v_untargeted$obs_mae_dist - clinical_v_untargeted$null_mae_dist
        ),
        fill = mae_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed MAE - Null MAE",
           title = "Clinical vs Untargeted WW",
           fill = "")
    
  }
  
  ## Targeted vs untargeted ----
  if (is.null(targeted_v_untargeted)) {
    mae_targeted_untargeted <- NULL
  }
  else{
    mae_targeted_untargeted <- ggplot() +
      geom_histogram(
        aes(
          targeted_v_untargeted$obs_mae_dist - targeted_v_untargeted$null_mae_dist
        ),
        fill = mae_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed MAE - Null MAE",
           title = "Targeted vs Unargeted WW",
           fill = "")
    
  }
  
  ## Combine ----
  mae_plots <- list(mae_clinical_targeted,
                    mae_clinical_untargeted,
                    mae_targeted_untargeted)
  mae_plots <- Filter(Negate(is.null), mae_plots)
  mae_plot <- Reduce(`+`, mae_plots) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "bottom")
  
  # ---- AOB ----
  ## Clinical vs targeted ----
  if (is.null(clinical_v_targeted)) {
    aob1_clinical_targeted <- NULL
  }
  else{
    aob1_clinical_targeted <- ggplot() +
      geom_histogram(
        aes(
          clinical_v_targeted$obs_aob_dist - clinical_v_targeted$null_aob_dist
        ),
        fill = aob_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed AOB - Null AOB",
           title = "Clinical vs Targeted WW",
           fill = "")
  }
  
  ## Clinical vs untargeted ----
  if (is.null(clinical_v_untargeted)) {
    aob1_clinical_untargeted <- NULL
  }
  else{
    aob1_clinical_untargeted <- ggplot() +
      geom_histogram(
        aes(
          clinical_v_untargeted$obs_aob_dist - clinical_v_untargeted$null_aob_dist
        ),
        fill = aob_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed AOB - Null AOB",
           title = "Clinical vs Untargeted WW",
           fill = "")
    
  }
  
  ## Targeted vs untargeted ----
  if (is.null(targeted_v_untargeted)) {
    aob1_targeted_untargeted <- NULL
  }
  else{
    aob1_targeted_untargeted <- ggplot() +
      geom_histogram(
        aes(
          targeted_v_untargeted$obs_aob_dist - targeted_v_untargeted$null_aob_dist
        ),
        fill = aob_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed AOB - Null AOB",
           title = "Targeted vs Unargeted WW",
           fill = "")
    
  }
  
  ## Combine ----
  aob1_plots <- list(aob1_clinical_targeted,
                     aob1_clinical_untargeted,
                     aob1_targeted_untargeted)
  aob1_plots <- Filter(Negate(is.null), aob1_plots)
  aob1_plot <- Reduce(`+`, aob1_plots) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "bottom")
  
  # ---- Cor ----
  ## Clinical vs targeted ----
  if (is.null(clinical_v_targeted)) {
    cor_clinical_targeted <- NULL
  }
  else{
    cor_clinical_targeted <- ggplot() +
      geom_histogram(
        aes(
          clinical_v_targeted$obs_cor_dist - clinical_v_targeted$null_cor_dist
        ),
        fill = cor_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed Tau - Null Tau",
           title = "Clinical vs Targeted WW",
           fill = "")
  }
  
  ## Clinical vs untargeted ----
  if (is.null(clinical_v_untargeted)) {
    cor_clinical_untargeted <- NULL
  }
  else{
    cor_clinical_untargeted <- ggplot() +
      geom_histogram(
        aes(
          clinical_v_untargeted$obs_cor_dist - clinical_v_untargeted$null_cor_dist
        ),
        fill = cor_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed Tau - Null Tau",
           title = "Clinical vs Untargeted WW",
           fill = "")
    
  }
  
  ## Targeted vs untargeted ----
  if (is.null(targeted_v_untargeted)) {
    cor_targeted_untargeted <- NULL
  }
  else{
    cor_targeted_untargeted <- ggplot() +
      geom_histogram(
        aes(
          targeted_v_untargeted$obs_cor_dist - targeted_v_untargeted$null_cor_dist
        ),
        fill = cor_col,
        col = "black"
      ) +
      geom_vline(
        aes(xintercept = 0),
        col = "#b3b3b3",
        linetype = "dashed",
        size = 1
      ) +
      labs(x = "Observed Tau - Null Tau",
           title = "Targeted vs Unargeted WW",
           fill = "")
    
  }
  
  ## Combine ----
  cor_plots <- list(cor_clinical_targeted,
                    cor_clinical_untargeted,
                    cor_targeted_untargeted)
  cor_plots <- Filter(Negate(is.null), cor_plots)
  cor_plot <- Reduce(`+`, cor_plots) +
    plot_layout(guides = "collect", axis_titles = "collect") &
    theme(legend.position = "bottom")
  
  
  # ---- ALL ----
  col_header <- function(label) {
    p <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = label,
        fontface = "bold",
        size = 5
      ) +
      theme_void()
    wrap_elements(full = p, ignore_tag = TRUE)
  }
  
  header_labels <- c(
    if (!is.null(clinical_v_targeted))
      "Clinical vs Targeted WW",
    if (!is.null(clinical_v_untargeted))
      "Clinical vs Untargeted WW",
    if (!is.null(targeted_v_untargeted))
      "Targeted WW vs Untargeted WW"
  )
  headers <- lapply(header_labels, col_header)
  
  all_plot <- Reduce(`|`, headers) /
    Reduce(`+`, mae_plots) / Reduce(`+`, aob1_plots) / Reduce(`+`, cor_plots) +
    plot_layout(heights = c(0.08, 1, 1, 1)) +
    plot_annotation(tag_levels = "a") &
    labs(title = "")
  
  return(all_plot)
}
