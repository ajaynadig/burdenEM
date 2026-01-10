#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))

# Parse command-line argument: results directory
parser <- ArgumentParser(description = "Plotting script for simulations")
parser$add_argument("results_dir", help = "Directory containing results tables used for plotting")
args <- parser$parse_args()

results_dir <- args$results_dir
tables_dir  <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Helper to find exactly one file matching the given pattern (anchored suffix)
find_one <- function(dir, pattern) {
  matches <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
  if (length(matches) != 1) return(NA_character_)
  matches[[1]]
}

# Plotting function (distribution): one page per abbreviation
plot_distribution <- function(df) {
  out_pdf <- file.path(figures_dir, "distribution.pdf")
  abbrs <- unique(df$abbreviation)
  pdf(out_pdf)
  for (abbr in abbrs) {
    s <- df[df$abbreviation == abbr, , drop = FALSE]
    title_txt <- gsub("_", " ", abbr)

    p <- ggplot(s, aes(x = needed_genes)) +
      geom_path(aes(y = `variance_mean.true`, color = "True")) +
      geom_path(aes(y = `variance_mean.estimated`, color = "Estimated")) +
      geom_errorbar(aes(ymin = `variance_mean.estimated` - `variance_sd.estimated`,
                        ymax = `variance_mean.estimated` + `variance_sd.estimated`),
                    width = 0) +
      scale_x_log10() +
      coord_cartesian(xlim = c(1, max(s$needed_genes)), ylim = c(0, 1)) +
      scale_color_manual(values = c("Estimated" = "black", "True" = "red"), name = NULL) +
      labs(x = "fraction of genes", y = "fraction of variance", title = title_txt) +
      theme_minimal()

    print(p)
  }
  dev.off()
  cat(out_pdf, "\n")
}

plot_heritability <- function(df) {
  p_total <- plot_prefix_estimates(df, prefix = "total_h2", axis_label = "Burden heritability", min_x = 0)
  p_pos   <- plot_prefix_estimates(df, prefix = "positive_h2", axis_label = "Positive-effect heritability", min_x = 0)
  list(total_h2 = p_total, positive_h2 = p_pos)
}

# Plot calibration for mean standard error of heritability
plot_heritability_se <- function(df) {
  p_se <- plot_prefix_estimates(
    df,
    prefix = "total_h2",
    axis_label = "Burden heritability SE",
    min_x = 0,
    col_est = "total_h2_se_rms.estimated",
    col_true = "total_h2_sd.estimated",
    col_sd = NULL
  )
  list(total_h2_se = p_se)
}

plot_polygenicity <- function(df) {
  metrics <- c(
    "softmax_polygenicity",
    "effective_polygenicity",
    "entropy_polygenicity",
    "effective_effect_var_polygenicity"
  )
  plots <- list()
  for (m in metrics) {
    tmp <- data.frame(abbreviation = df$abbreviation, stringsAsFactors = FALSE)
    tmp[[paste0(m, "_mean.estimated")]] <- df[[paste0(m, "_mean_log10.estimated")]]
    tmp[[paste0(m, "_mean.true")]]      <- df[[paste0(m, "_mean_log10.true")]]
    sd_name <- paste0(m, "_sd.estimated")
    if (sd_name %in% names(df)) {
      tmp[[sd_name]] <- df[[sd_name]]
    }
    max_log <- log10(20000)
    plots[[m]] <- plot_prefix_estimates(
      tmp,
      prefix = m,
      axis_label = "Polygenicity (genes)",
      min_x = 0,
      max_x = max_log,
      label_linear_from_log10 = TRUE
    )
  }
  plots
}

plot_calibration <- function(df) {
  out_pdf <- file.path(figures_dir, "calibration.pdf")
  abbrs <- unique(df$abbreviation)
  pdf(out_pdf)
  for (abbr in abbrs) {
    s <- df[df$abbreviation == abbr, , drop = FALSE]
    s <- s[order(s$`post_mean.mean`), , drop = FALSE]
    title_txt <- gsub("_", " ", abbr)

    min_x <- min(s$`post_mean.mean`)
    max_x <- max(s$`post_mean.mean`)

    p <- ggplot(s, aes(x = `post_mean.mean`, y = `true_mean.mean`)) +
      geom_point() +
      geom_path() +
      # y error bars (vertical)
      geom_errorbar(aes(ymin = `true_mean.mean` - `true_mean.sd`,
                        ymax = `true_mean.mean` + `true_mean.sd`),
                    width = 0) +
      # x error bars (horizontal)
      geom_segment(aes(x = `post_mean.mean` - `post_mean.sd`,
                       xend = `post_mean.mean` + `post_mean.sd`,
                       y = `true_mean.mean`, yend = `true_mean.mean`)) +
      # y = x reference line spanning x-range
      geom_segment(x = min_x, xend = max_x, y = min_x, yend = max_x, linetype = "dashed") +
      labs(x = "estimated posterior effect size",
           y = "actual effect size",
           title = title_txt) +
      theme_minimal() +
      coord_equal(xlim = c(min_x, max_x), ylim = c(min_x, max_x))

    print(p)
  }
  dev.off()
  cat(out_pdf, "\n")
}

# Shared helper for heritability/polygenicity style plots
# Inputs:
#   df: data.frame containing columns 'abbreviation',
#       paste0(prefix, "_mean.estimated"), paste0(prefix, "_mean.true"),
#       and optionally paste0(prefix, "_sd.estimated").
#   prefix: column name prefix (character)
#   axis_label: optional x-axis label (defaults to prefix)
#   min_x: optional minimum x-axis value (numeric)
plot_prefix_estimates <- function(df, prefix, axis_label = NULL, min_x = NULL, max_x = NULL, label_linear_from_log10 = FALSE, col_est = NULL, col_true = NULL, col_sd = NULL) {
  if (is.null(axis_label)) axis_label <- prefix

  col_est  <- if (is.null(col_est)) paste0(prefix, "_mean.estimated") else col_est
  col_true <- if (is.null(col_true)) paste0(prefix, "_mean.true") else col_true
  col_sd   <- if (is.null(col_sd))   paste0(prefix, "_sd.estimated")   else col_sd
  has_sd   <- !is.null(col_sd) && col_sd %in% names(df)

  n <- nrow(df)
  d <- data.frame(
    abbreviation = df$abbreviation,
    y = seq_len(n),
    x_est = df[[col_est]],
    x_true = df[[col_true]],
    stringsAsFactors = FALSE
  )
  if (has_sd) d$x_sd <- df[[col_sd]]

  p <- ggplot(d, aes(y = y)) +
    # vertical segment at the true value for each row (solid line)
    geom_segment(aes(x = x_true, xend = x_true, y = y - 0.35, yend = y + 0.35,
                     color = "true", linetype = "true"))

  if (has_sd) {
    # horizontal error bar centered at estimate
    p <- p + geom_segment(aes(x = x_est - x_sd, xend = x_est + x_sd, y = y, yend = y,
                              color = "estimated (s.d.)", linetype = "estimated (s.d.)"))
  }

  # y labels with underscores replaced by spaces
  y_labels <- gsub("_", " ", d$abbreviation)

  p <- p +
    geom_point(aes(x = x_est, color = "estimated (s.d.)"), size = 2) +
    scale_y_continuous(breaks = d$y, labels = y_labels) +
    labs(x = axis_label, y = NULL) +
    theme_minimal() +
    scale_color_manual(values = c("true" = "grey40", "estimated (s.d.)" = "black"), name = NULL) +
    scale_linetype_manual(values = c("true" = "solid", "estimated (s.d.)" = "solid"), name = NULL)

  # optional axis limits
  xlo <- if (is.null(min_x)) NA else min_x
  xhi <- if (is.null(max_x)) NA else max_x
  p <- p + coord_cartesian(xlim = c(xlo, xhi), ylim = c(1, n))

  # optionally label x-axis in linear space when data are log10
  if (label_linear_from_log10) {
    rng_min <- if (!is.null(min_x)) min_x else suppressWarnings(min(c(d$x_est, d$x_true), na.rm = TRUE))
    rng_max <- if (!is.null(max_x)) max_x else suppressWarnings(max(c(d$x_est, d$x_true), na.rm = TRUE))
    brks <- pretty(c(rng_min, rng_max), n = 5)
    p <- p + scale_x_continuous(breaks = brks, labels = function(b) signif(10^b, 3))
  }

  p
}

# Validate required files (exact suffix matches inside tables/)
distribution_file <- find_one(tables_dir, "meta\\.distribution\\.tsv$")
heritability_file <- find_one(tables_dir, "meta\\.heritability\\.tsv$")
polygenicity_file <- find_one(tables_dir, "meta\\.polygenicity\\.tsv$")
calibration_file  <- find_one(tables_dir, "meta\\.calibration\\.tsv$")

# Heritability figures
if (!is.na(heritability_file)) {
  h2_df <- read.delim(heritability_file, sep = "\t", header = TRUE, check.names = FALSE)
  h2_plots <- plot_heritability(h2_df)
  out_h2_total <- file.path(figures_dir, "heritability_total_h2.pdf")
  ggsave(out_h2_total, h2_plots$total_h2)
  cat(out_h2_total, "\n")
  out_h2_pos <- file.path(figures_dir, "heritability_positive_h2.pdf")
  ggsave(out_h2_pos, h2_plots$positive_h2)
  cat(out_h2_pos, "\n")
  se_plot <- plot_heritability_se(h2_df)
  out_h2_se <- file.path(figures_dir, "heritability_se.pdf")
  ggsave(out_h2_se, se_plot$total_h2_se)
  cat(out_h2_se, "\n")
}

# Distribution figure (multi-page PDF)
if (!is.na(distribution_file)) {
  dist_df <- read.delim(distribution_file, sep = "\t", header = TRUE, check.names = FALSE)
  plot_distribution(dist_df)
}

# Polygenicity figures
if (!is.na(polygenicity_file)) {
  poly_df <- read.delim(polygenicity_file, sep = "\t", header = TRUE, check.names = FALSE)
  poly_plots <- plot_polygenicity(poly_df)
  for (nm in names(poly_plots)) {
    out_poly <- file.path(figures_dir, paste0("polygenicity_", nm, ".pdf"))
    ggsave(out_poly, poly_plots[[nm]])
    cat(out_poly, "\n")
  }
}

# Calibration figure (multi-page PDF)
if (!is.na(calibration_file)) {
  cal_df <- read.delim(calibration_file, sep = "\t", header = TRUE, check.names = FALSE)
  plot_calibration(cal_df)
}
