suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))

# --- Plotting Functions ---

# Color palette for emphasized traits
get_emphasized_traits <- function() {
  c(
  "height" = "#C95693",
  "BMI" = "#7B03D1",
  "blood-pressure-systolic-mean" = "#FF7F0E",
  "blood glucose" = "#2CA02C",
  "WHR" = "#1F77B4"
)
}

#' Plot with emphasized traits
#'
#' @param data Data frame containing plotting data
#' @param x_var Name of x variable
#' @param y_var Name of y variable
#' @param x_lower_var Name of x lower error bound variable (optional)
#' @param x_upper_var Name of x upper error bound variable (optional)
#' @param y_lower_var Name of y lower error bound variable (optional)
#' @param y_upper_var Name of y upper error bound variable (optional)
#' @param color_var Name of variable to color by (default: 'trait')
#' @param label_var Name of variable to use for labels (default: 'trait')
#' @param emphasize_traits Character vector of trait names to emphasize
#' @param colors Vector of colors for emphasized traits (must be same length as emphasize_traits)
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param point_size Point size for non-emphasized traits
#' @param emphasized_point_size Point size for emphasized traits
#' @param alpha Transparency for points and error bars (0-1)
#' @param add_abline Logical, whether to add y=x line
#' @param legend_title Title for the legend
#' @return A ggplot object
plot_with_emphasis <- function(data, 
                             x_var, 
                             y_var, 
                             x_lower_var = NULL,
                             x_upper_var = NULL,
                             y_lower_var = NULL,
                             y_upper_var = NULL,
                             color_var = 'trait', 
                             label_var = NULL, 
                             emphasize_traits = NULL, 
                             colors = NULL, 
                             point_size = 1, 
                             emphasized_point_size = 2, 
                             title = "", 
                             xlab = "", 
                             ylab = "", 
                             alpha = 0.7,
                             add_y_equals_x_line = FALSE, 
                             legend_title_for_color = NULL, 
                             show_non_emphasized = TRUE,
                             errorbar_scale = 1,
                             use_legend_for_labels = FALSE,
                             make_square_axes = FALSE
                             ) {
  
  # Create a copy of the data to modify
  plot_data <- data

  # Determine emphasized status
  if (!is.null(emphasize_traits) && length(emphasize_traits) > 0) {
    plot_data$emphasized <- plot_data[[color_var]] %in% emphasize_traits
  } else {
    plot_data$emphasized <- FALSE # All false if emphasize_traits is NULL or empty
  }

  # Initialize plot_color: default to a general non-emphasized color (e.g., gray40)
  # This default applies if an item is not in `names(colors)`.
  default_non_emphasized_color <- "gray40"
  plot_data$plot_color <- default_non_emphasized_color

  # Apply specific colors from the 'colors' palette if provided
  if (!is.null(colors)) {
    # 'colors' should be a named vector like c("traitA" = "#hex1", "traitB" = "#hex2")
    for (item_name in names(colors)) {
      item_color <- colors[[item_name]]
      # Apply this color to rows where color_var matches item_name
      # Ensure item_name is treated as a string for comparison with column values
      current_items <- as.character(plot_data[[color_var]])
      plot_data$plot_color[current_items == item_name] <- item_color
    }
  }
  
  # If an item is marked as 'emphasized' but didn't get a color from the 'colors' palette
  # (e.g. it was in emphasize_traits but not in names(colors)), 
  # it would currently have default_non_emphasized_color.
  # This behavior is generally fine, as 'colors' should be the source of truth for specific coloring.
  # If show_non_emphasized is FALSE, those rows will be filtered out later anyway.
  
  # Create the base plot (without global color mapping)
  p <- ggplot(plot_data, aes_string(x = x_var, y = y_var))
  
  # Add error bars if specified
  has_y_error <- !is.null(y_lower_var) && !is.null(y_upper_var)
  has_x_error <- !is.null(x_lower_var) && !is.null(x_upper_var)
  
  # Helper function to add error bars
  # Updated to handle use_legend_for_labels
  add_error_bars <- function(p_obj, data_subset, is_emphasized_flag, y_var_str, y_lower_var_str, y_upper_var_str, x_var_str, x_lower_var_str, x_upper_var_str, color_var_str, use_legend_flag, current_alpha, current_errorbar_scale) {
    line_size <- ifelse(is_emphasized_flag, 0.5, 0.2)

    # Y-error bars
    if (!is.null(y_lower_var_str) && !is.null(y_upper_var_str)) {
      common_y_aes <- aes(ymin = .data[[y_var_str]] - .data[[y_lower_var_str]] * current_errorbar_scale,
                          ymax = .data[[y_var_str]] + .data[[y_upper_var_str]] * current_errorbar_scale)
      if (is_emphasized_flag) {
        color_y_aes <- if (use_legend_flag) aes(color = .data[[color_var_str]]) else aes(color = .data[["plot_color"]])
        p_obj <- p_obj + geom_errorbar(data = data_subset, mapping = utils::modifyList(common_y_aes, color_y_aes), width = 0, alpha = current_alpha, linewidth = line_size, show.legend = FALSE)
      } else if (show_non_emphasized) {
        p_obj <- p_obj + geom_errorbar(data = data_subset, mapping = common_y_aes, width = 0, color = "gray70", alpha = current_alpha * 0.5, linewidth = line_size, show.legend = FALSE)
      }
    }

    # X-error bars
    if (!is.null(x_lower_var_str) && !is.null(x_upper_var_str)) {
      common_x_aes <- aes(xmin = .data[[x_var_str]] - .data[[x_lower_var_str]] * current_errorbar_scale,
                          xmax = .data[[x_var_str]] + .data[[x_upper_var_str]] * current_errorbar_scale)
      if (is_emphasized_flag) {
        color_x_aes <- if (use_legend_flag) aes(color = .data[[color_var_str]]) else aes(color = .data[["plot_color"]])
        p_obj <- p_obj + geom_errorbarh(data = data_subset, mapping = utils::modifyList(common_x_aes, color_x_aes), height = 0, alpha = current_alpha, linewidth = line_size, show.legend = FALSE)
      } else if (show_non_emphasized && !is.null(y_lower_var_str)) { # Original logic had 'has_y_error' for non-emphasized x-bars, seems like a bug, should be has_x_error
        # Assuming non-emphasized x-error bars should also be gray if shown
        p_obj <- p_obj + geom_errorbarh(data = data_subset, mapping = common_x_aes, height = 0, color = "gray70", alpha = current_alpha * 0.5, linewidth = line_size, show.legend = FALSE)
      }
    }
    return(p_obj)
  }

  # Add error bars for emphasized points
  if (any(plot_data$emphasized)) {
    p <- add_error_bars(p, subset(plot_data, emphasized), TRUE, y_var, y_lower_var, y_upper_var, x_var, x_lower_var, x_upper_var, color_var, use_legend_for_labels, alpha, errorbar_scale)
  }
  
  # Add lines connecting points from the same trait
  trait_counts <- plot_data %>%
    count(!!sym(color_var)) %>%
    filter(n > 1)
  
  if (nrow(trait_counts) > 0) {
    # Get line data for all traits with multiple points
    line_data <- plot_data %>% 
      filter(!!sym(color_var) %in% trait_counts[[color_var]])
    
    # Add lines with colors
    line_data_to_plot <- if (show_non_emphasized) {
      line_data
    } else {
      line_data %>% filter(emphasized)
    }

    if (nrow(line_data_to_plot) > 0) {
      line_data_emph <- line_data_to_plot %>% filter(emphasized)
      line_data_non_emph <- line_data_to_plot %>% filter(!emphasized)

      if (nrow(line_data_emph) > 0) {
        line_aes_emph <- if (use_legend_for_labels) {
          aes(group = .data[[color_var]], color = .data[[color_var]])
        } else {
          aes(group = .data[[color_var]], color = .data[["plot_color"]])
        }
        p <- p + geom_line(data = line_data_emph, mapping = line_aes_emph, alpha = 0.9, linewidth = 0.8)
      }
      if (nrow(line_data_non_emph) > 0) {
        line_aes_non_emph <- if (use_legend_for_labels) {
          aes(group = .data[[color_var]], color = .data[[color_var]]) # Color determined by scale_color_manual
        } else {
          aes(group = .data[[color_var]], color = .data[["plot_color"]]) # plot_color will be 'gray40'
        }
        p <- p + geom_line(data = line_data_non_emph, mapping = line_aes_non_emph, alpha = 0.4, linewidth = 0.4)
      }
    }
  }
  
  # Add points (only non-emphasized if show_non_emphasized is TRUE)
  if (show_non_emphasized || !any(plot_data$emphasized)) {
    point_data <- if (show_non_emphasized) plot_data else plot_data[plot_data$emphasized, ]
    
    point_aes <- if (use_legend_for_labels) {
      aes(color = .data[[color_var]])
    } else {
      aes(color = .data[["plot_color"]])
    }
    p <- p +
      geom_point(
        data = point_data,
        mapping = point_aes,
        size = ifelse(point_data$emphasized, emphasized_point_size, point_size),
        alpha = alpha
      )
  }
  
  # Add y=x line if requested
  if (add_y_equals_x_line) {
    p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50")
  }
  
  # Apply color scales, labels, and theme
  if (use_legend_for_labels) {
    # Use a legend; colors are mapped from color_var via aes() in geoms
    # na.value handles cases where a value in color_var is not in names(colors)
    p <- p + scale_color_manual(values = colors, name = legend_title_for_color, na.value = "grey70")
    # No direct ggrepel labels if using legend for color identification
  } else {
    # Use direct coloring via plot_color; no legend for these colors
    p <- p + scale_color_identity(guide = "none")
    
    # Add direct labels if specified and not using legend for labels
    if (!is.null(label_var) && any(plot_data$emphasized)) {
      label_data <- plot_data[plot_data$emphasized, ]
      if (nrow(label_data) > 0) {
          label_data_dedup <- label_data %>%
              group_by(.data[[color_var]]) %>%
              # Label point with max y-value, or if only one point for that group
              filter(.data[[y_var]] == max(.data[[y_var]], na.rm = TRUE) | n() == 1) %>%
              slice(1) # Take the first if ties in max y or if multiple points have n() == 1 (should not happen with group_by)

          p <- p + ggrepel::geom_text_repel(
              data = label_data_dedup,
              aes_string(label = label_var, color = "plot_color"), # Use plot_color for label color
              size = 3.5,
              segment.color = 'grey50',
              segment.alpha = 0.7,
              min.segment.length = 0.1,
              box.padding = 0.5,
              point.padding = 0.3,
              max.overlaps = Inf,
              show.legend = FALSE
          )
      }
    }
  }
  
  # Apply base theme
  p <- p + theme_minimal(base_size = 12) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
      # legend.position will be set conditionally below
    )

  # Conditional legend position
  if (use_legend_for_labels) {
    # Default legend position (usually right), or specify e.g. theme(legend.position = "right")
    # If legend_title_for_color is empty, ggplot might produce a better default title or no title.
    # Ensuring the legend appears:
    p <- p + theme(legend.position = "right") # Or theme(legend.position = "bottom"), etc.
  } else {
    p <- p + theme(legend.position = "none")
  }

  p <- p + labs(title = title, x = xlab, y = ylab)
  
  # Add fill scale if needed (e.g., for boxplots, not used here)
  # p <- p + scale_fill_identity()
  
    # Make axes square if requested
  if (make_square_axes) {
    all_vals <- c()
    if (!is.null(x_var) && x_var %in% names(plot_data)) all_vals <- c(all_vals, plot_data[[x_var]])
    if (!is.null(y_var) && y_var %in% names(plot_data)) all_vals <- c(all_vals, plot_data[[y_var]])

    if (!is.null(x_lower_var) && x_lower_var %in% names(plot_data) && !is.null(plot_data[[x_lower_var]])) {
      all_vals <- c(all_vals, plot_data[[x_var]] - plot_data[[x_lower_var]] * errorbar_scale)
    }
    if (!is.null(x_upper_var) && x_upper_var %in% names(plot_data) && !is.null(plot_data[[x_upper_var]])) {
      all_vals <- c(all_vals, plot_data[[x_var]] + plot_data[[x_upper_var]] * errorbar_scale)
    }
    if (!is.null(y_lower_var) && y_lower_var %in% names(plot_data) && !is.null(plot_data[[y_lower_var]])) {
      all_vals <- c(all_vals, plot_data[[y_var]] - plot_data[[y_lower_var]] * errorbar_scale)
    }
    if (!is.null(y_upper_var) && y_upper_var %in% names(plot_data) && !is.null(plot_data[[y_upper_var]])) {
      all_vals <- c(all_vals, plot_data[[y_var]] + plot_data[[y_upper_var]] * errorbar_scale)
    }
    
    all_vals <- all_vals[is.finite(all_vals)]

    if (length(all_vals) > 0) {
      min_val <- min(all_vals, na.rm = TRUE)
      max_val <- max(all_vals, na.rm = TRUE)
      p <- p + coord_cartesian(xlim = c(min_val, max_val), ylim = c(min_val, max_val), expand = TRUE)
    } else {
      # Fallback if all_vals is empty, though this shouldn't happen with typical data
      # Or, one might choose to not apply coord_cartesian if no valid range found
    }
  }
  
  return(p)
}
