###############################################################################
############################# 5 : PLOTS #######################################
###############################################################################

###############################################################################
# 1. GRID SEARCH PLOT
###############################################################################

# Plot grid search results (convex hull in grey, best guess in blue, targets in red dashed)
# Arguments:
#   - grid_result: output from grid_search() => list(grid=..., best_guess=...)
#   - targets: named list with the target values for x_var and y_var
#   - x_var, y_var: column names in grid_result$grid (and in best_guess) to plot
#   - x_target, y_target: names in targets for the dashed target lines
#   - title, subtitle: plot labels
#   - x_lab, y_lab: axis labels (optional)
plot_grid_search_generic <- function(grid_result, targets,
                                     x_var, y_var, 
                                     x_target, y_target,
                                     title, subtitle,
                                     x_lab = x_var, y_lab = y_var) {
  
  # Extract data
  df <- grid_result$grid
  best <- grid_result$best_guess
  
  # Keep only finite points (avoid chull / ggplot errors)
  ok <- is.finite(df[[x_var]]) & is.finite(df[[y_var]])
  df_ok <- df[ok, , drop = FALSE]
  
  # If everything is NA/Inf, return an empty plot with message
  if (nrow(df_ok) == 0) {
    p <- ggplot() +
      theme_bw() +
      labs(title = title, subtitle = "No finite points to plot (all NA/Inf).")
    return(p)
  }
  
  # Targets
  tx <- targets[[x_target]]
  ty <- targets[[y_target]]
  
  # Axis ranges (for label placement)
  x_min <- min(df_ok[[x_var]])
  x_max <- max(df_ok[[x_var]])
  y_min <- min(df_ok[[y_var]])
  y_max <- max(df_ok[[y_var]])
  x_pad <- (x_max - x_min) * 0.06
  y_pad <- (y_max - y_min) * 0.08
  
  # Best guess label uses parameter names (first two grid columns)
  param_names <- names(grid_result$grid)[1:2]
  best_label <- paste0(param_names[1], "=", signif(best[[param_names[1]]], 4),
                       "\n", param_names[2], "=", signif(best[[param_names[2]]], 4))
  
  # Plot
  p <- ggplot(df_ok, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    
    # Grey hull
    ggforce::geom_mark_hull(
      concavity = 20,
      expand = grid::unit(0, "mm"),
      radius = grid::unit(0, "mm"),
      fill = "grey80", alpha = 0.5, color = NA
    ) +
    
    # Cloud of points
    geom_point(alpha = 0.3, size = 0.5) +
    
    # Target dashed lines
    geom_vline(xintercept = tx, linetype = "dashed", color = "red") +
    geom_hline(yintercept = ty, linetype = "dashed", color = "red") +
    
    # Target values on axes (red)
    annotate("text", x = tx, y = -Inf, label = sprintf("%.4g", tx),
             color = "red", vjust = 1.6, size = 3) +
    annotate("text", x = -Inf, y = ty, label = sprintf("%.4g", ty),
             color = "red", hjust = 1.6, size = 3) +
    
    # Best guess in blue
    geom_point(data = best, aes(x = .data[[x_var]], y = .data[[y_var]]),
               color = "blue", size = 3) +
    geom_label(data = best, aes(x = .data[[x_var]], y = .data[[y_var]], label = best_label),
               color = "blue", fill = NA, linewidth = 0.3, hjust = -0.1, vjust = 0.5, size = 3) +
    
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(plot.margin = margin(10, 10, 30, 40)) +
    labs(x = x_lab, y = y_lab, title = title, subtitle = subtitle)
  
  return(p)
}


# 1) beta grid
plot_grid_search_beta <- function(grid_result, targets) {
  plot_grid_search_generic(
    grid_result = grid_result,
    targets = targets,
    x_var = "prevalence_h", y_var = "prevalence_c",
    x_target = "prevalence_h", y_target = "prevalence_c",
    title = "Grid search: beta_h, beta_c",
    subtitle = "Feasible carriage space",
    x_lab = "Carriage prevalence hospital",
    y_lab = "Carriage prevalence community"
  )
}

# 2) sigma grid
plot_grid_search_sigma <- function(grid_result, targets) {
  
  # Convert from per person per day -> per 100k per year
  scale <- 365 * 1e5
  
  grid2 <- grid_result
  grid2$grid <- within(grid2$grid, {
    incidence_h <- incidence_h * scale
    incidence_c <- incidence_c * scale
  })
  grid2$best_guess <- within(grid2$best_guess, {
    incidence_h <- incidence_h * scale
    incidence_c <- incidence_c * scale
  })
  
  targets2 <- targets
  targets2$incidence_h <- targets2$incidence_h * scale
  targets2$incidence_c <- targets2$incidence_c * scale
  
  plot_grid_search_generic(
    grid_result = grid2,
    targets = targets2,
    x_var = "incidence_h", y_var = "incidence_c",
    x_target = "incidence_h", y_target = "incidence_c",
    title = "Grid search: sigma_h, sigma_c",
    subtitle = "Feasible incidence space",
    x_lab = "Incidence hospital (per 100k per year)",
    y_lab = "Incidence community (per 100k per year)"
  )
}


# 3) k grid
plot_grid_search_k <- function(grid_result, targets) {
  plot_grid_search_generic(
    grid_result = grid_result,
    targets = targets,
    x_var = "recid_1", y_var = "recid_2",
    x_target = "recid_1", y_target = "recid_2",
    title = "Grid search: k_II, k_III",
    subtitle = "Feasible recurrence space",
    x_lab = "Recurrence rate 1",
    y_lab = "Recurrence rate 2+"
  )
}


###############################################################################
# 2. CALIBRATION PLOT
###############################################################################

make_table_theme <- function(nrow_df) {
  gridExtra::ttheme_default(
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5, fontsize = 9),
      bg_params = list(fill = rep(c("white", "grey95"), length.out = nrow_df))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10, col = "white"),
      bg_params = list(fill = "steelblue", col = "white")
    )
  )
}

set_cell_fill <- function(grob, fill) {
  if (!is.null(grob$children) && length(grob$children) > 0) {
    rect_idx <- which(vapply(grob$children, function(ch) inherits(ch, "rect"), logical(1)))
    if (length(rect_idx) >= 1) {
      grob$children[[rect_idx[1]]]$gp$fill <- fill
      return(grob)
    }
  }
  grob$gp$fill <- fill
  grob
}

apply_row_fill <- function(tab, rows, fill, ncol_df) {
  for (r in rows) {
    for (c in seq_len(ncol_df)) {
      idx <- which(tab$layout$t == r + 1 & tab$layout$l == c)
      if (length(idx) == 1) {
        tab$grobs[[idx]] <- set_cell_fill(tab$grobs[[idx]], fill)
      }
    }
  }
  tab
}

apply_row_bold <- function(tab, rows, ncol_df) {
  for (r in rows) {
    for (c in seq_len(ncol_df)) {
      idx <- which(tab$layout$t == r + 1 & tab$layout$l == c)
      if (length(idx) == 1) {
        tab$grobs[[idx]]$gp$fontface <- "bold"
      }
    }
  }
  tab
}

make_pretty_table <- function(df) {
  gridExtra::tableGrob(df, rows = NULL, theme = make_table_theme(nrow(df)))
}

plot_calibration_table1 <- function(params_calib, metrics_calib, target_metrics, show = TRUE) {
  
  # Helper formatting
  fmt <- function(x) ifelse(is.na(x), "", sprintf("%.4g", x))
  fmt_pct <- function(x) ifelse(is.na(x), "", sprintf("%+.1f%%", x))
  
  # Outputs (last time point)
  prev_h <- tail(metrics_calib$carriage$prev_h, 1) * 100
  prev_c <- tail(metrics_calib$carriage$prev_c, 1) * 100
  inc_h  <- tail(metrics_calib$incidence_instant$inc_h_total_abs, 1) * 365 * 1e5
  inc_c  <- tail(metrics_calib$incidence_instant$inc_c_total_abs, 1) * 365 * 1e5
  rec1   <- tail(metrics_calib$recurrence$rec1, 1) * 100
  rec2   <- tail(metrics_calib$recurrence$rec2, 1) * 100
  
  # Targets
  t_prev_h <- target_metrics$prevalence_h * 100
  t_prev_c <- target_metrics$prevalence_c * 100
  t_inc_h  <- target_metrics$incidence_h * 365 * 1e5
  t_inc_c  <- target_metrics$incidence_c * 365 * 1e5
  t_rec1   <- target_metrics$recid_1 * 100
  t_rec2   <- target_metrics$recid_2 * 100
  
  # Relative error (%)
  rel <- function(val, target) {
    if (is.na(target) || target == 0) return(NA_real_)
    100 * (val / target - 1)
  }
  
  df <- data.frame(
    Item = c(
      "CALIBRATED PARAMETERS:",
      "  β_h (transmission hospital)",
      "  β_c (transmission community)",
      "  σ_h (progression hospital)",
      "  σ_c (progression community)",
      "  k_II (recurrence factor 1)",
      "  k_III (recurrence factor 2+)",
      "",
      "OUTPUTS:",
      "  Carriage hospital",
      "  Carriage community",
      "  Incidence hospital",
      "  Incidence community",
      "  Recurrence rate 1",
      "  Recurrence rate 2+"
    ),
    Value = c(
      "",
      fmt(params_calib["beta_h"]),
      fmt(params_calib["beta_c"]),
      fmt(params_calib["sigma_h"]),
      fmt(params_calib["sigma_c"]),
      fmt(params_calib["k_II"]),
      fmt(params_calib["k_III"]),
      "",
      "",
      fmt(prev_h),
      fmt(prev_c),
      fmt(inc_h),
      fmt(inc_c),
      fmt(rec1),
      fmt(rec2)
    ),
    Units = c(
      "",
      "day^-1", "day^-1", "day^-1", "day^-1", "-", "-",
      "",
      "",
      "%", "%", "/100k/year", "/100k/year", "%", "%"
    ),
    Target = c(
      "",
      "", "", "", "", "", "",
      "",
      "",
      fmt(t_prev_h),
      fmt(t_prev_c),
      fmt(t_inc_h),
      fmt(t_inc_c),
      fmt(t_rec1),
      fmt(t_rec2)
    ),
    RelError = c(
      "",
      "", "", "", "", "", "",
      "",
      "",
      fmt_pct(rel(prev_h, t_prev_h)),
      fmt_pct(rel(prev_c, t_prev_c)),
      fmt_pct(rel(inc_h, t_inc_h)),
      fmt_pct(rel(inc_c, t_inc_c)),
      fmt_pct(rel(rec1, t_rec1)),
      fmt_pct(rel(rec2, t_rec2))
    ),
    stringsAsFactors = FALSE
  )
  
  tab <- make_pretty_table(df)
  
  if (show) gridExtra::grid.arrange(tab)
  invisible(tab)
}

plot_calibration_table2 <- function(metrics_calib, target_metrics, show = TRUE) {
  
  fmt <- function(x) ifelse(is.na(x), "", sprintf("%.4g", x))
  
  # Basic values
  N_h <- metrics_calib$population$N_h[1]
  N_c <- metrics_calib$population$N_c[1]
  N_tot <- N_h + N_c
  
  prev_h <- tail(metrics_calib$carriage$prev_h, 1)
  prev_c <- tail(metrics_calib$carriage$prev_c, 1)
  prev_tot <- (prev_h * N_h + prev_c * N_c) / N_tot
  
  scale_inc <- 365 * 1e5
  
  inc_h_abs <- tail(metrics_calib$incidence_instant$inc_h_total_abs, 1) * scale_inc
  inc_c_abs <- tail(metrics_calib$incidence_instant$inc_c_total_abs, 1) * scale_inc
  inc_tot_abs <- inc_h_abs + inc_c_abs
  inc_h_rel <- tail(metrics_calib$incidence_instant$inc_h_total_rel, 1) * scale_inc
  inc_c_rel <- tail(metrics_calib$incidence_instant$inc_c_total_rel, 1) * scale_inc
  
  inc_h_primo_abs <- tail(metrics_calib$incidence_instant$inc_h_primo_abs, 1) * scale_inc
  inc_c_primo_abs <- tail(metrics_calib$incidence_instant$inc_c_primo_abs, 1) * scale_inc
  inc_primo_abs <- inc_h_primo_abs + inc_c_primo_abs
  
  inc_h_rec_abs <- tail(metrics_calib$incidence_instant$inc_h_rec_abs, 1) * scale_inc
  inc_c_rec_abs <- tail(metrics_calib$incidence_instant$inc_c_rec_abs, 1) * scale_inc
  inc_rec_abs <- inc_h_rec_abs + inc_c_rec_abs
  
  rec1 <- tail(metrics_calib$recurrence$rec1, 1) * 100
  rec2 <- tail(metrics_calib$recurrence$rec2, 1) * 100
  
  # R0h / R0c from NGM blocks if available
  R0h <- NA_real_
  R0c <- NA_real_
  if (!is.null(metrics_calib$R0$K)) {
    K <- metrics_calib$R0$K
    if (nrow(K) >= 14) {
      K_h <- K[1:7, 1:7]
      K_c <- K[8:14, 8:14]
      R0h <- max(Mod(eigen(K_h, only.values = TRUE)$values))
      R0c <- max(Mod(eigen(K_c, only.values = TRUE)$values))
    }
  }
  
  df <- data.frame(
    Output = c(
      "Prev total",
      "Prev hospital",
      "Prev community",
      "CDI incidence total (absolute)",
      "CDI incidence hospital (absolute)",
      "CDI incidence community (absolute)",
      "CDI incidence hospital (relative)",
      "CDI incidence community (relative)",
      "Primo CDI incidence total (absolute)",
      "Primo CDI incidence hospital (absolute)",
      "Primo CDI incidence community (absolute)",
      "Rec CDI incidence total (absolute)",
      "Rec CDI incidence hospital (absolute)",
      "Rec CDI incidence community (absolute)",
      "Recidive I_II / I",
      "Recidive I_III / I_II",
      "R0 hospital",
      "R0 community"
    ),
    Value = c(
      fmt(prev_tot * 100),
      fmt(prev_h * 100),
      fmt(prev_c * 100),
      fmt(inc_tot_abs),
      fmt(inc_h_abs),
      fmt(inc_c_abs),
      fmt(inc_h_rel),
      fmt(inc_c_rel),
      fmt(inc_primo_abs),
      fmt(inc_h_primo_abs),
      fmt(inc_c_primo_abs),
      fmt(inc_rec_abs),
      fmt(inc_h_rec_abs),
      fmt(inc_c_rec_abs),
      fmt(rec1),
      fmt(rec2),
      fmt(R0h),
      fmt(R0c)
    ),
    Units = c(
      "%", "%", "%", "/100k/year", "/100k/year", "/100k/year", "/100k/year", "/100k/year",
      "/100k/year", "/100k/year", "/100k/year", "/100k/year", "/100k/year", "/100k/year",
      "%", "%", "-", "-"
    ),
    stringsAsFactors = FALSE
  )
  
  tab <- make_pretty_table(df)
  if (show) gridExtra::grid.arrange(tab)
  invisible(tab)
}

plot_calibration_table3 <- function(params_calib, metrics_calib, alpha_end = NA_real_, show = TRUE) {
  
  fmt <- function(x) ifelse(is.na(x), "", sprintf("%.4g", x))
  same_val <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12))
  
  # Values
  N_h <- metrics_calib$population$N_h[1]
  N_c <- metrics_calib$population$N_c[1]
  if (is.na(alpha_end) && !is.null(metrics_calib$R0$alpha_end)) {
    alpha_end <- metrics_calib$R0$alpha_end
  }
  
  rows <- data.frame(
    Parameter = c(
      "N (population size)",
      "β (transmission rate)",
      "ν (relative infectiousness I vs C)",
      "σ (baseline progression C→I)",
      "τ (antibiotic exposure rate)",
      "ω (dysbiosis recovery rate)",
      "γ (natural clearance)",
      "ε (CDI resolution)",
      "p (post-CDI protection)",
      "φ (loss of period at risk)",
      "prop (discharge to S_II vs C_II)",
      "δ (discharge rate)",
      "α (admission at equilibrium)"
    ),
    Hospital = c(
      N_h,
      params_calib["beta_h"],
      params_calib["nu_h"],
      params_calib["sigma_h"],
      params_calib["tau_h"],
      params_calib["omega_h"],
      params_calib["gamma_h"],
      params_calib["epsilon_h"],
      params_calib["p_h"],
      params_calib["phi_h"],
      params_calib["prop"],
      params_calib["delta"],
      alpha_end
    ),
    Community = c(
      N_c,
      params_calib["beta_c"],
      params_calib["nu_c"],
      params_calib["sigma_c"],
      params_calib["tau_c"],
      params_calib["omega_c"],
      params_calib["gamma_c"],
      params_calib["epsilon_c"],
      params_calib["p_c"],
      params_calib["phi_c"],
      params_calib["prop"],
      params_calib["delta"],
      alpha_end
    ),
    Units = c(
      "persons",
      "day^-1", "–", "day^-1", "day^-1", "day^-1", "day^-1", "day^-1",
      "–", "day^-1", "–", "day^-1", "day^-1"
    ),
    Stratified = c(
      "no",
      "no",
      "no",
      "yes",
      "no",
      "no",
      "no",
      "no",
      "no",
      "no",
      "no",
      "yes",
      "yes"
    ),
    stringsAsFactors = FALSE
  )
  
  # Format values as strings
  df <- rows
  df$Hospital <- fmt(df$Hospital)
  df$Community <- fmt(df$Community)
  
  tab <- make_pretty_table(df)
  
  # Merge cells when Hospital == Community
  for (i in seq_len(nrow(rows))) {
    if (same_val(rows$Hospital[i], rows$Community[i])) {
      idx_h <- which(tab$layout$t == i + 1 & tab$layout$l == 2)
      idx_c <- which(tab$layout$t == i + 1 & tab$layout$l == 3)
      if (length(idx_h) == 1) {
        tab$layout$r[idx_h] <- 3
        tab$grobs[[idx_h]]$x <- grid::unit(0.5, "npc")
        tab$grobs[[idx_h]]$hjust <- 0.5
      }
      if (length(idx_c) == 1) {
        tab$grobs[[idx_c]] <- grid::nullGrob()
      }
    }
  }
  
  if (show) gridExtra::grid.arrange(tab)
  invisible(tab)
}

plot_calibration_table4 <- function(params_calib, metrics_calib, alpha_end = NA_real_, show = TRUE) {
  
  fmt <- function(x) ifelse(is.na(x), "", sprintf("%.4g", x))
  same_val <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12))
  
  if (is.na(alpha_end) && !is.null(metrics_calib$R0$alpha_end)) {
    alpha_end <- metrics_calib$R0$alpha_end
  }
  
  rows <- data.frame(
    Parameter = c(
      "σ (baseline progression C→I)",
      "k_A (progression factor for CA)",
      "k_II (recurrence factor 1)",
      "k_III (recurrence factor 2+)",
      "σ_A (CA → I)",
      "σ_II (C_II → I_II)",
      "σ_III (C_III → I_III)",
      "δ (discharge S, C)",
      "δ_I (discharge I)",
      "δ_II (discharge I_II)",
      "δ_III (discharge I_III)",
      "α (admission S, C)",
      "α_I (admission I)",
      "α_II (admission I_II)",
      "α_III (admission I_III)"
    ),
    Hospital = c(
      params_calib["sigma_h"],
      params_calib["k_A"],
      params_calib["k_II"],
      params_calib["k_III"],
      params_calib["k_A"] * params_calib["sigma_h"],
      params_calib["k_II"] * params_calib["sigma_h"],
      params_calib["k_III"] * params_calib["sigma_h"],
      params_calib["delta"],
      params_calib["delta_I"],
      params_calib["delta_II"],
      params_calib["delta_III"],
      alpha_end,
      alpha_end * params_calib["w_I"],
      alpha_end * params_calib["w_II"],
      alpha_end * params_calib["w_III"]
    ),
    Community = c(
      params_calib["sigma_c"],
      params_calib["k_A"],
      params_calib["k_II"],
      params_calib["k_III"],
      params_calib["k_A"] * params_calib["sigma_c"],
      params_calib["k_II"] * params_calib["sigma_c"],
      params_calib["k_III"] * params_calib["sigma_c"],
      params_calib["delta"],
      params_calib["delta_I"],
      params_calib["delta_II"],
      params_calib["delta_III"],
      alpha_end,
      alpha_end * params_calib["w_I"],
      alpha_end * params_calib["w_II"],
      alpha_end * params_calib["w_III"]
    ),
    Units = c(
      "day^-1",
      "–",
      "–",
      "–",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1",
      "day^-1"
    ),
    stringsAsFactors = FALSE
  )
  
  df <- rows
  df$Hospital <- fmt(df$Hospital)
  df$Community <- fmt(df$Community)
  
  tab <- make_pretty_table(df)
  
  # Merge cells when Hospital == Community
  for (i in seq_len(nrow(rows))) {
    if (same_val(rows$Hospital[i], rows$Community[i])) {
      idx_h <- which(tab$layout$t == i + 1 & tab$layout$l == 2)
      idx_c <- which(tab$layout$t == i + 1 & tab$layout$l == 3)
      if (length(idx_h) == 1) {
        tab$layout$r[idx_h] <- 3
        tab$grobs[[idx_h]]$x <- grid::unit(0.5, "npc")
        tab$grobs[[idx_h]]$hjust <- 0.5
      }
      if (length(idx_c) == 1) {
        tab$grobs[[idx_c]] <- grid::nullGrob()
      }
    }
  }
  
  if (show) gridExtra::grid.arrange(tab)
  invisible(tab)
}




###############################################################################
# 3. MODEL DYNAMICS PLOTS 
###############################################################################

# Plot model dynamics (hospital and community side by side)
plot_dynamics <- function(ode_result, targets, N_h, N_c) {
  
  # Convert target prevalence to absolute counts (for dashed horizontal line)
  target_C_h <- targets$prevalence_h * N_h
  target_C_c <- targets$prevalence_c * N_c
  
  # ---- Build aggregated time series for Hospital ----
  df_h <- data.frame(
    time = ode_result$time,
    Susceptible = ode_result$S0_h + ode_result$SA_h + ode_result$S_II_h + ode_result$S_III_h,
    Colonized   = ode_result$C0_h + ode_result$CA_h + ode_result$C_II_h + ode_result$C_III_h,
    Infected    = ode_result$I_h  + ode_result$I_II_h + ode_result$I_III_h,
    Recurrences = ode_result$S_II_h + ode_result$C_II_h + ode_result$I_II_h +
      ode_result$S_III_h + ode_result$C_III_h + ode_result$I_III_h
  )
  
  # ---- Build aggregated time series for Community ----
  df_c <- data.frame(
    time = ode_result$time,
    Susceptible = ode_result$S0_c + ode_result$SA_c + ode_result$S_II_c + ode_result$S_III_c,
    Colonized   = ode_result$C0_c + ode_result$CA_c + ode_result$C_II_c + ode_result$C_III_c,
    Infected    = ode_result$I_c  + ode_result$I_II_c + ode_result$I_III_c,
    Recurrences = ode_result$S_II_c + ode_result$C_II_c + ode_result$I_II_c +
      ode_result$S_III_c + ode_result$C_III_c + ode_result$I_III_c
  )
  
  # build on ggplot (hospital or community)
  build_plot <- function(df, target_C, title_txt) {
    
    df_long <- reshape2::melt(df, id.vars = "time")  # long format for ggplot
    
    ggplot(df_long, aes(x = time, y = value, color = variable)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = target_C, linetype = "dashed", color = "red", linewidth = 0.8) +
      theme_bw() +
      labs(
        x = "Time (days)",
        y = "Number of individuals",
        color = "Group",
        title = title_txt
      ) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 9),
        plot.title   = element_text(face = "bold")
      )
  }
  
  # Build the two panels
  p_h <- build_plot(df_h, target_C_h, "Model dynamics - Hospital")
  p_c <- build_plot(df_c, target_C_c, "Model dynamics - Community")
  
  # Combine side by side
  p_both <- p_h + p_c + patchwork::plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  # Return everything (useful if you want to save individual panels)
  return(list(hospital = p_h, community = p_c, both = p_both))
}


# Plot admission rate alpha(t) over time (computed from state + params)
# This is the same logic as in calibration model: alpha = out_hc / den_alpha
plot_alpha_dynamics <- function(ode_result, params_vec) {
  
  alpha_vec <- numeric(nrow(ode_result))  # Allocate alpha vector
  
  # Compute alpha at each time point
  for (i in seq_len(nrow(ode_result))) {
    
    state <- as.list(ode_result[i, ]) # Extract current row as a named list (includes "time", but we do not use it)
    
    # Compute totals (Hospital and Community)
    tot_h <- compute_totals(state$S0_h, state$C0_h, state$SA_h, state$CA_h, state$I_h, state$S_II_h, state$C_II_h, state$I_II_h, state$S_III_h, state$C_III_h, state$I_III_h)
    tot_c <- compute_totals(state$S0_c, state$C0_c, state$SA_c, state$CA_c, state$I_c, state$S_II_c, state$C_II_c, state$I_II_c, state$S_III_c, state$C_III_c, state$I_III_c)
    
    # Hospital -> Community outflows (numerator)
    out_hc <- params_vec["delta"] * (tot_h$S + tot_h$C) + params_vec["delta_I"] * state$I_h + params_vec["delta_II"] * state$I_II_h + params_vec["delta_III"] * state$I_III_h
    
    # Community admission pool (denominator)
    den_alpha <- params_vec["w"] * (tot_c$S + tot_c$C) + params_vec["w_I"] * state$I_c + params_vec["w_II"] * state$I_II_c + params_vec["w_III"] * state$I_III_c
    
    # Alpha (avoid division by zero)
    alpha_vec[i] <- if (den_alpha == 0) 0 else out_hc / den_alpha
  }
  
  # Build data frame for plotting
  alpha_df <- data.frame(time = ode_result$time, alpha = alpha_vec)
  
  # Plot alpha(t)
  p <- ggplot(alpha_df, aes(x = time, y = alpha)) +
    geom_line(linewidth = 0.9) +
    theme_bw() +
    labs(
      x = "Time (days)",
      y = "Admission rate α (day⁻¹)",
      title = "Dynamic admission rate α over time"
    ) +
    theme(
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}
