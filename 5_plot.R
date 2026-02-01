###############################################################################
############################# 5 : PLOTS #######################################
###############################################################################

###############################################################################
# 1. GRID SEARCH PLOT
###############################################################################

###############################################################################
# GENERIC GRID-SEARCH PLOT (works for beta / sigma / k)
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
  
  # Convex hull (need at least 3 points)
  hull_data <- NULL
  if (nrow(df_ok) >= 3) {
    hull_idx <- chull(df_ok[[x_var]], df_ok[[y_var]])
    hull_data <- df_ok[hull_idx, , drop = FALSE]
  }
  
  # Targets
  tx <- targets[[x_target]]
  ty <- targets[[y_target]]
  
  # Plot
  p <- ggplot(df_ok, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    
    # Grey hull (only if it exists)
    { if (!is.null(hull_data)) geom_polygon(data = hull_data, fill = "grey80", alpha = 0.5) } +
    
    # Cloud of points
    geom_point(alpha = 0.3, size = 0.5) +
    
    # Target dashed lines
    geom_vline(xintercept = tx, linetype = "dashed", color = "red") +
    geom_hline(yintercept = ty, linetype = "dashed", color = "red") +
    
    # Best guess in blue
    geom_point(data = best, aes(x = .data[[x_var]], y = .data[[y_var]]),
               color = "blue", size = 3) +
    
    theme_bw() +
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
  plot_grid_search_generic(
    grid_result = grid_result,
    targets = targets,
    x_var = "incidence_h", y_var = "incidence_c",
    x_target = "incidence_h", y_target = "incidence_c",
    title = "Grid search: sigma_h, sigma_c",
    subtitle = "Feasible incidence space",
    x_lab = "Incidence hospital (per day)",
    y_lab = "Incidence community (per day)"
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
  p_both <- p_h + p_c + patchwork::plot_layout(ncol = 2, guides = "collect")
  
  # Return everything (useful if you want to save individual panels)
  return(list(hospital = p_h, community = p_c, both = p_both))
}


###############################################################################
# 2. ALPHA DYNAMICS PLOT
###############################################################################

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
