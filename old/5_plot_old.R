###############################################################################
############################# 5 : PLOTS #######################################
###############################################################################

# This file contains visualization functions for:
# 1. Grid search results (beta, sigma, k calibration steps)
# 2. Model dynamics and equilibrium verification
# 3. Calibration outputs vs targets

###############################################################################
# 1. MODEL DYNAMICS PLOTS
###############################################################################

# Plot model dynamics until equilibrium (all compartments + aggregated views)
# Shows convergence to equilibrium and comparison with target values
# Arguments:
#   - ode_result: output from run_model_to_equilibrium() (data.frame with time + compartments)
#   - targets: named list with target values (portage_h, portage_c, etc.)
#   - N_h: hospital population size
#   - N_c: community population size
plot_dynamics <- function(ode_result, targets, N_h, N_c) {
  
  # Add aggregated compartments for visualization
  ode_result$C_h <- with(ode_result, C0_h + CA_h + C_II_h + C_III_h)
  ode_result$C_c <- with(ode_result, C0_c + CA_c + C_II_c + C_III_c)
  ode_result$I_tot_h <- with(ode_result, I_h + I_II_h + I_III_h)
  ode_result$I_tot_c <- with(ode_result, I_c + I_II_c + I_III_c)
  ode_result$R_h <- with(ode_result, S_II_h + S_III_h + C_II_h + C_III_h + I_II_h + I_III_h)
  ode_result$R_c <- with(ode_result, S_II_c + S_III_c + C_II_c + C_III_c + I_II_c + I_III_c)
  
  # Convert targets from prevalence to absolute numbers
  target_C_h <- targets$portage_h * N_h
  target_C_c <- targets$portage_c * N_c
  
  # ---- Hospital - all compartments ----
  p_h_all <- ode_result %>%
    select(time, S0_h, SA_h, C0_h, CA_h, I_h,
           S_II_h, C_II_h, I_II_h, S_III_h, C_III_h, I_III_h) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot(aes(time, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    theme_bw() +
    labs(x = "Time (days)", 
         y = "Number of individuals",
         color = "Compartment",
         title = "Model dynamics - Hospital (all compartments)")
  
  # ---- Hospital - aggregated ----
  p_h_agg <- ode_result %>%
    select(time, S0_h, SA_h, C_h, I_tot_h, R_h) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot(aes(time, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = target_C_h, linetype = "dashed", color = "red") +
    annotate("text", 
             x = max(ode_result$time) * 0.8, 
             y = target_C_h * 1.1,
             label = "Target carriage", 
             color = "red", 
             size = 3) +
    theme_bw() +
    labs(x = "Time (days)", 
         y = "Number of individuals",
         color = "Compartment",
         title = "Model dynamics - Hospital (aggregated)")
  
  # ---- Community - all compartments ----
  p_c_all <- ode_result %>%
    select(time, S0_c, SA_c, C0_c, CA_c, I_c,
           S_II_c, C_II_c, I_II_c, S_III_c, C_III_c, I_III_c) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot(aes(time, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    theme_bw() +
    labs(x = "Time (days)", 
         y = "Number of individuals",
         color = "Compartment",
         title = "Model dynamics - Community (all compartments)")
  
  # ---- Community - aggregated ----
  p_c_agg <- ode_result %>%
    select(time, S0_c, SA_c, C_c, I_tot_c, R_c) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot(aes(time, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = target_C_c, linetype = "dashed", color = "red") +
    annotate("text", 
             x = max(ode_result$time) * 0.8, 
             y = target_C_c * 1.1,
             label = "Target carriage", 
             color = "red", 
             size = 3) +
    theme_bw() +
    labs(x = "Time (days)", 
         y = "Number of individuals",
         color = "Compartment",
         title = "Model dynamics - Community (aggregated)")
  
  # Return list of plots
  return(list(
    hospital_all = p_h_all,
    hospital_agg = p_h_agg,
    community_all = p_c_all,
    community_agg = p_c_agg
  ))
}

# Plot dynamic admission rate (alpha) over time during calibration
# Shows how alpha adjusts to balance hospital discharge flows
# Arguments:
#   - ode_result: output from run_model_to_equilibrium() with alpha_mode = "dynamic"
#   - params_vec: named parameter vector used in the model
plot_alpha_dynamics <- function(ode_result, params_vec) {
  
  # Compute alpha at each time point
  alpha_vec <- numeric(nrow(ode_result))
  
  for (i in seq_len(nrow(ode_result))) {
    
    # Extract current state
    state <- as.list(ode_result[i, ])
    
    # Compute totals using compute_totals() from file 2
    tot_h <- compute_totals(state$S0_h, state$C0_h, state$SA_h, state$CA_h, 
                            state$I_h, state$S_II_h, state$C_II_h, state$I_II_h, 
                            state$S_III_h, state$C_III_h, state$I_III_h)
    
    tot_c <- compute_totals(state$S0_c, state$C0_c, state$SA_c, state$CA_c, 
                            state$I_c, state$S_II_c, state$C_II_c, state$I_II_c, 
                            state$S_III_c, state$C_III_c, state$I_III_c)
    
    # Hospital outflows
    out_hc <- params_vec["delta"] * (tot_h$S + tot_h$C) +
      params_vec["delta_I"] * state$I_h +
      params_vec["delta_II"] * state$I_II_h +
      params_vec["delta_III"] * state$I_III_h
    
    # Community admission pool
    den_alpha <- params_vec["w"] * (tot_c$S + tot_c$C) +
      params_vec["w_I"] * state$I_c +
      params_vec["w_II"] * state$I_II_c +
      params_vec["w_III"] * state$I_III_c
    
    # Alpha (avoid division by zero)
    alpha_vec[i] <- if (den_alpha == 0) 0 else out_hc / den_alpha
  }
  
  # Create data frame for plotting
  alpha_df <- data.frame(
    time = ode_result$time,
    alpha = alpha_vec
  )
  
  # Create plot
  p <- ggplot(alpha_df, aes(x = time, y = alpha)) +
    geom_line(linewidth = 0.8, color = "steelblue") +
    theme_bw() +
    labs(x = "Time (days)",
         y = "Admission rate α (day⁻¹)",
         title = "Dynamic admission rate (α) over time",
         subtitle = "α adjusts to balance hospital discharge flows")
  
  return(p)
}

###############################################################################
# 2. GRID SEARCH PLOTS
###############################################################################

# Plot grid search results for beta_h and beta_c calibration
# Shows the feasible space of carriage prevalence values
# Arguments:
#   - grid_result: output from grid_search() function
#   - targets: named list with portage_h and portage_c target values
plot_grid_search_beta <- function(grid_result, targets) {
  
  # Extract best guess from grid result
  best_guess <- grid_result$best_guess
  
  # Compute convex hull for feasible space visualization
  hull_idx <- chull(grid_result$grid$portage_h, grid_result$grid$portage_c)
  hull_data <- grid_result$grid[hull_idx, ]
  
  # Create plot
  p <- ggplot(grid_result$grid, aes(x = portage_h, y = portage_c)) +
    geom_polygon(data = hull_data, fill = "grey80", alpha = 0.5) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_vline(xintercept = targets$portage_h, linetype = "dashed", color = "red") +
    geom_hline(yintercept = targets$portage_c, linetype = "dashed", color = "red") +
    geom_point(data = best_guess, aes(x = portage_h, y = portage_c), 
               color = "blue", size = 3) +
    theme_bw() +
    labs(x = "Carriage prevalence hospital", 
         y = "Carriage prevalence community",
         title = "Grid search: beta_h, beta_c",
         subtitle = "Feasible carriage space")
  
  return(p)
}

# Plot grid search results for sigma_h and sigma_c calibration
# Shows the feasible space of CDI incidence values
# Arguments:
#   - grid_result: output from grid_search() function
#   - targets: named list with incidence_h and incidence_c target values
plot_grid_search_sigma <- function(grid_result, targets) {
  
  # Extract best guess from grid result
  best_guess <- grid_result$best_guess
  
  # Compute convex hull for feasible space visualization
  hull_idx <- chull(grid_result$grid$incidence_h, grid_result$grid$incidence_c)
  hull_data <- grid_result$grid[hull_idx, ]
  
  # Create plot
  p <- ggplot(grid_result$grid, aes(x = incidence_h, y = incidence_c)) +
    geom_polygon(data = hull_data, fill = "grey80", alpha = 0.5) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_vline(xintercept = targets$incidence_h, linetype = "dashed", color = "red") +
    geom_hline(yintercept = targets$incidence_c, linetype = "dashed", color = "red") +
    geom_point(data = best_guess, aes(x = incidence_h, y = incidence_c), 
               color = "blue", size = 3) +
    theme_bw() +
    labs(x = "Incidence hospital (per day)", 
         y = "Incidence community (per day)",
         title = "Grid search: sigma_h, sigma_c",
         subtitle = "Feasible incidence space")
  
  return(p)
}

# Plot grid search results for k_II and k_III calibration
# Shows the feasible space of recurrence rate values
# Arguments:
#   - grid_result: output from grid_search() function
#   - targets: named list with recid_1 and recid_2 target values
plot_grid_search_k <- function(grid_result, targets) {
  
  # Extract best guess from grid result
  best_guess <- grid_result$best_guess
  
  # Compute convex hull for feasible space visualization
  hull_idx <- chull(grid_result$grid$recid_1, grid_result$grid$recid_2)
  hull_data <- grid_result$grid[hull_idx, ]
  
  # Create plot
  p <- ggplot(grid_result$grid, aes(x = recid_1, y = recid_2)) +
    geom_polygon(data = hull_data, fill = "grey80", alpha = 0.5) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_vline(xintercept = targets$recid_1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = targets$recid_2, linetype = "dashed", color = "red") +
    geom_point(data = best_guess, aes(x = recid_1, y = recid_2),
               color = "blue", size = 3) +
    theme_bw() +
    labs(x = "Recurrence rate 1",
         y = "Recurrence rate 2+",
         title = "Grid search: k_II, k_III",
         subtitle = "Feasible recurrence space")
  
  return(p)
}

###############################################################################
# 3. PARAMETER AND OUTPUTS TABLES
###############################################################################

create_outputs_table <- function(last_state, params, N_h, N_c) {
  
  N_h_eq <- get_total_pop(last_state, "h", "both")
  N_c_eq <- get_total_pop(last_state, "c", "both")
  N_tot_eq <- N_h_eq + N_c_eq
  
  # Metrics
  carriage_h <- get_carriage(last_state, "total", "h", "both") / N_h_eq
  carriage_c <- get_carriage(last_state, "total", "c", "both") / N_c_eq
  
  incidence_h <- get_incidence(last_state, params, "total", "h", "both") / N_tot_eq * 100000 * 365
  incidence_c <- get_incidence(last_state, params, "total", "c", "both") / N_tot_eq * 100000 * 365
  
  recurrence_1 <- compute_recurrence_prevalence(last_state, "rec_1", "total", "both")
  recurrence_2 <- compute_recurrence_prevalence(last_state, "rec_2", "total", "both")
  
  dysbiosis_h <- get_dysbiosis(last_state, "h", "both") / N_h_eq
  dysbiosis_c <- get_dysbiosis(last_state, "c", "both") / N_c_eq
  
  out <- data.frame(
    Metric = c("N_h (eq)", "N_c (eq)",
               "Carriage_h", "Carriage_c",
               "Incidence_h", "Incidence_c",
               "Recurrence_1", "Recurrence_2+",
               "Dysbiosis_h", "Dysbiosis_c"),
    Value = c(
      sprintf("%.0f", N_h_eq),
      sprintf("%.0f", N_c_eq),
      sprintf("%.1f%%", carriage_h * 100),
      sprintf("%.1f%%", carriage_c * 100),
      sprintf("%.2f", incidence_h),
      sprintf("%.2f", incidence_c),
      sprintf("%.1f%%", recurrence_1 * 100),
      sprintf("%.1f%%", recurrence_2 * 100),
      sprintf("%.1f%%", dysbiosis_h * 100),
      sprintf("%.1f%%", dysbiosis_c * 100)
    ),
    Unit = c(
      "individuals", "individuals",
      "percent", "percent",
      "cases/100k/year (N_total)", "cases/100k/year (N_total)",
      "percent", "percent",
      "percent", "percent"
    ),
    stringsAsFactors = FALSE
  )
  
  tt <- ttheme_default(
    core = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 9)),
    colhead = list(fg_params = list(fontface = "bold", fontsize = 10),
                   bg_params = list(fill = "steelblue", col = "white"))
  )
  
  grid.arrange(
    textGrob("Model Outputs at Equilibrium",
             gp = gpar(fontsize = 14, fontface = "bold")),
    tableGrob(out, rows = NULL, theme = tt),
    ncol = 1,
    heights = c(0.05, 0.95)
  )
}





create_all_parameters_table <- function(params_final, N_h, N_c, alpha_eq) {
  
  # Build separate data for display
  param_names <- c(
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
    "k_A (progression factor CA)",
    "prop (discharge to S_II vs C_II)",
    "δ (discharge rate)",
    "α (admission at equilibrium)"
  )
  
  hosp_vals <- c(
    sprintf("%d", N_h),
    sprintf("%.6f", params_final["beta_h"]),
    sprintf("%.1f", params_final["nu_h"]),
    sprintf("%.6f", params_final["sigma_h"]),
    sprintf("%.5f", params_final["tau_h"]),
    sprintf("%.3f", params_final["omega_h"]),
    sprintf("%.3f", params_final["gamma_h"]),
    sprintf("%.3f", params_final["epsilon_h"]),
    sprintf("%.2f", params_final["p_h"]),
    sprintf("%.3f", params_final["phi_h"]),
    sprintf("%.1f", params_final["k_A"]),
    sprintf("%.2f", params_final["prop"]),
    sprintf("%.5f", params_final["delta"]),
    sprintf("%.7f", alpha_eq)
  )
  
  comm_vals <- c(
    sprintf("%d", N_c),
    sprintf("%.6f", params_final["beta_c"]),
    sprintf("%.1f", params_final["nu_c"]),
    sprintf("%.6f", params_final["sigma_c"]),
    sprintf("%.5f", params_final["tau_c"]),
    sprintf("%.3f", params_final["omega_c"]),
    sprintf("%.3f", params_final["gamma_c"]),
    sprintf("%.3f", params_final["epsilon_c"]),
    sprintf("%.2f", params_final["p_c"]),
    sprintf("%.3f", params_final["phi_c"]),
    sprintf("%.1f", params_final["k_A"]),
    sprintf("%.2f", params_final["prop"]),
    "",
    ""
  )
  
  units <- c(
    "individuals",
    "day⁻¹", "-", "day⁻¹", "day⁻¹", "day⁻¹",
    "day⁻¹", "day⁻¹", "-", "day⁻¹", "-", "-",
    "day⁻¹", "day⁻¹"
  )
  
  # Create display columns with merged values shown as "Both"
  hosp_display <- hosp_vals
  comm_display <- comm_vals
  
  for (i in seq_along(hosp_vals)) {
    if (hosp_vals[i] == comm_vals[i] && hosp_vals[i] != "") {
      hosp_display[i] <- "—"  # Em dash to indicate merged
      comm_display[i] <- hosp_vals[i]  # Show value only in community column
    }
  }
  
  # Build table data
  table_data <- data.frame(
    Parameter = param_names,
    Hospital = hosp_display,
    Community = comm_display,
    Unit = units,
    stringsAsFactors = FALSE
  )
  
  # Create custom theme with centered alignment for merged cells
  tt <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5, fontsize = 9),
      bg_params = list(fill = rep(c("grey95", "white"), length.out = nrow(table_data)))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10),
      bg_params = list(fill = "steelblue", col = "white")
    )
  )
  
  # Create table grob
  tg <- tableGrob(table_data, rows = NULL, theme = tt)
  
  # Add note about merged cells
  note <- textGrob("— indicates parameter shared between Hospital and Community",
                   x = 0.5, y = 0.5,
                   gp = gpar(fontsize = 8, fontface = "italic", col = "grey40"))
  
  # Add title and arrange
  grid.arrange(
    textGrob("All Model Parameters",
             gp = gpar(fontsize = 14, fontface = "bold")),
    tg,
    note,
    ncol = 1,
    heights = c(0.05, 0.90, 0.05)
  )
}


# Create table of stratified parameters (progression and flow rates)
# Shows derived parameters that vary by disease stage or severity
# Arguments:
#   - params_final: named parameter vector after calibration
#   - alpha_eq: equilibrium admission rate
create_stratified_parameters_table <- function(params_final, alpha_eq) {
  
  # Compute derived progression rates
  sigma_A_h <- params_final["k_A"] * params_final["sigma_h"]
  sigma_II_h <- params_final["k_II"] * params_final["sigma_h"]
  sigma_III_h <- params_final["k_III"] * params_final["sigma_h"]
  
  sigma_A_c <- params_final["k_A"] * params_final["sigma_c"]
  sigma_II_c <- params_final["k_II"] * params_final["sigma_c"]
  sigma_III_c <- params_final["k_III"] * params_final["sigma_c"]
  
  # Compute derived admission rates
  alpha_I <- alpha_eq * params_final["w_I"]
  alpha_II <- alpha_eq * params_final["w_II"]
  alpha_III <- alpha_eq * params_final["w_III"]
  
  # Build table data
  table_data <- data.frame(
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
      sprintf("%.6f", params_final["sigma_h"]),
      sprintf("%.3f", params_final["k_A"]),
      sprintf("%.3f", params_final["k_II"]),
      sprintf("%.3f", params_final["k_III"]),
      sprintf("%.6f", sigma_A_h),
      sprintf("%.6f", sigma_II_h),
      sprintf("%.6f", sigma_III_h),
      sprintf("%.5f", params_final["delta"]),
      sprintf("%.5f", params_final["delta_I"]),
      sprintf("%.5f", params_final["delta_II"]),
      sprintf("%.5f", params_final["delta_III"]),
      sprintf("%.7f", alpha_eq),
      sprintf("%.7f", alpha_I),
      sprintf("%.7f", alpha_II),
      sprintf("%.7f", alpha_III)
    ),
    
    Community = c(
      sprintf("%.6f", params_final["sigma_c"]),
      "",
      "",
      "",
      sprintf("%.6f", sigma_A_c),
      sprintf("%.6f", sigma_II_c),
      sprintf("%.6f", sigma_III_c),
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      ""
    ),
    
    Unit = c(
      "day⁻¹", "-", "-", "-",
      "day⁻¹", "day⁻¹", "day⁻¹",
      "day⁻¹", "day⁻¹", "day⁻¹", "day⁻¹",
      "day⁻¹", "day⁻¹", "day⁻¹", "day⁻¹"
    )
  )
  
  # Create table theme
  tt <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5, fontsize = 9),
      bg_params = list(fill = rep(c("grey95", "white"), length.out = nrow(table_data)))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10),
      bg_params = list(fill = "steelblue", col = "white")
    )
  )
  
  # Create table grob
  tg <- tableGrob(table_data, rows = NULL, theme = tt)
  
  # Add title and arrange
  grid.arrange(
    textGrob("Stratified Parameters (Progression & Flow Rates)",
             gp = gpar(fontsize = 14, fontface = "bold")),
    tg,
    ncol = 1,
    heights = c(0.05, 0.95)
  )
}

# Create table comparing model outputs with target metrics
# Shows calibration quality
# Arguments:
#   - outputs: list with portage_h, portage_c, incidence_h, incidence_c, recid_1, recid_2
#   - targets: list with same structure
#   - params_final: calibrated parameters (for displaying beta, sigma, k values)
create_calibration_comparison_table <- function(outputs, targets, params_final) {
  
  # Compute relative errors
  rel_err <- function(value, target) {
    if (target == 0) return(NA)
    ((value - target) / target) * 100
  }
  
  # Build table data
  table_data <- data.frame(
    Metric = c(
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
      sprintf("%.6f", params_final["beta_h"]),
      sprintf("%.6f", params_final["beta_c"]),
      sprintf("%.6f", params_final["sigma_h"]),
      sprintf("%.6f", params_final["sigma_c"]),
      sprintf("%.3f", params_final["k_II"]),
      sprintf("%.3f", params_final["k_III"]),
      "",
      "",
      sprintf("%.2f%%", outputs$portage_h * 100),
      sprintf("%.2f%%", outputs$portage_c * 100),
      sprintf("%.2f", outputs$incidence_h * 100000 * 365),
      sprintf("%.2f", outputs$incidence_c * 100000 * 365),
      sprintf("%.1f%%", outputs$recid_1 * 100),
      sprintf("%.1f%%", outputs$recid_2 * 100)
    ),
    
    Target = c(
      "", rep("-", 6), "", "",
      sprintf("%.2f%%", targets$portage_h * 100),
      sprintf("%.2f%%", targets$portage_c * 100),
      sprintf("%.2f", targets$incidence_h * 100000 * 365),
      sprintf("%.2f", targets$incidence_c * 100000 * 365),
      sprintf("%.1f%%", targets$recid_1 * 100),
      sprintf("%.1f%%", targets$recid_2 * 100)
    ),
    
    Relative_Error = c(
      "", rep("-", 6), "", "",
      sprintf("%.1f%%", rel_err(outputs$portage_h, targets$portage_h)),
      sprintf("%.1f%%", rel_err(outputs$portage_c, targets$portage_c)),
      sprintf("%.1f%%", rel_err(outputs$incidence_h, targets$incidence_h)),
      sprintf("%.1f%%", rel_err(outputs$incidence_c, targets$incidence_c)),
      sprintf("%.1f%%", rel_err(outputs$recid_1, targets$recid_1)),
      sprintf("%.1f%%", rel_err(outputs$recid_2, targets$recid_2))
    )
  )
  
  # Create table theme
  tt <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0, x = 0.05, fontsize = 9),
      bg_params = list(fill = c(
        "lightyellow",        # header
        rep("honeydew2", 6),  # parameters
        "white",              # blank
        "lightyellow",        # outputs header
        rep("aliceblue", 6)   # outputs
      ))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10),
      bg_params = list(fill = "steelblue", col = "white")
    )
  )
  
  # Add title and arrange
  grid.arrange(
    textGrob("Calibration Summary",
             gp = gpar(fontsize = 14, fontface = "bold")),
    tableGrob(table_data, rows = NULL, theme = tt),
    ncol = 1,
    heights = c(0.05, 0.95)
  )
}


# Create detailed table of incidence and recurrence metrics
# Shows stratified metrics by primary/recurrent and hospital/community
# Arguments:
#   - last_state: final equilibrium state from ODE
#   - params: parameter vector
create_detailed_metrics_table <- function(last_state, params) {
  
  N_h_eq <- get_total_pop(last_state, "h", "both")
  N_c_eq <- get_total_pop(last_state, "c", "both")
  N_tot  <- N_h_eq + N_c_eq
  
  to_100k_year <- function(x) x * 100000 * 365
  
  # Absolute incidences (per day)
  inc_h_abs <- get_incidence(last_state, params, "total", "h", "both")
  inc_c_abs <- get_incidence(last_state, params, "total", "c", "both")
  
  inc_primo_h_abs <- get_incidence(last_state, params, "primo", "h", "both")
  inc_primo_c_abs <- get_incidence(last_state, params, "primo", "c", "both")
  
  inc_rec_h_abs <- get_incidence(last_state, params, "rec", "h", "both")
  inc_rec_c_abs <- get_incidence(last_state, params, "rec", "c", "both")
  
  # ---- Per TOTAL population (/N_total) ----
  inc_tot     <- (inc_h_abs + inc_c_abs) / N_tot
  inc_h_per_N <- inc_h_abs / N_tot
  inc_c_per_N <- inc_c_abs / N_tot
  
  inc_primo_tot     <- (inc_primo_h_abs + inc_primo_c_abs) / N_tot
  inc_primo_h_per_N <- inc_primo_h_abs / N_tot
  inc_primo_c_per_N <- inc_primo_c_abs / N_tot
  
  inc_rec_tot     <- (inc_rec_h_abs + inc_rec_c_abs) / N_tot
  inc_rec_h_per_N <- inc_rec_h_abs / N_tot
  inc_rec_c_per_N <- inc_rec_c_abs / N_tot
  
  # ---- Relative incidences (/N_h and /N_c) ----
  inc_h_per_Nh <- inc_h_abs / N_h_eq
  inc_c_per_Nc <- inc_c_abs / N_c_eq
  
  inc_primo_h_per_Nh <- inc_primo_h_abs / N_h_eq
  inc_primo_c_per_Nc <- inc_primo_c_abs / N_c_eq
  
  inc_rec_h_per_Nh <- inc_rec_h_abs / N_h_eq
  inc_rec_c_per_Nc <- inc_rec_c_abs / N_c_eq
  
  recid_1 <- compute_recurrence_prevalence(last_state, "rec_1", "total", "both")
  recid_2 <- compute_recurrence_prevalence(last_state, "rec_2", "total", "both")
  
  table_data <- data.frame(
    Metric = c(
      "CDI incidence (total, /N_total)",
      "  - Hospital contribution (/N_total)",
      "  - Community contribution (/N_total)",
      "  - Hospital incidence (/N_h)",
      "  - Community incidence (/N_c)",
      
      "Primary CDI incidence (/N_total)",
      "  - Hospital contribution (/N_total)",
      "  - Community contribution (/N_total)",
      "  - Hospital incidence (/N_h)",
      "  - Community incidence (/N_c)",
      
      "Recurrent CDI incidence (/N_total)",
      "  - Hospital contribution (/N_total)",
      "  - Community contribution (/N_total)",
      "  - Hospital incidence (/N_h)",
      "  - Community incidence (/N_c)",
      
      "Recurrence rate 1 (I_II / I)",
      "Recurrence rate 2+ (I_III / I_II)"
    ),
    Value = c(
      sprintf("%.2f", to_100k_year(inc_tot)),
      sprintf("%.2f", to_100k_year(inc_h_per_N)),
      sprintf("%.2f", to_100k_year(inc_c_per_N)),
      sprintf("%.2f", to_100k_year(inc_h_per_Nh)),
      sprintf("%.2f", to_100k_year(inc_c_per_Nc)),
      
      sprintf("%.2f", to_100k_year(inc_primo_tot)),
      sprintf("%.2f", to_100k_year(inc_primo_h_per_N)),
      sprintf("%.2f", to_100k_year(inc_primo_c_per_N)),
      sprintf("%.2f", to_100k_year(inc_primo_h_per_Nh)),
      sprintf("%.2f", to_100k_year(inc_primo_c_per_Nc)),
      
      sprintf("%.2f", to_100k_year(inc_rec_tot)),
      sprintf("%.2f", to_100k_year(inc_rec_h_per_N)),
      sprintf("%.2f", to_100k_year(inc_rec_c_per_N)),
      sprintf("%.2f", to_100k_year(inc_rec_h_per_Nh)),
      sprintf("%.2f", to_100k_year(inc_rec_c_per_Nc)),
      
      sprintf("%.1f%%", recid_1 * 100),
      sprintf("%.1f%%", recid_2 * 100)
    ),
    Unit = c(
      rep("cases/100k pop/year", 15),
      "%", "%"
    ),
    stringsAsFactors = FALSE
  )
  
  tt <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0, x = 0.05, fontsize = 9),
      bg_params = list(fill = c(
        # Total block (5 rows)
        "lightyellow", "grey95", "white", "grey95", "white",
        # Primary block (5 rows)
        "lightyellow", "grey95", "white", "grey95", "white",
        # Recurrent block (5 rows)
        "lightyellow", "grey95", "white", "grey95", "white",
        # Recurrence rates (2 rows)
        "lightyellow", "lightyellow"
      ))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10),
      bg_params = list(fill = "steelblue", col = "white")
    )
  )
  
  grid.arrange(
    textGrob("Detailed Incidence & Recurrence Metrics",
             gp = gpar(fontsize = 14, fontface = "bold")),
    tableGrob(table_data, rows = NULL, theme = tt),
    ncol = 1,
    heights = c(0.05, 0.95)
  )
}



################################################################################
################################################################################



