###############################################################################
############################# PLOTS #######################################
###############################################################################

###############################################################################
# ---- 1: GRID SEARCH PLOTS ----
###############################################################################
plot_grid_search_beta <- function(grid_result, best_guess, targets) {
  hull_idx <- chull(grid_result$grid$portage_h, grid_result$grid$portage_c)
  hull_data <- grid_result$grid[hull_idx, ]
  
  ggplot(grid_result$grid, aes(x = portage_h, y = portage_c)) +
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
}

plot_grid_search_sigma <- function(grid_result, best_guess, targets) {
  hull_idx <- chull(grid_result$grid$incidence_h, grid_result$grid$incidence_c)
  hull_data <- grid_result$grid[hull_idx, ]
  
  ggplot(grid_result$grid, aes(x = incidence_h, y = incidence_c)) +
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
}

plot_grid_search_k <- function(grid_result, best_guess, targets) {
  hull_idx <- chull(grid_result$grid$recid_1, grid_result$grid$recid_2)
  hull_data <- grid_result$grid[hull_idx, ]
  
  ggplot(grid_result$grid, aes(x = recid_1, y = recid_2)) +
    geom_polygon(data = hull_data, fill = "grey80", alpha = 0.5) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_vline(xintercept = targets$recid_1, linetype = "dashed", color = "red") +
    geom_hline(yintercept = targets$recid_2, linetype = "dashed", color = "red") +
    geom_point(data = best_guess, aes(x = recid_1, y = recid_2),
               color = "blue", size = 3) +
    theme_bw() +
    labs(
      x = "Recurrence rate 1",
      y = "Recurrence rate 2+",
      title = "Grid search: k_II, k_III",
      subtitle = "Feasible recurrence space"
    )
}

###############################################################################
# ---- 2 : CALIBRATION ---- 
###############################################################################

# ---- DYNAMICS PLOTS (UNTIL EQUILIBRIUM) ----
plot_dynamics <- function(ode_result, params_vec, targets, N_h, N_c) {
  
  # Add aggregated compartments
  ode_result$C_h <- with(ode_result, C0_h + CA_h + C_II_h + C_III_h)
  ode_result$C_c <- with(ode_result, C0_c + CA_c + C_II_c + C_III_c)
  ode_result$I_tot_h <- with(ode_result, I_h + I_II_h + I_III_h)
  ode_result$I_tot_c <- with(ode_result, I_c + I_II_c + I_III_c)
  ode_result$R_h <- with(ode_result, S_II_h + S_III_h + C_II_h + C_III_h + I_II_h + I_III_h)
  ode_result$R_c <- with(ode_result, S_II_c + S_III_c + C_II_c + C_III_c + I_II_c + I_III_c)
  ode_result$N_h_sim <- with(ode_result, S0_h + SA_h + S_II_h + S_III_h + C_h + I_tot_h + R_h)
  ode_result$N_c_sim <- with(ode_result, S0_c + SA_c + S_II_c + S_III_c + C_c + I_tot_c + R_c)
  
  # Targets in absolute numbers
  target_C_h <- targets$portage_h * N_h
  target_C_c <- targets$portage_c * N_c
  
  # Hospital - all compartments
  p_h_all <- ode_result %>%
    dplyr::select(time, S0_h, SA_h, C0_h, CA_h, I_h,
                  S_II_h, C_II_h, I_II_h, S_III_h, C_III_h, I_III_h) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Time (days)", y = "Number of individuals",
      color = "Hospital:", title = "Model dynamics - Hospital (all compartments)"
    )
  
  # Hospital - aggregated with target and N_h line
  p_h_agg <- ode_result %>%
    dplyr::select(time, S0_h, SA_h, C_h, I_tot_h, R_h, N_h_sim) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = target_C_h, linetype = "dashed", color = "red") +
    ggplot2::annotate(
      "text", x = max(ode_result$time) * 0.8, y = target_C_h * 1.1,
      label = "Target carriage", color = "red", size = 3
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Time (days)", y = "Number of individuals",
      color = "Hospital:", title = "Model dynamics - Hospital (aggregated)"
    )
  
  # Community - all compartments
  p_c_all <- ode_result %>%
    dplyr::select(time, S0_c, SA_c, C0_c, CA_c, I_c,
                  S_II_c, C_II_c, I_II_c, S_III_c, C_III_c, I_III_c) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Time (days)", y = "Number of individuals",
      color = "Community:", title = "Model dynamics - Community (all compartments)"
    )
  
  # Community - aggregated with target and N_c line
  p_c_agg <- ode_result %>%
    dplyr::select(time, S0_c, SA_c, C_c, I_tot_c, R_c, N_c_sim) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot2::ggplot(ggplot2::aes(time, value, color = variable)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = target_C_c, linetype = "dashed", color = "red") +
    ggplot2::annotate(
      "text", x = max(ode_result$time) * 0.8, y = target_C_c * 1.1,
      label = "Target carriage", color = "red", size = 3
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Time (days)", y = "Number of individuals",
      color = "Community:", title = "Model dynamics - Community (aggregated)"
    )
  
  return(list(
    hospital_all = p_h_all,
    hospital_agg = p_h_agg,
    community_all = p_c_all,
    community_agg = p_c_agg
  ))
}

# ---- R0 vs NU PLOT ----
plot_R0_vs_nu <- function(params_vec, init_cond, time_vec, nu_range = c(1, 30)) {
  
  compute_R0_for_nu <- function(nu_val) {
    params_temp <- params_vec
    params_temp["nu"] <- nu_val
    
    tryCatch({
      result <- run_model_to_equilibrium(params_temp, init_cond, time_vec)
      
      last <- get_equilibrium_state(result)
      last <- as.list(last[1, ])  # <<< PATCH: cohérence avec compute_R0_approx
      
      c(
        R0_h = compute_R0_approx(last_state, params_temp, "h"),
        R0_c = compute_R0_approx(last_state, params_temp, "c")
      )
    }, error = function(e) {
      c(R0_h = NA, R0_c = NA)
    })
  }
  
  nu_values <- seq(nu_range[1], nu_range[2], by = 0.5)
  R0_grid <- do.call(rbind, lapply(nu_values, function(nu_val) {
    R0_pair <- compute_R0_for_nu(nu_val)
    data.frame(
      nu = nu_val,
      R0_h = unname(R0_pair["R0_h"]),
      R0_c = unname(R0_pair["R0_c"])
    )
  }))
  
  R0_grid <- R0_grid[complete.cases(R0_grid), ]
  
  candidates <- subset(R0_grid, R0_c > 1.001)
  best_nu <- if (nrow(candidates) > 0) {
    candidates[which.min(candidates$R0_h), ]
  } else {
    NULL
  }
  
  p <- ggplot2::ggplot(R0_grid, ggplot2::aes(x = R0_h, y = R0_c)) +
    ggplot2::geom_point(ggplot2::aes(color = nu), size = 3) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    ggplot2::scale_color_viridis_c(name = "ν") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = expression(R[0]~"Hospital"),
      y = expression(R[0]~"Community"),
      title = "Optimal choice of ν (relative infectiousness I vs C)",
      subtitle = "Constraint: R0 community > 1 | Objective: Minimize R0 hospital"
    )
  
  if (!is.null(best_nu)) {
    p <- p +
      ggplot2::geom_point(
        data = best_nu, ggplot2::aes(x = R0_h, y = R0_c),
        color = "black", fill = "gold", shape = 21, size = 6, stroke = 1.2
      ) +
      ggplot2::geom_label(
        data = best_nu, ggplot2::aes(label = sprintf("Best ν = %.1f", nu)),
        nudge_x = 0.05, nudge_y = 0.05, fontface = "bold",
        fill = "gold", label.size = 0.4
      )
  }
  
  list(plot = p, grid = R0_grid, best_nu = best_nu)
}

# ---- ALPHA DYNAMICS ----

plot_alpha_dynamics <- function(ode_result, params_vec) {
  
  # Calcul d'alpha à chaque pas de temps
  alpha_data <- lapply(seq_len(nrow(ode_result)), function(i) {
    st <- as.list(ode_result[i, ])
    
    tot_h <- compute_totals(
      st$S0_h, st$SA_h, st$S_II_h, st$S_III_h,
      st$C0_h, st$CA_h, st$C_II_h, st$C_III_h,
      st$I_h, st$I_II_h, st$I_III_h
    )
    
    tot_c <- compute_totals(
      st$S0_c, st$SA_c, st$S_II_c, st$S_III_c,
      st$C0_c, st$CA_c, st$C_II_c, st$C_III_c,
      st$I_c, st$I_II_c, st$I_III_c
    )
    
    out_hc <- params_vec["delta"]   * (tot_h$S + tot_h$C) +
      params_vec["delta_I"] * st$I_h +
      params_vec["delta_II"]* st$I_II_h +
      params_vec["delta_III"]*st$I_III_h
    
    den_alpha <- params_vec["w"]   * (tot_c$S + tot_c$C) +
      params_vec["w_I"] * st$I_c +
      params_vec["w_II"]* st$I_II_c +
      params_vec["w_III"]*st$I_III_c
    
    alpha <- out_hc / den_alpha
    
    data.frame(time = st$time, alpha = alpha)
  })
  
  alpha_df <- do.call(rbind, alpha_data)
  
  # Plot
  p <- ggplot2::ggplot(alpha_df, ggplot2::aes(x = time, y = alpha)) +
    ggplot2::geom_line(linewidth = 0.8, color = "steelblue") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Time (days)",
      y = "Admission rate α (day⁻¹)",
      title = "Dynamic admission rate (α) over time",
      subtitle = "α adjusts to balance hospital discharge flows"
    )
  
  return(p)
}

# ---- ALL PARAMETERS TABLE ----

create_all_parameters_table <- function(params_final, N_h, N_c, alpha_eq) {
  
  # ---- Table content ----
  table_data <- data.frame(
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
      "δ (discharge rate)",
      "α (admission at equilibrium)"
    ),
    
    Hospital = c(
      sprintf("%d", N_h),
      sprintf("%.6f", params_final["beta_h"]),
      sprintf("%.1f", params_final["nu"]),
      sprintf("%.6f", params_final["sigma_h"]),
      sprintf("%.5f", params_final["tau_h"]),
      sprintf("%.3f", params_final["omega_h"]),
      sprintf("%.3f", params_final["gamma"]),
      sprintf("%.3f", params_final["epsilon"]),
      sprintf("%.2f", params_final["p"]),
      sprintf("%.3f", params_final["phi"]),
      sprintf("%.5f", params_final["delta"]),
      sprintf("%.7f", alpha_eq)
    ),
    
    Community = c(
      sprintf("%d", N_c),
      sprintf("%.6f", params_final["beta_c"]),
      "",
      sprintf("%.6f", params_final["sigma_c"]),
      sprintf("%.5f", params_final["tau_c"]),
      sprintf("%.3f", params_final["omega_c"]),
      "",
      "",
      "",
      "",
      "",
      ""
    ),
    
    Strat = c(
      "no", "no", "no",
      "yes",
      "no", "no", "no", "no", "no", "no",
      "yes",
      "yes"
    ),
    
    Unit = c(
      "individuals",
      "day⁻¹", "-",
      "day⁻¹",
      "day⁻¹",
      "day⁻¹",
      "day⁻¹", "day⁻¹", "-", "day⁻¹",
      "day⁻¹",
      "day⁻¹"
    )
  )
  
  # ---- Table theme ----
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
  
  # Create base table
  tg <- tableGrob(table_data, rows = NULL, theme = tt)
  
  # Lignes à fusionner (indices dans le corps du tableau, pas dans table_data)
  # row 3 = ν, rows 7-10 = γ,ε,p,φ, row 11 = δ, row 12 = α
  rows_to_merge <- c(3, 7, 8, 9, 10, 11, 12)
  
  # Fusion : remplace les cellules Hospital et Community par une seule cellule centrée
  for (i in rows_to_merge) {
    row_in_grob <- i + 1  # +1 car la première ligne est le header
    
    # Récupère la valeur (depuis colonne Hospital)
    val <- table_data$Hospital[i]
    
    # Crée un textGrob centré sur 2 colonnes
    merged_cell <- textGrob(
      val,
      gp = gpar(fontsize = 9)
    )
    
    # Remplace les cellules des colonnes 2 (Hospital) et 3 (Community)
    tg$grobs[[which(tg$layout$t == row_in_grob & tg$layout$l == 2)]] <- merged_cell
    tg$grobs[[which(tg$layout$t == row_in_grob & tg$layout$l == 3)]] <- nullGrob()
    
    # Étend la cellule Hospital pour couvrir les 2 colonnes
    tg$layout$l[tg$layout$t == row_in_grob & tg$layout$l == 2] <- 2
    tg$layout$r[tg$layout$t == row_in_grob & tg$layout$r == 2] <- 3
  }
  
  # Arrange avec titre
  grid.arrange(
    textGrob(
      "All Model Parameters",
      gp = gpar(fontsize = 14, fontface = "bold")
    ),
    tg,
    ncol = 1,
    heights = c(0.05, 0.95)
  )
}

# ---- STRATIFIED PARAMETERS ----
create_stratified_parameters_table <- function(params_final, N_h, N_c, alpha_eq) {
  
  # ---- Derived progression rates ----
  sigma_A_h   <- params_final["k_A"]  * params_final["sigma_h"]
  sigma_II_h  <- params_final["k_II"] * params_final["sigma_h"]
  sigma_III_h <- params_final["k_III"] * params_final["sigma_h"]
  
  sigma_A_c   <- params_final["k_A"]  * params_final["sigma_c"]
  sigma_II_c  <- params_final["k_II"] * params_final["sigma_c"]
  sigma_III_c <- params_final["k_III"] * params_final["sigma_c"]
  
  # ---- Derived admission rates ----
  alpha_I   <- alpha_eq * params_final["w_I"]
  alpha_II  <- alpha_eq * params_final["w_II"]
  alpha_III <- alpha_eq * params_final["w_III"]
  
  # ---- Table content ----
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
  
  # ---- Table theme ----
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
  
  # Create base table
  tg <- tableGrob(table_data, rows = NULL, theme = tt)
  
  # Lignes à fusionner : k (2,3,4), delta (8,9,10,11), alpha (12,13,14,15)
  rows_to_merge <- c(2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 15)
  
  for (i in rows_to_merge) {
    row_in_grob <- i + 1
    val <- table_data$Hospital[i]
    
    # Centrer le texte
    merged_cell <- textGrob(val, gp = gpar(fontsize = 9, hjust = 0.5))
    
    tg$grobs[[which(tg$layout$t == row_in_grob & tg$layout$l == 2)]] <- merged_cell
    tg$grobs[[which(tg$layout$t == row_in_grob & tg$layout$l == 3)]] <- nullGrob()
    
    tg$layout$l[tg$layout$t == row_in_grob & tg$layout$l == 2] <- 2
    tg$layout$r[tg$layout$t == row_in_grob & tg$layout$r == 2] <- 3
  }
  
  # Arrange avec titre
  grid.arrange(
    textGrob(
      "Stratified Parameters (Progression & Flow Rates)",
      gp = gpar(fontsize = 14, fontface = "bold")
    ),
    tg,
    ncol = 1,
    heights = c(0.05, 0.95)
  )
}

# ---- CALIBRATED PARAMETERS & TARGETS TABLE ----

create_calibration_table <- function(params_final, metrics, targets) {
  
  table_data <- data.frame(
    Parameter = c(
      "β_h (transmission hospital)",
      "β_c (transmission community)",
      "σ_h (progression C→I hospital)",
      "σ_c (progression C→I community)",
      "k_II (recurrence factor 1)",
      "k_III (recurrence factor 2+)",
      "",
      "Asymptomatic carriage hospital",
      "Asymptomatic carriage community",
      "CDI incidence hospital (HOHA)",
      "CDI incidence community (CO)",
      "Recurrence rate 1",
      "Recurrence rate 2+"
    ),
    
    Value = c(
      sprintf("%.6f", params_final["beta_h"]),
      sprintf("%.6f", params_final["beta_c"]),
      sprintf("%.6f", params_final["sigma_h"]),
      sprintf("%.6f", params_final["sigma_c"]),
      sprintf("%.3f", params_final["k_II"]),
      sprintf("%.3f", params_final["k_III"]),
      "",
      sprintf("%.2f%%", metrics$portage_h * 100),
      sprintf("%.2f%%", metrics$portage_c * 100),
      sprintf("%.2f", metrics$inc_h * 100000 * 365),
      sprintf("%.2f", metrics$inc_c * 100000 * 365),
      sprintf("%.1f%%", metrics$recid_1_tot * 100),
      sprintf("%.1f%%", metrics$recid_2_tot * 100)
    ),
    
    Target = c(
      rep("-", 6),
      "",
      sprintf("%.2f%%", targets$portage_h * 100),
      sprintf("%.2f%%", targets$portage_c * 100),
      sprintf("%.2f", targets$incidence_h * 100000 * 365),
      sprintf("%.2f", targets$incidence_c * 100000 * 365),
      sprintf("%.1f%%", targets$recid_1 * 100),
      sprintf("%.1f%%", targets$recid_2 * 100)
    ),
    
    Relative_Error = c(
      rep("-", 6),
      "",
      sprintf("%.1f%%", metrics$errors$err_portage_h * 100),
      sprintf("%.1f%%", metrics$errors$err_portage_c * 100),
      sprintf("%.1f%%", metrics$errors$err_inc_h * 100),
      sprintf("%.1f%%", metrics$errors$err_inc_c * 100),
      sprintf("%.1f%%", metrics$errors$err_recid_1_tot * 100),
      sprintf("%.1f%%", metrics$errors$err_recid_2_tot * 100)
    )
  )
  
  tt <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0, x = 0.05, fontsize = 9),
      bg_params = list(fill = c(
        rep("honeydew2", 6),
        "white",
        rep("lightyellow", 6)
      ))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10),
      bg_params = list(fill = "steelblue", col = "white")
    )
  )
  
  grid.arrange(
    textGrob("Calibrated Parameters & Targets",
             gp = gpar(fontsize = 14, fontface = "bold")),
    tableGrob(table_data, rows = NULL, theme = tt),
    ncol = 1,
    heights = c(0.05, 0.95)
  )
}

# ---- MODEL OUTPUTS TABLE ----

create_outputs_table <- function(metrics, N_h, N_c) {
  
  table_data <- data.frame(
    Metric = c(
      "Asymptomatic carriage hospital",
      "Asymptomatic carriage community",
      "CDI incidence hospital",
      "CDI incidence community",
      "CDI incidence total",
      "Recurrence rate 1 (average H+C)",
      "Recurrence rate 2+ (average H+C)",
      "R0 hospital (approximation)",
      "R0 community (approximation)"
    ),
    
    Value = c(
      sprintf("%.2f%%", metrics$portage_h * 100),
      sprintf("%.2f%%", metrics$portage_c * 100),
      sprintf("%.2f", metrics$inc_h * 100000 * 365),
      sprintf("%.2f", metrics$inc_c * 100000 * 365),
      sprintf("%.2f", metrics$inc_tot * 100000* 365),
      sprintf("%.1f%%", metrics$recid_1_tot * 100),
      sprintf("%.1f%%", metrics$recid_2_tot * 100),
      sprintf("%.2f", metrics$R0_h),
      sprintf("%.2f", metrics$R0_c)
    ),
    
    Unit = c(
      "%", "%",
      "cases/100k pop/year", 
      "cases/100k pop/year",
      "cases/100k pop/year",
      "%", "%", "-", "-"
    )
  )
  
  tt <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0, x = 0.05, fontsize = 9),
      bg_params = list(fill = rep(c("grey95", "white"), length.out = nrow(table_data)))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10),
      bg_params = list(fill = "steelblue", col = "white")
    )
  )
  
  grid.arrange(
    textGrob("Model Outputs",
             gp = gpar(fontsize = 14, fontface = "bold")),
    tableGrob(table_data, rows = NULL, theme = tt),
    ncol = 1,
    heights = c(0.05, 0.95)
  )
}


# ---- DETAILED INCIDENCE & RECURRENCE TABLE ----

create_detailed_metrics_table <- function(metrics) {
  
  table_data <- data.frame(
    Metric = c(
      "CDI incidence (total)",
      "  - Hospital",
      "  - Community",
      "Primo-CDI incidence",
      "  - Hospital",
      "  - Community",
      "Recurrent CDI incidence",
      "  - Hospital",
      "  - Community",
      "Recurrence prevalence 1 (I_II/I)",
      "Recurrence prevalence 2 (I_III/I_II)"
    ),
    
    Value = c(
      sprintf("%.2f", metrics$inc_tot * 100000 * 365),
      sprintf("%.2f", metrics$inc_h * 100000 * 365),
      sprintf("%.2f", metrics$inc_c * 100000 * 365),
      sprintf("%.2f", metrics$inc_primo_tot * 100000 * 365),
      sprintf("%.2f", metrics$inc_primo_h * 100000 * 365),
      sprintf("%.2f", metrics$inc_primo_c * 100000 * 365),
      sprintf("%.2f", metrics$inc_rec_tot * 100000 * 365),
      sprintf("%.2f", metrics$inc_rec_h * 100000 * 365),
      sprintf("%.2f", metrics$inc_rec_c * 100000 * 365),
      sprintf("%.1f%%", metrics$recid_1_tot * 100),
      sprintf("%.1f%%", metrics$recid_2_tot * 100)
    ),
    
    Unit = c(
      "total pop", "per 100k pop/year", "per 100k pop/year",
      "total pop", "per 100k pop/year", "per 100k pop/year",
      "total pop", "per 100k pop/year", "per 100k pop/year",
      "%", "%"
    )
  )
  
  tt <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0, x = 0.05, fontsize = 9),
      bg_params = list(fill = c(
        "lightyellow", "grey95", "white",      # CDI total
        "lightyellow", "grey95", "white",      # Primo-CDI
        "lightyellow", "grey95", "white",      # Recurrent CDI
        "lightyellow", "lightyellow"           # Recurrence prevalence (2 lignes jaunes consécutives)
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

###############################################################################
# ---- 3 : PLOT ANTIBIOTIC REDUCTION ----
###############################################################################

# ---- CDI REDUCTION ----

plot_antibiotic_reduction_facets <- function(scenario_results) {

  # Extract data for all settings
  data_list <- list()

  for (setting in names(scenario_results$scenarios)) {
    for (reduction in names(scenario_results$scenarios[[setting]])) {

      scenario <- scenario_results$scenarios[[setting]][[reduction]]

      data_list[[length(data_list) + 1]] <- data.frame(
        setting = setting,
        reduction = as.numeric(reduction) * 100,
        inc_reduction_total = -scenario$inc_reduction_total,
        inc_reduction_h = -scenario$inc_reduction_h,
        inc_reduction_c = -scenario$inc_reduction_c
      )
    }
  }

  df <- do.call(rbind, data_list)

  # Create facet labels
  setting_labels <- c(
    hospital = "Hospital only",
    community = "Community only",
    both = "Both settings"
  )

  df$setting <- factor(df$setting, levels = names(setting_labels), labels = setting_labels)

  # Plot - Total population
  p_total <- ggplot(df, aes(x = reduction, y = inc_reduction_total)) +
    geom_line(linewidth = 1.2, color = "#8B4513") +
    geom_point(size = 3, color = "#8B4513") +
    facet_wrap(~setting, nrow = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "lightblue")
    ) +
    labs(
      x = "Antibiotic use reduction (%)",
      y = "Relative change in CDI incidence (%)",
      title = "Impact of antibiotic stewardship on CDI incidence - Total population",
      subtitle = "Reduction in antibiotic prescribing across different settings"
    )

  # Plot - By setting (hospital vs community separately)
  df_long <- reshape2::melt(df,
                            id.vars = c("setting", "reduction"),
                            measure.vars = c("inc_reduction_h", "inc_reduction_c"),
                            variable.name = "population",
                            value.name = "inc_reduction")

  df_long$population <- factor(df_long$population,
                               levels = c("inc_reduction_h", "inc_reduction_c"),
                               labels = c("Hospital", "Community"))

  p_separate <- ggplot(df_long, aes(x = reduction, y = inc_reduction, color = population)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~setting, nrow = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(
      name = "Population",
      values = c("Hospital" = "#E41A1C", "Community" = "#377EB8")
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "lightblue"),
      legend.position = "bottom"
    ) +
    labs(
      x = "Antibiotic use reduction (%)",
      y = "Relative change in CDI incidence (%)",
      title = "Impact of antibiotic stewardship on CDI incidence - By population",
      subtitle = "Hospital and community populations shown separately"
    )

  list(
    total = p_total,
    by_population = p_separate
  )
}

# ---- PORTAGE ASYMPTOMATIQUE REDUCTION ----

plot_antibiotic_reduction_carriage_facets <- function(scenario_results) {
  baseline <- scenario_results$baseline$metrics

  data_list <- list()
  for (setting in names(scenario_results$scenarios)) {
    for (reduction in names(scenario_results$scenarios[[setting]])) {
      sc <- scenario_results$scenarios[[setting]][[reduction]]

      # % change (négatif = baisse) comme tes plots CDI actuels
      car_change_h <- -sc$carriage_reduction_h
      car_change_c <- -sc$carriage_reduction_c
      car_change_total <- -((baseline$portage_tot - sc$metrics$portage_tot) / baseline$portage_tot * 100)

      data_list[[length(data_list) + 1]] <- data.frame(
        setting = setting,
        reduction = as.numeric(reduction) * 100,
        car_change_total = car_change_total,
        car_change_h = car_change_h,
        car_change_c = car_change_c
      )
    }
  }

  df <- do.call(rbind, data_list)

  setting_labels <- c(hospital = "Hospital only", community = "Community only", both = "Both settings")
  df$setting <- factor(df$setting, levels = names(setting_labels), labels = setting_labels)

  # Total
  p_total <- ggplot(df, aes(x = reduction, y = car_change_total)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~setting, nrow = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    theme(strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "lightblue")) +
    labs(
      x = "Antibiotic use reduction (%)",
      y = "Relative change in asymptomatic carriage (%)",
      title = "Impact of antibiotic stewardship on asymptomatic carriage - Total"
    )

  # Hospital vs Community
  df_long <- reshape2::melt(df,
                            id.vars = c("setting", "reduction"),
                            measure.vars = c("car_change_h", "car_change_c"),
                            variable.name = "population",
                            value.name = "car_change")
  df_long$population <- factor(df_long$population,
                               levels = c("car_change_h", "car_change_c"),
                               labels = c("Hospital", "Community"))

  p_separate <- ggplot(df_long, aes(x = reduction, y = car_change, color = population)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_wrap(~setting, nrow = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    theme(strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "lightblue"),
          legend.position = "bottom") +
    labs(
      x = "Antibiotic use reduction (%)",
      y = "Relative change in asymptomatic carriage (%)",
      color = "Population",
      title = "Impact of antibiotic stewardship on asymptomatic carriage - By population"
    )

  list(total = p_total, by_population = p_separate)
}


