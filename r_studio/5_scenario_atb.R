###############################################################################
#################### 5 : ATB REDUCTION SCENARIOS ##############################
###############################################################################

# 1) Compute incidence + carriage time series from an ODE output
compute_atb_timeseries_metrics <- function(ode_df, params_vec) {
  
  needed <- c(
    "time",
    "C0_h","CA_h","C_II_h","C_III_h","C0_c","CA_c","C_II_c","C_III_c",
    "S0_h","SA_h","S_II_h","S_III_h","I_h","I_II_h","I_III_h",
    "S0_c","SA_c","S_II_c","S_III_c","I_c","I_II_c","I_III_c"
  )
  
  miss <- setdiff(needed, colnames(ode_df))
  if (length(miss) > 0) {
    stop(
      "Colonnes manquantes dans la sortie ODE: ",
      paste(miss, collapse=", "),
      "\nVérifie les noms de compartiments dans cdiff_micro()."
    )
  }
  
  out_mat <- vapply(seq_len(nrow(ode_df)), function(i) {
    
    st <- as.list(ode_df[i, , drop = FALSE])
    
    Nh <- compute_population_totals(st, "h")
    Nc <- compute_population_totals(st, "c")
    N_tot <- Nh + Nc  # <-- CORRECTION: population totale
    
    # Incidence par jour, puis × 365 × 100k
    inc_h_primo <- compute_CDI_incidence(st, params_vec, "h", "primo") / N_tot * 365 * 1e5
    inc_h_rec   <- compute_CDI_incidence(st, params_vec, "h", "recurrent") / N_tot * 365 * 1e5
    
    inc_c_primo <- compute_CDI_incidence(st, params_vec, "c", "primo") / N_tot * 365 * 1e5
    inc_c_rec   <- compute_CDI_incidence(st, params_vec, "c", "recurrent") / N_tot * 365 * 1e5
    
    car_h   <- compute_carriage_prevalence(st, "h")
    car_c   <- compute_carriage_prevalence(st, "c")
    car_tot <- compute_carriage_prevalence(st, "both")
    
    c(
      inc_h_total_100k_year = inc_h_primo + inc_h_rec,
      inc_c_total_100k_year = inc_c_primo + inc_c_rec,
      inc_h_primo_100k_year = inc_h_primo,
      inc_h_rec_100k_year   = inc_h_rec,
      inc_c_primo_100k_year = inc_c_primo,
      inc_c_rec_100k_year   = inc_c_rec,
      car_h_pct = car_h * 100,
      car_c_pct = car_c * 100,
      car_total_pct = car_tot * 100
    )
    
  }, numeric(9))
  
  out <- as.data.frame(t(out_mat))
  out$time <- ode_df$time
  
  out
}

# 2) Run baseline + scenarios after equilibrium (step change at t=0)
run_atb_timecourse_scenarios <- function(params_calibrated, init_cond, time_to_eq,
                                         horizon_days = 5*365,
                                         scenarios = list(
                                           super  = list(red_c = 0.25, red_h = 0.10),
                                           moyen  = list(red_c = 0.20, red_h = 0.05),
                                           faible = list(red_c = 0.20, red_h = 0.00)
                                         )) {
  
  # Étape 1 : amener le système à l’équilibre avec les paramètres calibrés
  eq_df <- run_model_to_equilibrium(params_calibrated, init_cond, time_to_eq)
  
  # Dernier état de l’équilibre, utilisé comme condition initiale des projections
  y0 <- as.numeric(eq_df[nrow(eq_df), -1])
  names(y0) <- colnames(eq_df)[-1]
  
  # Grille temporelle pour la projection post-intervention
  times <- seq(0, horizon_days, by = 1)
  
  # ------------------
  # Baseline (aucune réduction ATB)
  # ------------------
  ode_base <- as.data.frame(
    deSolve::lsoda(y = y0, times = times, func = cdiff_micro, parms = params_calibrated)
  )
  
  # Calcul des métriques temporelles associées
  met_base <- compute_atb_timeseries_metrics(ode_base, params_calibrated)
  met_base$scenario <- "baseline"
  
  # Initialisation de la liste de résultats
  out_list <- list(
    baseline = list(
      params  = params_calibrated,
      ode     = ode_base,
      metrics = met_base
    )
  )
  
  # ------------------
  # Scénarios de réduction ATB
  # ------------------
  for (sc_name in names(scenarios)) {
    
    sc <- scenarios[[sc_name]]
    
    # Copie des paramètres calibrés
    p_sc <- params_calibrated
    
    # Application de la réduction d’exposition aux antibiotiques
    p_sc[["tau_c"]] <- params_calibrated[["tau_c"]] * (1 - sc$red_c)
    p_sc[["tau_h"]] <- params_calibrated[["tau_h"]] * (1 - sc$red_h)
    
    # Simulation ODE avec les nouveaux paramètres
    ode_sc <- as.data.frame(
      deSolve::lsoda(y = y0, times = times, func = cdiff_micro, parms = p_sc)
    )
    
    # Calcul des métriques temporelles pour ce scénario
    met_sc <- compute_atb_timeseries_metrics(ode_sc, p_sc)
    met_sc$scenario <- sc_name
    
    # Stockage des résultats
    out_list[[sc_name]] <- list(
      params  = p_sc,
      ode     = ode_sc,
      metrics = met_sc
    )
  }
  
  # Retourne l’ensemble des scénarios (baseline + interventions)
  out_list
}
















# 3) Plot incidence timecourses + annual comparison
plot_atb_timecourses_incidence <- function(tc_results, year_days = 365) {
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  
  df_all <- bind_rows(lapply(tc_results, function(x) x$metrics))
  
  # Force baseline to be flat (FIRST value)
  baseline_first <- df_all %>% 
    filter(scenario == "baseline", time == min(time)) %>%
    select(-time)
  
  df_all <- df_all %>%
    filter(scenario != "baseline") %>%
    bind_rows(
      baseline_first %>%
        slice(rep(1, length(unique(df_all$time)))) %>%
        mutate(time = sort(unique(df_all$time)), scenario = "baseline")
    )
  
  df_all <- df_all %>%
    mutate(
      scenario = recode(scenario, super = "fort"),
      scenario = factor(scenario, levels = c("baseline","faible","moyen","fort"))
    )
  
  pal <- c(
    "faible"   = "#C77DFF",
    "moyen"    = "#60A5FA",
    "fort"     = "#6EE7B7",
    "baseline" = "#808080"
  )
  
  Tmax <- max(df_all$time, na.rm = TRUE)
  
  # ---------- TIMECOURSES (Hospital + Community separate) ----------
  inc_long <- df_all %>%
    pivot_longer(
      cols = c(inc_h_total_100k_year, inc_c_total_100k_year),
      names_to = "population",
      values_to = "incidence"
    ) %>%
    mutate(
      population = factor(population, 
                          levels = c("inc_h_total_100k_year", "inc_c_total_100k_year"),
                          labels = c("Hospital", "Community"))
    )
  
  p_ts <- ggplot(inc_long, aes(x = time, y = incidence, color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1) +
    facet_wrap(~population, nrow = 1, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = c(baseline="dashed", fort="solid", moyen="solid", faible="solid")) +
    labs(
      x = "Time (days)", 
      y = "CDI incidence per 100k population / year",
      title = "CDI incidence after an antibiotic-use reduction (step change at t=0)"
    )
  
  # ---------- ANNUAL COMPARE (cumulative, with Total + Primo + Rec) ----------
  inc_last_year <- df_all %>%
    filter(time >= (Tmax - (year_days - 1))) %>%
    group_by(scenario) %>%
    summarise(
      total_h = sum(inc_h_total_100k_year, na.rm = TRUE),
      primo_h = sum(inc_h_primo_100k_year, na.rm = TRUE),
      rec_h   = sum(inc_h_rec_100k_year, na.rm = TRUE),
      total_c = sum(inc_c_total_100k_year, na.rm = TRUE),
      primo_c = sum(inc_c_primo_100k_year, na.rm = TRUE),
      rec_c   = sum(inc_c_rec_100k_year, na.rm = TRUE),
      .groups = "drop"
    )
  
  base_vals <- inc_last_year %>% filter(scenario == "baseline")
  
  # Hospital
  inc_long_h <- inc_last_year %>%
    select(scenario, total_h, primo_h, rec_h) %>%
    pivot_longer(-scenario, names_to = "type", values_to = "inc") %>%
    mutate(
      type = factor(type, levels = c("total_h", "primo_h", "rec_h"),
                    labels = c("Total CDI", "Primo CDI", "Recurrent CDI"))
    )
  
  base_h <- base_vals %>%
    select(total_h, primo_h, rec_h) %>%
    pivot_longer(everything(), names_to = "type", values_to = "base_inc") %>%
    mutate(type = factor(type, levels = c("total_h", "primo_h", "rec_h"),
                         labels = c("Total CDI", "Primo CDI", "Recurrent CDI")))
  
  inc_long_h <- inc_long_h %>%
    left_join(base_h, by = "type") %>%
    mutate(pct_change = (inc - base_inc) / base_inc * 100)
  
  p_annual_h <- ggplot(inc_long_h, aes(x = type, y = inc, fill = scenario)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = ifelse(scenario == "baseline", "", 
                                 sprintf("%.1f%%", pct_change))),
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    labs(
      x = NULL,
      y = "CDI incidence (per 100k pop/year)",
      fill = "Scenario",
      title = "Hospital"
    ) +
    theme(legend.position = "none")
  
  # Community
  inc_long_c <- inc_last_year %>%
    select(scenario, total_c, primo_c, rec_c) %>%
    pivot_longer(-scenario, names_to = "type", values_to = "inc") %>%
    mutate(
      type = factor(type, levels = c("total_c", "primo_c", "rec_c"),
                    labels = c("Total CDI", "Primo CDI", "Recurrent CDI"))
    )
  
  base_c <- base_vals %>%
    select(total_c, primo_c, rec_c) %>%
    pivot_longer(everything(), names_to = "type", values_to = "base_inc") %>%
    mutate(type = factor(type, levels = c("total_c", "primo_c", "rec_c"),
                         labels = c("Total CDI", "Primo CDI", "Recurrent CDI")))
  
  inc_long_c <- inc_long_c %>%
    left_join(base_c, by = "type") %>%
    mutate(pct_change = (inc - base_inc) / base_inc * 100)
  
  p_annual_c <- ggplot(inc_long_c, aes(x = type, y = inc, fill = scenario)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = ifelse(scenario == "baseline", "", 
                                 sprintf("%.1f%%", pct_change))),
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    labs(
      x = NULL,
      y = NULL,  # Pas d'axe Y à droite
      fill = "Scenario",
      title = "Community"
    )
  
  p_annual_combined <- p_annual_h + p_annual_c +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(title = "Annual CDI incidence by type (last year)")
  
  # ---------- INSTANT COMPARE (instantaneous at end) ----------
  inst_end <- df_all %>%
    filter(time == Tmax) %>%
    select(scenario, 
           inc_h_total_100k_year, inc_h_primo_100k_year, inc_h_rec_100k_year,
           inc_c_total_100k_year, inc_c_primo_100k_year, inc_c_rec_100k_year)
  
  base_inst <- inst_end %>% filter(scenario == "baseline")
  
  # Hospital
  inst_long_h <- inst_end %>%
    select(scenario, inc_h_total_100k_year, inc_h_primo_100k_year, inc_h_rec_100k_year) %>%
    pivot_longer(-scenario, names_to = "type", values_to = "inc") %>%
    mutate(
      type = factor(type, 
                    levels = c("inc_h_total_100k_year", "inc_h_primo_100k_year", "inc_h_rec_100k_year"),
                    labels = c("Total CDI", "Primo CDI", "Recurrent CDI"))
    )
  
  base_inst_h <- base_inst %>%
    select(inc_h_total_100k_year, inc_h_primo_100k_year, inc_h_rec_100k_year) %>%
    pivot_longer(everything(), names_to = "type", values_to = "base_inc") %>%
    mutate(type = factor(type,
                         levels = c("inc_h_total_100k_year", "inc_h_primo_100k_year", "inc_h_rec_100k_year"),
                         labels = c("Total CDI", "Primo CDI", "Recurrent CDI")))
  
  inst_long_h <- inst_long_h %>%
    left_join(base_inst_h, by = "type") %>%
    mutate(pct_change = (inc - base_inc) / base_inc * 100)
  
  p_instant_h <- ggplot(inst_long_h, aes(x = type, y = inc, fill = scenario)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = ifelse(scenario == "baseline", "", 
                                 sprintf("%.1f%%", pct_change))),
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    labs(
      x = NULL,
      y = "CDI incidence (per 100k pop/year)",
      fill = "Scenario",
      title = "Hospital"
    ) +
    theme(legend.position = "none")
  
  # Community
  inst_long_c <- inst_end %>%
    select(scenario, inc_c_total_100k_year, inc_c_primo_100k_year, inc_c_rec_100k_year) %>%
    pivot_longer(-scenario, names_to = "type", values_to = "inc") %>%
    mutate(
      type = factor(type,
                    levels = c("inc_c_total_100k_year", "inc_c_primo_100k_year", "inc_c_rec_100k_year"),
                    labels = c("Total CDI", "Primo CDI", "Recurrent CDI"))
    )
  
  base_inst_c <- base_inst %>%
    select(inc_c_total_100k_year, inc_c_primo_100k_year, inc_c_rec_100k_year) %>%
    pivot_longer(everything(), names_to = "type", values_to = "base_inc") %>%
    mutate(type = factor(type,
                         levels = c("inc_c_total_100k_year", "inc_c_primo_100k_year", "inc_c_rec_100k_year"),
                         labels = c("Total CDI", "Primo CDI", "Recurrent CDI")))
  
  inst_long_c <- inst_long_c %>%
    left_join(base_inst_c, by = "type") %>%
    mutate(pct_change = (inc - base_inc) / base_inc * 100)
  
  p_instant_c <- ggplot(inst_long_c, aes(x = type, y = inc, fill = scenario)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = ifelse(scenario == "baseline", "", 
                                 sprintf("%.1f%%", pct_change))),
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    labs(
      x = NULL,
      y = NULL,
      fill = "Scenario",
      title = "Community"
    )
  
  p_instant_combined <- p_instant_h + p_instant_c +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(title = "Instantaneous CDI incidence by type (end of simulation)")
  
  list(
    timecourses = p_ts,
    annual_compare = p_annual_combined,
    instant_compare = p_instant_combined,
    data = df_all
  )
}

# 4) Plot carriage timecourses + annual comparison
plot_atb_timecourses_carriage <- function(tc_results) {
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  
  df_all <- bind_rows(lapply(tc_results, function(x) x$metrics))
  
  # On a besoin aussi des compartiments individuels
  ode_all <- bind_rows(lapply(names(tc_results), function(sc_name) {
    tc_results[[sc_name]]$ode %>% mutate(scenario = sc_name)
  }))
  
  # Force baseline flat (FIRST value)
  baseline_first <- df_all %>% 
    filter(scenario == "baseline", time == min(time)) %>%
    select(-time)
  
  df_all <- df_all %>%
    filter(scenario != "baseline") %>%
    bind_rows(
      baseline_first %>%
        slice(rep(1, length(unique(df_all$time)))) %>%
        mutate(time = sort(unique(df_all$time)), scenario = "baseline")
    )
  
  df_all <- df_all %>%
    mutate(
      scenario = recode(scenario, super = "fort"),
      scenario = factor(scenario, levels = c("baseline","faible","moyen","fort"))
    )
  
  ode_all <- ode_all %>%
    mutate(
      scenario = recode(scenario, super = "fort"),
      scenario = factor(scenario, levels = c("baseline","faible","moyen","fort"))
    )
  
  pal <- c(
    "faible"   = "#C77DFF",
    "moyen"    = "#60A5FA",
    "fort"     = "#6EE7B7",
    "baseline" = "#808080"
  )
  
  # ---------- TIMECOURSES ----------
  car_long <- df_all %>%
    pivot_longer(
      cols = c(car_h_pct, car_c_pct),
      names_to = "population",
      values_to = "carriage_pct"
    ) %>%
    mutate(
      population = factor(population, levels = c("car_h_pct","car_c_pct"),
                          labels = c("Hospital","Community"))
    )
  
  p_ts <- ggplot(car_long, aes(x = time, y = carriage_pct, color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1) +
    facet_wrap(~population, nrow = 1, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = c(baseline="dashed", fort="solid", moyen="solid", faible="solid")) +
    labs(
      x = "Time (days)",
      y = "Asymptomatic carriage prevalence (%)",
      title = "Asymptomatic carriage over time after an antibiotic-use reduction (step change at t=0)"
    )
  
  # ---------- ANNUAL COMPARE (Total + Primo + Recurrent carriage) ----------
  Tmax <- max(ode_all$time, na.rm = TRUE)
  
  car_end <- ode_all %>%
    filter(time == Tmax) %>%
    mutate(
      # Total carriage
      car_total_h = (C0_h + CA_h + C_II_h + C_III_h) / (S0_h + SA_h + S_II_h + S_III_h + C0_h + CA_h + C_II_h + C_III_h + I_h + I_II_h + I_III_h) * 100,
      car_total_c = (C0_c + CA_c + C_II_c + C_III_c) / (S0_c + SA_c + S_II_c + S_III_c + C0_c + CA_c + C_II_c + C_III_c + I_c + I_II_c + I_III_c) * 100,
      # Primo carriage (C0 + CA)
      car_primo_h = (C0_h + CA_h) / (S0_h + SA_h + S_II_h + S_III_h + C0_h + CA_h + C_II_h + C_III_h + I_h + I_II_h + I_III_h) * 100,
      car_primo_c = (C0_c + CA_c) / (S0_c + SA_c + S_II_c + S_III_c + C0_c + CA_c + C_II_c + C_III_c + I_c + I_II_c + I_III_c) * 100,
      # Recurrent carriage (C_II + C_III)
      car_rec_h = (C_II_h + C_III_h) / (S0_h + SA_h + S_II_h + S_III_h + C0_h + CA_h + C_II_h + C_III_h + I_h + I_II_h + I_III_h) * 100,
      car_rec_c = (C_II_c + C_III_c) / (S0_c + SA_c + S_II_c + S_III_c + C0_c + CA_c + C_II_c + C_III_c + I_c + I_II_c + I_III_c) * 100
    ) %>%
    select(scenario, car_total_h, car_primo_h, car_rec_h, car_total_c, car_primo_c, car_rec_c)
  
  base_car <- car_end %>% filter(scenario == "baseline")
  
  # Hospital
  car_long_h <- car_end %>%
    select(scenario, car_total_h, car_primo_h, car_rec_h) %>%
    pivot_longer(-scenario, names_to = "type", values_to = "car_pct") %>%
    mutate(
      type = factor(type, levels = c("car_total_h", "car_primo_h", "car_rec_h"),
                    labels = c("Total carriage", "Primo carriage", "Recurrent carriage"))
    )
  
  base_car_h <- base_car %>%
    select(car_total_h, car_primo_h, car_rec_h) %>%
    pivot_longer(everything(), names_to = "type", values_to = "base_car") %>%
    mutate(type = factor(type, levels = c("car_total_h", "car_primo_h", "car_rec_h"),
                         labels = c("Total carriage", "Primo carriage", "Recurrent carriage")))
  
  car_long_h <- car_long_h %>%
    left_join(base_car_h, by = "type") %>%
    mutate(pct_change = (car_pct - base_car) / base_car * 100)
  
  p_car_h <- ggplot(car_long_h, aes(x = type, y = car_pct, fill = scenario)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = ifelse(scenario == "baseline", "", 
                                 sprintf("%.1f%%", pct_change))),
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Asymptomatic carriage prevalence (%)",
      fill = "Scenario",
      title = "Hospital"
    ) +
    theme(legend.position = "none")
  
  # Community
  car_long_c <- car_end %>%
    select(scenario, car_total_c, car_primo_c, car_rec_c) %>%
    pivot_longer(-scenario, names_to = "type", values_to = "car_pct") %>%
    mutate(
      type = factor(type, levels = c("car_total_c", "car_primo_c", "car_rec_c"),
                    labels = c("Total carriage", "Primo carriage", "Recurrent carriage"))
    )
  
  base_car_c <- base_car %>%
    select(car_total_c, car_primo_c, car_rec_c) %>%
    pivot_longer(everything(), names_to = "type", values_to = "base_car") %>%
    mutate(type = factor(type, levels = c("car_total_c", "car_primo_c", "car_rec_c"),
                         labels = c("Total carriage", "Primo carriage", "Recurrent carriage")))
  
  car_long_c <- car_long_c %>%
    left_join(base_car_c, by = "type") %>%
    mutate(pct_change = (car_pct - base_car) / base_car * 100)
  
  p_car_c <- ggplot(car_long_c, aes(x = type, y = car_pct, fill = scenario)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = ifelse(scenario == "baseline", "", 
                                 sprintf("%.1f%%", pct_change))),
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    labs(
      x = NULL,
      y = NULL,
      fill = "Scenario",
      title = "Community"
    )
  
  p_car_combined <- p_car_h + p_car_c +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(title = "Asymptomatic carriage prevalence by type (end of simulation)")
  
  list(
    timecourses = p_ts,
    annual_compare = p_car_combined,
    data = df_all
  )
}


























# 5) Créer un tableau récapitulatif des scénarios sous forme de plot (gridExtra)
create_atb_scenarios_summary_table <- function(tc_results, p_inc_results, p_car_results, year_days = 365) {
  library(dplyr); library(gridExtra); library(grid); library(patchwork)
  
  df_all <- p_inc_results$data
  Tmax <- max(df_all$time, na.rm = TRUE)
  
  # --- récupérer réductions ATB depuis les params ---
  base_tau_c <- tc_results$baseline$params[["tau_c"]]
  base_tau_h <- tc_results$baseline$params[["tau_h"]]
  
  scen_names_raw <- setdiff(names(tc_results), "baseline")
  scen_map <- function(x) ifelse(x == "super", "fort", x)
  
  scen_reductions <- lapply(scen_names_raw, function(nm) {
    nm2 <- scen_map(nm)
    tau_c_sc <- tc_results[[nm]]$params[["tau_c"]]
    tau_h_sc <- tc_results[[nm]]$params[["tau_h"]]
    data.frame(
      scenario = nm2,
      red_c = 1 - tau_c_sc / base_tau_c,
      red_h = 1 - tau_h_sc / base_tau_h
    )
  }) %>% bind_rows() %>%
    group_by(scenario) %>%
    summarise(red_c = max(red_c), red_h = max(red_h), .groups="drop") %>%
    mutate(
      scenario = factor(scenario, levels = c("faible","moyen","fort")),
      red_c_lbl = sprintf("%.0f%%", 100*red_c),
      red_h_lbl = sprintf("%.0f%%", 100*red_h)
    )
  
  # Helper to create styled table grob as ggplot
  .make_table <- function(df, title_text, note_text = NULL) {
    tg <- gridExtra::tableGrob(
      df, rows = NULL,
      theme = gridExtra::ttheme_default(
        base_size = 10,
        core = list(
          fg_params = list(hjust = 0.5, x = 0.5),
          bg_params = list(fill = c("white", "#F5F5F5"))
        ),
        colhead = list(
          fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5),
          bg_params = list(fill = "#E5E7EB")
        )
      )
    )
    title <- grid::textGrob(title_text, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    if (is.null(note_text)) {
      out <- gridExtra::arrangeGrob(title, tg, heights = c(0.15, 0.85), ncol = 1)
    } else {
      note <- grid::textGrob(note_text, gp = grid::gpar(fontsize = 9, col = "gray30"),
                             hjust = 0, x = 0.02)
      out <- gridExtra::arrangeGrob(title, tg, note, heights = c(0.12, 0.78, 0.10), ncol = 1)
    }
    # Convertir en ggplot pour compatibilité avec patchwork
    ggplot2::ggplot() + 
      annotation_custom(out) +
      theme_void()
  }
  
  # =========================
  # A) Incidence INSTANTANÉE (début vs fin)
  # =========================
  inst <- df_all %>%
    filter(time %in% c(0, Tmax)) %>%
    select(time, scenario, inc_h_total_100k_year, inc_c_total_100k_year) %>%
    mutate(scenario = recode(scenario, super = "fort"),
           scenario = factor(scenario, levels = c("baseline","faible","moyen","fort")))
  
  inst_w <- inst %>%
    mutate(tp = ifelse(time == 0, "start", "end")) %>%
    select(-time) %>%
    tidyr::pivot_wider(names_from = tp, values_from = c(inc_h_total_100k_year, inc_c_total_100k_year))
  
  base_inst <- inst_w %>% filter(scenario == "baseline")
  base_end_h <- base_inst$inc_h_total_100k_year_end
  base_end_c <- base_inst$inc_c_total_100k_year_end
  
  inst_tab <- inst_w %>%
    filter(scenario != "baseline") %>%
    left_join(scen_reductions, by = "scenario") %>%
    mutate(
      chg_end_h = (inc_h_total_100k_year_end - base_end_h) / base_end_h * 100,
      chg_end_c = (inc_c_total_100k_year_end - base_end_c) / base_end_c * 100,
      scenario = factor(scenario, levels = c("faible","moyen","fort"))  # <-- AJOUTER
    ) %>%
    arrange(scenario) %>%  # <-- AJOUTER
    transmute(
      Scenario = tools::toTitleCase(as.character(scenario)),
      `ATB red.\nCommunity` = red_c_lbl,
      `ATB red.\nHospital`  = red_h_lbl,
      `Inc start\nHosp` = round(inc_h_total_100k_year_start, 2),
      `Inc end\nHosp`   = round(inc_h_total_100k_year_end, 2),
      `Δ% end\nHosp`    = sprintf("%.1f%%", chg_end_h),
      `Inc start\nComm` = round(inc_c_total_100k_year_start, 2),
      `Inc end\nComm`   = round(inc_c_total_100k_year_end, 2),
      `Δ% end\nComm`    = sprintf("%.1f%%", chg_end_c)
    )
  
  tab_inst_plot <- .make_table(
    inst_tab,
    "Instantaneous CDI incidence (start vs end of simulation)",
    "Note: Δ% computed vs baseline at end. Units: /100k pop/year"
  )
  
  # =========================
  # B) Incidence CUMULÉE annuelle (Year1 vs LastYear)
  # =========================
  df_tag <- df_all %>%
    mutate(
      scenario = recode(scenario, super = "fort"),
      scenario = factor(scenario, levels = c("baseline","faible","moyen","fort")),
      period = case_when(
        time <= (year_days - 1)          ~ "Year1",
        time >= (Tmax - (year_days - 1)) ~ "LastYear",
        TRUE                             ~ NA_character_
      )
    ) %>%
    filter(!is.na(period)) %>%
    mutate(period = factor(period, levels = c("Year1","LastYear")))
  
  inc_annual <- df_tag %>%
    group_by(scenario, period) %>%
    summarise(
      inc_h = sum(inc_h_total_100k_year, na.rm = TRUE),
      inc_c = sum(inc_c_total_100k_year,     na.rm = TRUE),
      .groups="drop"
    ) %>%
    tidyr::pivot_wider(names_from = period, values_from = c(inc_h, inc_c))
  
  base_last_h <- inc_annual %>% filter(scenario=="baseline") %>% pull(inc_h_LastYear)
  base_last_c <- inc_annual %>% filter(scenario=="baseline") %>% pull(inc_c_LastYear)
  
  inc_cum_tab <- inc_annual %>%
    filter(scenario != "baseline") %>%
    left_join(scen_reductions, by = "scenario") %>%
    mutate(
      chg_last_h = (inc_h_LastYear - base_last_h) / base_last_h * 100,
      chg_last_c = (inc_c_LastYear - base_last_c) / base_last_c * 100,
      scenario = factor(scenario, levels = c("faible","moyen","fort"))  # <-- AJOUTER CETTE LIGNE
    ) %>%
    arrange(scenario) %>%  # <-- AJOUTER CETTE LIGNE
    transmute(
      Scenario = tools::toTitleCase(as.character(scenario)),
      `ATB red.\nCommunity` = red_c_lbl,
      `ATB red.\nHospital`  = red_h_lbl,
      `Cum Y1\nHosp`  = round(inc_h_Year1, 1),
      `Cum Last\nHosp` = round(inc_h_LastYear, 1),
      `Δ% last\nHosp` = sprintf("%.1f%%", chg_last_h),
      `Cum Y1\nComm`  = round(inc_c_Year1, 1),
      `Cum Last\nComm` = round(inc_c_LastYear, 1),
      `Δ% last\nComm` = sprintf("%.1f%%", chg_last_c)
    )
  
  tab_cum_plot <- .make_table(
    inc_cum_tab,
    "Annual cumulative CDI incidence (Year 1 vs Last year)",
    "Note: Δ% computed vs baseline in the last year."
  )
  
  # =========================
  # C) Portage asymptomatique (début vs fin)
  # =========================
  car <- p_car_results$data %>%
    mutate(
      scenario = recode(scenario, super = "fort"),
      scenario = factor(scenario, levels = c("baseline","faible","moyen","fort"))
    ) %>%
    filter(time %in% c(0, Tmax)) %>%
    select(time, scenario, car_h_pct, car_c_pct) %>%
    mutate(tp = ifelse(time == 0, "start", "end")) %>%
    select(-time) %>%
    tidyr::pivot_wider(names_from = tp, values_from = c(car_h_pct, car_c_pct))
  
  base_car <- car %>% filter(scenario=="baseline")
  base_end_h <- base_car$car_h_pct_end
  base_end_c <- base_car$car_c_pct_end
  
  car_tab <- car %>%
    filter(scenario != "baseline") %>%
    left_join(scen_reductions, by="scenario") %>%
    mutate(
      chg_end_h = (car_h_pct_end - base_end_h) / base_end_h * 100,
      chg_end_c = (car_c_pct_end - base_end_c) / base_end_c * 100,
      scenario = factor(scenario, levels = c("faible","moyen","fort"))  # <-- AJOUTER
    ) %>%
    arrange(scenario) %>%  # <-- AJOUTER
    transmute(
      Scenario = tools::toTitleCase(as.character(scenario)),
      `ATB red.\nCommunity` = red_c_lbl,
      `ATB red.\nHospital`  = red_h_lbl,
      `Car start\nHosp (%)` = round(car_h_pct_start, 2),
      `Car end\nHosp (%)`   = round(car_h_pct_end, 2),
      `Δ% end\nHosp`        = sprintf("%.1f%%", chg_end_h),
      `Car start\nComm (%)` = round(car_c_pct_start, 2),
      `Car end\nComm (%)`   = round(car_c_pct_end, 2),
      `Δ% end\nComm`        = sprintf("%.1f%%", chg_end_c)
    )
  
  tab_car_plot <- .make_table(
    car_tab,
    "Asymptomatic carriage (start vs end of simulation)",
    "Note: Δ% computed vs baseline at end of simulation."
  )
  
  # =========================
  # COMBINER LES 3 TABLEAUX AVEC PATCHWORK
  # =========================
  combined_plot <- tab_inst_plot / tab_cum_plot / tab_car_plot +
    plot_annotation(
      title = "Summary of antibiotic reduction scenarios",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  list(
    incidence_instant = tab_inst_plot,
    incidence_cum     = tab_cum_plot,
    carriage          = tab_car_plot,
    combined          = combined_plot,  # <-- NOUVEAU: tous les tableaux ensemble
    data = list(inst = inst_tab, cum = inc_cum_tab, car = car_tab)
  )
}