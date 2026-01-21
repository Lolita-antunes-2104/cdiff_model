###############################################################################
########################## 6 : VACCINATION SCENARIOS ##########################
###############################################################################

###############################################################################
# ---- 1) Metrics time series from ODE output ----
###############################################################################
# Cette fonction transforme la sortie brute de l’ODE (états compartimentaux)
# en séries temporelles de métriques épidémiologiques interprétables :
# - incidence CDI (hospitalière et communautaire)
# - prévalence de portage asymptomatique
#
# Elle est générique et utilisée par tous les scénarios vaccin / ATB.
###############################################################################
compute_vacc_timeseries_metrics <- function(ode_df,
                                            params_vec,
                                            params_baseline = NULL) {
  
  # Liste des colonnes nécessaires issues du modèle ODE
  # (compartiments colonisés, infectés et susceptibles)
  needed <- c(
    "time",
    "C0_h","CA_h","C_II_h","C_III_h","C0_c","CA_c","C_II_c","C_III_c",
    "S0_h","SA_h","S_II_h","S_III_h","I_h","I_II_h","I_III_h",
    "S0_c","SA_c","S_II_c","S_III_c","I_c","I_II_c","I_III_c"
  )
  
  # Vérification de cohérence : arrêt immédiat si des colonnes manquent
  miss <- setdiff(needed, colnames(ode_df))
  if (length(miss) > 0) {
    stop("Missing columns in ODE output: ", paste(miss, collapse = ", "))
  }
  
  # Boucle sur chaque pas de temps pour calculer les métriques
  out_mat <- vapply(seq_len(nrow(ode_df)), function(i) {
    
    # Extraction de l’état du système au temps i
    st <- as.list(ode_df[i, , drop = FALSE])
    
    # Tailles de population hospitalière et communautaire
    Nh <- compute_population_totals(st, "h")
    Nc <- compute_population_totals(st, "c")
    
    # ---- Incidence CDI ----
    # À t = 0, on force l’utilisation des paramètres baseline
    # afin de garantir la continuité mathématique au moment
    # de l’introduction de l’intervention (step change).
    if (!is.null(params_baseline) && ode_df$time[i] == 0) {
      inc_h <- compute_CDI_incidence(st, params_baseline, "h", "total") / Nh
      inc_c <- compute_CDI_incidence(st, params_baseline, "c", "total") / Nc
    } else {
      # Pour t > 0, on utilise les paramètres du scénario courant
      inc_h <- compute_CDI_incidence(st, params_vec, "h", "total") / Nh
      inc_c <- compute_CDI_incidence(st, params_vec, "c", "total") / Nc
    }
    
    # ---- Portage asymptomatique ----
    # Prévalence instantanée de portage (colonisation)
    car_h   <- compute_carriage_prevalence(st, "h")
    car_c   <- compute_carriage_prevalence(st, "c")
    car_tot <- compute_carriage_prevalence(st, "both")
    
    # Sortie des métriques, mises à l’échelle pour interprétation
    c(
      inc_h_100k_beddays_day = inc_h * 1e5,   # /100k bed-days / jour
      inc_c_100k_pop_day     = inc_c * 1e5,   # /100k pop / jour
      car_h_pct              = car_h * 100,   # %
      car_c_pct              = car_c * 100,   # %
      car_total_pct          = car_tot * 100  # %
    )
    
  }, numeric(5))
  
  # Mise en forme finale en data.frame
  out <- as.data.frame(t(out_mat))
  out$time <- ode_df$time
  
  return(out)
}

###############################################################################
# ---- 2A) Direct vaccine timecourse: sigma <- sigma * (1 - VE*VC) ----
###############################################################################
# Cette fonction simule un vaccin "direct" contre le CDI,
# modélisé comme une réduction du taux de progression
# colonisation -> infection (sigma).
#
# L’effet du vaccin est proportionnel :
#   - à la couverture vaccinale (VC)
#   - à l’efficacité vaccinale (VE)
#
# L’intervention est implémentée comme un step change
# appliqué à t = 0, après l’équilibre.
###############################################################################
run_direct_vaccine_timecourse <- function(params_calibrated,
                                          init_cond,
                                          time_to_eq,
                                          horizon_days = 5*365,
                                          VC_levels = c(0.2, 0.4, 0.6, 0.8),
                                          VE = c(0.3, 0.6, 0.9)) {
  
  # ---- 1) Calcul de l’état d’équilibre baseline ----
  # Le système est d’abord laissé évoluer jusqu’à l’équilibre
  # avant toute intervention.
  y0 <- .get_equilibrium_y0(params_calibrated, init_cond, time_to_eq)
  
  # ---- Baseline incidence at equilibrium (I0) ----
  st_eq <- as.list(y0)
  
  Nh_eq <- compute_population_totals(st_eq, "h")
  Nc_eq <- compute_population_totals(st_eq, "c")
  
  I0_h <- compute_CDI_incidence(st_eq, params_calibrated, "h", "total") / Nh_eq * 1e5
  I0_c <- compute_CDI_incidence(st_eq, params_calibrated, "c", "total") / Nc_eq * 1e5
  
  I0_cum_h <- 365 * I0_h
  I0_cum_c <- 365 * I0_c
  
  # Grille temporelle post-intervention
  times <- seq(0, horizon_days, by = 1)
  
  # Liste de sortie contenant baseline + scénarios vaccinaux
  out_list <- list()
  
  # ---- 2) Simulation baseline (sans vaccination) ----
  # Sert de référence pour toutes les comparaisons
  out_list$baseline <- .make_baseline_block(y0, times, params_calibrated)
  
  # ---- 3) Boucle sur les niveaux de couverture et d’efficacité ----
  for (vc in VC_levels) {
    for (ve in VE) {
      
      # Facteur d’efficacité globale du vaccin
      # (part de la population effectivement protégée)
      eff <- 1 - (ve * vc)
      eff <- max(0, eff)  # sécurité numérique
      
      # Copie des paramètres calibrés
      p_sc <- params_calibrated
      
      # Application du vaccin :
      # réduction du taux sigma (hospitalier et communautaire)
      p_sc[["sigma_h"]] <- params_calibrated[["sigma_h"]] * eff
      p_sc[["sigma_c"]] <- params_calibrated[["sigma_c"]] * eff
      
      # Nom du scénario (traçabilité)
      sc_name <- sprintf(
        "direct_VC%02d_VE%02d",
        round(100 * vc),
        round(100 * ve)
      )
      
      # ---- 4) Simulation ODE post-vaccination ----
      ode_sc <- as.data.frame(
        deSolve::lsoda(
          y     = y0,
          times = times,
          func  = cdiff_micro,
          parms = p_sc
        )
      )
      
      # ---- 5) Calcul des métriques ----
      # params_baseline est passé explicitement
      # pour assurer la continuité à t = 0
      met_sc <- compute_vacc_timeseries_metrics(
        ode_sc,
        params_vec      = p_sc,
        params_baseline = params_calibrated
      )
      
      # Métadonnées du scénario
      met_sc$vacc_type <- "direct"
      met_sc$scenario  <- sc_name
      met_sc$VC        <- vc
      met_sc$VE        <- ve
      met_sc$red_c     <- NA_real_
      met_sc$red_h     <- NA_real_
      
      # Stockage du scénario
      out_list[[sc_name]] <- list(
        params  = p_sc,
        ode     = ode_sc,
        metrics = met_sc
      )
    }
  }
  
  attr(out_list, "I0_cum_h") <- I0_cum_h
  attr(out_list, "I0_cum_c") <- I0_cum_c
  
  return(out_list)
}

###############################################################################
# ---- 2B) Indirect vaccine 
###############################################################################
# ---- 2C) Les deux en même temps
###############################################################################


###############################################################################
# ---- 3) Plot timecourses avec courbes distinctes ----
###############################################################################
plot_vacc_timecourses_incidence <- function(tc_results, title = "CDI incidence after a direct vaccination (step change at t=0)") {
  library(dplyr); library(ggplot2); library(patchwork)
  
  df_all <- bind_rows(lapply(tc_results, function(x) x$metrics))
  Tmax <- max(df_all$time)
  
  # Palette pour VC (similaire aux scénarios ATB)
  pal_vc <- c(
    "0.2" = "#C77DFF",  # violet (faible)
    "0.4" = "#7C3AED",  # violet foncé
    "0.6" = "#60A5FA",  # bleu (moyen)
    "0.8" = "#6EE7B7"   # vert (fort)
  )
  
  # baseline (flat equilibrium)
  base <- df_all %>% filter(scenario == "baseline", time == 0)
  base_h <- base$inc_h_100k_beddays_day
  base_c <- base$inc_c_100k_pop_day
  
  # Préparer les données
  df_plot <- df_all %>% 
    filter(scenario != "baseline") %>%
    mutate(
      VC_factor = factor(VC, levels = c(0.2, 0.4, 0.6, 0.8)),
      VC_label = paste0("VC = ", scales::percent(VC, accuracy = 1)),
      VE_label = paste0("VE = ", scales::percent(VE, accuracy = 1))
    )
  
  # ---------- TIMECOURSES ----------
  p_ts_h <- ggplot(df_plot, aes(x = time, y = inc_h_100k_beddays_day, color = VC_factor, linetype = VC_factor)) +
    geom_hline(yintercept = base_h, linetype = "dashed", color = "#DC2626", linewidth = 1) +
    geom_line(linewidth = 1) +
    facet_wrap(~VE_label, ncol = 3) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = pal_vc, name = "Vaccine Coverage") +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid"), name = "Vaccine Coverage") +
    labs(x = "Time (days)", y = "CDI per 100k bed-days / day", title = "Hospital")
  
  p_ts_c <- ggplot(df_plot, aes(x = time, y = inc_c_100k_pop_day, color = VC_factor, linetype = VC_factor)) +
    geom_hline(yintercept = base_c, linetype = "dashed", color = "#DC2626", linewidth = 1) +
    geom_line(linewidth = 1) +
    facet_wrap(~VE_label, ncol = 3) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = pal_vc, name = "Vaccine Coverage") +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid"), name = "Vaccine Coverage") +
    labs(x = "Time (days)", y = "CDI per 100k pop / day", title = "Community")
  
  p_ts <- p_ts_h + p_ts_c +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(title = title)
  
  # ---------- ANNUAL COMPARE (LAST YEAR, cumulative) ----------
  inc_last_year <- df_all %>%
    filter(time >= (Tmax - 364)) %>%
    group_by(scenario, VC, VE) %>%
    summarise(
      inc_cum_h = sum(inc_h_100k_beddays_day, na.rm = TRUE),
      inc_cum_c = sum(inc_c_100k_pop_day, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Baseline cumulative incidence = equilibrium incidence × 365
  base_end_h <- attr(tc_results, "I0_cum_h")
  base_end_c <- attr(tc_results, "I0_cum_c")
  
  recap_h <- inc_last_year %>%
    filter(scenario != "baseline") %>%
    mutate(
      rel_change_end = (inc_cum_h - base_end_h) / base_end_h * 100,
      VC_factor = factor(VC, levels = c(0.2, 0.4, 0.6, 0.8)),
      VE_factor = factor(VE, levels = c(0.3, 0.6, 0.9)),
      VC_label = paste0("VC=", scales::percent(VC, accuracy = 1)),
      VE_label = paste0("VE=", scales::percent(VE, accuracy = 1))
    )
  
  recap_c <- inc_last_year %>%
    filter(scenario != "baseline") %>%
    mutate(
      rel_change_end = (inc_cum_c - base_end_c) / base_end_c * 100,
      VC_factor = factor(VC, levels = c(0.2, 0.4, 0.6, 0.8)),
      VE_factor = factor(VE, levels = c(0.3, 0.6, 0.9)),
      VC_label = paste0("VC=", scales::percent(VC, accuracy = 1)),
      VE_label = paste0("VE=", scales::percent(VE, accuracy = 1))
    )
  
  # Style ATB inversé : barres groupées par VE, une barre par VC
  p_annual_h <- ggplot(recap_h, aes(x = VE_label, y = rel_change_end, fill = VC_factor)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(
      values = c("0.2" = "#C77DFF", "0.4" = "#7C3AED", "0.6" = "#60A5FA", "0.8" = "#6EE7B7"),
      name = "Vaccine Coverage",
      labels = c("VC=20%", "VC=40%", "VC=60%", "VC=80%")
    ) +
    labs(x = "Vaccine Efficacy", y = "Relative change vs baseline (%)", title = "Hospital") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  p_annual_c <- ggplot(recap_c, aes(x = VE_label, y = rel_change_end, fill = VC_factor)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(
      values = c("0.2" = "#C77DFF", "0.4" = "#7C3AED", "0.6" = "#60A5FA", "0.8" = "#6EE7B7"),
      name = "Vaccine Coverage",
      labels = c("VC=20%", "VC=40%", "VC=60%", "VC=80%")
    ) +
    labs(x = "Vaccine Efficacy", y = "Relative change vs baseline (%)", title = "Community") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  p_annual_compare <- p_annual_h + p_annual_c +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(title = "Annual CDI incidence change vs baseline (last year)")
  
  # ---------- INSTANT COMPARE (END OF SIM, instantaneous) ----------
  inst_end <- df_all %>%
    filter(time == Tmax) %>%
    select(scenario, VC, VE, inc_h_100k_beddays_day, inc_c_100k_pop_day)
  
  base_inst_h <- inst_end %>% filter(scenario == "baseline") %>% pull(inc_h_100k_beddays_day)
  base_inst_c <- inst_end %>% filter(scenario == "baseline") %>% pull(inc_c_100k_pop_day)
  
  inst_h <- inst_end %>%
    filter(scenario != "baseline") %>%
    mutate(
      rel_change_inst = (inc_h_100k_beddays_day - base_inst_h) / base_inst_h * 100,
      VC_factor = factor(VC, levels = c(0.2, 0.4, 0.6, 0.8)),
      VE_factor = factor(VE, levels = c(0.3, 0.6, 0.9)),
      VC_label = paste0("VC=", scales::percent(VC, accuracy = 1)),
      VE_label = paste0("VE=", scales::percent(VE, accuracy = 1))
    )
  
  inst_c <- inst_end %>%
    filter(scenario != "baseline") %>%
    mutate(
      rel_change_inst = (inc_c_100k_pop_day - base_inst_c) / base_inst_c * 100,
      VC_factor = factor(VC, levels = c(0.2, 0.4, 0.6, 0.8)),
      VE_factor = factor(VE, levels = c(0.3, 0.6, 0.9)),
      VC_label = paste0("VC=", scales::percent(VC, accuracy = 1)),
      VE_label = paste0("VE=", scales::percent(VE, accuracy = 1))
    )
  
  # Style ATB inversé : barres groupées par VE, une barre par VC
  p_inst_h <- ggplot(inst_h, aes(x = VE_label, y = rel_change_inst, fill = VC_factor)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(
      values = c("0.2" = "#C77DFF", "0.4" = "#7C3AED", "0.6" = "#60A5FA", "0.8" = "#6EE7B7"),
      name = "Vaccine Coverage",
      labels = c("VC=20%", "VC=40%", "VC=60%", "VC=80%")
    ) +
    labs(x = "Vaccine Efficacy", y = "Relative change vs baseline (%)", title = "Hospital (instant at end)") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  p_inst_c <- ggplot(inst_c, aes(x = VE_label, y = rel_change_inst, fill = VC_factor)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(
      values = c("0.2" = "#C77DFF", "0.4" = "#7C3AED", "0.6" = "#60A5FA", "0.8" = "#6EE7B7"),
      name = "Vaccine Coverage",
      labels = c("VC=20%", "VC=40%", "VC=60%", "VC=80%")
    ) +
    labs(x = "Vaccine Efficacy", y = "Relative change vs baseline (%)", title = "Community (instant at end)") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  p_inst_compare <- p_inst_h + p_inst_c +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(title = "Instantaneous CDI incidence change vs baseline (end of simulation)")
  
  list(
    timecourses = p_ts,
    annual_compare = p_annual_compare,
    instant_compare = p_inst_compare,
    data = df_all,
    recap_data_h = recap_h,
    recap_data_c = recap_c,
    inst_data_h = inst_h,
    inst_data_c = inst_c
  )
}

###############################################################################
# ---- 3bis) Plot carriage timecourses STYLE ATB ----
###############################################################################
plot_vacc_timecourses_carriage <- function(tc_results, title = "Asymptomatic carriage after a direct vaccination (step change at t=0)") {
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  
  df_all <- bind_rows(lapply(tc_results, function(x) x$metrics))
  Tmax <- max(df_all$time)
  
  # Palette pour VC
  pal_vc <- c(
    "0.2" = "#C77DFF",
    "0.4" = "#7C3AED",
    "0.6" = "#60A5FA",
    "0.8" = "#6EE7B7"
  )
  
  # baseline values
  base <- df_all %>% filter(scenario == "baseline", time == 0)
  base_h <- base$car_h_pct
  base_c <- base$car_c_pct
  
  # Préparer données
  df_plot <- df_all %>% 
    filter(scenario != "baseline") %>%
    mutate(
      VC_factor = factor(VC, levels = c(0.2, 0.4, 0.6, 0.8)),
      VC_label = paste0("VC = ", scales::percent(VC, accuracy = 1)),
      VE_label = paste0("VE = ", scales::percent(VE, accuracy = 1))
    )
  
  # ---------- TIMECOURSES: Hospital + Community ----------
  car_long <- df_plot %>%
    pivot_longer(c(car_h_pct, car_c_pct),
                 names_to = "pop", values_to = "car") %>%
    mutate(pop = recode(pop,
                        car_h_pct = "Hospital",
                        car_c_pct = "Community"))
  
  # Baseline pour chaque population
  base_long <- data.frame(
    pop = c("Hospital", "Community"),
    baseline_val = c(base_h, base_c)
  )
  
  p_ts <- ggplot(car_long, aes(time, car, color = VC_factor, linetype = VC_factor)) +
    geom_hline(data = base_long, aes(yintercept = baseline_val),
               linetype = "dashed", color = "#DC2626", linewidth = 1) +
    geom_line(linewidth = 1) +
    facet_grid(pop ~ VE_label, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = pal_vc, name = "Vaccine Coverage") +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid"), name = "Vaccine Coverage") +
    labs(x = "Time (days)", y = "Asymptomatic carriage prevalence (%)")
  
  p_ts <- p_ts + plot_annotation(title = title)
  
  # ---------- ANNUAL COMPARE ----------
  car_end <- df_all %>%
    filter(time == Tmax) %>%
    select(scenario, VC, VE, car_h_pct, car_c_pct)
  
  base_end_h <- car_end %>% filter(scenario == "baseline") %>% pull(car_h_pct)
  base_end_c <- car_end %>% filter(scenario == "baseline") %>% pull(car_c_pct)
  
  recap_h <- car_end %>%
    filter(scenario != "baseline") %>%
    mutate(
      rel_change = (car_h_pct - base_end_h) / base_end_h * 100,
      VC_factor = factor(VC, levels = c(0.2, 0.4, 0.6, 0.8)),
      VE_factor = factor(VE, levels = c(0.3, 0.6, 0.9)),
      VC_label = paste0("VC=", scales::percent(VC, accuracy = 1)),
      VE_label = paste0("VE=", scales::percent(VE, accuracy = 1))
    )
  
  recap_c <- car_end %>%
    filter(scenario != "baseline") %>%
    mutate(
      rel_change = (car_c_pct - base_end_c) / base_end_c * 100,
      VC_factor = factor(VC, levels = c(0.2, 0.4, 0.6, 0.8)),
      VE_factor = factor(VE, levels = c(0.3, 0.6, 0.9)),
      VC_label = paste0("VC=", scales::percent(VC, accuracy = 1)),
      VE_label = paste0("VE=", scales::percent(VE, accuracy = 1))
    )
  
  # Style ATB inversé : barres groupées par VE, une barre par VC
  p_annual_h <- ggplot(recap_h, aes(x = VE_label, y = rel_change, fill = VC_factor)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(
      values = c("0.2" = "#C77DFF", "0.4" = "#7C3AED", "0.6" = "#60A5FA", "0.8" = "#6EE7B7"),
      name = "Vaccine Coverage",
      labels = c("VC=20%", "VC=40%", "VC=60%", "VC=80%")
    ) +
    labs(x = "Vaccine Efficacy", y = "Relative change vs baseline (%)", title = "Hospital") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  p_annual_c <- ggplot(recap_c, aes(x = VE_label, y = rel_change, fill = VC_factor)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(
      values = c("0.2" = "#C77DFF", "0.4" = "#7C3AED", "0.6" = "#60A5FA", "0.8" = "#6EE7B7"),
      name = "Vaccine Coverage",
      labels = c("VC=20%", "VC=40%", "VC=60%", "VC=80%")
    ) +
    labs(x = "Vaccine Efficacy", y = "Relative change vs baseline (%)", title = "Community") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  p_annual_compare <- p_annual_h + p_annual_c +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(title = "Final asymptomatic carriage prevalence change vs baseline")
  
  list(
    timecourses = p_ts,
    annual_compare = p_annual_compare,
    data = df_all,
    recap_data_h = recap_h,
    recap_data_c = recap_c
  )
}

###############################################################################
# ---- 4B) Indirect vaccine: instantaneous + annual comparisons ----
###############################################################################
plot_indirect_vaccine_comparisons <- function(tc_indirect, year_days = 365) {
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  
  df_all <- bind_rows(lapply(tc_indirect, function(x) x$metrics))
  Tmax <- max(df_all$time, na.rm = TRUE)
  
  # Palette
  pal <- c(
    "faible" = "#C77DFF",
    "moyen"  = "#60A5FA",
    "fort"   = "#6EE7B7"
  )
  
  df_all <- df_all %>%
    mutate(scenario = factor(scenario, levels = c("baseline","faible","moyen","fort")))
  
  # ---------- INSTANTANEOUS (end) ----------
  inst_end <- df_all %>%
    filter(time == Tmax) %>%
    select(scenario, inc_h_100k_beddays_day, inc_c_100k_pop_day)
  
  base_inst_h <- inst_end %>% filter(scenario == "baseline") %>% pull(inc_h_100k_beddays_day)
  base_inst_c <- inst_end %>% filter(scenario == "baseline") %>% pull(inc_c_100k_pop_day)
  
  inst_h <- inst_end %>%
    filter(scenario != "baseline") %>%
    mutate(rel_change = (inc_h_100k_beddays_day - base_inst_h) / base_inst_h * 100)
  
  inst_c <- inst_end %>%
    filter(scenario != "baseline") %>%
    mutate(rel_change = (inc_c_100k_pop_day - base_inst_c) / base_inst_c * 100)
  
  p_inst_h <- ggplot(inst_h, aes(x = scenario, y = rel_change, fill = scenario)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = pal[c("faible","moyen","fort")]) +
    labs(x = NULL, y = "Relative change vs baseline (%)", title = "Hospital (instant at end)") +
    theme(legend.position = "none")
  
  p_inst_c <- ggplot(inst_c, aes(x = scenario, y = rel_change, fill = scenario)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = pal[c("faible","moyen","fort")]) +
    labs(x = NULL, y = "Relative change vs baseline (%)", title = "Community (instant at end)") +
    theme(legend.position = "none")
  
  p_inst <- p_inst_h + p_inst_c +
    plot_annotation(title = "Indirect vaccine: instantaneous CDI incidence change vs baseline (end of simulation)")
  
  # ---------- ANNUAL CUMULATIVE ----------
  inc_last_year <- df_all %>%
    filter(time >= (Tmax - (year_days - 1))) %>%
    group_by(scenario) %>%
    summarise(
      inc_cum_h = sum(inc_h_100k_beddays_day, na.rm = TRUE),
      inc_cum_c = sum(inc_c_100k_pop_day, na.rm = TRUE),
      .groups = "drop"
    )
  
  inc_last_year <- inc_last_year %>% filter(scenario != "baseline")
  
  base_end_h <- inc_last_year %>% filter(scenario == "baseline") %>% pull(inc_cum_h)
  base_end_c <- inc_last_year %>% filter(scenario == "baseline") %>% pull(inc_cum_c)
  
  recap_h <- inc_last_year %>%
    filter(scenario != "baseline") %>%
    mutate(rel_change = (inc_cum_h - base_end_h) / base_end_h * 100)
  
  recap_c <- inc_last_year %>%
    filter(scenario != "baseline") %>%
    mutate(rel_change = (inc_cum_c - base_end_c) / base_end_c * 100)
  
  p_annual_h <- ggplot(recap_h, aes(x = scenario, y = rel_change, fill = scenario)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = pal[c("faible","moyen","fort")]) +
    labs(x = NULL, y = "Relative change vs baseline (%)", title = "Hospital") +
    theme(legend.position = "none")
  
  p_annual_c <- ggplot(recap_c, aes(x = scenario, y = rel_change, fill = scenario)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_bw() +
    scale_fill_manual(values = pal[c("faible","moyen","fort")]) +
    labs(x = NULL, y = "Relative change vs baseline (%)", title = "Community") +
    theme(legend.position = "none")
  
  p_annual <- p_annual_h + p_annual_c +
    plot_annotation(title = "Indirect vaccine: annual CDI incidence change vs baseline (last year)")
  
  list(
    instant_compare = p_inst,
    annual_compare = p_annual,
    data_instant_h = inst_h,
    data_instant_c = inst_c,
    data_annual_h = recap_h,
    data_annual_c = recap_c
  )
}

###############################################################################
# ---- 5) Créer un tableau récapitulatif des scénarios sous forme de plot ----
###############################################################################
create_vacc_scenarios_summary_table <- function(tc_results, p_inc_results, p_car_results, year_days = 365) {
  library(dplyr); library(gridExtra); library(grid); library(patchwork)
  
  df_all <- p_inc_results$data
  Tmax <- max(df_all$time, na.rm = TRUE)
  
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
    select(time, scenario, VC, VE, inc_h_100k_beddays_day, inc_c_100k_pop_day) %>%
    filter(scenario != "baseline")
  
  inst_w <- inst %>%
    mutate(tp = ifelse(time == 0, "start", "end")) %>%
    select(-time) %>%
    tidyr::pivot_wider(names_from = tp, values_from = c(inc_h_100k_beddays_day, inc_c_100k_pop_day))
  
  # Baseline values
  base_inst <- df_all %>% 
    filter(scenario == "baseline", time %in% c(0, Tmax)) %>%
    select(time, inc_h_100k_beddays_day, inc_c_100k_pop_day) %>%
    mutate(tp = ifelse(time == 0, "start", "end")) %>%
    select(-time) %>%
    tidyr::pivot_wider(names_from = tp, values_from = c(inc_h_100k_beddays_day, inc_c_100k_pop_day))
  
  base_end_h <- base_inst$inc_h_100k_beddays_day_end
  base_end_c <- base_inst$inc_c_100k_pop_day_end
  
  inst_tab <- inst_w %>%
    mutate(
      chg_end_h = (inc_h_100k_beddays_day_end - base_end_h) / base_end_h * 100,
      chg_end_c = (inc_c_100k_pop_day_end - base_end_c) / base_end_c * 100,
      VC_pct = scales::percent(VC, accuracy = 1),
      VE_pct = scales::percent(VE, accuracy = 1)
    ) %>%
    arrange(VE, VC) %>%
    transmute(
      `VC` = VC_pct,
      `VE` = VE_pct,
      `Inc start\nHosp` = round(inc_h_100k_beddays_day_start, 2),
      `Inc end\nHosp`   = round(inc_h_100k_beddays_day_end, 2),
      `Δ% end\nHosp`    = sprintf("%.1f%%", chg_end_h),
      `Inc start\nComm` = round(inc_c_100k_pop_day_start, 2),
      `Inc end\nComm`   = round(inc_c_100k_pop_day_end, 2),
      `Δ% end\nComm`    = sprintf("%.1f%%", chg_end_c)
    )
  
  tab_inst_plot <- .make_table(
    inst_tab,
    "Instantaneous CDI incidence (start vs end of simulation)",
    "Note: Δ% computed vs baseline at end. Units: /100k bed-days/day (Hosp), /100k pop/day (Comm)."
  )
  
  # =========================
  # B) Incidence CUMULÉE annuelle
  # Baseline = équilibre analytique (365 × I0)
  # Post-intervention = somme dynamique sur la dernière année
  # =========================
  
  # Identifier la dernière année de simulation
  df_last_year <- df_all %>%
    filter(time >= (Tmax - (year_days - 1)))
  
  # Incidence cumulée sur la dernière année (scénarios uniquement)
  inc_annual <- df_last_year %>%
    group_by(scenario, VC, VE) %>%
    summarise(
      inc_h_LastYear = sum(inc_h_100k_beddays_day, na.rm = TRUE),
      inc_c_LastYear = sum(inc_c_100k_pop_day,     na.rm = TRUE),
      .groups = "drop"
    )
  
  # Baseline cumulée analytique (équilibre × 365)
  base_last_h <- attr(tc_results, "I0_cum_h")
  base_last_c <- attr(tc_results, "I0_cum_c")
  
  # Tableau récapitulatif (comparaison vs baseline)
  inc_cum_tab <- inc_annual %>%
    filter(scenario != "baseline") %>%
    mutate(
      chg_last_h = (inc_h_LastYear - base_last_h) / base_last_h * 100,
      chg_last_c = (inc_c_LastYear - base_last_c) / base_last_c * 100,
      VC_pct = scales::percent(VC, accuracy = 1),
      VE_pct = scales::percent(VE, accuracy = 1)
    ) %>%
    arrange(VE, VC) %>%
    transmute(
      `VC` = VC_pct,
      `VE` = VE_pct,
      `Cum Last\nHosp` = round(inc_h_LastYear, 1),
      `Δ% last\nHosp`  = sprintf("%.1f%%", chg_last_h),
      `Cum Last\nComm` = round(inc_c_LastYear, 1),
      `Δ% last\nComm`  = sprintf("%.1f%%", chg_last_c)
    )
  
  tab_cum_plot <- .make_table(
    inc_cum_tab,
    "Annual cumulative CDI incidence (last year vs baseline equilibrium)",
    "Baseline = 365 × incidence at equilibrium; scenarios = dynamic cumulative incidence over last year."
  )
  
  
  # =========================
  # C) Portage asymptomatique (début vs fin)
  # =========================
  car <- p_car_results$data %>%
    filter(time %in% c(0, Tmax)) %>%
    select(time, scenario, VC, VE, car_h_pct, car_c_pct) %>%
    filter(scenario != "baseline") %>%
    mutate(tp = ifelse(time == 0, "start", "end")) %>%
    select(-time) %>%
    tidyr::pivot_wider(names_from = tp, values_from = c(car_h_pct, car_c_pct))
  
  # Baseline carriage
  base_car <- p_car_results$data %>%
    filter(scenario == "baseline", time %in% c(0, Tmax)) %>%
    select(time, car_h_pct, car_c_pct) %>%
    mutate(tp = ifelse(time == 0, "start", "end")) %>%
    select(-time) %>%
    tidyr::pivot_wider(names_from = tp, values_from = c(car_h_pct, car_c_pct))
  
  base_end_h <- base_car$car_h_pct_end
  base_end_c <- base_car$car_c_pct_end
  
  car_tab <- car %>%
    mutate(
      chg_end_h = (car_h_pct_end - base_end_h) / base_end_h * 100,
      chg_end_c = (car_c_pct_end - base_end_c) / base_end_c * 100,
      VC_pct = scales::percent(VC, accuracy = 1),
      VE_pct = scales::percent(VE, accuracy = 1)
    ) %>%
    arrange(VE, VC) %>%
    transmute(
      `VC` = VC_pct,
      `VE` = VE_pct,
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
      title = "Summary of direct vaccination scenarios",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  list(
    incidence_instant = tab_inst_plot,
    incidence_cum     = tab_cum_plot,
    carriage          = tab_car_plot,
    combined          = combined_plot,
    data = list(inst = inst_tab, cum = inc_cum_tab, car = car_tab)
  )
}
plot_vacc_carriage_comparisons <- function(tc_results, year_days = 365, vacc_type = "direct") {
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  
  df_all <- bind_rows(lapply(tc_results, function(x) x$metrics))
  Tmax <- max(df_all$time, na.rm = TRUE)
  
  # ---------- END OF SIMULATION ----------
  car_end <- df_all %>%
    filter(time == Tmax) %>%
    select(scenario, VC, VE, car_h_pct, car_c_pct)
  
  base_end_h <- car_end %>% filter(scenario == "baseline") %>% pull(car_h_pct)
  base_end_c <- car_end %>% filter(scenario == "baseline") %>% pull(car_c_pct)
  
  if (vacc_type == "direct") {
    # Barres pour direct vaccine avec facet par VE
    car_rel <- car_end %>%
      filter(scenario != "baseline") %>%
      mutate(
        rel_change_h = (car_h_pct - base_end_h) / base_end_h * 100,
        rel_change_c = (car_c_pct - base_end_c) / base_end_c * 100,
        VC_label = paste0("VC = ", scales::percent(VC, accuracy = 1)),
        VE_label = paste0("VE = ", scales::percent(VE, accuracy = 1)),
        VC_label = factor(VC_label, levels = paste0("VC = ", scales::percent(sort(unique(VC)), accuracy = 1))),
        VE_label = factor(VE_label, levels = paste0("VE = ", scales::percent(sort(unique(VE)), accuracy = 1)))
      )
    
    # Palette de couleurs pour VC (cohérente avec les timecourses)
    vc_colors <- c(
      "VC = 20%" = "#C77DFF",
      "VC = 40%" = "#7C3AED",
      "VC = 60%" = "#60A5FA",
      "VC = 80%" = "#6EE7B7"
    )
    
    # Hospital
    p_h <- ggplot(car_rel, aes(x = VC_label, y = rel_change_h, fill = VC_label)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      facet_wrap(~VE_label, ncol = 1) +
      theme_bw() +
      scale_fill_manual(values = vc_colors) +
      labs(
        x = NULL, 
        y = "Relative change vs baseline (%)", 
        title = "Hospital"
      ) +
      theme(
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Community
    p_c <- ggplot(car_rel, aes(x = VC_label, y = rel_change_c, fill = VC_label)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      facet_wrap(~VE_label, ncol = 1) +
      theme_bw() +
      scale_fill_manual(values = vc_colors) +
      labs(
        x = NULL, 
        y = "Relative change vs baseline (%)", 
        title = "Community"
      ) +
      theme(
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    p_final <- p_h + p_c +
      plot_annotation(title = "Direct vaccine: final asymptomatic carriage change vs baseline")
    
  } else {
    # Bars for indirect vaccine (comme avant)
    pal <- c("faible" = "#C77DFF", "moyen" = "#60A5FA", "fort" = "#6EE7B7")
    
    car_rel_h <- car_end %>%
      filter(scenario != "baseline") %>%
      mutate(
        rel_change = (car_h_pct - base_end_h) / base_end_h * 100,
        scenario = factor(scenario, levels = c("faible","moyen","fort"))
      )
    
    car_rel_c <- car_end %>%
      filter(scenario != "baseline") %>%
      mutate(
        rel_change = (car_c_pct - base_end_c) / base_end_c * 100,
        scenario = factor(scenario, levels = c("faible","moyen","fort"))
      )
    
    p_h <- ggplot(car_rel_h, aes(x = scenario, y = rel_change, fill = scenario)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      theme_bw() +
      scale_fill_manual(values = pal) +
      labs(x = NULL, y = "Relative change vs baseline (%)", title = "Hospital") +
      theme(legend.position = "none")
    
    p_c <- ggplot(car_rel_c, aes(x = scenario, y = rel_change, fill = scenario)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      theme_bw() +
      scale_fill_manual(values = pal) +
      labs(x = NULL, y = "Relative change vs baseline (%)", title = "Community") +
      theme(legend.position = "none")
    
    p_final <- p_h + p_c +
      plot_annotation(title = "Indirect vaccine: final asymptomatic carriage change vs baseline")
  }
  
  list(final_carriage = p_final, data = car_end)
}