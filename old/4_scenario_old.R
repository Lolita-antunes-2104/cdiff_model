###############################################################################
# ATB SCENARIOS — METRICS
###############################################################################

# Compute instantaneous incidence + carriage time series from an ODE output
compute_atb_timeseries_metrics <- function(ode_df, params_vec) {
  
  needed <- c(
    "time",
    "C0_h","CA_h","C_II_h","C_III_h",
    "C0_c","CA_c","C_II_c","C_III_c",
    "S0_h","SA_h","S_II_h","S_III_h",
    "I_h","I_II_h","I_III_h",
    "S0_c","SA_c","S_II_c","S_III_c",
    "I_c","I_II_c","I_III_c"
  )
  
  miss <- setdiff(needed, colnames(ode_df))
  if (length(miss) > 0) {
    stop("Missing compartments in ODE output: ", paste(miss, collapse = ", "))
  }
  
  out <- lapply(seq_len(nrow(ode_df)), function(i) {
    
    st <- as.numeric(ode_df[i, -1])
    names(st) <- colnames(ode_df)[-1]
    
    Nh <- get_total_pop(st, "h", "both")
    Nc <- get_total_pop(st, "c", "both")
    Ntot <- Nh + Nc
    
    inc_h <- get_incidence(st, params_vec, "total", "h", "both") / Ntot * 365 * 1e5
    inc_c <- get_incidence(st, params_vec, "total", "c", "both") / Ntot * 365 * 1e5
    
    car_h <- get_carriage(st, "total", "h", "both") * 100
    car_c <- get_carriage(st, "total", "c", "both") * 100
    
    data.frame(
      time = ode_df$time[i],
      inc_h_100k_year = inc_h,
      inc_c_100k_year = inc_c,
      car_h_pct = car_h,
      car_c_pct = car_c
    )
  })
  
  do.call(rbind, out)
}



###############################################################################
# ATB SCENARIOS — RUNNER
###############################################################################

run_atb_timecourse_scenarios <- function(params_fixed_alpha,
                                         init_cond_scenarios,
                                         horizon_days = 5 * 365,
                                         scenarios = list(
                                           baseline = list(red_c = 0),
                                           min      = list(red_c = 0.01),
                                           max      = list(red_c = 1)
                                         )) {
  
  times <- seq(0, horizon_days, by = 1)
  out <- list()
  
  for (sc_name in names(scenarios)) {
    
    red_c <- scenarios[[sc_name]]$red_c
    
    ode <- run_model_to_equilibrium(
      params_fixed_alpha,
      init_cond_scenarios,
      times,
      alpha_mode = "fixed",
      atb_reduction_h = 0,
      atb_reduction_c = red_c
    )
    
    metrics <- compute_atb_timeseries_metrics(ode, params_fixed_alpha)
    metrics$scenario <- sc_name
    
    out[[sc_name]] <- list(
      ode = ode,
      metrics = metrics,
      red_c = red_c
    )
  }
  
  out
}

###############################################################################
# ATB SCENARIOS — INCIDENCE PLOTS
###############################################################################

label_atb_scenarios <- function(df, tc_results) {
  
  # extraire les réductions ATB depuis tc_results
  scen_info <- lapply(names(tc_results), function(sc) {
    data.frame(
      scenario = sc,
      red_c = tc_results[[sc]]$red_c
    )
  }) |> dplyr::bind_rows()
  
  scen_info <- scen_info |>
    dplyr::mutate(
      scenario = dplyr::recode(scenario,
                               baseline = "baseline",
                               min = "faible",
                               max = "fort"),
      label = dplyr::case_when(
        red_c == 0    ~ "Baseline (0% ATB reduction)",
        red_c == 0.01 ~ "ATB −1% (community)",
        red_c == 1    ~ "ATB −100% (community)",
        TRUE          ~ paste0("ATB −", round(100 * red_c), "% (community)")
      )
    )
  
  df |>
    dplyr::left_join(scen_info[, c("scenario", "label")], by = "scenario")
}


###############################################################################
# ATB SCENARIOS — INCIDENCE PLOTS
###############################################################################

plot_atb_incidence <- function(tc_results) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  
  # ------------------------------------------------------------------
  # Build dataframe
  # ------------------------------------------------------------------
  df <- bind_rows(lapply(tc_results, `[[`, "metrics"))
  
  df <- df %>%
    mutate(
      scenario = recode(scenario,
                        min = "faible",
                        max = "fort"),
      scenario = factor(scenario, levels = c("baseline", "faible", "fort")),
      label = case_when(
        scenario == "baseline" ~ "Baseline (0% ATB reduction)",
        scenario == "faible"   ~ "ATB −1% (community)",
        scenario == "fort"     ~ "ATB −100% (community)"
      )
    )
  
  pal <- c(
    "Baseline (0% ATB reduction)" = "#808080",
    "ATB −1% (community)"         = "#C77DFF",
    "ATB −100% (community)"       = "#6EE7B7"
  )
  
  # ------------------------------------------------------------------
  # Force baseline to be flat (value at t = 0)
  # ------------------------------------------------------------------
  baseline_t0 <- df %>%
    filter(scenario == "baseline", time == min(time)) %>%
    select(-time)
  
  df <- df %>%
    filter(scenario != "baseline") %>%
    bind_rows(
      baseline_t0 %>%
        slice(rep(1, length(unique(df$time)))) %>%
        mutate(
          time = sort(unique(df$time)),
          scenario = "baseline",
          label = "Baseline (0% ATB reduction)"
        )
    )
  
  # ------------------------------------------------------------------
  # TIMECOURSES (Hospital + Community)
  # ------------------------------------------------------------------
  inc_long <- df %>%
    pivot_longer(
      cols = c(inc_h_100k_year, inc_c_100k_year),
      names_to = "pop",
      values_to = "inc"
    ) %>%
    mutate(
      pop = factor(
        pop,
        levels = c("inc_h_100k_year", "inc_c_100k_year"),
        labels = c("Hospital", "Community")
      )
    )
  
  p_ts <- ggplot(
    inc_long,
    aes(time, inc, color = label, linetype = label)
  ) +
    geom_line(linewidth = 1) +
    facet_wrap(~pop, scales = "free_y") +
    theme_bw() +
    scale_color_manual(
      values = pal,
      breaks = c(
        "Baseline (0% ATB reduction)",
        "ATB −1% (community)",
        "ATB −100% (community)"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Baseline (0% ATB reduction)" = "dashed",
        "ATB −1% (community)"         = "solid",
        "ATB −100% (community)"       = "solid"
      ),
      breaks = c(
        "Baseline (0% ATB reduction)",
        "ATB −1% (community)",
        "ATB −100% (community)"
      )
    ) +
    labs(
      x = "Time (days)",
      y = "CDI incidence per 100k / year",
      color = "Scenario",
      linetype = "Scenario",
      title = "Instantaneous CDI incidence after antibiotic reduction"
    )
  
  # ------------------------------------------------------------------
  # INSTANT COMPARE (end of simulation) — Hospital + Community + Δ%
  # ------------------------------------------------------------------
  Tmax <- max(df$time, na.rm = TRUE)
  
  inst <- df %>% filter(time == Tmax)
  base_vals <- inst %>% filter(scenario == "baseline")
  
  inst_bar <- inst %>%
    pivot_longer(
      cols = c(inc_h_100k_year, inc_c_100k_year),
      names_to = "pop",
      values_to = "inc"
    ) %>%
    mutate(
      pop = factor(
        pop,
        levels = c("inc_h_100k_year", "inc_c_100k_year"),
        labels = c("Hospital", "Community")
      )
    ) %>%
    left_join(
      base_vals %>%
        pivot_longer(
          cols = c(inc_h_100k_year, inc_c_100k_year),
          names_to = "pop",
          values_to = "base_inc"
        ),
      by = "pop"
    ) %>%
    mutate(
      delta_pct = (inc - base_inc) / base_inc * 100
    )
  
  inst_bar <- inst_bar %>%
    left_join(
      inst %>% select(scenario, label) %>% distinct(),
      by = "scenario"
    )
  
  p_inst <- ggplot(
    inst_bar,
    aes(pop, inc, fill = label)
  ) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(
      aes(label = ifelse(scenario == "baseline",
                         "",
                         sprintf("%.1f%%", delta_pct))),
      position = position_dodge(width = 0.7),
      vjust = -0.5,
      size = 3
    ) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    labs(
      x = NULL,
      y = "CDI incidence per 100k / year",
      fill = "Scenario",
      title = "CDI incidence at end of simulation (Δ% vs baseline)"
    )
  
  # ------------------------------------------------------------------
  # Return
  # ------------------------------------------------------------------
  list(
    timecourses = p_ts,
    instant_compare = p_inst,
    data = df
  )
}



###############################################################################
# ATB SCENARIOS — CARRIAGE PLOTS
###############################################################################
plot_atb_carriage <- function(tc_results) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  
  # ------------------------------------------------------------------
  # Build dataframe
  # ------------------------------------------------------------------
  df <- bind_rows(lapply(tc_results, `[[`, "metrics"))
  
  df <- df %>%
    mutate(
      scenario = recode(scenario,
                        min = "faible",
                        max = "fort"),
      scenario = factor(scenario, levels = c("baseline", "faible", "fort")),
      label = case_when(
        scenario == "baseline" ~ "Baseline (0% ATB reduction)",
        scenario == "faible"   ~ "ATB −1% (community)",
        scenario == "fort"     ~ "ATB −100% (community)"
      )
    )
  
  pal <- c(
    "Baseline (0% ATB reduction)" = "#808080",
    "ATB −1% (community)"         = "#C77DFF",
    "ATB −100% (community)"       = "#6EE7B7"
  )
  
  # ------------------------------------------------------------------
  # Force baseline to be flat (value at t = 0)
  # ------------------------------------------------------------------
  baseline_t0 <- df %>%
    filter(scenario == "baseline", time == min(time)) %>%
    select(-time)
  
  df <- df %>%
    filter(scenario != "baseline") %>%
    bind_rows(
      baseline_t0 %>%
        slice(rep(1, length(unique(df$time)))) %>%
        mutate(
          time = sort(unique(df$time)),
          scenario = "baseline",
          label = "Baseline (0% ATB reduction)"
        )
    )
  
  # ------------------------------------------------------------------
  # TIMECOURSES (Hospital + Community)
  # ------------------------------------------------------------------
  car_long <- df %>%
    pivot_longer(
      cols = c(car_h_pct, car_c_pct),
      names_to = "pop",
      values_to = "car"
    ) %>%
    mutate(
      pop = factor(
        pop,
        levels = c("car_h_pct", "car_c_pct"),
        labels = c("Hospital", "Community")
      )
    )
  
  p_ts <- ggplot(
    car_long,
    aes(time, car, color = label, linetype = label)
  ) +
    geom_line(linewidth = 1) +
    facet_wrap(~pop, scales = "free_y") +
    theme_bw() +
    scale_color_manual(
      values = pal,
      breaks = c(
        "Baseline (0% ATB reduction)",
        "ATB −1% (community)",
        "ATB −100% (community)"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Baseline (0% ATB reduction)" = "dashed",
        "ATB −1% (community)"         = "solid",
        "ATB −100% (community)"       = "solid"
      ),
      breaks = c(
        "Baseline (0% ATB reduction)",
        "ATB −1% (community)",
        "ATB −100% (community)"
      )
    ) +
    labs(
      x = "Time (days)",
      y = "Asymptomatic carriage prevalence (%)",
      color = "Scenario",
      linetype = "Scenario",
      title = "Asymptomatic carriage after antibiotic reduction"
    )
  
  # ------------------------------------------------------------------
  # INSTANT COMPARE (end of simulation) — Hospital + Community + Δ%
  # ------------------------------------------------------------------
  Tmax <- max(df$time, na.rm = TRUE)
  
  inst <- df %>% filter(time == Tmax)
  base_vals <- inst %>% filter(scenario == "baseline")
  
  inst_bar <- inst %>%
    pivot_longer(
      cols = c(car_h_pct, car_c_pct),
      names_to = "pop",
      values_to = "car"
    ) %>%
    mutate(
      pop = factor(
        pop,
        levels = c("car_h_pct", "car_c_pct"),
        labels = c("Hospital", "Community")
      )
    ) %>%
    left_join(
      base_vals %>%
        pivot_longer(
          cols = c(car_h_pct, car_c_pct),
          names_to = "pop",
          values_to = "base_car"
        ),
      by = "pop"
    ) %>%
    mutate(
      delta_pct = (car - base_car) / base_car * 100
    )
  
  p_inst <- ggplot(
    inst_bar,
    aes(pop, car, fill = label)
  ) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(
      aes(label = ifelse(scenario == "baseline",
                         "",
                         sprintf("%.1f%%", delta_pct))),
      position = position_dodge(width = 0.7),
      vjust = -0.5,
      size = 3
    ) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Asymptomatic carriage prevalence (%)",
      fill = "Scenario",
      title = "Asymptomatic carriage at end of simulation (Δ% vs baseline)"
    )
  
  # ------------------------------------------------------------------
  # Return
  # ------------------------------------------------------------------
  list(
    timecourses = p_ts,
    instant_compare = p_inst,
    data = df
  )
}
