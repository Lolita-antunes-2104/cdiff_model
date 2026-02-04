###############################################################################
############################ MAIN_SCENARIO ####################################
###############################################################################

###############################################################################
# 0. SET UP
###############################################################################

# Clear workspace
rm(list = ls())

# Packages & sources
source("0_package.R")
source("1_model.R")
source("2_function.R")

# Load calibrated results
res <- readRDS("results_equilibrium.rds")

# Parameters
params <- res$params_calib

# Equilibrium state to be used as initial conditions for scenarios
state_eq <- res$state_eq

# Time
time <- seq(0, 365 * 5, by = 1)

###############################################################################
# 1. TEST VALIDITY OF THE STRUCTURE
###############################################################################

# VE = 0 (no vaccine effect) and no ATB reduction
params_test <- c(params, tau_mult_red = 1, sigma_mult_v = 1)

# Helper to build initial conditions with vaccine coverage vc
make_init <- function(vc) {
  nv <- 1 - vc
  v  <- vc
  g  <- function(x) as.numeric(state_eq[x])
  c(
    # Hospital non-vaccinated
    S0_h_nv = g("S0_h") * nv, SA_h_nv = g("SA_h") * nv,
    C0_h_nv = g("C0_h") * nv, CA_h_nv = g("CA_h") * nv,
    I_h_nv  = g("I_h")  * nv,
    S_II_h_nv = g("S_II_h") * nv, C_II_h_nv = g("C_II_h") * nv,
    I_II_h_nv = g("I_II_h") * nv,
    S_III_h_nv = g("S_III_h") * nv, C_III_h_nv = g("C_III_h") * nv,
    I_III_h_nv = g("I_III_h") * nv,
    CumI_h_nv = 0, primo_CumI_h_nv = 0, rec_CumI_h_nv = 0,

    # Community non-vaccinated
    S0_c_nv = g("S0_c") * nv, SA_c_nv = g("SA_c") * nv,
    C0_c_nv = g("C0_c") * nv, CA_c_nv = g("CA_c") * nv,
    I_c_nv  = g("I_c")  * nv,
    S_II_c_nv = g("S_II_c") * nv, C_II_c_nv = g("C_II_c") * nv,
    I_II_c_nv = g("I_II_c") * nv,
    S_III_c_nv = g("S_III_c") * nv, C_III_c_nv = g("C_III_c") * nv,
    I_III_c_nv = g("I_III_c") * nv,
    CumI_c_nv = 0, primo_CumI_c_nv = 0, rec_CumI_c_nv = 0,

    # Hospital vaccinated
    S0_h_v = g("S0_h") * v, SA_h_v = g("SA_h") * v,
    C0_h_v = g("C0_h") * v, CA_h_v = g("CA_h") * v,
    I_h_v  = g("I_h")  * v,
    S_II_h_v = g("S_II_h") * v, C_II_h_v = g("C_II_h") * v,
    I_II_h_v = g("I_II_h") * v,
    S_III_h_v = g("S_III_h") * v, C_III_h_v = g("C_III_h") * v,
    I_III_h_v = g("I_III_h") * v,
    CumI_h_v = 0, primo_CumI_h_v = 0, rec_CumI_h_v = 0,

    # Community vaccinated
    S0_c_v = g("S0_c") * v, SA_c_v = g("SA_c") * v,
    C0_c_v = g("C0_c") * v, CA_c_v = g("CA_c") * v,
    I_c_v  = g("I_c")  * v,
    S_II_c_v = g("S_II_c") * v, C_II_c_v = g("C_II_c") * v,
    I_II_c_v = g("I_II_c") * v,
    S_III_c_v = g("S_III_c") * v, C_III_c_v = g("C_III_c") * v,
    I_III_c_v = g("I_III_c") * v,
    CumI_c_v = 0, primo_CumI_c_v = 0, rec_CumI_c_v = 0
  )
}

# Baseline calibration model (no vaccination compartments)
out_base <- deSolve::lsoda(y = state_eq, times = time, func = cdiff_model_for_calibration, parms = params)
metrics_base <- compute_metrics_calib(out_base, params)

# Test 1: VC = 0, VE = 0
init_test1 <- make_init(0)
out_test1 <- deSolve::lsoda(y = init_test1, times = time, func = cdiff_model_for_scenario, parms = params_test)
metrics_test1 <- compute_metrics_scenario(out_test1, params_test)

# Test 2: VC = 100, VE = 0
init_test2 <- make_init(1)
out_test2 <- deSolve::lsoda(y = init_test2, times = time, func = cdiff_model_for_scenario, parms = params_test)
metrics_test2 <- compute_metrics_scenario(out_test2, params_test)

###############################################################################
# 2. TABLES (3 RESULTS SIDE BY SIDE)
###############################################################################

# Totals table
sum_calib <- function(out) {
  x <- as.list(out[nrow(out), ])
  th <- compute_totals(x$S0_h, x$C0_h, x$SA_h, x$CA_h, x$I_h,
                       x$S_II_h, x$C_II_h, x$I_II_h, x$S_III_h, x$C_III_h, x$I_III_h)
  tc <- compute_totals(x$S0_c, x$C0_c, x$SA_c, x$CA_c, x$I_c,
                       x$S_II_c, x$C_II_c, x$I_II_c, x$S_III_c, x$C_III_c, x$I_III_c)
  R <- x$S_II_h + x$C_II_h + x$I_II_h + x$S_III_h + x$C_III_h + x$I_III_h +
       x$S_II_c + x$C_II_c + x$I_II_c + x$S_III_c + x$C_III_c + x$I_III_c
  c(S_total = th$S + tc$S, C_total = th$C + tc$C, I_total = th$I_tot + tc$I_tot, Recurrences_total = R)
}

sum_scen <- function(out) {
  x <- as.list(out[nrow(out), ])
  th_nv <- compute_totals(x$S0_h_nv, x$C0_h_nv, x$SA_h_nv, x$CA_h_nv, x$I_h_nv,
                          x$S_II_h_nv, x$C_II_h_nv, x$I_II_h_nv, x$S_III_h_nv, x$C_III_h_nv, x$I_III_h_nv)
  tc_nv <- compute_totals(x$S0_c_nv, x$C0_c_nv, x$SA_c_nv, x$CA_c_nv, x$I_c_nv,
                          x$S_II_c_nv, x$C_II_c_nv, x$I_II_c_nv, x$S_III_c_nv, x$C_III_c_nv, x$I_III_c_nv)
  th_v <- compute_totals(x$S0_h_v, x$C0_h_v, x$SA_h_v, x$CA_h_v, x$I_h_v,
                         x$S_II_h_v, x$C_II_h_v, x$I_II_h_v, x$S_III_h_v, x$C_III_h_v, x$I_III_h_v)
  tc_v <- compute_totals(x$S0_c_v, x$C0_c_v, x$SA_c_v, x$CA_c_v, x$I_c_v,
                         x$S_II_c_v, x$C_II_c_v, x$I_II_c_v, x$S_III_c_v, x$C_III_c_v, x$I_III_c_v)
  R <- x$S_II_h_nv + x$C_II_h_nv + x$I_II_h_nv + x$S_III_h_nv + x$C_III_h_nv + x$I_III_h_nv +
       x$S_II_c_nv + x$C_II_c_nv + x$I_II_c_nv + x$S_III_c_nv + x$C_III_c_nv + x$I_III_c_nv +
       x$S_II_h_v + x$C_II_h_v + x$I_II_h_v + x$S_III_h_v + x$C_III_h_v + x$I_III_h_v +
       x$S_II_c_v + x$C_II_c_v + x$I_II_c_v + x$S_III_c_v + x$C_III_c_v + x$I_III_c_v
  c(S_total = th_nv$S + tc_nv$S + th_v$S + tc_v$S,
    C_total = th_nv$C + tc_nv$C + th_v$C + tc_v$C,
    I_total = th_nv$I_tot + tc_nv$I_tot + th_v$I_tot + tc_v$I_tot,
    Recurrences_total = R)
}

tab_totals <- rbind(
  Base = sum_calib(out_base),
  Test1_VC0_VE0 = sum_scen(out_test1),
  Test2_VC100_VE0 = sum_scen(out_test2)
)

print(tab_totals)

# Epidemiological table
scale_inc <- 365 * 1e5

tab_epi <- data.frame(
  prev_h = c(tail(metrics_base$carriage$prev_h, 1),
             tail(metrics_test1$carriage$prev_h, 1),
             tail(metrics_test2$carriage$prev_h, 1)) * 100,
  prev_c = c(tail(metrics_base$carriage$prev_c, 1),
             tail(metrics_test1$carriage$prev_c, 1),
             tail(metrics_test2$carriage$prev_c, 1)) * 100,
  inc_h = c(tail(metrics_base$incidence_instant$inc_h_total_abs, 1),
            tail(metrics_test1$incidence_instant$inc_h_total_abs, 1),
            tail(metrics_test2$incidence_instant$inc_h_total_abs, 1)) * scale_inc,
  inc_c = c(tail(metrics_base$incidence_instant$inc_c_total_abs, 1),
            tail(metrics_test1$incidence_instant$inc_c_total_abs, 1),
            tail(metrics_test2$incidence_instant$inc_c_total_abs, 1)) * scale_inc,
  recid_1 = c(tail(metrics_base$recurrence$rec1, 1),
              tail(metrics_test1$recurrence$rec1, 1),
              tail(metrics_test2$recurrence$rec1, 1)) * 100,
  recid_2 = c(tail(metrics_base$recurrence$rec2, 1),
              tail(metrics_test1$recurrence$rec2, 1),
              tail(metrics_test2$recurrence$rec2, 1)) * 100
)

rownames(tab_epi) <- c("Base", "Test1_VC0_VE0", "Test2_VC100_VE0")

print(tab_epi)




###############################################################################
# 3. ATB SCENARIOS (VC = 50%, VE = 0)
###############################################################################

vc <- 0.5
init_vc50 <- make_init(vc)

make_params_atb <- function(red) {
  p <- params
  p["tau_h"] <- p["tau_h"] / red
  c(p, tau_mult_red = red, sigma_mult_v = 1)
}

out_base50 <- deSolve::lsoda(y = init_vc50, times = time, func = cdiff_model_for_scenario, parms = make_params_atb(1))
out_atb05  <- deSolve::lsoda(y = init_vc50, times = time, func = cdiff_model_for_scenario, parms = make_params_atb(0.95))
out_atb10  <- deSolve::lsoda(y = init_vc50, times = time, func = cdiff_model_for_scenario, parms = make_params_atb(0.90))
out_atb20  <- deSolve::lsoda(y = init_vc50, times = time, func = cdiff_model_for_scenario, parms = make_params_atb(0.80))

metrics_base50 <- compute_metrics_scenario(out_base50, make_params_atb(1))
metrics_05     <- compute_metrics_scenario(out_atb05,  make_params_atb(0.95))
metrics_10     <- compute_metrics_scenario(out_atb10,  make_params_atb(0.90))
metrics_20     <- compute_metrics_scenario(out_atb20,  make_params_atb(0.80))

scenarios_bar <- c("ATB -5%", "ATB -10%", "ATB -20%")
cols <- c("ATB -5%" = "#E6F4FF", "ATB -10%" = "#A8D4FF", "ATB -20%" = "#5AA8E8")
cols_prev <- c("ATB -5%" = "#F2E8FF", "ATB -10%" = "#D8C2FF", "ATB -20%" = "#B18CFF")

prev_plot <- rbind(
  data.frame(setting = "Hospital", scenario = scenarios_bar,
             value = c(metrics_05$carriage$prev_h, metrics_10$carriage$prev_h, metrics_20$carriage$prev_h)),
  data.frame(setting = "Community", scenario = scenarios_bar,
             value = c(metrics_05$carriage$prev_c, metrics_10$carriage$prev_c, metrics_20$carriage$prev_c))
)
prev_plot$scenario <- factor(prev_plot$scenario, levels = scenarios_bar)

scale_inc <- 365 * 1e5

inc_plot <- rbind(
  data.frame(type = "Total CDI", scenario = scenarios_bar,
             value = c(
               (metrics_05$incidence_instant$inc_h_total_abs + metrics_05$incidence_instant$inc_c_total_abs) * scale_inc,
               (metrics_10$incidence_instant$inc_h_total_abs + metrics_10$incidence_instant$inc_c_total_abs) * scale_inc,
               (metrics_20$incidence_instant$inc_h_total_abs + metrics_20$incidence_instant$inc_c_total_abs) * scale_inc
             )),
  data.frame(type = "Primo CDI", scenario = scenarios_bar,
             value = c(
               (metrics_05$incidence_instant$inc_h_primo_abs + metrics_05$incidence_instant$inc_c_primo_abs) * scale_inc,
               (metrics_10$incidence_instant$inc_h_primo_abs + metrics_10$incidence_instant$inc_c_primo_abs) * scale_inc,
               (metrics_20$incidence_instant$inc_h_primo_abs + metrics_20$incidence_instant$inc_c_primo_abs) * scale_inc
             )),
  data.frame(type = "Recidive CDI", scenario = scenarios_bar,
             value = c(
               (metrics_05$incidence_instant$inc_h_rec_abs + metrics_05$incidence_instant$inc_c_rec_abs) * scale_inc,
               (metrics_10$incidence_instant$inc_h_rec_abs + metrics_10$incidence_instant$inc_c_rec_abs) * scale_inc,
               (metrics_20$incidence_instant$inc_h_rec_abs + metrics_20$incidence_instant$inc_c_rec_abs) * scale_inc
             ))
)
inc_plot$type <- factor(inc_plot$type, levels = c("Total CDI", "Primo CDI", "Recidive CDI"))
inc_plot$scenario <- factor(inc_plot$scenario, levels = scenarios_bar)

base_prev <- rbind(
  data.frame(setting = "Hospital", value = metrics_base50$carriage$prev_h),
  data.frame(setting = "Community", value = metrics_base50$carriage$prev_c)
)
base_inc <- data.frame(
  type = c("Total CDI", "Primo CDI", "Recidive CDI"),
  value = c(
    (metrics_base50$incidence_instant$inc_h_total_abs + metrics_base50$incidence_instant$inc_c_total_abs) * scale_inc,
    (metrics_base50$incidence_instant$inc_h_primo_abs + metrics_base50$incidence_instant$inc_c_primo_abs) * scale_inc,
    (metrics_base50$incidence_instant$inc_h_rec_abs + metrics_base50$incidence_instant$inc_c_rec_abs) * scale_inc
  )
)
base_inc$type <- factor(base_inc$type, levels = c("Total CDI", "Primo CDI", "Recidive CDI"))
base_inc$x <- as.numeric(base_inc$type)
base_inc$x_start <- base_inc$x - 0.45
base_inc$x_end <- base_inc$x + 0.45

label_prev <- merge(prev_plot, base_prev, by = "setting", suffixes = c("", "_base"))
label_prev$lab <- sprintf("%+.0f%%", 100 * (label_prev$value / label_prev$value_base - 1))

label_inc <- merge(inc_plot, base_inc[, c("type", "value")], by = "type", suffixes = c("", "_base"))
label_inc$lab <- sprintf("%+.0f%%", 100 * (label_inc$value / label_inc$value_base - 1))

p_prev <- ggplot(prev_plot, aes(x = scenario, y = value, fill = scenario)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_hline(data = base_prev, aes(yintercept = value), color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_text(data = label_prev, aes(y = value, label = lab), vjust = 1.2, size = 3) +
  facet_wrap(~ setting, ncol = 2) +
  scale_fill_manual(values = cols_prev) +
  theme_bw() +
  labs(x = NULL, y = "Colonization prevalence", title = "Colonization prevalence (relative change)")

p_inc <- ggplot(inc_plot, aes(x = type, y = value, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", linewidth = 0.3) +
  geom_segment(data = base_inc,
               aes(x = x_start, xend = x_end, y = value, yend = value),
               color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_text(data = label_inc, aes(x = type, y = value, label = lab),
            position = position_dodge(width = 0.8), vjust = 1.2, size = 3) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  labs(x = NULL, y = "CDI incidence (/100k/year)", title = "CDI incidence (relative change)")

print(p_prev)
print(p_inc)

p_all <- p_prev / p_inc
print(p_all)
ggplot2::ggsave("plots_atb.png", p_all, width = 12, height = 10, dpi = 300)







###############################################################################
# 4. VACCINATION SCENARIOS (VC 0->100%, VE = 0.3 / 0.5 / 0.7)
###############################################################################

coverages <- c(0, 0.25, 0.50, 0.75, 0.90)
ves <- c(0.3, 0.5, 0.6)

params_v0 <- c(params, tau_mult_red = 1, sigma_mult_v = 1)
out_v0 <- deSolve::lsoda(y = make_init(0), times = time, func = cdiff_model_for_scenario, parms = params_v0)
metrics_v0 <- compute_metrics_scenario(out_v0, params_v0)

base_prev_h <- metrics_v0$carriage$prev_h
base_prev_c <- metrics_v0$carriage$prev_c
base_inc_total_h <- metrics_v0$incidence_instant$inc_h_total_abs * scale_inc
base_inc_total_c <- metrics_v0$incidence_instant$inc_c_total_abs * scale_inc
base_inc_primo_h <- metrics_v0$incidence_instant$inc_h_primo_abs * scale_inc
base_inc_primo_c <- metrics_v0$incidence_instant$inc_c_primo_abs * scale_inc
base_inc_rec_h <- metrics_v0$incidence_instant$inc_h_rec_abs * scale_inc
base_inc_rec_c <- metrics_v0$incidence_instant$inc_c_rec_abs * scale_inc

rows <- list()
idx <- 1
  for (ve in ves) {
    params_ve <- c(params, tau_mult_red = 1, sigma_mult_v = 1 - ve)
    for (vc in coverages) {
      out <- deSolve::lsoda(y = make_init(vc), times = time, func = cdiff_model_for_scenario, parms = params_ve)
      m <- compute_metrics_scenario(out, params_ve)
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "Hospital", metric = "Prevalence",
      rel = 100 * (m$carriage$prev_h / base_prev_h - 1)
    ); idx <- idx + 1
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "Community", metric = "Prevalence",
      rel = 100 * (m$carriage$prev_c / base_prev_c - 1)
    ); idx <- idx + 1
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "All", metric = "Total CDI",
      rel = 100 * (((m$incidence_instant$inc_h_total_abs + m$incidence_instant$inc_c_total_abs) * scale_inc) / (base_inc_total_h + base_inc_total_c) - 1)
    ); idx <- idx + 1
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "Hospital", metric = "Primo CDI",
      rel = 100 * ((m$incidence_instant$inc_h_primo_abs * scale_inc) / base_inc_primo_h - 1)
    ); idx <- idx + 1
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "Community", metric = "Primo CDI",
      rel = 100 * ((m$incidence_instant$inc_c_primo_abs * scale_inc) / base_inc_primo_c - 1)
    ); idx <- idx + 1
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "All", metric = "Primo CDI",
      rel = 100 * (((m$incidence_instant$inc_h_primo_abs + m$incidence_instant$inc_c_primo_abs) * scale_inc) / (base_inc_primo_h + base_inc_primo_c) - 1)
    ); idx <- idx + 1
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "Hospital", metric = "Recidive CDI",
      rel = 100 * ((m$incidence_instant$inc_h_rec_abs * scale_inc) / base_inc_rec_h - 1)
    ); idx <- idx + 1
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "Community", metric = "Recidive CDI",
      rel = 100 * ((m$incidence_instant$inc_c_rec_abs * scale_inc) / base_inc_rec_c - 1)
    ); idx <- idx + 1
    rows[[idx]] <- data.frame(
      ve = ve, vc = vc, setting = "All", metric = "Recidive CDI",
      rel = 100 * (((m$incidence_instant$inc_h_rec_abs + m$incidence_instant$inc_c_rec_abs) * scale_inc) / (base_inc_rec_h + base_inc_rec_c) - 1)
    ); idx <- idx + 1
  }
}

vacc_df <- do.call(rbind, rows)
vacc_df$ve <- factor(vacc_df$ve, levels = ves)

ve_cols_prev <- c("0.3" = "#C79BFF", "0.5" = "#7B2CFF", "0.6" = "#4A00C9")
ve_cols_inc  <- c("0.3" = "#8CC7FF", "0.5" = "#2C7DFF", "0.6" = "#0047B8")

plot_vacc <- function(metric_name, title_txt, cols_use, facet_setting = TRUE) {
  df <- vacc_df[vacc_df$metric == metric_name, ]
  p <- ggplot(df, aes(x = vc * 100, y = rel, color = ve)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values = cols_use, name = "Vaccine efficacy") +
    coord_cartesian(ylim = c(-100, 0)) +
    theme_bw() +
    labs(x = "Vaccine coverage (%)", y = "Relative change (%)", title = title_txt)
  if (facet_setting) {
    p <- p + facet_wrap(~ setting, ncol = 2)
  }
  p
}

plot_vacc_cdi <- function(cols_use) {
  df <- vacc_df[vacc_df$metric %in% c("Total CDI", "Primo CDI", "Recidive CDI") & vacc_df$setting == "All", ]
  df$metric <- factor(df$metric, levels = c("Total CDI", "Primo CDI", "Recidive CDI"))
  ggplot(df, aes(x = vc * 100, y = rel, color = ve)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ metric, ncol = 3) +
    scale_color_manual(values = cols_use, name = "Vaccine efficacy") +
    coord_cartesian(ylim = c(-100, 0)) +
    theme_bw() +
    labs(x = "Vaccine coverage (%)", y = "Relative change (%)", title = "CDI incidence (relative change)")
}

p_v_prev  <- plot_vacc("Prevalence", "Colonization prevalence (relative change)", ve_cols_prev)
p_v_cdi   <- plot_vacc_cdi(ve_cols_inc)

print(p_v_prev)
print(p_v_cdi)
 
p_v_all <- p_v_prev / p_v_cdi

print(p_v_all)
ggplot2::ggsave("plots_vaccination.png", p_v_all, width = 12, height = 10, dpi = 300)

###############################################################################
# 5. ATB + VACCINATION COMBINED (relative change)
###############################################################################

atb_levels <- c(1.00, 0.95, 0.90, 0.80)
atb_labels <- c("ATB -0%", "ATB -5%", "ATB -10%", "ATB -20%")

rows2 <- list()
idx2 <- 1
for (a in seq_along(atb_levels)) {
  red <- atb_levels[a]
  for (ve in ves) {
    params_ve_atb <- c(params, tau_mult_red = red, sigma_mult_v = 1 - ve)
    for (vc in coverages) {
      out <- deSolve::lsoda(y = make_init(vc), times = time, func = cdiff_model_for_scenario, parms = params_ve_atb)
      m <- compute_metrics_scenario(out, params_ve_atb)
      rows2[[idx2]] <- data.frame(atb = atb_labels[a], ve = ve, vc = vc, setting = "Hospital", metric = "Prevalence",
                                  rel = 100 * (m$carriage$prev_h / base_prev_h - 1)); idx2 <- idx2 + 1
      rows2[[idx2]] <- data.frame(atb = atb_labels[a], ve = ve, vc = vc, setting = "Community", metric = "Prevalence",
                                  rel = 100 * (m$carriage$prev_c / base_prev_c - 1)); idx2 <- idx2 + 1
      rows2[[idx2]] <- data.frame(atb = atb_labels[a], ve = ve, vc = vc, setting = "Hospital", metric = "Total CDI",
                                  rel = 100 * ((m$incidence_instant$inc_h_total_abs * scale_inc) / base_inc_total_h - 1)); idx2 <- idx2 + 1
      rows2[[idx2]] <- data.frame(atb = atb_labels[a], ve = ve, vc = vc, setting = "Community", metric = "Total CDI",
                                  rel = 100 * ((m$incidence_instant$inc_c_total_abs * scale_inc) / base_inc_total_c - 1)); idx2 <- idx2 + 1
      rows2[[idx2]] <- data.frame(atb = atb_labels[a], ve = ve, vc = vc, setting = "Hospital", metric = "Primo CDI",
                                  rel = 100 * ((m$incidence_instant$inc_h_primo_abs * scale_inc) / base_inc_primo_h - 1)); idx2 <- idx2 + 1
      rows2[[idx2]] <- data.frame(atb = atb_labels[a], ve = ve, vc = vc, setting = "Community", metric = "Primo CDI",
                                  rel = 100 * ((m$incidence_instant$inc_c_primo_abs * scale_inc) / base_inc_primo_c - 1)); idx2 <- idx2 + 1
      rows2[[idx2]] <- data.frame(atb = atb_labels[a], ve = ve, vc = vc, setting = "Hospital", metric = "Recidive CDI",
                                  rel = 100 * ((m$incidence_instant$inc_h_rec_abs * scale_inc) / base_inc_rec_h - 1)); idx2 <- idx2 + 1
      rows2[[idx2]] <- data.frame(atb = atb_labels[a], ve = ve, vc = vc, setting = "Community", metric = "Recidive CDI",
                                  rel = 100 * ((m$incidence_instant$inc_c_rec_abs * scale_inc) / base_inc_rec_c - 1)); idx2 <- idx2 + 1
    }
  }
}

comb_df <- do.call(rbind, rows2)
comb_df$ve <- factor(comb_df$ve, levels = ves)
comb_df$atb <- factor(comb_df$atb, levels = atb_labels)

plot_comb <- function(metric_name, title_txt, cols_use) {
  df <- comb_df[comb_df$metric == metric_name, ]
  ggplot(df, aes(x = vc * 100, y = rel, color = ve)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_grid(atb ~ setting) +
    scale_color_manual(values = cols_use, name = "Vaccine efficacy") +
    coord_cartesian(ylim = c(-100, 0)) +
    theme_bw() +
    labs(x = "Vaccine coverage (%)", y = "Relative change (%)", title = title_txt)
}

p_c_prev  <- plot_comb("Prevalence", "Colonization prevalence (relative change)", ve_cols_prev)
p_c_total <- plot_comb("Total CDI", "Total CDI (relative change)", ve_cols_inc)
p_c_primo <- plot_comb("Primo CDI", "Primo CDI (relative change)", ve_cols_inc)
p_c_rec   <- plot_comb("Recidive CDI", "Recidive CDI (relative change)", ve_cols_inc)

print(p_c_prev)
print(p_c_total)
print(p_c_primo)
print(p_c_rec)

p_c_all <- p_c_prev / p_c_total
p_c_all <- p_c_primo / p_c_rec

print(p_c_all)
ggplot2::ggsave("plots_atb_vaccination.png", p_c_all, width = 12, height = 18, dpi = 300)
