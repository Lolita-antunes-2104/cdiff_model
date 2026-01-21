###############################################################################
#################################### MAIN #####################################
###############################################################################

# ---- Packages & sources ----
source("0_packages.R")
source("1_model_new.R")
source("2_metrics.R")
source("3_grid_search.R")
source("4_calibration.R")
source("5_scenario_atb.R")
source("plot.R")
source("6_scen_vacc.R")
source("6_scena_vacc_alphaS.R")

###############################################################################
# ---- PARAMETERS ----
###############################################################################

# ---- Population sizes ----
N_h <- 100
N_c <- 50000

# ---- Fixed parameters ----
nu <- 7
gamma <- 0.012
epsilon <- 0.055
p <- 0.5
phi <- 0.018
k_A <- 4
omega_h <- 0
omega_c <- 0.022

# ---- Initial guesses (to be calibrated) ----
beta_h <- 0.06
beta_c  <- 0.02
sigma_h <- 0.003
sigma_c <- 0.001
k_II <- 2
k_III <- 3

# ---- Flow rates (hospital ↔ community) ----
delta <- 1 / 7.5
delta_I <- 1 / 16
delta_II <- 1 / 20
delta_III <- 1 / 25

# ---- Severity weights for community → hospital flows ----
w <- 1
w_I <- 1.5
w_II <- 2
w_III <- 2.5

# ---- Antibiotic exposure rates ----
tau_h <- -log(1 - 0.234) / (1 / delta)
tau_c <- -log(1 - 0.0188) / 7

# ---- Parameter vector for ODE solver ----
params <- c(
  beta_h  = beta_h,  beta_c  = beta_c,
  sigma_h = sigma_h, sigma_c = sigma_c,
  tau_h   = tau_h,   tau_c   = tau_c,
  omega_h = omega_h, omega_c = omega_c,
  nu      = nu,      gamma   = gamma,
  epsilon = epsilon, p       = p,
  phi     = phi,
  k_A     = k_A,     k_II    = k_II,   k_III  = k_III,
  w       = w,
  w_I     = w_I,
  w_II    = w_II,
  w_III   = w_III,
  delta   = delta,   delta_I = delta_I, delta_II = delta_II, delta_III = delta_III
)

# ---- Initial conditions ----
init.cond <- create_initial_conditions(N_h, N_c)

# ---- Integration time ----
time <- seq(0, 4000, by=1)

# ---- Calibration targets ----
targets <- list(
  portage_h   = 0.069,
  portage_c   = 0.0137,
  incidence_h = 11.8 / (100000 * 365),
  incidence_c = 17.7 / (100000 * 365),
  recid_1     = 0.25,
  recid_2     = 0.50
)

###################### 
# verification dynamique correcte 
##############

res <- run_model_to_equilibrium(params, init.cond, time)

res$N_h <- with(res,
                S0_h + SA_h + S_II_h + S_III_h +
                  C0_h + CA_h + C_II_h + C_III_h +
                  I_h + I_II_h + I_III_h
)

res$N_c <- with(res,
                S0_c + SA_c + S_II_c + S_III_c +
                  C0_c + CA_c + C_II_c + C_III_c +
                  I_c + I_II_c + I_III_c
)

range(res$N_h)
range(res$N_c)


###############################################################################
# ---- 1 : GRID SEARCH ----
###############################################################################

cat("\n=== GRID SEARCH ===\n")

# Time grid (defined in 1_model_def.R)
time_vec <- time

# ---- (1) Grid search: beta_h / beta_c (carriage) ----
beta_ranges <- list(
  beta_h = c(0.001, 0.5, 30),
  beta_c = c(0.001, 0.5, 30)
)

beta_res <- grid_search(
  param_names     = c("beta_h", "beta_c"),
  param_ranges    = beta_ranges,
  metric_function = compute_metrics_beta,
  target_metrics  = list(portage_h = targets$portage_h, portage_c = targets$portage_c),
  params_base     = params,
  init_cond       = init.cond,
  time_vec        = time_vec,
  n_cores         = NULL
)

params["beta_h"] <- as.numeric(beta_res$best_guess$beta_h)
params["beta_c"] <- as.numeric(beta_res$best_guess$beta_c)

cat("\n--- BEST beta ---\n")
print(beta_res$best_guess)

p_beta <- plot_grid_search_beta(beta_res, beta_res$best_guess, targets)
print(p_beta)

# ---- (2) Grid search: sigma_h / sigma_c (incidence) ----
sigma_ranges <- list(
  sigma_h = c(0.000001, 0.01, 30),
  sigma_c = c(0.000001, 0.01, 30)
)

sigma_res <- grid_search(
  param_names     = c("sigma_h", "sigma_c"),
  param_ranges    = sigma_ranges,
  metric_function = compute_metrics_sigma,
  target_metrics  = list(incidence_h = targets$incidence_h, incidence_c = targets$incidence_c),
  params_base     = params,
  init_cond       = init.cond,
  time_vec        = time_vec,
  n_cores         = NULL
)

params["sigma_h"] <- as.numeric(sigma_res$best_guess$sigma_h)
params["sigma_c"] <- as.numeric(sigma_res$best_guess$sigma_c)

cat("\n--- BEST sigma ---\n")
print(sigma_res$best_guess)

p_sigma <- plot_grid_search_sigma(sigma_res, sigma_res$best_guess, targets)
print(p_sigma)

# ---- (3) Grid search: k_II / k_III (recurrence) ----
k_ranges <- list(
  k_II  = c(1, 1000, 30),
  k_III = c(1, 1000, 30)
)

k_res <- grid_search(
  param_names     = c("k_II", "k_III"),
  param_ranges    = k_ranges,
  metric_function = compute_metrics_k,
  target_metrics  = list(recid_1 = targets$recid_1, recid_2 = targets$recid_2),
  params_base     = params,
  init_cond       = init.cond,
  time_vec        = time_vec,
  n_cores         = NULL
)

params["k_II"]  <- as.numeric(k_res$best_guess$k_II)
params["k_III"] <- as.numeric(k_res$best_guess$k_III)

cat("\n--- BEST k ---\n")
print(k_res$best_guess)

p_k <- plot_grid_search_k(k_res, k_res$best_guess, targets)
print(p_k)

###############################################################################
# ---- 2 : CALIBRATION ----
###############################################################################

cat("\n=== TEST FONCTION OBJECTIF ===\n")

# Crée la fonction objectif
obj_test <- create_objective_function(targets, params, init.cond, time_vec)

# Teste avec les paramètres initiaux
test_log <- log(initial_params_calib)
cat("Paramètres en log:\n")
print(test_log)

# Appelle la fonction
cat("\nAppel de la fonction objectif...\n")
result_test <- obj_test(test_log)

cat("\nRésultat:", result_test, "\n")

if (result_test >= 1e10) {
  cat("ERREUR: La fonction objectif retourne 1e10 !\n")
  cat("Il y a un bug dans create_objective_function\n")
} else {
  cat("OK: La fonction objectif fonctionne, valeur =", result_test, "\n")
}

cat("\n=== MULTISTRAT CALIBRATION ===\n")

initial_params_calib <- c(
  beta_h  = params["beta_h"],    # ±50% autour du grid search
  beta_c  = params["beta_c"],
  sigma_h = params["sigma_h"],
  sigma_c = params["sigma_c"],
  k_II    = params["k_II"],
  k_III   = params["k_III"]
)
calibration_results <- run_calibration(
  initial_params = initial_params_calib,
  target_metrics = targets,
  params_base    = params,
  init_cond      = init.cond,
  time_vec       = time_vec,
  n_starts       = 10
)

# === DIAGNOSTIC ===
cat("\n=== DIAGNOSTIC CALIBRATION ===\n")
cat("Nombre de runs :", length(calibration_results$all_results), "\n")
cat("Valeurs objectif de tous les runs :\n")
print(sapply(calibration_results$all_results, function(x) x$value))
cat("\nConvergence codes :\n")
print(sapply(calibration_results$all_results, function(x) x$convergence))
cat("\nMeilleure valeur objectif :", calibration_results$best$value, "\n")

calibrated_params <- calibration_results$best$par_natural
cat("\n=== CALIBRATED PARAMETERS ===\n")
print(calibrated_params)

#Elle sert à ré-exécuter le modèle avec les paramètres calibrés, puis à calculer/retourner les métriques finales (portage, incidences, récidives, R0, erreurs vs cibles) + la trajectoire ODE
final_results <- compute_calibrated_metrics(
  calibrated_params = calibrated_params,
  params_base       = params,
  init_cond         = init.cond,
  time_vec          = time_vec,
  targets           = targets
)

# ---- Vérifier si équilibre ----
eq_df_test <- run_model_to_equilibrium(final_results$params_final, init.cond, time)
y0_test <- unlist(eq_df_test[nrow(eq_df_test), -1], use.names = TRUE)
dy0_test <- cdiff_micro(0, y0_test, final_results$params_final)[[1]]
cat(sprintf("\n[Equilibrium check] max(|dy/dt|) = %.3e\n", max(abs(dy0_test))))

# ---- Quasi-equilibrium diagnostic ----
eq_df_test <- run_model_to_equilibrium(final_results$params_final, init.cond, time)

tail_eq <- eq_df_test[(nrow(eq_df_test)-200):nrow(eq_df_test), ]

# populations
range(with(tail_eq, S0_h + SA_h + S_II_h + S_III_h + C0_h + CA_h + C_II_h + C_III_h + I_h + I_II_h + I_III_h))

range(with(tail_eq, S0_c + SA_c + S_II_c + S_III_c + C0_c + CA_c + C_II_c + C_III_c + I_c + I_II_c + I_III_c))

# epidemiological stability
max(abs(diff(tail_eq$I_h)))
max(abs(diff(tail_eq$I_c)))
max(abs(diff(tail_eq$C0_h + tail_eq$CA_h)))
max(abs(diff(tail_eq$C0_c + tail_eq$CA_c)))

# ---- Compute alpha at equilibrium (last state) ----
last_state <- as.list(final_results$ode_result[nrow(final_results$ode_result), ])

tot_h <- compute_totals(
  last_state$S0_h, last_state$SA_h, last_state$S_II_h, last_state$S_III_h,
  last_state$C0_h, last_state$CA_h, last_state$C_II_h, last_state$C_III_h,
  last_state$I_h, last_state$I_II_h, last_state$I_III_h
)

tot_c <- compute_totals(
  last_state$S0_c, last_state$SA_c, last_state$S_II_c, last_state$S_III_c,
  last_state$C0_c, last_state$CA_c, last_state$C_II_c, last_state$C_III_c,
  last_state$I_c, last_state$I_II_c, last_state$I_III_c
)

out_hc <- final_results$params_final["delta"]   * (tot_h$S + tot_h$C) +
  final_results$params_final["delta_I"] * last_state$I_h +
  final_results$params_final["delta_II"]* last_state$I_II_h +
  final_results$params_final["delta_III"]*last_state$I_III_h

den_alpha <- final_results$params_final["w"]   * (tot_c$S + tot_c$C) +
  final_results$params_final["w_I"] * last_state$I_c +
  final_results$params_final["w_II"]* last_state$I_II_c +
  final_results$params_final["w_III"]*last_state$I_III_c

alpha_eq <- out_hc / den_alpha


# ---- Plots ----
plots_dyn <- plot_dynamics(final_results$ode_result, final_results$params_final, targets, N_h, N_c)
print(plots_dyn$hospital_agg + plots_dyn$community_agg)

R0_analysis <- plot_R0_vs_nu(final_results$params_final, init.cond, time_vec)
print(R0_analysis$plot)

p_alpha <- plot_alpha_dynamics(final_results$ode_result, final_results$params_final)
print(p_alpha)


# ---- Tables ----
create_all_parameters_table(final_results$params_final, N_h, N_c, alpha_eq)
create_calibration_table(final_results$params_final, final_results$metrics, targets)
create_outputs_table(final_results$metrics, N_h, N_c)
create_detailed_metrics_table(final_results$metrics)
create_stratified_parameters_table(final_results$params_final, N_h, N_c, alpha_eq)


###############################################################################
# ---- 3 : ANTIBIOTIC REDUCTION SCENARIO ----
###############################################################################

cat("\n=== ANTIBIOTIC REDUCTION SCENARIO ===\n")

tc_atb <- run_atb_timecourse_scenarios(
  params_calibrated = final_results$params_final,
  init_cond = init.cond,
  time_to_eq = time,
  horizon_days = 5*365,
  scenarios = list(
    super  = list(red_c = 0.20, red_h = 0),
    moyen  = list(red_c = 0.10, red_h = 0),
    faible = list(red_c = 0.05, red_h = 0)
  )
)

# ---- Incidence plots ----
p_inc_atb <- plot_atb_timecourses_incidence(tc_atb)
print(p_inc_atb$timecourses)
print(p_inc_atb$annual_compare)
print(p_inc_atb$instant_compare)

# ---- Carriage plots ----
p_car_atb <- plot_atb_timecourses_carriage(tc_atb)
print(p_car_atb$timecourses)
print(p_car_atb$annual_compare)

# ---- Summary table ----
tabs <- create_atb_scenarios_summary_table(tc_atb, p_inc_atb, p_car_atb)
tabs$combined


###############################################################################
# ---- 4 : VACCINATION SCENARIO ----
###############################################################################

cat("\n=== VACCINATION SCENARIOS ===\n")

# Prepare equilibrium state from calibration
calibrated_eq <- as.list(final_results$ode_result[nrow(final_results$ode_result), ])


# Run vaccination scenarios
vacc_results <- run_vacc_scenarios(
  params_calibrated = final_results$params_final,
  calibrated_equilibrium = calibrated_eq,
  time_vec = seq(0, 3650, by = 1),  # 10 years simulation
  VE_values = c(0.3, 0.5, 0.7),
  VC_values = c(0.2, 0.4, 0.6, 0.8)
)

# Quick population check
cat("\n--- Quick population check ---\n")
test_sc <- vacc_results[["VE30_VC20"]]
first <- test_sc$ode[1, ]
last <- test_sc$ode[nrow(test_sc$ode), ]
N_h_start <- sum(first[grepl("_h_", names(first))])
N_h_end <- sum(last[grepl("_h_", names(last))])
N_c_start <- sum(first[grepl("_c_", names(first))])
N_c_end <- sum(last[grepl("_c_", names(last))])
cat(sprintf("N_h: %.2f -> %.2f (diff: %.6f)\n", N_h_start, N_h_end, abs(N_h_end - N_h_start)))
cat(sprintf("N_c: %.2f -> %.2f (diff: %.6f)\n", N_c_start, N_c_end, abs(N_c_end - N_c_start)))

# Verify test scenarios
verification <- verify_test_scenarios(vacc_results)

# Plot vaccination impact
p_vacc <- plot_vacc_impact(vacc_results, final_results$metrics)
print(p_vacc$combined)

# Individual plots if needed
print(p_vacc$inc_total)
print(p_vacc$inc_primo)
print(p_vacc$inc_rec)
print(p_vacc$portage)

# Summary table
tab_vacc <- create_vacc_summary_table(vacc_results, final_results$metrics)
grid.draw(tab_vacc)

print(calibrated_eq)
print(names(calibrated_eq))
















