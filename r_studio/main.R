###############################################################################
#################################### MAIN #####################################
###############################################################################

# ---- Packages & sources ----
source("0_packages.R")
source("1_model_def.R")
source("2_metrics.R")
source("3_grid_search.R")
source("4_calibration.R")
source("5_scenario_atb.R")
source("6_scenario_vacc.R")
source("plot.R")

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

flow_ratio <- N_h / N_c
alpha <- delta * flow_ratio
alpha_I <- delta_I * flow_ratio
alpha_II <- delta_II * flow_ratio
alpha_III <- delta_III * flow_ratio

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
  alpha   = alpha,   alpha_I = alpha_I, alpha_II = alpha_II, alpha_III = alpha_III,
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


###############################################################################
# ---- 1 : GRID SEARCH ----
###############################################################################

cat("\n=== GRID SEARCH ===\n")

# Time grid (defined in 1_model_def.R)
time_vec <- time

# ---- (1) Grid search: beta_h / beta_c (carriage) ----
beta_ranges <- list(
  beta_h = c(0.001, 0.1, 20),
  beta_c = c(0.001, 0.1, 20)
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
  sigma_h = c(0.00001, 0.01, 25),
  sigma_c = c(0.00001, 0.01, 25)
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
  k_II  = c(1, 200, 25),
  k_III = c(1, 200, 25)
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

cat("\n=== MULTISTRAT CALIBRATION ===\n")

initial_params_calib <- c(
  beta_h  = params["beta_h"],
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

# ---- Plots ----
plots_dyn <- plot_dynamics(final_results$ode_result, final_results$params_final, targets, N_h, N_c)
print(plots_dyn$hospital_agg + plots_dyn$community_agg)

R0_analysis <- plot_R0_vs_nu(final_results$params_final, init.cond, time_vec)
print(R0_analysis$plot)

# ---- Tables ----
create_all_parameters_table(final_results$params_final, N_h, N_c)
create_calibration_table(final_results$params_final, final_results$metrics, targets)
create_outputs_table(final_results$metrics, N_h, N_c)
create_detailed_metrics_table(final_results$metrics)

# ---- Vérifier si équilibre ----
eq_df_test <- run_model_to_equilibrium(final_results$params_final, init.cond, time)
y0_test <- unlist(eq_df_test[nrow(eq_df_test), -1], use.names = TRUE)
dy0_test <- cdiff_micro(0, y0_test, final_results$params_final)[[1]]
cat(sprintf("\n[Equilibrium check] max(|dy/dt|) = %.3e\n", max(abs(dy0_test))))

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
    super  = list(red_c = 0.25, red_h = 0.10),
    moyen  = list(red_c = 0.20, red_h = 0.05),
    faible = list(red_c = 0.20, red_h = 0.00)
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





















