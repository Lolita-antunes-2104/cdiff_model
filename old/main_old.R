###############################################################################
#################################### MAIN #####################################
###############################################################################

###############################################################################
# 0 : SET UP
###############################################################################

# Clear workspace
rm(list = ls())

# Packages & sources 
source("0_packages.R")
source("1_model.R")
source("2_functions.R")
source("3_calibration.R")
source("4_scenario.R") 
source("5_plot.R") 

###############################################################################
# 1 : MODEL CONSTRUCTION -------------------------------------------------------
###############################################################################

# ---- Population sizes ----
N_h <- 100
N_c <- 50000

# ---- Fixed parameters ----
nu_h <- 7
nu_c <- 7

gamma_h <- 0.012
gamma_c <- 0.012

epsilon_h <- 0
epsilon_c <- 0.06

p_h <- 0.5
p_c <- 0.5

phi_h <- 0.018
phi_c <- 0.018

k_A <- 4

omega_h <- 0
omega_c <- 0.022

prop <- 0.5  

# ---- Initial guesses (to be calibrated) ----
beta_h  <- 0.06
beta_c  <- 0.02
sigma_h <- 0.003
sigma_c <- 0.001
k_II    <- 2
k_III   <- 3

# ---- Flow rates (hospital ↔ community) ----
delta     <- 1 / 7.5
delta_I   <- 1 / 16.6
delta_II  <- 1 / 20
delta_III <- 1 / 25

# ---- Severity weights for community → hospital flows ----
w     <- 1
w_I   <- 1.5
w_II  <- 2
w_III <- 2.5

# ---- Antibiotic exposure rates ----
tau_h <- -log(1 - 0.234) / (1 / delta)
tau_c <- -log(1 - 0.0188) / 7

# ---- Parameter vector for ODE solver ----
# NOTE: keep names consistent with file 1 + file 2 expectations
params <- c(
  beta_h  = beta_h,  beta_c  = beta_c,
  sigma_h = sigma_h, sigma_c = sigma_c,
  tau_h   = tau_h,   tau_c   = tau_c,
  omega_h = omega_h, omega_c = omega_c,
  nu_h    = nu_h,    nu_c    = nu_c,
  gamma_h = gamma_h, gamma_c = gamma_c,
  epsilon_h = epsilon_h, epsilon_c = epsilon_c,
  p_h     = p_h,     p_c     = p_c,
  phi_h   = phi_h,   phi_c   = phi_c,
  k_A     = k_A,     k_II    = k_II,   k_III = k_III,
  w       = w,       w_I     = w_I,    w_II  = w_II,  w_III = w_III,
  delta   = delta,   delta_I = delta_I, delta_II = delta_II, delta_III = delta_III,
  prop    = prop
)

# ---- Initial conditions (pre-calibration) ----
prev_primo <- 0.01
prev_rec   <- 0.005
init_cond <- create_initial_conditions_precalibration(N_h, N_c, prev_primo, prev_rec)

# ---- Time vector for equilibration ----
time_burnin <- 10000
time_vec <- seq(0, time_burnin, by = 1)

# ---- Target metrics from literature ----
targets <- list(
  portage_h   = 0.07,
  portage_c   = 0.015,
  incidence_h = 12 / 100000 / 365,
  incidence_c = 18 / 100000 / 365,
  recid_1     = 0.25,
  recid_2     = 0.50
)

# ---- Quick model sanity check ----
print(body(cdiff_hc_model)) 
alpha_mode <- "dynamic"
res_eq  <- run_model_to_equilibrium(params, init_cond, time_vec, alpha_mode = alpha_mode)
last_eq <- extract_last_state(res_eq)   # excludes "time"

# ---- Equilibrium diagnostics (gradient / RHS) ----
res_eq <- run_model_to_equilibrium(params, init_cond, time_vec, alpha_mode = alpha_mode)
check_equilibrium_rhs(res_eq, params, alpha_mode = "dynamic")

# ---- Visualize dynamics BEFORE calibration ----
dynamics_pre_calib <- plot_dynamics(res_eq, targets, N_h, N_c)

# Optionally display them
dynamics_pre_calib$hospital_agg
dynamics_pre_calib$community_agg

# Visualize alpha dynamics before calibration
plot_alpha_dynamics(res_eq, params)

# ---- Compute initial metrics (before calibration) ----
N_h_init <- get_total_pop(last_eq, "h", "both")
N_c_init <- get_total_pop(last_eq, "c", "both")
N_tot_init <- N_h_init + N_c_init

outputs_init <- list(
  portage_h = get_carriage(last_eq, "total", "h", "both") / N_h_init,
  portage_c = get_carriage(last_eq, "total", "c", "both") / N_c_init,
  incidence_h = get_incidence(last_eq, params, "total", "h", "both") / N_tot_init,
  incidence_c = get_incidence(last_eq, params, "total", "c", "both") / N_tot_init,
  recid_1 = compute_recurrence_prevalence(last_eq, "rec_1", "total", "both"),
  recid_2 = compute_recurrence_prevalence(last_eq, "rec_2", "total", "both")
)

cat("\n=== BEFORE CALIBRATION ===\n")
print(outputs_init)
cat("\n=== TARGETS ===\n")
print(targets)

###############################################################################
# 2 : GRID SEARCH --------------------------------------------------------------
###############################################################################

# ---- (1) Grid search: beta_h / beta_c (carriage) ----
beta_grid <- grid_search(
  param_names   = c("beta_h", "beta_c"),
  param_ranges  = list(beta_h = c(0.01, 0.15, 30),
                       beta_c = c(0.001, 0.08, 30)),
  metric_function = compute_metrics_beta,
  target_metrics  = list(portage_h = targets$portage_h, portage_c = targets$portage_c),
  params_base   = params,
  init_cond     = init_cond,
  time_vec      = time_vec,
  n_cores       = NULL
)

best_beta <- beta_grid$best_guess[, c("beta_h", "beta_c")]
params[c("beta_h", "beta_c")] <- as.numeric(best_beta)

plot_grid_search_beta(beta_grid, targets)

# ---- (2) Grid search: sigma_h / sigma_c (incidence) ----
sigma_grid <- grid_search(
  param_names   = c("sigma_h", "sigma_c"),
  param_ranges  = list(sigma_h = c(1e-4, 1e-2, 30),
                       sigma_c = c(1e-4, 1e-2, 30)),
  metric_function = compute_metrics_sigma,
  target_metrics  = list(incidence_h = targets$incidence_h, incidence_c = targets$incidence_c),
  params_base   = params,
  init_cond     = init_cond,
  time_vec      = time_vec,
  n_cores       = NULL
)

best_sigma <- sigma_grid$best_guess[, c("sigma_h", "sigma_c")]
params[c("sigma_h", "sigma_c")] <- as.numeric(best_sigma)

plot_grid_search_sigma(sigma_grid, targets)

# ---- (3) Grid search: k_II / k_III (recurrence) ----
k_grid <- grid_search(
  param_names   = c("k_II", "k_III"),
  param_ranges  = list(k_II  = c(1, 100, 30),
                       k_III = c(1, 500, 30)),
  metric_function = compute_metrics_k,
  target_metrics  = list(recid_1 = targets$recid_1, recid_2 = targets$recid_2),
  params_base   = params,
  init_cond     = init_cond,
  time_vec      = time_vec,
  n_cores       = NULL
)

best_k <- k_grid$best_guess[, c("k_II", "k_III")]
params[c("k_II", "k_III")] <- as.numeric(best_k)

plot_grid_search_k(k_grid, targets)

###############################################################################
# 3 : CALIBRATION --------------------------------------------------------------
###############################################################################

# ---- Initial point for optimization (from grids) ----
initial_params <- params[c("beta_h", "beta_c", "sigma_h", "sigma_c", "k_II", "k_III")]

# ---- Multi-start calibration ----
calib <- run_calibration(
  initial_params = initial_params,
  target_metrics = targets,
  params_base    = params,
  init_cond      = init_cond,
  time_vec       = time_vec,
  n_starts       = 10,
  seed           = 123,
  n_cores        = NULL
)

# ---- Update params with calibrated values ----
params[names(calib$best$par)] <- calib$best$par

# ---- Equilibrium verification after calibration ----
res_eq_cal <- run_model_to_equilibrium(params, init_cond, time_vec, alpha_mode = alpha_mode)

# IMPORTANT: exclude 'time' to avoid polluting totals/tables
last_cal <- extract_last_state(res_eq_cal)

# ---- After calibration outputs ----
N_h_eq <- get_total_pop(last_cal, "h", "both")
N_c_eq <- get_total_pop(last_cal, "c", "both")
N_tot_eq <- N_h_eq + N_c_eq

outputs <- list(
  portage_h = get_carriage(last_cal, "total", "h", "both") / N_h_eq,
  portage_c = get_carriage(last_cal, "total", "c", "both") / N_c_eq,
  incidence_h = get_incidence(last_cal, params, "total", "h", "both") / N_tot_eq,
  incidence_c = get_incidence(last_cal, params, "total", "c", "both") / N_tot_eq,
  recid_1 = compute_recurrence_prevalence(last_cal, "rec_1", "total", "both"),
  recid_2 = compute_recurrence_prevalence(last_cal, "rec_2", "total", "both")
)

print(outputs)

# ---- Visualize dynamics AFTER calibration ----
dynamics_post_calib <- plot_dynamics(res_eq_cal, targets, N_h, N_c)
dynamics_post_calib$hospital_agg
dynamics_post_calib$community_agg

# ---- Equilibrium alpha (robust, uses the same logic as the model) ----
last_cal_list <- as.list(last_cal)

tot_h_eq <- compute_totals(
  last_cal_list$S0_h, last_cal_list$C0_h, last_cal_list$SA_h, last_cal_list$CA_h, last_cal_list$I_h,
  last_cal_list$S_II_h, last_cal_list$C_II_h, last_cal_list$I_II_h,
  last_cal_list$S_III_h, last_cal_list$C_III_h, last_cal_list$I_III_h
)

tot_c_eq <- compute_totals(
  last_cal_list$S0_c, last_cal_list$C0_c, last_cal_list$SA_c, last_cal_list$CA_c, last_cal_list$I_c,
  last_cal_list$S_II_c, last_cal_list$C_II_c, last_cal_list$I_II_c,
  last_cal_list$S_III_c, last_cal_list$C_III_c, last_cal_list$I_III_c
)

alpha_eq <- compute_alpha(
  alpha_mode = "dynamic",
  tot_h = tot_h_eq,
  tot_c = tot_c_eq,
  I_h = last_cal_list$I_h,
  I_II_h = last_cal_list$I_II_h,
  I_III_h = last_cal_list$I_III_h,
  I_c = last_cal_list$I_c,
  I_II_c = last_cal_list$I_II_c,
  I_III_c = last_cal_list$I_III_c,
  params = params
)$alpha

if (!is.finite(alpha_eq)) alpha_eq <- 0

# ---- Tables ----
create_all_parameters_table(params, N_h, N_c, alpha_eq)
create_stratified_parameters_table(params, alpha_eq)
create_calibration_comparison_table(outputs, targets, params)
create_outputs_table(last_cal, params, N_h, N_c)
create_detailed_metrics_table(last_cal, params)


saveRDS(
  list(
    params   = params,
    state_eq = last_cal,
    outputs  = outputs,
    alpha_eq = alpha_eq
  ),
  "results_calibration.rds"
)


res <- readRDS("results_calibration.rds")
names(res)
res$params
res$state_eq
res$outputs
res$alpha_eq


###############################################################################
###############################################################################
# 1. EQUILIBRIUM WITH FIXED ALPHA
###############################################################################

# Run to equilibrium with dynamic alpha
res_eq <- run_model_to_equilibrium(
  params,
  init_cond,
  time_vec,
  alpha_mode = "dynamic"
)

last_eq <- extract_last_state(res_eq)

# Compute equilibrium alpha
last <- as.list(last_eq)

tot_h <- compute_totals(
  last$S0_h, last$C0_h, last$SA_h, last$CA_h, last$I_h,
  last$S_II_h, last$C_II_h, last$I_II_h,
  last$S_III_h, last$C_III_h, last$I_III_h
)

tot_c <- compute_totals(
  last$S0_c, last$C0_c, last$SA_c, last$CA_c, last$I_c,
  last$S_II_c, last$C_II_c, last$I_II_c,
  last$S_III_c, last$C_III_c, last$I_III_c
)

alpha_eq <- compute_alpha(
  alpha_mode = "dynamic",
  tot_h = tot_h,
  tot_c = tot_c,
  I_h = last$I_h, I_II_h = last$I_II_h, I_III_h = last$I_III_h,
  I_c = last$I_c, I_II_c = last$I_II_c, I_III_c = last$I_III_c,
  params = params
)

# Fixed-alpha parameters
params_fixed <- params
params_fixed["alpha"]     <- alpha_eq$alpha
params_fixed["alpha_I"]   <- alpha_eq$alpha_I
params_fixed["alpha_II"]  <- alpha_eq$alpha_II
params_fixed["alpha_III"] <- alpha_eq$alpha_III

###############################################################################
# 2. SUMMARY FUNCTION
###############################################################################

summarize_state <- function(state, params) {
  
  N_h <- get_total_pop(state, "h", "both")
  N_c <- get_total_pop(state, "c", "both")
  N_tot <- N_h + N_c
  
  C_h <- get_carriage(state, "total", "h", "both")
  C_c <- get_carriage(state, "total", "c", "both")
  
  I_h <- get_incidence(state, params, "total", "h", "both")
  I_c <- get_incidence(state, params, "total", "c", "both")
  
  A_h <- get_dysbiosis(state, "h", "both")
  A_c <- get_dysbiosis(state, "c", "both")
  
  data.frame(
    carriage_h = C_h / N_h,
    carriage_c = C_c / N_c,
    incidence_h = I_h / N_tot,
    incidence_c = I_c / N_tot,
    C_h = C_h,
    C_c = C_c,
    I_h = I_h,
    I_c = I_c,
    A_h = A_h,
    A_c = A_c
  )
}

###############################################################################
# 3. ATB SCENARIOS
###############################################################################

horizon <- 5 * 365
times   <- seq(0, horizon, by = 1)

# Run baseline with fixed alpha
ode_baseline <- run_model_to_equilibrium(
  params_fixed,
  last_eq,
  times,
  alpha_mode = "fixed",
  atb_reduction_h = 0,
  atb_reduction_c = 0
)

last_baseline <- extract_last_state(ode_baseline)
baseline_summary <- summarize_state(last_baseline, params_fixed)

# ATB reduction scenarios
atb_sets <- list(
  set1 = c(0.00, 0.01, 1.00),
  set2 = c(0.05, 0.10, 0.20)
)

atb_results <- list()

for (set in names(atb_sets)) {
  for (red in atb_sets[[set]]) {
    
    ode <- run_model_to_equilibrium(
      params_fixed,
      last_eq,
      times,
      alpha_mode = "fixed",
      atb_reduction_h = 0,
      atb_reduction_c = red
    )
    
    final <- extract_last_state(ode)
    sc_sum <- summarize_state(final, params_fixed)
    
    out <- data.frame(
      scenario = paste0("ATB_c_", red),
      reduction = red,
      delta_carriage_h = sc_sum$carriage_h - baseline_summary$carriage_h,
      delta_carriage_c = sc_sum$carriage_c - baseline_summary$carriage_c,
      delta_incidence_h = sc_sum$incidence_h - baseline_summary$incidence_h,
      delta_incidence_c = sc_sum$incidence_c - baseline_summary$incidence_c,
      delta_C_h = sc_sum$C_h - baseline_summary$C_h,
      delta_C_c = sc_sum$C_c - baseline_summary$C_c,
      delta_I_h = sc_sum$I_h - baseline_summary$I_h,
      delta_I_c = sc_sum$I_c - baseline_summary$I_c,
      delta_A_h = sc_sum$A_h - baseline_summary$A_h,
      delta_A_c = sc_sum$A_c - baseline_summary$A_c
    )
    
    atb_results[[paste(set, red)]] <- out
  }
}

atb_table <- do.call(rbind, atb_results)
print(atb_table)

###############################################################################
# 4. VACCINATION SCENARIOS
###############################################################################

VC_vals <- c(0.25, 0.50, 0.75)
VE_vals <- c(0.25, 0.50, 0.75)

init_vacc <- create_initial_conditions_vaccination(
  last_state_without_vacc = last_eq,
  coverage = 0.25
)

# Vérifications
print(length(init_vacc))  # Devrait être 66 (11 compartiments × 2 settings × 2 vacc status + 12 CumI)
print(sum(is.na(init_vacc)))  # Devrait être 0
print(names(init_vacc))  # Vérifie que tous les noms se terminent par _v ou _nv


for (VC in VC_vals) {
  
  # Initial conditions with fixed vaccination coverage
  init_vacc <- create_initial_conditions_vaccination(
    last_state_without_vacc = last_eq,
    coverage = VC
  )
  
  # Baseline vaccinated (VE = 0)
  ode_base_vacc <- run_model_to_equilibrium(
    params_fixed,
    init_vacc,
    times,
    model_fn = cdiff_hc_vacc_model,  # ← IMPORTANT
    alpha_mode = "fixed",
    VE = 0,
    vacc_type = "both"
  )
  
  last_base_vacc <- extract_last_state(ode_base_vacc)
  base_vacc_summary <- summarize_state(last_base_vacc, params_fixed)
  
  # Vaccination scenarios (VE > 0)
  for (VE in VE_vals) {
    
    ode <- run_model_to_equilibrium(
      params_fixed,
      init_vacc,
      times,
      model_fn = cdiff_hc_vacc_model,  # ← IMPORTANT
      alpha_mode = "fixed",
      VE = VE,
      vacc_type = "both"
    )
    
    final <- extract_last_state(ode)
    sc_sum <- summarize_state(final, params_fixed)
    
    out <- data.frame(
      VC = VC,
      VE = VE,
      delta_carriage_h  = sc_sum$carriage_h  - base_vacc_summary$carriage_h,
      delta_carriage_c  = sc_sum$carriage_c  - base_vacc_summary$carriage_c,
      delta_incidence_h = sc_sum$incidence_h - base_vacc_summary$incidence_h,
      delta_incidence_c = sc_sum$incidence_c - base_vacc_summary$incidence_c,
      delta_C_h = sc_sum$C_h - base_vacc_summary$C_h,
      delta_C_c = sc_sum$C_c - base_vacc_summary$C_c,
      delta_I_h = sc_sum$I_h - base_vacc_summary$I_h,
      delta_I_c = sc_sum$I_c - base_vacc_summary$I_c,
      delta_A_h = sc_sum$A_h - base_vacc_summary$A_h,
      delta_A_c = sc_sum$A_c - base_vacc_summary$A_c
    )
    
    vacc_results[[paste0("VC_", VC, "_VE_", VE)]] <- out
  }
}

vacc_table <- do.call(rbind, vacc_results)
print(vacc_table)










# 
# ###############################################################################
# # Build equilibrium-based initial conditions for scenarios
# ###############################################################################
# 
# # --- 1) Run model to equilibrium with dynamic alpha ---
# alpha_mode <- "dynamic"
# 
# res_eq_cal <- run_model_to_equilibrium(
#   params,
#   init_cond,
#   time_vec,
#   alpha_mode = alpha_mode
# )
# 
# # Extract last state (remove time)
# last_eq <- extract_last_state(res_eq_cal)
# 
# # --- 2) Compute equilibrium alpha (same logic as in the model) ---
# last_eq_list <- as.list(last_eq)
# 
# tot_h_eq <- compute_totals(
#   last_eq_list$S0_h, last_eq_list$C0_h, last_eq_list$SA_h, last_eq_list$CA_h,
#   last_eq_list$I_h,
#   last_eq_list$S_II_h, last_eq_list$C_II_h, last_eq_list$I_II_h,
#   last_eq_list$S_III_h, last_eq_list$C_III_h, last_eq_list$I_III_h
# )
# 
# tot_c_eq <- compute_totals(
#   last_eq_list$S0_c, last_eq_list$C0_c, last_eq_list$SA_c, last_eq_list$CA_c,
#   last_eq_list$I_c,
#   last_eq_list$S_II_c, last_eq_list$C_II_c, last_eq_list$I_II_c,
#   last_eq_list$S_III_c, last_eq_list$C_III_c, last_eq_list$I_III_c
# )
# 
# alpha_eq <- compute_alpha(
#   alpha_mode = "dynamic",
#   tot_h = tot_h_eq,
#   tot_c = tot_c_eq,
#   I_h = last_eq_list$I_h,
#   I_II_h = last_eq_list$I_II_h,
#   I_III_h = last_eq_list$I_III_h,
#   I_c = last_eq_list$I_c,
#   I_II_c = last_eq_list$I_II_c,
#   I_III_c = last_eq_list$I_III_c,
#   params = params
# )
# 
# # Safety check
# if (!is.finite(alpha_eq$alpha)) alpha_eq$alpha <- 0
# 
# # --- 3) Build fixed-alpha parameter vector ---
# params_fixed_alpha <- params
# params_fixed_alpha["alpha"]      <- alpha_eq$alpha
# params_fixed_alpha["alpha_I"]    <- alpha_eq$alpha_I
# params_fixed_alpha["alpha_II"]   <- alpha_eq$alpha_II
# params_fixed_alpha["alpha_III"]  <- alpha_eq$alpha_III
# 
# # --- 4) Rerun model with FIXED alpha to get clean equilibrium ICs ---
# res_eq_fixed <- run_model_to_equilibrium(
#   params_fixed_alpha,
#   last_eq,
#   time_vec,
#   alpha_mode = "fixed"
# )
# 
# init_cond_scenarios <- extract_last_state(res_eq_fixed)
# 
# # --- 5) (Optional sanity check) ---
# check_equilibrium_rhs(res_eq_fixed, params_fixed_alpha, alpha_mode = "fixed")
# 
# 
# ###############################################################################
# # ATB TEST SCENARIOS — STEP CHANGE AT t = 0
# ###############################################################################
# 
# # Horizon temporel du scénario
# horizon_days <- 5 * 365
# times_scen   <- seq(0, horizon_days, by = 1)
# 
# # Définition des scénarios (réduction ATB communauté uniquement)
# atb_test_scenarios <- list(
#   baseline = list(red_c = 0.00),
#   min      = list(red_c = 0.01),
#   max      = list(red_c = 1.00)
# )
# 
# # Conteneur résultats
# atb_test_results <- list()
# 
# for (sc_name in names(atb_test_scenarios)) {
#   
#   red_c <- atb_test_scenarios[[sc_name]]$red_c
#   
#   ode_sc <- run_model_to_equilibrium(
#     params_fixed_alpha,
#     init_cond_scenarios,
#     times_scen,
#     alpha_mode = "fixed",
#     atb_reduction_h = 0,
#     atb_reduction_c = red_c
#   )
#   
#   atb_test_results[[sc_name]] <- list(
#     scenario = sc_name,
#     red_c    = red_c,
#     ode      = ode_sc
#   )
# }
# 
# 
# # Vérifier que baseline est bien stable
# check_equilibrium_rhs(
#   atb_test_results$baseline$ode,
#   params_fixed_alpha,
#   alpha_mode = "fixed"
# )
# 
# tc_atb <- run_atb_timecourse_scenarios(
#   params_fixed_alpha,
#   init_cond_scenarios
# )
# 
# p_inc_atb <- plot_atb_incidence(tc_atb)
# p_car_atb <- plot_atb_carriage(tc_atb)
# 
# print(p_inc_atb$timecourses)
# print(p_inc_atb$instant_compare)
# 
# print(p_car_atb$timecourses)
# print(p_car_atb$instant_compare)
# 
# 
# ###############################################################################
# # VACCINATIONS TEST SCENARIOS — STEP CHANGE AT t = 0
# ###############################################################################
# 
# last_cal <- extract_last_state(res_eq_cal)
# 
# init_eq <- create_initial_conditions_from_equilibrium(last_cal)
# 
# 
# # Vaccination scenario settings
# 
# vacc_coverage <- 0.70   # 70% vaccinated at start AND end
# 
# ###############################################################################
# # Initial conditions for vaccination scenarios
# ###############################################################################
# 
# init_vacc <- create_initial_conditions_vaccination(
#   last_state_without_vacc = init_eq,
#   coverage = vacc_coverage
# )
# 
# ###############################################################################
# # Vaccination + ATB scenarios
# ###############################################################################
# 
# vacc_scenarios <- list(
#   vacc1        = list(vacc_type = "vacc1", atb_red = 0),
#   vacc2        = list(vacc_type = "vacc2", atb_red = 0),
#   vacc3        = list(vacc_type = "vacc3", atb_red = 0),
#   vacc1_atb    = list(vacc_type = "vacc1", atb_red = 0.5),
#   vacc2_atb    = list(vacc_type = "vacc2", atb_red = 0.5),
#   vacc3_atb    = list(vacc_type = "vacc3", atb_red = 0.5)
# )
# 
# ###############################################################################
# # VACCINATION SCENARIOS — RUNNER
# ###############################################################################
# 
# run_vaccination_scenarios <- function(params_fixed_alpha,
#                                       init_vacc,
#                                       horizon_days = 5 * 365,
#                                       scenarios) {
#   
#   times <- seq(0, horizon_days, by = 1)
#   out <- list()
#   
#   for (sc in names(scenarios)) {
#     
#     vacc_type <- scenarios[[sc]]$vacc_type
#     atb_red   <- scenarios[[sc]]$atb_red
#     
#     ode <- run_model_to_equilibrium(
#       params_fixed_alpha,
#       init_vacc,
#       times,
#       alpha_mode = "fixed",
#       vacc_type = vacc_type,
#       atb_reduction_h = 0,
#       atb_reduction_c = atb_red
#     )
#     
#     out[[sc]] <- ode
#   }
#   
#   out
# }
# 
# ###############################################################################
# # Run vaccination scenarios
# ###############################################################################
# 
# vacc_results <- run_vaccination_scenarios(
#   params_fixed_alpha = params_fixed_alpha,
#   init_vacc          = init_vacc,
#   horizon_days       = 5 * 365,
#   scenarios          = vacc_scenarios
# )
