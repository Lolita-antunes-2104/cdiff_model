###############################################################################
################################ MAIN #########################################
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
source("3_calibration.R")
source("5_plot.R")


###############################################################################
# 1. PARAMETRES, POPULATIONS and TARGETS 
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
k_A <- 5
omega_h <- 0
omega_c <- 0.022
prop <- 0.5

# ---- Initial guesses (to be calibrated) ----
beta_h <- 0.05
beta_c <- 0.05
sigma_h <- 0.0005
sigma_c <- 0.0005
k_II <- 100
k_III <- 100

# ---- Flow rates (hospital ↔ community) ----
delta <- 1 / 7.5
delta_I <- 1 / 15
delta_II <- 1 / 18
delta_III <- 1 / 20

# ---- Severity weights for community → hospital flows ----
w <- 1
w_I <- 1.5
w_II <- 3
w_III <- 4

# ---- Antibiotic exposure rates ----
tau_h <- -log(1 - 0.234) / (1 / delta)
tau_c <- -log(1 - 0.0188) / 7

# ---- Parameter vector for ODE solver ----
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

alpha_fixed <- if ((w * N_c) == 0) 0 else delta * N_h / (w * N_c)
params <- c(params, alpha_fixed = alpha_fixed)





###############################################################################
# 2. INITIAL CONDITIONS & TIME FOR CALIBRATION
###############################################################################

# Initial conditions for the calibration model (hospital + community, primary + recurrences)
init_cond <- create_initial_conditions_precalibration(
  N_h = N_h,
  N_c = N_c,
  prev_primo = 0.01,   # initial primary colonization/infection "seed" (choose small)
  prev_rec   = 0.005   # initial recurrence "seed" (choose small)
)

# Time grids (days)
# Calibration runs (grid + optim)
time_vec_calib <- seq(0, 50000, by = 1)
# Equilibrium run (can be longer if needed)
time_vec_eq <- seq(0, 100000, by = 1)
threshold_grid  <- 1e-6   # faster 
threshold_calib <- 1e-10  # stricter

# ---- Target metrics from literature ----
target_metrics <- list(
  prevalence_h   = 0.07,
  prevalence_c   = 0.015,
  incidence_h = 12 / 100000 / 365, # unités
  incidence_c = 18 / 100000 / 365, # unités
  recid_1     = 0.25,
  recid_2     = 0.50
)





###############################################################################
# 3. GRID SEARCH (1 -> beta, 2 -> sigma, 3 -> k)
###############################################################################

# ---- Grid 1: beta_h, beta_c (match carriage prevalence) ----
beta_grid <- grid_search(
  param_names = c("beta_h", "beta_c"),
  param_ranges = list(
    beta_h = c(0.001, 0.5, 25),   # from, to, number of points
    beta_c = c(0.001, 0.5, 25)
  ),
  metric_function = compute_metrics_beta,
  target_metrics  = list(prevalence_h = target_metrics$prevalence_h,
                         prevalence_c = target_metrics$prevalence_c),
  params_base = params,
  init_cond   = init_cond,
  time_vec    = time_vec_calib,
  n_cores     = 4,
  threshold = threshold_grid)

best_beta <- beta_grid$best_guess[, c("beta_h", "beta_c")]
params[c("beta_h", "beta_c")] <- as.numeric(best_beta)

p_beta <- plot_grid_search_beta(beta_grid,  target_metrics)
print(p_beta)


# ---- Grid 2: sigma_h, sigma_c (match incidence per N_tot) ----
sigma_grid <- grid_search(
  param_names = c("sigma_h", "sigma_c"),
  param_ranges = list(
    sigma_h = c(0.000001, 0.001, 25),
    sigma_c = c(0.000001, 0.001, 25)
  ),
  metric_function = compute_metrics_sigma,
  target_metrics  = list(incidence_h = target_metrics$incidence_h,
                         incidence_c = target_metrics$incidence_c),
  params_base = params,     # params already contains best beta from Grid 1
  init_cond   = init_cond,
  time_vec    = time_vec_calib,
  n_cores     = 4, 
  threshold = threshold_grid)

best_sigma <- sigma_grid$best_guess[, c("sigma_h", "sigma_c")]
params[c("sigma_h", "sigma_c")] <- as.numeric(best_sigma)

p_sigma <- plot_grid_search_sigma(sigma_grid, target_metrics)
print(p_sigma)


# ---- Grid 3: k_II, k_III (match recurrence prevalence) ----
k_grid <- grid_search(
  param_names = c("k_II", "k_III"),
  param_ranges = list(
    k_II  = c(1, 1000, 25),
    k_III = c(1, 1000, 25)
  ),
  metric_function = compute_metrics_k,
  target_metrics  = list(recid_1 = target_metrics$recid_1,
                         recid_2 = target_metrics$recid_2),
  params_base = params,     # params already contains best beta + best sigma
  init_cond   = init_cond,
  time_vec    = time_vec_calib,
  n_cores     = 4,
  threshold = threshold_grid)

best_k <- k_grid$best_guess[, c("k_II", "k_III")]
params[c("k_II", "k_III")] <- as.numeric(best_k)

p_k <- plot_grid_search_k(k_grid, target_metrics)
print(p_k)





###############################################################################
# 4. CALIBRATION (MULTI-START OPTIMIZATION)
###############################################################################

# Settings for multi-start Nelder-Mead
n_starts <- 10
n_cores  <- 4
maxit    <- 500

# Run calibration (optimizes beta_h, beta_c, sigma_h, sigma_c, k_II, k_III)
calib_res <- run_calibration(
  initial_params = params,         # start near the grid-search best values
  target_metrics = target_metrics, # prevalence_h/c, incidence_h/c, recid_1/2
  params_base    = params,         # baseline parameter vector (will be overwritten for the 6 calibrated params)
  init_cond      = init_cond,
  time_vec       = time_vec_calib,
  n_starts       = n_starts,
  n_cores        = n_cores,
  maxit          = maxit,
  threshold      = threshold_calib
)

# Build the final calibrated parameter vector
params_calib <- params
params_calib[names(calib_res$best$par)] <- calib_res$best$par

# Simulate once with calibrated parameters and compute metrics
out_calib <- run_calib_model_to_equilibrium(params_calib, init_cond, time_vec_calib, threshold_calib)
metrics_calib <- compute_metrics_calib(out_calib, params_calib)

# Quick check: print last values vs targets
cat("\n--- CALIBRATION CHECK (last time point) ---\n")
cat("prevalence_h:", tail(metrics_calib$carriage$prev_h, 1), " target:", target_metrics$prevalence_h, "\n")
cat("prevalence_c:", tail(metrics_calib$carriage$prev_c, 1), " target:", target_metrics$prevalence_c, "\n")
cat("incidence_h :", tail(metrics_calib$incidence_instant$inc_h_total_abs, 1), " target:", target_metrics$incidence_h, "\n")
cat("incidence_c :", tail(metrics_calib$incidence_instant$inc_c_total_abs, 1), " target:", target_metrics$incidence_c, "\n")
cat("recid_1     :", tail(metrics_calib$recurrence$rec1, 1), " target:", target_metrics$recid_1, "\n")
cat("recid_2     :", tail(metrics_calib$recurrence$rec2, 1), " target:", target_metrics$recid_2, "\n")





###############################################################################
# 5. RUN UNTIL EQUILIBRIUM WITH CALIBRATED PARAMETERS TO FIND INITIAL CONDITIONS
###############################################################################

# Run the calibration model with the final calibrated parameters (stop early if equilibrium)
out_eq <- run_calib_model_to_equilibrium(
  params   = params_calib,   # calibrated parameters
  init_cond = init_cond,     # initial conditions
  time_vec = time_vec_eq,        # long horizon, but will stop early if equilibrium reached
  threshold = threshold_calib)

# Get the last state (equilibrium state) as a named numeric vector : use as initial conditions for scenario runs !!
state_eq <- get_last_state(out_eq)

# Dynamics plots (hospital + community side by side)
dyn_plots <- plot_dynamics(
  ode_result = out_eq,
  targets = list(prevalence_h = target_metrics$prevalence_h, prevalence_c = target_metrics$prevalence_c),
  N_h = N_h,
  N_c = N_c)

print(dyn_plots$both)  # side-by-side plot

# Alpha plot
p_alpha <- plot_alpha_dynamics(out_eq, params_calib)
print(p_alpha)

# SAVE EQUILIBRIUM RUN OUTPUTS AS RDS (Bundle everything in a single object)
res_eq <- list(
  out_eq = out_eq,
  state_eq = state_eq,
  params_calib = params_calib,
  threshold_calib = threshold_calib,
  target_metrics = target_metrics,
  N_h = N_h,
  N_c = N_c
)
saveRDS(res_eq, "results_equilibrium.rds")





















