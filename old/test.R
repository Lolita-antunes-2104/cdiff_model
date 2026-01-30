###############################################################################
# QUICK TEST - SCENARIOS ATB & VACCINATION
###############################################################################

# Clear workspace
rm(list = ls())

# Load required functions (assuming you have sourced them)
# source("1_model.R")
# source("2_functions.R")
# ... etc

###############################################################################
# 1. FAKE CALIBRATED PARAMETERS & EQUILIBRIUM STATE
###############################################################################

# Population sizes
N_h <- 100
N_c <- 50000

# Fake calibrated parameters
params_fixed <- c(
  beta_h  = 0.08,
  beta_c  = 0.025,
  sigma_h = 0.005,
  sigma_c = 0.002,
  tau_h   = 0.035,
  tau_c   = 0.0027,
  omega_h = 0.05,
  omega_c = 0.022,
  nu_h    = 7,
  nu_c    = 7,
  gamma_h = 0.012,
  gamma_c = 0.012,
  epsilon_h = 0.05,
  epsilon_c = 0.06,
  p_h     = 0.5,
  p_c     = 0.5,
  phi_h   = 0.018,
  phi_c   = 0.018,
  k_A     = 4,
  k_II    = 2.5,
  k_III   = 3.5,
  w       = 1,
  w_I     = 1.5,
  w_II    = 2,
  w_III   = 2.5,
  delta   = 1/7.5,
  delta_I = 1/16.6,
  delta_II = 1/20,
  delta_III = 1/25,
  prop    = 0.5,
  alpha   = 0.0002,
  alpha_I = 0.0003,
  alpha_II = 0.0004,
  alpha_III = 0.0005
)

# Fake equilibrium state (invented compartment values)
last_eq <- c(
  # Hospital
  S0_h = 65, SA_h = 10, C0_h = 5, CA_h = 3, I_h = 2,
  S_II_h = 4, C_II_h = 3, I_II_h = 2,
  S_III_h = 3, C_III_h = 2, I_III_h = 1,
  CumI_primo_h = 0, CumI_rec_h = 0, CumI_h = 0,
  
  # Community
  S0_c = 42000, SA_c = 3000, C0_c = 2500, CA_c = 1000, I_c = 500,
  S_II_c = 300, C_II_c = 250, I_II_c = 150,
  S_III_c = 150, C_III_c = 100, I_III_c = 50,
  CumI_primo_c = 0, CumI_rec_c = 0, CumI_c = 0
)

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

cat("\n=== RUNNING BASELINE (ATB) ===\n")

# Baseline with fixed alpha
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

cat("Baseline computed.\n")

# ATB reduction scenarios
atb_sets <- list(
  set1 = c(0.00, 0.01, 1.00),
  set2 = c(0.05, 0.10, 0.20)
)

atb_results <- list()

cat("\n=== RUNNING ATB SCENARIOS ===\n")

for (set in names(atb_sets)) {
  for (red in atb_sets[[set]]) {
    
    cat(sprintf("Running ATB scenario: %s, reduction = %.2f\n", set, red))
    
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

cat("\n=== ATB RESULTS ===\n")
print(atb_table)

###############################################################################
# 4. VACCINATION SCENARIOS
###############################################################################

VC_vals <- c(0.25, 0.50, 0.75)
VE_vals <- c(0.25, 0.50, 0.75)

vacc_results <- list()

cat("\n=== RUNNING VACCINATION SCENARIOS ===\n")

for (VC in VC_vals) {
  
  cat(sprintf("\n--- Vaccination Coverage = %.0f%% ---\n", VC * 100))
  
  # Initial conditions with fixed vaccination coverage
  init_vacc <- create_initial_conditions_vaccination(
    last_state_without_vacc = last_eq,
    coverage = VC
  )
  
  # Baseline vaccinated (VE = 0)
  cat(sprintf("  Running baseline (VE = 0)...\n"))
  
  ode_base_vacc <- run_model_to_equilibrium(
    params_fixed,
    init_vacc,
    times,
    model_fn = cdiff_hc_vacc_model,
    alpha_mode = "fixed",
    VE = 0,
    vacc_type = "both"
  )
  
  last_base_vacc <- extract_last_state(ode_base_vacc)
  base_vacc_summary <- summarize_state(last_base_vacc, params_fixed)
  
  # Vaccination scenarios (VE > 0)
  for (VE in VE_vals) {
    
    cat(sprintf("  Running VE = %.0f%%...\n", VE * 100))
    
    ode <- run_model_to_equilibrium(
      params_fixed,
      init_vacc,
      times,
      model_fn = cdiff_hc_vacc_model,
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

cat("\n=== VACCINATION RESULTS ===\n")
print(vacc_table)

cat("\n=== TEST COMPLETED ===\n")