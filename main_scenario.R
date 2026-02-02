###############################################################################
############################ MAIN_SCENARIO ####################################
###############################################################################

###############################################################################
# 0. SET UP
###############################################################################

# Clear workspace
rm(list = ls())

# Packages & sources
source("0_packages.R")
source("1_model.R")
source("2_function.R")   
source("5_plot.R")       

# Load calibrated results
res <- readRDS("results_calibration.rds")

# Build scenario parameter vector (ensure numeric)
params <- c(
  as.numeric(res$params),
  alpha     = as.numeric(res$alpha_eq$alpha),
  alpha_I   = as.numeric(res$alpha_eq$alpha_I),
  alpha_II  = as.numeric(res$alpha_eq$alpha_II),
  alpha_III = as.numeric(res$alpha_eq$alpha_III)
)

names(params) <- c(names(res$params), "alpha", "alpha_I", "alpha_II", "alpha_III")

print(params)

# Equilibrium state to be used as initial conditions for scenarios
state_eq <- res$state_eq

# 
time <- seq(0, 365*5, by = 1)


###############################################################################
# 1. TEST VALIDITY OF THE STRUCTURE 
###############################################################################

# Paramètres : VE = 0 (pas d'effet vaccinal)
params_test <- c(
  params,
  tau_mult_red = 1, # 
  sigma_mult_v = 1  # VE = 0 sur infection
)

# ============================================================================
# SCÉNARIO TEST 1 : VC = 0, VE = 0
# ============================================================================

init_test1 <- c(
  # Hospital non-vaccinated
  S0_h_nv = as.numeric(state_eq["S0_h"]), SA_h_nv = as.numeric(state_eq["SA_h"]),
  C0_h_nv = as.numeric(state_eq["C0_h"]), CA_h_nv = as.numeric(state_eq["CA_h"]),
  I_h_nv = as.numeric(state_eq["I_h"]),
  S_II_h_nv = as.numeric(state_eq["S_II_h"]), C_II_h_nv = as.numeric(state_eq["C_II_h"]),
  I_II_h_nv = as.numeric(state_eq["I_II_h"]),
  S_III_h_nv = as.numeric(state_eq["S_III_h"]), C_III_h_nv = as.numeric(state_eq["C_III_h"]),
  I_III_h_nv = as.numeric(state_eq["I_III_h"]),
  
  # Community non-vaccinated
  S0_c_nv = as.numeric(state_eq["S0_c"]), SA_c_nv = as.numeric(state_eq["SA_c"]),
  C0_c_nv = as.numeric(state_eq["C0_c"]), CA_c_nv = as.numeric(state_eq["CA_c"]),
  I_c_nv = as.numeric(state_eq["I_c"]),
  S_II_c_nv = as.numeric(state_eq["S_II_c"]), C_II_c_nv = as.numeric(state_eq["C_II_c"]),
  I_II_c_nv = as.numeric(state_eq["I_II_c"]),
  S_III_c_nv = as.numeric(state_eq["S_III_c"]), C_III_c_nv = as.numeric(state_eq["C_III_c"]),
  I_III_c_nv = as.numeric(state_eq["I_III_c"]),
  
  # Hospital vaccinated (initialement 0)
  S0_h_v = 0, SA_h_v = 0, C0_h_v = 0, CA_h_v = 0, I_h_v = 0,
  S_II_h_v = 0, C_II_h_v = 0, I_II_h_v = 0,
  S_III_h_v = 0, C_III_h_v = 0, I_III_h_v = 0,
  
  # Community vaccinated (initialement 0)
  S0_c_v = 0, SA_c_v = 0, C0_c_v = 0, CA_c_v = 0, I_c_v = 0,
  S_II_c_v = 0, C_II_c_v = 0, I_II_c_v = 0,
  S_III_c_v = 0, C_III_c_v = 0, I_III_c_v = 0
)

out_test1 <- lsoda(y = init_test0, times = time, func = cdiff_model, parms = params_test)
metrics_test1 <- compute_metrics_vacc(out_test1, params_test)

# ============================================================================
# SCÉNARIO TEST 2 : VC = 100, VE = 0
# ============================================================================

init_test2 <- c(
  # Hospital non-vaccinated (initialement 0)
  S0_h_nv = 0, SA_h_nv = 0, C0_h_nv = 0, CA_h_nv = 0, I_h_nv = 0,
  S_II_h_nv = 0, C_II_h_nv = 0, I_II_h_nv = 0,
  S_III_h_nv = 0, C_III_h_nv = 0, I_III_h_nv = 0,
  
  # Community non-vaccinated (initialement 0)
  S0_c_nv = 0, SA_c_nv = 0, C0_c_nv = 0, CA_c_nv = 0, I_c_nv = 0,
  S_II_c_nv = 0, C_II_c_nv = 0, I_II_c_nv = 0,
  S_III_c_nv = 0, C_III_c_nv = 0, I_III_c_nv = 0,
  
  # Hospital vaccinated (tout le monde)
  S0_h_v = as.numeric(state_eq["S0_h"]), SA_h_v = as.numeric(state_eq["SA_h"]),
  C0_h_v = as.numeric(state_eq["C0_h"]), CA_h_v = as.numeric(state_eq["CA_h"]),
  I_h_v = as.numeric(state_eq["I_h"]),
  S_II_h_v = as.numeric(state_eq["S_II_h"]), C_II_h_v = as.numeric(state_eq["C_II_h"]),
  I_II_h_v = as.numeric(state_eq["I_II_h"]),
  S_III_h_v = as.numeric(state_eq["S_III_h"]), C_III_h_v = as.numeric(state_eq["C_III_h"]),
  I_III_h_v = as.numeric(state_eq["I_III_h"]),
  
  # Community vaccinated (tout le monde)
  S0_c_v = as.numeric(state_eq["S0_c"]), SA_c_v = as.numeric(state_eq["SA_c"]),
  C0_c_v = as.numeric(state_eq["C0_c"]), CA_c_v = as.numeric(state_eq["CA_c"]),
  I_c_v = as.numeric(state_eq["I_c"]),
  S_II_c_v = as.numeric(state_eq["S_II_c"]), C_II_c_v = as.numeric(state_eq["C_II_c"]),
  I_II_c_v = as.numeric(state_eq["I_II_c"]),
  S_III_c_v = as.numeric(state_eq["S_III_c"]), C_III_c_v = as.numeric(state_eq["C_III_c"]),
  I_III_c_v = as.numeric(state_eq["I_III_c"])
)

out_test2 <- lsoda(y = init_test2, times = time, func = cdiff_model, parms = params_test)
metrics_test2 <- compute_metrics_vacc(out_test2, params_test)

###############################################################################
# 2. ATB  
###############################################################################


###############################################################################
# 3. VACCINAION  
###############################################################################


