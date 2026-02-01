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
source("5_plot.R")

###############################################################################
# 0. PARAMETRES & POPULATIONS
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
delta_I   <- 1 / 15
delta_II  <- 1 / 20
delta_III <- 1 / 25

# ---- Severity weights for community → hospital flows ----
w     <- 1
w_I   <- 2
w_II  <- 4
w_III <- 8

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

# ---- Target metrics from literature ----
targets <- list(
  portage_h   = 0.07,
  portage_c   = 0.015,
  incidence_h = 12 / 100000 / 365, # unités
  incidence_c = 18 / 100000 / 365, # unités
  recid_1     = 0.25,
  recid_2     = 0.50
)

###############################################################################
# 1. Faire tourner modèle jusqu'à l'équilibre
###############################################################################

# Conditions initiales (pre-calibration)
# attention les conditions initiales doivent 
y0 <- create_initial_conditions_precalibration(N_h, N_c, prev_primo = 0.01, prev_rec = 0.005)

# Grille de temps pour le burn-in
t_burnin <- seq(0, 10000, by = 1)

# Résolution ODE
sol_burnin <- lsoda(y = y0,
                    times = t_burnin,
                    func  = cdiff_model_for_calibration,
                    parms = params)

# Extraire l'état final comme named vector
pop_eq <- setNames(sol_burnin[nrow(sol_burnin), -1], names(y0))

# Vérifier l'équilibre au dernier timestep
res_eq <- check_equilibrium(pop_eq, params, cdiff_model_for_calibration, threshold = 1e-6)

if (!res_eq$at_equilibrium) {
  cat("Pas d'équilibre sur la fenêtre du burn-in\n")
  cat("Compartiment le plus loin :", res_eq$worst_comp, " |dX/dt| =", res_eq$max_deriv, "\n")
  print(head(res_eq$details, 5))
} else {
  # Chercher le premier timestep où l'équilibre est atteint (binary search)
  lo <- 1; hi <- nrow(sol_burnin)
  while (lo < hi) {
    mid <- (lo + hi) %/% 2
    pop_mid <- setNames(sol_burnin[mid, -1], names(y0))
    if (check_equilibrium(pop_mid, params, cdiff_model_for_calibration, threshold = 1e-6)$at_equilibrium) {
      hi <- mid
    } else {
      lo <- mid + 1
    }
  }
  t_eq <- sol_burnin[lo, 1]
  cat("Équilibre atteint à t =", t_eq, "jours (max |dX/dt| =", res_eq$max_deriv, ")\n")
}

# Plot de la dynamique 
plots <- plot_dynamics(sol_dyn_df, targets, N_h, N_c)














###############################################################################
# 4. TABLEAU RÉSUMÉ : DÉBUT vs FIN DE LA DYNAMIQUE
###############################################################################
###############################################################################
# 4. TABLEAU RÉSUMÉ : DÉBUT vs FIN DE LA DYNAMIQUE
###############################################################################
###############################################################################
# 4. TABLEAU RÉSUMÉ : DÉBUT vs FIN DE LA DYNAMIQUE
###############################################################################
###############################################################################
# 4. TABLEAU RÉSUMÉ : DÉBUT vs FIN DE LA DYNAMIQUE
###############################################################################
###############################################################################
# 4. TABLEAU RÉSUMÉ : DÉBUT vs FIN DE LA DYNAMIQUE
###############################################################################
###############################################################################
# 4. TABLEAU RÉSUMÉ : DÉBUT vs FIN DE LA DYNAMIQUE
###############################################################################
###############################################################################
# 4. TABLEAU RÉSUMÉ : DÉBUT vs FIN DE LA DYNAMIQUE
###############################################################################

plots <- plot_dynamics(sol_dyn_df, targets, N_h, N_c)

# Afficher les 4 plots (2x2)
gridExtra::grid.arrange(plots$hospital_all,
                        plots$hospital_agg,
                        plots$community_all,
                        plots$community_agg,
                        nrow = 2, ncol = 2)


# faire un tableau plot dans lequel il y a le nombre de gens dans tous les compartiments

# main 

# Create initial conditions from a calibrated equilibrium state for intervention scenarios
# prendre les dernières valeurs de la résolution de l'équation lsoa



