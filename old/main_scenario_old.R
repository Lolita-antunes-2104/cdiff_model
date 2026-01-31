# Clear workspace
rm(list = ls())

# Packages & sources 
source("0_packages.R")
source("11_model_vac.R")

# Load calibrated results 
res <- readRDS("results_calibration.rds")

# Préparer les paramètres
params <- c(
  res$params,
  alpha     = as.numeric(res$alpha_eq$alpha),
  alpha_I   = as.numeric(res$alpha_eq$alpha_I),
  alpha_II  = as.numeric(res$alpha_eq$alpha_II),
  alpha_III = as.numeric(res$alpha_eq$alpha_III)
)

state_eq <- res$state_eq

# ============================================================================
# SCÉNARIO TEST 0 : Modèle sans vaccination pendant 5 ans
# ============================================================================

init_test0 <- c(
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

# Paramètres : VE = 0 (pas d'effet vaccinal)
params_test0 <- c(
  params,
  sigma_mult_v = 1,  # VE = 0 sur infection
  gamma_mult_v = 1   # VE = 0 sur portage
)

time <- seq(0, 365*5, by = 1)
out_test0 <- lsoda(y = init_test0, times = time, func = cdiff_model, parms = params_test0)
metrics_test0 <- compute_metrics_vacc(out_test0, params_test0)

# ============================================================================
# SCÉNARIO TEST 1 : 100% non-vaccinés, VE = 0
# ============================================================================

init_test1 <- init_test0  # Même initialisation
params_test1 <- params_test0  # Mêmes paramètres

out_test1 <- lsoda(y = init_test1, times = time, func = cdiff_model, parms = params_test1)
metrics_test1 <- compute_metrics_vacc(out_test1, params_test1)

# ============================================================================
# SCÉNARIO TEST 2 : 100% vaccinés, VE = 0
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

params_test2 <- params_test0  # VE = 0

out_test2 <- lsoda(y = init_test2, times = time, func = cdiff_model, parms = params_test2)
metrics_test2 <- compute_metrics_vacc(out_test2, params_test2)

# ============================================================================
# CRÉER TABLEAUX COMPARATIFS SÉPARÉS
# ============================================================================

# --- Clean tiny numerical noise (from ODE solver) ---
clean_mat <- function(x, eps = 1e-12) {
  x2 <- x
  x2[abs(x2) < eps] <- 0
  x2
}

# Fonction pour extraire les compartiments
extraire_compartiments <- function(out, label) {
  
  out_clean <- clean_mat(out, eps = 1e-12)
  
  t0 <- 1  # Première ligne
  tf <- nrow(out_clean)  # Dernière ligne
  
  # Populations totales
  N_h_nv_0 <- sum(out_clean[t0, c("S0_h_nv", "SA_h_nv", "C0_h_nv", "CA_h_nv", "I_h_nv",
                            "S_II_h_nv", "C_II_h_nv", "I_II_h_nv",
                            "S_III_h_nv", "C_III_h_nv", "I_III_h_nv")])
  N_c_nv_0 <- sum(out_clean[t0, c("S0_c_nv", "SA_c_nv", "C0_c_nv", "CA_c_nv", "I_c_nv",
                            "S_II_c_nv", "C_II_c_nv", "I_II_c_nv",
                            "S_III_c_nv", "C_III_c_nv", "I_III_c_nv")])
  N_h_v_0 <- sum(out_clean[t0, c("S0_h_v", "SA_h_v", "C0_h_v", "CA_h_v", "I_h_v",
                           "S_II_h_v", "C_II_h_v", "I_II_h_v",
                           "S_III_h_v", "C_III_h_v", "I_III_h_v")])
  N_c_v_0 <- sum(out_clean[t0, c("S0_c_v", "SA_c_v", "C0_c_v", "CA_c_v", "I_c_v",
                           "S_II_c_v", "C_II_c_v", "I_II_c_v",
                           "S_III_c_v", "C_III_c_v", "I_III_c_v")])
  
  N_h_nv_f <- sum(out_clean[tf, c("S0_h_nv", "SA_h_nv", "C0_h_nv", "CA_h_nv", "I_h_nv",
                            "S_II_h_nv", "C_II_h_nv", "I_II_h_nv",
                            "S_III_h_nv", "C_III_h_nv", "I_III_h_nv")])
  N_c_nv_f <- sum(out_clean[tf, c("S0_c_nv", "SA_c_nv", "C0_c_nv", "CA_c_nv", "I_c_nv",
                            "S_II_c_nv", "C_II_c_nv", "I_II_c_nv",
                            "S_III_c_nv", "C_III_c_nv", "I_III_c_nv")])
  N_h_v_f <- sum(out_clean[tf, c("S0_h_v", "SA_h_v", "C0_h_v", "CA_h_v", "I_h_v",
                           "S_II_h_v", "C_II_h_v", "I_II_h_v",
                           "S_III_h_v", "C_III_h_v", "I_III_h_v")])
  N_c_v_f <- sum(out_clean[tf, c("S0_c_v", "SA_c_v", "C0_c_v", "CA_c_v", "I_c_v",
                           "S_II_c_v", "C_II_c_v", "I_II_c_v",
                           "S_III_c_v", "C_III_c_v", "I_III_c_v")])
  
  # Porteurs C
  C_h_nv_0 <- sum(out_clean[t0, c("C0_h_nv", "CA_h_nv", "C_II_h_nv", "C_III_h_nv")])
  C_c_nv_0 <- sum(out_clean[t0, c("C0_c_nv", "CA_c_nv", "C_II_c_nv", "C_III_c_nv")])
  C_h_v_0 <- sum(out_clean[t0, c("C0_h_v", "CA_h_v", "C_II_h_v", "C_III_h_v")])
  C_c_v_0 <- sum(out_clean[t0, c("C0_c_v", "CA_c_v", "C_II_c_v", "C_III_c_v")])
  
  C_h_nv_f <- sum(out_clean[tf, c("C0_h_nv", "CA_h_nv", "C_II_h_nv", "C_III_h_nv")])
  C_c_nv_f <- sum(out_clean[tf, c("C0_c_nv", "CA_c_nv", "C_II_c_nv", "C_III_c_nv")])
  C_h_v_f <- sum(out_clean[tf, c("C0_h_v", "CA_h_v", "C_II_h_v", "C_III_h_v")])
  C_c_v_f <- sum(out_clean[tf, c("C0_c_v", "CA_c_v", "C_II_c_v", "C_III_c_v")])
  
  # Infectés I
  I_h_nv_0 <- sum(out_clean[t0, c("I_h_nv", "I_II_h_nv", "I_III_h_nv")])
  I_c_nv_0 <- sum(out_clean[t0, c("I_c_nv", "I_II_c_nv", "I_III_c_nv")])
  I_h_v_0 <- sum(out_clean[t0, c("I_h_v", "I_II_h_v", "I_III_h_v")])
  I_c_v_0 <- sum(out_clean[t0, c("I_c_v", "I_II_c_v", "I_III_c_v")])
  
  I_h_nv_f <- sum(out_clean[tf, c("I_h_nv", "I_II_h_nv", "I_III_h_nv")])
  I_c_nv_f <- sum(out_clean[tf, c("I_c_nv", "I_II_c_nv", "I_III_c_nv")])
  I_h_v_f <- sum(out_clean[tf, c("I_h_v", "I_II_h_v", "I_III_h_v")])
  I_c_v_f <- sum(out_clean[tf, c("I_c_v", "I_II_c_v", "I_III_c_v")])
  
  data.frame(
    Scenario = label,
    Temps = c("t=0", "t=5ans"),
    N_h_nv = c(N_h_nv_0, N_h_nv_f),
    N_c_nv = c(N_c_nv_0, N_c_nv_f),
    N_h_v = c(N_h_v_0, N_h_v_f),
    N_c_v = c(N_c_v_0, N_c_v_f),
    C_h_nv = c(C_h_nv_0, C_h_nv_f),
    C_c_nv = c(C_c_nv_0, C_c_nv_f),
    C_h_v = c(C_h_v_0, C_h_v_f),
    C_c_v = c(C_c_v_0, C_c_v_f),
    I_h_nv = c(I_h_nv_0, I_h_nv_f),
    I_c_nv = c(I_c_nv_0, I_c_nv_f),
    I_h_v = c(I_h_v_0, I_h_v_f),
    I_c_v = c(I_c_v_0, I_c_v_f)
  )
}

# Fonction pour extraire les outputs (prévalences et incidences)
extraire_outputs <- function(metrics, label) {
  t0 <- 1  # Première ligne
  tf <- nrow(metrics$carriage)  # Dernière ligne
  
  # Prévalences
  prev_h_nv_0 <- metrics$carriage$prev_h_nv[t0]
  prev_c_nv_0 <- metrics$carriage$prev_c_nv[t0]
  prev_h_v_0 <- metrics$carriage$prev_h_v[t0]
  prev_c_v_0 <- metrics$carriage$prev_c_v[t0]
  
  prev_h_nv_f <- metrics$carriage$prev_h_nv[tf]
  prev_c_nv_f <- metrics$carriage$prev_c_nv[tf]
  prev_h_v_f <- metrics$carriage$prev_h_v[tf]
  prev_c_v_f <- metrics$carriage$prev_c_v[tf]
  
  # Incidences
  inc_h_nv_0 <- metrics$incidence$inc_h_nv_n[t0]
  inc_c_nv_0 <- metrics$incidence$inc_c_nv_n[t0]
  inc_h_v_0 <- metrics$incidence$inc_h_v_n[t0]
  inc_c_v_0 <- metrics$incidence$inc_c_v_n[t0]
  
  inc_h_nv_f <- metrics$incidence$inc_h_nv_n[tf]
  inc_c_nv_f <- metrics$incidence$inc_c_nv_n[tf]
  inc_h_v_f <- metrics$incidence$inc_h_v_n[tf]
  inc_c_v_f <- metrics$incidence$inc_c_v_n[tf]
  
  data.frame(
    Scenario = label,
    Temps = c("t=0", "t=5ans"),
    Prev_h_nv = c(prev_h_nv_0, prev_h_nv_f),
    Prev_c_nv = c(prev_c_nv_0, prev_c_nv_f),
    Prev_h_v = c(prev_h_v_0, prev_h_v_f),
    Prev_c_v = c(prev_c_v_0, prev_c_v_f),
    Inc_h_nv = c(inc_h_nv_0, inc_h_nv_f),
    Inc_c_nv = c(inc_c_nv_0, inc_c_nv_f),
    Inc_h_v = c(inc_h_v_0, inc_h_v_f),
    Inc_c_v = c(inc_c_v_0, inc_c_v_f)
  )
}

# Créer les tableaux séparés
tab_comp_0 <- extraire_compartiments(out_test0, "Test 0: Sans vaccination")
tab_comp_1 <- extraire_compartiments(out_test1, "Test 1: 100% non-vaccinés, VE=0")
tab_comp_2 <- extraire_compartiments(out_test2, "Test 2: 100% vaccinés, VE=0")

tab_out_0 <- extraire_outputs(metrics_test0, "Test 0: Sans vaccination")
tab_out_1 <- extraire_outputs(metrics_test1, "Test 1: 100% non-vaccinés, VE=0")
tab_out_2 <- extraire_outputs(metrics_test2, "Test 2: 100% vaccinés, VE=0")

# Tableaux finaux
tableau_compartiments <- rbind(tab_comp_0, tab_comp_1, tab_comp_2)
tableau_outputs <- rbind(tab_out_0, tab_out_1, tab_out_2)



# ============================================================================
# TABLEAU COMPARTIMENTS ET OUTPUTS
# Créer les tableaux sous forme de plots
library(gridExtra)
library(grid)
library(dplyr)

options(scipen = 999)  # évite l'affichage en 5e+04

tab_comp_format <- tableau_compartiments %>%
  mutate(across(where(is.numeric), ~format(signif(., 3), scientific = FALSE, trim = TRUE)))

tab_out_format <- tableau_outputs %>%
  mutate(
    Prev_h_nv = 100 * Prev_h_nv,
    Prev_c_nv = 100 * Prev_c_nv,
    Prev_h_v  = 100 * Prev_h_v,
    Prev_c_v  = 100 * Prev_c_v,
    Inc_h_nv  = Inc_h_nv * 365 * 1e5,
    Inc_c_nv  = Inc_c_nv * 365 * 1e5,
    Inc_h_v   = Inc_h_v  * 365 * 1e5,
    Inc_c_v   = Inc_c_v  * 365 * 1e5
  ) %>%
  mutate(across(where(is.numeric), ~format(signif(., 3), scientific = FALSE, trim = TRUE)))

# Créer les grobs pour chaque tableau
grob_comp <- tableGrob(tab_comp_format, rows = NULL,
                       theme = ttheme_default(
                         base_size = 8,
                         core = list(fg_params = list(hjust = 1, x = 0.95)),
                         colhead = list(fg_params = list(fontface = "bold"))
                       ))

grob_out <- tableGrob(tab_out_format, rows = NULL,
                      theme = ttheme_default(
                        base_size = 8,
                        core = list(fg_params = list(hjust = 1, x = 0.95)),
                        colhead = list(fg_params = list(fontface = "bold"))
                      ))

# Ajouter des titres
titre_comp <- textGrob("COMPARTIMENTS (Populations N, Porteurs C, Infectés I)", 
                       gp = gpar(fontsize = 12, fontface = "bold"))
titre_out <- textGrob("OUTPUTS (Prévalences en %, Incidences pour 100,000 PA)", 
                      gp = gpar(fontsize = 12, fontface = "bold"))

# AFFICHER DIRECTEMENT DANS PLOTS (sans sauvegarder)
grid.arrange(
  titre_comp,
  grob_comp,
  titre_out,
  grob_out,
  ncol = 1,
  heights = c(0.5, 4, 0.5, 4)
)















###############################################################################
###############################################################################
# 27 vaccination scenarios -> 3 "plots" (each plot = 2 tables), 9 rows each
# Assumes you already did:
#   source("0_packages.R")
#   source("11_model_vac.R")
#   res <- readRDS("results_calibration.rds")
#   params <- c(res$params, alpha=..., alpha_I=..., alpha_II=..., alpha_III=...)
#   state_eq <- res$state_eq
###############################################################################

# ---- User inputs ----
years <- 5
time  <- seq(0, 365*years, by = 1)

# ---- Helpers ----

# Split equilibrium state into vaccinated / non-vaccinated at t=0 using coverage cov in [0,1]
build_init_from_eq <- function(state_eq, cov) {
  split_nv <- function(x) (1 - cov) * as.numeric(x)
  split_v  <- function(x) cov * as.numeric(x)
  
  c(
    # Hospital non-vaccinated
    S0_h_nv    = split_nv(state_eq["S0_h"]),    SA_h_nv    = split_nv(state_eq["SA_h"]),
    C0_h_nv    = split_nv(state_eq["C0_h"]),    CA_h_nv    = split_nv(state_eq["CA_h"]),
    I_h_nv     = split_nv(state_eq["I_h"]),
    S_II_h_nv  = split_nv(state_eq["S_II_h"]),  C_II_h_nv  = split_nv(state_eq["C_II_h"]),
    I_II_h_nv  = split_nv(state_eq["I_II_h"]),
    S_III_h_nv = split_nv(state_eq["S_III_h"]), C_III_h_nv = split_nv(state_eq["C_III_h"]),
    I_III_h_nv = split_nv(state_eq["I_III_h"]),
    
    # Community non-vaccinated
    S0_c_nv    = split_nv(state_eq["S0_c"]),    SA_c_nv    = split_nv(state_eq["SA_c"]),
    C0_c_nv    = split_nv(state_eq["C0_c"]),    CA_c_nv    = split_nv(state_eq["CA_c"]),
    I_c_nv     = split_nv(state_eq["I_c"]),
    S_II_c_nv  = split_nv(state_eq["S_II_c"]),  C_II_c_nv  = split_nv(state_eq["C_II_c"]),
    I_II_c_nv  = split_nv(state_eq["I_II_c"]),
    S_III_c_nv = split_nv(state_eq["S_III_c"]), C_III_c_nv = split_nv(state_eq["C_III_c"]),
    I_III_c_nv = split_nv(state_eq["I_III_c"]),
    
    # Hospital vaccinated
    S0_h_v     = split_v(state_eq["S0_h"]),     SA_h_v     = split_v(state_eq["SA_h"]),
    C0_h_v     = split_v(state_eq["C0_h"]),     CA_h_v     = split_v(state_eq["CA_h"]),
    I_h_v      = split_v(state_eq["I_h"]),
    S_II_h_v   = split_v(state_eq["S_II_h"]),   C_II_h_v   = split_v(state_eq["C_II_h"]),
    I_II_h_v   = split_v(state_eq["I_II_h"]),
    S_III_h_v  = split_v(state_eq["S_III_h"]),  C_III_h_v  = split_v(state_eq["C_III_h"]),
    I_III_h_v  = split_v(state_eq["I_III_h"]),
    
    # Community vaccinated
    S0_c_v     = split_v(state_eq["S0_c"]),     SA_c_v     = split_v(state_eq["SA_c"]),
    C0_c_v     = split_v(state_eq["C0_c"]),     CA_c_v     = split_v(state_eq["CA_c"]),
    I_c_v      = split_v(state_eq["I_c"]),
    S_II_c_v   = split_v(state_eq["S_II_c"]),   C_II_c_v   = split_v(state_eq["C_II_c"]),
    I_II_c_v   = split_v(state_eq["I_II_c"]),
    S_III_c_v  = split_v(state_eq["S_III_c"]),  C_III_c_v  = split_v(state_eq["C_III_c"]),
    I_III_c_v  = split_v(state_eq["I_III_c"])
  )
}

clean_mat <- function(x, eps = 1e-12) {
  x2 <- x
  x2[abs(x2) < eps] <- 0
  x2
}

# Extract compartments at final time (one row = one scenario)
extract_compartments_tf <- function(out, label) {
  outc <- clean_mat(out, eps = 1e-12)
  tf <- nrow(outc)
  
  # Totals
  N_h_nv <- sum(outc[tf, c("S0_h_nv","SA_h_nv","C0_h_nv","CA_h_nv","I_h_nv","S_II_h_nv","C_II_h_nv","I_II_h_nv","S_III_h_nv","C_III_h_nv","I_III_h_nv")])
  N_c_nv <- sum(outc[tf, c("S0_c_nv","SA_c_nv","C0_c_nv","CA_c_nv","I_c_nv","S_II_c_nv","C_II_c_nv","I_II_c_nv","S_III_c_nv","C_III_c_nv","I_III_c_nv")])
  N_h_v  <- sum(outc[tf, c("S0_h_v","SA_h_v","C0_h_v","CA_h_v","I_h_v","S_II_h_v","C_II_h_v","I_II_h_v","S_III_h_v","C_III_h_v","I_III_h_v")])
  N_c_v  <- sum(outc[tf, c("S0_c_v","SA_c_v","C0_c_v","CA_c_v","I_c_v","S_II_c_v","C_II_c_v","I_II_c_v","S_III_c_v","C_III_c_v","I_III_c_v")])
  
  # Carriers (C)
  C_h_nv <- sum(outc[tf, c("C0_h_nv","CA_h_nv","C_II_h_nv","C_III_h_nv")])
  C_c_nv <- sum(outc[tf, c("C0_c_nv","CA_c_nv","C_II_c_nv","C_III_c_nv")])
  C_h_v  <- sum(outc[tf, c("C0_h_v","CA_h_v","C_II_h_v","C_III_h_v")])
  C_c_v  <- sum(outc[tf, c("C0_c_v","CA_c_v","C_II_c_v","C_III_c_v")])
  
  # Infected (I)
  I_h_nv <- sum(outc[tf, c("I_h_nv","I_II_h_nv","I_III_h_nv")])
  I_c_nv <- sum(outc[tf, c("I_c_nv","I_II_c_nv","I_III_c_nv")])
  I_h_v  <- sum(outc[tf, c("I_h_v","I_II_h_v","I_III_h_v")])
  I_c_v  <- sum(outc[tf, c("I_c_v","I_II_c_v","I_III_c_v")])
  
  data.frame(
    Scenario = label,
    N_h_nv = N_h_nv, N_c_nv = N_c_nv, N_h_v = N_h_v, N_c_v = N_c_v,
    C_h_nv = C_h_nv, C_c_nv = C_c_nv, C_h_v = C_h_v, C_c_v = C_c_v,
    I_h_nv = I_h_nv, I_c_nv = I_c_nv, I_h_v = I_h_v, I_c_v = I_c_v
  )
}

# Extract outputs at final time (prevalence in %, incidence per 100k PY)
extract_outputs_tf <- function(metrics, label) {
  tf <- nrow(metrics$carriage)
  
  prev_h_nv <- 100 * metrics$carriage$prev_h_nv[tf]
  prev_c_nv <- 100 * metrics$carriage$prev_c_nv[tf]
  prev_h_v  <- 100 * metrics$carriage$prev_h_v[tf]
  prev_c_v  <- 100 * metrics$carriage$prev_c_v[tf]
  
  # Annualised per 100k person-years (same transform you used)
  inc_h_nv <- metrics$incidence$inc_h_nv_n[tf] * 365 * 1e5
  inc_c_nv <- metrics$incidence$inc_c_nv_n[tf] * 365 * 1e5
  inc_h_v  <- metrics$incidence$inc_h_v_n[tf]  * 365 * 1e5
  inc_c_v  <- metrics$incidence$inc_c_v_n[tf]  * 365 * 1e5
  
  data.frame(
    Scenario = label,
    Prev_h_nv = prev_h_nv, Prev_c_nv = prev_c_nv, Prev_h_v = prev_h_v, Prev_c_v = prev_c_v,
    Inc_h_nv  = inc_h_nv,  Inc_c_nv  = inc_c_nv,  Inc_h_v  = inc_h_v,  Inc_c_v  = inc_c_v
  )
}

# Run one scenario
run_one_scenario <- function(params_base, state_eq, cov, VE, vacc_type = c("sigma","gamma","both"), atb_red_c = 0) {
  vacc_type <- match.arg(vacc_type)
  
  # Map VE to multipliers
  sigma_mult_v <- 1
  gamma_mult_v <- 1
  if (vacc_type == "sigma") {
    sigma_mult_v <- 1 - VE
    gamma_mult_v <- 1
  } else if (vacc_type == "gamma") {
    sigma_mult_v <- 1
    gamma_mult_v <- 1 + VE
  } else if (vacc_type == "both") {
    sigma_mult_v <- 1 - VE
    gamma_mult_v <- 1 + VE
  }
  
  # Initial state
  init <- build_init_from_eq(state_eq, cov = cov)
  
  # Scenario params: only tau_c reduced in community
  params_i <- c(
    params_base,
    tau_c        = as.numeric(params_base[["tau_c"]]) * (1 - atb_red_c),
    sigma_mult_v = sigma_mult_v,
    gamma_mult_v = gamma_mult_v
  )
  
  out <- lsoda(y = init, times = time, func = cdiff_model, parms = params_i)
  met <- compute_metrics_vacc(out, params_i)
  
  list(out = out, metrics = met, sigma_mult_v = sigma_mult_v, gamma_mult_v = gamma_mult_v)
}

# Make one "plot" = titles + 2 tables (9 rows)
plot_group_tables <- function(df_comp, df_out, group_title) {
  
  library(gridExtra)
  library(grid)
  library(dplyr)
  
  options(scipen = 999)
  
  # Format numbers
  dfc <- df_comp %>%
    mutate(across(where(is.numeric), ~format(signif(., 2), scientific = FALSE, trim = TRUE)))
  
  dfo <- df_out %>%
    mutate(across(where(is.numeric), ~format(signif(., 2), scientific = FALSE, trim = TRUE)))
  
  grob_comp <- tableGrob(dfc, rows = NULL,
                         theme = ttheme_default(
                           base_size = 8,
                           core = list(fg_params = list(hjust = 1, x = 0.95)),
                           colhead = list(fg_params = list(fontface = "bold"))
                         ))
  
  grob_out <- tableGrob(dfo, rows = NULL,
                        theme = ttheme_default(
                          base_size = 8,
                          core = list(fg_params = list(hjust = 1, x = 0.95)),
                          colhead = list(fg_params = list(fontface = "bold"))
                        ))
  
  titre_group <- textGrob(group_title, gp = gpar(fontsize = 13, fontface = "bold"))
  titre_comp  <- textGrob("COMPARTIMENTS (Populations N, Porteurs C, Infectés I) — t = 5 ans",
                          gp = gpar(fontsize = 11, fontface = "bold"))
  titre_out   <- textGrob("OUTPUTS (Prévalences en %, Incidences pour 100,000 PA) — t = 5 ans",
                          gp = gpar(fontsize = 11, fontface = "bold"))
  
  grid.arrange(
    titre_group,
    titre_comp, grob_comp,
    titre_out,  grob_out,
    ncol = 1,
    heights = c(0.6, 0.5, 3.8, 0.5, 3.2)
  )
}

###############################################################################
# DEFINE THE 3 GROUPS (each group = 9 scenarios)
###############################################################################

vacc_types <- c("sigma","gamma","both")

# ---- Group 1: vary VC (coverage) = 30/60/90, with VE = 40%, ATB red = 0% ----
VC_levels_g1 <- c(0.30, 0.60, 0.90)
VE_g1 <- 0.40
ATB_g1 <- 0.00

# ---- Group 2: vary VE = 20/40/60, with VC = 60%, ATB red = 0% ----
VE_levels_g2 <- c(0.20, 0.40, 0.60)
VC_g2 <- 0.60
ATB_g2 <- 0.00

# ---- Group 3: vary ATB red in community = 5/10/20, with VC = 60%, VE = 40% ----
ATB_levels_g3 <- c(0.05, 0.10, 0.20)
VC_g3 <- 0.60
VE_g3 <- 0.40

###############################################################################
# RUN ALL 27 SCENARIOS, BUILD 3x (comp table + out table) with 9 rows each
###############################################################################

# ----- Group 1 runs -----
rows_comp_g1 <- list()
rows_out_g1  <- list()

k <- 1
for (vc in VC_levels_g1) {
  for (vt in vacc_types) {
    
    label <- sprintf("G1-%d: VC=%d%% | VE=40%% | ATB=0%% | vacc=%s",
                     k, round(100*vc), vt)
    
    res_i <- run_one_scenario(params, state_eq, cov = vc, VE = VE_g1, vacc_type = vt, atb_red_c = ATB_g1)
    
    rows_comp_g1[[k]] <- extract_compartments_tf(res_i$out, label)
    rows_out_g1[[k]]  <- extract_outputs_tf(res_i$metrics, label)
    
    k <- k + 1
  }
}

df_comp_g1 <- do.call(rbind, rows_comp_g1)
df_out_g1  <- do.call(rbind, rows_out_g1)

# ----- Group 2 runs -----
rows_comp_g2 <- list()
rows_out_g2  <- list()

k <- 1
for (ve in VE_levels_g2) {
  for (vt in vacc_types) {
    
    label <- sprintf("G2-%d: VC=60%% | VE=%d%% | ATB=0%% | vacc=%s",
                     k, round(100*ve), vt)
    
    res_i <- run_one_scenario(params, state_eq, cov = VC_g2, VE = ve, vacc_type = vt, atb_red_c = ATB_g2)
    
    rows_comp_g2[[k]] <- extract_compartments_tf(res_i$out, label)
    rows_out_g2[[k]]  <- extract_outputs_tf(res_i$metrics, label)
    
    k <- k + 1
  }
}

df_comp_g2 <- do.call(rbind, rows_comp_g2)
df_out_g2  <- do.call(rbind, rows_out_g2)

# ----- Group 3 runs -----
rows_comp_g3 <- list()
rows_out_g3  <- list()

k <- 1
for (atb in ATB_levels_g3) {
  for (vt in vacc_types) {
    
    label <- sprintf("G3-%d: VC=60%% | VE=40%% | ATB=%d%% (comm) | vacc=%s",
                     k, round(100*atb), vt)
    
    res_i <- run_one_scenario(params, state_eq, cov = VC_g3, VE = VE_g3, vacc_type = vt, atb_red_c = atb)
    
    rows_comp_g3[[k]] <- extract_compartments_tf(res_i$out, label)
    rows_out_g3[[k]]  <- extract_outputs_tf(res_i$metrics, label)
    
    k <- k + 1
  }
}

df_comp_g3 <- do.call(rbind, rows_comp_g3)
df_out_g3  <- do.call(rbind, rows_out_g3)

###############################################################################
# PLOT (3 pages): each page = 2 tables, 9 rows each
###############################################################################

plot_group_tables(df_comp_g1, df_out_g1, "GROUPE 1 (9 scénarios) — VE=40%, ATB=0%, VC = 30/60/90 × 3 types vaccins")
plot_group_tables(df_comp_g2, df_out_g2, "GROUPE 2 (9 scénarios) — VC=60%, ATB=0%, VE = 20/40/60 × 3 types vaccins")
plot_group_tables(df_comp_g3, df_out_g3, "GROUPE 3 (9 scénarios) — VC=60%, VE=40%, ATB(comm) = 5/10/20 × 3 types vaccins")


###############################################################################
# EXTRA PLOTS: same 3 pages, but OUTPUTS collapsed to totals (nv+v)
# Columns: Scenario, Prev_h, Prev_c, Inc_h, Inc_c
###############################################################################

collapse_outputs_total <- function(df_out, df_comp) {
  # df_out : has Prev_h_nv, Prev_h_v, Prev_c_nv, Prev_c_v, Inc_*_nv, Inc_*_v
  # df_comp: has N_h_nv, N_h_v, N_c_nv, N_c_v (same Scenario order)
  
  stopifnot(all(df_out$Scenario == df_comp$Scenario))
  
  Nh_nv <- as.numeric(df_comp$N_h_nv); Nh_v <- as.numeric(df_comp$N_h_v)
  Nc_nv <- as.numeric(df_comp$N_c_nv); Nc_v <- as.numeric(df_comp$N_c_v)
  
  # Prevalence totals (%) = weighted average
  Prev_h <- (as.numeric(df_out$Prev_h_nv) * Nh_nv + as.numeric(df_out$Prev_h_v) * Nh_v) / (Nh_nv + Nh_v)
  Prev_c <- (as.numeric(df_out$Prev_c_nv) * Nc_nv + as.numeric(df_out$Prev_c_v) * Nc_v) / (Nc_nv + Nc_v)
  
  # Incidence totals (per 100k PY) = weighted average
  Inc_h  <- (as.numeric(df_out$Inc_h_nv)  * Nh_nv + as.numeric(df_out$Inc_h_v)  * Nh_v) / (Nh_nv + Nh_v)
  Inc_c  <- (as.numeric(df_out$Inc_c_nv)  * Nc_nv + as.numeric(df_out$Inc_c_v)  * Nc_v) / (Nc_nv + Nc_v)
  
  data.frame(
    Scenario = df_out$Scenario,
    Prev_h = Prev_h,
    Prev_c = Prev_c,
    Inc_h  = Inc_h,
    Inc_c  = Inc_c
  )
}

df_out_g1_tot <- collapse_outputs_total(df_out_g1, df_comp_g1)
df_out_g2_tot <- collapse_outputs_total(df_out_g2, df_comp_g2)
df_out_g3_tot <- collapse_outputs_total(df_out_g3, df_comp_g3)

# --- Plot the same 3 pages, with the new "total" outputs table ---
plot_group_tables(df_comp_g1, df_out_g1_tot, "GROUPE 1 (9 scénarios) — OUTPUTS totaux (nv+v)")
plot_group_tables(df_comp_g2, df_out_g2_tot, "GROUPE 2 (9 scénarios) — OUTPUTS totaux (nv+v)")
plot_group_tables(df_comp_g3, df_out_g3_tot, "GROUPE 3 (9 scénarios) — OUTPUTS totaux (nv+v)")


