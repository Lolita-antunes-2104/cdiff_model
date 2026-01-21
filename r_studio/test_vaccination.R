# ---- TEST DE VALIDATION SIMPLE ----
source("0_packages.R")
source("1_model_new.R")         
source("2_metrics.R")         
source("6_scena_vacc_alphaS.R") 

# Vérifier si la correction est en place
cat("\n=== VÉRIFICATION DU CODE ===\n")
test_code <- deparse(cdiff_vacc)
if (any(grepl("ifelse\\(den_alpha_v == 0, 0, out_hc_v / den_alpha_v\\)", test_code, fixed = FALSE))) {
  cat("✓ Correction présente dans le code!\n")
} else {
  cat("✗ Correction ABSENTE - le fichier n'a pas été sauvegardé correctement\n")
  cat("Lignes contenant 'alpha_v':\n")
  print(grep("alpha_v", test_code, value = TRUE))
}

# ---- PARAMÈTRES DE BASE ----
params_base <- c(
  beta_h = 0.116404, beta_c = 0.010793, nu = 7.0,
  sigma_h = 0.000820, sigma_c = 0.000020,
  tau_h = 0.03554, tau_c = 0.00271,
  omega_h = 0.000, omega_c = 0.022, gamma = 0.012,
  epsilon = 0.055, p = 0.50, phi = 0.018,
  k_A = 4.0, k_II = 291.099, k_III = 668.060,
  w = 1, w_I = 1.5, w_II = 2, w_III = 2.5,
  delta = 0.13333,
  delta_I = 0.0625,
  delta_II = 0.05,
  delta_III = 0.04
)

# ---- ÉQUILIBRE ----
N_h <- 100
N_c <- 50000
init_cond <- create_initial_conditions(N_h, N_c)
time_to_eq <- seq(0, 4000, by = 1)

eq_result <- run_model_to_equilibrium(params_base, init_cond, time_to_eq)
calibrated_eq <- as.list(eq_result[nrow(eq_result), ])

# ---- PARAMÈTRES VACCINAUX ----
params_vacc <- c(
  beta_h = 0.116404, beta_c = 0.010793, nu = 7.0,
  sigma_h = 0.000820, sigma_c = 0.000020,
  tau_h = 0.03554, tau_c = 0.00271,
  omega_h = 0.000, omega_c = 0.022, gamma = 0.012,
  epsilon = 0.055, p = 0.50, phi = 0.018,
  k_A = 4.0, k_II = 291.099, k_III = 668.060,
  w = 1, w_I = 1.5, w_II = 2, w_III = 2.5,
  delta_nv = 0.13333, delta_I_nv = 0.0625,
  delta_II_nv = 0.05, delta_III_nv = 0.04,
  delta_v = 0.13333, delta_I_v = 0.0625,
  delta_II_v = 0.05, delta_III_v = 0.04,
  VE = 0
)

# ---- TEST 1 : 100% NON-VACCINÉS (VC = 0) ----
cat("\n=== TEST 1: 100% non-vaccinés ===\n")
init_test1 <- create_vacc_initial_conditions(calibrated_eq, VC = 0)
params_test1 <- params_vacc
params_test1["VE"] <- 0

cat("\n=== VÉRIFICATION DES CONDITIONS INITIALES ===\n")

# Afficher l'équilibre de base
cat("\nÉquilibre du modèle de base (premiers compartiments):\n")
print(unlist(calibrated_eq)[1:11])

# Vérifier les conditions initiales test1
cat("\n\nConditions initiales TEST 1 (VC=0):\n")
cat("Premiers compartiments _nv:\n")
print(init_test1[1:11])
cat("\nPremiers compartiments _v (devraient être 0):\n")
print(init_test1[23:33])

cat("\n\nSomme totale _nv:", sum(init_test1[grepl("_nv", names(init_test1))]), "\n")
cat("Somme totale _v:", sum(init_test1[grepl("_v", names(init_test1))]), "\n")

# Vérifier s'il y a des valeurs négatives ou NaN dans les conditions initiales
cat("\nVérification de validité:\n")
cat("  Valeurs négatives dans init_test1?", any(init_test1 < 0, na.rm = TRUE), "\n")
cat("  Valeurs NaN dans init_test1?", any(is.nan(init_test1)), "\n")
cat("  Valeurs infinies dans init_test1?", any(is.infinite(init_test1)), "\n")

# Tester la fonction cdiff_vacc directement à t=0
cat("\n=== TEST DIRECT DE cdiff_vacc à t=0 ===\n")
derivs_t0 <- cdiff_vacc(t = 0, pop = init_test1, params = params_test1)
cat("Nombre de dérivées:", length(derivs_t0[[1]]), "\n")
cat("Dérivées NaN?", any(is.nan(unlist(derivs_t0))), "\n")
cat("Dérivées infinies?", any(is.infinite(unlist(derivs_t0))), "\n")

# Afficher les premières dérivées
cat("\nPremières dérivées:\n")
print(unlist(derivs_t0)[1:11])

# Simuler juste les 5 premiers pas de temps pour voir où ça casse
cat("\n=== SIMULATION DES 5 PREMIERS PAS ===\n")
time_short <- seq(0, 5, by = 1)

ode_short <- as.data.frame(
  deSolve::lsoda(
    y = init_test1, 
    times = time_short, 
    func = cdiff_vacc, 
    parms = params_test1,
    maxsteps = 50000,
    rtol = 1e-6,
    atol = 1e-8
  )
)

print(ode_short[, 1:6])  # Afficher temps + premiers compartiments

# Vérifier à quel moment apparaissent les NaN
for (i in 1:nrow(ode_short)) {
  n_nan <- sum(is.nan(as.numeric(ode_short[i, -1])))
  if (n_nan > 0) {
    cat(sprintf("\n*** ALERTE: Temps t=%d: %d valeurs NaN détectées ***\n", 
                ode_short[i, "time"], n_nan))
    break
  }
}




# ---- DIAGNOSTIC DÉTAILLÉ DES DÉRIVÉES ----
cat("\n=== DIAGNOSTIC DÉTAILLÉ ===\n")

# Afficher TOUTES les dérivées pour voir lesquelles sont NaN
all_derivs <- unlist(derivs_t0)
nan_derivs <- names(all_derivs)[is.nan(all_derivs)]

cat("\nDérivées qui sont NaN:\n")
print(nan_derivs)

cat("\nValeurs des dérivées NaN:\n")
print(all_derivs[is.nan(all_derivs)])

# Tester manuellement le calcul des totaux
cat("\n=== TEST MANUEL DES CALCULS ===\n")

# Extraire les valeurs pour le calcul manuel
S0_h_nv <- init_test1["S0_h_nv"]
SA_h_nv <- init_test1["SA_h_nv"]
C0_h_nv <- init_test1["C0_h_nv"]
# etc.

# Test: population vaccinée est vide (VC=0)
cat("\nPopulations vaccinées (devraient être 0):\n")
vacc_pops <- init_test1[grepl("_v$", names(init_test1))]
cat("  Somme des compartiments _v:", sum(vacc_pops), "\n")
cat("  Min:", min(vacc_pops), "Max:", max(vacc_pops), "\n")

# Le problème probable: division par zéro dans le calcul de alpha_v
cat("\n=== CALCUL DE alpha_v (pour population vaccinée vide) ===\n")

# Simuler le calcul
w <- params_test1["w"]
w_I <- params_test1["w_I"]
w_II <- params_test1["w_II"]
w_III <- params_test1["w_III"]

# Pour la population vaccinée (qui est vide)
S0_c_v <- init_test1["S0_c_v"]
SA_c_v <- init_test1["SA_c_v"]
C0_c_v <- init_test1["C0_c_v"]
CA_c_v <- init_test1["CA_c_v"]
I_c_v <- init_test1["I_c_v"]
I_II_c_v <- init_test1["I_II_c_v"]
I_III_c_v <- init_test1["I_III_c_v"]

# Calcul du dénominateur
den_alpha_v_manual <- w * (S0_c_v + SA_c_v + C0_c_v + CA_c_v) + 
  w_I * I_c_v + 
  w_II * I_II_c_v + 
  w_III * I_III_c_v

cat("Dénominateur alpha_v:", den_alpha_v_manual, "\n")

# Calcul du numérateur
delta_v <- params_test1["delta_v"]
delta_I_v <- params_test1["delta_I_v"]

S0_h_v <- init_test1["S0_h_v"]
SA_h_v <- init_test1["SA_h_v"]
C0_h_v <- init_test1["C0_h_v"]
CA_h_v <- init_test1["CA_h_v"]
I_h_v <- init_test1["I_h_v"]

num_alpha_v_manual <- delta_v * (S0_h_v + SA_h_v + C0_h_v + CA_h_v) + 
  delta_I_v * I_h_v

cat("Numérateur alpha_v:", num_alpha_v_manual, "\n")

# Division
if (den_alpha_v_manual == 0) {
  cat("\n*** PROBLÈME: Division par zéro! ***\n")
  cat("Le dénominateur est 0 car la population vaccinée en communauté est vide\n")
  cat("Cela crée alpha_v = 0/0 = NaN\n")
} else {
  alpha_v_manual <- num_alpha_v_manual / den_alpha_v_manual
  cat("alpha_v:", alpha_v_manual, "\n")
}

