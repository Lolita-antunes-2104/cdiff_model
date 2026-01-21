###############################################################################
######################## 2 : METRICS FUNCTIONS ################################
###############################################################################

# Helper : récupère un paramètre dans un vecteur nommé
get_param <- function(params_vec, name) { # Déclare la fonction (vecteur de params + nom à extraire)
  if (is.null(names(params_vec)) || !(name %in% names(params_vec))) {  # Vérifie que params_vec est nommé et contient `name`
    stop(sprintf("Paramètre manquant ou params_vec non nommé : '%s'", name)) # Stoppe avec un message explicite si absent
  }
  as.numeric(params_vec[[name]]) # Renvoie la valeur du paramètre convertie en numérique
}

# Helper : division sûre (évite NaN/Inf)
safe_div <- function(num, den, default = 0) { # Déclare la fonction (num/den + valeur par défaut)
  ifelse(den == 0, default, num / den) # Si den==0 → default, sinon → num/den (vectorisé)
}

# Calcule les totaux de population N à partir d'un état
compute_population_totals <- function(last_state, setting = c("h", "c", "both")) {
  setting <- match.arg(setting)
  
  N_h <- last_state$S0_h + last_state$SA_h + last_state$S_II_h + last_state$S_III_h + 
    last_state$C0_h + last_state$CA_h + last_state$C_II_h + last_state$C_III_h + 
    last_state$I_h + last_state$I_II_h + last_state$I_III_h
  
  N_c <- last_state$S0_c + last_state$SA_c + last_state$S_II_c + last_state$S_III_c + 
    last_state$C0_c + last_state$CA_c + last_state$C_II_c + last_state$C_III_c + 
    last_state$I_c + last_state$I_II_c + last_state$I_III_c
  
  if (setting == "h") return(N_h)
  if (setting == "c") return(N_c)
  if (setting == "both") return(N_h + N_c)}

# Calcule prévalence de portage asymptomatique
compute_carriage_prevalence <- function(last_state, setting = c("h","c","both")) {
  setting <- match.arg(setting) # Valide/normalise `setting` pour qu’il soit uniquement "h", "c" ou "both" (sinon erreur)
  
  C_h <- last_state$C0_h + last_state$CA_h + last_state$C_II_h + last_state$C_III_h
  C_c <- last_state$C0_c + last_state$CA_c + last_state$C_II_c + last_state$C_III_c
  
  N_h <- compute_population_totals(last_state, "h")
  N_c <- compute_population_totals(last_state, "c")
  
  if (setting == "h") return(C_h / N_h)
  if (setting == "c") return(C_c / N_c)
  if (setting == "both") return((C_h + C_c) / (N_h + N_c))}

# Incidence CDI (total / primo / récidive) par setting
compute_CDI_incidence <- function(last_state, params_vec, setting = c("h", "c", "both"), type = c("total", "primo", "recurrent")) {
  setting <- match.arg(setting) # Force setting à être "h", "c" ou "both" (sinon erreur)
  type <- match.arg(type) # Force type à être "total", "primo" ou "recurrent" (sinon erreur)
  
  # Ici on récupère σ et les k car l’incidence est un flux dépendant des taux de transition, contrairement à la prévalence (stock C/N) qui n’en a pas besoin.
  k_A   <- get_param(params_vec, "k_A")
  k_II  <- get_param(params_vec, "k_II")
  k_III <- get_param(params_vec, "k_III")
  sigma_h <- get_param(params_vec, "sigma_h")
  sigma_c <- get_param(params_vec, "sigma_c")
  
  calc_inc <- function(sigma, C0, CA, C_II, C_III) {
    primo <- as.numeric(sigma * C0 + k_A * sigma * CA)
    rec   <- as.numeric(k_II * sigma * C_II + k_III * sigma * C_III)
    c(primo = primo, recurrent = rec, total = primo + rec)}
  
  inc_h <- calc_inc(sigma_h, last_state$C0_h, last_state$CA_h, last_state$C_II_h, last_state$C_III_h)
  inc_c <- calc_inc(sigma_c, last_state$C0_c, last_state$CA_c, last_state$C_II_c, last_state$C_III_c)
  
  if (setting == "h") return(as.numeric(inc_h[type]))
  if (setting == "c") return(as.numeric(inc_c[type]))
  if (setting == "both") return(as.numeric((inc_h + inc_c)[type]))}

# Cumul d’incidence CDI sur une période (total / primo / récidive) et par setting (par défaut = 365 jours)
compute_cumulative_CDI_incidence <- function(ode_result, params_vec, period_days = 365, setting = c("h","c","both"), type = c("total","primo","recurrent")) {  # Déclare la fonction + options (période, setting, type)
  setting <- match.arg(setting)  # Valide/normalise setting ("h","c","both")
  type <- match.arg(type) # Valide/normalise type ("total","primo","recurrent")
  
  n_rows <- nrow(ode_result) # Nombre total de points (lignes) dans la simulation
  start_idx <- max(1, n_rows - period_days + 1) # Index de début : les `period_days` derniers jours (sans dépasser 1)
  period_data <- ode_result[start_idx:n_rows, ] # Sous-ensemble : uniquement la période sélectionnée
  
  daily_incidence <- sapply(1:nrow(period_data), function(i) { # Calcule l’incidence pour chaque jour de la période…
    compute_CDI_incidence(period_data[i, ], params_vec, # …en appelant la fonction d’incidence sur la ligne i
                          setting = setting, type = type)
  })
  
  sum(daily_incidence) # Somme sur la période = incidence cumulée
}

# Prévalence de récidive (1ère ou 2ème+) par setting
compute_recurrence_prevalence <- function(last_state, setting = c("h", "c", "both"), type = c("rec_1", "rec_2"), default = 0) {
  setting <- match.arg(setting)
  type <- match.arg(type)
  
  calc_prev_rec <- function(I, I_II, I_III) { 
    c(rec_1 = safe_div(I_II, I, default = default), # safe_div évite NaN/Inf quand le dénominateur (I) vaut 0 
      rec_2 = safe_div(I_III, I_II, default = default))
    }
  
  prev_h <- calc_prev_rec(last_state$I_h, last_state$I_II_h, last_state$I_III_h)
  prev_c <- calc_prev_rec(last_state$I_c, last_state$I_II_c, last_state$I_III_c)
  prev_both <- calc_prev_rec(last_state$I_h + last_state$I_c, last_state$I_II_h + last_state$I_II_c, last_state$I_III_h + last_state$I_III_c)
  
  if (setting == "h") return(prev_h[type])
  if (setting == "c") return(prev_c[type])
  if (setting == "both") return(prev_both[type])
}

# Calcule la proportion de personnes en dysbiose post-ATB (SA + CA) / N
compute_dysbiosis_count <- function(last_state, setting = c("h", "c", "both"), default = 0) {
  setting <- match.arg(setting)
  
  dysb_h <- last_state$SA_h + last_state$CA_h
  dysb_c <- last_state$SA_c + last_state$CA_c
  
  N_h <- compute_population_totals(last_state, "h")
  N_c <- compute_population_totals(last_state, "c")
  
  if (setting == "h") return(dysb_h / N_h)
  if (setting == "c") return(dysb_c / N_c)
  if (setting == "both") return((dysb_h + dysb_c) / (N_h + N_c))}

# # Calcule R0 = formule quentin
# compute_R0_approx <- function(params_vec, setting = c("h", "c"), default = 0) {
#   setting <- match.arg(setting)
#   
#   gamma <- get_param(params_vec, "gamma")
#   epsilon <- get_param(params_vec, "epsilon")
#   nu <- get_param(params_vec, "nu")
#   
#   if (setting == "h") {
#     beta_h  <- get_param(params_vec, "beta_h")
#     sigma_h <- get_param(params_vec, "sigma_h")
#     delta   <- get_param(params_vec, "delta")
#     delta_I <- get_param(params_vec, "delta_I")
#     
#     den_C <- gamma + sigma_h + delta
#     den_I <- epsilon + delta_I
#     if (den_C <= 0 || den_I <= 0) return(default)
#     
#     return(beta_h / den_C + (nu * beta_h / den_I) * (sigma_h / den_C))
#   }
#   
#   if (setting == "c") {
#     beta_c  <- get_param(params_vec, "beta_c")
#     sigma_c <- get_param(params_vec, "sigma_c")
#     # alpha   <- get_param(params_vec, "alpha")
#     # alpha_I <- get_param(params_vec, "alpha_I")
#     alpha   <- get_param(params_vec, "delta") * get_param(params_vec, "w")
#     alpha_I <- get_param(params_vec, "delta_I") * get_param(params_vec, "w_I")
#     
#     den_C <- gamma + sigma_c + alpha
#     den_I <- epsilon + alpha_I
#     if (den_C <= 0 || den_I <= 0) return(default)
#     
#     return(beta_c / den_C + (nu * beta_c / den_I) * (sigma_c / den_C))}
# }

# Calcule R0 approximatif
compute_alpha_eq <- function(last_state, params_vec) {
  
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
  
  out_hc <- params_vec["delta"]   * (tot_h$S + tot_h$C) +
    params_vec["delta_I"] * last_state$I_h +
    params_vec["delta_II"]* last_state$I_II_h +
    params_vec["delta_III"]*last_state$I_III_h
  
  den_alpha <- params_vec["w"]   * (tot_c$S + tot_c$C) +
    params_vec["w_I"] * last_state$I_c +
    params_vec["w_II"]* last_state$I_II_c +
    params_vec["w_III"]*last_state$I_III_c
  
  alpha <- out_hc / den_alpha
  
  list(
    alpha   = alpha,
    alpha_I = alpha * params_vec["w_I"]
  )
}

compute_R0_approx <- function(last_state, params_vec,
                              setting = c("h", "c"), default = 0) {
  
  setting <- match.arg(setting)
  
  gamma   <- get_param(params_vec, "gamma")
  epsilon <- get_param(params_vec, "epsilon")
  nu      <- get_param(params_vec, "nu")
  k_A     <- get_param(params_vec, "k_A")
  k_II    <- get_param(params_vec, "k_II")
  k_III   <- get_param(params_vec, "k_III")
  
  calc_R0 <- function(suff, sigma, tau, omega, flow_out, flow_out_I) {
    
    C0   <- last_state[[paste0("C0_", suff)]]
    CA   <- last_state[[paste0("CA_", suff)]]
    CII  <- last_state[[paste0("C_II_", suff)]]
    CIII <- last_state[[paste0("C_III_", suff)]]
    
    C_tot <- C0 + CA + CII + CIII
    if (!is.finite(C_tot) || C_tot <= 0) return(default)
    
    T_C0   <- 1 / (gamma + tau + sigma + omega + flow_out)
    T_CA   <- 1 / (gamma + k_A * sigma + omega + flow_out)
    T_CII  <- 1 / (gamma + k_II * sigma + flow_out)
    T_CIII <- 1 / (gamma + k_III * sigma + flow_out)
    
    w_C0   <- C0   / C_tot
    w_CA   <- CA   / C_tot
    w_CII  <- CII  / C_tot
    w_CIII <- CIII / C_tot
    
    T_C <- w_C0*T_C0 + w_CA*T_CA + w_CII*T_CII + w_CIII*T_CIII
    T_I <- 1 / (epsilon + flow_out_I)
    
    beta <- get_param(params_vec, paste0("beta_", suff))
    beta * (T_C + nu * T_I)
  }
  
  if (setting == "h") {
    return(calc_R0(
      "h",
      sigma = get_param(params_vec, "sigma_h"),
      tau   = get_param(params_vec, "tau_h"),
      omega = get_param(params_vec, "omega_h"),
      flow_out   = get_param(params_vec, "delta"),
      flow_out_I = get_param(params_vec, "delta_I")
    ))
  }
  
  if (setting == "c") {
    alpha_eq <- compute_alpha_eq(last_state, params_vec)
    return(calc_R0(
      "c",
      sigma = get_param(params_vec, "sigma_c"),
      tau   = get_param(params_vec, "tau_c"),
      omega = get_param(params_vec, "omega_c"),
      flow_out   = alpha_eq$alpha,
      flow_out_I = alpha_eq$alpha_I
    ))
  }
}


###############################################################################
# ---- complete metrics function ----
compute_all_metrics <- function(last_state, params_vec, ode_result = NULL, targets = NULL) {
  
  # Conversion automatique si nécessaire
  if (!is.list(last_state)) {
    last_state <- as.list(last_state[1, ])
  }
  
  N_h <- compute_population_totals(last_state, "h")
  N_c <- compute_population_totals(last_state, "c")
  
  # Carriage
  portage_tot <- compute_carriage_prevalence(last_state, "both")
  portage_h <- compute_carriage_prevalence(last_state, "h")
  portage_c <- compute_carriage_prevalence(last_state, "c")
  
  # Incidence CDI instantanée (total = primo + recidive)
  inc_h <- unname(compute_CDI_incidence(last_state, params_vec, "h", "total")) / (N_h+N_c) 
  inc_c <- unname(compute_CDI_incidence(last_state, params_vec, "c", "total")) / (N_h+N_c) 
  inc_tot <- inc_h + inc_c
  
  # Incidence primo-CDI instantanée
  inc_primo_h <- unname(compute_CDI_incidence(last_state, params_vec, "h", "primo")) / (N_h+N_c) 
  inc_primo_c <- unname(compute_CDI_incidence(last_state, params_vec, "c", "primo")) / (N_h+N_c) 
  inc_primo_tot <- inc_primo_h + inc_primo_c
  
  # Incidence rCDI instantanée
  inc_rec_h <- unname(compute_CDI_incidence(last_state, params_vec, "h", "recurrent")) / (N_h+N_c)
  inc_rec_c <- unname(compute_CDI_incidence(last_state, params_vec, "c", "recurrent")) / (N_h+N_c) 
  inc_rec_tot   <- inc_rec_h   + inc_rec_c
  
  # Recurrence prevalence
  recid_1_tot <- compute_recurrence_prevalence(last_state, "both", "rec_1")
  recid_1_h <- compute_recurrence_prevalence(last_state, "h", "rec_1")
  recid_1_c <- compute_recurrence_prevalence(last_state, "c", "rec_1")
  
  recid_2_tot <- compute_recurrence_prevalence(last_state, "both", "rec_2")
  recid_2_h <- compute_recurrence_prevalence(last_state, "h", "rec_2")
  recid_2_c <- compute_recurrence_prevalence(last_state, "c", "rec_2")
  
  # R0
  R0_h <- compute_R0_approx(last_state, params_vec, "h")
  R0_c <- compute_R0_approx(last_state, params_vec, "c")
  
  metrics <- list(
    portage_tot = portage_tot,
    portage_h = portage_h,
    portage_c = portage_c,
    inc_tot = inc_tot,
    inc_h = inc_h,
    inc_c = inc_c,
    inc_primo_tot = inc_primo_tot,
    inc_primo_h = inc_primo_h,
    inc_primo_c = inc_primo_c,
    inc_rec_tot = inc_rec_tot,
    inc_rec_h = inc_rec_h,
    inc_rec_c = inc_rec_c,
    recid_1_tot = recid_1_tot,
    recid_1_h = recid_1_h,
    recid_1_c = recid_1_c,
    recid_2_tot = recid_2_tot,
    recid_2_h = recid_2_h,
    recid_2_c = recid_2_c,
    R0_h = R0_h,
    R0_c = R0_c,
    N_h = N_h,
    N_c = N_c
  )
  
  # Incidences cumulées (si ode_result fourni)
  if (!is.null(ode_result)) {
    metrics$inc_cum_tot = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "both", "total")
    metrics$inc_cum_h = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "h", "total")
    metrics$inc_cum_c = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "c", "total")
    
    metrics$inc_cum_primo_tot = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "both", "primo")
    metrics$inc_cum_primo_h = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "h", "primo")
    metrics$inc_cum_primo_c = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "c", "primo")
    
    metrics$inc_cum_rec_tot = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "both", "recurrent")
    metrics$inc_cum_rec_h = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "h", "recurrent")
    metrics$inc_cum_rec_c = compute_cumulative_CDI_incidence(ode_result, params_vec, 365, "c", "recurrent")
  }
  
  # Add errors if targets provided
  if (!is.null(targets)) {
    metrics$errors <- list(
      err_portage_h = abs(metrics$portage_h - targets$portage_h) / targets$portage_h,
      err_portage_c = abs(metrics$portage_c - targets$portage_c) / targets$portage_c,
      err_inc_h = abs(metrics$inc_h - targets$incidence_h) / targets$incidence_h,
      err_inc_c = abs(metrics$inc_c - targets$incidence_c) / targets$incidence_c,
      err_recid_1_tot = abs(metrics$recid_1_tot - targets$recid_1) / targets$recid_1,
      err_recid_2_tot = abs(metrics$recid_2_tot - targets$recid_2) / targets$recid_2
    )
  }
  
  return(metrics)
}

