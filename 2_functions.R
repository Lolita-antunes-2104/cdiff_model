###############################################################################
############################# 2 : FUNCTIONS ###################################
###############################################################################

###############################################################################
#  1. MODEL FUNCTIONS FOR CALIBRATION MODEL
###############################################################################

# Create initial conditions for calibration model
create_initial_conditions_precalibration <- function(N_h, N_c, prev_primo = 0.01, prev_rec = 0.005) {
  # Safety check to not exceed 100%
  if (4*prev_primo + 6*prev_rec >= 1) {
    stop("Initial prevalences are too high and lead to an invalid initial state (total population exceeds 100%).")}
  # Hospital
  SA_h <- prev_primo * N_h
  C0_h <- prev_primo * N_h
  CA_h <- prev_primo * N_h
  I_h <- prev_primo * N_h
  S_II_h  <- prev_rec * N_h
  C_II_h <- prev_rec * N_h
  I_II_h <- prev_rec * N_h
  S_III_h <- prev_rec * N_h
  C_III_h <- prev_rec * N_h
  I_III_h <- prev_rec * N_h
  S0_h <- N_h - (C0_h + CA_h + C_II_h + C_III_h + I_h + I_II_h + I_III_h + S_II_h + S_III_h + SA_h)
  # Community
  SA_c <- prev_primo * N_c
  C0_c <- prev_primo * N_c
  CA_c <- prev_primo * N_c
  I_c <- prev_primo * N_c
  S_II_c  <- prev_rec * N_c
  C_II_c <- prev_rec * N_c
  I_II_c <- prev_rec * N_c
  S_III_c <- prev_rec * N_c
  C_III_c <- prev_rec * N_c
  I_III_c <- prev_rec * N_c
  S0_c <- N_c - (C0_c + CA_c + C_II_c + C_III_c + I_c + I_II_c + I_III_c + S_II_c + S_III_c + SA_c)
  # Cumulative incidence (reset)
  CumI_h       <- 0
  primo_CumI_h <- 0 
  rec_CumI_h   <- 0
  CumI_c       <- 0
  primo_CumI_c <- 0 
  rec_CumI_c   <- 0
  
  # Return state vector
  return(c(
    # Hospital
    S0_h = S0_h, SA_h = SA_h, C0_h = C0_h, CA_h = CA_h, I_h = I_h,
    S_II_h = S_II_h, C_II_h = C_II_h, I_II_h = I_II_h,
    S_III_h = S_III_h, C_III_h = C_III_h, I_III_h = I_III_h,
    CumI_primo_h = CumI_primo_h, CumI_rec_h = CumI_rec_h, CumI_h = CumI_h, 
    # Community
    S0_c = S0_c, SA_c = SA_c, C0_c = C0_c, CA_c = CA_c, I_c = I_c,
    S_II_c = S_II_c, C_II_c = C_II_c, I_II_c = I_II_c,
    S_III_c = S_III_c, C_III_c = C_III_c, I_III_c = I_III_c,
    CumI_primo_c = CumI_primo_c, CumI_rec_c = CumI_rec_c, CumI_c = CumI_c))
}

# Check if model has reached equilibrium by evalutating the derivative





###############################################################################
# 2. METRICS FUNCTIONS FOR CALIBRATION MODEL
###############################################################################

# Simple R0 proxy (time series)
compute_R0 <- function(params) {
  with(as.list(params), {
    
    # Au DFE, la population est entierement susceptible.
    # lambda_h = beta_h * (C_h + nu_h * I_tot_h) / N_h
    # On linearise : lambda_h ~= beta_h / N_h * (C_h + nu_h * I_tot_h)
    # Le coefficient devant chaque C ou I dans lambda est donc beta/N.
    # Au DFE, N_h et N_c sont les populations totales initiales.
    
    # On suppose qu'on a N_h_0 et N_c_0 comme populations au DFE
    # (ou on les passe en parametre)
    
    bh <- beta_h / N_h_0   # coefficient de transmission par individu, hôpital
    bc <- beta_c / N_c_0   # même, communauté
    
    # Ordre des compartiments infectieux :
    # [C0_h, CA_h, I_h, C_II_h, I_II_h, C_III_h, I_III_h,
    #  C0_c, CA_c, I_c, C_II_c, I_II_c, C_III_c, I_III_c]
    #   1      2    3     4       5       6        7
    #   8      9   10    11      12      13       14
    
    # =========================================================
    # MATRICE F : nouveaux infectés seulement
    # =========================================================
    # Les nouveaux infectés viennent de lambda * S.
    # Au DFE, S0 = N (tout le monde est susceptible S0_h = N_h_0, etc.)
    # lambda_h = bh * (C0_h + CA_h + C_II_h + C_III_h + nu_h*(I_h + I_II_h + I_III_h))
    #
    # Nouveaux infectés au DFE :
    #   dC0_h : +lambda_h * S0_h  -> S0_h = N_h_0 -> coeff = bh * N_h_0 = beta_h
    #   dCA_h : +lambda_h * SA_h  -> SA_h = 0 au DFE -> 0
    #   dC_II_h: +lambda_h * S_II_h -> S_II_h = 0 au DFE -> 0
    #   dC_III_h: même -> 0
    # Même raisonnement pour communauté.
    # Donc seuls C0_h et C0_c reçoivent des nouveaux infectés au DFE.
    
    F <- matrix(0, nrow=14, ncol=14)
    
    # Ligne 1 (C0_h) : nouveaux infectés = beta_h * (C + nu*I) au DFE
    # Colonnes correspondant à C0_h(1), CA_h(2), I_h(3), C_II_h(4), I_II_h(5), C_III_h(6), I_III_h(7)
    F[1, 1] <- beta_h          # C0_h
    F[1, 2] <- beta_h          # CA_h
    F[1, 3] <- beta_h * nu_h   # I_h
    F[1, 4] <- beta_h          # C_II_h
    F[1, 5] <- beta_h * nu_h   # I_II_h
    F[1, 6] <- beta_h          # C_III_h
    F[1, 7] <- beta_h * nu_h   # I_III_h
    
    # Ligne 8 (C0_c) : même pour communauté
    F[8,  8] <- beta_c
    F[8,  9] <- beta_c
    F[8, 10] <- beta_c * nu_c
    F[8, 11] <- beta_c
    F[8, 12] <- beta_c * nu_c
    F[8, 13] <- beta_c
    F[8, 14] <- beta_c * nu_c
    
    # =========================================================
    # MATRICE V : transitions de sortie des compartiments infectieux
    # (termes diagonaux = sortie, termes hors-diagonaux = flux vers d'autres compartiments infectieux)
    # V est la matrice de transition ENTRE compartiments infectieux.
    # Convention : V[i,j] = flux de j vers i (comme pour F).
    # La diagonale V[i,i] = somme des taux de sortie de i.
    # =========================================================
    
    # --- Hôpital ---
    # C0_h (1) : sorties = gamma_h + tau_h + sigma_h + delta
    V[1,1] <- gamma_h + tau_h + sigma_h + delta
    
    # CA_h (2) : sorties = gamma_h + omega_h + k_A*sigma_h + delta
    V[2,2] <- gamma_h + omega_h + k_A*sigma_h + delta
    
    # I_h (3) : sorties = epsilon_h + delta_I
    V[3,3] <- epsilon_h + delta_I
    # I_h reçoit de C0_h et CA_h
    V[3,1] <- -sigma_h          # flux C0_h -> I_h
    V[3,2] <- -k_A*sigma_h      # flux CA_h -> I_h
    
    # C_II_h (4) : sorties = gamma_h + k_II*sigma_h + delta
    V[4,4] <- gamma_h + k_II*sigma_h + delta
    # C_II_h reçoit (1-p_h)*epsilon_h*I_h
    V[4,3] <- -(1-p_h)*epsilon_h
    
    # I_II_h (5) : sorties = epsilon_h + delta_II
    V[5,5] <- epsilon_h + delta_II
    # I_II_h reçoit de C_II_h
    V[5,4] <- -k_II*sigma_h
    
    # C_III_h (6) : sorties = gamma_h + k_III*sigma_h + delta
    V[6,6] <- gamma_h + k_III*sigma_h + delta
    # C_III_h reçoit (1-p_h)*epsilon_h*(I_II_h + I_III_h)
    V[6,5] <- -(1-p_h)*epsilon_h   # depuis I_II_h
    V[6,7] <- -(1-p_h)*epsilon_h   # depuis I_III_h
    
    # I_III_h (7) : sorties = epsilon_h + delta_III
    V[7,7] <- epsilon_h + delta_III
    # I_III_h reçoit de C_III_h
    V[7,6] <- -k_III*sigma_h
    
    # --- Communauté (même structure, indices 8-14) ---
    # C0_c (8)
    V[8,8] <- gamma_c + tau_c + sigma_c + alpha
    
    # CA_c (9)
    V[9,9] <- gamma_c + omega_c + k_A*sigma_c + alpha
    
    # I_c (10)
    V[10,10] <- epsilon_c + alpha_I
    V[10,8]  <- -sigma_c
    V[10,9]  <- -k_A*sigma_c
    
    # C_II_c (11)
    V[11,11] <- gamma_c + k_II*sigma_c + alpha
    V[11,10] <- -(1-p_c)*epsilon_c
    
    # I_II_c (12)
    V[12,12] <- epsilon_c + alpha_II
    V[12,11] <- -k_II*sigma_c
    
    # C_III_c (13)
    V[13,13] <- gamma_c + k_III*sigma_c + alpha
    V[13,12] <- -(1-p_c)*epsilon_c   # depuis I_II_c
    V[13,14] <- -(1-p_c)*epsilon_c   # depuis I_III_c
    
    # I_III_c (14)
    V[14,14] <- epsilon_c + alpha_III
    V[14,13] <- -k_III*sigma_c
    
    # =========================================================
    # R0 = rayon spectral de F %*% solve(V)
    # =========================================================
    NGM <- F %*% solve(V)
    R0 <- max(Mod(eigen(NGM)$values))
    
    return(R0)
  })
}

# Main metrics function for the calibration model
compute_metrics_calib <- function(out, params) {
  
  out <- as.data.frame(out)
  
  needed <- c(
    "time",
    # Hospital compartments
    "S0_h","SA_h","C0_h","CA_h","I_h","S_II_h","C_II_h","I_II_h","S_III_h","C_III_h","I_III_h",
    # Community compartments
    "S0_c","SA_c","C0_c","CA_c","I_c","S_II_c","C_II_c","I_II_c","S_III_c","C_III_c","I_III_c",
    # Cumulative counters
    "CumI_h","primo_CumI_h","rec_CumI_h",
    "CumI_c","primo_CumI_c","rec_CumI_c"
  )

  with(as.list(params), {
    
    # ---- Populations (Nh, Nc, Ntot) ----
    final_tot_h <- compute_totals(out$S0_h, out$C0_h, out$SA_h, out$CA_h, out$I_h, out$S_II_h, out$C_II_h, out$I_II_h, out$S_III_h, out$C_III_h, out$I_III_h)
    final_tot_c <- compute_totals(out$S0_c, out$C0_c, out$SA_c, out$CA_c, out$I_c, out$S_II_c, out$C_II_c, out$I_II_c, out$S_III_c, out$C_III_c, out$I_III_c)
    N_h <- final_tot_h$N
    N_c <- final_tot_c$N
    N_tot <- N_h + N_c
    
    # ---- Carriage prevalence ----
    C_h <- final_tot_h$C
    C_c <- final_tot_c$C
    prev_h <- C_h / N_h
    prev_c <- C_c / N_c
    
    # ---- Instantaneous incidence (total / primary / recurrent) ----
    inc_h_total <- sigma_h*out$C0_h + k_A*sigma_h*out$CA_h + k_II*sigma_h*out$C_II_h + k_III*sigma_h*out$C_III_h
    inc_h_primo <- sigma_h*out$C0_h + k_A*sigma_h*out$CA_h
    inc_h_rec   <- k_II*sigma_h*out$C_II_h + k_III*sigma_h*out$C_III_h
    
    inc_c_total <- sigma_c*out$C0_c + k_A*sigma_c*out$CA_c + k_II*sigma_c*out$C_II_c + k_III*sigma_c*out$C_III_c
    inc_c_primo <- sigma_c*out$C0_c + k_A*sigma_c*out$CA_c
    inc_c_rec   <- k_II*sigma_c*out$C_II_c + k_III*sigma_c*out$C_III_c
    
    # Normalized instantaneous incidence
    # - absolute: per N_tot
    # - relative: per setting population (N_h or N_c)
    inc_h_total_abs <- inc_h_total / N_tot
    inc_h_primo_abs <- inc_h_primo / N_tot
    inc_h_rec_abs   <- inc_h_rec / N_tot
    
    inc_c_total_abs <- inc_c_total / N_tot
    inc_c_primo_abs <- inc_c_primo / N_tot
    inc_c_rec_abs   <- inc_c_rec / N_tot
    
    inc_h_total_abs <- inc_h_total / N_h
    inc_h_primo_abs <- inc_h_primo / N_h
    inc_h_rec_abs   <- inc_h_rec / N_h
    
    inc_h_total_abs <- inc_h_total / N_c
    inc_h_primo_abs <- inc_h_primo / N_c
    inc_h_rec_abs   <- inc_h_rec / N_c
    
    # ---- Cumulative incidence from CumI_* compartments ----
    # Make cumulative start at 0 (difference from t=0)
    cum_h_total <- out$CumI_h - out$CumI_h[1]
    cum_h_primo <- out$primo_CumI_h - out$primo_CumI_h[1]
    cum_h_rec   <- out$rec_CumI_h - out$rec_CumI_h[1]
    
    cum_c_total <- out$CumI_c - out$CumI_c[1]
    cum_c_primo <- out$primo_CumI_c - out$primo_CumI_c[1]
    cum_c_rec   <- out$rec_CumI_c - out$rec_CumI_c[1]
    
    # Normalized cumulative incidence
    cum_h_total_abs <- cum_h_total / N_tot
    cum_h_primo_abs <- cum_h_primo / N_tot
    cum_h_rec_abs   <- cum_h_rec / N_tot
    
    cum_c_total_abs <- cum_c_total / N_tot
    cum_c_primo_abs <- cum_c_primo / N_tot
    cum_c_rec_abs   <- cum_c_rec / N_tot
    
    cum_h_total_rel <- cum_h_total / N_h
    cum_h_primo_rel <- cum_h_primo / N_h
    cum_h_rec_rel   <- cum_h_rec / N_h
    
    cum_c_total_rel <- cum_c_total / N_c
    cum_c_primo_rel <- cum_c_primo / N_c
    cum_c_rec_rel   <- cum_c_rec / N_c
    
    # ---- Recurrence prevalence ----
    rec1_h <- out$I_II_h / out$I_h
    rec2_h <- out$I_III_h / out$I_II_h
    
    rec1_c <- out$I_II_c / out$I_c
    rec2_c <- out$I_III_c / out$I_II_c
    
    # ---- R0 (simple proxy time series) ----
    R0_df <- compute_R0(out, params)
    
    return(list(
      population = data.frame(time = out$time, N_h = N_h, N_c = N_c, N_tot = N_tot),
      
      carriage = data.frame(time = out$time, prev_h = prev_h, prev_c = prev_c),
      
      incidence_instant = data.frame(time = out$time,
        inc_h_total, inc_h_primo, inc_h_rec,
        inc_c_total, inc_c_primo, inc_c_rec,
        inc_h_total_abs, inc_h_primo_abs, inc_h_rec_abs,
        inc_c_total_abs, inc_c_primo_abs, inc_c_rec_abs,
        inc_h_total_rel, inc_h_primo_rel, inc_h_rec_rel,
        inc_c_total_rel, inc_c_primo_rel, inc_c_rec_rel
      ),
      
      incidence_cumulative = data.frame(time = out$time,
        cum_h_total, cum_h_primo, cum_h_rec,
        cum_c_total, cum_c_primo, cum_c_rec,
        cum_h_total_abs, cum_h_primo_abs, cum_h_rec_abs,
        cum_c_total_abs, cum_c_primo_abs, cum_c_rec_abs,
        cum_h_total_rel, cum_h_primo_rel, cum_h_rec_rel,
        cum_c_total_rel, cum_c_primo_rel, cum_c_rec_rel
      ),
      
      recurrence = data.frame(time = out$time, rec1_h = rec1_h, rec2_h = rec2_h, rec1_c = rec1_c, rec2_c = rec2_c),
      
      R0 = R0_df
    ))
  })
}

###############################################################################
# 3. METRICS FUNCTIONS FOR VACCINATION MODEL
###############################################################################

# 1) prevalence
# prendre à la fin de résoudre cdiff_model_for_calibration par desolve lsoda 
# prendre le C de compute_metrics 
# et faire C_h /N_h et C_c / N_c

# 2) Compute instantaneous CDI incidence (primary, recurrent, or total)
# prendre à la fin de résoudre cdiff_model_for_calibration par desolve lsoda 
# inc_h <- sigma_h*C0_h + k_A*sigma_h*CA_h + k_II*sigma_h*C_II_h + k_III*sigma_h*C_III_h
# primo_inc_h <- sigma_h*C0_h + k_A*sigma_h*CA_h 
# rec_inc_h <- k_II*sigma_h*C_II_h + k_III*sigma_h*C_III_h
# inc_c <- sigma_c*C0_c + k_A*sigma_c*CA_c + k_II*sigma_c*C_II_c + k_III*sigma_c*C_III_c
# primo_inc_c <- sigma_c*C0_c + k_A*sigma_c*CA_c 
# rec_inc_c <- k_II*sigma_c*C_II_c + k_III*sigma_c*C_III_c


# Normalised instantaneous CDI incidence (absolue (/N_tot) / relative (_N_h or N_c))

# 3) Compute cumulative CDI incidence from CumI_* compartments (primary, recurrent or total)

# Normalised(absolue Compute cumulative CDI incidence  (/N_tot) / relative (_N_h ou N_c))

# 4) Compute recurrence prevalence from infected compartments
# prendre à la fin de résoudre cdiff_model_for_calibration par desolve lsoda 
# rec_1 = I_II / I , rec_2 = I_III / I_II

# 5) Faire fonction simple pour R0

# FOR SCENARIO MODEL



