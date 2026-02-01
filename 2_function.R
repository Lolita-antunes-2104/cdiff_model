###############################################################################
############################# 2 : FUNCTIONS ###################################
###############################################################################

###############################################################################
#  1. MODEL FUNCTIONS
###############################################################################

# Compute total population counts for a given setting
compute_totals <- function(S0, C0, SA, CA, I, S_II, C_II, I_II, S_III, C_III, I_III) {
  return(list(S = S0 + SA + S_II + S_III,
              C = C0 + CA + C_II + C_III,
              I_tot = I + I_II + I_III,
              N = S0 + SA + S_II + S_III + C0 + CA + C_II + C_III + I + I_II + I_III))
}

# Compute force of infection (transmission rate)
compute_lambda <- function(beta, C, I_tot, N, nu) {
  if (N == 0) return(0)
  return(beta * (C + nu * I_tot) / N)
}

# Create initial conditions for CALIBRATION model
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
  
  # Return state vector
  return(c(
    # Hospital
    S0_h = S0_h, SA_h = SA_h, C0_h = C0_h, CA_h = CA_h, I_h = I_h,
    S_II_h = S_II_h, C_II_h = C_II_h, I_II_h = I_II_h,
    S_III_h = S_III_h, C_III_h = C_III_h, I_III_h = I_III_h,
    # Community
    S0_c = S0_c, SA_c = SA_c, C0_c = C0_c, CA_c = CA_c, I_c = I_c,
    S_II_c = S_II_c, C_II_c = C_II_c, I_II_c = I_II_c,
    S_III_c = S_III_c, C_III_c = C_III_c, I_III_c = I_III_c))
}

# Create initial conditions for vaccination model (add CumI = 0)




###############################################################################
# 2. METRICS FUNCTIONS FOR CALIBRATION MODEL
###############################################################################

# compute R0 with Next Generation Matrix
compute_R0 <- function(params, N_h0, N_c0) {
  params <- as.list(params)

  with(params, {
    
    # ---- DFE susceptible split due to ATB dysbiosis (S0 <-> SA) ----
    # At DFE (no colonization/infection), S0 and SA are at equilibrium under tau/omega.
    den_h <- tau_h + omega_h
    den_c <- tau_c + omega_c
    
    f0_h <- if (den_h > 0) omega_h / den_h else 1
    fA_h <- if (den_h > 0) tau_h   / den_h else 0
    f0_c <- if (den_c > 0) omega_c / den_c else 1
    fA_c <- if (den_c > 0) tau_c   / den_c else 0
    
    # ---- Alpha at DFE: must be a fixed number ----
    # If not provided, use the DFE value implied by your dynamic formula with I=0:
    # out = delta * N_h0 ; den = w * N_c0  => alpha* = delta*N_h0/(w*N_c0)
    alpha_I   <- alpha_fixed * w_I
    alpha_II  <- alpha_fixed * w_II
    alpha_III <- alpha_fixed * w_III
    
    # ---- Infectious-state ordering (fixed) ----
    # 1 C0_h, 2 CA_h, 3 I_h, 4 CII_h, 5 I2_h, 6 CIII_h, 7 I3_h,
    # 8 C0_c, 9 CA_c,10 I_c,11 CII_c,12 I2_c,13 CIII_c,14 I3_c
    idx <- list(C0_h=1, CA_h=2, I_h=3, CII_h=4, I2_h=5, CIII_h=6, I3_h=7,
                C0_c=8, CA_c=9, I_c=10, CII_c=11, I2_c=12, CIII_c=13, I3_c=14)
    
    n <- 14
    F <- matrix(0, n, n)
    V <- matrix(0, n, n)
    
    # Give names (helps debugging)
    nm <- names(idx)
    rownames(F) <- colnames(F) <- nm
    rownames(V) <- colnames(V) <- nm
    
    # =========================================================
    # F: new infections only
    # New infections enter C0 via lambda*S0 and CA via lambda*SA.
    # lambda_h uses (C0+CA+CII+CIII + nu*(I+I2+I3)) / N
    # =========================================================
    
    hosp_C <- c(idx$C0_h, idx$CA_h, idx$CII_h, idx$CIII_h)
    hosp_I <- c(idx$I_h,  idx$I2_h,  idx$I3_h)
    comm_C <- c(idx$C0_c, idx$CA_c, idx$CII_c, idx$CIII_c)
    comm_I <- c(idx$I_c,  idx$I2_c,  idx$I3_c)
    
    # Hospital: new infections into C0_h (row 1) and CA_h (row 2)
    F[idx$C0_h, hosp_C] <- beta_h * f0_h
    F[idx$C0_h, hosp_I] <- beta_h * f0_h * nu_h
    F[idx$CA_h, hosp_C] <- beta_h * fA_h
    F[idx$CA_h, hosp_I] <- beta_h * fA_h * nu_h
    
    # Community: new infections into C0_c (row 8) and CA_c (row 9)
    F[idx$C0_c, comm_C] <- beta_c * f0_c
    F[idx$C0_c, comm_I] <- beta_c * f0_c * nu_c
    F[idx$CA_c, comm_C] <- beta_c * fA_c
    F[idx$CA_c, comm_I] <- beta_c * fA_c * nu_c
    
    # =========================================================
    # V: everything else (transfers + outflows)
    # Convention: x' = F x - V x
    # Diagonal: total outflow rate from compartment
    # Off-diagonal: - (rate) for inflow from another infected compartment
    # =========================================================
    
    # ---- Hospital block ----
    # C0_h: - (gamma + tau + sigma + delta) * C0_h + omega*CA_h + alpha*C0_c
    V[idx$C0_h, idx$C0_h] <- gamma_h + tau_h + sigma_h + delta
    V[idx$C0_h, idx$CA_h] <- -omega_h
    V[idx$C0_h, idx$C0_c] <- -alpha_fixed
    
    # CA_h: + tau*C0_h - (gamma + omega + kA*sigma + delta)*CA_h + alpha*CA_c
    V[idx$CA_h, idx$CA_h] <- gamma_h + omega_h + k_A * sigma_h + delta
    V[idx$CA_h, idx$C0_h] <- -tau_h
    V[idx$CA_h, idx$CA_c] <- -alpha_fixed
    
    # I_h: + sigma*C0_h + kA*sigma*CA_h - (epsilon + delta_I)*I_h + alpha_I*I_c
    V[idx$I_h, idx$I_h]   <- epsilon_h + delta_I
    V[idx$I_h, idx$C0_h]  <- -sigma_h
    V[idx$I_h, idx$CA_h]  <- -k_A * sigma_h
    V[idx$I_h, idx$I_c]   <- -alpha_I
    
    # CII_h: + (1-p)*epsilon*I_h - (gamma + kII*sigma + delta)*CII_h + alpha*CII_c
    V[idx$CII_h, idx$CII_h] <- gamma_h + k_II * sigma_h + delta
    V[idx$CII_h, idx$I_h]   <- -(1 - p_h) * epsilon_h
    V[idx$CII_h, idx$CII_c] <- -alpha_fixed
    
    # I2_h: + kII*sigma*CII_h - (epsilon + delta_II)*I2_h + alpha_II*I2_c
    V[idx$I2_h, idx$I2_h]   <- epsilon_h + delta_II
    V[idx$I2_h, idx$CII_h]  <- -k_II * sigma_h
    V[idx$I2_h, idx$I2_c]   <- -alpha_II
    
    # CIII_h: + (1-p)*epsilon*(I2_h + I3_h) - (gamma + kIII*sigma + delta)*CIII_h + alpha*CIII_c
    V[idx$CIII_h, idx$CIII_h] <- gamma_h + k_III * sigma_h + delta
    V[idx$CIII_h, idx$I2_h]   <- -(1 - p_h) * epsilon_h
    V[idx$CIII_h, idx$I3_h]   <- -(1 - p_h) * epsilon_h
    V[idx$CIII_h, idx$CIII_c] <- -alpha_fixed
    
    # I3_h: + kIII*sigma*CIII_h - (epsilon + delta_III)*I3_h + alpha_III*I3_c
    V[idx$I3_h, idx$I3_h]   <- epsilon_h + delta_III
    V[idx$I3_h, idx$CIII_h] <- -k_III * sigma_h
    V[idx$I3_h, idx$I3_c]   <- -alpha_III
    
    # ---- Community block ----
    # C0_c: - (gamma + tau + sigma + alpha)*C0_c + omega*CA_c + delta*C0_h
    V[idx$C0_c, idx$C0_c] <- gamma_c + tau_c + sigma_c + alpha_fixed
    V[idx$C0_c, idx$CA_c] <- -omega_c
    V[idx$C0_c, idx$C0_h] <- -delta
    
    # CA_c: + tau*C0_c - (gamma + omega + kA*sigma + alpha)*CA_c + delta*CA_h
    V[idx$CA_c, idx$CA_c] <- gamma_c + omega_c + k_A * sigma_c + alpha_fixed
    V[idx$CA_c, idx$C0_c] <- -tau_c
    V[idx$CA_c, idx$CA_h] <- -delta
    
    # I_c: + sigma*C0_c + kA*sigma*CA_c - (epsilon + alpha_I)*I_c
    V[idx$I_c, idx$I_c]   <- epsilon_c + alpha_I
    V[idx$I_c, idx$C0_c]  <- -sigma_c
    V[idx$I_c, idx$CA_c]  <- -k_A * sigma_c
    
    # CII_c: + (1-p)*epsilon*I_c - (gamma + kII*sigma + alpha)*CII_c + delta*CII_h + (1-prop)*delta_I*I_h
    V[idx$CII_c, idx$CII_c] <- gamma_c + k_II * sigma_c + alpha_fixed
    V[idx$CII_c, idx$I_c]   <- -(1 - p_c) * epsilon_c
    V[idx$CII_c, idx$CII_h] <- -delta
    V[idx$CII_c, idx$I_h]   <- -(1 - prop) * delta_I
    
    # I2_c: + kII*sigma*CII_c - (epsilon + alpha_II)*I2_c
    V[idx$I2_c, idx$I2_c]   <- epsilon_c + alpha_II
    V[idx$I2_c, idx$CII_c]  <- -k_II * sigma_c
    
    # CIII_c: + (1-p)*epsilon*(I2_c + I3_c) - (gamma + kIII*sigma + alpha)*CIII_c + delta*CIII_h
    #         + (1-prop)*(delta_II*I2_h + delta_III*I3_h)
    V[idx$CIII_c, idx$CIII_c] <- gamma_c + k_III * sigma_c + alpha_fixed
    V[idx$CIII_c, idx$I2_c]   <- -(1 - p_c) * epsilon_c
    V[idx$CIII_c, idx$I3_c]   <- -(1 - p_c) * epsilon_c
    V[idx$CIII_c, idx$CIII_h] <- -delta
    V[idx$CIII_c, idx$I2_h]   <- -(1 - prop) * delta_II
    V[idx$CIII_c, idx$I3_h]   <- -(1 - prop) * delta_III
    
    # I3_c: + kIII*sigma*CIII_c - (epsilon + alpha_III)*I3_c
    V[idx$I3_c, idx$I3_c]   <- epsilon_c + alpha_III
    V[idx$I3_c, idx$CIII_c] <- -k_III * sigma_c
    
    # ---- Next-generation matrix ----
    K <- F %*% solve(V)
    ev <- eigen(K, only.values = TRUE)$values
    R0 <- max(Mod(ev))
    
    list(R0 = R0, K = K, F = F, V = V,
         dfe_fractions = c(f0_h=f0_h, fA_h=fA_h, f0_c=f0_c, fA_c=fA_c),
         alpha_fixed = alpha_fixed)
  })
}


# Main metrics function for the calibration model
compute_metrics_calib <- function(out, params) {
  
  out <- as.data.frame(out) 
  last <- out[nrow(out), ]   # last row of out, at equilibrium 

  with(as.list(params), {
    
    # ---- Populations (Nh, Nc, Ntot) ----
    final_tot_h <- compute_totals(last$S0_h, last$C0_h, last$SA_h, last$CA_h, last$I_h, last$S_II_h, last$C_II_h, last$I_II_h, last$S_III_h, last$C_III_h, last$I_III_h)
    final_tot_c <- compute_totals(last$S0_c, last$C0_c, last$SA_c, last$CA_c, last$I_c, last$S_II_c, last$C_II_c, last$I_II_c, last$S_III_c, last$C_III_c, last$I_III_c)
    N_h <- final_tot_h$N
    N_c <- final_tot_c$N
    N_tot <- N_h + N_c
    
    # ---- Carriage prevalence ----
    C_h <- final_tot_h$C
    C_c <- final_tot_c$C
    prev_h <- C_h / N_h
    prev_c <- C_c / N_c
    
    # ---- Instantaneous incidence (total / primary / recurrent) ----
    inc_h_total <- sigma_h*last$C0_h + k_A*sigma_h*last$CA_h + k_II*sigma_h*last$C_II_h + k_III*sigma_h*last$C_III_h
    inc_h_primo <- sigma_h*last$C0_h + k_A*sigma_h*last$CA_h
    inc_h_rec   <- k_II*sigma_h*last$C_II_h + k_III*sigma_h*last$C_III_h
    
    inc_c_total <- sigma_c*last$C0_c + k_A*sigma_c*last$CA_c + k_II*sigma_c*last$C_II_c + k_III*sigma_c*last$C_III_c
    inc_c_primo <- sigma_c*last$C0_c + k_A*sigma_c*last$CA_c
    inc_c_rec   <- k_II*sigma_c*last$C_II_c + k_III*sigma_c*last$C_III_c
    
    # Normalized instantaneous incidence
    # - absolute: per N_tot
    # - relative: per setting population (N_h or N_c)
    inc_h_total_abs <- inc_h_total / N_tot
    inc_h_primo_abs <- inc_h_primo / N_tot
    inc_h_rec_abs   <- inc_h_rec / N_tot
    
    inc_c_total_abs <- inc_c_total / N_tot
    inc_c_primo_abs <- inc_c_primo / N_tot
    inc_c_rec_abs   <- inc_c_rec / N_tot
    
    inc_h_total_rel <- inc_h_total / N_h
    inc_h_primo_rel <- inc_h_primo / N_h
    inc_h_rec_rel   <- inc_h_rec / N_h
    
    inc_c_total_rel <- inc_c_total / N_c
    inc_c_primo_rel <- inc_c_primo / N_c
    inc_c_rec_rel   <- inc_c_rec / N_c
    
    # ---- Recurrence prevalence ----
    I1 <- last$I_h + last$I_c                                       # I primary (h+c)
    I2 <- last$I_II_h + last$I_II_c                                 # I recurrence 1 (h+c)
    I3 <- last$I_III_h + last$I_III_c                               # I recurrence 2+ (h+c)
    rec1 <- ifelse(I1 == 0, 0, I2 / I1)                             # recid_1
    rec2 <- ifelse(I2 == 0, 0, I3 / I2)                             # recid_2
    
    # ---- R0 (simple proxy time series) ----
    R0_df <- compute_R0(params, N_h0 = N_h, N_c0 = N_c)
    
    return(list(
      population = data.frame(time = last$time, N_h = N_h, N_c = N_c, N_tot = N_tot),
      
      carriage = data.frame(time = last$time, prev_h = prev_h, prev_c = prev_c),
      
      incidence_instant = data.frame(time = last$time,
        inc_h_total, inc_h_primo, inc_h_rec,
        inc_c_total, inc_c_primo, inc_c_rec,
        inc_h_total_abs, inc_h_primo_abs, inc_h_rec_abs,
        inc_c_total_abs, inc_c_primo_abs, inc_c_rec_abs,
        inc_h_total_rel, inc_h_primo_rel, inc_h_rec_rel,
        inc_c_total_rel, inc_c_primo_rel, inc_c_rec_rel),
      
      recurrence = data.frame(time = last$time, rec1 = rec1, rec2 = rec2),
      
      R0 = R0_df
    ))
  })
}





###############################################################################
# 3. METRICS FUNCTIONS FOR VACCINATION and ATB REDUCTION MODEL
###############################################################################

compute_metrics_scenario <- function(out, params) {
  
  out <- as.data.frame(out)
  last <- out[nrow(out), ]   # last row of out, at equilibrium 
  
  with(as.list(params), {
    
    # ---- Vaccinated parameters (same as in the ODE) ----
    sigma_h_v <- sigma_h * sigma_mult_v
    sigma_c_v <- sigma_c * sigma_mult_v
    
    # ---- Populations (N_h, N_c, N_tot)  ----
    final_tot_h_nv <- compute_totals(last$S0_h_nv, last$C0_h_nv, last$SA_h_nv, last$CA_h_nv, last$I_h_nv, last$S_II_h_nv, last$C_II_h_nv, last$I_II_h_nv, last$S_III_h_nv, last$C_III_h_nv, last$I_III_h_nv)
    final_tot_c_nv <- compute_totals(last$S0_c_nv, last$C0_c_nv, last$SA_c_nv, last$CA_c_nv, last$I_c_nv, last$S_II_c_nv, last$C_II_c_nv, last$I_II_c_nv, last$S_III_c_nv, last$C_III_c_nv, last$I_III_c_nv)
    final_tot_h_v  <- compute_totals(last$S0_h_v,  last$C0_h_v,  last$SA_h_v,  last$CA_h_v,  last$I_h_v,  last$S_II_h_v,  last$C_II_h_v,  last$I_II_h_v,  last$S_III_h_v,  last$C_III_h_v,  last$I_III_h_v)
    final_tot_c_v  <- compute_totals(last$S0_c_v,  last$C0_c_v,  last$SA_c_v,  last$CA_c_v,  last$I_c_v,  last$S_II_c_v,  last$C_II_c_v,  last$I_II_c_v,  last$S_III_c_v,  last$C_III_c_v,  last$I_III_c_v)

    N_h_nv  <- final_tot_h_nv$N
    N_c_nv  <- final_tot_c_nv$N
    N_h_v   <- final_tot_h_v$N
    N_c_v   <- final_tot_c_v$N
    N_h     <- N_h_nv + N_h_v
    N_c     <- N_c_nv + N_c_v
    N_tot   <- N_h + N_c
    N_nv    <- N_h_nv + N_c_nv
    N_v     <- N_h_v + N_c_v
    
    # ---- Carriage prevalence ----
    C_h <- final_tot_h_nv$C + final_tot_h_v$C
    C_c <- final_tot_c_nv$C + final_tot_c_v$C
    prev_h  <- C_h / N_h
    prev_c  <- C_c / N_c
    
    # ---- Instantaneous incidence ----
    # Non-vaccinated
    inc_h_nv_total <- sigma_h*last$C0_h_nv + k_A*sigma_h*last$CA_h_nv + k_II*sigma_h*last$C_II_h_nv + k_III*sigma_h*last$C_III_h_nv
    inc_h_nv_primo <- sigma_h*last$C0_h_nv + k_A*sigma_h*last$CA_h_nv
    inc_h_nv_rec   <- k_II*sigma_h*last$C_II_h_nv + k_III*sigma_h*last$C_III_h_nv
    
    inc_c_nv_total <- sigma_c*last$C0_c_nv + k_A*sigma_c*last$CA_c_nv + k_II*sigma_c*last$C_II_c_nv + k_III*sigma_c*last$C_III_c_nv
    inc_c_nv_primo <- sigma_c*last$C0_c_nv + k_A*sigma_c*last$CA_c_nv
    inc_c_nv_rec   <- k_II*sigma_c*last$C_II_c_nv + k_III*sigma_c*last$C_III_c_nv
    
    # Vaccinated (use sigma_h_v / sigma_c_v)
    inc_h_v_total <- sigma_h_v*last$C0_h_v + k_A*sigma_h_v*last$CA_h_v + k_II*sigma_h_v*last$C_II_h_v + k_III*sigma_h_v*last$C_III_h_v
    inc_h_v_primo <- sigma_h_v*last$C0_h_v + k_A*sigma_h_v*last$CA_h_v
    inc_h_v_rec   <- k_II*sigma_h_v*last$C_II_h_v + k_III*sigma_h_v*last$C_III_h_v
    
    inc_c_v_total <- sigma_c_v*last$C0_c_v + k_A*sigma_c_v*last$CA_c_v + k_II*sigma_c_v*last$C_II_c_v + k_III*sigma_c_v*last$C_III_c_v
    inc_c_v_primo <- sigma_c_v*last$C0_c_v + k_A*sigma_c_v*last$CA_c_v
    inc_c_v_rec   <- k_II*sigma_c_v*last$C_II_c_v + k_III*sigma_c_v*last$C_III_c_v
    
    # Aggregated by setting
    inc_h_total <- inc_h_nv_total + inc_h_v_total
    inc_h_primo <- inc_h_nv_primo + inc_h_v_primo
    inc_h_rec   <- inc_h_nv_rec   + inc_h_v_rec
    
    inc_c_total <- inc_c_nv_total + inc_c_v_total
    inc_c_primo <- inc_c_nv_primo + inc_c_v_primo
    inc_c_rec   <- inc_c_nv_rec   + inc_c_v_rec
    
    # Normalized instantaneous incidence
    # Absolute (per N_tot)
    inc_h_total_abs <- inc_h_total / N_tot
    inc_h_primo_abs <- inc_h_primo / N_tot
    inc_h_rec_abs   <- inc_h_rec   / N_tot
    inc_c_total_abs <- inc_c_total / N_tot
    inc_c_primo_abs <- inc_c_primo / N_tot
    inc_c_rec_abs   <- inc_c_rec   / N_tot
    # Relative (per setting population)
    inc_h_total_rel <- inc_h_total / N_h
    inc_h_primo_rel <- inc_h_primo / N_h
    inc_h_rec_rel   <- inc_h_rec   / N_h
    inc_c_total_rel <- inc_c_total / N_c
    inc_c_primo_rel <- inc_c_primo / N_c
    inc_c_rec_rel   <- inc_c_rec   / N_c
    
    # ---- Cumulative incidence (from CumI compartments) ----
    # Make cumulative start at 0 (difference from t=0)
    
    # Non-vaccinated
    cum_h_nv_total <- last$CumI_h_nv       - out$CumI_h_nv[1] # final value - initial value = cumul between the start and the end
    cum_h_nv_primo <- last$primo_CumI_h_nv - out$primo_CumI_h_nv[1]
    cum_h_nv_rec   <- last$rec_CumI_h_nv   - out$rec_CumI_h_nv[1]
    
    cum_c_nv_total <- last$CumI_c_nv       - out$CumI_c_nv[1]
    cum_c_nv_primo <- last$primo_CumI_c_nv - out$primo_CumI_c_nv[1]
    cum_c_nv_rec   <- last$rec_CumI_c_nv   - out$rec_CumI_c_nv[1]
    
    # Vaccinated
    cum_h_v_total <- last$CumI_h_v       - out$CumI_h_v[1]
    cum_h_v_primo <- last$primo_CumI_h_v - out$primo_CumI_h_v[1]
    cum_h_v_rec   <- last$rec_CumI_h_v   - out$rec_CumI_h_v[1]
    
    cum_c_v_total <- last$CumI_c_v       - out$CumI_c_v[1]
    cum_c_v_primo <- last$primo_CumI_c_v - out$primo_CumI_c_v[1]
    cum_c_v_rec   <- last$rec_CumI_c_v   - out$rec_CumI_c_v[1]
    
    # Aggregated by setting
    cum_h_total <- cum_h_nv_total + cum_h_v_total
    cum_h_primo <- cum_h_nv_primo + cum_h_v_primo
    cum_h_rec   <- cum_h_nv_rec   + cum_h_v_rec
    
    cum_c_total <- cum_c_nv_total + cum_c_v_total
    cum_c_primo <- cum_c_nv_primo + cum_c_v_primo
    cum_c_rec   <- cum_c_nv_rec   + cum_c_v_rec
    
    # Normalized cumulative incidence
    # Absolute (per N_tot)
    cum_h_total_abs <- cum_h_total / N_tot
    cum_h_primo_abs <- cum_h_primo / N_tot
    cum_h_rec_abs   <- cum_h_rec   / N_tot
    cum_c_total_abs <- cum_c_total / N_tot
    cum_c_primo_abs <- cum_c_primo / N_tot
    cum_c_rec_abs   <- cum_c_rec   / N_tot
    # Relative (per setting population)
    cum_h_total_rel <- cum_h_total / N_h
    cum_h_primo_rel <- cum_h_primo / N_h
    cum_h_rec_rel   <- cum_h_rec   / N_h
    cum_c_total_rel <- cum_c_total / N_c
    cum_c_primo_rel <- cum_c_primo / N_c
    cum_c_rec_rel   <- cum_c_rec   / N_c
    
    # ---- Recurrence prevalence (protected against div by 0) ----
    # Aggregated I across v/nv
    I1 <- last$I_h_nv + last$I_h_v + last$I_c_nv + last$I_c_v                    # I primary (h+c)
    I2 <- last$I_II_h_nv + last$I_II_h_v + last$I_II_c_nv  + last$I_II_c_v       # I recurrence 1 (h+c)
    I3 <- last$I_III_h_nv + last$I_III_h_v + last$I_III_c_nv + last$I_III_c_v    # I recurrence 2+ (h+c)
    rec1 <- ifelse(I1 == 0, 0, I2 / I1)                                      # recid_1
    rec2 <- ifelse(I2 == 0, 0, I3 / I2)                                      # recid_2
    
    # ---- Return ----
    return(list(
      population = data.frame(time = last$time, N_h = N_h, N_c = N_c, N_tot = N_tot, N_nv = N_nv, N_v = N_v),
      
      carriage = data.frame(time = last$time, prev_h = prev_h, prev_c = prev_c),
      
      incidence_instant = data.frame(time = last$time, 
                                     inc_h_total, inc_h_primo, inc_h_rec,
                                     inc_c_total, inc_c_primo, inc_c_rec,
                                     inc_h_total_abs, inc_h_primo_abs, inc_h_rec_abs,
                                     inc_c_total_abs, inc_c_primo_abs, inc_c_rec_abs,
                                     inc_h_total_rel, inc_h_primo_rel, inc_h_rec_rel,
                                     inc_c_total_rel, inc_c_primo_rel, inc_c_rec_rel),
      
      incidence_cumulative = data.frame(time = last$time,
                                        cum_h_total, cum_h_primo, cum_h_rec,
                                        cum_c_total, cum_c_primo, cum_c_rec,
                                        cum_h_total_abs, cum_h_primo_abs, cum_h_rec_abs,
                                        cum_c_total_abs, cum_c_primo_abs, cum_c_rec_abs,
                                        cum_h_total_rel, cum_h_primo_rel, cum_h_rec_rel,
                                        cum_c_total_rel, cum_c_primo_rel, cum_c_rec_rel),
      
      recurrence = data.frame(time = last$time, rec1 = rec1, rec2 = rec2)
      
    ))
  })
}






