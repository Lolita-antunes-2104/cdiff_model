###############################################################################
################################ 1 : MODEL ####################################
###############################################################################

###############################################################################
# 1. HELPERS FUNCTIONS
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

###############################################################################
# 2. STRATIFIED MODEL HOSPITAL/COMMUNITY FOR CALIBRATION
###############################################################################

# This function defines the system of ODEs to be solved by lsoda() 
# specific function signature (t, pop, params) is required by lsoda
cdiff_model_for_calibration <- function(t, pop, params) {
  
  # Combine state variables (pop) and parameters (params) into a single list so that all compartments and parameters can be accessed by name directly
  with(as.list(c(pop, params)), {
    
    # ---- Totals (hospital/community) ----
    tot_h <- compute_totals(S0_h, C0_h, SA_h, CA_h, I_h, S_II_h, C_II_h, I_II_h, S_III_h, C_III_h, I_III_h)
    tot_c <- compute_totals(S0_c, C0_c, SA_c, CA_c, I_c, S_II_c, C_II_c, I_II_c, S_III_c, C_III_c, I_III_c)
    
    # ---- Force of infection (transmission rate S -> C) ----
    lambda_h <- compute_lambda(beta_h, tot_h$C, tot_h$I_tot, tot_h$N, nu_h)
    lambda_c <- compute_lambda(beta_c, tot_c$C, tot_c$I_tot, tot_c$N, nu_c)
    
    # ---- Dynamic alpha (hospital admissions) ----
    out_hc <- delta * (tot_h$S + tot_h$C) + delta_I * I_h + delta_II * I_II_h + delta_III * I_III_h
    den_alpha <- w * (tot_c$S + tot_c$C) + w_I * I_c + w_II * I_II_c + w_III * I_III_c
    alpha <- if (is.na(den_alpha) || den_alpha <= 0) 0 else out_hc / den_alpha
    
    # admissions rates by infection state
    alpha_I   <- alpha * w_I
    alpha_II  <- alpha * w_II
    alpha_III <- alpha * w_III
    
    #  ---- ODE HOSPITAL ----
    # Primary
    dS0_h <- -lambda_h*S0_h  + gamma_h*C0_h  - tau_h*S0_h  + omega_h*SA_h  + phi_h*(S_II_h + S_III_h)                 + alpha*S0_c - delta*S0_h
    dSA_h <- -lambda_h*SA_h  + gamma_h*CA_h  + tau_h*S0_h  - omega_h*SA_h                                             + alpha*SA_c - delta*SA_h
    dC0_h <-  lambda_h*S0_h  - gamma_h*C0_h  - tau_h*C0_h  + omega_h*CA_h  - sigma_h*C0_h                             + alpha*C0_c - delta*C0_h
    dCA_h <-  lambda_h*SA_h  - gamma_h*CA_h  + tau_h*C0_h  - omega_h*CA_h  - k_A*sigma_h*CA_h                         + alpha*CA_c - delta*CA_h
    dI_h  <-  sigma_h*C0_h   + k_A*sigma_h*CA_h  - epsilon_h*I_h                                                      + alpha_I*I_c - delta_I*I_h
    # First recurrence
    dS_II_h <-  p_h*epsilon_h*I_h     - lambda_h*S_II_h  + gamma_h*C_II_h  - phi_h*S_II_h                             + alpha*S_II_c - delta*S_II_h
    dC_II_h <- (1-p_h)*epsilon_h*I_h  + lambda_h*S_II_h  - gamma_h*C_II_h  - k_II*sigma_h*C_II_h                      + alpha*C_II_c - delta*C_II_h
    dI_II_h <-  k_II*sigma_h*C_II_h   - epsilon_h*I_II_h                                                              + alpha_II*I_II_c - delta_II*I_II_h
    # Second+ recurrence
    dS_III_h <-  p_h*epsilon_h*(I_II_h + I_III_h)     - lambda_h*S_III_h  + gamma_h*C_III_h  - phi_h*S_III_h          + alpha*S_III_c - delta*S_III_h
    dC_III_h <- (1-p_h)*epsilon_h*(I_II_h + I_III_h)  + lambda_h*S_III_h  - gamma_h*C_III_h  - k_III*sigma_h*C_III_h  + alpha*C_III_c - delta*C_III_h
    dI_III_h <-  k_III*sigma_h*C_III_h  - epsilon_h*I_III_h                                                           + alpha_III*I_III_c - delta_III*I_III_h
    
    #  ---- ODE COMMUNITY ----
    # Primary 
    dS0_c <- -lambda_c*S0_c  + gamma_c*C0_c  - tau_c*S0_c  + omega_c*SA_c  + phi_c*(S_II_c + S_III_c)                 - alpha*S0_c + delta*S0_h
    dSA_c <- -lambda_c*SA_c  + gamma_c*CA_c  + tau_c*S0_c  - omega_c*SA_c                                             - alpha*SA_c + delta*SA_h
    dC0_c <-  lambda_c*S0_c  - gamma_c*C0_c  - tau_c*C0_c  + omega_c*CA_c  - sigma_c*C0_c                             - alpha*C0_c + delta*C0_h
    dCA_c <-  lambda_c*SA_c  - gamma_c*CA_c  + tau_c*C0_c  - omega_c*CA_c  - k_A*sigma_c*CA_c                         - alpha*CA_c + delta*CA_h
    dI_c  <-  sigma_c*C0_c  + k_A*sigma_c*CA_c   - epsilon_c*I_c                                                      - alpha_I*I_c
    # First recurrence
    dS_II_c <-  p_c*epsilon_c*I_c     - lambda_c*S_II_c  + gamma_c*C_II_c  - phi_c*S_II_c                             - alpha*S_II_c + delta*S_II_h + prop*delta_I*I_h
    dC_II_c <- (1-p_c)*epsilon_c*I_c  + lambda_c*S_II_c  - gamma_c*C_II_c  - k_II*sigma_c*C_II_c                      - alpha*C_II_c + delta*C_II_h + (1 - prop)*delta_I*I_h
    dI_II_c <-  k_II*sigma_c*C_II_c   - epsilon_c*I_II_c                                                              - alpha_II*I_II_c
    # Second+ recurrence
    dS_III_c <-  p_c*epsilon_c*(I_II_c + I_III_c)     - lambda_c*S_III_c  + gamma_c*C_III_c  - phi_c*S_III_c          - alpha*S_III_c + delta*S_III_h + prop*(delta_II*I_II_h + delta_III*I_III_h)
    dC_III_c <- (1-p_c)*epsilon_c*(I_II_c + I_III_c)  + lambda_c*S_III_c  - gamma_c*C_III_c  - k_III*sigma_c*C_III_c  - alpha*C_III_c + delta*C_III_h + (1-prop)*(delta_II*I_II_h + delta_III*I_III_h)
    dI_III_c <-  k_III*sigma_c*C_III_c  - epsilon_c*I_III_c                                                           - alpha_III*I_III_c
    
    return(list(
      c(dS0_h, dSA_h, dC0_h, dCA_h, dI_h,
        dS_II_h, dC_II_h, dI_II_h,
        dS_III_h, dC_III_h, dI_III_h,
        
        dS0_c, dSA_c, dC0_c, dCA_c, dI_c,
        dS_II_c, dC_II_c, dI_II_c,
        dS_III_c, dC_III_c, dI_III_c)))
  })
}

###############################################################################
# 3. STRATIFIED MODEL HOSPITAL/COMMUNITY/VACCINATION FOR SCENARIO
###############################################################################

cdiff_model_for_scenario <- function(t, pop, params) {
  
  # Combine state variables (pop) and parameters (params) into a single list so that all compartments and parameters can be accessed by name directly
  with(as.list(c(pop, params)), {
    
    # ---- Totals (hospital/community and vaccinated/non-vaccinated) ----
    tot_h_nv <- compute_totals(S0_h_nv, C0_h_nv, SA_h_nv, CA_h_nv, I_h_nv, S_II_h_nv, C_II_h_nv, I_II_h_nv, S_III_h_nv, C_III_h_nv, I_III_h_nv)
    tot_c_nv <- compute_totals(S0_c_nv, C0_c_nv, SA_c_nv, CA_c_nv, I_c_nv, S_II_c_nv, C_II_c_nv, I_II_c_nv, S_III_c_nv, C_III_c_nv, I_III_c_nv)
    tot_h_v <- compute_totals(S0_h_v, C0_h_v, SA_h_v, CA_h_v, I_h_v, S_II_h_v, C_II_h_v, I_II_h_v, S_III_h_v, C_III_h_v, I_III_h_v)
    tot_c_v <- compute_totals(S0_c_v, C0_c_v, SA_c_v, CA_c_v, I_c_v, S_II_c_v, C_II_c_v, I_II_c_v, S_III_c_v, C_III_c_v, I_III_c_v)
    
    C_h <- tot_h_nv$C + tot_h_v$C
    C_c <- tot_c_nv$C + tot_c_v$C
    
    I_tot_h <- tot_h_nv$I + tot_h_v$I
    I_tot_c <- tot_c_nv$I + tot_c_v$I
    
    N_h <- tot_h_nv$N + tot_h_v$N 
    N_c <- tot_c_nv$N + tot_c_v$N 
  
    # ---- Force of infection (shared across vaccinated and non-vaccinated, but separated for hospital and community) ----
    lambda_h <- compute_lambda(beta_h, C_h, I_tot_h, N_tot_h, nu_h)
    lambda_c <- compute_lambda(beta_c, C_c, I_tot_c, N_tot_c, nu_c)
    
    # ---- Fixed alpha (hospital admissions) ----
    
    # ---- Vaccinated et ATB parameters ----
    sigma_h_v <- sigma_h * sigma_mult_v
    sigma_c_v <- sigma_c * sigma_mult_v
    
    tau_h_red <- tau_h * tau_mult_red
    tau_c_red <- tau_c * tau_mult_red
    
    # ---- ODE HOSPITAL / NON VACCINATED ----
    # Primary
    dS0_h_nv <- -lambda_h*S0_h_nv  + gamma_h*C0_h_nv  - tau_h*S0_h_nv  + omega_h*SA_h_nv  + phi_h*(S_II_h_nv + S_III_h_nv)             + alpha*S0_c_nv - delta*S0_h_nv
    dSA_h_nv <- -lambda_h*SA_h_nv  + gamma_h*CA_h_nv  + tau_h*S0_h_nv  - omega_h*SA_h_nv                                               + alpha*SA_c_nv - delta*SA_h_nv
    dC0_h_nv <-  lambda_h*S0_h_nv  - gamma_h*C0_h_nv  - tau_h*C0_h_nv  + omega_h*CA_h_nv  - sigma_h*C0_h_nv                            + alpha*C0_c_nv - delta*C0_h_nv
    dCA_h_nv <-  lambda_h*SA_h_nv  - gamma_h*CA_h_nv  + tau_h*C0_h_nv  - omega_h*CA_h_nv  - k_A*sigma_h*CA_h_nv                        + alpha*CA_c_nv - delta*CA_h_nv
    dI_h_nv  <-  sigma_h*C0_h_nv  + k_A*sigma_h*CA_h_nv  - epsilon_h*I_h_nv                                                            + alpha_I*I_c_nv - delta_I*I_h_nv
    # First recurrence
    dS_II_h_nv <-  p_h*epsilon_h*I_h_nv     - lambda_h*S_II_h_nv  + gamma_h*C_II_h_nv  - phi_h*S_II_h_nv                               + alpha*S_II_c_nv - delta*S_II_h_nv
    dC_II_h_nv <- (1-p_h)*epsilon_h*I_h_nv  + lambda_h*S_II_h_nv  - gamma_h*C_II_h_nv  - k_II*sigma_h*C_II_h_nv                        + alpha*C_II_c_nv - delta*C_II_h_nv
    dI_II_h_nv <-  k_II*sigma_h*C_II_h_nv   - epsilon_h*I_II_h_nv                                                                      + alpha_II*I_II_c_nv - delta_II*I_II_h_nv
    # Second+ recurrence
    dS_III_h_nv <-  p_h*epsilon_h*(I_II_h_nv + I_III_h_nv)    - lambda_h*S_III_h_nv  + gamma_h*C_III_h_nv  - phi_h*S_III_h_nv          + alpha*S_III_c_nv - delta*S_III_h_nv
    dC_III_h_nv <- (1-p_h)*epsilon_h*(I_II_h_nv + I_III_h_nv) + lambda_h*S_III_h_nv  - gamma_h*C_III_h_nv  - k_III*sigma_h*C_III_h_nv  + alpha*C_III_c_nv - delta*C_III_h_nv
    dI_III_h_nv <-  k_III*sigma_h*C_III_h_nv  - epsilon_h*I_III_h_nv                                                                   + alpha_III*I_III_c_nv - delta_III*I_III_h_nv

    # ---- ODE COMMUNITY / NON VACCINATED ----
    # Primary
    dS0_c_nv <- -lambda_c*S0_c_nv  + gamma_c*C0_c_nv  - tau_c*S0_c_nv  + omega_c*SA_c_nv  + phi_c*(S_II_c_nv + S_III_c_nv)             - alpha*S0_c_nv + delta*S0_h_nv
    dSA_c_nv <- -lambda_c*SA_c_nv  + gamma_c*CA_c_nv  + tau_c*S0_c_nv  - omega_c*SA_c_nv                                               - alpha*SA_c_nv + delta*SA_h_nv
    dC0_c_nv <-  lambda_c*S0_c_nv  - gamma_c*C0_c_nv  - tau_c*C0_c_nv  + omega_c*CA_c_nv  - sigma_c*C0_c_nv                            - alpha*C0_c_nv + delta*C0_h_nv
    dCA_c_nv <-  lambda_c*SA_c_nv  - gamma_c*CA_c_nv  + tau_c*C0_c_nv  - omega_c*CA_c_nv  - k_A*sigma_c*CA_c_nv                        - alpha*CA_c_nv + delta*CA_h_nv
    dI_c_nv  <-  sigma_c*C0_c_nv  + k_A*sigma_c*CA_c_nv - epsilon_c*I_c_nv                                                             - alpha_I*I_c_nv
    # First recurrence
    dS_II_c_nv <-  p_c*epsilon_c*I_c_nv     - lambda_c*S_II_c_nv  + gamma_c*C_II_c_nv  - phi_c*S_II_c_nv                               - alpha*S_II_c_nv + delta*S_II_h_nv + prop*delta_I*I_h_nv
    dC_II_c_nv <- (1-p_c)*epsilon_c*I_c_nv  + lambda_c*S_II_c_nv  - gamma_c*C_II_c_nv  - k_II*sigma_c*C_II_c_nv                        - alpha*C_II_c_nv + delta*C_II_h_nv + (1 - prop)*delta_I*I_h_nv
    dI_II_c_nv <-  k_II*sigma_c*C_II_c_nv   - epsilon_c*I_II_c_nv                                                                      - alpha_II*I_II_c_nv
    # Second+ recurrence
    dS_III_c_nv <-  p_c*epsilon_c*(I_II_c_nv + I_III_c_nv)    - lambda_c*S_III_c_nv  + gamma_c*C_III_c_nv  - phi_c*S_III_c_nv          - alpha*S_III_c_nv + delta*S_III_h_nv + prop*(delta_II*I_II_h_nv + delta_III*I_III_h_nv)
    dC_III_c_nv <- (1-p_c)*epsilon_c*(I_II_c_nv + I_III_c_nv) + lambda_c*S_III_c_nv  - gamma_c*C_III_c_nv  - k_III*sigma_c*C_III_c_nv  - alpha*C_III_c_nv + delta*C_III_h_nv + (1-prop)*(delta_II*I_II_h_nv + delta_III*I_III_h_nv)
    dI_III_c_nv <-  k_III*sigma_c*C_III_c_nv  - epsilon_c*I_III_c_nv                                                                   - alpha_III*I_III_c_nv
    
    # ---- HOSPITAL / VACCINATED ----
    # Primary
    dS0_h_v <- -lambda_h*S0_h_v  + gamma_h*C0_h_v  - tau_h*S0_h_v  + omega_h*SA_h_v  + phi_h*(S_II_h_v + S_III_h_v)              + alpha*S0_c_v - delta*S0_h_v
    dSA_h_v <- -lambda_h*SA_h_v  + gamma_h*CA_h_v  + tau_h*S0_h_v  - omega_h*SA_h_v                                              + alpha*SA_c_v - delta*SA_h_v
    dC0_h_v <-  lambda_h*S0_h_v  - gamma_h*C0_h_v  - tau_h*C0_h_v  + omega_h*CA_h_v  - sigma_h*C0_h_v                            + alpha*C0_c_v - delta*C0_h_v
    dCA_h_v <-  lambda_h*SA_h_v  - gamma_h*CA_h_v  + tau_h*C0_h_v  - omega_h*CA_h_v  - k_A*sigma_h*CA_h_v                        + alpha*CA_c_v - delta*CA_h_v
    dI_h_v  <-  sigma_h*C0_h_v  + k_A*sigma_h*CA_h_v  - epsilon_h*I_h_v                                                          + alpha_I*I_c_v - delta_I*I_h_v
    # First recurrence
    dS_II_h_v <-  p_h*epsilon_h*I_h_v     - lambda_h*S_II_h_v  + gamma_h*C_II_h_v  - phi_h*S_II_h_v                              + alpha*S_II_c_v - delta*S_II_h_v
    dC_II_h_v <- (1-p_h)*epsilon_h*I_h_v  + lambda_h*S_II_h_v  - gamma_h*C_II_h_v  - k_II*sigma_h*C_II_h_v                       + alpha*C_II_c_v - delta*C_II_h_v
    dI_II_h_v <-  k_II*sigma_h*C_II_h_v   - epsilon_h*I_II_h_v                                                                   + alpha_II*I_II_c_v - delta_II*I_II_h_v
    # Second+ recurrence
    dS_III_h_v <-  p_h*epsilon_h*(I_II_h_v + I_III_h_v)    - lambda_h*S_III_h_v  + gamma_h*C_III_h_v  - phi_h*S_III_h_v          + alpha*S_III_c_v - delta*S_III_h_v
    dC_III_h_v <- (1-p_h)*epsilon_h*(I_II_h_v + I_III_h_v) + lambda_h*S_III_h_v  - gamma_h*C_III_h_v  - k_III*sigma_h*C_III_h_v  + alpha*C_III_c_v - delta*C_III_h_v
    dI_III_h_v <-  k_III*sigma_h*C_III_h_v  - epsilon_h*I_III_h_v                                                                + alpha_III*I_III_c_v - delta_III*I_III_h_v
    
    # ---- COMMUNITY / VACCINATED ----
    # Primary compartments
    dS0_c_v <- -lambda_c*S0_c_v  + gamma_c*C0_c_v  - tau_c*S0_c_v  + omega_c*SA_c_v  + phi_c*(S_II_c_v + S_III_c_v)              - alpha*S0_c_v + delta*S0_h_v
    dSA_c_v <- -lambda_c*SA_c_v  + gamma_c*CA_c_v  + tau_c*S0_c_v  - omega_c*SA_c_v                                              - alpha*SA_c_v + delta*SA_h_v
    dC0_c_v <-  lambda_c*S0_c_v  - gamma_c*C0_c_v  - tau_c*C0_c_v  + omega_c*CA_c_v  - sigma_c*C0_c_v                            - alpha*C0_c_v + delta*C0_h_v
    dCA_c_v <-  lambda_c*SA_c_v  - gamma_c*CA_c_v  + tau_c*C0_c_v  - omega_c*CA_c_v  - k_A*sigma_c*CA_c_v                        - alpha*CA_c_v + delta*CA_h_v
    dI_c_v  <-  sigma_c*C0_c_v  + k_A*sigma_c*CA_c_v - epsilon_c*I_c_v                                                           - alpha_I*I_c_v
    # First recurrence
    dS_II_c_v <-  p_c*epsilon_c*I_c_v     - lambda_c*S_II_c_v  + gamma_c*C_II_c_v  - phi_c*S_II_c_v                              - alpha*S_II_c_v + delta*S_II_h_v + prop*delta_I*I_h_v
    dC_II_c_v <- (1-p_c)*epsilon_c*I_c_v  + lambda_c*S_II_c_v  - gamma_c*C_II_c_v  - k_II*sigma_c*C_II_c_v                       - alpha*C_II_c_v + delta*C_II_h_v + (1 - prop)*delta_I*I_h_v
    dI_II_c_v <-  k_II*sigma_c*C_II_c_v   - epsilon_c*I_II_c_v                                                                   - alpha_II*I_II_c_v
    # Second+ recurrence
    dS_III_c_v <-  p_c*epsilon_c*(I_II_c_v + I_III_c_v)    - lambda_c*S_III_c_v  + gamma_c*C_III_c_v  - phi_c*S_III_c_v          - alpha*S_III_c_v + delta*S_III_h_v + prop*(delta_II*I_II_h_v + delta_III*I_III_h_v)
    dC_III_c_v <- (1-p_c)*epsilon_c*(I_II_c_v + I_III_c_v) + lambda_c*S_III_c_v  - gamma_c*C_III_c_v  - k_III*sigma_c*C_III_c_v  - alpha*C_III_c_v + delta*C_III_h_v + (1-prop)*(delta_II*I_II_h_v + delta_III*I_III_h_v)
    dI_III_c_v <-  k_III*sigma_c*C_III_c_v  - epsilon_c*I_III_c_v                                                                - alpha_III*I_III_c_v
    
    return(list(
      c(
        dS0_h_nv, dSA_h_nv, dC0_h_nv, dCA_h_nv, dI_h_nv,
        dS_II_h_nv, dC_II_h_nv, dI_II_h_nv,
        dS_III_h_nv, dC_III_h_nv, dI_III_h_nv,
        
        dS0_c_nv, dSA_c_nv, dC0_c_nv, dCA_c_nv, dI_c_nv,
        dS_II_c_nv, dC_II_c_nv, dI_II_c_nv,
        dS_III_c_nv, dC_III_c_nv, dI_III_c_nv,
        
        dS0_h_v, dSA_h_v, dC0_h_v, dCA_h_v, dI_h_v,
        dS_II_h_v, dC_II_h_v, dI_II_h_v,
        dS_III_h_v, dC_III_h_v, dI_III_h_v,
        
        dS0_c_v, dSA_c_v, dC0_c_v, dCA_c_v, dI_c_v,
        dS_II_c_v, dC_II_c_v, dI_II_c_v,
        dS_III_c_v, dC_III_c_v, dI_III_c_v)
    ))
  })
}



