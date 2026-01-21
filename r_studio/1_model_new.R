###############################################################################
###################### 1 : MODEL DEFINITION  ##################################
###############################################################################

###############################################################################
# ---- HELPER FUNCTIONS ----
###############################################################################

# Compute total population counts for a given setting (hospital or community)
compute_totals <- function(S0, SA, S_II, S_III, C0, CA, C_II, C_III, I, I_II, I_III) {
  list(S = S0 + SA + S_II + S_III, 
       C = C0 + CA + C_II + C_III, 
       I_tot = I + I_II + I_III,
       N = S0 + SA + S_II + S_III + C0 + CA + C_II + C_III + I + I_II + I_III)}

# Compute the force of infection (lambda)
compute_lambda <- function(beta, C, I_tot, N, nu) { beta * (C + nu * I_tot) / N }

# Compute stage-specific progression rates from sigmas
compute_sigmas <- function(sigma_base, k_A, k_II, k_III) {
  list(sigma_A = k_A * sigma_base, 
       sigma_II = k_II * sigma_base, 
       sigma_III = k_III * sigma_base)}

# Run the ODE model over a time grid (used to reach equilibrium) -> returns the lsoda output as a data.frame
run_model_to_equilibrium <- function(params_vec, init_cond, time_vec) {
  out <- as.data.frame(
    lsoda(
      y = init_cond,
      times = time_vec,
      func = cdiff_micro,
      parms = params_vec,
      atol = 1e-10,
      rtol = 1e-10,
      maxsteps = 50000
    )
  )
  return(out)
}

# Extract the final state (last row) from an ODE output (equilibrium state when the simulation is long enough)
get_equilibrium_state <- function(ode_result) {tail(ode_result, 1)}

###############################################################################
# ---- CLOSTRIDIUM DIFFICILE TRANSMISSION MODEL ----
###############################################################################

cdiff_micro <- function(t, pop, params) {
  with(as.list(c(pop, params)), {
    
    # Totals hospital / community 
    tot_h <- compute_totals(S0_h, SA_h, S_II_h, S_III_h, C0_h, CA_h, C_II_h, C_III_h, I_h, I_II_h, I_III_h)
    tot_c <- compute_totals(S0_c, SA_c, S_II_c, S_III_c, C0_c, CA_c, C_II_c, C_III_c, I_c, I_II_c, I_III_c)
    
    # ---- Dynamic renormalisation of alpha rates ----
    
    # Total hospital -> community outflow
    out_hc <- delta * (tot_h$S + tot_h$C) + delta_I * I_h + delta_II * I_II_h + delta_III * I_III_h
    
    # Denominator for renormalisation
    den_alpha <- w * (tot_c$S + tot_c$C) + w_I   * I_c + w_II  * I_II_c + w_III * I_III_c
    
    # Safety
    if (den_alpha <= 0) den_alpha <- 1
    
    # Base alpha
    alpha <- out_hc / den_alpha
    
    # Ordered admission rates
    alpha_I   <- alpha * w_I
    alpha_II  <- alpha * w_II
    alpha_III <- alpha * w_III
    
    # Forces of infection
    lambda_h <- compute_lambda(beta_h, tot_h$C, tot_h$I_tot, tot_h$N, nu)
    lambda_c <- compute_lambda(beta_c, tot_c$C, tot_c$I_tot, tot_c$N, nu)
    
    # Sigmas
    sig_h <- compute_sigmas(sigma_h, k_A, k_II, k_III)
    sig_c <- compute_sigmas(sigma_c, k_A, k_II, k_III)
    
    # EDO - Hospital
    dS0_h <- -lambda_h*S0_h + gamma*C0_h - tau_h*S0_h + omega_h*SA_h + phi*(S_II_h + S_III_h) + alpha*S0_c - delta*S0_h
    
    dSA_h <- -lambda_h*SA_h + gamma*CA_h + tau_h*S0_h - omega_h*SA_h +  alpha*SA_c - delta*SA_h
    
    dC0_h <- lambda_h*S0_h - gamma*C0_h - tau_h*C0_h + omega_h*CA_h - sigma_h*C0_h + alpha*C0_c - delta*C0_h
    
    dCA_h <- lambda_h*SA_h - gamma*CA_h + tau_h*C0_h - omega_h*CA_h - sig_h$sigma_A*CA_h + alpha*CA_c - delta*CA_h
    
    dI_h <- sigma_h*C0_h + sig_h$sigma_A*CA_h - epsilon*I_h + alpha_I*I_c - delta_I*I_h
    
    dS_II_h <- p*epsilon*I_h + gamma*C_II_h - lambda_h*S_II_h - phi*S_II_h + alpha*S_II_c - delta*S_II_h
    
    dC_II_h <- (1-p)*epsilon*I_h + lambda_h*S_II_h - gamma*C_II_h - sig_h$sigma_II*C_II_h + alpha*C_II_c - delta*C_II_h
    
    dI_II_h <- sig_h$sigma_II*C_II_h - epsilon*I_II_h + alpha_II*I_II_c - delta_II*I_II_h
    
    dS_III_h <- p*epsilon*(I_II_h + I_III_h) + gamma*C_III_h - lambda_h*S_III_h - phi*S_III_h + alpha*S_III_c - delta*S_III_h
    
    dC_III_h <- (1-p)*epsilon*(I_II_h + I_III_h) + lambda_h*S_III_h - gamma*C_III_h - sig_h$sigma_III*C_III_h + alpha*C_III_c - delta*C_III_h
    
    dI_III_h <- sig_h$sigma_III*C_III_h - epsilon*I_III_h + alpha_III*I_III_c - delta_III*I_III_h
    
    # EDO - Community
    dS0_c <- -lambda_c*S0_c + gamma*C0_c - tau_c*S0_c + omega_c*SA_c + phi*(S_II_c + S_III_c) - alpha*S0_c + delta*S0_h
    
    dSA_c <- -lambda_c*SA_c + gamma*CA_c + tau_c*S0_c - omega_c*SA_c - alpha*SA_c + delta*SA_h
    
    dC0_c <- lambda_c*S0_c - gamma*C0_c - tau_c*C0_c + omega_c*CA_c - sigma_c*C0_c - alpha*C0_c + delta*C0_h
    
    dCA_c <- lambda_c*SA_c - gamma*CA_c + tau_c*C0_c - omega_c*CA_c - sig_c$sigma_A*CA_c - alpha*CA_c + delta*CA_h
    
    dI_c <- sigma_c*C0_c + sig_c$sigma_A*CA_c - epsilon*I_c - alpha_I*I_c 
    
    dS_II_c <- p*epsilon*I_c + gamma*C_II_c - lambda_c*S_II_c - phi*S_II_c - alpha*S_II_c + delta*S_II_h + p*delta_I*I_h
    
    dC_II_c <- (1-p)*epsilon*I_c + lambda_c*S_II_c - gamma*C_II_c - sig_c$sigma_II*C_II_c - alpha*C_II_c + delta*C_II_h + (1-p)*delta_I*I_h
    
    dI_II_c <- sig_c$sigma_II*C_II_c - epsilon*I_II_c - alpha_II*I_II_c 
    
    dS_III_c <- p*epsilon*(I_II_c + I_III_c) + gamma*C_III_c - lambda_c*S_III_c - phi*S_III_c - alpha*S_III_c + delta*S_III_h + p*delta_II*I_II_h + p*delta_III*I_III_h
    
    dC_III_c <- (1-p)*epsilon*(I_II_c + I_III_c) + lambda_c*S_III_c - gamma*C_III_c - sig_c$sigma_III*C_III_c - alpha*C_III_c + delta*C_III_h  + (1-p)*delta_II*I_II_h + (1-p)*delta_III*I_III_h
    
    dI_III_c <- sig_c$sigma_III*C_III_c - epsilon*I_III_c - alpha_III*I_III_c 
    
    list(c(
      S0_h=dS0_h, SA_h=dSA_h, C0_h=dC0_h, CA_h=dCA_h, I_h=dI_h,
      S_II_h=dS_II_h, C_II_h=dC_II_h, I_II_h=dI_II_h,
      S_III_h=dS_III_h, C_III_h=dC_III_h, I_III_h=dI_III_h,
      S0_c=dS0_c, SA_c=dSA_c, C0_c=dC0_c, CA_c=dCA_c, I_c=dI_c,
      S_II_c=dS_II_c, C_II_c=dC_II_c, I_II_c=dI_II_c,
      S_III_c=dS_III_c, C_III_c=dC_III_c, I_III_c=dI_III_c))
  })
}

###############################################################################
# ---- PARAMETERS ----
###############################################################################

# see main

###############################################################################
# ---- INITIAL CONDITIONS (fonction) ----
###############################################################################

create_initial_conditions <- function(N_h, N_c, prev = 0.01, prev_rec = 0.005) {
  
  stopifnot(4*prev + 6*prev_rec < 1) # vérifie que la somme des prévalences initiales < 100% (S0 > 0)
  
  init_h <- list(SA_h = N_h * prev, C0_h = N_h * prev, CA_h = N_h * prev, I_h = N_h * prev,
                 S_II_h = N_h * prev_rec, C_II_h = N_h * prev_rec, I_II_h = N_h * prev_rec,
                 S_III_h = N_h * prev_rec, C_III_h = N_h * prev_rec, I_III_h = N_h * prev_rec)
  init_h$S0_h <- N_h - sum(unlist(init_h))
  
  init_c <- list(SA_c = N_c * prev, C0_c = N_c * prev, CA_c = N_c * prev, I_c = N_c * prev,
                 S_II_c = N_c * prev_rec, C_II_c = N_c * prev_rec, I_II_c = N_c * prev_rec,
                 S_III_c = N_c * prev_rec, C_III_c = N_c * prev_rec, I_III_c = N_c * prev_rec)
  init_c$S0_c <- N_c - sum(unlist(init_c))
  
  c(S0_h=init_h$S0_h, SA_h=init_h$SA_h, C0_h=init_h$C0_h, CA_h=init_h$CA_h, I_h=init_h$I_h,
    S_II_h=init_h$S_II_h, C_II_h=init_h$C_II_h, I_II_h=init_h$I_II_h,
    S_III_h=init_h$S_III_h, C_III_h=init_h$C_III_h, I_III_h=init_h$I_III_h,
    S0_c=init_c$S0_c, SA_c=init_c$SA_c, C0_c=init_c$C0_c, CA_c=init_c$CA_c, I_c=init_c$I_c,
    S_II_c=init_c$S_II_c, C_II_c=init_c$C_II_c, I_II_c=init_c$I_II_c,
    S_III_c=init_c$S_III_c, C_III_c=init_c$C_III_c, I_III_c=init_c$I_III_c)
}





