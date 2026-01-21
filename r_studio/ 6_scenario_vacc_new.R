###############################################################################
#################### 6 : VACCINATION SCENARIOS ################################
###############################################################################

###############################################################################
# ---- VACCINATION MODEL (stratified by vaccination status) ----
###############################################################################

cdiff_vacc <- function(t, pop, params) {
  with(as.list(c(pop, params)), {
    
    # ---- Totals for non-vaccinated (nv) ----
    tot_h_nv <- compute_totals(
      S0_h_nv, SA_h_nv, S_II_h_nv, S_III_h_nv,
      C0_h_nv, CA_h_nv, C_II_h_nv, C_III_h_nv,
      I_h_nv, I_II_h_nv, I_III_h_nv
    )
    
    tot_c_nv <- compute_totals(
      S0_c_nv, SA_c_nv, S_II_c_nv, S_III_c_nv,
      C0_c_nv, CA_c_nv, C_II_c_nv, C_III_c_nv,
      I_c_nv, I_II_c_nv, I_III_c_nv
    )
    
    # ---- Totals for vaccinated (v) ----
    tot_h_v <- compute_totals(
      S0_h_v, SA_h_v, S_II_h_v, S_III_h_v,
      C0_h_v, CA_h_v, C_II_h_v, C_III_h_v,
      I_h_v, I_II_h_v, I_III_h_v
    )
    
    tot_c_v <- compute_totals(
      S0_c_v, SA_c_v, S_II_c_v, S_III_c_v,
      C0_c_v, CA_c_v, C_II_c_v, C_III_c_v,
      I_c_v, I_II_c_v, I_III_c_v
    )
    
    # ---- Dynamic alpha for NON-VACCINATED ----
    out_hc_nv <- delta * (tot_h_nv$S + tot_h_nv$C) +
      delta_I   * I_h_nv +
      delta_II  * I_II_h_nv +
      delta_III * I_III_h_nv
    
    den_alpha_nv <- w * (tot_c_nv$S + tot_c_nv$C) +
      w_I   * I_c_nv +
      w_II  * I_II_c_nv +
      w_III * I_III_c_nv
    
    if (den_alpha_nv <= 0) den_alpha_nv <- 1
    alpha_nv <- out_hc_nv / den_alpha_nv
    alpha_I_nv   <- alpha_nv * w_I
    alpha_II_nv  <- alpha_nv * w_II
    alpha_III_nv <- alpha_nv * w_III
    
    # ---- Dynamic alpha for VACCINATED ----
    out_hc_v <- delta * (tot_h_v$S + tot_h_v$C) +
      delta_I   * I_h_v +
      delta_II  * I_II_h_v +
      delta_III * I_III_h_v
    
    den_alpha_v <- w * (tot_c_v$S + tot_c_v$C) +
      w_I   * I_c_v +
      w_II  * I_II_c_v +
      w_III * I_III_c_v
    
    if (den_alpha_v <= 0) den_alpha_v <- 1
    alpha_v <- out_hc_v / den_alpha_v
    alpha_I_v   <- alpha_v * w_I
    alpha_II_v  <- alpha_v * w_II
    alpha_III_v <- alpha_v * w_III
    
    # ---- Common force of infection (FOI) ----
    # Hospital: total C and I across v + nv
    C_h_total <- tot_h_nv$C + tot_h_v$C
    I_h_total <- tot_h_nv$I_tot + tot_h_v$I_tot
    N_h_total <- tot_h_nv$N + tot_h_v$N
    
    lambda_h <- compute_lambda(beta_h, C_h_total, I_h_total, N_h_total, nu)
    
    # Community: total C and I across v + nv
    C_c_total <- tot_c_nv$C + tot_c_v$C
    I_c_total <- tot_c_nv$I_tot + tot_c_v$I_tot
    N_c_total <- tot_c_nv$N + tot_c_v$N
    
    lambda_c <- compute_lambda(beta_c, C_c_total, I_c_total, N_c_total, nu)
    
    # ---- Sigmas (with vaccine effect for vaccinated) ----
    # Non-vaccinated
    sig_h_nv <- compute_sigmas(sigma_h, k_A, k_II, k_III)
    sig_c_nv <- compute_sigmas(sigma_c, k_A, k_II, k_III)
    
    # Vaccinated (reduced by VE)
    sigma_h_v <- sigma_h * (1 - VE)
    sigma_c_v <- sigma_c * (1 - VE)
    sig_h_v <- compute_sigmas(sigma_h_v, k_A, k_II, k_III)
    sig_c_v <- compute_sigmas(sigma_c_v, k_A, k_II, k_III)
    
    # ===============================================
    # EDO - NON-VACCINATED POPULATION
    # ===============================================
    
    # Hospital - non-vaccinated
    dS0_h_nv <- -lambda_h*S0_h_nv + gamma*C0_h_nv - tau_h*S0_h_nv + omega_h*SA_h_nv + 
      phi*(S_II_h_nv + S_III_h_nv) + alpha_nv*S0_c_nv - delta*S0_h_nv
    
    dSA_h_nv <- -lambda_h*SA_h_nv + gamma*CA_h_nv + tau_h*S0_h_nv - omega_h*SA_h_nv + 
      alpha_nv*SA_c_nv - delta*SA_h_nv
    
    dC0_h_nv <- lambda_h*S0_h_nv - gamma*C0_h_nv - tau_h*C0_h_nv + omega_h*CA_h_nv - 
      sigma_h*C0_h_nv + alpha_nv*C0_c_nv - delta*C0_h_nv
    
    dCA_h_nv <- lambda_h*SA_h_nv - gamma*CA_h_nv + tau_h*C0_h_nv - omega_h*CA_h_nv - 
      sig_h_nv$sigma_A*CA_h_nv + alpha_nv*CA_c_nv - delta*CA_h_nv
    
    dI_h_nv <- sigma_h*C0_h_nv + sig_h_nv$sigma_A*CA_h_nv - epsilon*I_h_nv + 
      alpha_I_nv*I_c_nv - delta_I*I_h_nv
    
    dS_II_h_nv <- p*epsilon*I_h_nv + gamma*C_II_h_nv - lambda_h*S_II_h_nv - phi*S_II_h_nv + 
      alpha_nv*S_II_c_nv - delta*S_II_h_nv
    
    dC_II_h_nv <- (1-p)*epsilon*I_h_nv + lambda_h*S_II_h_nv - gamma*C_II_h_nv - 
      sig_h_nv$sigma_II*C_II_h_nv + alpha_nv*C_II_c_nv - delta*C_II_h_nv
    
    dI_II_h_nv <- sig_h_nv$sigma_II*C_II_h_nv - epsilon*I_II_h_nv + 
      alpha_II_nv*I_II_c_nv - delta_II*I_II_h_nv
    
    dS_III_h_nv <- p*epsilon*(I_II_h_nv + I_III_h_nv) + gamma*C_III_h_nv - lambda_h*S_III_h_nv - 
      phi*S_III_h_nv + alpha_nv*S_III_c_nv - delta*S_III_h_nv
    
    dC_III_h_nv <- (1-p)*epsilon*(I_II_h_nv + I_III_h_nv) + lambda_h*S_III_h_nv - gamma*C_III_h_nv - 
      sig_h_nv$sigma_III*C_III_h_nv + alpha_nv*C_III_c_nv - delta*C_III_h_nv
    
    dI_III_h_nv <- sig_h_nv$sigma_III*C_III_h_nv - epsilon*I_III_h_nv + 
      alpha_III_nv*I_III_c_nv - delta_III*I_III_h_nv
    
    # Community - non-vaccinated
    dS0_c_nv <- -lambda_c*S0_c_nv + gamma*C0_c_nv - tau_c*S0_c_nv + omega_c*SA_c_nv + 
      phi*(S_II_c_nv + S_III_c_nv) - alpha_nv*S0_c_nv + delta*S0_h_nv
    
    dSA_c_nv <- -lambda_c*SA_c_nv + gamma*CA_c_nv + tau_c*S0_c_nv - omega_c*SA_c_nv - 
      alpha_nv*SA_c_nv + delta*SA_h_nv
    
    dC0_c_nv <- lambda_c*S0_c_nv - gamma*C0_c_nv - tau_c*C0_c_nv + omega_c*CA_c_nv - 
      sigma_c*C0_c_nv - alpha_nv*C0_c_nv + delta*C0_h_nv
    
    dCA_c_nv <- lambda_c*SA_c_nv - gamma*CA_c_nv + tau_c*C0_c_nv - omega_c*CA_c_nv - 
      sig_c_nv$sigma_A*CA_c_nv - alpha_nv*CA_c_nv + delta*CA_h_nv
    
    dI_c_nv <- sigma_c*C0_c_nv + sig_c_nv$sigma_A*CA_c_nv - epsilon*I_c_nv - alpha_I_nv*I_c_nv
    
    dS_II_c_nv <- p*epsilon*I_c_nv + gamma*C_II_c_nv - lambda_c*S_II_c_nv - phi*S_II_c_nv - 
      alpha_nv*S_II_c_nv + delta*S_II_h_nv + p*delta_I*I_h_nv
    
    dC_II_c_nv <- (1-p)*epsilon*I_c_nv + lambda_c*S_II_c_nv - gamma*C_II_c_nv - 
      sig_c_nv$sigma_II*C_II_c_nv - alpha_nv*C_II_c_nv + delta*C_II_h_nv + (1-p)*delta_I*I_h_nv
    
    dI_II_c_nv <- sig_c_nv$sigma_II*C_II_c_nv - epsilon*I_II_c_nv - alpha_II_nv*I_II_c_nv
    
    dS_III_c_nv <- p*epsilon*(I_II_c_nv + I_III_c_nv) + gamma*C_III_c_nv - lambda_c*S_III_c_nv - 
      phi*S_III_c_nv - alpha_nv*S_III_c_nv + delta*S_III_h_nv + p*delta_II*I_II_h_nv + p*delta_III*I_III_h_nv
    
    dC_III_c_nv <- (1-p)*epsilon*(I_II_c_nv + I_III_c_nv) + lambda_c*S_III_c_nv - gamma*C_III_c_nv - 
      sig_c_nv$sigma_III*C_III_c_nv - alpha_nv*C_III_c_nv + delta*C_III_h_nv + 
      (1-p)*delta_II*I_II_h_nv + (1-p)*delta_III*I_III_h_nv
    
    dI_III_c_nv <- sig_c_nv$sigma_III*C_III_c_nv - epsilon*I_III_c_nv - alpha_III_nv*I_III_c_nv
    
    # ===============================================
    # EDO - VACCINATED POPULATION
    # ===============================================
    
    # Hospital - vaccinated
    dS0_h_v <- -lambda_h*S0_h_v + gamma*C0_h_v - tau_h*S0_h_v + omega_h*SA_h_v + 
      phi*(S_II_h_v + S_III_h_v) + alpha_v*S0_c_v - delta*S0_h_v
    
    dSA_h_v <- -lambda_h*SA_h_v + gamma*CA_h_v + tau_h*S0_h_v - omega_h*SA_h_v + 
      alpha_v*SA_c_v - delta*SA_h_v
    
    dC0_h_v <- lambda_h*S0_h_v - gamma*C0_h_v - tau_h*C0_h_v + omega_h*CA_h_v - 
      sigma_h_v*C0_h_v + alpha_v*C0_c_v - delta*C0_h_v
    
    dCA_h_v <- lambda_h*SA_h_v - gamma*CA_h_v + tau_h*C0_h_v - omega_h*CA_h_v - 
      sig_h_v$sigma_A*CA_h_v + alpha_v*CA_c_v - delta*CA_h_v
    
    dI_h_v <- sigma_h_v*C0_h_v + sig_h_v$sigma_A*CA_h_v - epsilon*I_h_v + 
      alpha_I_v*I_c_v - delta_I*I_h_v
    
    dS_II_h_v <- p*epsilon*I_h_v + gamma*C_II_h_v - lambda_h*S_II_h_v - phi*S_II_h_v + 
      alpha_v*S_II_c_v - delta*S_II_h_v
    
    dC_II_h_v <- (1-p)*epsilon*I_h_v + lambda_h*S_II_h_v - gamma*C_II_h_v - 
      sig_h_v$sigma_II*C_II_h_v + alpha_v*C_II_c_v - delta*C_II_h_v
    
    dI_II_h_v <- sig_h_v$sigma_II*C_II_h_v - epsilon*I_II_h_v + 
      alpha_II_v*I_II_c_v - delta_II*I_II_h_v
    
    dS_III_h_v <- p*epsilon*(I_II_h_v + I_III_h_v) + gamma*C_III_h_v - lambda_h*S_III_h_v - 
      phi*S_III_h_v + alpha_v*S_III_c_v - delta*S_III_h_v
    
    dC_III_h_v <- (1-p)*epsilon*(I_II_h_v + I_III_h_v) + lambda_h*S_III_h_v - gamma*C_III_h_v - 
      sig_h_v$sigma_III*C_III_h_v + alpha_v*C_III_c_v - delta*C_III_h_v
    
    dI_III_h_v <- sig_h_v$sigma_III*C_III_h_v - epsilon*I_III_h_v + 
      alpha_III_v*I_III_c_v - delta_III*I_III_h_v
    
    # Community - vaccinated
    dS0_c_v <- -lambda_c*S0_c_v + gamma*C0_c_v - tau_c*S0_c_v + omega_c*SA_c_v + 
      phi*(S_II_c_v + S_III_c_v) - alpha_v*S0_c_v + delta*S0_h_v
    
    dSA_c_v <- -lambda_c*SA_c_v + gamma*CA_c_v + tau_c*S0_c_v - omega_c*SA_c_v - 
      alpha_v*SA_c_v + delta*SA_h_v
    
    dC0_c_v <- lambda_c*S0_c_v - gamma*C0_c_v - tau_c*C0_c_v + omega_c*CA_c_v - 
      sigma_c_v*C0_c_v - alpha_v*C0_c_v + delta*C0_h_v
    
    dCA_c_v <- lambda_c*SA_c_v - gamma*CA_c_v + tau_c*C0_c_v - omega_c*CA_c_v - 
      sig_c_v$sigma_A*CA_c_v - alpha_v*CA_c_v + delta*CA_h_v
    
    dI_c_v <- sigma_c_v*C0_c_v + sig_c_v$sigma_A*CA_c_v - epsilon*I_c_v - alpha_I_v*I_c_v
    
    dS_II_c_v <- p*epsilon*I_c_v + gamma*C_II_c_v - lambda_c*S_II_c_v - phi*S_II_c_v - 
      alpha_v*S_II_c_v + delta*S_II_h_v + p*delta_I*I_h_v
    
    dC_II_c_v <- (1-p)*epsilon*I_c_v + lambda_c*S_II_c_v - gamma*C_II_c_v - 
      sig_c_v$sigma_II*C_II_c_v - alpha_v*C_II_c_v + delta*C_II_h_v + (1-p)*delta_I*I_h_v
    
    dI_II_c_v <- sig_c_v$sigma_II*C_II_c_v - epsilon*I_II_c_v - alpha_II_v*I_II_c_v
    
    dS_III_c_v <- p*epsilon*(I_II_c_v + I_III_c_v) + gamma*C_III_c_v - lambda_c*S_III_c_v - 
      phi*S_III_c_v - alpha_v*S_III_c_v + delta*S_III_h_v + p*delta_II*I_II_h_v + p*delta_III*I_III_h_v
    
    dC_III_c_v <- (1-p)*epsilon*(I_II_c_v + I_III_c_v) + lambda_c*S_III_c_v - gamma*C_III_c_v - 
      sig_c_v$sigma_III*C_III_c_v - alpha_v*C_III_c_v + delta*C_III_h_v + 
      (1-p)*delta_II*I_II_h_v + (1-p)*delta_III*I_III_h_v
    
    dI_III_c_v <- sig_c_v$sigma_III*C_III_c_v - epsilon*I_III_c_v - alpha_III_v*I_III_c_v
    
    # Return all derivatives
    list(c(
      # Non-vaccinated
      S0_h_nv=dS0_h_nv, SA_h_nv=dSA_h_nv, C0_h_nv=dC0_h_nv, CA_h_nv=dCA_h_nv, I_h_nv=dI_h_nv,
      S_II_h_nv=dS_II_h_nv, C_II_h_nv=dC_II_h_nv, I_II_h_nv=dI_II_h_nv,
      S_III_h_nv=dS_III_h_nv, C_III_h_nv=dC_III_h_nv, I_III_h_nv=dI_III_h_nv,
      S0_c_nv=dS0_c_nv, SA_c_nv=dSA_c_nv, C0_c_nv=dC0_c_nv, CA_c_nv=dCA_c_nv, I_c_nv=dI_c_nv,
      S_II_c_nv=dS_II_c_nv, C_II_c_nv=dC_II_c_nv, I_II_c_nv=dI_II_c_nv,
      S_III_c_nv=dS_III_c_nv, C_III_c_nv=dC_III_c_nv, I_III_c_nv=dI_III_c_nv,
      # Vaccinated
      S0_h_v=dS0_h_v, SA_h_v=dSA_h_v, C0_h_v=dC0_h_v, CA_h_v=dCA_h_v, I_h_v=dI_h_v,
      S_II_h_v=dS_II_h_v, C_II_h_v=dC_II_h_v, I_II_h_v=dI_II_h_v,
      S_III_h_v=dS_III_h_v, C_III_h_v=dC_III_h_v, I_III_h_v=dI_III_h_v,
      S0_c_v=dS0_c_v, SA_c_v=dSA_c_v, C0_c_v=dC0_c_v, CA_c_v=dCA_c_v, I_c_v=dI_c_v,
      S_II_c_v=dS_II_c_v, C_II_c_v=dC_II_c_v, I_II_c_v=dI_II_c_v,
      S_III_c_v=dS_III_c_v, C_III_c_v=dC_III_c_v, I_III_c_v=dI_III_c_v
    ))
  })
}


###############################################################################
# ---- CREATE INITIAL CONDITIONS FROM CALIBRATED EQUILIBRIUM ----
###############################################################################

create_vacc_initial_conditions <- function(calibrated_equilibrium, VC) {
  # calibrated_equilibrium = last row of ODE result from calibration
  # VC = vaccination coverage (proportion vaccinated, between 0 and 1)
  
  # Extract compartment names (without _h or _c suffix for base model)
  base_compartments <- c(
    "S0_h", "SA_h", "C0_h", "CA_h", "I_h",
    "S_II_h", "C_II_h", "I_II_h",
    "S_III_h", "C_III_h", "I_III_h",
    "S0_c", "SA_c", "C0_c", "CA_c", "I_c",
    "S_II_c", "C_II_c", "I_II_c",
    "S_III_c", "C_III_c", "I_III_c"
  )
  
  init_vacc <- c()
  
  # Split each compartment according to VC
  for (comp in base_compartments) {
    base_value <- as.numeric(calibrated_equilibrium[[comp]])
    
    # Non-vaccinated gets (1 - VC) of the population
    init_vacc[[paste0(comp, "_nv")]] <- base_value * (1 - VC)
    
    # Vaccinated gets VC of the population
    init_vacc[[paste0(comp, "_v")]] <- base_value * VC
  }
  
  return(init_vacc)
}

###############################################################################
# ---- RUN VACCINATION SCENARIOS ----
###############################################################################

run_vacc_scenarios <- function(params_calibrated, calibrated_equilibrium, 
                               time_vec, 
                               VE_values = c(0.3, 0.5, 0.7),
                               VC_values = c(0.2, 0.4, 0.6, 0.8)) {
  
  results <- list()
  
  # Test scenarios
  cat("\n=== Test Scenario 1: 100% non-vaccinated ===\n")
  init_test1 <- create_vacc_initial_conditions(calibrated_equilibrium, VC = 0)
  params_test1 <- params_calibrated
  params_test1["VE"] <- 0
  
  ode_test1 <- as.data.frame(
    deSolve::lsoda(y = init_test1, times = time_vec, func = cdiff_vacc, parms = params_test1)
  )
  
  results[["test1_VC0"]] <- list(
    params = params_test1,
    ode = ode_test1,
    VC = 0,
    VE = 0
  )
  
  cat("\n=== Test Scenario 2: 100% vaccinated with VE=0 ===\n")
  init_test2 <- create_vacc_initial_conditions(calibrated_equilibrium, VC = 1)
  params_test2 <- params_calibrated
  params_test2["VE"] <- 0
  
  ode_test2 <- as.data.frame(
    deSolve::lsoda(y = init_test2, times = time_vec, func = cdiff_vacc, parms = params_test2)
  )
  
  results[["test2_VC1_VE0"]] <- list(
    params = params_test2,
    ode = ode_test2,
    VC = 1,
    VE = 0
  )
  
  # Main scenarios: VE x VC combinations
  cat("\n=== Main Vaccination Scenarios ===\n")
  
  for (VE in VE_values) {
    for (VC in VC_values) {
      
      scenario_name <- sprintf("VE%.0f_VC%.0f", VE*100, VC*100)
      cat(sprintf("Running %s (VE=%.0f%%, VC=%.0f%%)...\n", scenario_name, VE*100, VC*100))
      
      # Create initial conditions
      init_scenario <- create_vacc_initial_conditions(calibrated_equilibrium, VC)
      
      # Set parameters
      params_scenario <- params_calibrated
      params_scenario["VE"] <- VE
      
      # Run ODE
      ode_scenario <- as.data.frame(
        deSolve::lsoda(y = init_scenario, times = time_vec, func = cdiff_vacc, parms = params_scenario)
      )
      
      # Store results
      results[[scenario_name]] <- list(
        params = params_scenario,
        ode = ode_scenario,
        VC = VC,
        VE = VE
      )
    }
  }
  
  return(results)
}

###############################################################################
# ---- COMPUTE VACCINATION METRICS ----
###############################################################################

compute_vacc_metrics <- function(ode_result, params_vec) {
  # Take last state (equilibrium or end of simulation)
  last <- as.list(ode_result[nrow(ode_result), ])
  
  # Total populations
  N_h_nv <- with(last, S0_h_nv + SA_h_nv + S_II_h_nv + S_III_h_nv + 
                   C0_h_nv + CA_h_nv + C_II_h_nv + C_III_h_nv + 
                   I_h_nv + I_II_h_nv + I_III_h_nv)
  
  N_c_nv <- with(last, S0_c_nv + SA_c_nv + S_II_c_nv + S_III_c_nv + 
                   C0_c_nv + CA_c_nv + C_II_c_nv + C_III_c_nv + 
                   I_c_nv + I_II_c_nv + I_III_c_nv)
  
  N_h_v <- with(last, S0_h_v + SA_h_v + S_II_h_v + S_III_h_v + 
                  C0_h_v + CA_h_v + C_II_h_v + C_III_h_v + 
                  I_h_v + I_II_h_v + I_III_h_v)
  
  N_c_v <- with(last, S0_c_v + SA_c_v + S_II_c_v + S_III_c_v + 
                  C0_c_v + CA_c_v + C_II_c_v + C_III_c_v + 
                  I_c_v + I_II_c_v + I_III_c_v)
  
  N_total <- N_h_nv + N_c_nv + N_h_v + N_c_v
  
  # Carriage prevalence
  C_h_nv <- with(last, C0_h_nv + CA_h_nv + C_II_h_nv + C_III_h_nv)
  C_c_nv <- with(last, C0_c_nv + CA_c_nv + C_II_c_nv + C_III_c_nv)
  C_h_v  <- with(last, C0_h_v + CA_h_v + C_II_h_v + C_III_h_v)
  C_c_v  <- with(last, C0_c_v + CA_c_v + C_II_c_v + C_III_c_v)
  
  portage_h_nv <- C_h_nv / N_h_nv
  portage_c_nv <- C_c_nv / N_c_nv
  portage_h_v  <- C_h_v / N_h_v
  portage_c_v  <- C_c_v / N_c_v
  portage_total <- (C_h_nv + C_c_nv + C_h_v + C_c_v) / N_total
  
  # CDI incidence (primo + recurrent)
  # Non-vaccinated
  sigma_h <- params_vec["sigma_h"]
  sigma_c <- params_vec["sigma_c"]
  k_A   <- params_vec["k_A"]
  k_II  <- params_vec["k_II"]
  k_III <- params_vec["k_III"]
  
  inc_primo_h_nv <- (sigma_h * last$C0_h_nv + k_A * sigma_h * last$CA_h_nv) / N_total * 365 * 1e5
  inc_rec_h_nv   <- (k_II * sigma_h * last$C_II_h_nv + k_III * sigma_h * last$C_III_h_nv) / N_total * 365 * 1e5
  
  inc_primo_c_nv <- (sigma_c * last$C0_c_nv + k_A * sigma_c * last$CA_c_nv) / N_total * 365 * 1e5
  inc_rec_c_nv   <- (k_II * sigma_c * last$C_II_c_nv + k_III * sigma_c * last$C_III_c_nv) / N_total * 365 * 1e5
  
  # Vaccinated (with VE)
  VE <- params_vec["VE"]
  sigma_h_v <- sigma_h * (1 - VE)
  sigma_c_v <- sigma_c * (1 - VE)
  
  inc_primo_h_v <- (sigma_h_v * last$C0_h_v + k_A * sigma_h_v * last$CA_h_v) / N_total * 365 * 1e5
  inc_rec_h_v   <- (k_II * sigma_h_v * last$C_II_h_v + k_III * sigma_h_v * last$C_III_h_v) / N_total * 365 * 1e5
  
  inc_primo_c_v <- (sigma_c_v * last$C0_c_v + k_A * sigma_c_v * last$CA_c_v) / N_total * 365 * 1e5
  inc_rec_c_v   <- (k_II * sigma_c_v * last$C_II_c_v + k_III * sigma_c_v * last$C_III_c_v) / N_total * 365 * 1e5
  
  # Total incidence (all populations combined)
  inc_primo_total <- inc_primo_h_nv + inc_primo_c_nv + inc_primo_h_v + inc_primo_c_v
  inc_rec_total   <- inc_rec_h_nv + inc_rec_c_nv + inc_rec_h_v + inc_rec_c_v
  inc_total       <- inc_primo_total + inc_rec_total
  
  list(
    portage_total = portage_total,
    portage_h_nv = portage_h_nv,
    portage_c_nv = portage_c_nv,
    portage_h_v = portage_h_v,
    portage_c_v = portage_c_v,
    inc_total = inc_total,
    inc_primo_total = inc_primo_total,
    inc_rec_total = inc_rec_total,
    N_total = N_total,
    N_h_nv = N_h_nv,
    N_c_nv = N_c_nv,
    N_h_v = N_h_v,
    N_c_v = N_c_v
  )
}



###############################################################################
# ---- PLOT VACCINATION SCENARIOS ----
###############################################################################

plot_vacc_impact <- function(vacc_results, baseline_metrics) {
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  
  # Extract metrics for all scenarios (excluding tests)
  scenario_names <- names(vacc_results)
  main_scenarios <- scenario_names[grepl("^VE", scenario_names)]
  
  metrics_list <- lapply(main_scenarios, function(sc_name) {
    sc <- vacc_results[[sc_name]]
    metrics <- compute_vacc_metrics(sc$ode, sc$params)
    
    data.frame(
      scenario = sc_name,
      VE = sc$VE,
      VC = sc$VC,
      inc_total = metrics$inc_total,
      inc_primo = metrics$inc_primo_total,
      inc_rec = metrics$inc_rec_total,
      portage = metrics$portage_total * 100
    )
  })
  
  df <- do.call(rbind, metrics_list)
  
  # Baseline values (no vaccination)
  base_inc_total <- baseline_metrics$inc_tot * 100000 * 365
  base_inc_primo <- baseline_metrics$inc_primo_tot * 100000 * 365
  base_inc_rec <- baseline_metrics$inc_rec_tot * 100000 * 365
  base_portage <- baseline_metrics$portage_tot * 100
  
  # Compute relative reductions
  df <- df %>%
    mutate(
      red_inc_total = (base_inc_total - inc_total) / base_inc_total * 100,
      red_inc_primo = (base_inc_primo - inc_primo) / base_inc_primo * 100,
      red_inc_rec = (base_inc_rec - inc_rec) / base_inc_rec * 100,
      red_portage = (base_portage - portage) / base_portage * 100,
      VE_label = factor(sprintf("VE = %.0f%%", VE*100)),
      VC_num = VC * 100
    )
  
  # Color palette
  pal_VE <- c(
    "VE = 30%" = "#C77DFF",
    "VE = 50%" = "#60A5FA",
    "VE = 70%" = "#6EE7B7"
  )
  
  # ---- PLOT 1: Total CDI incidence reduction ----
  p_inc_total <- ggplot(df, aes(x = VC_num, y = red_inc_total, color = VE_label, group = VE_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = pal_VE) +
    theme_bw() +
    labs(
      x = "Vaccination coverage (%)",
      y = "Reduction in total CDI incidence (%)",
      color = "Vaccine efficacy",
      title = "Impact on total CDI incidence"
    ) +
    theme(legend.position = "bottom")
  
  # ---- PLOT 2: Primo-CDI incidence reduction ----
  p_inc_primo <- ggplot(df, aes(x = VC_num, y = red_inc_primo, color = VE_label, group = VE_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = pal_VE) +
    theme_bw() +
    labs(
      x = "Vaccination coverage (%)",
      y = "Reduction in primo-CDI incidence (%)",
      color = "Vaccine efficacy",
      title = "Impact on primo-CDI incidence"
    ) +
    theme(legend.position = "bottom")
  
  # ---- PLOT 3: Recurrent CDI incidence reduction ----
  p_inc_rec <- ggplot(df, aes(x = VC_num, y = red_inc_rec, color = VE_label, group = VE_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = pal_VE) +
    theme_bw() +
    labs(
      x = "Vaccination coverage (%)",
      y = "Reduction in recurrent CDI incidence (%)",
      color = "Vaccine efficacy",
      title = "Impact on recurrent CDI incidence"
    ) +
    theme(legend.position = "bottom")
  
  # ---- PLOT 4: Asymptomatic carriage reduction ----
  p_portage <- ggplot(df, aes(x = VC_num, y = red_portage, color = VE_label, group = VE_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = pal_VE) +
    theme_bw() +
    labs(
      x = "Vaccination coverage (%)",
      y = "Reduction in asymptomatic carriage (%)",
      color = "Vaccine efficacy",
      title = "Impact on asymptomatic carriage"
    ) +
    theme(legend.position = "bottom")
  
  # ---- COMBINED PLOT ----
  p_combined <- (p_inc_total + p_inc_primo) / (p_inc_rec + p_portage) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(
      title = "Impact of vaccination on CDI burden",
      subtitle = "Relative reduction compared to no vaccination scenario"
    )
  
  list(
    inc_total = p_inc_total,
    inc_primo = p_inc_primo,
    inc_rec = p_inc_rec,
    portage = p_portage,
    combined = p_combined,
    data = df
  )
}

###############################################################################
# ---- VERIFY TEST SCENARIOS ----
###############################################################################

verify_test_scenarios <- function(vacc_results) {
  cat("\n=== VERIFICATION OF TEST SCENARIOS ===\n")
  
  test1 <- vacc_results$test1_VC0
  test2 <- vacc_results$test2_VC1_VE0
  
  metrics1 <- compute_vacc_metrics(test1$ode, test1$params)
  metrics2 <- compute_vacc_metrics(test2$ode, test2$params)
  
  cat("\nTest 1 (100% non-vaccinated):\n")
  cat(sprintf("  Total CDI incidence: %.2f per 100k/year\n", metrics1$inc_total))
  cat(sprintf("  Asymptomatic carriage: %.2f%%\n", metrics1$portage_total * 100))
  
  cat("\nTest 2 (100% vaccinated, VE=0):\n")
  cat(sprintf("  Total CDI incidence: %.2f per 100k/year\n", metrics2$inc_total))
  cat(sprintf("  Asymptomatic carriage: %.2f%%\n", metrics2$portage_total * 100))
  
  cat("\nDifferences (should be ~0):\n")
  cat(sprintf("  Incidence difference: %.6f\n", abs(metrics1$inc_total - metrics2$inc_total)))
  cat(sprintf("  Carriage difference: %.6f\n", abs(metrics1$portage_total - metrics2$portage_total)))
  
  # Check conservation of N_h and N_c
  cat("\n=== VERIFICATION OF POPULATION CONSERVATION ===\n")
  
  # For test1 (all non-vaccinated)
  last1 <- test1$ode[nrow(test1$ode), ]
  N_h_nv_test1 <- metrics1$N_h_nv
  N_c_nv_test1 <- metrics1$N_c_nv
  
  cat("\nTest 1 populations:\n")
  cat(sprintf("  N_h (non-vaccinated): %.2f\n", N_h_nv_test1))
  cat(sprintf("  N_c (non-vaccinated): %.2f\n", N_c_nv_test1))
  
  # Check a main scenario
  main_scenario_name <- names(vacc_results)[grepl("^VE", names(vacc_results))][1]
  main_scenario <- vacc_results[[main_scenario_name]]
  metrics_main <- compute_vacc_metrics(main_scenario$ode, main_scenario$params)
  
  VC_expected <- main_scenario$VC
  VC_observed_h <- metrics_main$N_h_v / (metrics_main$N_h_v + metrics_main$N_h_nv)
  VC_observed_c <- metrics_main$N_c_v / (metrics_main$N_c_v + metrics_main$N_c_nv)
  
  cat(sprintf("\n%s (VC=%.0f%%):\n", main_scenario_name, VC_expected*100))
  cat(sprintf("  Expected VC: %.2f%%\n", VC_expected*100))
  cat(sprintf("  Observed VC hospital: %.2f%%\n", VC_observed_h*100))
  cat(sprintf("  Observed VC community: %.2f%%\n", VC_observed_c*100))
  cat(sprintf("  Total N_h: %.2f (nv: %.2f, v: %.2f)\n", 
              metrics_main$N_h_nv + metrics_main$N_h_v, metrics_main$N_h_nv, metrics_main$N_h_v))
  cat(sprintf("  Total N_c: %.2f (nv: %.2f, v: %.2f)\n", 
              metrics_main$N_c_nv + metrics_main$N_c_v, metrics_main$N_c_nv, metrics_main$N_c_v))
  
  invisible(list(
    test1_metrics = metrics1,
    test2_metrics = metrics2,
    main_scenario_metrics = metrics_main
  ))
}

###############################################################################
# ---- CREATE SUMMARY TABLE ----
###############################################################################

create_vacc_summary_table <- function(vacc_results, baseline_metrics) {
  library(dplyr); library(gridExtra); library(grid)
  
  scenario_names <- names(vacc_results)
  main_scenarios <- scenario_names[grepl("^VE", scenario_names)]
  
  summary_list <- lapply(main_scenarios, function(sc_name) {
    sc <- vacc_results[[sc_name]]
    metrics <- compute_vacc_metrics(sc$ode, sc$params)
    
    data.frame(
      Scenario = sc_name,
      VE = sprintf("%.0f%%", sc$VE * 100),
      VC = sprintf("%.0f%%", sc$VC * 100),
      `Total CDI` = sprintf("%.2f", metrics$inc_total),
      `Primo CDI` = sprintf("%.2f", metrics$inc_primo_total),
      `Recurrent CDI` = sprintf("%.2f", metrics$inc_rec_total),
      `Carriage (%)` = sprintf("%.2f", metrics$portage_total * 100),
      check.names = FALSE
    )
  })
  
  summary_df <- do.call(rbind, summary_list)
  
  # Add baseline row
  baseline_row <- data.frame(
    Scenario = "Baseline (no vaccination)",
    VE = "-",
    VC = "0%",
    `Total CDI` = sprintf("%.2f", baseline_metrics$inc_tot * 100000 * 365),
    `Primo CDI` = sprintf("%.2f", baseline_metrics$inc_primo_tot * 100000 * 365),
    `Recurrent CDI` = sprintf("%.2f", baseline_metrics$inc_rec_tot * 100000 * 365),
    `Carriage (%)` = sprintf("%.2f", baseline_metrics$portage_tot * 100),
    check.names = FALSE
  )
  
  summary_df <- rbind(baseline_row, summary_df)
  
  # Create table
  tt <- ttheme_default(
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5, fontsize = 9),
      bg_params = list(fill = c("lightyellow", rep(c("white", "#F5F5F5"), length.out = nrow(summary_df)-1)))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", fontsize = 10),
      bg_params = list(fill = "steelblue", col = "white")
    )
  )
  
  tg <- tableGrob(summary_df, rows = NULL, theme = tt)
  
  title <- textGrob(
    "Vaccination Scenarios Summary",
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  note <- textGrob(
    "CDI incidence in cases per 100,000 population per year",
    gp = gpar(fontsize = 9, col = "gray30"),
    hjust = 0, x = 0.02
  )
  
  arrangeGrob(title, tg, note, heights = c(0.1, 0.8, 0.1), ncol = 1)
}

