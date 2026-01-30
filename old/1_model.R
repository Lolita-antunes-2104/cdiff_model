###############################################################################
############################# 1 : MODEL #######################################
###############################################################################

# This file defines the core ODE models used in the analysis.
# 1. cdiff_core() implements the within-setting CDI transmission and progression dynamics.
# 2. cdiff_hc_model() extends this core model to a hospital–community structure with optional antibiotic interventions.
# 3. cdiff_hc_vacc_model() further extends the model by stratifying the population by vaccination status and vaccine type.

###############################################################################
# 1. CORE SETTING DYNAMICS
###############################################################################

# This function defines the core CDI transmission and progression dynamics within a single setting (hospital or community).
# It includes infection, clearance, antibiotic-induced dysbiosis, and recurrence processes, but does NOT include any demographic exchanges.
# Those processes are added at a higher level.

cdiff_core <- function(S0, SA, C0, CA, I, S_II, C_II, I_II, S_III, C_III, I_III,
                       lambda, gamma, epsilon, p, phi, omega, tau, 
                       sigma, sigma_A, sigma_II, sigma_III) {
  
  # Primary compartments
  dS0 <- -lambda*S0 + gamma*C0 - tau*S0 + omega*SA + phi*(S_II + S_III) # Recovered recurrent cases return to S0
  dSA <- -lambda*SA + gamma*CA + tau*S0 - omega*SA
  dC0 <- lambda*S0 - gamma*C0 - tau*C0 + omega*CA - sigma*C0
  dCA <- lambda*SA - gamma*CA + tau*C0 - omega*CA - sigma_A*CA
  dI  <- sigma*C0 + sigma_A*CA - epsilon*I
  
  # First recurrence
  dS_II  <- p*epsilon*I - lambda*S_II + gamma*C_II - phi*S_II
  dC_II  <- (1-p)*epsilon*I + lambda*S_II - gamma*C_II - sigma_II*C_II
  dI_II  <- sigma_II*C_II - epsilon*I_II
  
  # Second+ recurrence
  dS_III <- p*epsilon*(I_II + I_III) - lambda*S_III + gamma*C_III - phi*S_III
  dC_III <- (1-p)*epsilon*(I_II + I_III) + lambda*S_III - gamma*C_III - sigma_III*C_III
  dI_III <- sigma_III*C_III - epsilon*I_III
  
  list(dS0 = dS0, dSA = dSA, dC0 = dC0, dCA = dCA, dI = dI,
    dS_II = dS_II, dC_II = dC_II, dI_II = dI_II,
    dS_III = dS_III, dC_III = dC_III, dI_III = dI_III)
}

###############################################################################
# 2. STRATIFIED MODEL HOSPITAL/COMMUNITY WITH ATB REDUCTION
###############################################################################

# This function defines the system of ODEs to be solved by lsoda()
# for a model stratified into hospital (h) and community (c) settings.
# The specific function signature (t, pop, params) is required by lsoda:
#  - t     : current time (provided automatically by the solver)
#  - pop   : current state vector (all compartments), updated at each time step
#  - params: fixed model parameters (not integrated over time)
# The function must return a vector of derivatives d(pop)/dt in the same order as pop.

cdiff_hc_model <- function(t, pop, params, alpha_mode, 
                           atb_reduction_h = 0, atb_reduction_c = 0, 
                           VE = 0, vacc_type = "both") { # ATB reduction scenario parameters
  
  # Combine state variables (pop) and parameters (params) into a single list so that all compartments and parameters can be accessed by name directly
  with(as.list(c(pop, params)), {
    
    # ---- Totals (hospital and community) ----
    tot_h <- compute_totals(S0_h, C0_h, SA_h, CA_h, I_h, S_II_h, C_II_h, I_II_h, S_III_h, C_III_h, I_III_h)
    tot_c <- compute_totals(S0_c, C0_c, SA_c, CA_c, I_c, S_II_c, C_II_c, I_II_c, S_III_c, C_III_c, I_III_c)
    
    # ---- ATB reduction effect ----
    tau_h_eff <- tau_h * (1 - atb_reduction_h)
    tau_c_eff <- tau_c * (1 - atb_reduction_c)
    
    # ---- Force of infection (transmission S -> C) ----
    lambda_h <- compute_lambda(beta_h, tot_h$C, tot_h$I_tot, tot_h$N, nu_h)
    lambda_c <- compute_lambda(beta_c, tot_c$C, tot_c$I_tot, tot_c$N, nu_c)
    
    # ---- Progression infection rates (C -> I) ----
    sig_h <- compute_sigmas(sigma_h, k_A, k_II, k_III)
    sig_c <- compute_sigmas(sigma_c, k_A, k_II, k_III)
    
    # ---- Instantaneous incidence (for cumulative incidence) ----
    # hospital
    inc_primo_h <- sigma_h * C0_h + sig_h$sigma_A * CA_h
    inc_rec_h   <- sig_h$sigma_II * C_II_h + sig_h$sigma_III * C_III_h
    inc_h <- inc_primo_h + inc_rec_h
    # community
    inc_primo_c <- sigma_c * C0_c + sig_c$sigma_A * CA_c
    inc_rec_c   <- sig_c$sigma_II * C_II_c + sig_c$sigma_III * C_III_c
    inc_c <- inc_primo_c + inc_rec_c
    # cumulative incidence counter (no outflow)
    dCumI_primo_h <- inc_primo_h
    dCumI_rec_h   <- inc_rec_h
    dCumI_h       <- inc_primo_h + inc_rec_h
    dCumI_primo_c <- inc_primo_c
    dCumI_rec_c   <- inc_rec_c
    dCumI_c       <- inc_primo_c + inc_rec_c
    
    # ---- Alpha (hospital admissions) ----
    # alpha_mode : dynamic : calibration / fixed : atb interventions scenarios
    alphas <- compute_alpha(alpha_mode, tot_h, tot_c, I_h, I_II_h, I_III_h, I_c, I_II_c, I_III_c, params) 
    alpha <- alphas$alpha
    alpha_I <- alphas$alpha_I
    alpha_II <- alphas$alpha_II
    alpha_III <- alphas$alpha_III 
    
    # ---- Core dynamics (no flows) ----
    hosp <-  cdiff_core(S0_h, SA_h, C0_h, CA_h, I_h, S_II_h, C_II_h, I_II_h, S_III_h, C_III_h, I_III_h,
                        lambda_h, gamma_h, epsilon_h, p_h, phi_h, omega_h, tau_h_eff, # effective atb exposure rate (atb intervention)
                        sigma_h, sig_h$sigma_A, sig_h$sigma_II, sig_h$sigma_III)
    comm <-  cdiff_core(S0_c, SA_c, C0_c, CA_c, I_c, S_II_c, C_II_c, I_II_c, S_III_c, C_III_c, I_III_c,
                        lambda_c, gamma_c, epsilon_c, p_c, phi_c, omega_c, tau_c_eff, # effective atb exposure rate (atb intervention)
                        sigma_c, sig_c$sigma_A, sig_c$sigma_II, sig_c$sigma_III)
    
    # ---- Add flows between hospital <-> community ----
    compartments_h <- list(S0 = S0_h, SA = SA_h, C0 = C0_h, CA = CA_h, I = I_h, S_II = S_II_h, C_II = C_II_h, I_II = I_II_h, S_III = S_III_h, C_III = C_III_h, I_III = I_III_h)
    compartments_c <- list(S0 = S0_c, SA = SA_c, C0 = C0_c, CA = CA_c, I = I_c, S_II = S_II_c, C_II = C_II_c, I_II = I_II_c, S_III = S_III_c, C_III = C_III_c, I_III = I_III_c)
    result <- add_hc_flows(hosp, comm, compartments_h, compartments_c, alpha, alpha_I, alpha_II, alpha_III, delta, delta_I, delta_II, delta_III, prop)
    hosp <- result$hospital
    comm <- result$community
    
    # ---- Return (order must match pop) ----
    list(c(S0_h = hosp$dS0, SA_h = hosp$dSA, C0_h = hosp$dC0, CA_h = hosp$dCA, I_h = hosp$dI,
           S_II_h = hosp$dS_II, C_II_h = hosp$dC_II, I_II_h = hosp$dI_II,
           S_III_h = hosp$dS_III, C_III_h = hosp$dC_III, I_III_h = hosp$dI_III, 
           CumI_primo_h = dCumI_primo_h, CumI_rec_h = dCumI_rec_h, CumI_h = dCumI_h, # cumulative CDI incidence (hospital)
           
           S0_c = comm$dS0, SA_c = comm$dSA, C0_c = comm$dC0, CA_c = comm$dCA, I_c = comm$dI,
           S_II_c = comm$dS_II, C_II_c = comm$dC_II, I_II_c = comm$dI_II,
           S_III_c = comm$dS_III, C_III_c = comm$dC_III, I_III_c = comm$dI_III, 
           CumI_primo_c = dCumI_primo_c, CumI_rec_c = dCumI_rec_c, CumI_c = dCumI_c # cumulative CDI incidence (community)
    ))
  })
}

###############################################################################
# 3. STRATIFIED MODEL HOSPITAL/COMMUNITY + VACCINATION WITH ATB REDUCTION
###############################################################################

# This function defines the system of ODEs to be solved by lsoda() (same structure as before)
# stratified by setting (hospital / community) and vaccination status.
# It extends the baseline model by incorporating vaccine effects on infection progression and/or colonization clearance, as well as antibiotic reduction scenarios.

cdiff_hc_vacc_model <- function(t, pop, params, alpha_mode = "fixed",
                                atb_reduction_h = 0, atb_reduction_c = 0, # ATB reduction scenario parameters
                                VE = 0, vacc_type = c("inf", "car", "both")) { # vaccination scenario parameters
  
  vacc_type <- match.arg(vacc_type)  # ensure vacc_type is a single valid option ("inf", "car", or "both")
  
  # Combine state variables (pop) and parameters (params) into a single list so that all compartments and parameters can be accessed by name directly
  with(as.list(c(pop, params)), {
    
    # ---- Totals (hospital/community and non-vaccinated/vaccinated) ----
    tot_h_nv <- compute_totals(S0_h_nv, C0_h_nv, SA_h_nv, CA_h_nv, I_h_nv, S_II_h_nv, C_II_h_nv, I_II_h_nv, S_III_h_nv, C_III_h_nv, I_III_h_nv)
    tot_c_nv <- compute_totals(S0_c_nv, C0_c_nv, SA_c_nv, CA_c_nv, I_c_nv, S_II_c_nv, C_II_c_nv, I_II_c_nv, S_III_c_nv, C_III_c_nv, I_III_c_nv)
    tot_h_v <- compute_totals(S0_h_v, C0_h_v, SA_h_v, CA_h_v, I_h_v, S_II_h_v, C_II_h_v, I_II_h_v, S_III_h_v, C_III_h_v, I_III_h_v)
    tot_c_v <- compute_totals(S0_c_v, C0_c_v, SA_c_v, CA_c_v, I_c_v, S_II_c_v, C_II_c_v, I_II_c_v, S_III_c_v, C_III_c_v, I_III_c_v)
    # ---- Aggregated totals by setting (used for lambda and alpha) ----
    tot_h <- compute_totals(S0 = S0_h_nv + S0_h_v, SA = SA_h_nv + SA_h_v, C0 = C0_h_nv + C0_h_v, CA = CA_h_nv + CA_h_v, I = I_h_nv + I_h_v, S_II = S_II_h_nv + S_II_h_v, C_II = C_II_h_nv + C_II_h_v, I_II = I_II_h_nv + I_II_h_v, S_III = S_III_h_nv + S_III_h_v, C_III = C_III_h_nv + C_III_h_v, I_III = I_III_h_nv + I_III_h_v)
    tot_c <- compute_totals(S0 = S0_c_nv + S0_c_v, SA = SA_c_nv + SA_c_v, C0 = C0_c_nv + C0_c_v, CA = CA_c_nv + CA_c_v, I = I_c_nv + I_c_v, S_II = S_II_c_nv + S_II_c_v, C_II = C_II_c_nv + C_II_c_v, I_II = I_II_c_nv + I_II_c_v, S_III = S_III_c_nv + S_III_c_v, C_III = C_III_c_nv + C_III_c_v, I_III = I_III_c_nv + I_III_c_v)
    tot_nv <- compute_totals(S0 = S0_h_nv + S0_c_nv, SA = SA_h_nv + SA_c_nv, C0 = C0_h_nv + C0_c_nv, CA = CA_h_nv + CA_c_nv, I = I_h_nv + I_c_nv, S_II = S_II_h_nv + S_II_c_nv, C_II = C_II_h_nv + C_II_c_nv, I_II = I_II_h_nv + I_II_c_nv, S_III = S_III_h_nv + S_III_c_nv, C_III = C_III_h_nv + C_III_c_nv, I_III = I_III_h_nv + I_III_c_nv)
    tot_v <- compute_totals(S0 = S0_h_v + S0_c_v, SA = SA_h_v + SA_c_v, C0 = C0_h_v + C0_c_v, CA = CA_h_v + CA_c_v, I = I_h_v + I_c_v, S_II = S_II_h_v + S_II_c_v, C_II = C_II_h_v + C_II_c_v, I_II = I_II_h_v + I_II_c_v, S_III = S_III_h_v + S_III_c_v, C_III = C_III_h_v + C_III_c_v, I_III = I_III_h_v + I_III_c_v)
    
    # ---- ATB reduction effect ----
    tau_h_eff <- tau_h * (1 - atb_reduction_h)
    tau_c_eff <- tau_c * (1 - atb_reduction_c)
    
    # ---- Vaccination effects ----
    # Reduction in infection progression (sigma) and/or increase in clearance (gamma)
    sigma_mult_v <- if (vacc_type %in% c("inf", "both")) (1 - VE) else 1
    gamma_mult_v <- if (vacc_type %in% c("car", "both")) (1 + VE) else 1
    
    # ---- Force of infection (shared across vaccinated and non-vaccinated, but separated for hospital and community) ----
    lambda_h <- compute_lambda(beta_h, tot_h$C, tot_h$I_tot, tot_h$N, nu_h)
    lambda_c <- compute_lambda(beta_c, tot_c$C, tot_c$I_tot, tot_c$N, nu_c)
    
    # ---- Sigmas with VE (hospital/community and non-vaccinated/vaccinated) ----
    # Vaccine effect on infection progression applies to primary and recurrent infections
    sig_h_nv <- compute_sigmas(sigma_h, k_A, k_II, k_III)
    sig_c_nv <- compute_sigmas(sigma_c, k_A, k_II, k_III) 
    sig_h_v <- compute_sigmas(sigma_h * sigma_mult_v, k_A, k_II, k_III)
    sig_c_v <- compute_sigmas(sigma_c * sigma_mult_v, k_A, k_II, k_III)
    
    # ---- Instantaneous incidence (for cumulative incidence) ----
    # Non-vaccinated
    inc_primo_h_nv <- sigma_h * C0_h_nv + sig_h_nv$sigma_A * CA_h_nv
    inc_rec_h_nv   <- sig_h_nv$sigma_II * C_II_h_nv + sig_h_nv$sigma_III * C_III_h_nv
    inc_h_nv       <- inc_primo_h_nv + inc_rec_h_nv
    inc_primo_c_nv <- sigma_c * C0_c_nv + sig_c_nv$sigma_A * CA_c_nv
    inc_rec_c_nv   <- sig_c_nv$sigma_II * C_II_c_nv + sig_c_nv$sigma_III * C_III_c_nv
    inc_c_nv       <- inc_primo_c_nv + inc_rec_c_nv
    # Vaccinated
    inc_primo_h_v <- (sigma_h * sigma_mult_v) * C0_h_v + sig_h_v$sigma_A * CA_h_v
    inc_rec_h_v   <- sig_h_v$sigma_II * C_II_h_v + sig_h_v$sigma_III * C_III_h_v
    inc_h_v       <- inc_primo_h_v + inc_rec_h_v
    inc_primo_c_v <- (sigma_c * sigma_mult_v) * C0_c_v + sig_c_v$sigma_A * CA_c_v
    inc_rec_c_v   <- sig_c_v$sigma_II * C_II_c_v + sig_c_v$sigma_III * C_III_c_v
    inc_c_v       <- inc_primo_c_v + inc_rec_c_v
    # Cumulative incidence counter (no outflow)
    # hospital 
    dCumI_primo_h_nv <- inc_primo_h_nv
    dCumI_rec_h_nv   <- inc_rec_h_nv
    dCumI_h_nv       <- inc_primo_h_nv + inc_rec_h_nv
    dCumI_primo_h_v  <- inc_primo_h_v
    dCumI_rec_h_v    <- inc_rec_h_v
    dCumI_h_v        <- inc_primo_h_v + inc_rec_h_v
    # community 
    dCumI_primo_c_nv <- inc_primo_c_nv
    dCumI_rec_c_nv   <- inc_rec_c_nv
    dCumI_c_nv       <- inc_primo_c_nv + inc_rec_c_nv
    dCumI_primo_c_v  <- inc_primo_c_v
    dCumI_rec_c_v    <- inc_rec_c_v
    dCumI_c_v        <- inc_primo_c_v + inc_rec_c_v
    
    # ---- Alpha (ONLY fixed) ----
    # alpha_mode is set to "fixed" by default, as vaccination scenarios are evaluated using hospital admission rates calibrated beforehand
    alphas <- compute_alpha(alpha_mode, tot_h, tot_c, I_h = I_h_nv + I_h_v, I_II_h = I_II_h_nv + I_II_h_v, I_III_h = I_III_h_nv + I_III_h_v,
                            I_c = I_c_nv + I_c_v, I_II_c = I_II_c_nv + I_II_c_v, I_III_c = I_III_c_nv + I_III_c_v, params)
    alpha <- alphas$alpha
    alpha_I <- alphas$alpha_I
    alpha_II <- alphas$alpha_II
    alpha_III <- alphas$alpha_III
    
    # ---- Core dynamics (no hospital–community flows) ---- 
    # Non-vaccinated individuals: baseline progression rates, with ATB-modified tau
    h_nv <-  cdiff_core(S0_h_nv, SA_h_nv, C0_h_nv, CA_h_nv, I_h_nv, S_II_h_nv, C_II_h_nv, I_II_h_nv, S_III_h_nv, C_III_h_nv, I_III_h_nv,
                        lambda_h, gamma_h, # baseline clearance rate
                        epsilon_h, p_h, phi_h, omega_h, tau_h_eff, # effective ATB exposure rate (ATB intervention)
                        sigma_h, sig_h_nv$sigma_A, sig_h_nv$sigma_II, sig_h_nv$sigma_III) # baseline progression rate
    c_nv <-  cdiff_core(S0_c_nv, SA_c_nv, C0_c_nv, CA_c_nv, I_c_nv, S_II_c_nv, C_II_c_nv, I_II_c_nv, S_III_c_nv, C_III_c_nv, I_III_c_nv,
                        lambda_c, gamma_c, # baseline clearance rate
                        epsilon_c, p_c, phi_c, omega_c, tau_c_eff, # effective ATB exposure rate (ATB intervention)
                        sigma_c, sig_c_nv$sigma_A, sig_c_nv$sigma_II, sig_c_nv$sigma_III) # baseline progression rate 
    # Vaccinated individuals: modified progression and/or clearance rates depending on vaccine type
    h_v <-  cdiff_core(S0_h_v, SA_h_v, C0_h_v, CA_h_v, I_h_v, S_II_h_v, C_II_h_v, I_II_h_v, S_III_h_v, C_III_h_v, I_III_h_v,
                       lambda_h, gamma_h * gamma_mult_v, # increased clearance if vaccine reduces carriage
                       epsilon_h, p_h, phi_h, omega_h, tau_h_eff, # effective ATB exposure rate (same as non-vaccinated)
                       sigma_h * sigma_mult_v, sig_h_v$sigma_A, sig_h_v$sigma_II, sig_h_v$sigma_III) # reduced progression to infection if vaccine prevents infection
    c_v <-  cdiff_core(S0_c_v, SA_c_v, C0_c_v, CA_c_v, I_c_v, S_II_c_v, C_II_c_v, I_II_c_v, S_III_c_v, C_III_c_v, I_III_c_v,
                       lambda_c, gamma_c * gamma_mult_v,  # increased clearance if vaccine reduces carriage
                       epsilon_c, p_c, phi_c, omega_c, tau_c_eff, # effective ATB exposure rate (same as non-vaccinated)
                       sigma_c * sigma_mult_v, sig_c_v$sigma_A, sig_c_v$sigma_II, sig_c_v$sigma_III) # reduced progression to infection if vaccine prevents infection
    
    # ---- Add flows between hospital <-> community ----
    # non-vaccinated 
    comp_h_nv <- list(S0 = S0_h_nv, SA = SA_h_nv, C0 = C0_h_nv, CA = CA_h_nv, S_II = S_II_h_nv, C_II = C_II_h_nv, S_III = S_III_h_nv, C_III = C_III_h_nv, I = I_h_nv, I_II = I_II_h_nv, I_III = I_III_h_nv)
    comp_c_nv <- list(S0 = S0_c_nv, SA = SA_c_nv, C0 = C0_c_nv, CA = CA_c_nv, S_II = S_II_c_nv, C_II = C_II_c_nv, S_III = S_III_c_nv, C_III = C_III_c_nv, I = I_c_nv, I_II = I_II_c_nv, I_III = I_III_c_nv)
    result_nv <- add_hc_flows(h_nv, c_nv, comp_h_nv, comp_c_nv, alpha, alpha_I, alpha_II, alpha_III, delta, delta_I, delta_II, delta_III, prop)
    h_nv <- result_nv$hospital
    c_nv <- result_nv$community
    # vaccinated
    comp_h_v <- list(S0 = S0_h_v, SA = SA_h_v, C0 = C0_h_v, CA = CA_h_v, S_II = S_II_h_v, C_II = C_II_h_v, S_III = S_III_h_v, C_III = C_III_h_v, I = I_h_v, I_II = I_II_h_v, I_III = I_III_h_v)
    comp_c_v <- list(S0 = S0_c_v, SA = SA_c_v, C0 = C0_c_v, CA = CA_c_v, S_II = S_II_c_v, C_II = C_II_c_v, S_III = S_III_c_v, C_III = C_III_c_v, I = I_c_v, I_II = I_II_c_v, I_III = I_III_c_v)
    result_v <- add_hc_flows(h_v, c_v, comp_h_v, comp_c_v, alpha, alpha_I, alpha_II, alpha_III, delta, delta_I, delta_II, delta_III, prop)
    h_v <- result_v$hospital
    c_v <- result_v$community
    
    # ---- Return ----
    list(c(
      # Hospital non-vaccinated
      S0_h_nv = h_nv$dS0, SA_h_nv = h_nv$dSA, C0_h_nv = h_nv$dC0, CA_h_nv = h_nv$dCA, I_h_nv = h_nv$dI,
      S_II_h_nv = h_nv$dS_II, C_II_h_nv = h_nv$dC_II, I_II_h_nv = h_nv$dI_II,
      S_III_h_nv = h_nv$dS_III, C_III_h_nv = h_nv$dC_III, I_III_h_nv = h_nv$dI_III,
      CumI_primo_h_nv = dCumI_primo_h_nv, CumI_rec_h_nv = dCumI_rec_h_nv, CumI_h_nv = dCumI_h_nv, # cumulative CDI incidence (non-vaccinated / hospital)
      # Community non-vaccinated
      S0_c_nv = c_nv$dS0, SA_c_nv = c_nv$dSA, C0_c_nv = c_nv$dC0, CA_c_nv = c_nv$dCA, I_c_nv = c_nv$dI,
      S_II_c_nv = c_nv$dS_II, C_II_c_nv = c_nv$dC_II, I_II_c_nv = c_nv$dI_II,
      S_III_c_nv = c_nv$dS_III, C_III_c_nv = c_nv$dC_III, I_III_c_nv = c_nv$dI_III,
      CumI_primo_c_nv = dCumI_primo_c_nv, CumI_rec_c_nv = dCumI_rec_c_nv, CumI_c_nv = dCumI_c_nv, # cumulative CDI incidence (non-vaccinated / community)
      # Hospital vaccinated
      S0_h_v = h_v$dS0, SA_h_v = h_v$dSA, C0_h_v = h_v$dC0, CA_h_v = h_v$dCA, I_h_v = h_v$dI,
      S_II_h_v = h_v$dS_II, C_II_h_v = h_v$dC_II, I_II_h_v = h_v$dI_II,
      S_III_h_v = h_v$dS_III, C_III_h_v = h_v$dC_III, I_III_h_v = h_v$dI_III,
      CumI_primo_h_v = dCumI_primo_h_v, CumI_rec_h_v = dCumI_rec_h_v, CumI_h_v = dCumI_h_v, # cumulative CDI incidence (vaccinated / hospital)
      # Community vaccinated
      S0_c_v = c_v$dS0, SA_c_v = c_v$dSA, C0_c_v = c_v$dC0, CA_c_v = c_v$dCA, I_c_v = c_v$dI,
      S_II_c_v = c_v$dS_II, C_II_c_v = c_v$dC_II, I_II_c_v = c_v$dI_II,
      S_III_c_v = c_v$dS_III, C_III_c_v = c_v$dC_III, I_III_c_v = c_v$dI_III,
      CumI_primo_c_v = dCumI_primo_c_v, CumI_rec_c_v = dCumI_rec_c_v, CumI_c_v = dCumI_c_v # cumulative CDI incidence (vaccinated / community)
    ))
  })
}




