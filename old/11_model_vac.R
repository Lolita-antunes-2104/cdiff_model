###############################################################################
###############################################################################

cdiff_model <- function(t, pop, params) {
  
  # Combine state variables (pop) and parameters (params) into a single list so that all compartments and parameters can be accessed by name directly
  with(as.list(c(pop, params)), {
    
  # ---- Totals (hospital/community and non-vaccinated/vaccinated) ----
    N_h_nv <- S0_h_nv + C0_h_nv + SA_h_nv + CA_h_nv + I_h_nv + S_II_h_nv + C_II_h_nv + I_II_h_nv + S_III_h_nv + C_III_h_nv + I_III_h_nv
    N_c_nv <- S0_c_nv + C0_c_nv + SA_c_nv + CA_c_nv + I_c_nv + S_II_c_nv + C_II_c_nv + I_II_c_nv + S_III_c_nv + C_III_c_nv + I_III_c_nv
    N_h_v <- S0_h_v + C0_h_v + SA_h_v + CA_h_v + I_h_v + S_II_h_v + C_II_h_v + I_II_h_v + S_III_h_v + C_III_h_v + I_III_h_v
    N_c_v <- S0_c_v + C0_c_v + SA_c_v + CA_c_v + I_c_v + S_II_c_v + C_II_c_v + I_II_c_v + S_III_c_v + C_III_c_v + I_III_c_v
    
    C_h_nv <- C0_h_nv + CA_h_nv + C_II_h_nv + C_III_h_nv 
    C_c_nv <- C0_c_nv + CA_c_nv + C_II_c_nv + C_III_c_nv 
    C_h_v <- C0_h_v + CA_h_v + C_II_h_v + C_III_h_v
    C_c_v <- C0_c_v + CA_c_v + C_II_c_v + C_III_c_v 
    
    I_tot_h_nv <- I_h_nv + I_II_h_nv + I_III_h_nv
    I_tot_c_nv <- I_c_nv + I_II_c_nv + I_III_c_nv
    I_tot_h_v <- I_h_v + I_II_h_v + I_III_h_v
    I_tot_c_v <- I_c_v + I_II_c_v + I_III_c_v
    
    N_h <- N_h_nv + N_h_v
    N_c <- N_c_nv + N_c_v 
    
    N_nv <- N_h_nv + N_c_nv
    N_v <- N_h_v + N_c_v
    
    C_h <- C_h_nv + C_h_v
    C_c <- C_c_nv + C_c_v 
    
    I_tot_h <- I_tot_h_nv + I_tot_h_v
    I_tot_c <- I_tot_c_nv + I_tot_c_v 
  
    # ---- Force of infection (shared across vaccinated and non-vaccinated, but separated for hospital and community) ----
    lambda_h <- (beta_h * ( C_h + nu_h * I_tot_h) / N_h )
    lambda_c <- (beta_c * ( C_c + nu_c * I_tot_c) / N_c )
    
    # Vaccinated et atb parameters
    sigma_h_v <- sigma_h * sigma_mult_v
    sigma_c_v <- sigma_c * sigma_mult_v
    
    gamma_h_v <- gamma_h * gamma_mult_v
    gamma_c_v <- gamma_c * gamma_mult_v
    
    tau_h_eff <- tau_h
    tau_c_eff <- tau_c
    
    
    # HOSPITAL / NON VACCINATED 
    # PRIMARY
    dS0_h_nv <- -lambda_h*S0_h_nv + gamma_h*C0_h_nv - tau_h*S0_h_nv + omega_h*SA_h_nv + phi_h*(S_II_h_nv + S_III_h_nv) + alpha*S0_c_nv - delta*S0_h_nv
    dSA_h_nv <- -lambda_h*SA_h_nv + gamma_h*CA_h_nv + tau_h*S0_h_nv - omega_h*SA_h_nv + alpha*SA_c_nv - delta*SA_h_nv
    dC0_h_nv <-  lambda_h*S0_h_nv - gamma_h*C0_h_nv - tau_h*C0_h_nv + omega_h*CA_h_nv - sigma_h*C0_h_nv + alpha*C0_c_nv - delta*C0_h_nv
    dCA_h_nv <-  lambda_h*SA_h_nv - gamma_h*CA_h_nv + tau_h*C0_h_nv - omega_h*CA_h_nv - k_A*sigma_h*CA_h_nv + alpha*CA_c_nv - delta*CA_h_nv
    dI_h_nv  <-  sigma_h*C0_h_nv + k_A*sigma_h*CA_h_nv - epsilon_h*I_h_nv + alpha_I*I_c_nv - delta_I*I_h_nv
    # FIRST RECURRENCE
    dS_II_h_nv <-  p_h*epsilon_h*I_h_nv - lambda_h*S_II_h_nv + gamma_h*C_II_h_nv - phi_h*S_II_h_nv + alpha*S_II_c_nv - delta*S_II_h_nv
    dC_II_h_nv <- (1-p_h)*epsilon_h*I_h_nv + lambda_h*S_II_h_nv - gamma_h*C_II_h_nv - k_II*sigma_h*C_II_h_nv + alpha*C_II_c_nv - delta*C_II_h_nv
    dI_II_h_nv <-  k_II*sigma_h*C_II_h_nv - epsilon_h*I_II_h_nv + alpha_II*I_II_c_nv - delta_II*I_II_h_nv
    # SECOND+ RECURRENCE
    dS_III_h_nv <-  p_h*epsilon_h*(I_II_h_nv + I_III_h_nv) - lambda_h*S_III_h_nv + gamma_h*C_III_h_nv - phi_h*S_III_h_nv + alpha*S_III_c_nv - delta*S_III_h_nv
    dC_III_h_nv <- (1-p_h)*epsilon_h*(I_II_h_nv + I_III_h_nv) + lambda_h*S_III_h_nv - gamma_h*C_III_h_nv - k_III*sigma_h*C_III_h_nv + alpha*C_III_c_nv - delta*C_III_h_nv
    dI_III_h_nv <-  k_III*sigma_h*C_III_h_nv - epsilon_h*I_III_h_nv + alpha_III*I_III_c_nv - delta_III*I_III_h_nv

    # COMMUNITY / NON VACCINATED 
    # Primary compartments
    dS0_c_nv <- -lambda_c*S0_c_nv + gamma_c*C0_c_nv - tau_c*S0_c_nv + omega_c*SA_c_nv + phi_c*(S_II_c_nv + S_III_c_nv) - alpha*S0_c_nv + delta*S0_h_nv
    dSA_c_nv <- -lambda_c*SA_c_nv + gamma_c*CA_c_nv + tau_c*S0_c_nv - omega_c*SA_c_nv - alpha*SA_c_nv + delta*SA_h_nv
    dC0_c_nv <-  lambda_c*S0_c_nv - gamma_c*C0_c_nv - tau_c*C0_c_nv + omega_c*CA_c_nv - sigma_c*C0_c_nv - alpha*C0_c_nv + delta*C0_h_nv
    dCA_c_nv <-  lambda_c*SA_c_nv - gamma_c*CA_c_nv + tau_c*C0_c_nv - omega_c*CA_c_nv - k_A*sigma_c*CA_c_nv - alpha*CA_c_nv + delta*CA_h_nv
    dI_c_nv  <-  sigma_c*C0_c_nv + k_A*sigma_c*CA_c_nv - epsilon_c*I_c_nv - alpha_I*I_c_nv
    
    dS_II_c_nv <-  p_c*epsilon_c*I_c_nv - lambda_c*S_II_c_nv + gamma_c*C_II_c_nv - phi_c*S_II_c_nv - alpha*S_II_c_nv + delta*S_II_h_nv + prop*delta_I*I_h_nv
    dC_II_c_nv <- (1-p_c)*epsilon_c*I_c_nv + lambda_c*S_II_c_nv - gamma_c*C_II_c_nv - k_II*sigma_c*C_II_c_nv - alpha*C_II_c_nv + delta*C_II_h_nv + (1 - prop)*delta_I*I_h_nv
    dI_II_c_nv <-  k_II*sigma_c*C_II_c_nv - epsilon_c*I_II_c_nv - alpha_II*I_II_c_nv
    
    dS_III_c_nv <-  p_c*epsilon_c*(I_II_c_nv + I_III_c_nv) - lambda_c*S_III_c_nv + gamma_c*C_III_c_nv - phi_c*S_III_c_nv - alpha*S_III_c_nv + delta*S_III_h_nv + prop*(delta_II*I_II_h_nv + delta_III*I_III_h_nv)
    dC_III_c_nv <- (1-p_c)*epsilon_c*(I_II_c_nv + I_III_c_nv) + lambda_c*S_III_c_nv - gamma_c*C_III_c_nv - k_III*sigma_c*C_III_c_nv + (1-prop)*(delta_II*I_II_h_nv + delta_III*I_III_h_nv)
    dI_III_c_nv <-  k_III*sigma_c*C_III_c_nv - epsilon_c*I_III_c_nv - alpha_III*I_III_c_nv
    
    # HOSPITAL / VACCINATED 
    # PRIMARY
    dS0_h_v <- -lambda_h*S0_h_v + gamma_h_v*C0_h_v - tau_h_eff*S0_h_v + omega_h*SA_h_v + phi_h*(S_II_h_v + S_III_h_v) + alpha*S0_c_v - delta*S0_h_v
    dSA_h_v <- -lambda_h*SA_h_v + gamma_h_v*CA_h_v + tau_h_eff*S0_h_v - omega_h*SA_h_v + alpha*SA_c_v - delta*SA_h_v
    dC0_h_v <-  lambda_h*S0_h_v - gamma_h_v*C0_h_v - tau_h_eff*C0_h_v + omega_h*CA_h_v - sigma_h_v*C0_h_v + alpha*C0_c_v - delta*C0_h_v
    dCA_h_v <-  lambda_h*SA_h_v - gamma_h_v*CA_h_v + tau_h_eff*C0_h_v - omega_h*CA_h_v - k_A*sigma_h_v*CA_h_v + alpha*CA_c_v - delta*CA_h_v
    dI_h_v  <-  sigma_h_v*C0_h_v + k_A*sigma_h_v*CA_h_v - epsilon_h*I_h_v + alpha_I*I_c_v - delta_I*I_h_v
    # FIRST RECURRENCE
    dS_II_h_v <-  p_h*epsilon_h*I_h_v - lambda_h*S_II_h_v + gamma_h_v*C_II_h_v - phi_h*S_II_h_v + alpha*S_II_c_v - delta*S_II_h_v
    dC_II_h_v <- (1-p_h)*epsilon_h*I_h_v + lambda_h*S_II_h_v - gamma_h_v*C_II_h_v - k_II*sigma_h_v*C_II_h_v + alpha*C_II_c_v - delta*C_II_h_v
    dI_II_h_v <-  k_II*sigma_h_v*C_II_h_v - epsilon_h*I_II_h_v + alpha_II*I_II_c_v - delta_II*I_II_h_v
    # SECOND+ RECURRENCE
    dS_III_h_v <-  p_h*epsilon_h*(I_II_h_v + I_III_h_v) - lambda_h*S_III_h_v + gamma_h_v*C_III_h_v - phi_h*S_III_h_v + alpha*S_III_c_v - delta*S_III_h_v
    dC_III_h_v <- (1-p_h)*epsilon_h*(I_II_h_v + I_III_h_v) + lambda_h*S_III_h_v - gamma_h_v*C_III_h_v - k_III*sigma_h_v*C_III_h_v + alpha*C_III_c_v - delta*C_III_h_v
    dI_III_h_v <-  k_III*sigma_h_v*C_III_h_v - epsilon_h*I_III_h_v + alpha_III*I_III_c_v - delta_III*I_III_h_v
    
    # COMMUNITY / VACCINATED 
    # Primary compartments
    dS0_c_v <- -lambda_c*S0_c_v + gamma_c_v*C0_c_v - tau_c_eff*S0_c_v + omega_c*SA_c_v + phi_c*(S_II_c_v + S_III_c_v) - alpha*S0_c_v + delta*S0_h_v
    dSA_c_v <- -lambda_c*SA_c_v + gamma_c_v*CA_c_v + tau_c_eff*S0_c_v - omega_c*SA_c_v - alpha*SA_c_v + delta*SA_h_v
    dC0_c_v <-  lambda_c*S0_c_v - gamma_c_v*C0_c_v - tau_c_eff*C0_c_v + omega_c*CA_c_v - sigma_c_v*C0_c_v - alpha*C0_c_v + delta*C0_h_v
    dCA_c_v <-  lambda_c*SA_c_v - gamma_c_v*CA_c_v + tau_c_eff*C0_c_v - omega_c*CA_c_v - k_A*sigma_c_v*CA_c_v - alpha*CA_c_v + delta*CA_h_v
    dI_c_v  <-  sigma_c_v*C0_c_v + k_A*sigma_c_v*CA_c_v - epsilon_c*I_c_v - alpha_I*I_c_v
    
    dS_II_c_v <-  p_c*epsilon_c*I_c_v - lambda_c*S_II_c_v + gamma_c_v*C_II_c_v - phi_c*S_II_c_v - alpha*S_II_c_v + delta*S_II_h_v + prop*delta_I*I_h_v
    dC_II_c_v <- (1-p_c)*epsilon_c*I_c_v + lambda_c*S_II_c_v - gamma_c_v*C_II_c_v - k_II*sigma_c_v*C_II_c_v - alpha*C_II_c_v + delta*C_II_h_v + (1 - prop)*delta_I*I_h_v
    dI_II_c_v <-  k_II*sigma_c_v*C_II_c_v - epsilon_c*I_II_c_v - alpha_II*I_II_c_v
    
    dS_III_c_v <-  p_c*epsilon_c*(I_II_c_v + I_III_c_v) - lambda_c*S_III_c_v + gamma_c_v*C_III_c_v - phi_c*S_III_c_v - alpha*S_III_c_v + delta*S_III_h_v + prop*(delta_II*I_II_h_v + delta_III*I_III_h_v)
    dC_III_c_v <- (1-p_c)*epsilon_c*(I_II_c_v + I_III_c_v) + lambda_c*S_III_c_v - gamma_c_v*C_III_c_v - k_III*sigma_c_v*C_III_c_v + (1-prop)*(delta_II*I_II_h_v + delta_III*I_III_h_v)
    dI_III_c_v <-  k_III*sigma_c_v*C_III_c_v - epsilon_c*I_III_c_v - alpha_III*I_III_c_v
    
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


###############################################################################
# autres
###############################################################################
compute_metrics_vacc <- function(out, params) {
  
  with(as.list(params), {
    
    sigma_h_v <- sigma_h * sigma_mult_v
    sigma_c_v <- sigma_c * sigma_mult_v
    
    # populations
    N_h_nv <- out[,"S0_h_nv"] + out[,"SA_h_nv"] + out[,"C0_h_nv"] + out[,"CA_h_nv"] +
      out[,"I_h_nv"] + out[,"S_II_h_nv"] + out[,"C_II_h_nv"] + out[,"I_II_h_nv"] +
      out[,"S_III_h_nv"] + out[,"C_III_h_nv"] + out[,"I_III_h_nv"]
    
    N_c_nv <- out[,"S0_c_nv"] + out[,"SA_c_nv"] + out[,"C0_c_nv"] + out[,"CA_c_nv"] +
      out[,"I_c_nv"] + out[,"S_II_c_nv"] + out[,"C_II_c_nv"] + out[,"I_II_c_nv"] +
      out[,"S_III_c_nv"] + out[,"C_III_c_nv"] + out[,"I_III_c_nv"]
    
    N_h_v <- out[,"S0_h_v"] + out[,"SA_h_v"] + out[,"C0_h_v"] + out[,"CA_h_v"] +
      out[,"I_h_v"] + out[,"S_II_h_v"] + out[,"C_II_h_v"] + out[,"I_II_h_v"] +
      out[,"S_III_h_v"] + out[,"C_III_h_v"] + out[,"I_III_h_v"]
    
    N_c_v <- out[,"S0_c_v"] + out[,"SA_c_v"] + out[,"C0_c_v"] + out[,"CA_c_v"] +
      out[,"I_c_v"] + out[,"S_II_c_v"] + out[,"C_II_c_v"] + out[,"I_II_c_v"] +
      out[,"S_III_c_v"] + out[,"C_III_c_v"] + out[,"I_III_c_v"]
    
    # incidence instantanée
    inc_h_nv <- sigma_h*out[,"C0_h_nv"] + k_A*sigma_h*out[,"CA_h_nv"] +
      k_II*sigma_h*out[,"C_II_h_nv"] + k_III*sigma_h*out[,"C_III_h_nv"]
    
    inc_c_nv <- sigma_c*out[,"C0_c_nv"] + k_A*sigma_c*out[,"CA_c_nv"] +
      k_II*sigma_c*out[,"C_II_c_nv"] + k_III*sigma_c*out[,"C_III_c_nv"]
    
    inc_h_v <- sigma_h_v*out[,"C0_h_v"] + k_A*sigma_h_v*out[,"CA_h_v"] +
      k_II*sigma_h_v*out[,"C_II_h_v"] + k_III*sigma_h_v*out[,"C_III_h_v"]
    
    inc_c_v <- sigma_c_v*out[,"C0_c_v"] + k_A*sigma_c_v*out[,"CA_c_v"] +
      k_II*sigma_c_v*out[,"C_II_c_v"] + k_III*sigma_c_v*out[,"C_III_c_v"]
    
    # incidence normalisée
    inc_h_nv_n <- inc_h_nv / (N_h_nv + N_c_nv)
    inc_c_nv_n <- inc_c_nv / (N_h_nv + N_c_nv)
    inc_h_v_n  <- inc_h_v  / (N_h_v + N_c_v)
    inc_c_v_n  <- inc_c_v  / (N_h_v + N_c_v)
    
    # prévalence de portage
    prev_h_nv <- (out[,"C0_h_nv"] + out[,"CA_h_nv"] + out[,"C_II_h_nv"] + out[,"C_III_h_nv"]) / N_h_nv
    prev_c_nv <- (out[,"C0_c_nv"] + out[,"CA_c_nv"] + out[,"C_II_c_nv"] + out[,"C_III_c_nv"]) / N_c_nv
    prev_h_v  <- (out[,"C0_h_v"]  + out[,"CA_h_v"]  + out[,"C_II_h_v"]  + out[,"C_III_h_v"])  / N_h_v
    prev_c_v  <- (out[,"C0_c_v"]  + out[,"CA_c_v"]  + out[,"C_II_c_v"]  + out[,"C_III_c_v"])  / N_c_v
    
    return(list(
      incidence = data.frame(
        time = out[,"time"],
        inc_h_nv, inc_c_nv, inc_h_v, inc_c_v,
        inc_h_nv_n, inc_c_nv_n, inc_h_v_n, inc_c_v_n
      ),
      carriage = data.frame(
        time = out[,"time"],
        prev_h_nv, prev_c_nv, prev_h_v, prev_c_v
      )
    ))
  })
}


# # ---- Instantaneous incidence ----
# # Non-vaccinated
# ####hospital
# inc_primo_h_nv <- sigma_h*C0_h_nv + k_A*sigma_h*CA_h_nv
# inc_rec_h_nv   <- k_II*sigma_h*C_II_h_nv + k_II*sigma_h*C_III_h_nv
# inc_h_nv       <- inc_primo_h_nv + inc_rec_h_nv
# ####community
# inc_primo_c_nv <- sigma_c*C0_c_nv + k_A*sigma_c*CA_c_nv
# inc_rec_c_nv   <- k_II*sigma_c*C_II_c_nv + k_II*sigma_c*C_III_c_nv
# inc_c_nv       <- inc_primo_c_nv + inc_rec_c_nv
# # Vaccinated
# ####hospital
# inc_primo_h_v <- sigma_h*C0_h_v + k_A*sigma_h*CA_h_v
# inc_rec_h_v   <- k_II*sigma_h*C_II_h_v + k_II*sigma_h*C_III_h_v
# inc_h_v       <- inc_primo_h_v + inc_rec_h_v
# ####community
# inc_primo_c_v <- sigma_c*C0_c_v + k_A*sigma_c*CA_c_v
# inc_rec_c_v   <- k_II*sigma_c*C_II_c_v + k_II*sigma_c*C_III_c_v
# inc_c_v       <- inc_primo_c_v + inc_rec_c_v
# 
# # ---- Normalised Instantaneous incidence ----
# # Non-vaccinated
# abs_inc_primo_h_nv <- inc_primo_h_nv / N_h_nv
# abs_inc_rec_h_nv   <- inc_rec_h_nv / N_h_nv
# abs_inc_h_nv       <- inc_h_nv / N_h_nv
# abs_inc_primo_c_nv <- inc_primo_c_nv / N_c_nv
# abs_inc_rec_c_nv   <- inc_rec_c_nv / N_c_nv
# abs_inc_c_nv       <- inc_c_nv / N_c_nv
# # Vaccinated
# abs_inc_primo_h_v <- inc_primo_h_v / N_h_v
# abs_inc_rec_h_v   <- inc_rec_h_v / N_h_v
# abs_inc_h_v       <- inc_h_v / N_h_v
# abs_inc_primo_c_v <- inc_primo_c_v / N_c_v
# abs_inc_rec_c_v   <- inc_rec_c_v / N_c_v
# abs_inc_c_v       <- inc_c_v / N_c_v
# 
# # ---- Prevalence of colonisée ----
# # Non-vaccinated
# prev_h_nv      <- (C0_h_nv + CA_h_nv + C_II_h_nv + C_III_h_nv) / N_h_nv
# prev_c_nv      <- (C0_c_nv + CA_c_nv + C_II_c_nv + C_III_c_nv) / N_c_nv
# # Vaccinated
# prev_h_v       <- (C0_h_v + CA_h_v + C_II_h_v + C_III_h_v) / N_h_v
# prev_c_v       <- (C0_c_v + CA_c_v + C_II_c_v + C_III_c_v) / N_c_v
# 
# 
