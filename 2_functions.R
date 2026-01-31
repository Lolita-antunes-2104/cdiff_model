###############################################################################
############################# 2 : FUNCTIONS ###################################
###############################################################################

###############################################################################
#  1. MODEL FUNCTIONS
###############################################################################

# ################### FOR CALIBRATION MODEL ################### 

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
  CumI_primo_h <- 0
  CumI_rec_h   <- 0
  CumI_h       <- 0
  CumI_primo_c <- 0
  CumI_rec_c   <- 0
  CumI_c       <- 0
  
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
# 2. METRICS FUNCTIONS
###############################################################################

# FOR CALIBRATION MODEL 

# prevalence
# prendre le C de compute metrics et faire C_h /N_h et C_c / N_c

# Compute instantaneous CDI incidence (primary, recurrent, or total)
# prendre 

# and normalised (absolue (/N_tot) / relative (_N_h or N_c))

# Compute cumulative CDI incidence from CumI_* compartments (primary, recurrent or total)

# and normalised(absolue (/N_tot) / relative (_N_h ou N_c))

# Compute recurrence prevalence from infected compartments
# rec_1 = I_II / I , rec_2 = I_III / I_II

# Faire fonction pour R0

# FOR SCENARIO MODEL



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




