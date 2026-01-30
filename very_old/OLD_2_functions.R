###############################################################################
############################## 2 : FUNCTIONS ##################################
###############################################################################
#
# This file centralizes all functions used throughout the project.
# It is organized into two sections:
# 1. Model functions (transmission, progression, flows)
# 2. Metrics functions (epidemiological outputs)
#
###############################################################################
#
# R base functions used in this file include:
#
# - grepl()     : test whether a character string contains a pattern
# - names()     : retrieve names of vectors, lists, or data.frames
# - paste0()    : concatenate character strings without separator
# - sum()       : sum numeric values
# - if / else   : conditional logic
# - for()       : explicit loops
# - return()    : explicit return values from functions
# - match.arg() : argument validation for function inputs
#
###############################################################################

###############################################################################
# 1. MODEL FUNCTIONS
###############################################################################

# Compute total population counts for a given setting
compute_totals <- function(S0, C0, SA, CA, I, S_II, C_II, I_II, S_III, C_III, I_III) {
  list(S = S0 + SA + S_II + S_III,
       C = C0 + CA + C_II + C_III,
       I_tot = I + I_II + I_III,
       N = S0 + SA + S_II + S_III + C0 + CA + C_II + C_III + I + I_II + I_III)
}

# Compute force of infection (lambda)
compute_lambda <- function(beta, C, I_tot, N, nu) {
  if (N == 0) return(0)
  beta * (C + nu * I_tot) / N
}

# Compute stage-specific progression rates
compute_sigmas <- function(sigma_base, k_A, k_II, k_III) {
  list(sigma_A = k_A * sigma_base,
    sigma_II = k_II * sigma_base,
    sigma_III = k_III * sigma_base)
}

# Compute dynamic hospital admission rates (used during calibration)
compute_alpha <- function(alpha_mode, tot_h, tot_c, I_h, I_II_h, I_III_h, I_c, I_II_c, I_III_c, params) {
  if (alpha_mode == "dynamic") {
    out_hc <- params["delta"]*(tot_h$S + tot_h$C) + params["delta_I"]*I_h + params["delta_II"]*I_II_h + params["delta_III"]*I_III_h
    den_alpha <- params["w"]*(tot_c$S + tot_c$C) + params["w_I"]*I_c + params["w_II"]*I_II_c + params["w_III"]*I_III_c
    alpha <- if (den_alpha == 0) 0 
             else out_hc / den_alpha
    return(list(alpha = alpha, alpha_I = alpha*params["w_I"], alpha_II = alpha*params["w_II"], alpha_III = alpha * params["w_III"]))
  }
  if (alpha_mode == "fixed") {
    return(list(alpha = params["alpha"], alpha_I = params["alpha_I"], alpha_II = params["alpha_II"], alpha_III = params["alpha_III"]))
  }
  stop("alpha_mode must be either 'dynamic' or 'fixed' (only fixed for intervention scenarios)")
}

# Apply hospital–community flows to ODE derivatives
add_hc_flows <- function(dh, dc, compartments_h, compartments_c, 
                         alpha, alpha_I, alpha_II, alpha_III, delta, delta_I, delta_II, delta_III, p) {
  
  dh$dS0    <- dh$dS0    + alpha * compartments_c$S0        - delta * compartments_h$S0
  dh$dSA    <- dh$dSA    + alpha * compartments_c$SA        - delta * compartments_h$SA
  dh$dC0    <- dh$dC0    + alpha * compartments_c$C0        - delta * compartments_h$C0
  dh$dCA    <- dh$dCA    + alpha * compartments_c$CA        - delta * compartments_h$CA
  dh$dI     <- dh$dI     + alpha_I * compartments_c$I       - delta_I * compartments_h$I # Infected compartments (severity-specific)
  dh$dS_II  <- dh$dS_II  + alpha * compartments_c$S_II      - delta * compartments_h$S_II
  dh$dC_II  <- dh$dC_II  + alpha * compartments_c$C_II      - delta * compartments_h$C_II
  dh$dI_II  <- dh$dI_II  + alpha_II * compartments_c$I_II   - delta_II * compartments_h$I_II # Infected compartments (severity-specific)
  dh$dS_III <- dh$dS_III + alpha * compartments_c$S_III     - delta * compartments_h$S_III
  dh$dC_III <- dh$dC_III + alpha * compartments_c$C_III     - delta * compartments_h$C_III
  dh$dI_III <- dh$dI_III + alpha_III * compartments_c$I_III - delta_III * compartments_h$I_III # Infected compartments (severity-specific)
  
  dc$dS0    <- dc$dS0    - alpha * compartments_c$S0        + delta * compartments_h$S0
  dc$dSA    <- dc$dSA    - alpha * compartments_c$SA        + delta * compartments_h$SA
  dc$dC0    <- dc$dC0    - alpha * compartments_c$C0        + delta * compartments_h$C0
  dc$dCA    <- dc$dCA    - alpha * compartments_c$CA        + delta * compartments_h$CA
  dc$dI     <- dc$dI     - alpha_I   * compartments_c$I     + delta_I * compartments_h$I # Infected compartments (severity-specific)
  dc$dS_II  <- dc$dS_II                                     + p * delta_I * compartments_h$I # Discharge of infected patients to community recurrence states
  dc$dC_II  <- dc$dC_II                                     + (1 - p) * delta_I  * compartments_h$I # Discharge of infected patients to community recurrence states
  dc$dI_II  <- dc$dI_II  - alpha_II * compartments_c$I_II   + delta_II * compartments_h$I_II # Infected compartments (severity-specific)
  dc$dS_III <- dc$dS_III                                    + p * (delta_II * compartments_h$I_II + delta_III * compartments_h$I_III) # Discharge of infected patients to community recurrence states
  dc$dC_III <- dc$dC_III                                    + (1 - p) * (delta_II * compartments_h$I_II + delta_III * compartments_h$I_III) # Discharge of infected patients to community recurrence states
  dc$dI_III <- dc$dI_III - alpha_III * compartments_c$I_III + delta_III * compartments_h$I_III # Infected compartments (severity-specific)
  
  list(hospital = dh, community = dc)
}

# Create initial conditions before calibration (with vaccination model)
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
  CI_h <- 0
  CI_c <- 0
  # Return state vector
  c(
    # Hospital
    S0_h = S0_h, SA_h = SA_h, C0_h = C0_h, CA_h = CA_h, I_h = I_h,
    S_II_h = S_II_h, C_II_h = C_II_h, I_II_h = I_II_h,
    S_III_h = S_III_h, C_III_h = C_III_h, I_III_h = I_III_h,
    CI_h = CI_h,
    # Community
    S0_c = S0_c, SA_c = SA_c, C0_c = C0_c, CA_c = CA_c, I_c = I_c,
    S_II_c = S_II_c, C_II_c = C_II_c, I_II_c = I_II_c,
    S_III_c = S_III_c, C_III_c = C_III_c, I_III_c = I_III_c,
    CI_c = CI_c
  )
}

# Create initial conditions from calibrated equilibrium (without vaccination model)
create_initial_conditions_from_equilibrium <- function(last_state) {
  init <- last_state # copy the equilibrium state
  for (name in names(init)) { # loop over all compartment names
    if (grepl("^CI_", name)) {init[name] <- 0}} # reset cumultative incidence compartments 
  init # return initial conditions for the simulation 
}

# Create initial conditions for vaccination model from an equilibrium without vaccination
# With coverage vaccination !!!
create_initial_conditions_vaccination <- function(last_state_without_vacc, coverage = 0) {
  # Initialize empty state vector for the model with vaccination because vaccinated (_v) and non-vaccinated (_nv) compartments don't exist in the non-vaccinated equilibrium state
  init <- c() 
  # loop over all compartment names
  for (name in names(last_state_without_vacc)) { 
    
    # Skip cumulative incidence here because the vaccinated model is built from scratch (init <- c())
    # and CI_* compartments do not exist yet and will be added and reset explicitly at the end of the function
    if (grepl("^CI_", name)) next
    
    # Extract the value of the current compartment
    value <- last_state_without_vacc[name]
    
    # Hospital compartments 
    if (grepl("_h$", name)) {
      init[paste0(name, "_v")]  <- coverage * value # vaccinated individuals
      init[paste0(name, "_nv")] <- (1 - coverage) * value # non-vaccinated individuals
    }
    
    # Community compartment 
    if (grepl("_c$", name)) {
      init[paste0(name, "_v")]  <- coverage * value # vaccinated individuals
      init[paste0(name, "_nv")] <- (1 - coverage) * value # non-vaccinated individuals
    }
  }
  
  # Add cumulative incidence compartments and initialize them to zero.
  init["CI_h_v"]  <- 0
  init["CI_h_nv"] <- 0
  init["CI_c_v"]  <- 0
  init["CI_c_nv"] <- 0
  
  # Return the complete initial state for the vaccinated model
  init
}


###############################################################################
# 3. METRICS FUNCTIONS
###############################################################################
# Metrics are computed from the last state of the ODE (equilibrium or end
# of intervention). Vaccination stratification is handled transparently.
###############################################################################

###############################################################################
# GLOBAL DEFINITIONS
###############################################################################

# Colonization compartments by type
carriage_types <- list(
  C_primo = c("C0", "CA"),
  C_rec   = c("C_II", "C_III"),
  C_tot   = c("C0", "CA", "C_II", "C_III")
)

###############################################################################
# HELPERS
###############################################################################

# Check whether the model includes vaccination stratification
has_vaccination <- function(state) {
  return(any(grepl("_(v|nv)$", names(state))))
}

# Total population size
get_population <- function(state,
                           setting = c("h", "c", "total"),
                           vacc = c("v", "nv", "both")) {
  
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  vars <- names(state)
  
  # Setting
  if (setting == "h") vars <- vars[grepl("_h", vars)]
  if (setting == "c") vars <- vars[grepl("_c", vars)]
  
  # Vaccination
  if (has_vaccination(state)) {
    if (vacc == "v")  vars <- vars[grepl("_v$", vars)]
    if (vacc == "nv") vars <- vars[grepl("_nv$", vars)]
  } else {
    if (vacc != "both") return(NA_real_)
  }
  
  return(sum(state[vars]))
}

###############################################################################
# CARRIAGE PREVALENCE
###############################################################################

compute_carriage_prevalence <- function(last_state,
                                        type = c("C_primo", "C_rec", "C_tot"),
                                        setting = c("h", "c", "total"),
                                        vacc = c("v", "nv", "both")) {
  
  type <- match.arg(type)
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  settings_sel <- if (setting == "total") c("h", "c") else setting
  
  if (has_vaccination(last_state)) {
    vacc_sel <- if (vacc == "both") c("v", "nv") else vacc
  } else {
    if (vacc != "both") return(NA_real_)
    vacc_sel <- ""
  }
  
  # Numerator: carriage count
  C <- 0
  for (s in settings_sel) {
    for (v in vacc_sel) {
      if (v == "") {
        C <- C + sum(last_state[paste0(carriage_types[[type]], "_", s)])
      } else {
        C <- C + sum(last_state[paste0(carriage_types[[type]], "_", s, "_", v)])
      }
    }
  }
  
  # Denominator: population size
  N <- get_population(last_state, setting, vacc)
  
  if (N == 0) return(NA_real_)
  return(C / N)
}

###############################################################################
# CDI INCIDENCE – INSTANTANEOUS (FROM I COMPARTMENTS)
###############################################################################

get_CDI_incidence <- function(last_state,
                              type = c("primo", "recurrent", "total"),
                              setting = c("h", "c", "total"),
                              vacc = c("v", "nv", "both")) {
  
  type <- match.arg(type)
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  I_vars <- switch(type,
                   primo     = "I",
                   recurrent = c("I_II", "I_III"),
                   total     = c("I", "I_II", "I_III"))
  
  settings_sel <- if (setting == "total") c("h", "c") else setting
  
  if (has_vaccination(last_state)) {
    vacc_sel <- if (vacc == "both") c("v", "nv") else vacc
  } else {
    if (vacc != "both") return(NA_real_)
    vacc_sel <- ""
  }
  
  I <- 0
  for (s in settings_sel) {
    for (v in vacc_sel) {
      if (v == "") {
        I <- I + sum(last_state[paste0(I_vars, "_", s)])
      } else {
        I <- I + sum(last_state[paste0(I_vars, "_", s, "_", v)])
      }
    }
  }
  
  return(I)
}

compute_CDI_incidence <- function(last_state,
                                  type = c("primo", "recurrent", "total"),
                                  setting = c("h", "c", "total"),
                                  vacc = c("v", "nv", "both")) {
  
  type <- match.arg(type)
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  I <- get_CDI_incidence(last_state, type, setting, vacc)
  N <- get_population(last_state, setting, vacc)
  
  if (N == 0) return(NA_real_)
  return(I / N)
}

###############################################################################
# CDI INCIDENCE – CUMULATIVE (FROM CI_* IN ODE OUTPUT)
###############################################################################

compute_cumulative_CDI_incidence <- function(ode_df,
                                             setting = c("h", "c", "total"),
                                             vacc = c("v", "nv", "both"),
                                             t_start = 0,
                                             t_end = NULL) {
  
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  if (is.null(t_end)) t_end <- max(ode_df$time)
  if (!has_vaccination(ode_df)) vacc <- "both"
  
  CI_vars <- c()
  
  if (setting %in% c("h", "total")) {
    CI_vars <- c(CI_vars,
                 if (vacc == "both") grep("^CI_h", names(ode_df), value = TRUE)
                 else paste0("CI_h_", vacc))
  }
  
  if (setting %in% c("c", "total")) {
    CI_vars <- c(CI_vars,
                 if (vacc == "both") grep("^CI_c", names(ode_df), value = TRUE)
                 else paste0("CI_c_", vacc))
  }
  
  CI_end   <- sum(ode_df[ode_df$time == t_end, CI_vars])
  CI_start <- sum(ode_df[ode_df$time == t_start, CI_vars])
  
  return(CI_end - CI_start)
}

###############################################################################
# RECURRENCE PREVALENCE
###############################################################################

compute_recurrence_prevalence <- function(last_state,
                                          level = c("rec_1", "rec_2"),
                                          setting = c("h", "c", "total"),
                                          vacc = c("v", "nv", "both")) {
  
  level <- match.arg(level)
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  I     <- get_CDI_incidence(last_state, "primo", setting, vacc)
  I_II  <- get_CDI_incidence(last_state, "recurrent", setting, vacc)
  I_III <- get_CDI_incidence(last_state, "total", setting, vacc) - I - I_II
  
  if (level == "rec_1") {
    if (I == 0) return(NA_real_)
    return(I_II / I)
  }
  
  if (level == "rec_2") {
    if (I_II == 0) return(NA_real_)
    return(I_III / I_II)
  }
}

###############################################################################
# DYSBIOSIS
###############################################################################

compute_dysbiosis_count <- function(last_state,
                                    setting = c("h", "c", "total"),
                                    vacc = c("v", "nv", "both")) {
  
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  settings_sel <- if (setting == "total") c("h", "c") else setting
  
  if (has_vaccination(last_state)) {
    vacc_sel <- if (vacc == "both") c("v", "nv") else vacc
  } else {
    if (vacc != "both") return(NA_real_)
    vacc_sel <- ""
  }
  
  D <- 0
  for (s in settings_sel) {
    for (v in vacc_sel) {
      if (v == "") {
        D <- D + sum(last_state[paste0(c("SA", "CA"), "_", s)])
      } else {
        D <- D + sum(last_state[paste0(c("SA", "CA"), "_", s, "_", v)])
      }
    }
  }
  
  return(D)
}

compute_dysbiosis_prevalence <- function(last_state,
                                         setting = c("h", "c", "total"),
                                         vacc = c("v", "nv", "both")) {
  
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  D <- compute_dysbiosis_count(last_state, setting, vacc)
  N <- get_population(last_state, setting, vacc)
  
  if (N == 0) return(NA_real_)
  return(D / N)
}




