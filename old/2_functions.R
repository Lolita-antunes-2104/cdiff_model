###############################################################################
####################### 2 : MODEL and METRICS FUNCTIONS #######################
###############################################################################

# This file centralizes all helper functions used throughout the project.
# It is organized into two sections:
# 1. Model helper functions: totals, force of infection, progression rates, hospital–community flows, and initial conditions (baseline and vaccination stratification).
# 2. Metric functions: population totals, carriage, incidence (instantaneous and cumulative),
#    recurrence prevalence, and dysbiosis burden, all supporting hospital/community/total and vaccinated/non-vaccinated/both outputs.

###############################################################################
# R base functions used in this file include:
# - grepl()     : test whether a character string contains a pattern
# - names()     : retrieve names of vectors, lists, or data.frames
# - paste0()    : concatenate character strings without separator
# - sum()       : sum numeric values
# - if / else   : conditional logic
# - for()       : explicit loops
# - return()    : explicit return values from functions
# - match.arg() : argument validation for function inputs
# - stop()      : stop execution with an error message (input validation / safety checks)
###############################################################################

###############################################################################
# 1. MODEL FUNCTIONS
###############################################################################

# Compute total population counts for a given setting
# Returns: list with S (susceptibles), C (colonized), I_tot (all infected), N (total population)
compute_totals <- function(S0, C0, SA, CA, I, S_II, C_II, I_II, S_III, C_III, I_III) {
  return(list(S = S0 + SA + S_II + S_III,
              C = C0 + CA + C_II + C_III,
              I_tot = I + I_II + I_III,
              N = S0 + SA + S_II + S_III + C0 + CA + C_II + C_III + I + I_II + I_III))
}

# Compute force of infection (transmission rate from colonized/infected to susceptibles)
# Returns: lambda (per capita transmission rate)
compute_lambda <- function(beta, C, I_tot, N, nu) {
  if (N == 0) return(0)
  return(beta * (C + nu * I_tot) / N)
}

# Compute stage-specific progression rates from colonization to infection
# Returns: list with sigma_A (dysbiotic), sigma_II (first recurrence), sigma_III (second+ recurrence)
compute_sigmas <- function(sigma_base, k_A, k_II, k_III) {
  return(list(sigma_A   = k_A  * sigma_base, 
              sigma_II  = k_II * sigma_base,
              sigma_III = k_III * sigma_base))
}

# Compute hospital admission rates (alpha) based on mode: dynamic for calibration, fixed for intervention scenarios
# Returns: list with alpha (base rate) and alpha_I, alpha_II, alpha_III (infection severity-specific rates)
# ATTENTION : 
# - tot_h and tot_c must be lists returned by compute_totals() (contain S, C, I_tot, N)
# - w, w_I, w_II, w_III are relative admission weights by infection state
# Dynamic mode is used for calibration, fixed mode for intervention scenarios
compute_alpha <- function(alpha_mode, tot_h, tot_c, I_h, I_II_h, I_III_h, I_c, I_II_c, I_III_c, params) {
  
  if (alpha_mode == "dynamic") {
    
    out_hc <- params["delta"] * (tot_h$S + tot_h$C)  + params["delta_I"] * I_h  + params["delta_II"] * I_II_h  + params["delta_III"] * I_III_h
    
    den_alpha <- params["w"] * (tot_c$S + tot_c$C)   + params["w_I"] * I_c      + params["w_II"] * I_II_c      + params["w_III"] * I_III_c
    
    alpha <- if (den_alpha == 0) 0 else out_hc / den_alpha
    
    return(list(alpha = alpha, alpha_I = alpha * params["w_I"], alpha_II = alpha * params["w_II"], alpha_III = alpha * params["w_III"]))
  }
  
  if (alpha_mode == "fixed") {
    return(list(alpha = params["alpha"], alpha_I = params["alpha_I"], alpha_II = params["alpha_II"],alpha_III  = params["alpha_III"]))
  }
  
  stop("alpha_mode must be either 'dynamic' or 'fixed'")
}

# Apply bidirectional flows between hospital and community to ODE derivatives
# Modifies hospital (dh) and community (dc) derivatives in place, accounting for admissions (alpha) and discharges (delta)
# Special handling: infected patients discharged from hospital transition to recurrence states in community
# Returns: list with updated hospital and community derivative lists
add_hc_flows <- function(dh, dc, compartments_h, compartments_c, 
                         alpha, alpha_I, alpha_II, alpha_III, delta, delta_I, delta_II, delta_III, prop) {
  
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
  dc$dI     <- dc$dI     - alpha_I * compartments_c$I # Infected compartments (severity-specific)
  dc$dS_II  <- dc$dS_II  - alpha * compartments_c$S_II      + delta * compartments_h$S_II     + prop * delta_I * compartments_h$I # Discharge of infected patients to community recurrence states
  dc$dC_II  <- dc$dC_II  - alpha * compartments_c$C_II      + delta * compartments_h$C_II     + (1 - prop) * delta_I * compartments_h$I # Discharge of infected patients to community recurrence states
  dc$dI_II  <- dc$dI_II  - alpha_II * compartments_c$I_II # Infected compartments (severity-specific)
  dc$dS_III <- dc$dS_III - alpha * compartments_c$S_III     + delta * compartments_h$S_III    + prop * (delta_II * compartments_h$I_II + delta_III * compartments_h$I_III) # Discharge of infected patients to community recurrence states
  dc$dC_III <- dc$dC_III - alpha * compartments_c$C_III     + delta * compartments_h$C_III    + (1 - prop) * (delta_II * compartments_h$I_II + delta_III * compartments_h$I_III) # Discharge of infected patients to community recurrence states
  dc$dI_III <- dc$dI_III - alpha_III * compartments_c$I_III # Infected compartments (severity-specific)
  
  return(list(hospital = dh, community = dc))
}

# Create initial conditions for model calibration with user-specified prevalence levels
# Distributes initial infections across primary and recurrent compartments in both hospital and community
# Returns: named vector of initial compartment values for all state variables
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

# Create initial conditions from a calibrated equilibrium state for intervention scenarios
# Resets cumulative incidence (CumI_*) compartments to zero while preserving all other compartment values
# Returns: named vector of initial conditions ready for intervention scenario simulations
create_initial_conditions_from_equilibrium <- function(last_state) {
  init <- last_state # copy the equilibrium state
  for (name in names(init)) { # loop over all compartment names
    if (grepl("^CumI_", name)) {init[name] <- 0}} # reset cumulative incidence compartments 
  return(init) # return initial conditions for the simulation 
}

# Create initial conditions for vaccination model by stratifying an equilibrium state without vaccination
# Splits each compartment into vaccinated (_v) and non-vaccinated (_nv) based on coverage proportion
# Resets all cumulative incidence compartments to zero for the new simulation
# Returns: named vector with doubled compartments (_h_v, _h_nv, _c_v, _c_nv) ready for vaccination scenarios
create_initial_conditions_vaccination <- function(last_state_without_vacc, coverage = 0) {
  
  if (coverage < 0 || coverage > 1) stop("coverage must be between 0 and 1")
  
  # Initialize empty state vector for the model with vaccination because vaccinated (_v) and non-vaccinated (_nv) compartments don't exist in the non-vaccinated equilibrium state
  init <- c() 
  
  # loop over all compartment names
  for (name in names(last_state_without_vacc)) { 
    
    # Skip cumulative incidence here because the vaccinated model is built from scratch (init <- c())
    # and CumI_* compartments do not exist yet and will be added and reset explicitly at the end of the function
    if (grepl("^CumI_", name)) next
    
    # Extract the value of the current compartment
    value <- last_state_without_vacc[name]
    
    # Hospital compartments 
    if (grepl("_h$", name)) {
      init[paste0(name, "_v")]  <- coverage * value # vaccinated individuals
      init[paste0(name, "_nv")] <- (1 - coverage) * value # non-vaccinated individuals
    }
    
    # Community compartment 
    else if (grepl("_c$", name)) {
      init[paste0(name, "_v")]  <- coverage * value # vaccinated individuals
      init[paste0(name, "_nv")] <- (1 - coverage) * value # non-vaccinated individuals
    }
  }
  
  # Add cumulative incidence counters (must match the ODE state names)
  cum_names <- c("CumI_primo_h_nv","CumI_rec_h_nv","CumI_h_nv",
                 "CumI_primo_c_nv","CumI_rec_c_nv","CumI_c_nv",
                 "CumI_primo_h_v","CumI_rec_h_v","CumI_h_v",
                 "CumI_primo_c_v","CumI_rec_c_v","CumI_c_v")
  init[cum_names] <- 0
  
  return(init)
}

###############################################################################
#  MODEL FUNCTIONS
###############################################################################

# ---- Model execution and equilibrium detection ----

# Run ODE model until equilibrium using lsoda solver
# Returns: data.frame with time column and all state variables
# run_model_to_equilibrium <- function(params, init_cond, time_vec, 
#                                      model_fn = cdiff_hc_model, 
#                                      alpha_mode = "dynamic",
#                                      atb_reduction_h = 0, 
#                                      atb_reduction_c = 0) {
#   
#   res <- deSolve::lsoda(y = init_cond, times = time_vec, func = model_fn, parms = params,
#                         alpha_mode = alpha_mode,
#                         atb_reduction_h = atb_reduction_h,
#                         atb_reduction_c = atb_reduction_c)
#   
#   as.data.frame(res)
# }

run_model_to_equilibrium <- function(params, init_cond, time_vec, 
                                     model_fn = cdiff_hc_model, 
                                     alpha_mode = "dynamic",
                                     atb_reduction_h = 0, 
                                     atb_reduction_c = 0,
                                     VE = 0,
                                     vacc_type = "both") {
  
  res <- deSolve::lsoda(
    y = init_cond, 
    times = time_vec, 
    func = model_fn, 
    parms = params,
    alpha_mode = alpha_mode,
    atb_reduction_h = atb_reduction_h,
    atb_reduction_c = atb_reduction_c,
    VE = VE,
    vacc_type = vacc_type
  )
  
  as.data.frame(res)
}

# Extract the final equilibrium state from ODE results
# Returns: single-row data.frame with all compartment values at final time point
get_equilibrium_state <- function(ode_result) {
  ode_result[nrow(ode_result), , drop = FALSE]
}

# Extract last state as a named numeric vector (used in calibration)
# Returns: named numeric vector of final compartment values
get_last_state <- function(res_eq) {
  last_df <- get_equilibrium_state(res_eq)
  x <- as.numeric(last_df[1, ])
  names(x) <- colnames(last_df)
  x
}

# Extract last state from raw ODE data.frame (alternative helper)
# Returns: named numeric vector excluding time column
extract_last_state <- function(ode_df) {
  vars <- setdiff(names(ode_df), "time")
  y <- as.numeric(ode_df[nrow(ode_df), vars])
  names(y) <- vars
  y
}

# Check if model has reached equilibrium by evaluating RHS derivatives
# Returns: list with equilibrium status (ok), maximum absolute/relative derivatives, and worst compartment name
check_equilibrium_rhs <- function(ode_df, params, alpha_mode = "dynamic",
                                  tol_abs = 1e-4, tol_rel = 1e-8,
                                  ignore_regex = "^CumI_") {
  
  last_row <- ode_df[nrow(ode_df), , drop = FALSE]
  if ("time" %in% names(last_row)) last_row$time <- NULL
  
  state <- as.numeric(last_row[1, ])
  names(state) <- colnames(last_row)
  
  t_last <- if ("time" %in% names(ode_df)) ode_df$time[nrow(ode_df)] else 0
  rhs <- unlist(cdiff_hc_model(t_last, state, params, alpha_mode = alpha_mode))
  
  # Clean names like "C_III_c.delta" -> "C_III_c"
  names(rhs) <- sub("\\..*$", "", names(rhs))
  rhs <- rhs[!duplicated(names(rhs))]
  
  # Ignore cumulative vars (never equilibrium)
  rhs_core <- rhs[!grepl(ignore_regex, names(rhs))]
  
  # Align state
  state_core <- state[names(rhs_core)]
  if (anyNA(state_core)) {
    missing <- names(rhs_core)[is.na(state_core)]
    stop("State missing: ", paste(missing, collapse = ", "))
  }
  
  abs_err <- abs(rhs_core)
  rel_err <- abs_err / (abs(state_core) + 1e-12)
  
  max_abs <- max(abs_err)
  max_rel <- max(rel_err)
  worst_abs <- names(abs_err)[which.max(abs_err)]
  worst_rel <- names(rel_err)[which.max(rel_err)]
  
  ok <- (max_abs < tol_abs) && (max_rel < tol_rel)
  
  list(
    ok = ok,
    max_abs = max_abs,
    max_rel = max_rel,
    worst_abs = worst_abs,
    worst_rel = worst_rel
  )
}






###############################################################################
# 2. METRICS FUNCTIONS
###############################################################################
# Epidemiological outputs computed from the ODE state.
# All metrics support:
#  - hospital / community / total populations
#  - vaccinated / non-vaccinated / combined populations
###############################################################################

# Check whether the model state includes vaccination stratification (_v/_nv suffixes)
# Returns: TRUE if vaccination compartments detected, FALSE otherwise
has_vaccination <- function(state) {
  return(any(grepl("_(v|nv)$", names(state))))
}

# Compute total population size for a given setting and vaccination status
# Returns: total number of individuals (sum of all non-CI compartments in the specified stratum)
get_total_pop <- function(last_state, setting = c("h", "c", "total"), vacc = c("v", "nv", "both")) {
  
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  vars <- names(last_state) # Start from all state variables
  vars <- vars[!grepl("^CumI_", vars)] # Remove cumulative incidence counters (not real population compartments)
  
  # Select setting (hospital / community / total) 
  if (setting == "h") {vars <- vars[grepl("_h", vars)]}
  if (setting == "c") {vars <- vars[grepl("_c", vars)]}
  
  # Select setting (hospital / community / both)
  if (has_vaccination(last_state)) {
    if (vacc == "v")  {vars <- vars[grepl("_v$", vars)]}
    if (vacc == "nv") {vars <- vars[grepl("_nv$", vars)]}
  } else {   
    # No vaccination structure in the model: only 'both' is meaningful (entire population)
    if (vacc != "both") return(NA_real_)
  }

  return(sum(last_state[vars])) # Sum selected compartments
}

# Compute absolute number of colonized individuals (primary, recurrent, or total)
# Returns: total count of C compartments for the selected setting and vaccination status
get_carriage <- function(last_state, type = c("primo", "rec", "total"), setting = c("h", "c", "total"), vacc = c("v", "nv", "both")) {
  
  type <- match.arg(type)
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  vars <- names(last_state)
  
  # Keep only colonization compartments
  if (type == "primo") {vars <- vars[grepl("^C0_|^CA_", vars)]}
  if (type == "rec")   {vars <- vars[grepl("^C_II_|^C_III_", vars)]}
  if (type == "total") {vars <- vars[grepl("^C0_|^CA_|^C_II_|^C_III_", vars)]}
  
  # Select setting
  if (setting == "h") {vars <- vars[grepl("_h", vars)]}
  if (setting == "c") {vars <- vars[grepl("_c", vars)]}
  
  # Select vaccination status
  if (has_vaccination(last_state)) {
    if (vacc == "v")  {vars <- vars[grepl("_v$", vars)]}
    if (vacc == "nv") {vars <- vars[grepl("_nv$", vars)]}
  } else {
    if (vacc != "both") return(NA_real_)
  }
  
  return(sum(last_state[vars]))
}

# Compute instantaneous CDI incidence (primary, recurrent, or total)
# Returns: absolute number of new infections per unit time (not divided by N)
get_incidence <- function(last_state, params,
                          type = c("primo", "rec", "total"),
                          setting = c("h", "c", "total"),
                          vacc = c("v", "nv", "both"),
                          VE = 0,
                          vacc_type = c("inf", "car", "both")) {
  
  type <- match.arg(type)
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  vacc_type <- match.arg(vacc_type)
  
  # Vaccine effect on progression (sigma) only for "inf" or "both"
  sigma_mult_v <- if (vacc_type %in% c("inf", "both")) (1 - VE) else 1
  
  # Select relevant colonization compartments
  pat <- switch(type,
                primo = "^(C0|CA)_",
                rec   = "^(C_II|C_III)_",
                total = "^(C0|CA|C_II|C_III)_")
  
  vars <- names(last_state)
  vars <- vars[!grepl("^CumI_", vars)]  # remove cumulative counters
  vars <- vars[grepl(pat, vars)]
  
  # Setting filter
  if (setting == "h") vars <- vars[grepl("_h", vars)]
  if (setting == "c") vars <- vars[grepl("_c", vars)]
  
  # Vaccination filter
  if (has_vaccination(last_state)) {
    if (vacc == "v")  vars <- vars[grepl("_v$", vars)]
    if (vacc == "nv") vars <- vars[grepl("_nv$", vars)]
  } else {
    if (vacc != "both") return(NA_real_)
  }
  
  # Pull numeric params WITHOUT names
  k_A   <- params[["k_A"]]
  k_II  <- params[["k_II"]]
  k_III <- params[["k_III"]]
  sigma_h <- params[["sigma_h"]]
  sigma_c <- params[["sigma_c"]]
  
  inc <- 0.0
  
  for (vname in vars) {
    
    # Pick sigma by setting (numeric, unnamed)
    sigma <- if (grepl("_h", vname)) sigma_h else sigma_c
    
    # Apply VE on sigma if vaccinated
    if (has_vaccination(last_state) && grepl("_v$", vname)) {
      sigma <- sigma * sigma_mult_v
    }
    
    # Apply stage multipliers
    if (grepl("^C0_", vname))    inc <- inc + sigma * last_state[[vname]]
    if (grepl("^CA_", vname))    inc <- inc + (k_A   * sigma) * last_state[[vname]]
    if (grepl("^C_II_", vname))  inc <- inc + (k_II  * sigma) * last_state[[vname]]
    if (grepl("^C_III_", vname)) inc <- inc + (k_III * sigma) * last_state[[vname]]
  }
  
  as.numeric(inc)
}




# Compute cumulative CDI incidence from CumI_* compartments
# Returns: total cumulative number of CDI cases over the simulation period
get_cumulative_incidence <- function(last_state, type = c("total", "primo", "rec"),
                                     setting = c("h", "c", "total"), vacc = c("v", "nv", "both")) {
  
  type <- match.arg(type)
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  vars <- names(last_state)
  
  # Choose ONE counter family to avoid double counting
  if (type == "primo") {
    vars <- vars[grepl("^CumI_primo_", vars)]
  }
  if (type == "rec") {
    vars <- vars[grepl("^CumI_rec_", vars)]
  }
  if (type == "total") {
    # Keep only total counters (exclude primo/rec)
    vars <- vars[grepl("^CumI_", vars) & !grepl("^CumI_(primo|rec)_", vars)]
  }
  
  # Setting filter
  if (setting == "h") vars <- vars[grepl("_h", vars)]
  if (setting == "c") vars <- vars[grepl("_c", vars)]
  
  # Vaccination filter
  if (has_vaccination(last_state)) {
    if (vacc == "v")  vars <- vars[grepl("_v$", vars)]
    if (vacc == "nv") vars <- vars[grepl("_nv$", vars)]
  } else {
    if (vacc != "both") return(NA_real_)
  }
  
  return(sum(last_state[vars]))
}

# Compute recurrence prevalence from infected compartments
# rec_1 = I_II / I , rec_2 = I_III / I_II
compute_recurrence_prevalence <- function(last_state,
                                          level = c("rec_1", "rec_2"),
                                          setting = c("h", "c", "total"),
                                          vacc = c("v", "nv", "both")) {
  
  level <- match.arg(level)
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  vars <- names(last_state)
  vars <- vars[!grepl("^CumI_", vars)] # Remove fake compartiments (security)
  
  # Select numerator and denominator compartments
  if (level == "rec_1") {
    num_vars <- vars[grepl("^I_II_", vars)]
    den_vars <- vars[grepl("^I_(h|c)(_|$)", vars)]}  # I_h, I_c, I_h_v, I_c_nv, etc.
  
  if (level == "rec_2") {
    num_vars <- vars[grepl("^I_III_", vars)]
    den_vars <- vars[grepl("^I_II_", vars)]}
  
  # Select setting
  if (setting == "h") {
    num_vars <- num_vars[grepl("_h", num_vars)]
    den_vars <- den_vars[grepl("_h", den_vars)]}
  if (setting == "c") {
    num_vars <- num_vars[grepl("_c", num_vars)]
    den_vars <- den_vars[grepl("_c", den_vars)]}
  
  # Select vaccination status
  if (has_vaccination(last_state)) {
    if (vacc == "v") {
      num_vars <- num_vars[grepl("_v$", num_vars)]
      den_vars <- den_vars[grepl("_v$", den_vars)]}
    if (vacc == "nv") {
      num_vars <- num_vars[grepl("_nv$", num_vars)]
      den_vars <- den_vars[grepl("_nv$", den_vars)]}
  } else {
    if (vacc != "both") return(NA_real_)
  }
  
  num <- sum(last_state[num_vars])
  den <- sum(last_state[den_vars])
  
  if (den == 0) return(NA_real_)
  
  return(num / den)
}


# Get absolute count of individuals with antibiotic-induced dysbiosis
# Returns: total number of individuals in SA and CA compartments
get_dysbiosis <- function(last_state, setting = c("h", "c", "total"), vacc = c("v", "nv", "both")) {
  
  setting <- match.arg(setting)
  vacc <- match.arg(vacc)
  
  vars <- names(last_state)
  vars <- vars[!grepl("^CumI_", vars)] # Remove fake compartiments (security)
  vars <- vars[grepl("^SA_|^CA_", vars)] # Keep only dysbiosis compartments
  
  # Select setting
  if (setting == "h") {vars <- vars[grepl("_h", vars)]}
  if (setting == "c") {vars <- vars[grepl("_c", vars)]}
  
  # Select vaccination status
  if (has_vaccination(last_state)) {
    if (vacc == "v")  {vars <- vars[grepl("_v$", vars)]}
    if (vacc == "nv") {vars <- vars[grepl("_nv$", vars)]}
  } else {
    if (vacc != "both") return(NA_real_)
  }
  
  return(sum(last_state[vars]))
}




