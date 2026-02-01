###############################################################################
################### 3 : CALIBRATION AND GRID SEARCH ###########################
###############################################################################

###############################################################################
# 1. HELPER FUNCTIONS
###############################################################################

# Run the model and stop early if equilibrium is reached (ignoring Cum* compartments)
run_model_to_equilibrium <- function(params, init_cond, time_vec) {
  
  out_all <- NULL                                 # we will store all simulation rows here
  pop <- init_cond                                # current state (starts at initial conditions)
  keep_idx <- !grepl("Cum", names(pop))           # TRUE = real compartments, FALSE = Cum* counters
  
  step <- 50                                      # simulate 50 days at a time
  
  t_start <- time_vec[1]                          # first time
  t_end <- time_vec[length(time_vec)]             # last time
  
  t <- t_start                                    # current time
  
  while (t < t_end) {                             # loop until we reach the end
    
    t_next <- min(t + step, t_end)                # next stop time
    
    times_chunk <- seq(t, t_next, by = 1)         # times for this chunk (daily)
    
    out <- deSolve::lsoda(y = pop, times = times_chunk,
                          func = cdiff_model_for_calibration, parms = params) # solve ODE on this chunk
    out <- as.data.frame(out)                     # convert to data.frame
    
    if (is.null(out_all)) {                       # first chunk
      out_all <- out                              # store everything
    } else {                                      # next chunks
      out_all <- rbind(out_all, out[-1, ])        # add rows, remove first row (duplicate time)
    }
    
    pop <- as.numeric(out[nrow(out), names(init_cond)])  # take last state as new initial state
    names(pop) <- names(init_cond)                       # keep names
    
    derivs_end <- cdiff_model_for_calibration(t = t_next, pop = pop, params = params)[[1]] # compute derivatives
    names(derivs_end) <- names(pop)                   # name derivatives
    
    abs_end <- abs(derivs_end[keep_idx])              # |dX/dt| without Cum*
    at_eq_end <- all(abs_end < 1e-6)                  # TRUE if equilibrium reached
    
    if (at_eq_end) {                                  # if equilibrium reached
      break                                           # stop early
    }
    
    t <- t_next                                       # move time forward
  }
  
  return(out_all)                                     # return the simulated trajectory
}


# Extract the last state (last row) as a named numeric vector (without "time")
get_last_state <- function(out) {
  last_row <- out[nrow(out), ]                 # take last row
  last_row <- last_row[, names(last_row) != "time"] # remove time column
  x <- as.numeric(last_row)                    # convert to numeric vector
  names(x) <- names(last_row)                  # keep compartment names
  return(x)                                    # return named state vector
}





###############################################################################
# 2. GENERIC GRID SEARCH (PARALLEL)
###############################################################################

grid_search <- function(param_names, param_ranges, metric_function,
                        target_metrics, params_base, init_cond, time_vec, n_cores) {
  
  # Build sequences for each parameter
  grid_list <- list()                                                             # empty list for sequences
  grid_list[[param_names[1]]] <- seq(param_ranges[[1]][1], param_ranges[[1]][2], length.out = param_ranges[[1]][3]) # seq 1
  grid_list[[param_names[2]]] <- seq(param_ranges[[2]][1], param_ranges[[2]][2], length.out = param_ranges[[2]][3]) # seq 2
  
  # Create the full grid
  param_grid <- expand.grid(grid_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) # all combinations
  
  # Create a parallel cluster
  cl <- parallel::makeCluster(n_cores)                                      # start workers
  parallel::clusterEvalQ(cl, library(deSolve))                              # load deSolve on workers
  parallel::clusterExport(cl, c("param_grid", "param_names", "metric_function",
                                "target_metrics", "params_base", "init_cond", "time_vec",
                                "run_model_to_equilibrium", "get_last_state",
                                "compute_totals", "compute_lambda",
                                "cdiff_model_for_calibration", "compute_metrics_calib"),
                          envir = environment())                             # send objects/functions to workers
  
  # Run metrics in parallel for each grid row
  res_list <- parallel::parLapply(cl, X = seq_len(nrow(param_grid)), fun = function(i) { # loop in parallel
    
    row <- as.numeric(param_grid[i, param_names])                            # take the two parameters
    names(row) <- param_names                                                # set names
    
    m <- metric_function(row, params_base, init_cond, time_vec)              # compute metrics
    return(m)                                                                # return named numeric vector
  })
  
  parallel::stopCluster(cl)                                                  # stop workers
  
  # Put results into a data.frame
  res_mat <- do.call(rbind, res_list)                                        # bind results row-wise
  res_df <- as.data.frame(res_mat)                                           # convert to data.frame
  
  # Combine grid + metrics
  out_grid <- cbind(param_grid, res_df)                                      # final grid table
  
  # Compute distance to targets (Euclidean)
  metric_names <- names(target_metrics)                                      # metric columns to compare
  target_vec <- unlist(target_metrics)                                       # targets as numeric vector
  dist_vec <- rep(NA_real_, nrow(out_grid))                                  # allocate distance vector
  
  for (i in seq_len(nrow(out_grid))) {                                       # loop over grid rows
    diff <- as.numeric(out_grid[i, metric_names]) - as.numeric(target_vec)   # difference to targets
    dist_vec[i] <- sqrt(sum(diff * diff))                                    # Euclidean distance
  }
  
  out_grid$distance <- dist_vec                                              # store distance
  
  best_idx <- which.min(out_grid$distance)                                   # index of best row
  best_guess <- out_grid[best_idx, , drop = FALSE]                           # best row
  
  return(list(grid = out_grid, best_guess = best_guess))                     # return everything
}

###############################################################################
# 3. METRIC FUNCTIONS FOR THE 3 GRID SEARCH
###############################################################################

# Grid 1: calibrate beta_h, beta_c using carriage prevalence (prev_h, prev_c)
compute_metrics_beta <- function(param_row, params_base, init_cond, time_vec) {
  
  params_tmp <- params_base                         # copy base params
  params_tmp["beta_h"] <- param_row["beta_h"]       # set beta_h
  params_tmp["beta_c"] <- param_row["beta_c"]       # set beta_c
  
  out <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)            # simulate
  metrics <- compute_metrics_calib(out, params_tmp)                           # compute metrics
  
  prev_h <- tail(metrics$carriage$prev_h, 1)                                  # last prev_h
  prev_c <- tail(metrics$carriage$prev_c, 1)                                  # last prev_c
  
  m <- c(portage_h = prev_h, portage_c = prev_c)                              # named vector
  return(m)                                                                   # return
}

# Grid 2: calibrate sigma_h, sigma_c using incidence rates per N_tot (abs)
compute_metrics_sigma <- function(param_row, params_base, init_cond, time_vec) {
  
  params_tmp <- params_base                           # copy base params
  params_tmp["sigma_h"] <- param_row["sigma_h"]       # set sigma_h
  params_tmp["sigma_c"] <- param_row["sigma_c"]       # set sigma_c
  
  out <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)            # simulate
  metrics <- compute_metrics_calib(out, params_tmp)                           # compute metrics
  
  inc_h <- tail(metrics$incidence_instant$inc_h_total_abs, 1)                 # last hospital incidence / N_tot
  inc_c <- tail(metrics$incidence_instant$inc_c_total_abs, 1)                 # last community incidence / N_tot
  
  m <- c(incidence_h = inc_h, incidence_c = inc_c)                            # named vector
  return(m)                                                                   # return
}

# Grid 3: calibrate k_II, k_III using recurrence prevalence (overall, h+c)
compute_metrics_k <- function(param_row, params_base, init_cond, time_vec) {
  
  params_tmp <- params_base                           # copy base params
  params_tmp["k_II"] <- param_row["k_II"]             # set k_II
  params_tmp["k_III"] <- param_row["k_III"]           # set k_III
  
  out <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)            # simulate
  last <- get_last_state(out)                                                 # last state
  
  I1 <- last["I_h"] + last["I_c"]                                             # primary infected (h+c)
  I2 <- last["I_II_h"] + last["I_II_c"]                                       # first recurrence infected (h+c)
  I3 <- last["I_III_h"] + last["I_III_c"]                                     # second+ recurrence infected (h+c)
  
  recid_1 <- ifelse(I1 == 0, 0, I2 / I1)                                      # recurrence 1 ratio
  recid_2 <- ifelse(I2 == 0, 0, I3 / I2)                                      # recurrence 2 ratio
  
  m <- c(recid_1 = recid_1, recid_2 = recid_2)                                # named vector
  return(m)                                                                   # return
}





###############################################################################
# 4. OBJECTIVE FUNCTION (FOR MULTI-START OPTIM)
###############################################################################

create_objective_function <- function(target_metrics, params_base, init_cond, time_vec) {
  
  par_names <- c("beta_h","beta_c","sigma_h","sigma_c","k_II","k_III")         # parameters to estimate
  
  objective_fn <- function(par_log) {                                         # par_log is log(parameters)
    
    par <- exp(par_log)                                                       # back to natural scale
    names(par) <- par_names                                                   # name parameters
    
    params_tmp <- params_base                                                 # copy base params
    params_tmp[par_names] <- par                                              # put calibrated params into params
    
    out <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)          # simulate
    metrics <- compute_metrics_calib(out, params_tmp)                         # compute metrics
    last <- get_last_state(out)                                               # last state (for recurrence totals)
    
    # Targets: carriage
    port_h <- tail(metrics$carriage$prev_h, 1)                                # portage_h
    port_c <- tail(metrics$carriage$prev_c, 1)                                # portage_c
    
    # Targets: incidence (per N_tot)
    inc_h <- tail(metrics$incidence_instant$inc_h_total_abs, 1)               # incidence_h
    inc_c <- tail(metrics$incidence_instant$inc_c_total_abs, 1)               # incidence_c
    
    # Targets: recurrence (overall)
    I1 <- last["I_h"] + last["I_c"]                                           # I primary (h+c)
    I2 <- last["I_II_h"] + last["I_II_c"]                                     # I recurrence 1 (h+c)
    I3 <- last["I_III_h"] + last["I_III_c"]                                   # I recurrence 2+ (h+c)
    rec1 <- ifelse(I1 == 0, 0, I2 / I1)                                       # recid_1
    rec2 <- ifelse(I2 == 0, 0, I3 / I2)                                       # recid_2
    
    # Relative errors (simple)
    e1 <- (port_h - target_metrics$portage_h) / target_metrics$portage_h      # error portage_h
    e2 <- (port_c - target_metrics$portage_c) / target_metrics$portage_c      # error portage_c
    e3 <- (inc_h  - target_metrics$incidence_h) / target_metrics$incidence_h  # error incidence_h
    e4 <- (inc_c  - target_metrics$incidence_c) / target_metrics$incidence_c  # error incidence_c
    e5 <- (rec1   - target_metrics$recid_1) / target_metrics$recid_1          # error recid_1
    e6 <- (rec2   - target_metrics$recid_2) / target_metrics$recid_2          # error recid_2
    
    sse <- e1*e1 + e2*e2 + e3*e3 + e4*e4 + e5*e5 + e6*e6                      # sum of squared relative errors
    return(sse)                                                               # return objective value
  }
  
  return(objective_fn)                                                        # return the objective function
}

###############################################################################
# 4. MULTI-START OPTIMIZATION (PARALLEL)
###############################################################################

run_multistart_optimization <- function(objective_fn, initial_log, n_starts, n_cores, maxit) {
  
  # Create starting points in log-space
  start_points <- list()                                                      # list of starts
  for (i in seq_len(n_starts)) {                                              # loop over starts
    start_points[[i]] <- initial_log + runif(length(initial_log), -0.3, 0.3)  # small random jitter
  }
  start_points[[1]] <- initial_log                                            # first start = exact initial guess
  
  # Create a parallel cluster
  cl <- parallel::makeCluster(n_cores)                                        # start workers
  parallel::clusterExport(cl, c("objective_fn", "start_points", "maxit",
                                "run_model_to_equilibrium", "get_last_state",
                                "compute_totals", "compute_lambda",
                                "cdiff_model_for_calibration", "compute_metrics_calib"),
                          envir = environment())                               # send needed objects
  parallel::clusterEvalQ(cl, library(deSolve))                                 # load deSolve on workers
  
  # Run optim in parallel (one optim per start)
  res_list <- parallel::parLapply(cl, X = seq_len(n_starts), fun = function(i) { # parallel loop
    
    opt <- optim(par = start_points[[i]], fn = objective_fn, method = "Nelder-Mead",
                 control = list(maxit = maxit))                               # run Nelder-Mead
    
    out <- list()                                                             # create output list
    out$par_log <- opt$par                                                    # best log-params
    out$par <- exp(opt$par)                                                   # best params (natural scale)
    out$value <- opt$value                                                    # objective value
    out$convergence <- opt$convergence                                        # convergence code
    return(out)                                                               # return
  })
  
  parallel::stopCluster(cl)                                                   # stop workers
  
  # Pick the best result
  values <- rep(NA_real_, length(res_list))                                   # allocate values
  for (i in seq_along(res_list)) {                                            # loop over results
    values[i] <- res_list[[i]]$value                                          # store value
  }
  best_idx <- which.min(values)                                               # best index
  
  best <- res_list[[best_idx]]                                                # best run
  
  return(list(best = best, all_results = res_list))                           # return everything
}

###############################################################################
# 5) CALIBRATION WRAPPER (ONE FUNCTION TO RUN EVERYTHING)
###############################################################################

run_calibration <- function(initial_params, target_metrics, params_base,
                            init_cond, time_vec, n_starts, n_cores, maxit) {
  
  par_names <- c("beta_h","beta_c","sigma_h","sigma_c","k_II","k_III")         # calibrated parameter names
  initial_log <- log(initial_params[par_names])                               # initial guess in log-space
  
  objective_fn <- create_objective_function(target_metrics, params_base, init_cond, time_vec) # build objective
  
  res <- run_multistart_optimization(objective_fn, initial_log, n_starts, n_cores, maxit)     # run multi-start
  
  return(res)                                                                 # return results
}
