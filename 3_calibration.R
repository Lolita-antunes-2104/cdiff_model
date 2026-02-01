###############################################################################
################### 3 : CALIBRATION AND GRID SEARCH ###########################
###############################################################################

###############################################################################
# 1. HELPER FUNCTIONS
###############################################################################

# Run the CALIBRATION model and stop early if equilibrium is reached 
run_calib_model_to_equilibrium <- function(params, init_cond, time_vec) {
  
  out_all <- NULL                                 # create an empty object that will progressively store the full simulated trajectory (all time rows) as we simulate chunk by chunk
  pop <- init_cond                                # copy the initial state vector into 'pop'; this is the current state that will be used as the initial condition for the next ODE solve
  step <- 100                                     # set the chunk length (in days) so we solve the ODE in chunk instead of the full horizon at once
  
  t <- time_vec[1]                                # Read the first element of time_vec: this is the simulation start time
  t_end <- time_vec[length(time_vec)]             # Read the last element of time_vec using length(): this is the simulation end time
  
  while (t < t_end) {                             # Start a loop that continues until the current time reaches the final time
    
    t_next <- min(t + step, t_end)                # Compute the next stopping time: min() returns the smallest of its arguments so we never go beyond t_end
    times_chunk <- seq(t, t_next, by = 1)         # Create a sequence of times for this chunk using seq(): here we build daily times from t to t_next with step size 1
    
    out <- deSolve::lsoda(                        # Call lsoda() from the deSolve package: it numerically solves a system of ODEs and returns the state at each requested time
      y = pop,                                                      # 'y' is the initial state vector at the start of this chunk (the current compartments values)
      times = times_chunk,                                          # 'times' is the vector of time points where lsoda will output the solution (here, each day)
      func = cdiff_model_for_calibration,                           # 'func' is the ODE right-hand-side function that returns dX/dt (your calibration hospital/community model)
      parms = params                                                # 'parms' is the parameters object passed into the ODE function (so the model can read beta, sigma, delta, etc.)
    )
    
    out <- as.data.frame(out)                      # Convert the solver output to a data.frame (tabular format) using as.data.frame() so we can easily row-bind chunks and index columns
    out_all <- if (is.null(out_all)) out else rbind(out_all, out[-1, ]) # If this is the first chunk (out_all is NULL), store out directly; otherwise rbind() appends rows and out[-1,] removes the first row to avoid duplicating the boundary time
    
    pop <- as.numeric(out[nrow(out), names(init_cond)])         # Extract the last row using nrow() (number of rows) and take the columns in the same order as init_cond; as.numeric() forces a plain numeric vector (drop data.frame structure)
    names(pop) <- names(init_cond)                              # Restore names on the state vector so compartments keep their labels (important because the ODE function expects named states like S0_h, C0_h, ...).
    
    derivs_end <- cdiff_model_for_calibration(t_next, pop, params)[[1]] # Compute the derivatives dX/dt at the end of the chunk by calling the ODE function directly; [[1]] extracts the numeric derivative vector from the list returned by the ODE function.
    names(derivs_end) <- names(pop)                             # Name the derivative entries so we can safely subset by compartment names and keep alignment with 'pop'.
    
    if (all(abs(derivs_end) < 1e-8)) break            # abs() takes absolute values; we subset to non-Cum compartments; all() returns TRUE only if every value is < 1e-6, and if so we break out of the loop (equilibrium reached).
    t <- t_next                                                 # Update the current time to the end of this chunk so the next loop iteration starts where we left off.
  }
  return(out_all)                                  # Return the full concatenated trajectory (a data.frame with one row per time and one column per compartment, including "time").
}


# Extract the last state (last row) as a named numeric vector (without "time")
get_last_state <- function(out) {   
  last <- out[nrow(out), ]                   # Take the last row of the output using nrow() (the final time point) and keep it as a 1-row data.frame
  last <- last[, names(last) != "time"]      # Remove the "time" column by selecting all columns whose name is not exactly "time" (keep only compartments)
  x <- as.numeric(last)                      # Convert the 1-row data.frame into a plain numeric vector (dropping data.frame attributes) using as.numeric()
  names(x) <- names(last)                    # Set the vector names to the compartment names so the returned vector is "named" (easier to pass back into lsoda as initial conditions)
  return(x)                                  # Return the named numeric vector containing the last values of all compartments (no time column)
}





###############################################################################
# 2. GENERIC GRID SEARCH (PARALLEL)
###############################################################################

grid_search <- function(param_names, param_ranges, metric_function,           # define a function that explores a grid of 2 parameters (param_names), runs a metric_function, and finds the best match to target_metrics using parallel computing
                        target_metrics, params_base, init_cond, time_vec, n_cores) { # function arguments: targets, baseline params, initial conditions, time vector for ODE, and number of CPU cores to use
  
  # Build sequences for each parameter
  grid_list <- list()                                                         # create an empty list that will store the sequences (vectors) of values to try for each parameter
  grid_list[[param_names[1]]] <- seq(param_ranges[[1]][1], param_ranges[[1]][2], length.out = param_ranges[[1]][3]) # create a sequence for parameter 1: seq(from, to, length.out = number_of_points)
  grid_list[[param_names[2]]] <- seq(param_ranges[[2]][1], param_ranges[[2]][2], length.out = param_ranges[[2]][3]) # create a sequence for parameter 2: seq(from, to, length.out = number_of_points)
  
  # Create the full grid
  param_grid <- expand.grid(grid_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) # expand.grid() creates all combinations of the sequences; KEEP.OUT.ATTRS=FALSE keeps the result simple; stringsAsFactors=FALSE prevents automatic factor conversion
  
  # Create a parallel cluster
  cl <- parallel::makeCluster(n_cores)                                        # makeCluster() starts a parallel cluster with n_cores worker R processes
  parallel::clusterEvalQ(cl, library(deSolve))                                # clusterEvalQ() runs the given expression on each worker; here it loads the deSolve package so lsoda() is available on workers
  parallel::clusterExport(cl, c("param_grid", "param_names", "metric_function",# clusterExport() sends objects/functions from the current R session to all workers so they can use them
                                "target_metrics", "params_base", "init_cond", "time_vec", # export the inputs needed inside the parallel loop (baseline params, initial state, time vector, etc.)
                                "run_calib_model_to_equilibrium", "get_last_state", # export helper functions used by metric_function (or by the model-run pipeline inside it)
                                "compute_totals", "compute_lambda", "compute_R0",   # export model helper functions used by the ODE function
                                "cdiff_model_for_calibration", "compute_metrics_calib"), # export the ODE function and the metrics function used to evaluate each simulation
                          envir = environment())                               # envir=environment() tells clusterExport() to look for these objects in the current function environment (not only in .GlobalEnv)
  
  # Run metrics in parallel for each grid row
  res_list <- parallel::parLapply(cl, X = seq_len(nrow(param_grid)), fun = function(i) { # parLapply() applies a function in parallel over X; here X is 1..number_of_grid_rows, so each worker processes different grid rows
    
    row <- as.numeric(param_grid[i, param_names])                              # take the i-th row of the grid, select the 2 parameter columns, and convert to a numeric vector
    names(row) <- param_names                                                  # assign names to the vector so metric_function can refer to parameters by name (e.g., row["beta_h"])
    
    m <- metric_function(row, params_base, init_cond, time_vec)                # call the user-provided metric_function: it usually runs the ODE with these parameters and returns computed metrics (as a named numeric vector)
    return(m)                                                                  # return the metrics for this grid row back to the master process
  })
  
  parallel::stopCluster(cl)                                                    # stopCluster() shuts down the worker processes and releases CPU/memory resources
  
  # Put results into a data.frame
  res_mat <- do.call(rbind, res_list)                                          # do.call(rbind, ...) stacks (row-binds) the list of metric vectors into a matrix (one row per grid point)
  res_df <- as.data.frame(res_mat)                                             # convert the matrix to a data.frame so it can be combined with the parameter grid easily
  
  # Combine grid + metrics
  out_grid <- cbind(param_grid, res_df)                                        # cbind() column-binds the parameter values and the computed metrics into one final results table
  
  # Compute distance to targets (Euclidean)
  metric_names <- names(target_metrics)                                        # extract the metric names (column names) we want to compare to targets (e.g., prevalence_h, incidence_c, ...)
  target_vec <- unlist(target_metrics)                                         # unlist() converts the target_metrics list into a plain numeric vector with the same order as metric_names
  dist_vec <- rep(NA_real_, nrow(out_grid))                                    # create a numeric vector of length nrow(out_grid) filled with NA_real_ to store distances for each grid row
  
  for (i in seq_len(nrow(out_grid))) {                                         # loop over all grid rows (sequentially on the master process) to compute the distance to the targets
    diff <- as.numeric(out_grid[i, metric_names]) - as.numeric(target_vec)     # compute the difference between simulated metrics (row i) and target values (same metric order)
    dist_vec[i] <- sqrt(sum(diff * diff))                                      # compute Euclidean distance: sqrt(sum((diff)^2))
  }
  
  out_grid$distance <- dist_vec                                                # add the computed distance vector as a new column named "distance" in the results table
  
  best_idx <- which.min(out_grid$distance)                                     # which.min() returns the row index of the smallest distance (best match to targets)
  best_guess <- out_grid[best_idx, , drop = FALSE]                             # subset the best row; drop=FALSE keeps it as a data.frame (not simplified to a vector)
  
  return(list(grid = out_grid, best_guess = best_guess))                       # return a list containing the full grid results table and the single best row
}





###############################################################################
# 3. METRIC FUNCTIONS FOR THE 3 GRID SEARCH
###############################################################################

# Grid 1: calibrate beta_h, beta_c using carriage prevalence (prev_h, prev_c)
compute_metrics_beta <- function(param_row, params_base, init_cond, time_vec) { # define the metric function used in Grid Search 1 (it receives 1 grid row + the baseline params + init conditions + time vector)
  
  params_tmp <- params_base                                                    # copy the baseline parameter vector into a temporary one (so we don't overwrite the original)
  params_tmp["beta_h"] <- param_row["beta_h"]                                  # replace beta_h in the temporary params by the grid value beta_h from this row
  params_tmp["beta_c"] <- param_row["beta_c"]                                  # replace beta_c in the temporary params by the grid value beta_c from this row
  
  out <- run_calib_model_to_equilibrium(params_tmp, init_cond, time_vec)             # run the ODE model (lsoda) with these temporary params, starting from init_cond, over time_vec
  metrics <- compute_metrics_calib(out, params_tmp)                            # compute all calibration metrics at the end of the simulation using the same params
  
  prev_h <- tail(metrics$carriage$prev_h, 1)                                   # take the last hospital carriage prevalence value (tail(...,1) = last element)
  prev_c <- tail(metrics$carriage$prev_c, 1)                                   # take the last community carriage prevalence value
  
  m <- c(prevalence_h = prev_h, prevalence_c = prev_c)                         # build a named numeric vector with exactly the metrics we want Grid 1 to match
  return(m)                                                                    # return this named numeric vector to grid_search()
}

# Grid 2: calibrate sigma_h, sigma_c using incidence rates per N_tot (abs)
compute_metrics_sigma <- function(param_row, params_base, init_cond, time_vec) { # define the metric function used in Grid Search 2 (it calibrates sigma parameters)
  
  params_tmp <- params_base                                                      # copy the baseline parameters into a temporary vector
  # IMPORTANT (logic): here, params_base is assumed to ALREADY contain the best beta_h and beta_c from Grid 1 (i.e., you updated params_base with best_guess beta before calling Grid 2)
  params_tmp["sigma_h"] <- param_row["sigma_h"]                                  # replace sigma_h by the grid value sigma_h from this row
  params_tmp["sigma_c"] <- param_row["sigma_c"]                                  # replace sigma_c by the grid value sigma_c from this row
  
  out <- run_calib_model_to_equilibrium(params_tmp, init_cond, time_vec)               # run the ODE model with (best betas + candidate sigmas)
  metrics <- compute_metrics_calib(out, params_tmp)                              # compute metrics at equilibrium / last time
  
  inc_h <- tail(metrics$incidence_instant$inc_h_total_abs, 1)                    # take the last hospital total incidence normalized by N_tot (absolute normalization)
  inc_c <- tail(metrics$incidence_instant$inc_c_total_abs, 1)                    # take the last community total incidence normalized by N_tot
  
  m <- c(incidence_h = inc_h, incidence_c = inc_c)                               # build a named numeric vector with exactly the metrics we want Grid 2 to match
  return(m)                                                                      # return this named numeric vector to grid_search()
}

# Grid 3: calibrate k_II, k_III using recurrence prevalence (overall, h+c)
compute_metrics_k <- function(param_row, params_base, init_cond, time_vec) {     # define the metric function used in Grid Search 3 (it calibrates recurrence multipliers)
  
  params_tmp <- params_base                                                      # copy the baseline parameters into a temporary vector
  # IMPORTANT (logic): here, params_base is assumed to ALREADY contain the best beta_h/beta_c from Grid 1 AND the best sigma_h/sigma_c from Grid 2 (i.e., you updated params_base with both best_guess beta and best_guess sigma before calling Grid 3)
  params_tmp["k_II"] <- param_row["k_II"]                                        # replace k_II by the grid value k_II from this row
  params_tmp["k_III"] <- param_row["k_III"]                                      # replace k_III by the grid value k_III from this row
  
  out <- run_calib_model_to_equilibrium(params_tmp, init_cond, time_vec)               # run the ODE model with (best betas + best sigmas + candidate k values)
  metrics <- compute_metrics_calib(out, params_tmp)                              # compute metrics at equilibrium / last time (includes recurrence prevalence)
  
  rec1 <- tail(metrics$recurrence$rec1, 1)                                       # take the last recurrence prevalence rec1 (I_II / I) from the metrics output
  rec2 <- tail(metrics$recurrence$rec2, 1)                                       # take the last recurrence prevalence rec2 (I_III / I_II) from the metrics output
  
  m <- c(recid_1 = rec1, recid_2 = rec2)                                         # build a named numeric vector with exactly the metrics we want Grid 3 to match
  return(m)                                                                      # return this named numeric vector to grid_search()
}





###############################################################################
# 4. OBJECTIVE FUNCTION (FOR MULTI-START OPTIM)
###############################################################################

create_objective_function <- function(target_metrics, params_base, init_cond, time_vec) { # define a function that BUILDS and RETURNS an objective function, using fixed inputs (targets + base params + init state + time grid)
  
  par_names <- c("beta_h","beta_c","sigma_h","sigma_c","k_II","k_III")         # define the names of the parameters we want to estimate during calibration (these will be optimized)
  
  objective_fn <- function(par_log) {                                         # define the objective function used by optim(); input par_log is the parameters on the log scale (so they can be any real number)
    
    par <- exp(par_log)                                                       # exp() transforms log-parameters back to the natural scale (so parameters are positive)
    names(par) <- par_names                                                   # assign the correct parameter names to the vector so we can insert them into params_tmp by name
    
    params_tmp <- params_base                                                 # copy the baseline parameter vector so we do not modify params_base directly
    params_tmp[par_names] <- par                                              # replace the 6 parameters in params_tmp by the candidate values proposed by the optimizer
    
    out <- run_calib_model_to_equilibrium(params_tmp, init_cond, time_vec)          # run the ODE model with these candidate parameters to reach equilibrium (or the end of time_vec)
    metrics <- compute_metrics_calib(out, params_tmp)                         # compute the calibration metrics from the simulated trajectory (carriage, incidence, recurrence, etc.)

    # Targets: carriage
    port_h <- tail(metrics$carriage$prev_h, 1)                                # take the last hospital carriage prevalence value (tail(...,1) = last element)
    port_c <- tail(metrics$carriage$prev_c, 1)                                # take the last community carriage prevalence value
    
    # Targets: incidence (per N_tot)
    inc_h <- tail(metrics$incidence_instant$inc_h_total_abs, 1)               # take the last hospital total incidence normalized by N_tot (absolute normalization)
    inc_c <- tail(metrics$incidence_instant$inc_c_total_abs, 1)               # take the last community total incidence normalized by N_tot
    
    # Targets: recurrence (overall)
    rec1 <- tail(metrics$recurrence$rec1, 1)                                  # take the last recurrence metric rec1 from metrics (WARNING: 'recurrence' looks misspelled; if your metrics output uses 'recurrence', this will error)
    rec2 <- tail(metrics$recurrence$rec2, 1)                                  # take the last recurrence metric rec2 (same warning about spelling)
    
    # Relative errors (simple)
    e1 <- (port_h - target_metrics$prevalence_h) / target_metrics$prevalence_h      # compute relative error for hospital carriage: (model - target) / target
    e2 <- (port_c - target_metrics$prevalence_c) / target_metrics$prevalence_c      # compute relative error for community carriage
    e3 <- (inc_h  - target_metrics$incidence_h) / target_metrics$incidence_h  # compute relative error for hospital incidence (normalized by N_tot)
    e4 <- (inc_c  - target_metrics$incidence_c) / target_metrics$incidence_c  # compute relative error for community incidence (normalized by N_tot)
    e5 <- (rec1   - target_metrics$recid_1) / target_metrics$recid_1          # compute relative error for recurrence metric 1
    e6 <- (rec2   - target_metrics$recid_2) / target_metrics$recid_2          # compute relative error for recurrence metric 2
    
    sse <- e1*e1 + e2*e2 + e3*e3 + e4*e4 + e5*e5 + e6*e6                      # compute the objective value as the sum of squared relative errors (SSE) across all 6 targets
    return(sse)                                                               # return the objective value to optim(); optim() will try to MINIMIZE this number
  }
  
  return(objective_fn)                                                        # return the objective function (a function you can pass into optim())
}





###############################################################################
# 5. MULTI-START OPTIMIZATION (PARALLEL)
###############################################################################

run_multistart_optimization <- function(objective_fn, initial_log, n_starts, n_cores, maxit) { # define a function that runs many optimizations (multi-start) in parallel, each starting from a different initial point in log-space
  
  # Create starting points in log-space
  start_points <- list()                                                      # create an empty list that will contain the starting vectors (one per optimization run)
  for (i in seq_len(n_starts)) {                                              # seq_len(n_starts) creates the integer sequence 1,2,...,n_starts; we loop over each start index i
    start_points[[i]] <- initial_log + runif(length(initial_log), -0.3, 0.3)  # runif(n, min, max) generates n random numbers uniformly; here we add small random noise to initial_log to create different starting points
  }
  start_points[[1]] <- initial_log                                            # force the first start to be exactly the initial guess (no randomness), so we always try the user's best guess
  
  # Create a parallel cluster
  cl <- parallel::makeCluster(n_cores)                                        # makeCluster(n_cores) starts n_cores worker R processes for parallel computation
  parallel::clusterExport(cl, c("objective_fn", "start_points", "maxit",       # clusterExport() copies these objects from the current environment to every worker so they can be used inside the parallel function
                                "run_calib_model_to_equilibrium", "get_last_state", # export helper functions that objective_fn might call (directly or indirectly)
                                "compute_totals", "compute_lambda", "compute_R0",   # export model helper functions used by the ODE system
                                "cdiff_model_for_calibration", "compute_metrics_calib"), # export the ODE function and metric function used during simulation
                          envir = environment())                               # envir=environment() means: export from the current function environment (not only from .GlobalEnv)
  parallel::clusterEvalQ(cl, library(deSolve))                                 # clusterEvalQ() runs code on each worker; here it loads deSolve so lsoda() is available if the objective calls it
  
  # Run optim in parallel (one optim per start)
  res_list <- parallel::parLapply(cl, X = seq_len(n_starts), fun = function(i) { # parLapply() runs the function on workers for each i in 1..n_starts; each worker will do one optimization run
    
    opt <- optim(par = start_points[[i]], fn = objective_fn, method = "Nelder-Mead", # optim() minimizes fn starting from par; Nelder-Mead is a derivative-free simplex method (no gradients needed)
                 control = list(maxit = maxit))                               # control=list(maxit=...) sets the maximum number of iterations for optim so it does not run forever
    
    out <- list()                                                             # create an empty list to store the results from this optimization run in a clean format
    out$par_log <- opt$par                                                    # store the best parameters found by optim on the log scale (the final parameter vector)
    out$par <- exp(opt$par)                                                   # convert the best log-parameters back to natural scale using exp() (so parameters are positive)
    out$value <- opt$value                                                    # store the best objective value (the minimized SSE returned by objective_fn)
    out$convergence <- opt$convergence                                        # store the convergence code (0 usually means success; other values mean issues)
    return(out)                                                               # return this run's result back to the master process
  })
  
  parallel::stopCluster(cl)                                                   # stopCluster() shuts down all worker processes to free CPU and memory
  
  # Pick the best result
  values <- rep(NA_real_, length(res_list))                                   # create a numeric vector filled with NA_real_ to store objective values from each run (one value per start)
  for (i in seq_along(res_list)) {                                            # seq_along(res_list) creates 1..length(res_list); we loop over all optimization outputs
    values[i] <- res_list[[i]]$value                                          # extract the objective value of run i and store it in the vector values
  }
  best_idx <- which.min(values)                                               # which.min() returns the index of the smallest value in values (best optimization run)
  
  best <- res_list[[best_idx]]                                                # extract the best run result (the one with smallest objective value)
  
  return(list(best = best, all_results = res_list))                           # return a list with (1) the best run and (2) all runs so you can inspect variability
}





###############################################################################
# 6. CALIBRATION WRAPPER (ONE FUNCTION TO RUN EVERYTHING)
###############################################################################

run_calibration <- function(initial_params, target_metrics, params_base,       # define a wrapper function that prepares inputs (log initial guess + objective function) and runs the multi-start optimization
                            init_cond, time_vec, n_starts, n_cores, maxit) {  # arguments: initial params, targets, baseline params, initial state, time vector, number of starts, number of cores, and max iterations
  
  par_names <- c("beta_h","beta_c","sigma_h","sigma_c","k_II","k_III")         # define the exact parameter names that are calibrated (must match names in initial_params and params_base)
  initial_log <- log(initial_params[par_names])                               # log() transforms the initial guesses to log-space so optim can search over all real numbers while ensuring positivity after exp()
  
  objective_fn <- create_objective_function(target_metrics, params_base, init_cond, time_vec) # build the objective function by "freezing" targets + base params + init conditions + time grid inside it
  
  res <- run_multistart_optimization(objective_fn, initial_log, n_starts, n_cores, maxit)     # run the parallel multi-start optimizations and get back the best run + all runs
  
  return(res)                                                                 # return the full results list (best + all_results)
}
