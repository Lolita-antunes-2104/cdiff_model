###############################################################################
################### 3 : CALIBRATION AND GRID SEARCH ###########################
###############################################################################
#
# This file contains all functions required to:
# 1) Explore parameter space using grid search
# 2) Calibrate the model using numerical optimization (multi-start)
#
# The calibration relies exclusively on:
# - the ODE model defined in file 1
# - the metric functions defined in file 2
#
###############################################################################

###############################################################################
# 1. GENERIC GRID SEARCH FUNCTION
###############################################################################

grid_search <- function(param_names, param_ranges, metric_function,
                        target_metrics, params_base,
                        init_cond, time_vec, n_cores = NULL) {
  
  # Number of cores (default: all minus one)
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # ---- Build parameter grid ----
  grid_list <- lapply(param_ranges, function(r) {
    seq(r[1], r[2], length.out = r[3])
  })
  names(grid_list) <- param_names
  
  param_grid <- do.call(expand.grid, grid_list)
  
  cat(sprintf(
    "Running grid search with %d parameter combinations (%d cores)\n",
    nrow(param_grid), n_cores
  ))
  
  # ---- Parallel backend ----
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    library(deSolve)
  })
  
  # ---- Export required objects to workers ----
  parallel::clusterExport(
    cl,
    c(
      # ---- Generic objects ----
      "metric_function",
      "params_base",
      "init_cond",
      "time_vec",
      
      # ---- Model runners ----
      "run_model_to_equilibrium",
      "get_equilibrium_state",
      
      # ---- Core model helpers (FILE 2 – MODEL FUNCTIONS) ----
      "compute_totals",
      "compute_lambda",
      "compute_sigmas",
      "compute_alpha",
      "add_hc_flows",
      
      # ---- Initial conditions helpers ----
      "create_initial_conditions_precalibration",
      "create_initial_conditions_from_equilibrium",
      "create_initial_conditions_vaccination",
      
      # ---- Metrics helpers (FILE 2 – METRICS FUNCTIONS) ----
      "has_vaccination",
      "get_population",
      "carriage_types",
      
      # ---- Epidemiological metrics ----
      "compute_carriage_prevalence",
      "get_CDI_count",
      "compute_CDI_rate",
      "compute_cumulative_CDI_incidence",
      "compute_recurrence_prevalence",
      "compute_dysbiosis_count",
      "compute_dysbiosis_prevalence"
    ),
    envir = environment()
  )
  
  # ---- Compute metrics for each grid point ----
  results <- parallel::parLapply(
    cl,
    seq_len(nrow(param_grid)),
    function(i) {
      row <- as.numeric(param_grid[i, ])
      names(row) <- colnames(param_grid)
      metric_function(row, params_base, init_cond, time_vec)
    }
  )
  
  results <- t(sapply(results, function(x) x))
  
  # ---- Attach metrics to grid ----
  for (m in colnames(results)) {
    param_grid[[m]] <- results[, m]
  }
  
  # ---- Distance to targets ----
  dist_sq <- rowSums(
    sapply(colnames(results), function(m) {
      (param_grid[[m]] - target_metrics[[m]])^2
    })
  )
  
  param_grid$distance <- sqrt(dist_sq)
  
  best_idx <- which.min(param_grid$distance)
  
  return(list(
    grid = param_grid,
    best_guess = param_grid[best_idx, ]
  ))
}

###############################################################################
# 2. METRIC FUNCTIONS USED FOR GRID SEARCH
###############################################################################

# --- 1) Calibration of beta_h / beta_c using carriage prevalence ---
compute_metrics_beta <- function(param_row, params_base,
                                 init_cond, time_vec) {
  
  params_tmp <- params_base
  params_tmp[c("beta_h", "beta_c")] <- param_row[c("beta_h", "beta_c")]
  
  res <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)
  last <- as.list(get_equilibrium_state(res)[1, ])
  
  return(c(
    portage_h = compute_carriage_prevalence(last, "C_tot", "h", "both"),
    portage_c = compute_carriage_prevalence(last, "C_tot", "c", "both")
  ))
}

# --- 2) Calibration of sigma_h / sigma_c using CDI incidence ---
compute_metrics_sigma <- function(param_row, params_base,
                                  init_cond, time_vec) {
  
  params_tmp <- params_base
  params_tmp[c("sigma_h", "sigma_c")] <- param_row[c("sigma_h", "sigma_c")]
  
  res <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)
  last <- as.list(get_equilibrium_state(res)[1, ])
  
  return(c(
    incidence_h = compute_CDI_rate(last, "total", "h", "both"),
    incidence_c = compute_CDI_rate(last, "total", "c", "both")
  ))
}

# --- 3) Calibration of k_II / k_III using recurrence proportions ---
compute_metrics_k <- function(param_row, params_base,
                              init_cond, time_vec) {
  
  params_tmp <- params_base
  params_tmp[c("k_II", "k_III")] <- param_row[c("k_II", "k_III")]
  
  res <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)
  last <- as.list(get_equilibrium_state(res)[1, ])
  
  return(c(
    recid_1 = compute_recurrence_prevalence(last, "rec_1", "total", "both"),
    recid_2 = compute_recurrence_prevalence(last, "rec_2", "total", "both")
  ))
}

###############################################################################
# 3. OBJECTIVE FUNCTION FOR OPTIMIZATION
###############################################################################

create_objective_function <- function(target_metrics, params_base, init_cond, time_vec,
                                      penalty = 1e10) {
  
  # Freeze objects inside the closure (important for parallel)
  force(target_metrics)
  force(params_base)
  force(init_cond)
  force(time_vec)
  force(penalty)
  
  # Expected calibrated parameter names (log-space)
  par_names <- c("beta_h", "beta_c", "sigma_h", "sigma_c", "k_II", "k_III")
  
  # Helper: safe relative error
  rel_err <- function(value, target) {
    if (!is.finite(value) || !is.finite(target) || target <= 0) return(NA_real_)
    (value - target) / target
  }
  
  # Helper: extract last equilibrium state as named numeric vector
  get_last_state_num <- function(res_eq) {
    last_df <- get_equilibrium_state(res_eq)
    x <- as.numeric(last_df[1, ])
    names(x) <- colnames(last_df)
    x
  }
  
  function(par_log) {
    
    # Map log-params to named natural-scale params
    par <- exp(par_log)
    names(par) <- par_names
    
    params_tmp <- params_base
    params_tmp[par_names] <- par
    
    out <- tryCatch({
      
      res  <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)
      last <- get_last_state_num(res)
      
      # Basic state validity checks
      if (any(!is.finite(last)) || any(last < 0)) return(penalty)
      
      # --- Metrics computed using FILE 2 functions ---
      N_h <- get_total_pop(last, setting = "h", vacc = "both")
      N_c <- get_total_pop(last, setting = "c", vacc = "both")
      
      port_h <- get_carriage(last, type = "total", setting = "h", vacc = "both")
      port_c <- get_carriage(last, type = "total", setting = "c", vacc = "both")
      
      inc_h_abs <- get_incidence(last, params_tmp, type = "total", setting = "h", vacc = "both")
      inc_c_abs <- get_incidence(last, params_tmp, type = "total", setting = "c", vacc = "both")
      
      # Convert to prevalence/rates (per capita per unit time)
      portage_h <- if (is.finite(N_h) && N_h > 0) port_h / N_h else NA_real_
      portage_c <- if (is.finite(N_c) && N_c > 0) port_c / N_c else NA_real_
      incidence_h <- if (is.finite(N_h) && N_h > 0) inc_h_abs / N_h else NA_real_
      incidence_c <- if (is.finite(N_c) && N_c > 0) inc_c_abs / N_c else NA_real_
      
      recid_1 <- compute_recurrence_prevalence(last, level = "rec_1", setting = "total", vacc = "both")
      recid_2 <- compute_recurrence_prevalence(last, level = "rec_2", setting = "total", vacc = "both")
      
      # --- Relative errors vs targets ---
      err <- c(
        rel_err(portage_h,   target_metrics$portage_h),
        rel_err(portage_c,   target_metrics$portage_c),
        rel_err(incidence_h, target_metrics$incidence_h),
        rel_err(incidence_c, target_metrics$incidence_c),
        rel_err(recid_1,     target_metrics$recid_1),
        rel_err(recid_2,     target_metrics$recid_2)
      )
      
      if (any(!is.finite(err))) return(penalty)
      sum(err^2)
      
    }, error = function(e) penalty)
    
    out
  }
}

###############################################################################
# 4. MULTI-START OPTIMIZATION
###############################################################################

run_multistart_optimization <- function(objective_fn,
                                        initial_log,
                                        n_starts = 10,
                                        seed = 123,
                                        n_cores = NULL,
                                        jitter = 0.3,
                                        maxit = 10000,
                                        verbose = TRUE) {
  
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  n_cores <- max(1, n_cores)
  
  set.seed(seed)
  
  # Build start points in log-space
  start_points <- lapply(seq_len(n_starts), function(i) {
    if (i == 1) initial_log else initial_log + runif(length(initial_log), -jitter, jitter)
  })
  
  # Single-core fallback (simpler + no export headaches)
  if (n_cores == 1) {
    if (verbose) cat(sprintf("Multi-start optimization: %d starts (1 core)\n", n_starts))
    
    one_run <- function(i) {
      opt <- optim(start_points[[i]], objective_fn, method = "Nelder-Mead",
                   control = list(maxit = maxit))
      list(par_log = opt$par, par = exp(opt$par), value = opt$value, convergence = opt$convergence)
    }
    
    res <- if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply::pblapply(seq_len(n_starts), one_run)
    } else {
      lapply(seq_len(n_starts), one_run)
    }
    
  } else {
    if (verbose) cat(sprintf("Multi-start optimization: %d starts (%d cores)\n", n_starts, n_cores))
    
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterEvalQ(cl, { library(deSolve); NULL })
    
    # Export everything needed by objective_fn on workers
    export_names <- unique(c(
      # Local objects
      "objective_fn", "start_points", "maxit",
      # File 1
      "run_model_to_equilibrium", "get_equilibrium_state",
      "cdiff_core", "cdiff_hc_model", "cdiff_hc_vacc_model",
      # File 2 helpers used by ODE / metrics
      "compute_totals", "compute_lambda", "compute_sigmas", "compute_alpha", "add_hc_flows",
      "has_vaccination", "get_total_pop", "get_carriage", "get_incidence", "get_cumulative_incidence",
      "compute_recurrence_prevalence", "get_dysbiosis"
    ))
    
    # Export only objects that exist (prevents clusterExport errors)
    export_names <- export_names[vapply(
      export_names,
      function(nm) exists(nm, envir = environment(), inherits = TRUE),
      logical(1)
    )]
    
    parallel::clusterExport(cl, export_names, envir = environment())
    
    one_run <- function(i) {
      opt <- optim(start_points[[i]], objective_fn, method = "Nelder-Mead",
                   control = list(maxit = maxit))
      list(par_log = opt$par, par = exp(opt$par), value = opt$value, convergence = opt$convergence)
    }
    
    res <- if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply::pblapply(seq_len(n_starts), one_run, cl = cl)
    } else {
      parallel::parLapply(cl, seq_len(n_starts), one_run)
    }
  }
  
  best_idx <- which.min(vapply(res, function(x) x$value, numeric(1)))
  
  list(
    best = res[[best_idx]],
    all_results = res
  )
}

###############################################################################
# 5. CALIBRATION WRAPPER
###############################################################################

run_calibration <- function(initial_params,
                            target_metrics,
                            params_base,
                            init_cond,
                            time_vec,
                            n_starts = 10,
                            seed = 123,
                            n_cores = NULL) {
  
  initial_log <- log(initial_params)
  
  obj_fn <- create_objective_function(
    target_metrics = target_metrics,
    params_base    = params_base,
    init_cond      = init_cond,
    time_vec       = time_vec
  )
  
  run_multistart_optimization(
    objective_fn = obj_fn,
    initial_log  = initial_log,
    n_starts     = n_starts,
    seed         = seed,
    n_cores      = n_cores
  )
}



