###############################################################################
################### 3 : CALIBRATION AND GRID SEARCH ###########################
###############################################################################
#
# This file contains the functions required to:
# 1) Explore parameter space using grid search (parallel)
# 2) Provide metric functions used during calibration steps
#
# Dependencies:
# - ODE model + solvers from file 1 (run_model_to_equilibrium, get_equilibrium_state, cdiff_* models)
# - Helper + metric functions from file 2 (compute_* helpers, get_* metrics)
#
###############################################################################

###############################################################################
# 1. GENERIC GRID SEARCH FUNCTION
###############################################################################

grid_search <- function(param_names,
                        param_ranges,
                        metric_function,
                        target_metrics,
                        params_base,
                        init_cond,
                        time_vec,
                        n_cores = NULL,
                        extra_exports = character(0),
                        verbose = TRUE) {
  
  # ---- Input checks ----
  if (length(param_names) == 0) stop("param_names is empty.")
  if (length(param_names) != length(param_ranges)) stop("param_names and param_ranges must have the same length.")
  if (!is.list(target_metrics) || is.null(names(target_metrics)) || any(names(target_metrics) == "")) {
    stop("target_metrics must be a *named* list (e.g., list(portage_h=..., portage_c=...)).")
  }
  
  metric_names <- names(target_metrics)
  
  # ---- Number of cores (default: all minus one) ----
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  
  # ---- Build parameter grid ----
  grid_list <- lapply(param_ranges, function(r) {
    if (length(r) != 3) stop("Each element of param_ranges must be c(min, max, n_points).")
    seq(r[1], r[2], length.out = r[3])
  })
  names(grid_list) <- param_names
  
  param_grid <- expand.grid(grid_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  if (verbose) {
    cat(sprintf("Running grid search with %d combinations (%d cores)\n", nrow(param_grid), n_cores))
  }
  
  # ---- Parallel backend ----
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    library(deSolve)
    NULL
  })
  
  # ---- Export required objects to workers ----
  export_names <- unique(c(
    # Generic objects
    "metric_function", "params_base", "init_cond", "time_vec",
    
    # File 1
    "run_model_to_equilibrium", "get_equilibrium_state",
    "cdiff_core", "cdiff_hc_model", "cdiff_hc_vacc_model",
    
    # File 2 - model helpers
    "compute_totals", "compute_lambda", "compute_sigmas", "compute_alpha", "add_hc_flows",
    "create_initial_conditions_precalibration",
    "create_initial_conditions_from_equilibrium",
    "create_initial_conditions_vaccination",
    
    # File 2 - metrics
    "has_vaccination", "get_total_pop", "get_carriage", "get_incidence", "get_cumulative_incidence",
    "compute_recurrence_prevalence", "get_dysbiosis",
    
    # Local helper used by compute_metrics_*
    "get_last_state",
    
    # Optional extras
    extra_exports
  ))
  
  # Export only objects that exist (robust to small name changes)
  export_env <- environment()  # <- IMPORTANT: contient params_base/init_cond/time_vec/metric_function
  export_names <- export_names[vapply(
    export_names,
    function(nm) exists(nm, envir = export_env, inherits = TRUE),
    logical(1)
  )]
  parallel::clusterExport(cl, export_names, envir = export_env)
  
  # ---- Compute metrics for each grid point ----
  results_list <- parallel::parLapply(
    cl,
    X = seq_len(nrow(param_grid)),
    fun = function(i) {
      row <- as.numeric(param_grid[i, param_names, drop = FALSE])
      names(row) <- param_names
      
      out <- tryCatch(
        metric_function(row, params_base, init_cond, time_vec),
        error = function(e) NULL
      )
      
      # Ensure a consistent named numeric vector output
      if (is.null(out)) {
        out <- setNames(rep(NA_real_, length(metric_names)), metric_names)
      } else {
        out <- out[metric_names]
        out <- suppressWarnings(as.numeric(out))
        names(out) <- metric_names
      }
      
      out
    }
  )
  
  results_mat <- do.call(rbind, results_list)
  results_df  <- as.data.frame(results_mat, stringsAsFactors = FALSE)
  
  # ---- Attach metrics to grid ----
  out_grid <- cbind(param_grid, results_df)
  
  # ---- Distance to targets (NA-safe: NA -> Inf) ----
  target_vec <- unlist(target_metrics[metric_names])
  diff_mat <- as.matrix(out_grid[, metric_names, drop = FALSE]) - matrix(target_vec, nrow = nrow(out_grid), ncol = length(target_vec), byrow = TRUE)
  dist_sq <- rowSums(diff_mat^2)
  dist_sq[!is.finite(dist_sq)] <- Inf
  
  out_grid$distance <- sqrt(dist_sq)
  
  best_idx <- which.min(out_grid$distance)
  
  list(
    grid = out_grid,
    best_guess = out_grid[best_idx, , drop = FALSE]
  )
}

###############################################################################
# 2. METRIC FUNCTIONS USED FOR GRID SEARCH
###############################################################################

# --- 1) Calibration of beta_h / beta_c using carriage prevalence ---
compute_metrics_beta <- function(param_row, params_base, init_cond, time_vec) {
  
  params_tmp <- params_base
  params_tmp[c("beta_h", "beta_c")] <- param_row[c("beta_h", "beta_c")]
  
  res  <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)
  last <- get_last_state(res)
  
  N_h <- get_total_pop(last, setting = "h", vacc = "both")
  N_c <- get_total_pop(last, setting = "c", vacc = "both")
  
  port_h <- get_carriage(last, type = "total", setting = "h", vacc = "both")
  port_c <- get_carriage(last, type = "total", setting = "c", vacc = "both")
  
  c(
    portage_h = if (is.finite(N_h) && N_h > 0) port_h / N_h else NA_real_,
    portage_c = if (is.finite(N_c) && N_c > 0) port_c / N_c else NA_real_
  )
}

# --- 2) Calibration of sigma_h / sigma_c using CDI incidence rate ---
# Returns an incidence *rate* (incidence / population), consistent with your older targets like 18.8/100000, etc.
# --- 2) Calibration of sigma_h / sigma_c using CDI incidence rate ---
# Returns incidence rates per TOTAL population (N_tot) to match targets like 18/100k/year.
compute_metrics_sigma <- function(param_row, params_base, init_cond, time_vec) {
  
  params_tmp <- params_base
  params_tmp[c("sigma_h", "sigma_c")] <- param_row[c("sigma_h", "sigma_c")]
  
  res  <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)
  last <- get_last_state(res)
  
  N_h <- get_total_pop(last, setting = "h", vacc = "both")
  N_c <- get_total_pop(last, setting = "c", vacc = "both")
  N_tot <- N_h + N_c
  
  inc_h <- get_incidence(last, params_tmp, type = "total", setting = "h", vacc = "both")
  inc_c <- get_incidence(last, params_tmp, type = "total", setting = "c", vacc = "both")
  
  c(
    incidence_h = if (is.finite(N_tot) && N_tot > 0) inc_h / N_tot else NA_real_,
    incidence_c = if (is.finite(N_tot) && N_tot > 0) inc_c / N_tot else NA_real_
  )
}


compute_metrics_k <- function(param_row, params_base, init_cond, time_vec) {
  
  params_tmp <- params_base
  params_tmp[c("k_II", "k_III")] <- param_row[c("k_II", "k_III")]
  
  res  <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)
  last <- get_last_state(res)
  
  c(
    recid_1 = compute_recurrence_prevalence(last, "rec_1", "total", "both"),
    recid_2 = compute_recurrence_prevalence(last, "rec_2", "total", "both")
  )
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
  
  function(par_log) {
    
    # Map log-params to named natural-scale params
    par <- exp(par_log)
    names(par) <- par_names
    
    params_tmp <- params_base
    params_tmp[par_names] <- par
    
    out <- tryCatch({
      
      res  <- run_model_to_equilibrium(params_tmp, init_cond, time_vec)
      last <- get_last_state(res)
      
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
      # Incidence targets are per TOTAL population
      N_tot <- N_h + N_c
      incidence_h <- if (is.finite(N_tot) && N_tot > 0) inc_h_abs / N_tot else NA_real_
      incidence_c <- if (is.finite(N_tot) && N_tot > 0) inc_c_abs / N_tot else NA_real_
      
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
      "compute_recurrence_prevalence", "get_dysbiosis",
      # File 3 helper
      "get_last_state"
    ))
    
    # Export only objects that exist (robust to small name changes)
    export_names <- export_names[vapply(
      export_names,
      function(nm) exists(nm, envir = .GlobalEnv, inherits = TRUE),
      logical(1)
    )]
    
    parallel::clusterExport(cl, export_names, envir = .GlobalEnv)
    
    
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

run_calibration <- function(initial_params, target_metrics, params_base, init_cond, time_vec,
                            n_starts = 10, seed = 123, n_cores = NULL) {
  
  par_names <- c("beta_h","beta_c","sigma_h","sigma_c","k_II","k_III")
  stopifnot(all(par_names %in% names(initial_params)))
  initial_log <- log(initial_params[par_names])
  
  obj_fn <- create_objective_function(target_metrics, params_base, init_cond, time_vec)
  
  run_multistart_optimization(obj_fn, initial_log, n_starts = n_starts, seed = seed, n_cores = n_cores)
}





