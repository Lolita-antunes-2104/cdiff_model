###############################################################################
############################# 4 : CALIBRATION #################################
###############################################################################

###############################################################################
# ---- objective function ----
###############################################################################
create_objective_function <- function(target_metrics, params_base, init_cond, time_vec) {
  
  force(target_metrics)
  force(params_base)
  force(init_cond)
  force(time_vec)
  
  function(par_log) {
    
    par <- exp(par_log)
    names(par) <- c("beta_h", "beta_c", "sigma_h", "sigma_c", "k_II", "k_III")
    
    params_temp <- params_base
    params_temp[names(par)] <- par
    
    tryCatch({
      
      result <- run_model_to_equilibrium(params_temp, init_cond, time_vec)
      
      last <- get_equilibrium_state(result)
      last <- as.list(last[1, ])
      
      if (any(is.na(unlist(last))) || any(unlist(last) < 0)) return(1e10)
      
      metrics <- compute_all_metrics(
        last_state = last, 
        params_vec = params_temp, 
        ode_result = NULL,
        targets = target_metrics
      )
      
      sum(unlist(metrics$errors)^2)
      
    }, error = function(e) {
      cat("Erreur dans objectif:", conditionMessage(e), "\n")  # Debug
      1e10
    })
  }
}


###############################################################################
# ---- multi-start optimization ----
###############################################################################
run_multistart_optimization <- function(objective_fn, initial_log,          # Lance plusieurs optimisations depuis plusieurs points de départ
                                        n_starts = 10, seed = 123, n_cores = NULL) {
  
  # n_cores robuste
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)     # Par défaut : tous les coeurs moins 1, au minimum 1
  
  set.seed(seed)                                                          # Fixe le seed pour reproductibilité des points de départ aléatoires
  cat(sprintf("\n=== Multi-start optimization (parallel: %d cores) ===\n", n_cores))  # Message console
  
  # Start points (log-espace)
  start_points <- lapply(seq_len(n_starts), function(i) {                  # Génère n_starts points initiaux en log-espace
    if (i == 1) initial_log else initial_log + runif(length(initial_log), -0.3, 0.3) 
    # Le 1er run démarre exactement à initial_log
    # Les autres : petite perturbation uniforme [-0.3, +0.3] sur chaque paramètre log
  })
  
  # Cluster
  cl <- parallel::makeCluster(n_cores)                                     # Crée un cluster de workers (process R séparés)
  on.exit(parallel::stopCluster(cl), add = TRUE)                           # Arrête le cluster à la fin, même si erreur
  
  parallel::clusterEvalQ(cl, {                                             # Sur chaque worker : charge deSolve (lsoda)
    library(deSolve)
  })
  
  # Exporte l'objectif + start_points (environnement local)
  parallel::clusterExport(
    cl,
    c("objective_fn", "start_points"),                                     # Envoie la fonction objectif et les points de départ aux workers
    envir = environment()                                                  # Export depuis l’environnement local de run_multistart_optimization
  )
  
  # Exporte le modèle + métriques (global)
  export_vars <- c(
    "run_model_to_equilibrium", "get_equilibrium_state",
    "compute_all_metrics",
    "cdiff_micro", "compute_totals", "compute_lambda", "compute_sigmas",
    "compute_carriage_prevalence", "compute_CDI_incidence",
    "compute_recurrence_prevalence", "compute_population_totals",
    "compute_R0_approx", "compute_cumulative_CDI_incidence",
    "compute_alpha_eq",
    "get_param", "safe_div",
    "params"
  )
  
  export_vars <- export_vars[export_vars %in% ls(.GlobalEnv)]               # Ne garde que ceux qui existent vraiment dans l’environnement global
  
  parallel::clusterExport(cl, export_vars, envir = .GlobalEnv)             # Export des fonctions/objets depuis le GlobalEnv vers les workers
  
  # Worker function
  one_run <- function(i) {                                                 # Fonction exécutée pour un départ i (1 optimisation)
    opt <- optim(                                                          # Lance optim (minimisation)
      start_points[[i]], objective_fn,                                      # Point initial + fonction objectif
      method  = "Nelder-Mead",                                             # Méthode robuste sans gradient (simplex)
      control = list(maxit = 10000, reltol = 1e-8)                          # Paramètres d'arrêt (beaucoup d'itérations, tolérance stricte)
    )
    
    opt_natural <- exp(opt$par)                                             # Convertit les paramètres optimisés (log) en échelle naturelle
    names(opt_natural) <- c("beta_h", "beta_c", "sigma_h", "sigma_c", "k_II", "k_III")  # Remet les noms
    
    list(                                                                  # Renvoie un objet structuré pour ce run
      par_log      = opt$par,                                              # Paramètres optimisés en log
      par_natural  = opt_natural,                                          # Paramètres optimisés en échelle naturelle
      value        = opt$value,                                            # Valeur finale de l'objectif (SSE des erreurs)
      convergence  = opt$convergence                                       # Code convergence (0 = OK en général)
    )
  }
  
  # Avec barre de progression si pbapply dispo
  if (requireNamespace("pbapply", quietly = TRUE)) {                        # Si pbapply installé, on affiche une barre de progression
    pbapply::pboptions(type = "txt")                                        # Barre de progression en mode texte
    all_results <- pbapply::pblapply(seq_len(n_starts), one_run, cl = cl)   # Lance n_starts runs avec progress bar (en utilisant le cluster)
  } else {
    all_results <- parallel::parLapply(cl, seq_len(n_starts), one_run)      # Sinon : parLapply classique (sans barre)
  }
  
  # Meilleur run
  best_idx <- which.min(vapply(all_results, function(x) x$value, numeric(1))) # Cherche l’index du run avec la plus petite valeur objectif
  
  list(
    best = all_results[[best_idx]],                                         # Le meilleur run
    all_results = all_results                                               # Tous les runs (utile pour diagnostic)
  )
}

###############################################################################
# ---- run calibration ----
###############################################################################
run_calibration <- function(initial_params, target_metrics, params_base,     # Wrapper “simple” : construit l’objectif + lance multi-start
                            init_cond, time_vec, n_starts = 10) {
  
  initial_log <- log(initial_params)                                         # Passe les paramètres initiaux en log (pour optimiser en positif)
  
  obj_fn <- create_objective_function(target_metrics, params_base, init_cond, time_vec) 
  # Construit la fonction objectif (fermeture) avec les bons objets capturés
  
  run_multistart_optimization(obj_fn, initial_log, n_starts)                 # Lance l’optimisation multi-start
}

###############################################################################
# ---- compute and display calibrated results ----
###############################################################################
compute_calibrated_metrics <- function(calibrated_params, params_base,       # Après calibration : simule et renvoie params + métriques + trajectoires ODE
                                       init_cond, time_vec, targets) {
  
  params_final <- params_base                                                # Copie des paramètres de base
  params_final[names(calibrated_params)] <- calibrated_params                # Remplace par les paramètres calibrés (beta/sigma/k)
  
  result <- run_model_to_equilibrium(params_final, init_cond, time_vec)      # Relance une simulation complète avec les params calibrés
  
  last <- get_equilibrium_state(result)                                      # Etat final
  last <- as.list(last[1, ])                                                 # format cohérent (liste nommée)
  
  # pas de return(1e10) ici (sinon tu casses le type de retour)
  if (any(is.na(unlist(last))) || any(unlist(last) < 0)) {                   # Vérifie rapidement si la solution est invalide
    warning("Solution finale invalide (NA ou négatif) dans compute_calibrated_metrics().") # Avertissement (mais on continue)
  }
  
  metrics <- compute_all_metrics(last, params_final, targets = targets)      # Calcule toutes les métriques + erreurs vs targets
  
  list(
    params_final = params_final,                                             # Paramètres finaux complets (base + calibrés)
    metrics = metrics,                                                       # Liste de métriques calculées
    ode_result = result                                                      # Trajectoire ODE complète (utile pour plots)
  )
}




