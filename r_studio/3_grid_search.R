###############################################################################
############################## 3 : GRID SEARCH ################################
###############################################################################

###############################################################################
# ---- GENERIC grid search function ----
###############################################################################
grid_search <- function(param_names, param_ranges, metric_function,
                        target_metrics, params_base, init_cond, time_vec, n_cores = NULL) {
  
  # Si l’utilisateur ne précise pas le nombre de coeurs, on prend (nb_coeurs - 1) pour éviter de saturer la machine
  # max(1, ...) garantit qu’on a au moins 1 coeur (utile sur petites machines / environnements limités)
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  
  # ---- Construction de la grille de paramètres ----
  # Pour chaque paramètre, on crée une séquence (min -> max) avec un certain nombre de points (range[3])
  grid_list <- lapply(param_ranges, function(range) {
    seq(range[1], range[2], length.out = range[3])
  })
  
  # On nomme chaque vecteur de valeurs avec le nom du paramètre correspondant (beta_h, beta_c, etc.)
  names(grid_list) <- param_names
  
  # expand.grid fabrique toutes les combinaisons possibles de la grille (cartésien)
  # ex : si 10 valeurs pour beta_h et 10 pour beta_c => 100 combinaisons
  param_grid <- do.call(expand.grid, grid_list)
  
  # Message informatif : combien de simulations on va lancer et sur combien de coeurs
  cat(sprintf("Computing metrics for %d parameter combinations (parallel: %d cores)...\n",
              nrow(param_grid), n_cores))
  
  # ---- Initialisation du calcul parallèle ----
  # On crée un cluster de workers (process R séparés) pour paralléliser les simulations
  cl <- parallel::makeCluster(n_cores)
  
  # Sécurité : même si la fonction plante, on arrête proprement le cluster à la fin
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # Sur chaque worker : on charge deSolve car lsoda est utilisé dans le modèle
  parallel::clusterEvalQ(cl, {
    library(deSolve)
  })
  
  # On “envoie” aux workers tout ce dont ils ont besoin :
  # - la grille et les objets passés à metric_function
  # - le modèle ODE et ses helpers
  # - les fonctions de métriques (car metric_function les appelle)
  # Important : envir = environment() permet d’exporter depuis l’environnement de grid_search (pas juste le global)
  parallel::clusterExport(
    cl,
    c(
      # Objets liés à grid_search / au calcul des métriques
      "metric_function", "params_base", "init_cond", "time_vec", "param_grid",
      
      # Modèle (simulation ODE jusqu'à l'équilibre)
      "run_model_to_equilibrium", "get_equilibrium_state", "cdiff_micro",
      "compute_totals", "compute_lambda", "compute_sigmas",
      
      # Helpers des métriques
      "get_param", "safe_div",
      
      # Fonctions de métriques exportées (utilisées par compute_all_metrics ou metric_function)
      "compute_population_totals", "compute_carriage_prevalence",
      "compute_CDI_incidence", "compute_recurrence_prevalence",
      "compute_R0_approx", "compute_cumulative_CDI_incidence",
      "compute_all_metrics"
    ),
    envir = environment()
  )
  # ---- Fin du setup parallèle ----
  
  # ---- Calcul des métriques pour chaque combinaison de la grille ----
  # parLapply répartit les indices i sur les workers
  results <- parallel::parLapply(cl, seq_len(nrow(param_grid)), function(i) {
    
    # On récupère la i-ème ligne de la grille (1 combinaison de paramètres)
    # as.numeric : on veut un vecteur numérique “propre”
    row <- as.numeric(param_grid[i, , drop = TRUE])
    
    # On redonne les noms (beta_h, beta_c, ...) au vecteur, pour pouvoir faire param_row["beta_h"]
    names(row) <- colnames(param_grid)
    
    # On appelle la fonction métrique spécifique (beta / sigma / k...)
    # Elle renvoie un petit vecteur nommé (ex: portage_h, portage_c)
    metric_function(row, params_base, init_cond, time_vec)
  })
  
  # ---- Mise en forme des résultats ----
  # results est une liste de vecteurs (1 vecteur par point de la grille)
  # sapply(...): empile en matrice, puis t() pour avoir 1 ligne = 1 combinaison de paramètres
  results <- t(sapply(results, function(x) x))
  
  # ---- Ajout des métriques à la grille ----
  # Pour chaque colonne métrique (portage_h, incidence_h, etc.), on ajoute la colonne au data.frame param_grid
  for (col_name in colnames(results)) {
    param_grid[[col_name]] <- results[, col_name]
  }
  
  # ---- Calcul d’une distance aux cibles ----
  # On prend les colonnes métriques et on calcule une distance euclidienne aux valeurs cibles
  metric_cols <- colnames(results)
  
  # distance_sq = somme des (écarts)^2 sur les métriques utilisées
  distance_sq <- rowSums(sapply(seq_along(metric_cols), function(i) {
    (param_grid[[metric_cols[i]]] - target_metrics[[metric_cols[i]]])^2
  }))
  
  # Distance finale = racine carrée (distance euclidienne)
  param_grid$distance <- sqrt(distance_sq)
  
  # On prend le meilleur point = celui avec distance minimale
  best_idx <- which.min(param_grid$distance)
  best_guess <- param_grid[best_idx, ]
  
  # On renvoie la grille complète (utile pour plots) + le “meilleur” point de la grille
  return(list(grid = param_grid, best_guess = best_guess))
}

###############################################################################
# ---- metric functions for grid search ----
###############################################################################

# ---- 1) Métriques pour calibrer beta_h et beta_c via le portage ----
compute_metrics_beta <- function(param_row, params_base, init_cond, time_vec) { 
  # On copie les paramètres de base (sinon on écraserait params_base)
  params_temp <- params_base
  
  # On remplace les 2 paramètres testés par la paire (beta_h, beta_c) de la grille
  params_temp["beta_h"] <- param_row["beta_h"]
  params_temp["beta_c"] <- param_row["beta_c"]
  
  # On simule le modèle jusqu’à équilibre (ou pseudo-équilibre selon time_vec)
  result <- run_model_to_equilibrium(params_temp, init_cond, time_vec)
  
  # On récupère l’état final (dernière ligne)
  last <- get_equilibrium_state(result)
  
  # On convertit la ligne en liste pour pouvoir accéder avec last$C0_h, last$I_h, etc.
  last <- as.list(last[1, ])
  
  # On renvoie les métriques de portage (prévalence colonisés / N) dans chaque setting
  c(
    portage_h = unname(compute_carriage_prevalence(last, "h")),
    portage_c = unname(compute_carriage_prevalence(last, "c"))
  )
}

# ---- 2) Métriques pour calibrer sigma_h et sigma_c via l’incidence ----
compute_metrics_sigma <- function(param_row, params_base, init_cond, time_vec) {
  # Copie des paramètres de base
  params_temp <- params_base
  
  # On remplace sigma_h et sigma_c par les valeurs testées
  params_temp["sigma_h"] <- param_row["sigma_h"]
  params_temp["sigma_c"] <- param_row["sigma_c"]
  
  # Simulation jusqu’à l’équilibre
  result <- run_model_to_equilibrium(params_temp, init_cond, time_vec)
  
  # Etat final
  last <- get_equilibrium_state(result)
  last <- as.list(last[1, ])
  
  # Populations totales (pour normaliser l'incidence en “par personne” / jour)
  N_h <- compute_population_totals(last, "h")
  N_c <- compute_population_totals(last, "c")
  
  # Incidence CDI totale (primo + récidives) / N, dans chaque setting
  inc_h <- compute_CDI_incidence(last, params_temp, "h", "total") / (N_h+N_c)
  inc_c <- compute_CDI_incidence(last, params_temp, "c", "total") / (N_h+N_c) 
  
  # On renvoie les deux métriques cibles pour la grille
  c(
    incidence_h = unname(inc_h),
    incidence_c = unname(inc_c)
  )
}

# ---- 3) Métriques pour calibrer k_II et k_III via les proportions de récidives ----
compute_metrics_k <- function(param_row, params_base, init_cond, time_vec) {
  # Copie des paramètres de base
  params_temp <- params_base
  
  # On remplace k_II et k_III par les valeurs testées
  params_temp["k_II"] <- param_row["k_II"]
  params_temp["k_III"] <- param_row["k_III"]
  
  # Simulation jusqu’à l’équilibre
  result <- run_model_to_equilibrium(params_temp, init_cond, time_vec)
  
  # Etat final
  last <- get_equilibrium_state(result)
  last <- as.list(last[1, ])
  
  # On renvoie :
  # recid_1 = proportion I_II / I (1ère récidive parmi primo)
  # recid_2 = proportion I_III / I_II (2e+ récidive parmi 1ères récidives)
  c(
    recid_1 = unname(compute_recurrence_prevalence(last, "both", "rec_1")),
    recid_2 = unname(compute_recurrence_prevalence(last, "both", "rec_2"))
  )
}







