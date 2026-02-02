###############################################################################
# 5. RUN UNTIL EQUILIBRIUM WITH CALIBRATED PARAMETERS TO FIND INITIAL CONDITIONS
###############################################################################

# Run the calibration model with the final calibrated parameters (stop early if equilibrium)
out_eq <- run_calib_model_to_equilibrium(
  params   = params_calib,   # calibrated parameters
  init_cond = init_cond,     # initial conditions
  time_vec = time_vec_eq,        # long horizon, but will stop early if equilibrium reached
  threshold = threshold_calib)

# Get the last state (equilibrium state) as a named numeric vector : use as initial conditions for scenario runs !!
state_eq <- get_last_state(out_eq)

# Dynamics plots (hospital + community side by side)
dyn_plots <- plot_dynamics(
  ode_result = out_eq,
  targets = list(prevalence_h = target_metrics$prevalence_h, prevalence_c = target_metrics$prevalence_c),
  N_h = N_h,
  N_c = N_c)

print(dyn_plots$both)  # side-by-side plot

# Alpha plot
p_alpha <- plot_alpha_dynamics(out_eq, params_calib)
print(p_alpha)

# SAVE EQUILIBRIUM RUN OUTPUTS AS RDS (Bundle everything in a single object)
res_eq <- list(
  out_eq = out_eq,
  state_eq = state_eq,
  params_calib = params_calib,
  threshold_calib = threshold_calib,
  target_metrics = target_metrics,
  N_h = N_h,
  N_c = N_c
)
saveRDS(res_eq, "results_equilibrium.rds")