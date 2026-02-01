plot_dynamics <- function(ode_result, targets, N_h, N_c) {
  
  # Add aggregated compartments for visualization
  ode_result$C_h <- with(ode_result, C0_h + CA_h + C_II_h + C_III_h)
  ode_result$C_c <- with(ode_result, C0_c + CA_c + C_II_c + C_III_c)
  ode_result$I_tot_h <- with(ode_result, I_h + I_II_h + I_III_h)
  ode_result$I_tot_c <- with(ode_result, I_c + I_II_c + I_III_c)
  ode_result$R_h <- with(ode_result, S_II_h + S_III_h + C_II_h + C_III_h + I_II_h + I_III_h)
  ode_result$R_c <- with(ode_result, S_II_c + S_III_c + C_II_c + C_III_c + I_II_c + I_III_c)
  
  # Convert targets from prevalence to absolute numbers
  target_C_h <- targets$portage_h * N_h
  target_C_c <- targets$portage_c * N_c
  
  # ---- Hospital - all compartments ----
  p_h_all <- ode_result %>%
    select(time, S0_h, SA_h, C0_h, CA_h, I_h,
           S_II_h, C_II_h, I_II_h, S_III_h, C_III_h, I_III_h) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot(aes(time, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    theme_bw() +
    labs(x = "Time (days)", 
         y = "Number of individuals",
         color = "Compartment",
         title = "Model dynamics - Hospital (all compartments)")
  
  # ---- Hospital - aggregated ----
  p_h_agg <- ode_result %>%
    select(time, S0_h, SA_h, C_h, I_tot_h, R_h) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot(aes(time, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = target_C_h, linetype = "dashed", color = "red") +
    annotate("text", 
             x = max(ode_result$time) * 0.8, 
             y = target_C_h * 1.1,
             label = "Target carriage", 
             color = "red", 
             size = 3) +
    theme_bw() +
    labs(x = "Time (days)", 
         y = "Number of individuals",
         color = "Compartment",
         title = "Model dynamics - Hospital (aggregated)")
  
  # ---- Community - all compartments ----
  p_c_all <- ode_result %>%
    select(time, S0_c, SA_c, C0_c, CA_c, I_c,
           S_II_c, C_II_c, I_II_c, S_III_c, C_III_c, I_III_c) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot(aes(time, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    theme_bw() +
    labs(x = "Time (days)", 
         y = "Number of individuals",
         color = "Compartment",
         title = "Model dynamics - Community (all compartments)")
  
  # ---- Community - aggregated ----
  p_c_agg <- ode_result %>%
    select(time, S0_c, SA_c, C_c, I_tot_c, R_c) %>%
    reshape2::melt(id.vars = "time") %>%
    ggplot(aes(time, value, color = variable)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = target_C_c, linetype = "dashed", color = "red") +
    annotate("text", 
             x = max(ode_result$time) * 0.8, 
             y = target_C_c * 1.1,
             label = "Target carriage", 
             color = "red", 
             size = 3) +
    theme_bw() +
    labs(x = "Time (days)", 
         y = "Number of individuals",
         color = "Compartment",
         title = "Model dynamics - Community (aggregated)")
  
  # Return list of plots
  return(list(
    hospital_all = p_h_all,
    hospital_agg = p_h_agg,
    community_all = p_c_all,
    community_agg = p_c_agg
  ))
}