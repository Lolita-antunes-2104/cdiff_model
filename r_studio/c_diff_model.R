#package ----
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

#c. difficile transmission model ----
cdiff_model <- function(t, pop, params) {
  
  with(as.list(c(pop, params)), {
    
    #total pop 
    N <- S0 + SA + C0 + CA + I + S_II + C_II + I_II + S_III + C_III + I_III
    Sh <- S0 + SA + S_II + S_III
    Ch <- C0 + CA + C_II + C_III 
    Ih <- I + I_II + I_III
    
    #infection force 
    lambda <- beta * (Ch + nu*Ih)/N
    
    #edo system
    dS0 <- -lambda*S0 - tau*S0 + gamma*C0 + phi*S_II + phi*S_III #+ u0 - u1*S0 #attention communauté : pas admission mais naissance donc u0 uniquement dans S0 !!!
    dSA <- -lambda*SA + gamma*CA + tau*S0 #+ u0 - u1*SA #différent u0 pour le prochain code 
    dC0 <- lambda*S0 - gamma*C0 - tau*C0 - sigma*C0 #+ u0 - u1*C0
    dCA <- lambda*SA + tau*C0 - gamma*CA - sigma_A*CA #+ u0 - u1*CA
    dI <- sigma*C0  + sigma_A*CA - epsilon*I #+ u0 - mu*I
    dS_II <- p*epsilon*I + gamma*C_II - lambda*S_II - phi*S_II #+ u0 - u1*S_II
    dC_II <- (1 - p)*epsilon*I + lambda*S_II - gamma*C_II - sigma_II*C_II #+ u0 - u1*C_II
    dI_II <- sigma_II*C_II - epsilon*I_II #+ u0 - mu_II*I_II
    dS_III <- p*epsilon*I_II + p*epsilon*I_III + gamma*C_III - lambda*S_III - phi*S_III #+ u0 - u1*S_III
    dC_III <- (1-p)*epsilon*I_II + (1-p)*epsilon*I_III + lambda*S_III - gamma*C_III - sigma_III*C_III #+ u0 - u1*C_III
    dI_III <- sigma_III*C_III - epsilon*I_III #+ u0 - mu_III*I_III

    list(c(dS0, dSA, dC0, dCA, dI, dS_II, dC_II, dI_II, dS_III, dC_III, dI_III))
  })
} 

# parameters ----
beta=0.021
nu=36
gamma=0.03
sigma=0.000445
sigma_A=sigma*6.67
sigma_II=sigma
sigma_III=sigma
tau=0.044
p=0.5
epsilon=0.07
phi=0.018
#mu=0.0036
#mu_II=1.88*0.0036
#mu_III=2.38*0.0036
#u0=0
#u1=0
params <- c(beta=beta, nu=nu, gamma=gamma, sigma=sigma, sigma_A=sigma_A, 
            sigma_II=sigma_II, sigma_III=sigma_III, 
            tau=tau, p=p, epsilon=epsilon, phi=phi) #vecteur nommé pour lsoda

# time for simulation ----
dt=1
Tmax=140
times <- seq(from=0, to=Tmax, by=dt)

# initial conditions ----
N = 1000
S0 = N*0.66 # the majority
SA = N*0.12 # % of people taking ATB /2
C0 = N*0.03# % of people colonised /2
CA = N*0.15# ATB+colonised
I = N*0.01
S_II = N*0.005
C_II = N*0.005
I_II = N*0.005
S_III = N*0.005
C_III = N*0.005
I_III = N*0.005
init.cond <- c(S0=S0, SA=SA, C0=C0, CA=CA, I=I,
               S_II=S_II, C_II=C_II, I_II=I_II, 
               S_III=S_III, C_III=C_III,I_III=I_III) 

# integration without calibration ----
edo_without_cal <- as.data.frame(lsoda(init.cond, times, cdiff_model, params))
# View(edo_without_cal)

# plot ----
edo_without_cal %>%
  melt(id = "time") %>%
  #dplyr::filter(variable %in% c("I","I_II","I_III")) %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Value", colour = "Population:")


# calibration function with least square ----


#APRES AJUSTEMENT DES PARAMETRES 
# simulating the endemic equilibrium for initials conditions ----



#fonction pour extraire métriques d'ajustement pour beta_try et nu_try (eq_days=j pr atteindre eq, sim_days=simuler incidence annuelle)
compute_metrics <- function(beta_try, nu_try, eq_days = 2000, sim_days = 365) {
  params_try <- c(beta = beta_try, 
                  nu = nu_try, 
                  gamma = gamma, 
                  sigma = sigma,
                  sigma_A = sigma_A, 
                  tau = tau, 
                  epsilon = epsilon, 
                  phi = phi)
  
  # intégration pour laisser converger système vers son équilibre
  eq_out <- as.data.frame(lsoda(init.cond, seq(0, eq_days, by = 1), cdiff_model, params_try))
  state_eq <- eq_out[nrow(eq_out), c("S0","SA","C0","CA","I","R")] #extraction dernière ligne tableau
  state_eq_num <- as.numeric(state_eq) #converti en vecteur nums
  names(state_eq_num) <- names(state_eq) #reassigne nom composantes``
  print(eq_out)
  
  # prévalence colonisée calculée avec éq
  preval_C <- (state_eq_num["C0"] + state_eq_num["CA"]) / sum(state_eq_num)
  
  # part de la FOI provenant de I calculée avec éq
  contrib_I <- nu_try * state_eq_num["I"]
  contrib_C <- state_eq_num["C0"] + state_eq_num["CA"]
  frac_from_I <- as.numeric(contrib_I / (contrib_I + contrib_C))
  
  # incidence annuelle calculée avec éq
  sim <- as.data.frame(lsoda(state_eq_num, seq(0, sim_days, by = 1), cdiff_model, params_try)) #simule à p de l'état d'équilibre pdt 1 an
  annual_incidence <- sum(params_try["sigma"] * sim$C0 + params_try["sigma_A"] * sim$CA) #somme brut des nouveaux cas : besoin conversion en patient/lits-jours
  return(list(prevalence_C = preval_C, frac_from_I = frac_from_I, annual_inc = annual_incidence)) 
}
#fin compute metric

#données ode
ode <- 

# appel test (verif metriques retournées)
m <- compute_metrics(beta_try=beta,nu_try=nu)

# affichage synthétique
print(m)

# intégration : fonction lsoda (pour résoudre edo) avec as.data.frame() autour pour un output plus pratique
result_ode <- as.data.frame(lsoda(init.cond, times, cdiff_model, params))

# méthode ggplot
result_ode %>%
  melt(id = "time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Value", colour = "Population:")
