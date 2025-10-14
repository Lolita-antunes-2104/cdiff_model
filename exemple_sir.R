###Paquets
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

###Model (method2)
cdiff_model <- function(t, pop, params) {
  
  with(as.list(c(pop, params)), {
  
    #population totale
    N = S0 + SA + C0 + CA + I + R
    
    #force d'infection
    lambda = beta * (C0 + CA + nu*I)/N
    
    #equations diff
    dS0 <- -lambda*S0 + gamma*C0 - tau*S0 
    dSA <- -lambda*SA + gamma*CA + tau*S0 + phi*R
    dC0 <- lambda*S0 - gamma*C0 - sigma*C0 - tau*C0
    dCA <- lambda*SA - gamma*CA - sigma_A*CA + tau*C0
    dI <- sigma*C0  + sigma_A*CA - epsilon*I + sigma_A*R
    dR <-  epsilon*I - sigma_A*R - phi*R
    
    res<-c(dS0, dSA, dC0, dCA, dI, dR)
    
    list(res)
    
  })

}

# Parametres
beta=0.021
nu=36
gamma=0.03
sigma=0.000445
sigma_A=sigma*6.67
tau=0.044
epsilon=0.07
phi=0.018

#Simulation
dt=1
Tmax=140

# Conditions initiales
N=1000
I_0=10
R_0=0
C0_0=40
CA_0=10
S0_0=N-C0_0-CA_0-I_0
SA_0=0


# Création des vecteurs à partir des différentes valeurs
times=seq(from=0, to=Tmax, by=dt)
init.cond=c(S0=S0_0, SA=SA_0, C0=C0_0, CA=CA_0, I=I_0, R=R_0) 
params=c(beta=beta, nu=nu, gamma=gamma, sigma=sigma, sigma_A=sigma_A, tau=tau, epsilon=epsilon, phi=phi)

# Fonction pour extraire métriques d'ajustement pour beta_try et nu_try (eq_days=j pr atteindre eq, sim_days=simuler incidence annuelle)
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
  names(state_eq_num) <- names(state_eq) #reassigne nom composantes
  
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

#Appel test (verif metriques retournées)
m <- compute_metrics(beta_try=beta,nu_try=nu)

#Affichage synthétique
print(m)

# Intégration : fonction lsoda (pour résoudre edo) avec as.data.frame() autour pour un output plus pratique
result <- as.data.frame(lsoda(init.cond, times, cdiff_model, params))

# Méthode ggplot
result %>%
  melt(id = "time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Value", colour = "Population:")

# fonction pour enregistrer le plot (enregistre par défaut le dernier plot affiché)
# ggsave("nom_du_plot.png")
