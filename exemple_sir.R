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

# Paramètres
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
N=30
I_0=1
R_0=0
C0_0=4
CA_0=1
S0_0=N-C0_0-CA_0-I_0
SA_0=0


# Création des vecteurs à partir des différentes valeurs
times=seq(from=0, to=Tmax, by=dt)
init.cond=c(S0=S0_0, SA=SA_0, C0=C0_0, CA=CA_0, I=I_0, R=R_0) 
params=c(beta, nu, gamma, sigma, sigma_A, tau, epsilon, phi)

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
