###Paquets
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

###Model (method2)
SIR.model <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
  
    N=S+I+R
    
    dS <- -beta*S*I/N + u0*N - u1*S
    dI <- beta*S*I/N  - gamma*I - u1*I
    dR <-  gamma*I - u1*R
    
    res<-c(dS, dI, dR)
    
    list(res)
    
  })

}

# Paramètres
beta=0.6
gamma=1/3
dt=1
Tmax=14
N=30
u0=0
u1=0

# Conditions initiales
I0=1
S0=N-I0
R0=0

# Création des vecteurs à partir des différentes valeurs
Time=seq(from=0,to=Tmax,by=dt)
Init.cond=c(S=S0,I=I0, R=R0) 
param=c(beta,gamma,u0,u1)

# Intégration : fonction lsoda (pour résoudre edo) avec as.data.frame() autour pour un output plus pratique
result <- as.data.frame(lsoda(Init.cond, Time, SIR.model, param))

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
