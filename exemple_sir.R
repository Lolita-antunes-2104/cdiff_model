
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
# si les packages ne chargent pas, c'est qu'ils ne sont pas installés! Donc il faut faire install.packages("nomdupackage")

# principe: on créé une fonction pour le modèle, puis on utilise la fonction lsoda() du package deSolve et on lui donne les conditions initiales,
# les timesteps désirés, les paramètres, et la fonction du modèle. la fonction lsoda() retourne le résultat du modèle à partir de ces infos.
 

# méthode 1 pour créer la fonction du modèle
# on extrait chaque paramètre et population à partir des vecteurs pop et param qu'on donne à la fonction
# attention aux indices, sinon on se retrouve à extraire la mauvaise valeur pour le paramètre/population! 
SIR.model <- function(t, pop, param) {
  
  S <- pop[1]
  I <- pop[2]
  R <- pop[3]
  
  beta=param[1]
  gamma=param[2]
  u0=param[3]
  u1=param[4]
  
  N=S+I+R
  
  dS <- -beta*S*I/N + u0*N - u1*S
  dI <- beta*S*I/N  - gamma*I - u1*I
  dR <-  gamma*I - u1*R
  
  res<-c(dS, dI, dR)
  
  list(res)
  
}

# méthode 2 pour créer la fonction du modèle
# with(as.list()) ça charge les vecteurs de population et paramètres comme ça on a pas besoin
# d'extraire chaque paramètre comme au dessus (S = pop[1] etc...)
# attention à bien nommer tous les éléments correctement quand on créé les vecteurs en dessous!
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

# on défini les paramètres
beta=0.6
gamma=1/3
dt=1
Tmax=14
N=30
u0=0
u1=0

# on défini les conditions initiales
I0=1
S0=N-I0
R0=0

# on créé les vecteurs à partir des différentes valeurs
Time=seq(from=0,to=Tmax,by=dt)
Init.cond=c(S=S0,I=I0, R=R0) 
param=c(beta,gamma,u0,u1)

# on utilise la fonction lsoda, avec as.data.frame() autour pour un output plus pratique
result <- as.data.frame(lsoda(Init.cond, Time, SIR.model, param))

# méthode "basique" pour plot
plot(Time,result$S,type="l",col="blue",xlab="Time (days)",ylab="Prevalence",ylim=c(0,N),bty="n", lwd=2)
lines(Time,result$I,type="l",col="orange", lwd=2)
lines(Time,result$R,type="l",col="green", lwd=2)
legend("topright",c("S","I","R"),col=c("blue","orange","green"),lty=1,bty="n", lwd=2)

# méthode ggplot
ggplot(result) +
  geom_line(aes(time, S, colour = "S")) +
  geom_line(aes(time, I, colour = "I")) +
  geom_line(aes(time, R, colour = "R")) +
  theme_bw()

# méthode ggplot + reshape (pour réorganiser le dataframe, plus pratique)
result %>%
  melt(id = "time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw()

# méthode ggplot + reshape + joli
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
