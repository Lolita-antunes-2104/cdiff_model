###############################################################################
############################### 0 : PACKAGES ###################################
###############################################################################

if (!require(deSolve)) install.packages("deSolve");     library(deSolve)
if (!require(ggplot2)) install.packages("ggplot2");     library(ggplot2)
if (!require(reshape2)) install.packages("reshape2");   library(reshape2)
if (!require(dplyr)) install.packages("dplyr");         library(dplyr)
if (!require(patchwork)) install.packages("patchwork"); library(patchwork)
if (!require(gtable)) install.packages("gtable");       library(gtable)
if (!require(gridExtra)) install.packages("gridExtra"); library(gridExtra)
if (!require(pbapply)) install.packages("pbapply");     library(pbapply)
if (!require(tidyr)) install.packages("tidyr");         library(tidyr)

# Packages de base (pas d'installation nécessaire)
library(grid)
library(parallel)


