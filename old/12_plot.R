# Créer les tableaux sous forme de plots
library(gridExtra)
library(grid)
library(dplyr)

# ============================================================================
# TABLEAU 1 : COMPARTIMENTS (Populations N, Porteurs C, Infectés I)
# ============================================================================

tab1 <- tableau_compartiments

# Formater les valeurs en notation scientifique si nécessaire
tab1_format <- tab1 %>%
  mutate(across(where(is.numeric), ~sprintf("%.2e", .)))

# Créer le plot du tableau 1
plot1 <- tableGrob(tab1_format, rows = NULL, 
                   theme = ttheme_default(
                     base_size = 8,
                     core = list(fg_params = list(hjust = 1, x = 0.95)),
                     colhead = list(fg_params = list(fontface = "bold"))
                   ))

grid.newpage()
grid.draw(plot1)

png("tableau1_compartiments.png", width = 14, height = 4, units = "in", res = 300)
grid.draw(plot1)
dev.off()

# ============================================================================
# TABLEAU 2 : OUTPUTS (Prévalences en % et Incidences pour 100,000 PA)
# ============================================================================

tab2 <- tableau_outputs %>%
  mutate(
    # Prévalences en %
    Prev_h_nv = 100 * Prev_h_nv,
    Prev_c_nv = 100 * Prev_c_nv,
    Prev_h_v  = 100 * Prev_h_v,
    Prev_c_v  = 100 * Prev_c_v,
    # Incidences pour 100,000 personnes-années
    Inc_h_nv  = Inc_h_nv * 365 * 1e5,
    Inc_c_nv  = Inc_c_nv * 365 * 1e5,
    Inc_h_v   = Inc_h_v  * 365 * 1e5,
    Inc_c_v   = Inc_c_v  * 365 * 1e5
  )

# Formater les valeurs
tab2_format <- tab2 %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), "NA", sprintf("%.4f", .))))

# Créer le plot du tableau 2
plot2 <- tableGrob(tab2_format, rows = NULL,
                   theme = ttheme_default(
                     base_size = 9,
                     core = list(fg_params = list(hjust = 1, x = 0.95)),
                     colhead = list(fg_params = list(fontface = "bold"))
                   ))

grid.newpage()
grid.draw(plot2)

png("tableau2_outputs.png", width = 12, height = 4, units = "in", res = 300)
grid.draw(plot2)
dev.off()

# ============================================================================
# TABLEAU 3 : VERSION COMPACTE - Uniquement prévalences et incidences clés
# ============================================================================

tab3 <- tableau_outputs %>%
  mutate(
    # Prévalences en %
    Prev_h_nv = 100 * Prev_h_nv,
    Prev_c_nv = 100 * Prev_c_nv,
    Prev_h_v  = 100 * Prev_h_v,
    Prev_c_v  = 100 * Prev_c_v,
    # Incidences pour 100,000 personnes-années
    Inc_h_nv  = Inc_h_nv * 365 * 1e5,
    Inc_c_nv  = Inc_c_nv * 365 * 1e5,
    Inc_h_v   = Inc_h_v  * 365 * 1e5,
    Inc_c_v   = Inc_c_v  * 365 * 1e5
  ) %>%
  select(Scenario, Temps, 
         Prev_h_nv, Prev_c_nv, 
         Inc_h_nv, Inc_c_nv)

# Renommer les colonnes pour plus de clarté
colnames(tab3) <- c("Scenario", "Temps", 
                    "Prev_hop(%)", "Prev_com(%)", 
                    "Inc_hop(100kPA)", "Inc_com(100kPA)")

tab3_format <- tab3 %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), "NA", sprintf("%.4f", .))))

plot3 <- tableGrob(tab3_format, rows = NULL,
                   theme = ttheme_default(
                     base_size = 10,
                     core = list(fg_params = list(hjust = 1, x = 0.95)),
                     colhead = list(fg_params = list(fontface = "bold"))
                   ))

grid.newpage()
grid.draw(plot3)

png("tableau3_outputs_compact.png", width = 10, height = 4, units = "in", res = 300)
grid.draw(plot3)
dev.off()


# ============================================================================
# TABLEAU 4 : COMPARTIMENTS ET OUTPUTS - L'un au-dessus de l'autre
# ============================================================================

# Préparer le tableau compartiments
tab_comp_format <- tableau_compartiments %>%
  mutate(across(where(is.numeric), ~sprintf("%.2e", .)))

# Préparer le tableau outputs avec conversions
tab_out_format <- tableau_outputs %>%
  mutate(
    # Prévalences en %
    Prev_h_nv = 100 * Prev_h_nv,
    Prev_c_nv = 100 * Prev_c_nv,
    Prev_h_v  = 100 * Prev_h_v,
    Prev_c_v  = 100 * Prev_c_v,
    # Incidences pour 100,000 PA
    Inc_h_nv  = Inc_h_nv * 365 * 1e5,
    Inc_c_nv  = Inc_c_nv * 365 * 1e5,
    Inc_h_v   = Inc_h_v  * 365 * 1e5,
    Inc_c_v   = Inc_c_v  * 365 * 1e5
  ) %>%
  mutate(across(where(is.numeric), ~ifelse(is.nan(.), "NA", sprintf("%.4f", .))))

# Créer les grobs pour chaque tableau
grob_comp <- tableGrob(tab_comp_format, rows = NULL,
                       theme = ttheme_default(
                         base_size = 8,
                         core = list(fg_params = list(hjust = 1, x = 0.95)),
                         colhead = list(fg_params = list(fontface = "bold"))
                       ))

grob_out <- tableGrob(tab_out_format, rows = NULL,
                      theme = ttheme_default(
                        base_size = 8,
                        core = list(fg_params = list(hjust = 1, x = 0.95)),
                        colhead = list(fg_params = list(fontface = "bold"))
                      ))

# Ajouter des titres
titre_comp <- textGrob("COMPARTIMENTS (Populations N, Porteurs C, Infectés I)", 
                       gp = gpar(fontsize = 12, fontface = "bold"))
titre_out <- textGrob("OUTPUTS (Prévalences en %, Incidences pour 100,000 PA)", 
                      gp = gpar(fontsize = 12, fontface = "bold"))

# Combiner les tableaux verticalement avec leurs titres
plot4 <- arrangeGrob(
  titre_comp,
  grob_comp,
  titre_out,
  grob_out,
  ncol = 1,
  heights = c(0.5, 4, 0.5, 4)
)

grid.newpage()
grid.draw(plot4)

png("tableau4_complet_vertical.png", width = 14, height = 10, units = "in", res = 300)
grid.draw(plot4)
dev.off()

