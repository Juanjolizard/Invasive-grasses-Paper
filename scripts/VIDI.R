####variable importance index###
###ver si multiplicar por 100 o no##
###especies##

library(dplyr)
library(tidyr)

setwd('C:/Niche_overlap/Modelling_results/Cenchrus-results/Results/Global/Values/')
csv_global <- read.csv('Cenchrus_setaseus_indvar.csv')
csv_global$perc <- csv_global$var.imp*100

###agrupamiento global##

contri <- csv_global %>%
  group_by(expl.var) %>%
  summarise(mean_contribution = mean(perc, na.rm = TRUE))

contri$Modelo <- 'Global'

###local##
setwd('C:/Niche_overlap/Modelling_results/Cenchrus-results/Results/Regional/Values/')
csv_regional <- read.csv('Cenchrus_setaseus_indvar.csv')
csv_regional$perc <- csv_regional$var.imp*100

contri_1 <- csv_regional %>%
  group_by(expl.var) %>%
  summarise(mean_contribution = mean(perc, na.rm = TRUE))

contri_1$Modelo <- 'Regional'

###pivot wider##
pivotado_global <- contri %>%
  pivot_wider(names_from = expl.var, values_from = mean_contribution)

pivotado_regional <- contri_1 %>%
  pivot_wider(names_from = expl.var, values_from = mean_contribution)

####diferencias de cuadrados##
df_final <- rbind(pivotado_global, pivotado_regional)

####suma de diferencias cuadráticas##
df_final[2,2]
VIDI <- (df_final[1,2]-df_final[2,2])^2 + (df_final[1,3]-df_final[2,3])^2 + (df_final[1,4]-df_final[2,4])^2
