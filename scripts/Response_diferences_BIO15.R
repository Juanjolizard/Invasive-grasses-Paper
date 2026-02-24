####analysing variable responses for modelling###
###BIO 15##

library(raster)
library(ggplot2)
library(scales)
library(dplyr)
library(metrica)
library(mgcv)

###modelo global y predictor global###
modelo_global <- raster('C:/Niche_overlap/outputs/Setaria_parviflora//current/Global.tif')
modelo_global <- as.data.frame(modelo_global, xy = T)
temp_global <- raster('C:/Niche_overlap/Data/global_bioclimate/CHELSA_bio10_15.tif')

extraccion <- raster::extract(temp_global, modelo_global[,c(1:2)])

df_modelo_global <- data.frame(modelo_global, extraccion)
df_modelo_global <- na.omit(df_modelo_global)
colnames(df_modelo_global) <- c('x','y','Suit','BIO15')
df_modelo_global$modelo <- 'Global'

# ggplot(df_modelo_global, aes(x = BIO15, y = Suit))+
#   geom_point(alpha = 0.5)+
#   geom_smooth()


####modelo regional y predictor regional###
modelo_regional <- raster('C:/Niche_overlap/outputs/Setaria_parviflora/current/Regional.tif')
modelo_regional <- as.data.frame(modelo_regional, xy = T)
temp_regional <- raster('C:/Niche_overlap/Data/regional_bioclimate/CHELSA_bio10_15.tif')

extraccion <- raster::extract(temp_regional, modelo_regional[,c(1:2)])

df_modelo_regional <- data.frame(modelo_regional, extraccion)
df_modelo_regional <- na.omit(df_modelo_regional)
colnames(df_modelo_regional) <- c('x','y','Suit','BIO15')
df_modelo_regional$modelo <- 'Regional'



###ajustar rango de la variable global al rango nativo para comparar### Será ajustar rango del local al global no? ATENCION
rango <- range(df_modelo_regional$BIO15)

df_modelo_global_filtrado <- dplyr::filter(df_modelo_global, BIO15 >= rango[1] & BIO15 <= rango[2])

###escalar respuestas y unir
df_modelo_global_filtrado$Suit <- scale(df_modelo_global_filtrado$Suit)
df_modelo_regional$Suit <- scale(df_modelo_regional$Suit)

df_final <- rbind(df_modelo_global_filtrado, df_modelo_regional)



###plot final###
ggplot(df_final, aes(x = BIO15, y = Suit, colour = modelo))+
  geom_smooth()+
  scale_colour_manual(name = "Model", values = c("gold", "#27F58B"))+
  xlab ("BIO 15")+
  ylab ("Scaled Habitat Suitability")+
  theme_bw()+
  theme(legend.position = "bottom")

range(df_modelo_global$BIO15)
###grdi##
gam_global <- gam(Suit ~ s(BIO15, k = 10), data = df_modelo_global_filtrado)
# gam_global1 <- gam(Suit ~ s(BIO1, k = 6), data = df_modelo_global_filtrado)
# gam_global2 <- gam(Suit ~ s(BIO1, bs = "cr"), data = df_modelo_global_filtrado)
# gam_global3 <- gam(Suit ~ te(BIO1), data = df_modelo_global_filtrado)
# gam_global4 <- gam(Suit ~ s(BIO1, k = 3), data = df_modelo_global_filtrado)
# 
# AIC(gam_global1, gam_global2, gam_global3, gam_global4) ####gam global 2##
# gam_global2


gam_regional <- gam(Suit ~ s(BIO15, k = 10), data = df_modelo_regional)


# ###analizar gams###
# concurvidad_global <- mgcv::concurvity(gam_global1)
# concurvidad_regional <- mgcv::concurvity(gam_regional)
# 
# concurvidad_regional
# Secuencia común de BIO1
bio_seq <- data.frame(
  BIO15 = seq(rango[1], rango[2], length.out = 1000)
)

# Predicciones
pred_g <- predict(gam_global, newdata = bio_seq, type = "response")
pred_r <- predict(gam_regional, newdata = bio_seq, type = "response")

# Escalar predicciones
pred_g_scaled <- scale(pred_g)[,1]
pred_r_scaled <- scale(pred_r)[,1]

df_predicciones <- data.frame(pred_g_scaled, pred_r_scaled)

cc <- metrica::CCC(data = df_predicciones, pred_g_scaled, pred_r_scaled)
cc ###concordance correlation coefficient de las predicciones 0.6846868 -0.808: mal##

# Calcular suma de cuadrados de diferencias
grdi <- sum((pred_g_scaled - pred_r_scaled)^2)

cat("GRDI (Gradient Response Deviation Index) para BIO15:", round(grdi, 4), "\n") ##3612.94