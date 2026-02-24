library(ecospat)
library(raster)
library(dismo)
library(viridis)
library(ggplot2)
library(gridExtra)
library(ade4)
library(ggExtra)


setwd('C:/Niche_overlap/Data/regional_bioclimate/')

### Read csv archive of distributional data.
sp1 <- read.csv("Cenchrus_canarias.csv", sep = ';', dec =',')
sp2 <- read.csv("Melinis_canarias.csv", sep = ";", dec = ',')
sp3 <- read.csv("Setaria_canarias.csv", sep = ";", dec = ',')

colnames(sp1) <- c("y", "x")
sp1 <- sp1[,c(2,1)]


### Load previously cropped layers that are relevant for comparison of ecological niches.
varclim1 <- stack(list.files("C:/Niche_overlap/Data/regional_bioclimate/", pattern=".tif", full.names=T))

### Create a table for each pixel (x, y coordinates) and extract environmental data.
clim1 <- as.data.frame(varclim1, xy = TRUE)
clim1 <- na.omit(clim1)

### Extract climatic variables for species occurrences.
ex1 <- raster::extract(varclim1, sp1[,c(1,2)])
occ_sp1 <- na.omit(data.frame(sp1, ex1))

ex2 <- raster::extract(varclim1, sp2[,c(1,2)])
occ_sp2 <- na.omit(data.frame(sp2, ex2))

ex3 <- raster::extract(varclim1, sp3[,c(1,2)])
occ_sp3 <- na.omit(data.frame(sp3, ex3))


occ_sp1$group <- "sp1"
occ_sp2$group <- "sp2"
occ_sp3$group <- "sp3"
clim1$group <- "background"

# Combina los datos de forma segura
# Selecciona solo las columnas de variables y la nueva columna de grupo
data <- rbind(
  clim1[, c(3:5, 6)],
  occ_sp1[, c(3:5, 6)],
  occ_sp2[, c(3:5, 6)],
  occ_sp3[, c(3:5, 6)]
)

# Elimina los NA de la tabla final si es necesario (ya deberían estar limpios)
data <- na.omit(data)

#######################################################################
######################## Environmental PCA ########################
#######################################################################
### Data for PCA. All sites for the study area (clim1) and species presences.
# Definir los pesos para el PCA
# 1 para el fondo, 0 para las ocurrencias
row.w.env <- ifelse(data$group == "background", 1, 0)

# Realizar el PCA
pca.cal <- dudi.pca(data[, 1:3], row.w = row.w.env, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

# Ahora, asigna las puntuaciones del PCA de forma segura usando la columna 'group'
scores.clim <- pca.cal$li[data$group == "background", ]
scores.sp1 <- pca.cal$li[data$group == "sp1", ]
scores.sp2 <- pca.cal$li[data$group == "sp2", ]
scores.sp3 <- pca.cal$li[data$group == "sp3", ]

# Ya puedes continuar con la creación de las cuadrículas de nicho y los análisis
z1_grid <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp1, R = 100)
z2_grid <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp2, R = 100)
z3_grid <- ecospat.grid.clim.dyn(scores.clim, scores.clim, scores.sp3, R = 100)
z <- list(z1_grid, z2_grid, z3_grid)

### Metrics for niche overlap.
ecospat.niche.overlap(z1 = z1_grid, z2 = z2_grid, cor = TRUE)
ecospat.niche.overlap(z1 = z1_grid, z2 = z3_grid, cor = TRUE)
ecospat.niche.overlap(z1 = z2_grid, z2 = z3_grid, cor = TRUE)

### The rest of the plotting and equivalency test code should work fine with the corrected variables.
# You can uncomment and run the rest of the script.

# Combina los datos de PCA en un solo dataframe
data_plot <- as.data.frame(rbind(scores.clim, scores.sp1, scores.sp2, scores.sp3))
colnames(data_plot) <- c("PC1", "PC2")

# Añade una columna para identificar el grupo (fondo o especie)
data_plot$group <- c(
  rep("Fondo", nrow(scores.clim)),
  rep("Cenchrus", nrow(scores.sp1)),
  rep("Melinis", nrow(scores.sp2)),
  rep("Setaria", nrow(scores.sp3))
)

###ploting###

ggplot(data = data_plot, aes(x = PC1, y = PC2)) +
  
  # 3. Pinta los puntos de fondo con un color neutro
  # geom_point(aes(color = "Fondo"),
  #            data = subset(data_plot, group == "Fondo"),
  #            size = 0.5,
  #            alpha = 0.1) +
  
  # 1. Pinta los puntos de todas las especies
  geom_point(aes(color = group),
             data = subset(data_plot, group != "Fondo"),
             size = 1.5,
             alpha = 0.5) +
  
  # 2. Pinta el contorno de densidad de cada especie
  geom_density_2d(aes(color = group),
                  data = subset(data_plot, group != "Fondo"),
                  linewidth = 1) +
  
  # Configuración de colores para puntos y contornos
  scale_color_manual(name = "Species",
                     values = c("Cenchrus" = "#d799f5", "Melinis" = "red", "Setaria" = "#10b67c")) +
  
  # Etiquetas y título
  labs(
       x = "PC1",
       y = "PC2") +
  
  # Tema del gráfico
  theme_bw() +
  theme(legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8),
        legend.position = "bottom")


### Are two niches identical? 
### Null hypothesis is that niches are identical and you reject it if the observed value is lower than expected under a null model.
a.dyn<-ecospat.niche.equivalency.test(z1=z1_grid , z2=z2_grid, rep=100, ncores = 8)
ecospat.plot.overlap.test(a.dyn,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn,"I","Equivalency")
a.dyn2<-ecospat.niche.equivalency.test(z1=z1_grid, z2=z3_grid, rep=100, ncores = 8)
ecospat.plot.overlap.test(a.dyn2,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn2,"I","Equivalency")
a.dyn3<-ecospat.niche.equivalency.test(z1=z2_grid, z2=z3_grid, rep=100, ncores = 8)
ecospat.plot.overlap.test(a.dyn3,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn3,"I","Equivalency")


a.dyn4<-ecospat.niche.equivalency.test(z1=z2_grid , z2=z1_grid, rep=100, ncores = 8)
ecospat.plot.overlap.test(a.dyn,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn,"I","Equivalency")
a.dyn5<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z1_grid, rep=100, ncores = 8)
ecospat.plot.overlap.test(a.dyn2,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn2,"I","Equivalency")
a.dyn6<-ecospat.niche.equivalency.test(z1=z3_grid, z2=z2_grid, rep=100, ncores = 8)
ecospat.plot.overlap.test(a.dyn3,"D","Equivalency")
ecospat.plot.overlap.test(a.dyn3,"I","Equivalency")



###alternative plot###

ecospat.plot.niche (z1_grid, title="Cenchrus", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z2_grid, title="Melinis", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z3_grid, title="Setaria", name.axis1="PC1", name.axis2="PC2", cor=F)


r1 <- raster(z[[1]]$z.cor)

r2 <- raster(z[[2]]$z.cor)

r3 <- raster(z[[3]]$z.cor)

####añadiendo los rasters


raster_to_df <- function(raster, name){
  # raster: lista con $x, $y, $z.uncor (como antes)
  mat <- raster$z.uncor
  if(is.list(mat)) mat <- mat[[1]]  # extraer matriz si viene en lista
  if(length(dim(mat)) == 3 && dim(mat)[3] == 1) mat <- mat[,,1]
  
  df <- expand.grid(PC1 = raster$x, PC2 = raster$y)
  df$Density <- as.vector(mat)
  df$Species <- name
  return(df)
}

df1 <- raster_to_df(z[[1]], "Cenchrus")
df2 <- raster_to_df(z[[2]], "Melinis")
df3 <- raster_to_df(z[[3]], "Setaria")

raster_df <- rbind(df1, df2, df3)
raster_df$Species <- as.factor(raster_df$Species)

####Plotting##

plot_1 <- ggplot(data = subset(raster_df, Density > 0), aes(x = PC1, y = PC2, fill = Species)) +
  # Raster con densidades
  geom_tile(aes(alpha = Density)) +
  
  # Contornos de densidad para cada especie
  geom_contour(aes(z = Density, color = Species), bins = 5, linewidth = 0.7) +
  labs(title = "Baseline (1979-2013)")+
  
  # Puntos invisibles para ggMarginal, mapeando color a Species
  geom_point(aes(colour = Species), alpha = 0) +
  
  # Colores
  scale_fill_manual(values = c("Cenchrus" = "#d799f5",
                               "Melinis" = "red",
                               "Setaria" = "#10b67c"), labels = c("C. setaceus", "M. repens", "S. parviflora")) +
  scale_color_manual(values = c("Cenchrus" = "#d799f5",
                                "Melinis" = "red",
                                "Setaria" = "#10b67c"), labels = c("C. setaceus", "M. repens", "S. parviflora")) +
  scale_alpha_continuous(range = c(0,0.5), guide = "none") +
  
  # Etiquetas y tema
  labs(x = "PC1", y = "PC2") +
  theme_bw() +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9, face = "italic"),
        legend.position = "bottom")
  

# Añadir curvas marginales separadas por especie
ggMarginal(plot_1, type = "density", groupColour = TRUE, groupFill = TRUE)







