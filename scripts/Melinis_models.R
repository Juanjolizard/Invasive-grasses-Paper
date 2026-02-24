library(terra)
library(covsel)
library(biomod2)
library(ecospat)
library(fs)
library(sgsR)

# Load the sabinaNSDM package
library(sabinaNSDM)

setwd('C:/Niche_overlap/Data/global_bioclimate/')

##prepare occurence data##
##global##
especie <- 'Melinis_repens'

sp_global <- read.csv('Mel_rep_global.csv')
sp_global1 <- read.csv('Mel_rep_global_gbif.csv')

sp_global <- sp_global[,c(2,1)]
colnames(sp_global) <- c('x','y')

sp_global_1 <- sp_global1[,c(2,1)]
colnames(sp_global_1) <- c('x','y')
sp_global_final <- rbind(sp_global, sp_global_1)
##local

sp_local <- read.csv('C:/Niche_overlap/Data/regional_bioclimate/Mel_rep_Canary.csv')
sp_local <- sp_local[,c(2,1)]
colnames(sp_local) <- c('x','y')

sp_local_1 <- read.csv('C:/Niche_overlap/Data/regional_bioclimate/Mel_rep_Canary_gbif.csv')
sp_local_1 <- sp_local_1[,c(2,1)]
colnames(sp_local_1) <- c('x','y')
sp_local_2 <- read.csv('C:/Niche_overlap/Data/regional_bioclimate/Melinis_propias.csv', sep = ';')
sp_local_2 <- sp_local_2[,c(2,3)]

sp_local_final <- rbind(sp_local, sp_local_1, sp_local_2)

###covariables: Mismos nombres, misma proyección, formato spatRaster##
lista_rasters <- list.files(pattern = '*.tif')
var_global <- terra::rast(lista_rasters)
names(var_global)

setwd('../regional_bioclimate/')
lista_rasters <- list.files(pattern = '*.tif')
var_regional <- terra::rast(lista_rasters)
names(var_regional)

##data
###cambiar el directorio de resultados##
setwd('C:/Niche_overlap/Modelling_results/Mellinis/')

nsdm_input <- NSDM.InputData(SpeciesName = especie,
                             spp.data.global = sp_global_final, 
                             spp.data.regional = sp_local_final, 
                             expl.var.global = var_global, 
                             expl.var.regional = var_regional,
                             # new.env = new.env,
                             # new.env.names = "scenario1",
                             Background.Global = NULL, 
                             Background.Regional = NULL,
                             Absences.Global = NULL,
                             Absences.Regional = NULL)

###formato de entrada##
nsdm_finput <- NSDM.FormattingData(nsdm_input,
                                   nPoints = 10000, # number of background points
                                   Min.Dist.Global = "resolution",
                                   Min.Dist.Regional = "resolution",
                                   Background.method = "stratified", # method “random" or "stratified” to generate background points 
                                   save.output = F) #save outputs locally
#####de momento que no cree carpeta con resultados, puntos globales griddeados a la resolución global de chelsa##


##VARIABLE SELECTION## 3 y 3, normalmente.
nsdm_selvars <- NSDM.SelectCovariates(nsdm_finput,
                                      maxncov.Global = 5,   # Max number of covariates to be selected at the global scale
                                      maxncov.Regional = 7, # Max number of covariates to be selected at the regional scale
                                      corcut = 0.7, #  correlation threshold
                                      algorithms = c("glm","gam","rf"),
                                      ClimaticVariablesBands = NULL,# covariate bands to be excluded in the covariate selection at the regional scale
                                      save.output = F)


summary(nsdm_selvars) ###see details
###MODELO GLOBAL##
nsdm_global <- NSDM.Global(nsdm_selvars,
                           algorithms = c("GAM","GBM", "RF", "MAXNET","GLM"),# Statistical algorithms used for modelling
                           CV.nb.rep = 10, # number of cross-validation repetitions
                           CV.perc = 0.8, # percentage of the data will be used for training in each cross-validation fold
                           metric.select.thresh = 0.8, #  AUC threshold to include replicates in the final ensemble model
                           CustomModelOptions = NULL, # Allows users to apply custom modelling options. 
                           save.output = TRUE, 
                           rm.biomod.folder = TRUE) # Remove the temporary folders created by `biomod2`
###MODELO REGIONAL###
nsdm_regional <- NSDM.Regional(nsdm_selvars,
                               algorithms = c("GAM","GBM", "RF", "MAXNET","GLM"),
                               CV.nb.rep = 10,
                               CV.perc = 0.8,
                               metric.select.thresh = 0.7,
                               CustomModelOptions = NULL, 
                               save.output = TRUE,
                               rm.biomod.folder = TRUE
)

###COGIENDO MODELO GLOBAL COMO COVARIABLE, probar AUC > 0.7 y 0.8.##
nsdm_covariate <- NSDM.Covariate(nsdm_global,
                                 algorithms = c("GAM","GBM", "RF", "MAXNET","GLM"),
                                 rm.corr=TRUE,
                                 CV.nb.rep = 10,
                                 CV.perc = 0.8,
                                 metric.select.thresh = 0.7,
                                 CustomModelOptions = NULL,
                                 save.output = TRUE,
                                 rm.biomod.folder = TRUE)

###MULTIPLY: ANIDADO DE MODELOS##

nsdm_multiply <- NSDM.Multiply(nsdm_global,
                               nsdm_regional,
                               method = "Arithmetic", # Method for averate model outputs: "Arithmetic" or "Geometric" mean
                               rescale = FALSE,
                               save.output=TRUE)