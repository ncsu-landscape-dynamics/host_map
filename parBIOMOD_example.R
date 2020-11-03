require(biomod2) 
require(abind)
require(ade4)
require(caret)
require(checkmate)
require(dismo)
require(doParallel)
require(dplyr)
require(earth)
require(ecospat)
require(ENMeval)
require(foreach)
require(foreign)
require(gam)
require(gbm)
require(ggplot2)
require(Hmisc)
require(lattice)
require(MASS)
require(maxnet)
require(mda)
require(mgcv)
require(methods)
require(nnet)
require(parallel)
require(PresenceAbsence)
require(pROC)
require(purrr)
require(randomForest)
require(raster)
require(rasterVis)
require(reshape)
require(rlang)
require(rpart)
require(sp)
require(stats)
require(testthat)
require(tidyr)
require(utils)
require(rgdal)

#usa <- readOGR('C:\\Users\\bjselige\\Downloads\\us_lower_48_states.shp')
usa <- readOGR('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\usa_boundaries\\us_lower_48_states.shp')
usa <- spTransform(usa, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

# 1. load our species data
#aa.bien <- read.csv('C:\\Users\\bjselige\\Documents\\tree_of_heaven\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]  
#aa.gbif <- read.csv('C:\\Users\\bjselige\\Documents\\tree_of_heaven\\Ailanthus.GBIF.csv')[, c(2,3)]; names(aa.gbif) <- c('Longitude', 'Latitude')
aa.gbif <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.GBIF.csv')[, c(2,3)]; names(aa.gbif) <- c('Longitude', 'Latitude')
aa.bien <- read.csv('H:\\Shared drives\\Data\\Table\\Global\\Ailanthus.BIEN.csv')[, c('Longitude', 'Latitude')]
aa.pts <- SpatialPoints(coords = unique(rbind(aa.gbif, aa.bien)))
aa.pts <- crop(aa.pts, usa)

# load the environmental raster layers (could be any supported format by the raster package)
# Environmental variables extracted from Worldclim 
#myExpl <- raster::getData('worldclim', download=T, var='bio', res=10)
biodir <- 'H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\worldclim1k\\US\\'
biovars <- stack(lapply(X=list.files(biodir), FUN=function(X){raster(paste(biodir, X, sep=''))}))
myExpl <- crop(stack(biovars[[1]], biovars[[6]], biovars[[12]]), extent(usa))
myExpl <- stack(aggregate(myExpl, 50))
myExpl <- crop(myExpl, extent(usa))
myExpl <- stack(raster::mask(myExpl, usa))
aa.ras <- rasterize(x=aa.pts, y=myExpl[[1]], fun='count', background=0); aa.ras <- (aa.ras*(myExpl[[1]]*0+1))>0
a2.pts <- rasterToPoints(aa.ras)
myRespName <- 'A_altissima'
myResp <- a2.pts[, 3] # the presence/absences data for our species
myRespXY <- a2.pts[, c(1,2)] # the XY coordinates of species data

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
myBiomodData
plot(myBiomodData)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Computing the models
t1 <- system.time(BIOMOD_Modeling(myBiomodData, 
                                  models = c('SRE',
                                             'GLM',
                                             'CTA',
                                             'FDA',
                                             'RF',
                                             'MARS',
                                             'GAM',
                                             'ANN',
                                             'GBM',
                                             'MAXENT.Phillips'),
                                  models.options = myBiomodOption,
                                  NbRunEval = 1, #5
                                  DataSplit = 80,
                                  Prevalence = 0.5,
                                  VarImport = 3,
                                  models.eval.meth = 'TSS',
                                  SaveObj = TRUE,
                                  rescal.all.models = FALSE,
                                  do.full.models = FALSE,
                                  modeling.id = paste(myRespName, "FirstModeling", sep="")))

source('C:\\Users\\bjselige\\host_map\\parBIOMOD_Modeling.R')
environment(parBIOMOD_Modeling) <- asNamespace('biomod2')
assignInNamespace("BIOMOD_Modeling", parBIOMOD_Modeling, ns = "biomod2")

# 4. Computing the parallel model
t2 <- system.time(biomod2::BIOMOD_Modeling(myBiomodData,
                                           models = c('SRE',
                                                      'GLM',
                                                      'CTA',
                                                      'FDA',
                                                      'RF',
                                                      'MARS',
                                                      'GAM',
                                                      'ANN',
                                                      'GBM',
                                                      'MAXENT.Phillips'),
                                           models.options = myBiomodOption,
                                           NbRunEval = 1, #5
                                           DataSplit = 80,
                                           Prevalence = 0.5,
                                           VarImport = 3,
                                           models.eval.meth = 'TSS',
                                           SaveObj = TRUE,
                                           rescal.all.models = FALSE,
                                           do.full.models = FALSE,
                                           modeling.id=paste(myRespName,"FirstModeling",sep="")))
