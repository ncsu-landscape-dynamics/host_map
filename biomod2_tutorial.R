# load the library
require(biomod2)
# Load dependencies
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

# 1. load our species data
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv", package="biomod2"))
head(DataSpecies)

myRespName <- 'GuloGulo' # the name of studied species
myResp <- as.numeric(DataSpecies[,myRespName]) # the presence/absences data for our species
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")] # the XY coordinates of species data
# load the environmental raster layers (could be any supported format by the raster package)
# Environmental variables extracted from Worldclim 
myExpl = stack(system.file( "external/bioclim/current/bio3.grd", package="biomod2"),
               system.file( "external/bioclim/current/bio4.grd", package="biomod2"),
               system.file( "external/bioclim/current/bio7.grd", package="biomod2"),
               system.file( "external/bioclim/current/bio11.grd", package="biomod2"),
               system.file( "external/bioclim/current/bio12.grd", package="biomod2"))

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
myBiomodData
plot(myBiomodData)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Computing the models
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c('SRE',
                                               'GLM',
                                               'GAM',
                                               'MARS',
                                               'FDA',
                                               'CTA',
                                               'GBM',
                                               'RF',
                                               'ANN',
                                               'MAXENT.Phillips'),
                                    models.options = myBiomodOption,
                                    NbRunEval=1,
                                    DataSplit=80,
                                    Prevalence=0.5,
                                    VarImport=3,
                                    models.eval.meth = c('TSS'),
                                    SaveObj = TRUE,
                                    rescal.all.models = FALSE, #TRUE
                                    do.full.models = FALSE,
                                    modeling.id = paste(myRespName,"FirstModeling",sep=""))


myBiomodModelEval <- get_evaluations(myBiomodModelOut) # get all models evaluation
dimnames(myBiomodModelEval) # print the dimnames of this object
myBiomodModelEval[c('TSS'),"Testing.data",,,] # print the eval scores of all selected models
get_variables_importance(myBiomodModelOut) # print variable importances

# 3.2 Ensembling the models
myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      chosen.models = 'all',
                                      em.by='all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.7),
                                      prob.mean = F,
                                      prob.cv = F, #don't use
                                      prob.ci = F,
                                      #prob.ci.alpha = 0.05,
                                      prob.median = F,
                                      committee.averaging = F,
                                      prob.mean.weight = T,
                                      prob.mean.weight.decay = 'proportional')

#myBiomodEM # print summary
get_evaluations(myBiomodEM) # get evaluation scores

### 4. projection over the globe under current conditions
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = myExpl,
                                  proj.name = 'current',
                                  selected.models = 'all',
                                  binary.meth = 'TSS',
                                  compress = 'xz',
                                  clamping.mask = F,
                                  output.format = '.grd')

myBiomodProj # summary of crated oject
list.files("GuloGulo/proj_current/") # files created on hard drive
plot(myBiomodProj) # make some plots sub-selected by str.grep argument
myCurrentProj <- get_predictions(myBiomodProj) # if you want to make custom plots, you can also get the projected map

myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                         projection.output = myBiomodProj)

#myBiomodEF # print summary
plot(myBiomodEF) # reduce layer names for plotting convegences


##### This code creates a prediction for future climate scenarios
# ### load environmental variables for the future.
# myExplFuture = stack( system.file( "external/bioclim/future/bio3.grd",
#                                    package="biomod2"),
#                       system.file( "external/bioclim/future/bio4.grd",
#                                    package="biomod2"),
#                       system.file( "external/bioclim/future/bio7.grd",
#                                    package="biomod2"),
#                       system.file( "external/bioclim/future/bio11.grd",
#                                    package="biomod2"),
#                       system.file( "external/bioclim/future/bio12.grd",
#                                    package="biomod2"))
# myBiomodProjFuture <- BIOMOD_Projection(
#   modeling.output = myBiomodModelOut,
#   new.env = myExplFuture,
#   proj.name = 'future',
#   selected.models = 'all',
#   binary.meth = 'TSS',
#   compress = 'xz',
#   clamping.mask = T,
#   output.format = '.grd')