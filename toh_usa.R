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
roads.d <- raster('H:\\Shared drives\\APHIS  Projects\\shared resources\\data\\Rails_Roads\\Products_generated_from_Rails_Roads\\roads.distance.tif')
roads.d <- resample(roads.d, biovars[[1]], method='bilinear')
myExpl <- crop(stack(biovars[[1]], biovars[[6]], biovars[[12]], roads.d), extent(usa))
#myExpl <- stack(aggregate(myExpl, 100))
myExpl <- crop(myExpl, extent(usa))
myExpl <- stack(raster::mask(myExpl, usa))
aa.ras <- rasterize(x=aa.pts, y=myExpl[[1]], fun='count', background=0); aa.ras <- (aa.ras*(myExpl[[1]]*0+1))>0
a2.pts <- rasterToPoints(aa.ras)
myRespName <- 'A_altissima'
myResp <- a2.pts[, 3] # the presence/absences data for our species
myResp[myResp==0] <- NA # setting 'true absences' to undefined
myRespXY <- a2.pts[, c(1,2)] # the XY coordinates of species data

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.strategy = 'random',
                                     PA.nb.absences = sum(myResp, na.rm=T))
myBiomodData
plot(myBiomodData)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Computing the models
myBiomodModelOut <- biomod2::BIOMOD_Modeling(myBiomodData,
                                             models = c(#'CTA', 'SRE',  'MAXENT.Phillips','MAXENT.Phillips.2'
                                               'GLM',
                                               'GAM',
                                               'MARS',
                                               'FDA',
                                               'GBM',
                                               'RF',
                                               'ANN'),
                                             models.options = myBiomodOption,
                                             NbRunEval = 5,
                                             DataSplit = 80,
                                             Prevalence = 0.5,
                                             VarImport = 3,
                                             models.eval.meth = 'TSS',
                                             SaveObj = TRUE,
                                             rescal.all.models = FALSE,
                                             do.full.models = FALSE,
                                             modeling.id=paste(myRespName,"FirstModeling",sep=""))


myBiomodModelEval <- get_evaluations(myBiomodModelOut) # get all models evaluation
dimnames(myBiomodModelEval) # print the dimnames of this object
myBiomodModelEval[c('TSS'),"Testing.data",,,] # print the eval scores of all selected models
vars_importance <- data.frame(get_variables_importance(myBiomodModelOut)) # print variable importances
# vars_ranked <- data.frame('SRE'=as.integer(rank(vars_importance[,1])),
#                           'GLM'=as.integer(rank(vars_importance[,2])),
#                           'GAM'=as.integer(rank(vars_importance[,3])),
#                           'MARS'=as.integer(rank(vars_importance[,4])),
#                           'FDA'=as.integer(rank(vars_importance[,5])),
#                           'CTA'=as.integer(rank(vars_importance[,6])),
#                           'GBM'=as.integer(rank(vars_importance[,7])),
#                           'RF'=as.integer(rank(vars_importance[,8])),
#                           'ANN'=as.integer(rank(vars_importance[,9])),
#                           'MAXENT.Phillips'=as.integer(rank(vars_importance[,10])))
# vars_ranked[,11] <- rowSums(vars_ranked)

# 3.2 Ensembling the models
myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      chosen.models = 'all',
                                      em.by='all',
                                      eval.metric = c('TSS'),
                                      eval.metric.quality.threshold = c(0.5),
                                      prob.mean = F,
                                      prob.cv = F, #don't use
                                      prob.ci = F, #prob.ci.alpha = 0.05,
                                      prob.median = F,
                                      committee.averaging = F,
                                      prob.mean.weight = T,
                                      prob.mean.weight.decay = 'proportional' )

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
plot(myBiomodProj)
myCurrentProj <- get_predictions(myBiomodProj) # if you want to make custom plots, you can also get the projected map

myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                         projection.output = myBiomodProj)

plot(myBiomodEF) # reduce layer names for plotting convegences

pred.out <- myBiomodEF@proj@val[[1]]
values(pred.out)  <- rescale(values(pred.out), to=c(0,1),
                             from=c(min(values(pred.out), na.rm=T), max(values(pred.out), na.rm=T)))
writeRaster(pred.out, filename = 'C:\\Users\\bjselige\\Desktop\\toh.conus_ensemble.tif', format="GTiff")

# #### output plot
# borders <- usa
# rast <- myBiomodEF@proj@val
# rpts <- rasterToPoints(rast)
# rdf <- as.data.frame(rpts)
# ggsdm <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
#   geom_path(data=borders, aes(x=long, y=lat, group=group), col='white', lwd=1.1, alpha=.3) + 
#   scale_fill_continuous(type='viridis') + 
#   theme_void() + theme(legend.position='none')
# 
# png(paste('C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\Figures\\usa.', 
#           gsub(':', '', substr(Sys.time(), 12, 19)), '.png', sep=''), 
#     height=1080, width=2160); plot(ggsdm); dev.off()