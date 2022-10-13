# function(spname, #name of species. this is the only required parameter, all else are optional
#          data, merge, #data will allow user to supply own points. merge defines whether or not they should be merged with latest BIEN points
#          threshold, f.score, #threshold creates a thresholded version of the output. f.score allows user to tune this threshold to their use
#          envi, var.select, #envi allows user to supply own environmental data. var.select allows for turning off variable selection algorithm
#          bbox, output) #bbox allows user to set own extent, default is sized by extent of the data. output allows for saving outputs locally
EM_mod <- function(spname){
  require(biomod2)
  require(plyr)
  require(dismo)
  require(raster)
  require(rgdal)
  require(scales)
  require(stringr)
  require(klaR)
  
  
  #update w world borders
  borders <- vect('C:\\Users\\bjselige\\Downloads\\ne_10m_admin_0_countries_lakes\\ne_10m_admin_0_countries_lakes.shp')
  
  
  
  # bbox <- c(24.5, -125, 49.5, -66.5)
  # us_can <- borders[borders$admin%in%c('United States of America', 'Canada'),]#, 'Mexico'),]
  # us_can <- us_can[which(!us_can$name%in%c('Alaska', 'Hawaii', 'Nunavut', 'Yukon', 'Northwest Territories')),]
  # l48 <- us_can[which(us_can$admin=='United States of America'),]
  # us_can <- crop(us_can, extent(l48))
  
  # bbox <- c(49.5, -8, 59.5, 2)
  # uk_main <- borders[borders$admin%in%c('United Kingdom'),]
  # uk_main <- uk_main[which(!uk_main$region%in%c('Northern Ireland')),]
  # uk_main <- crop(uk_main, extent(c(bbox[2], bbox[4], bbox[1], bbox[3])))
  
  # 
  source('C:\\Users\\bjselige\\host_map\\get_Envi.R')
  envi <- get_Envi(rnr=F)

  
  
  # Create clusters of related variables among variables in the environmental layer
  envi.cc <- corclust(x=data.frame(values(envi)))
  envi.cp <- plot(envi.cc, mincor=.7, selection='numeric')
  envi.cv <- cvtree(envi.cc, k=6) #update k so its dynamic
  
  source('C:\\Users\\bjselige\\host_map\\get_pts.1.R')
  pts <- get_pts.1(spname)
  pts <- data.frame('longitude'=as.numeric(pts$longitude), 'latitude'=as.numeric(pts$latitude))
  pts <- pts[which(!is.na(pts$latitude)),]
  pts <- SpatialPoints(pts, CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
  pts <- crop(pts, extent(envi))
  pts.ras <- rasterize(x=pts, y=envi, fun='count', background=0)
  pts.ras <- (pts.ras*(envi[[1]]*0+1))>0
  pts.2 <- rasterToPoints(pts.ras)
  
  myRespName <- str_replace(tolower(spname),' ', '_')
  myResp <- pts.2[, 3] # the presence/absences data for our species
  myResp[myResp==0] <- NA # setting 'true absences' to NA
  myRespXY <- pts.2[, c(1,2)] # the XY coordinates of species data
  
  #var impo goes here
  source('C:\\Users\\bjselige\\host_map\\var.imp.R')
  myExpl.sel <- var.select()
  
  myOptions <- BIOMOD_ModelingOptions()
  
  myData <- BIOMOD_FormatingData(resp.var = myResp,
                                 expl.var = myExpl.sel,
                                 resp.xy = myRespXY,
                                 resp.name = myRespName,
                                 PA.nb.rep = 1,
                                 PA.strategy = 'random',
                                 PA.nb.absences = sum(myResp, na.rm=T))
  
  # Notes on algorithm choice; CTA is redundant with Random Forest, FDA and SRE have relatively low performance
  myAlgos <- c('GLM', 'GAM', 'MARS', 'GBM',  'RF', 'ANN')#'CTA', 'FDA', 'SRE',  'MAXENT.Phillips','MAXENT.Phillips.2'
  # Notes on evaluation methods : POD/SR/FR is not useful, KAPPA, ACCURACY, TSS, and ETS all get about the same results.
  myEvals <- c('CSI', 'ROC', 'TSS', 'KAPPA', 'ACCURACY', 'ETS', 'BIAS') #'POD', 'FAR', 'SR'
  # 3. Computing the models
  myModels <- biomod2::BIOMOD_Modeling(data=myData, 
                                       models.options = myOptions,
                                       models = myAlgos,
                                       NbRunEval = 1, #3, #number of runs
                                       DataSplit = 80,
                                       Prevalence = 0.5,
                                       VarImport = 1,
                                       models.eval.meth = myEvals,
                                       SaveObj = T, # recommended to leave true
                                       rescal.all.models = F, #experimental don't use
                                       do.full.models = F, # use this option if you don't want to use a datasplit
                                       modeling.id=paste(myRespName,"FirstModeling",sep=""))
  
  # Evaluation
  myEval <- get_evaluations(myModels) # get all models evaluation
  # dimnames(myEval) # print the dimnames of this object
  # tss.df <- data.frame(myEval[c('TSS'),"Testing.data",,,])
  # tss.df <- rbind(cbind(tss.df, rowMeans(tss.df)), colMeans(tss.df))
  # colnames(tss.df)[length(tss.df)] <- c('Algo.avg')
  # rownames(tss.df)[nrow(tss.df)] <- c('Run.avg')
  # tss.df$Algo.avg[nrow(tss.df)] <- NA
  # tss.df$Algo.avg[nrow(tss.df)] <- mean(unlist(tss.df['Run.avg',], tss.df$Algo.avg), na.rm=T)
  # tss.df <- signif(tss.df,3)
  
  # 3.2 Ensembling the models
  myEnsemble <- BIOMOD_EnsembleModeling(modeling.output = myModels,
                                        chosen.models = 'all',
                                        em.by = 'PA_dataset+repet',
                                        eval.metric = NULL,
                                        #eval.metric.quality.threshold = .5,
                                        models.eval.meth= myEvals,
                                        prob.mean = T,
                                        prob.cv = F, #don't use'
                                        prob.ci = F, #prob.ci.alpha = 0.05,
                                        prob.median = F,
                                        committee.averaging = F,
                                        prob.mean.weight = F, # leave on true
                                        prob.mean.weight.decay = 'proportional')
  myEvalEM <- get_evaluations(myEnsemble)[[1]] # get evaluation scores
  
  ### 4. projection over the globe under current conditions
  myProj <- BIOMOD_Projection(modeling.output = myModels,
                              new.env = myExpl.sel,
                              proj.name = 'current',
                              selected.models = 'all',
                              binary.meth = myEvals,
                              compress = 'xz',
                              clamping.mask = F,
                              output.format = '.grd',
                              do.stack=T)
  myProj2 <- get_predictions(myProj) # if you want to make custom plots, you can also get the projected map
  
  myBinary <- llply(.data=c(1:length(myEvals)),
                    .fun=function(X){X.Eval <- myEvals[[X]]
                    X.Binary <- stack(paste(tolower(str_replace(spname, ' ', '.')), '\\proj_current\\proj_current_',
                                            tolower(str_replace(spname, ' ', '.')), '_', X.Eval, 'bin.grd', sep=''))
                    names(X.Binary) <- paste(X.Eval, substr(myProj@models.projected, nchar(spname)+1, max(nchar(myProj@models.projected))), sep='')
                    return(X.Binary)}); myBinary <- stack(unlist(myBinary))
  
  # myBinary <-stack(paste(tolower(str_replace(spname, ' ', '.')), '\\proj_current\\proj_current_',
  #                        tolower(str_replace(spname, ' ', '.')), '_', myEvals, 'bin.grd', sep=''))
  # names(myBinary) <- rep(paste(myEvals, substr(myModels@models.computed, 15, 27), sep='_'), length(myEvals))
  # myBinary <- stack(list.files('C:\\Users\\bjselige\\host_map\\juglans.nigra\\proj_current\\individual_projections\\', 'bin.grd', full=T))
  # names(myBinary) <- rep(paste(myEvals, substr(myModels@models.computed, 15, 27), sep='_'), length(myEvals))
  
  myProjEM <- BIOMOD_EnsembleForecasting(EM.output = myEnsemble,
                                         projection.output = myProj,
                                         selected.models = 'all',
                                         binary.meth = myEvals,
                                         compress = 'xz',
                                         clamping.mask = F,
                                         output.format = '.grd')
  
  myBinaryEM <- stack(paste(tolower(str_replace(spname, ' ', '.')), '\\proj_current\\proj_current_',
                            tolower(str_replace(spname, ' ', '.')), '_ensemble_', myEvals[[1]], 'bin.grd', sep=''))
  names(myBinaryEM) <- myEvals
  
  p.out <- myProjEM@proj@val[[1]]
  
  ##### 5. Thresholding
  p2 <- p.out
  trs <- seq(min(values(p2), na.rm=T), max(values(p2), na.rm=T),
             by=(max(values(p2), na.rm=T) - min(values(p2), na.rm=T))/100)
  trs.metrics <- ldply(.data=c(1:length(trs)), .fun=function(X){
    p2.tr <- p2>trs[X]
    x.zero <- sum(raster::extract(x=p2.tr, y=pts)==0, na.rm=T)/length(pts)
    x.area <- sum(values(p2.tr), na.rm=T)/length(na.omit(p2[]))
    x.score <- 1 - x.zero - x.area
    df.out <- data.frame('trs'=trs[X], 'zero'=x.zero, 'area'=x.area, 'score'=x.score, row.names=trs[X])
    return(df.out)}, .progress = 'text')
  
  tr.best <- trs.metrics[which.max(trs.metrics$score),]
  p.tr <- p.out>tr.best$trs
  p3 <- stack(myProjEM@proj@val[[1]], p.tr); names(p3) <- c('Raw', 'Binary')
  p3 <- stack(p3, myBinaryEM)
  
  myEval2 <- list(Eval=myEval, Ensemble=myEvalEM, BJS.metrics=list('best'=tr.best, 'full'=trs.metrics))
  
  outlist <- list(myData, myModels, myEval2, myProj2, p3)
  names(outlist) <- c('Data', 'Model', 'Evaluation', 'Projections', 'Ensemble')
  return(outlist)
}

require(plyr)

splist <- c(#'Ailanthus altissima', #treeofheaven 
  'Buxus', #Boxwood
  'Juglans nigra',# Black walmnut
  # 'Lonicera hispidula', # Honey suckle
  # 'Notholithocarpus densiflorus', # Tanoak
  # 'Pseudotsuga menziesii',
  # # 'Tamarix chinensis',
  # 'Tsuga canadensis',
  # #  'Tsuga caroliniana',
  # 'Umbellularia californica' # Bay laurel
)

sp_ems <- llply(.data=c(1:length(splist)),
                .fun=function(x){EM_mod(splist[[x]])},
                .progress='text'); names(sp_ems) <- splist


# pts <- BIEN::BIEN_occurrence_species(species=spname)
# pts <-  read.csv('C:\\Users\\bjselige\\Documents\\tree_of_heaven\\Ailanthus.BIEN.csv')  
# pts <- gbif(genus=strsplit(x=spname, split=' ')[[1]][1], species=strsplit(x=spname, split=' ')[[1]][2],
#             geo=T, sp=T, removeZeros=T, download=T, ext=extent(biovars[[1]]))
#source('C:\\Users\\bjselige\\host_map\\get_pts.1.R')
#pts <- get_pts.1(spname)

# #### output plot
# borders <- usa
# rast <- myProjEM@proj@val
# rpts <- rasterToPoints(rast)
# rdf <- as.data.frame(rpts)
# ggsdm <- ggplot() + geom_raster(data=rdf, aes(x=x, y=y, fill=rdf[,3])) + 
#   geom_path(data=borders, aes(x=long, y=lat, group=group), col='white', lwd=1.1, alpha=.3) + 
#   scale_fill_continuous(type='viridis') + 
#   theme_void() + theme(legend.position='none')
# png(paste('C:\\Users\\bjselige\\Documents\\Tree_of_Heaven\\Figures\\usa.', 
#           gsub(':', '', substr(Sys.time(), 12, 19)), '.png', sep=''), 
#     height=1080, width=2160); plot(ggsdm); dev.off()