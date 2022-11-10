get_Envi <- function(bio=T, lc=T, rnr=T, soil=T, res=1){#borders=NULL, pop=T, crop=T,
  require(terra)
  require(geodata)
  
  if(bio==T){
    biodir <- 'C:\\Users\\bjselige\\Documents\\localdata\\wc2.0_30s_bio\\'
    #biodir <- 'Q:\\Shared drives\\APHIS  Projects\\shared resources\\data\\worldclim1k\\US\\'
    biovar <- rast(list.files(biodir, full.names=T))
    names(biovar) <-  c('Mean Annual Temp', 'Mean Diurnal Range', 'Isothermality', 'Temp Seasonality',
                        'Max Temp Warmest Month', 'Min Temp Coldest Month', 'Temp Annual Range',
                        'Mean Temp Wettest Quarter', 'Mean Temp Driest Quarter', 'Mean Temp Warmest Quarter',
                        'Mean Temp Coldest Quarter', 'Annual Precip', 'Precip Wettest Month',
                        'Precip Driest Month', 'Precip Seasonality', 'Precip Wettest Quarter',
                        'Precip Driest Quarter', 'Precip Warmest Quarter', 'Precip Coldest Quarter')
    biocl <- data.frame(var=names(biovar), cluster=NA)
    biocl$cluster[which(biocl$var%in%c('Mean Annual Temp',
                                       'Max Temp Warmest Month',
                                       'Mean Temp Wettest Quarter',
                                       'Mean Temp Warmest Quarter'))] <- 'Temp 1'
    biocl$cluster[which(biocl$var%in%c('Mean Diurnal Range',
                                       'Isothermality',
                                       'Temp Seasonality',
                                       'Min Temp Coldest Month',
                                       'Temp Annual Range',
                                       'Mean Temp Driest Quarter',
                                       'Mean Temp Coldest Quarter'))] <- 'Temp 2'
    biocl$cluster[which(biocl$var%in%c('Annual Precip',
                                       'Precip Wettest Month',
                                       'Precip Wettest Quarter',
                                       'Precip Warmest Quarter'))] <- 'Precip 1'
    biocl$cluster[which(biocl$var%in%c('Precip Seasonality',
                                       'Precip Driest Month',
                                       'Precip Driest Quarter',
                                       'Precip Coldest Quarter'))] <- 'Precip 2'
  }
  
  if(rnr==T){# roaddir <- 'C:\\Users\\bjselige\\Documents\\localdata\\Products_generated_from_Rails_Roads\\'
    # # roaddir <- 'Q:\\Shared drives\\APHIS  Projects\\shared resources\\data\\Rails_Roads\\Products_generated_from_Rails_Roads\\'
    rnrdir <- 'C:\\Users\\bjselige\\Desktop\\ne_rnr\\distance\\'
    rnrvar <- rast(list.files(rnrdir, '.dist.wrld', full.names = T))
    names(rnrvar) <- c('Road Dist', 'Rail Dist')
    rnrcl <- data.frame(var=names(rnrvar), cluster='Roadsrails 1')
  }
  
  # if(pop==T){
  #   popvar <- rast('C:\\Users\\bjselige\\Downloads\\gpw-v4-population-density-rev11_2020_30_sec_tif\\gpw_v4_population_density_rev11_2020_30_sec.tif')
  #   # ghm <- rast('C:\\Users\\bjselige\\Downloads\\lulc-human-modification-terrestrial-systems-geographic-geotiff\\lulc-human-modification-terrestrial-systems_geographic.tif')
  #   # cpl <- rast('C:\\Users\\bjselige\\Downloads\\ucl_2014_v23.tif')
  #   names(popvar) <- 'Population Density'
  #   popcl <- data.frame(var=names(popvar), cluster='Population 1')
  # }
  # 
  # if(crop==T){
  #   cpldir <- 'C:\\Users\\bjselige\\Desktop\\cropland\\'
  #   cpl <- cropland("WorldCover", path=cpldir)
  # }
  
  if(lc==T){
    lcdir <- 'C:\\Users\\bjselige\\Desktop\\landcover\\'
    cropl <- landcover(var="cropland", path=lcdir)
    built <- landcover(var='built', path=lcdir)
    # trees <- landcover(var='trees', path=lcdir)
    # shrub <- landcover(var='shrubs', path=lcdir)
    # grass <- landcover(var='grassland', path=lcdir)
    # wetld <- landcover(var='wetland', path=lcdir)
    # snowy <- landcover(var='snow', path=lcdir)
    # baren <- landcover(var='bare', path=lcdir)
    lcvar <- c(cropl, built); names(lcvar) <- c('Cropland', 'Built')
    lcvar <- extend(lcvar, ext(biovar))
    lccl <- data.frame(var=names(lcvar), cluster='Landcover 1')
  }
  
  if(soil==T){
    soildir <- 'C:\\Users\\bjselige\\Desktop\\Soil_OpenLandMap\\resample\\'
    soil.files <- c('Soil_pH_0cm.tif', 'Soil_pH_mean.1km.tif', 'Soil_pH_200cm.tif',
                    'Soil_h2o_33kpa_0cm.tif', 'Soil_h2o_33kpa_mean.1km.tif',
                    'Soil_h2o_33kpa_200cm.tif', 'Soil_h2o_1500kpa_0cm.tif',
                    'Soil_h2o_1500kpa_mean.1km.tif', 'Soil_h2o_1500kpa_200cm.tif')
    solvar <- rast(paste(soildir, soil.files, sep=''))
    names(solvar) <- c('Soil_pH_0cm', 'Soil_pH_mean', 'Soil_pH_200cm',
                       'Soil_h2o_33_0cm', 'Soil_h2o_33_mean', 'Soil_h2o_33_200cm',
                       'Soil_h2o_1500_0cm', 'Soil_h2o_1500_mean', 'Soil_h2o_1500_200cm')
    solcl <- data.frame(var=names(solvar), cluster='Soil 1')
  }
  
  if(bio==F){biovar <- NULL; biocl <- NULL}
  if(lc==F){lcvar <- NULL; lccl <- NULL}
  if(rnr==F){rnrvar <- NULL; rnrcl <- NULL}
  #if(pop==F){popvar <- NULL; popcl <- NULL}
  if(soil==F){solvar <- NULL; solcl <- NULL}
  
  envi <- c(biovar, lcvar, rnrvar, solvar)
  clst <- rbind(biocl, lccl, rnrcl, solcl)
  clst$cluster <- as.integer(as.factor(clst$cluster))
  cl2 <- clst$cluster; names(cl2) <- gsub(' ', '.', clst$var)
  if(res>1){envi <- terra::aggregate(envi, fact=res)}
  # if(!is.null(borders)){envi <- terra::crop(x=envi, y=borders, mask=T)}
  return(list('rast'=envi, 'clust'=cl2))
}