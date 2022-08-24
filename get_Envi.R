get_Envi <- function(bio=T, pop=T, rnr=T, soil=T){
  # 1. load data
  if(bio==F){bio <- NULL}
  if(bio==T){
  biovars <- raster::getData('worldclim', download=T, var='bio', res=10)
  biovars <- crop(biovars, extent(usa))
  #biodir <- 'C:\\Users\\bjselige\\Documents\\localdata\\US\\'
  #biodir <- 'Q:\\Shared drives\\APHIS  Projects\\shared resources\\data\\worldclim1k\\US\\'
  #biovars <- stack(lapply(X=list.files(biodir), FUN=function(X){suppressWarnings(raster(paste(biodir, X, sep='')))}))
  names(biovars) <-  c('Mean Annual Temp', 'Mean Diurnal Range', 'Isothermality',
                       'Temp Seasonality', 'Max Temp Warmest Month', 'Min Temp Coldest Month',
                       'Temp Annual Range', 'Mean Temp Wettest Quarter', 'Mean Temp Driest Quarter',
                       'Mean Temp Warmest Quarter', 'Mean Temp Coldest Quarter', 'Annual Precip',
                       'Precip Wettest Month', 'Precip Driest Month', 'Precip Seasonality',
                       'Precip Wettest Quarter', 'Precip Driest Quarter', 'Precip Warmest Quarter', 'Precip Coldest Quarter')
  }
  
  if(rnr==F){rnr <- NULL}
  if(rnr==T){#roaddir <- 'Q:\\Shared drives\\APHIS  Projects\\shared resources\\data\\Rails_Roads\\Products_generated_from_Rails_Roads\\'
  roaddir <- 'C:\\Users\\bjselige\\Documents\\localdata\\Products_generated_from_Rails_Roads\\'
  roads.d <- suppressWarnings(raster(paste(roaddir, 'roads.distance.tif', sep='')))
  rails.d <- suppressWarnings(raster(paste(roaddir, 'rails.distance.1km.tif', sep='')))
  roads.d <- resample(roads.d, biovars[[1]], method='bilinear')
  rails.d <- resample(rails.d, biovars[[1]], method='bilinear')
  rnr.d <- min(stack(roads.d, rails.d))
  rnr.2 <- stack(roads.d, rails.d, rnr.d); names(rnr.2) <- c('Roads Dist.', 'Rails Dist.', 'Roads & Rails Dist.')
  }
  
  if(pop==F){pop.den <- NULL}
  if(pop==T){popdir <- 'C:\\Users\\bjselige\\Downloads\\gpw-v4-population-density-rev11_2020_30_sec_asc\\'
  if(file.exists(paste(popdir, 'gpw_v4_population_density_rev11_2020_30_sec_12.asc', sep=''))==F){
    popdir <- 'C:\\Users\\bjselige\\Downloads\\gpw-v4-population-density-rev11_2020_30_sec_asc\\'
    pop.den <- merge(raster(paste(popdir, 'gpw_v4_population_density_rev11_2020_30_sec_1.asc', sep='')),
                     raster(paste(popdir, 'gpw_v4_population_density_rev11_2020_30_sec_2.asc', sep='')))
    names(pop.den) <- 'popl.density'
    pop.den <- crop(pop.den, extent(biovars[[1]]))
    pop.den <- resample(pop.den, biovars[[1]], method='bilinear')
    pop.den[pop.den<0] <- 0
    pop.den2 <- pop.den
    pop.den2[pop.den<1] <- 0.1
    writeRaster(pop.den, paste(popdir, 'gpw_v4_population_density_rev11_2020_30_sec_12.asc', sep=''), overwrite=T)
  }
  pop.den <- raster(paste(popdir, 'gpw_v4_population_density_rev11_2020_30_sec_12.asc', sep=''))
  pop.den <- resample(pop.den, biovars[[1]], method='bilinear'); names(pop.den) <- 'Population Density'
  }

  if(soil==F){soil <- NULL}
  if(soil==T){#soildir <- 'Q:\\Shared drives\\Data\\Raster\\USA\\soils\\'
  soildir <- 'C:\\Users\\bjselige\\Documents\\localdata\\soils\\'
  sl <- stack(lapply(X=list.files(soildir), FUN=function(X){suppressWarnings(raster(paste(soildir, X, sep='')))}))
  names(sl) <- c('Soil Density', 'Soil Clay', 'Soil Organic C', 'Soil pH', 'Soil Sand', 'Soil Texture', 'Soil Moisture')
  sl <- sl[[c(1:5,7)]]
  }
  
  envi <- stack(biovars, pop.den, rnr.2, soil)
  return(envi)
}
# cl <- raster('C:\\Users\\bjselige\\Desktop\\nlcd\\cultivated_xra30.tif',
#              crs='+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
# cl.2 <- projectRaster(cl, crs='+proj=longlat +datum=WGS84 +no_defs ')
# cl.2 <- resample(cl.2, biovars[[1]]); names(cl.2) <- 'Cultivated'