get_Envi <- function(borders=NULL, bio=T, pop=T, rnr=T, soil=T, res=1){
  require(terra)
  
  if(bio==T){
    biodir <- 'C:\\Users\\bjselige\\Documents\\localdata\\wc2.0_30s_bio\\'
    #biodir <- 'Q:\\Shared drives\\APHIS  Projects\\shared resources\\data\\worldclim1k\\US\\'
    #biovars <- suppressWarnings(raster::getData('worldclim', download=T, var='bio', res=5))
    #biovars <- stack(lapply(X=list.files(biodir, full.names=T), FUN=function(X){suppressWarnings(raster(X))}))
    biovars <- rast(list.files(biodir, full.names=T))
    #biovars <- terra::resample(biovars, base, method='bilinear', threads=T)
    # biovars <- terra::crop(biovars, borders)
    # biovars <- terra::mask(biovars, borders)
    names(biovars) <-  c('Mean Annual Temp', 'Mean Diurnal Range', 'Isothermality', 'Temp Seasonality',
                         'Max Temp Warmest Month', 'Min Temp Coldest Month', 'Temp Annual Range',
                         'Mean Temp Wettest Quarter', 'Mean Temp Driest Quarter', 'Mean Temp Warmest Quarter',
                         'Mean Temp Coldest Quarter', 'Annual Precip', 'Precip Wettest Month',
                         'Precip Driest Month', 'Precip Seasonality', 'Precip Wettest Quarter',
                         'Precip Driest Quarter', 'Precip Warmest Quarter', 'Precip Coldest Quarter')
  }
  
  if(rnr==T){
    roaddir <- 'C:\\Users\\bjselige\\Documents\\localdata\\Products_generated_from_Rails_Roads\\'
  # roaddir <- 'Q:\\Shared drives\\APHIS  Projects\\shared resources\\data\\Rails_Roads\\Products_generated_from_Rails_Roads\\'
  roads.d <- rast(paste(roaddir, 'roads.dist.NE2.tif', sep=''))
  rails.d <- rast(paste(roaddir, 'rails.dist.NE2.tif', sep=''))
  roads.uk <- rast(paste(roaddir, 'roads.dist.UK.tif', sep=''))
  rails.uk <- rast(paste(roaddir, 'rails.dist.UK.tif', sep=''))
  roads.d <- terra::merge(roads.d, roads.uk)
  rails.d <- terra::merge(rails.d, rails.uk)
  roads.d <- terra::crop(roads.d, extent(borders)) 
  rails.d <- terra::crop(rails.d, extent(borders))
  roads.d <- terra::resample(roads.d, biovars[[1]], method='bilinear', threads=T)
  rails.d <- terra::resample(rails.d, biovars[[1]], method='bilinear', threads=T)
  rnr.2 <- stack(roads.d, rails.d); names(rnr.2) <- c('Roads Dist', 'Rails Dist')
  # rnr.d <- min(stack(roads.d, rails.d)); rnr.2 <- stack(roads.d, rails.d, rnr.d); names(rnr.2) <- c('Roads Dist.', 'Rails Dist.', 'Roads & Rails Dist.')
  }
  
  if(pop==T){
    pop.den <- rast('C:\\Users\\bjselige\\Downloads\\gpw-v4-population-density-rev11_2020_30_sec_tif\\gpw_v4_population_density_rev11_2020_30_sec.tif')
    # pop.den <- terra::crop(pop.den, extent(borders))
    # pop.den <- terra::mask(pop.den, borders)
    # pop.den <- terra::resample(pop.den, biovars[[1]], method='bilinear', threads=T)
    # pop.den <- stack(pop.den)
    names(pop.den) <- 'Population Density'
  }
  
  if(soil==T){soildir <- 'C:\\Users\\bjselige\\Desktop\\Soil_OpenLandMap\\resample\\'
  soil.files <- c('Soil_pH_0cm.tif', 'Soil_pH_mean.1km.tif', 'Soil_pH_200cm.tif',
                  'Soil_h2o_33kpa_0cm.tif', 'Soil_h2o_33kpa_mean.1km.tif', 'Soil_h2o_33kpa_200cm.tif',
                  'Soil_h2o_1500kpa_0cm.tif', 'Soil_h2o_1500kpa_mean.1km.tif', 'Soil_h2o_1500kpa_200cm.tif')
  sl <- rast(paste(soildir, soil.files, sep=''))
  names(sl) <- c('pH_0cm', 'pH_mean', 'pH_200cm', 'h2o_33_0cm', 'h2o_33_mean',
                 'h2o_33_200cm', 'h2o_15k_0cm', 'h2o_15k_mean', 'h2o_15k_200cm')
  }
  
  if(bio==F){biovars <- NULL}
  if(rnr==F){rnr.2 <- NULL}
  if(pop==F){pop.den <- NULL}
  if(soil==F){sl <- NULL}
  
  # if(is.null(borders)){
  #   borders <- vect('C:\\Users\\bjselige\\Downloads\\ne_10m_admin_0_countries_lakes\\ne_10m_admin_0_countries_lakes.shp')
  # }
  # if(file.exists('C:\\Users\\bjselige\\Documents\\localdata\\base.tif')==F){
  #   wrld <- vect('C:\\Users\\bjselige\\Downloads\\ne_10m_admin_0_countries_lakes\\ne_10m_admin_0_countries_lakes.shp')
  #   biodir <- 'C:\\Users\\bjselige\\Documents\\localdata\\wc2.0_30s_bio\\'
  #   biovars <- rast(list.files(biodir, full.names=T))
  #   base <- biovars[[1]]/biovars[[1]]
  #   base <- terra::mask(base, wrld, filename='C:\\Users\\bjselige\\Documents\\localdata\\base.tif')
  # }
  # if(file.exists('C:\\Users\\bjselige\\Documents\\localdata\\base.tif')){
  #   base <- rast('C:\\Users\\bjselige\\Documents\\localdata\\base.tif')
  #   }
  # if(ext(borders)!=ext(base)){}
  
  envi <- c(biovars, pop.den, rnr.2, sl)
  
  if(res>1){envi <- terra::aggregate(envi, fact=res)}
  if(!is.null(borders)){
    envi <- crop(envi, ext(borders))
    envi <- mask(envi, borders, threads=T)
  }
  
  return(envi)
}
# cl <- raster('C:\\Users\\bjselige\\Desktop\\nlcd\\cultivated_xra30.tif',
#              crs='+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
# cl.2 <- projectRaster(cl, crs='+proj=longlat +datum=WGS84 +no_defs ')
# cl.2 <- resample(cl.2, biovars[[1]]); names(cl.2) <- 'Cultivated'