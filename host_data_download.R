# Load libraries ####
library(folderfun)
library(rgdal)
library(raster)
library(rgbif)
library(data.table)
library(foreach)
library(rFIA)
library(dplyr)
library(neonUtilities)
library(sp)
library(sf)
library(rinat)
# Set working directories####
setff("wd", "C:/Users/blaginh/Desktop/host_dl_data/")#personal wd
setff("sd", "Q:/Shared drives/Data/Vector/USA")#shared wd
# Setup reference area (study extent)####
l48<-readOGR(ffsd("us_lower_48_states.gpkg"))
ref_area_geom<-l48[l48$STATE_NAME=="Connecticut"|l48$STATE_NAME=="Rhode Island",]
ref_area_xy<-extent(ref_area_geom)
#Setup host species of interest####
sn <- "Ailanthus altissima"
# Import occurrence records from major species occurrence databases####
## GBIF ####
### Download using gbif package####
gbif.out <- occ_data(scientificName = sn, hasCoordinate=TRUE,geometry=
                       c(ref_area_xy[1],ref_area_xy[3],ref_area_xy[2],ref_area_xy[4]),
                     limit=100000, hasGeospatialIssue = FALSE)
gbif.out<-as.data.table(gbif.out[["data"]])
gbif.out<-gbif.out[ ,`:=`(networkKeys = NULL)]

dir.create(ffwd("gbif"))
write.csv(gbif.out,ffwd("gbif/gbif.out.csv"))#export raw datasheet
rm(gbif.out)
### Clean####
gbif<-as.data.table(read.csv(ffwd("gbif/gbif.out.csv")))#import raw datasheet
gbif<-gbif[grep("inaturalist",gbif$references, invert=T, ignore.case = T)]#remove inat (dl separately)
gbif<-gbif[basisOfRecord=="HUMAN_OBSERVATION",]#retain only human observations

# gbif.flags<-unique(gbif$issues)#review gbif issues flagged & remove inaccurate points
# gbif.flags<-unlist(stringr::str_split(gbif.flags,","))
# gbif.flags<-foreach(i=1:length(gbif.flags), .combine="rbind") %do% {
#   x<-grepl(gbif.flags[i],gbif.issues$code)
#   gbif.issues[x,]
# }

gbif_t<-gbif[,.(year=year,p_a=occurrenceStatus,sciname=sn,
                lat=decimalLatitude, lon=decimalLongitude, db="gbif")]
gbif_t$p_a[gbif_t$p_a=="PRESENT"]<-"present"
gbif_t$p_a[gbif_t$p_a=="ABSENT"]<-"absent"
#EddMaps -invasive hosts only####
## Download from https://www.eddmaps.org/distribution/, no api####
unzip(zipfile = ffwd("eddmaps/35328.zip"), file="mappings.csv", 
      exdir = ffwd("eddmaps/"))#path to zip file downloaded from eddmaps
eddmaps<-as.data.table(read.csv(ffwd("eddmaps/mappings.csv")))#import raw datasheet

## Clean####
eddmaps_t<-eddmaps[!is.na(Latitude),.(year=ObsDate, p_a=OccStatus,sciname=sn,
                                      lat=Latitude,lon=Longitude, db="eddmaps")]
eddmaps_t$p_a<-ifelse(eddmaps_t$p_a=="Detected",yes="present",no="absent")
eddmaps_t$year<-strptime(eddmaps_t$year, '%m-%d-%y')$year+1900
# FIA####
## Download using the FIA package####
l48_text<-c('CT','RI')#study extent
options(timeout=900)#increase timeout to accomodate file size
foreach(i=l48_text) %do% getFIA(states = i, nCores=(detectCores()-1),dir = ffwd("fia/"), 
                                tables=c("TREE","PLOT"))

rbindlist_fread <- function(path, pattern = "*.csv") 
  {files = list.files(path, pattern, full.names = TRUE)
  rbindlist(lapply(files, function(x) fread(x)))}#reads & binds multiple csvs w/ same column headers

l48_fia_tree<-rbindlist_fread(path=ffwd("fia/"),pattern="TREE.csv")
l48_fia_plot<-rbindlist_fread(path=ffwd("fia/"),pattern="PLOT.csv")
l48_fia_t<-unique(merge(l48_fia_tree[,c("PLT_CN","SPCD")],l48_fia_plot[,c("CN","MEASYEAR","LAT","LON","PLOT","UNITCD","COUNTYCD","STATECD","INVYR")], by.x="PLT_CN",by.y="CN", all.x=TRUE))
write.csv(l48_fia_t,ffwd("/fia/fia.raw.csv"))
rm(l48_fia_tree)
rm(l48_fia_t)

## Clean####
fia<-as.data.table(read.csv(ffwd("/fia/fia.raw.csv")))[,-1]
fia$PLT_CN<-as.factor(fia$PLT_CN)
fia.aa<-fia[SPCD==341,]#fia species code for species of interest
fia.aa<-unique(fia.aa[,.(p_a="present",year=max(INVYR),db="FIA",sciname=sn), by=.(PLOT,UNITCD,COUNTYCD,STATECD,LAT,LON)])
fia.aa[,uid:=as.character(paste0(PLOT,UNITCD,COUNTYCD,STATECD))]
names(fia.aa)[5:6]<-c("lat","lon")
fia.naa<-unique(fia[SPCD!=341,.(p_a="absent",year=max(INVYR),db="FIA",sciname=sn), by=.(PLOT,UNITCD,COUNTYCD,STATECD,LAT,LON)])
fia.naa[,uid:=as.character(paste0(PLOT,UNITCD,COUNTYCD,STATECD))]
names(fia.naa)[5:6]<-c("lat","lon")
aasites<-unique(fia.aa[,.(uid,lat,lon)])
fia.naa.t<-as.data.table(anti_join(fia.naa, aasites, by=c("uid","lat","lon")))
fia.t<-unique(rbind(fia.aa[,-c("uid","PLOT","UNITCD","COUNTYCD","STATECD")],fia.naa.t[,-c("uid","PLOT","UNITCD","COUNTYCD","STATECD")]))

# NEON####
## Download using neon package
neon.out<-loadByProduct(dpID="DP1.10058.001", site="all", nCores = (detectCores()-1))
dir.create(ffwd("neon"))
write.csv(neon.out$div_10m2Data100m2Data,ffwd("neon/neon.div10m2and100m2.csv"))#trees & shrubs
write.csv(neon.out$div_1m2Data,ffwd("neon/neon.div_1m2.csv"))#shorter plants (e.g., ferns, forbs, grasses)

## Clean ####
neon<-as.data.table(read.csv(ffwd("neon/neon.div_1m2.csv")))# if host is trees & shrubs
#neon<-as.data.table(read.csv(ffwd("neon/neon.div_1m2.csv")))#if host is not tree or shrub

neon.aa<-neon[grep(sn, neon$scientificName, ignore.case = T)]
unique(neon.aa$scientificName)#previous will pull all species, subsp., & var.; confirm records pulled 
neon.aa.t<-neon.aa[,.(p_a="present",date=endDate,sciname=sn, lat=decimalLatitude, lon=decimalLongitude, db="neon", plotID, subplotID)]
neon.aa.t<-unique(neon.aa.t[,.(p_a,year=max(date),lat,lon,db,sciname,plotID), by=.(plotID, subplotID)])

neon.naa<-neon[grep(sn, neon$scientificName, ignore.case = T, invert = T)]
neon.naa.t<-unique(neon.naa[,.(p_a="absent",date=endDate,sciname=sn, lat=decimalLatitude, lon=decimalLongitude, db="neon", plotID, subplotID)])
aasites<-unique(neon.aa[,.(plotID,subplotID)])
neon.naa.t<-as.data.table(anti_join(neon.naa.t, aasites, by=c("plotID","subplotID")))
neon.naa.t<-unique(neon.naa.t[,.(p_a,year=max(date),lat,lon,db,sciname,plotID), by=.(plotID, subplotID)])
neon.t<-unique(rbind(neon.aa.t[,-c("plotID","subplotID")],neon.naa.t[,-c("plotID","subplotID")]))
neon.t$year<-strptime(neon.t$year, '%Y-%m-%d')$year+1900

#BLM ####
## Download using online link, no api####
dir.create(ffwd("blm"))
download.file(url="https://gis.blm.gov/AIMDownload/LayerPackages/BLM_AIM_TerrADat.zip",destfile = ffwd("/blm/BLM_AIM_TerrADat.zip"))
unzip(zipfile = ffwd("/blm/BLM_AIM_TerrADat.zip"), exdir = ffwd("blm/"))

## Clean

blm<-readOGR(ffwd("blm/BLM_AIM_TerrADat/v107/TerrADat.gdb"))
blmWGS84<-spTransform(blm,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
blm<-as.data.table(blmWGS84)
blm.aa<-blm[grep("AIAL",blm$Spp_Nox)]#USDA Plants Code for Ailanthus altissima 
blm.aa_t<-blm.aa[,.(date=DateVisited,p_a="present",lat=coords.x2,lon=coords.x1,db="blm", sciname=sn, PlotID=PlotID)]
blm.aa_t<-unique(blm.aa_t[,.(p_a,year=max(date),lat,lon,db,sciname,PlotID), by=.(PlotID)])

blm.naa<-blm[!grep("AIAL",blm$Spp_Nox)]
aasites<-as.data.table(blm.aa_t[,.(PlotID,lat,lon)])
blm.naa<-as.data.table(anti_join(blm.naa, aasites, by="PlotID"))
blm.naa_t<-unique(blm.naa[,.(p_a="absent",year=max(DateVisited),lat=coords.x2,lon=coords.x1,db="blm",sciname=sn,PlotID), by=.(PlotID)])
blm_t<-unique(rbind(blm.aa_t[,-("PlotID")],blm.naa_t[,-("PlotID")]))

#iNaturalist####
inat<-get_inat_obs(query = sn, bounds = c(ref_area_xy[3],ref_area_xy[1],ref_area_xy[4],ref_area_xy[2]),maxresults = 10000)
dir.create(ffwd("inat"))
write.csv(inat,ffwd("/inat/inat.raw.csv"))

#Clean####
inat<-as.data.table(read.csv(ffwd("inat/inat.raw.csv")))
inat_t<-inat[coordinates_obscured=="false" & quality_grade=="research" ,.(year=datetime, p_a="present",sciname=sn,
                                             lat=latitude,lon=longitude, db="inaturalist")]
inat_t$year<-as.Date.character(inat_t$year)
inat_t$year<-format(inat_t$year,format="%Y")
inat_t
#Combine data, convert to pts, crop to ref area, export cropped pts as geopackage & csv####
aa.occ<-rbind(gbif_t,eddmaps_t,fia.t,neon.t,blm_t,inat_t)
aa.pts <- st_as_sf(x = aa.occ,
                        coords = c("lon", "lat"),
                        crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

aa.pts.c<-st_crop(aa.pts,ref_area_xy)
aa.pts.c.df<-as.data.table(aa.pts.c)
aa.pts.c.df<-cbind(aa.pts.c.df[,1:4],st_coordinates(aa.pts.c))
names(aa.pts.c.df)[5:6]<-c("lon","lat")
dir.create(ffwd("final_outputs"))
write.csv(aa.pts.c.df,ffwd("final_outputs/aa.pts.csv"))
write_sf(aa.pts.c,ffwd("final_outputs/aa.pts.gpkg"),driver="GPKG")

#To do ####
##1. add VMI NPS 
##2. add Carolina Vegetation Survey)
##3. double-check that BISON is captured by GBIF
##4. for inat - split ref area into grid when max results are reached, re-run, and combine results
##5. eddmap added manually - need to use api