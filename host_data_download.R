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
library(doParallel)
library(BIEN)

# Set working directories####
setff("wd", "C:/Users/blaginh/Desktop/host_dl_data/")#personal wd
setff("sd", "Q:/Shared drives/Data/Vector/USA")#shared wd

# Setup reference area (study extent)####
l48<-readOGR(ffsd("us_lower_48_states.gpkg"))
ref_area_geom<-l48[l48$STATE_NAME=="Connecticut"|l48$STATE_NAME=="Rhode Island",]
ref_area_xy<-extent(ref_area_geom)

#Setup host species of interest####
sn <- "Ailanthus altissima"
usdacodes<-as.data.table(read.csv("/Users/blaginh/Desktop/host_dl_data/usdaplantscodes.txt"))
sn_code<-usdacodes[grep("Ailanthus altissima",usdacodes$Scientific.Name.with.Author),]
sn_code<-sn_code$Symbol

# Import occurrence records from major species occurrence databases####
## GBIF #### - bucket
### Download using gbif package####
gbif.out <- occ_data(scientificName = sn, hasCoordinate=TRUE,geometry=
                       c(ref_area_xy[1],ref_area_xy[3],ref_area_xy[2],ref_area_xy[4]),
                     limit=100000, hasGeospatialIssue = FALSE)
gbif.out<-as.data.table(gbif.out[["data"]])
gbif.out<-gbif.out[ ,`:=`(networkKeys = NULL)]

dir.create(ffwd("gbif"))
write.csv(gbif.out,ffwd("gbif/gbif.out.csv"))#export raw datasheet

### Clean####
gbif<-as.data.table(read.csv(ffwd("gbif/gbif.out.csv")))#import raw datasheet

#unique(gbif$references)#review unique sources & remove inat or other db dups
gbif<-gbif[grep("inaturalist",gbif$references, invert=T, ignore.case = T)]#remove inat (dl separately)
gbif<-gbif[basisOfRecord=="HUMAN_OBSERVATION",]#retain only human observations

# gbif.flags<-unique(gbif$issues)#review gbif issues flagged & remove inaccurate points
# gbif.flags<-unlist(stringr::str_split(gbif.flags,","))
# gbif.flags<-foreach(i=1:length(gbif.flags), .combine="rbind") %do% {
#   x<-grepl(gbif.flags[i],gbif.flags$code)
#   gbif.flags[x,]
#   }

gbif_t<-gbif[,.(year=year,p_a=occurrenceStatus,sciname=sn,
                lat=decimalLatitude, lon=decimalLongitude, db="gbif")]
gbif_t$p_a[gbif_t$p_a=="PRESENT"]<-"present"
gbif_t$p_a[gbif_t$p_a=="ABSENT"]<-"absent"

# #EddMaps -invasive hosts only####
# ## Download from https://www.eddmaps.org/distribution/, no api####
unzip(zipfile = ffwd("eddmaps/35328.zip"), file="mappings.csv",
      exdir = ffwd("eddmaps/"))#path to zip file downloaded from eddmaps
eddmaps<-as.data.table(read.csv(ffwd("eddmaps/mappings.csv")))#import raw datasheet

# ## Clean####
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

# #Retain all points
# fia$PLT_CN<-as.factor(fia$PLT_CN)
# fia.aa<-fia[SPCD==341,]#fia species code for species of interest
# fia.aa<-unique(fia.aa[,.(p_a="present",year=max(INVYR),db="FIA",sciname=sn), by=.(PLOT,UNITCD,COUNTYCD,STATECD,LAT,LON)])
# fia.aa[,uid:=as.character(paste0(PLOT,UNITCD,COUNTYCD,STATECD))]
# names(fia.aa)[5:6]<-c("lat","lon")
# fia.naa<-unique(fia[SPCD!=341,.(p_a="absent",year=max(INVYR),db="FIA",sciname=sn), by=.(PLOT,UNITCD,COUNTYCD,STATECD,LAT,LON)])
# fia.naa[,uid:=as.character(paste0(PLOT,UNITCD,COUNTYCD,STATECD))]
# names(fia.naa)[5:6]<-c("lat","lon")
# aasites<-unique(fia.aa[,.(uid,lat,lon)])
# fia.naa.t<-as.data.table(anti_join(fia.naa, aasites, by=c("uid","lat","lon")))
# fia.t<-unique(rbind(fia.aa[,-c("uid","PLOT","UNITCD","COUNTYCD","STATECD")],fia.naa.t[,-c("uid","PLOT","UNITCD","COUNTYCD","STATECD")]))

#Retain one point per plot; retain most recent year surveyed per plot.
#If focal sp present, retain most recent year the plot was surveyed & the focal sp was id at present.
#Pull out presences & most recent date surveyed per plot
fia_lookup<-read.csv("/Users/blaginh/Desktop/host_dl_data/fia_lookup.csv")
fia_lookup$sciname<-paste0(fia_lookup$Genus," ", fia_lookup$Species)
fia_cd<-as.integer(fia_lookup$FIA[grep(sn,fia_lookup$sciname)])
fia$PLT_CN<-as.factor(fia$PLT_CN)
fia.aa<-fia[SPCD==fia_cd,]
fia.aa<-unique(fia.aa[,.(p_a="present",date=max(MEASYEAR),db="FIA",sciname=sn), by=.(PLOT,UNITCD,COUNTYCD,STATECD,LAT,LON)])
fia.aa[,uid:=as.character(paste0(PLOT,UNITCD,COUNTYCD,STATECD))]
names(fia.aa)[5:6]<-c("lat","lon")

#Pull out absences & most recent date surveyed per plot
fia.naa<-unique(fia[SPCD!=fia_cd,.(p_a="absent",date=max(MEASYEAR),db="FIA",sciname=sn), by=.(PLOT,UNITCD,COUNTYCD,STATECD,LAT,LON)])
fia.naa[,uid:=as.character(paste0(PLOT,UNITCD,COUNTYCD,STATECD))]
names(fia.naa)[5:6]<-c("lat","lon")
aasites<-unique(fia.aa[,.(uid,lat,lon)])
fia.naa.t<-as.data.table(anti_join(fia.naa, aasites, by=c("uid","lat","lon")))#only keep absences per plot where a present *wasn't* recorded

fia.t<-unique(rbind(fia.aa[,-c("uid","PLOT","UNITCD","COUNTYCD","STATECD")],fia.naa.t[,-c("uid","PLOT","UNITCD","COUNTYCD","STATECD")]))
fia.t<-fia.t[,.(year=date,p_a,sciname=sn,lat, lon,db)]

# NEON####
## Download using neon package
neon.out<-loadByProduct(dpID="DP1.10058.001", site="all", nCores = (detectCores()-1))
dir.create(ffwd("neon"))
write.csv(neon.out$div_10m2Data100m2Data,ffwd("neon/neon.div10m2and100m2.csv"))#trees & shrubs
write.csv(neon.out$div_1m2Data,ffwd("neon/neon.div_1m2.csv"))#shorter plants (e.g., ferns, forbs, grasses)

## Clean ####
neon<-as.data.table(read.csv(ffwd("neon/neon.div10m2and100m2.csv")))# if host is trees & shrubs
#neon<-as.data.table(read.csv(ffwd("neon/neon.div_1m2.csv")))#if host is not tree or shrub

# #Retain all
# neon.aa<-neon[grep(sn, neon$scientificName, ignore.case = T)]
# unique(neon.aa$scientificName)#previous will pull all species, subsp., & var.; confirm records pulled 
# neon.aa.t<-neon.aa[,.(p_a="present",date=endDate,sciname=sn, lat=decimalLatitude, lon=decimalLongitude, db="neon", plotID, subplotID)]
# neon.aa.t<-unique(neon.aa.t[,.(p_a,year=max(date),lat,lon,db,sciname,plotID), by=.(plotID, subplotID)])
# 
# neon.naa<-neon[grep(sn, neon$scientificName, ignore.case = T, invert = T)]
# neon.naa.t<-unique(neon.naa[,.(p_a="absent",date=endDate,sciname=sn, lat=decimalLatitude, lon=decimalLongitude, db="neon", plotID, subplotID)])
# aasites<-unique(neon.aa[,.(plotID,subplotID)])
# neon.naa.t<-as.data.table(anti_join(neon.naa.t, aasites, by=c("plotID","subplotID")))
# neon.naa.t<-unique(neon.naa.t[,.(p_a,year=max(date),lat,lon,db,sciname,plotID), by=.(plotID, subplotID)])
# neon.t<-unique(rbind(neon.aa.t[,-c("plotID","subplotID")],neon.naa.t[,-c("plotID","subplotID")]))
# neon.t$year<-strptime(neon.t$year, '%Y-%m-%d')$year+1900

#Retain one point per plot; retain most recent year surveyed per plot.
#If focal sp present, retain most recent year the plot was surveyed & the focal sp was id at present.
#Pull out presences & most recent date surveyed per plot
neon10m100m.aa<-neon[grep(sn_code, neon$taxonID, ignore.case = T)]
neon10m100m.aa.t<-neon10m100m.aa[,.(p_a="present",date=endDate,sciname=taxonID, lat=decimalLatitude, lon=decimalLongitude, db="neon10m100m", plotID, subplotID)]
neon10m100m.aa.t<-unique(neon10m100m.aa.t[,.(p_a,date=max(date),lat,lon,db,sciname,plotID), by=.(plotID, subplotID)])

neon10m100m.naa<-neon[grep(sn_code, neon$taxonID, ignore.case = T, invert = T)]
neon10m100m.naa.t<-unique(neon10m100m.naa[,.(p_a="absent",date=endDate,sciname=sn_code, lat=decimalLatitude, lon=decimalLongitude, db="neon10m100m", plotID, subplotID)])
aasites<-unique(neon10m100m.aa[,.(plotID,subplotID)])
neon10m100m.naa.t<-as.data.table(anti_join(neon10m100m.naa.t, aasites, by=c("plotID","subplotID")))
neon10m100m.naa.t<-unique(neon10m100m.naa.t[,.(p_a,date=max(date),lat,lon,db,sciname,plotID), by=.(plotID, subplotID)])

neon10m100m.t<-unique(rbind(neon10m100m.aa.t[,-c("plotID","subplotID")],neon10m100m.naa.t[,-c("plotID","subplotID")]))

neon10m100m.t$date<-as.POSIXct(neon10m100m.t$date, '%Y-%m-%d')
neon10m100m.t$date<-format(neon10m100m.t$date, '%Y')
neon10m100m.t<-neon10m100m.t[,.(year=date,p_a,sciname=sn,lat, lon,db="neon")]

#BLM ####
## Download using online link, no api####
dir.create(ffwd("blm"))
download.file(url="https://gis.blm.gov/AIMDownload/LayerPackages/BLM_AIM_TerrADat.zip",destfile = ffwd("/blm/BLM_AIM_TerrADat.zip"))
unzip(zipfile = ffwd("/blm/BLM_AIM_TerrADat.zip"), exdir = ffwd("blm/"))

## Clean
# #Keep all
# blm<-readOGR(ffwd("blm/BLM_AIM_TerrADat/v107/TerrADat.gdb"))
# blmWGS84<-spTransform(blm,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# blm<-as.data.table(blmWGS84)
# blm.aa<-blm[grep("AIAL",blm$Spp_Nox)]#USDA Plants Code for Ailanthus altissima 
# blm.aa_t<-blm.aa[,.(date=DateVisited,p_a="present",lat=coords.x2,lon=coords.x1,db="blm", sciname=sn, PlotID=PlotID)]
# blm.aa_t<-unique(blm.aa_t[,.(p_a,year=max(date),lat,lon,db,sciname,PlotID), by=.(PlotID)])
# 
# blm.naa<-blm[!grep("AIAL",blm$Spp_Nox)]
# aasites<-as.data.table(blm.aa_t[,.(PlotID,lat,lon)])
# blm.naa<-as.data.table(anti_join(blm.naa, aasites, by="PlotID"))
# blm.naa_t<-unique(blm.naa[,.(p_a="absent",year=max(DateVisited),lat=coords.x2,lon=coords.x1,db="blm",sciname=sn,PlotID), by=.(PlotID)])
# blm_t<-unique(rbind(blm.aa_t[,-("PlotID")],blm.naa_t[,-("PlotID")]))

#Remove dups

blm<-readOGR(ffwd("blm/BLM_AIM_TerrADat/v107/TerrADat.gdb"))
blmWGS84<-spTransform(blm,"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
blm<-as.data.table(blmWGS84)
blm.aa<-blm[grep(sn_code,blm)]
blm.aa_t<-blm.aa[,.(date=DateVisited,p_a="present",lat=coords.x2,lon=coords.x1,db="blm", sciname="AIAL", PlotID=PlotID)]
blm.aa_t<-unique(blm.aa_t[,.(p_a,date=max(date),lat,lon,db,sciname,PlotID), by=.(PlotID)])

blm.naa<-blm[!grep("JUNI",blm$Spp_Nox)]
aasites<-as.data.table(blm.aa_t[,.(PlotID,lat,lon)])
blm.naa<-as.data.table(anti_join(blm.naa, aasites, by="PlotID"))
blm.naa_t<-unique(blm.naa[,.(p_a="absent",date=max(DateVisited),lat=coords.x2,lon=coords.x1,db="blm",sciname="AIAL",PlotID), by=.(PlotID)])
blm_t<-unique(rbind(blm.aa_t[,-("PlotID")],blm.naa_t[,-("PlotID")]))

blm_t$date<-as.POSIXct(blm_t$date, '%Y/%m/%d')
blm_t$date<-format(blm_t$date, '%Y')
blm_t<-blm_t[,.(lat,lon,year=date,sciname=sn,p_a,db)]

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

#BIEN####
bien.out <- BIEN_occurrence_species(species = sn)
dir.create(ffwd("bien"))
write.csv(bien.out,ffwd("bien/bien.raw.csv"))

#Clean####
bien<-as.data.table(read.csv(ffwd("bien/bien.raw.csv")))
#unique(bien$datasource) #remove db duplicates - gbif, inaturalist, fia, neon, etc.
bien<-bien[grep("gbif", bien$datasource, ignore.case=T, invert=T),]#only gbif present, so gbif removed
bien_t<-bien[ ,.(year=date_collected, p_a="present",sciname=sn,lat=latitude,lon=longitude, db="bien")]
bien_t$year<-as.Date.character(bien_t$year)
bien_t$year<-format(bien_t$year,format="%Y")

#Combine data, convert to pts, crop to ref area, export cropped pts as geopackage & csv####
aa.occ<-rbind(gbif_t,eddmaps_t,fia.t,neon10m100m.t,blm_t,inat_t,bien_t)

#Exclude absences where presence is recorded at same xy location####
aa.occ.p<-aa.occ[p_a=="present",]#presences
aa.occ.a<-aa.occ[p_a!="present",]#absences
aa.occ.nodups.a<-as.data.table(anti_join(aa.occ.a, aa.occ.p, by=c("lat","lon")))#remove absences at xy locs where presence is recorded
aa.occ.nodups<-rbind(aa.occ.nodups.a,aa.occ.p)#combine xy locs w/ absence-only points with presence points


aa.pts <- st_as_sf(x = aa.occ.nodups,
                        coords = c("lon", "lat"),
                        crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

aa.pts.c<-st_crop(aa.pts,ref_area_xy)
aa.pts.c.df<-as.data.table(aa.pts.c)
aa.pts.c.df<-cbind(aa.pts.c.df[,1:4],st_coordinates(aa.pts.c))
names(aa.pts.c.df)[5:6]<-c("lon","lat")
dir.create(ffwd("final_outputs"))
dir.create(ffwd("final_outputs/allpts"))
write.csv(aa.pts.c.df,ffwd("final_outputs/aa.pts.csv"))
write_sf(aa.pts.c,ffwd("final_outputs/aa.pts.gpkg"),driver="GPKG")
aa.pts.c.df$p_a<-as.factor(aa.pts.c.df$p_a)
summary(aa.pts.c.df$p_a)
plot(aa.pts.c.df$p_a, main="Total TOH Occurrences")
dbs<-unique(aa.pts.c.df$db)
foreach(i=1:length(dbs)) %do% {
  x<-aa.pts.c.df[db == dbs[i],]
  plot(x$p_a, main=paste0("Total ",dbs[i], " Occurrences"))
  
}

#Remove year & retain unique instances
aa.occ.noyr<-unique(aa.occ.nodups[,.(p_a,sciname,lat,lon,dbs)])
aa.pts.noyr <- st_as_sf(x = aa.occ.noyr,
                   coords = c("lon", "lat"),
                   crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

aa.pts.c.noyr<-st_crop(aa.pts.noyr,ref_area_xy)
plot(aa.pts.c.noyr[1])
aa.pts.c.noyr.df<-as.data.table(aa.pts.c.noyr)
aa.pts.c.noyr.df<-cbind(aa.pts.c.noyr.df[,1:3],st_coordinates(aa.pts.c.noyr))
names(aa.pts.c.noyr.df)[4:5]<-c("lon","lat")
dir.create(ffwd("final_outputs/noyr"))
write.csv(aa.pts.c.noyr.df,ffwd("final_outputs/aa.pts.noyr.csv"))
write_sf(aa.pts.c.noyr,ffwd("final_outputs/aa.pts.noyr.gpkg"),driver="GPKG")
aa.pts.c.noyr.df$p_a<-as.factor(aa.pts.c.noyr.df$p_a)
summary(aa.pts.c.noyr.df$p_a)
plot(aa.pts.c.noyr.df$p_a, main="Total TOH Occurrences")
dbs<-unique(aa.pts.c.noyr.df$db)
foreach(i=1:length(dbs)) %do% {
  x<-aa.pts.c.noyr.df[dbs == dbs[i],]
  plot(x$p_a, main=paste0("Total ",dbs[i], " Occurrences"))
  
}

#Remove year & dbs, & retain unique instances
aa.occ.noyrdbs<-unique(aa.occ.nodups[,.(p_a,sciname,lat,lon)])
aa.pts.noyrdbs <- st_as_sf(x = aa.occ.noyrdbs,
                        coords = c("lon", "lat"),
                        crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

aa.pts.c.noyrdbs<-st_crop(aa.pts.noyrdbs,ref_area_xy)
plot(aa.pts.c.noyrdbs[1])
aa.pts.c.noyrdbs.df<-as.data.table(aa.pts.c.noyrdbs)
aa.pts.c.noyrdbs.df<-cbind(aa.pts.c.noyrdbs.df[,1:2],st_coordinates(aa.pts.c.noyrdbs))
names(aa.pts.c.noyrdbs.df)[3:4]<-c("lon","lat")

write.csv(aa.pts.c.noyrdbs.df,ffwd("final_outputs/aa.pts.noyrdbs.csv"))
write_sf(aa.pts.c.noyrdbs,ffwd("final_outputs/aa.pts.noyrdbs.gpkg"),driver="GPKG")

aa.pts.c.noyrdbs.df$p_a<-as.factor(aa.pts.c.noyrdbs.df$p_a)
summary(aa.pts.c.noyrdbs.df$p_a)

plot(aa.pts.c.noyrdbs.df$p_a, main="Total TOH Occurrences")


#To do ####
##1. add additional US/US regional databses: VMI NPS, BISON, Carolina Veg. Survey
##2. for inat - split ref area into grid when max results are reached, re-run, and combine results
##3. eddmap added manually - need to use api
##4. automatically find and append usda plants codes to sn
##5. for db that are downloaded in their entirety - fia, neon, blm, place in aws3 bucket, dc that latest version has been dl, if not, re-dl