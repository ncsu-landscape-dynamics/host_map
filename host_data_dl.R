#GLOBAL: host data downloaders####
#GBIF####

#' @description function to retrieve GBIF occurrence data using rgbif::occ_download. NOTE: The
#' cleanest way to get around the API limits that works for small and large calls is
#' to use occ_download.Occ_download requires an env file with login information.
#' @param folder path to where your occurrence data will be downloaded to; if not already
#' created, two subfolders will be created within this path: a subfolder for the host species
#' called for rgbif::occ_data or *sn* and within this folder, another subfolder for the raw and
#' cleaned data from gbif. 
#' @param studyext a path to raster file that encompasses the area where you would
#' like to search for occurrence data
#' @param sn the scientific name of the host species or genus
#' @param tax.rank the taxonomic rank of the scientific name provided in *sn*; if genus, 
#' then all *species* occurrences within the extent and genus will be returned; if species, only
#' that species occurrences within the *studyext* will be returned. Tax.ranks above genus or below
#' species can be called using this function, but may not be properly "prepped" for SDM since cleaning
#' focuses on using a genus or species key (i.e.*tax.key*) to check the gbif data.
#' @param tax.key you must supply the taxonKey of the *sn* parameter. The *tax.key* 
#' can be identified by using rgbif::name_backbone(name=sn,rank=tax.rank). Use the 
#' usageKey from this output as the *tax.key*. 

gbif_occ_data<-function(folder, studyext, sn, tax.rank, tax.key){
  owd<-getwd()
  readRenviron(paste0(folder,".env"))#provide .env file containing GBIF login information
  dir.create(paste0(folder,sn))
  dir.create(paste0(folder,sn,"/gbif"))
  folderfun::setff("wd",paste0(folder,sn,"/gbif"))
  setwd(ffwd())
  study.area<-terra::rast(studyext)
  study.area<-terra::project(study.area, "+proj=longlat +datum=WGS84 +no_defs +type=crs")
  study.area<-terra::ext(study.area)
  if(study.area[1]<=-180){study.area[1]<--179.999999999}
  if(study.area[2]>=180){study.area[2]<-179.999999999}
  if(study.area[3]<=-90){study.area[3]<--89.999999999}
  if(study.area[4]>=90){study.area[4]<-89.999999999}
  study.area<-gbif_bbox2wkt(minx = study.area[1],miny =study.area[3],
                            maxx=study.area[2],maxy=study.area[4])
  n<-occ_data(taxonKey = tax.key, hasCoordinate = T, hasGeospatialIssue = F,
              geometry=study.area)$meta$count
  if(n>0){
    inat<-"28eb1a3f-1c15-4a95-931a-4af90ecb574d"
    status<-occ_download(pred("taxonKey", tax.key), pred("hasCoordinate", TRUE), 
                         pred_within(study.area),pred("hasGeospatialIssue", FALSE),
                         pred_not(pred("publishingOrg", inat)),
                         pred_not(pred("basisOfRecord", "MACHINE_OBSERVATION")))
    occ_download_wait(status[[1]])
    occ_download_get(status[[1]],)
    unzip(zipfile = paste0(ffwd(),"/",status[[1]],".zip"),
          files = "occurrence.txt",overwrite = T)
    gbif.out<-as.data.table(read.csv(file = paste0(ffwd(),"/","occurrence.txt"),
                                     sep='\t',header=T))
    gbif.out<-gbif.out[taxonomicStatus=="ACCEPTED" & taxonRank == "SPECIES",]
    gbif.out$dateretrieved<-Sys.Date()
    if(tax.rank=="genus"){gbif.out<-gbif.out[genusKey==tax.key,]}else if(tax.rank=="species"){
      gbif.out<-gbif.out[speciesKey==tax.key,]}else{
        print("Warning: GBIF tax.rank NOT genus or species; this may affect host SDM data prep")}
    gbif.out<-unique(gbif.out)
    gbif.out[occurrenceStatus=="PRESENT",]
    write.csv(gbif.out,paste0(ffwd(),"/", sn,".gbif.out.raw.csv"),row.names = F)
    gbif.out$occurrenceStatus<-ifelse(gbif.out$occurrenceStatus=="PRESENT",1,0)
    gbif.out.c<-unique(gbif.out[,.(year=year,p_a=occurrenceStatus,species=species,
                                   lat=decimalLatitude, lon=decimalLongitude,
                                   db="gbif",dbname=acceptedScientificName,fkey=gbifID,
                                   retrieved=dateretrieved)])
    sp<-strsplit(gbif.out.c$species," ")
    sp<-lapply(sp,function(x) cbind(genus=x[1],species=x[2]))
    sp<-as.data.table(do.call("rbind",sp))
    gbif.out.c$species<-sp$species
    gbif.out.c$genus<-sp$genus
    write.csv(gbif.out.c,paste0(ffwd(),"/",sn,".gbif.out.clean.csv"),
              row.names = F)}else{print("No GBIF occurrences found")}
  setwd(owd)}

#INAT####

#' @description function to get occurrences of species per year & month & date 
#' and ignore error when 0 results found.
#' @param sn the scientific name of your target species
#' @param date the start date (yyyy-mm-dd) of the year you want to return
#' records for.
#' @param max total number of results to return; must be between 1-10000
get_inat_obs_ymd <- function(sn, date, max, year, month, day) {
  year_d <- lubridate::year(date)
  month_d <- lubridate::month(date)
  day_d <- lubridate::day(date)
  tryCatch({
    if(year==T & month==F & day==F){
      x <- rinat::get_inat_obs(query = sn, maxresults = max, year = year_d, 
                               meta = TRUE, bounds=study.area, quality="research")} 
    else if(month==T & year==F & day==F){
      x <- rinat::get_inat_obs(query = sn, maxresults = max, year = year_d, 
                               month=month_d, meta = TRUE, bounds=study.area, quality="research")
    }
    else if(month==F & year==F & day==T){
      x <- rinat::get_inat_obs(query = sn, maxresults = max, year = year_d, 
                               month=month_d, day=day_d, meta = TRUE, bounds=study.area, 
                               quality="research")
    }
  }, warning = function(w) {
  }, error = function(e) {x$meta$found<-0
  })
}

#' @description this function will determine how to best split up naturalist 
#' api calls (temporally) to override the 10K record limit.Then, it will export
#' the combined inaturalist data from the api which will be exported as a csv. This
#' data will be cleaned to prep for host mapping. The cleaned data will also
#' be exported to a csv.
#' @param studyext a path to the raster file covering the study extent.
#' @param sn the scientific name of the host species or genus
#' @param date.range the start date & end dates to retrieve occurrence data.
#' @param max the number of records for the api to return; should be value
#' between 1-10000.
#' @param folder path to where your occurrence data will be downloaded to; if not already
#' created, two subfolders will be created within this path: a subfolder for the host species
#' called for inat::get_inat_obs or *sn* and within this folder, another subfolder for the raw and
#' cleaned data from inat.

get_inat_occ<-function(sn,studyext,date.range, max, folder) {
  study.area<-terra::rast(studyext)
  study.area<-terra::project(study.area, "+proj=longlat +datum=WGS84 +no_defs +type=crs")
  study.area<-terra::ext(study.area)
  study.area<-c(study.area[3],study.area[1],study.area[4],study.area[2])
  date.range<-lubridate::ymd(date.range)
  yr<-seq.Date(date.range[1], date.range[2], by="year")
  if(tail(lubridate::year(yr))[1]< lubridate::year(date.range[2])){yr<-c(yr[-1],date.range[2])}
  inat.out<-foreach(i=1:(length(yr)),.combine="rbind") %do% {
    x<-get_inat_obs_ymd(sn,date=yr[i],max=max,year=T, month=F,day=F)
    n<-x$meta$found
    t1<-yr[i]
    Sys.sleep(10)
    if(n>=10000){
      t1<-lubridate::ymd(paste0(lubridate::year(t1),"-01","-01"))
      t2<-lubridate::ymd(paste0(lubridate::year(t1)+1,"-01","-01"))
      mo<-seq.Date(t1,t2,by="month")
      mo<-mo[mo<=lubridate::ymd(date.range[2])]
      if(lubridate::year(tail(mo,n=1))==lubridate::year(date.range[2])){
        mo<-unique(c(mo,lubridate::ymd(date.range[2])))}
      foreach(j=1:(length(mo)),.combine="rbind") %do% {
        x<-get_inat_obs_ymd(sn,date=mo[j],max=max,month=T, year=F,day=F)
        n<-x$meta$found
        t1<-mo[j]
        t2<-mo[j+1]
        Sys.sleep(10)
        if(n>=10000){
          d<-seq.Date(t1,t2,by="day")
          foreach(k=1:(length(d)),.combine="rbind") %do% {
            x<-get_inat_obs_ymd(sn,date=d[k], max=max,day=T, year=F, month=F)
            n<-x$meta$found
            Sys.sleep(3)
            x$data}}else if(n<10000){x$data}else{}}}else if(n<10000){x$data}else{}}
  dir.create(paste0(folder,sn))
  dir.create(paste0(folder,sn,"/inaturalist"))
  write.csv(inat.out,paste0(folder,sn,"/inaturalist/",sn,".inat.out.raw.csv"), row.names = F)
  inat.out<-data.table::as.data.table(inat.out)
  inat.out<-inat.out[grep(sn, inat.out$scientific_name, ignore.case=T)]
  inat.out<-inat.out[geoprivacy!="obscured",]
  inat.out<-inat.out[taxon_geoprivacy!="obscured",]
  inat.out<-inat.out[quality_grade!="needs_id"]
  inat.out$datetime<-lubridate::ymd_hms(inat.out$datetime, tz="UTC")
  inat.out.c<-unique(inat.out[,.(year=lubridate::year(datetime),p_a=1,species=scientific_name,
                                 lat=latitude, lon=longitude,
                                 db="inat",dbname=scientific_name,fkey=id,
                                 retrieved=Sys.Date())])
  inat.out.c$species<-gsub("× ", "×", inat.out.c$species,fixed=T)
  unique(inat.out.c$species)
  sp<-strsplit(inat.out.c$species," ")
  sp<-lapply(sp,function(x) cbind(genus=x[1],species=x[2], belowspecies=x[3]))
  sp<-data.table::as.data.table(do.call("rbind",sp))
  inat.out.c$species<-sp$species
  inat.out.c$genus<-sp$genus
  inat.out.c$belowspecies<-sp$belowspecies
  write.csv(inat.out.c,paste0(folder,sn,"/inaturalist/",sn,".inat.out.clean.csv"), row.names = F)}
#BIEN####
#' @description this function will use the BIEN package to download occurrences
#' for the BIEN database. The returned data will be exported as a csv and the data
#' will also be cleaned for use in host SDM. The cleaned data will also be exported
#' to an CSV.
#' @param studyext a path to the raster file covering the study extent. If left blank,
#' the search will be conducted globally
#' @param sn the scientific name of the host species or genus
#' @param folder path to where your occurrence data will be downloaded to; if not already
#' created, two subfolders will be created within this path: a subfolder for the host species
#' called for BIEN or *sn* and within this folder, another subfolder for the raw and
#' cleaned data from BIEN.
#' @param tax.rank the taxonomic rank of the scientific name provided in *sn*; if genus, 
#' then all *species* occurrences within the extent and genus will be returned; if species, only
#' that species occurrences within the *studyext* will be returned.

get_bien_obs<-function(studyext,sn,date.range,folder,tax.rank){
  sn<-str_to_title(sn)
  if(tax.rank=='genus' & nchar(studyext)>1){
    study.area<-terra::rast(studyext)
    study.area<-terra::project(study.area, "+proj=longlat +datum=WGS84 +no_defs +type=crs")
    study.area<-terra::ext(study.area)
    bien.out<-data.table::as.data.table(BIEN::BIEN_occurrence_box(
      min.lat=study.area[3], max.lat = study.area[4], min.long = study.area[1], 
      max.long=study.area[2], genus=sn,natives.only = F))}
  else if(tax.rank=='species' & nchar(studyext)>1){
    study.area<-terra::rast(studyext)
    study.area<-terra::project(study.area, "+proj=longlat +datum=WGS84 +no_defs +type=crs")
    study.area<-terra::ext(study.area)
    bien.out<-data.table::as.data.table(BIEN::BIEN_occurrence_box(
      min.lat=study.area[3], max.lat = study.area[4], min.long = study.area[1], 
      max.long=study.area[2], species=sn,natives.only = F))}
  else if(tax.rank=='genus' & nchar(studyext)<1){
    bien.out<-data.table::as.data.table(BIEN::BIEN_occurrence_box(
      min.lat=-90, max.lat = 90, min.long = -180, max.long=180, 
      genus=sn,natives.only = F))}
  else if(tax.rank=='species' & nchar(studyext)<1){
    bien.out<-data.table::as.data.table(BIEN::BIEN_occurrence_box(
      min.lat=-90, max.lat = 90, min.long = -180, max.long=180, species=sn, 
      natives.only = F))}
  else{'tax.rank must be genus or species'}
  
  bien.out<-bien.out[datasource!='iNaturalist' & datasource!='GBIF' & datasource!= 'FIA',]
  bien.out<-bien.out[!is.na(latitude) | !is.na(longitude),]
  bien.out$uid<-paste0('bien', seq(1,length(bien.out$scrubbed_species_binomial),by=1))
  bien.out$retrieved<-Sys.Date()
  dir.create(paste0(folder,sn))
  dir.create(paste0(folder,sn,"/bien"))
  write.csv(bien.out,paste0(folder,sn,"/bien/",sn,".bien.out.raw.csv"), row.names = F)
  bien.out.c<-unique(bien.out[,.(year=lubridate::year(date_collected),p_a=1,
                               species=scrubbed_species_binomial,lat=latitude,
                               lon=longitude, db="bien",
                               dbname=scrubbed_species_binomial,fkey=uid,
                               retrieved=retrieved)])
  sp<-strsplit(bien.out.c$species," ")
  sp<-lapply(sp,function(x) cbind(genus=x[1],species=x[2], belowspecies=x[3]))
  sp<-data.table::as.data.table(do.call("rbind",sp))
  bien.out.c$species<-sp$species
  bien.out.c$genus<-sp$genus
  bien.out.c$belowspecies<-sp$belowspecies
  write.csv(bien.out.c,paste0(folder,sn,"/bien/",sn,".bien.out.clean.csv"), row.names = F)}

#L48 USA: host data downloaders####
#' @description this function will pull occurrences from the current FIA SQlite
#' database housed in the Data folder in the Shared Google Drive. The occurrences
#' in the FIA database have been sampled multiple times. An occurrence will be marked
#' as a presence for the *sn* if that species has been reported at any time at a location
#' over the course of the FIA sampling.
#' @param studyext a path to the raster file covering the study extent. If left blank,
#' the search will be conducted globally
#' @param sn the scientific name of the host species or genus
#' @param folder path to where your occurrence data will be downloaded to; if not already
#' created, two subfolders will be created within this path: a subfolder for the host species
#' called for FIA or *sn* and within this folder, another subfolder for the raw and
#' cleaned data from FIA.
studyext<-"~/Google Drive/My Drive/PhD/RA/ecff/data/ecff_rasters/ecff_2017.tif"
folder<-"~/Google Drive/My Drive/PhD/RA/ecff/data/hosts/"
sn<-"Prunus"

study.area<-terra::rast(studyext)
study.area<-terra::project(study.area, "+proj=longlat +datum=WGS84 +no_defs +type=crs")
us<-terra::crop(terra::vect(spData::us_states),study.area)
us.dt<-data.table::as.data.table(cbind(abb=datasets::state.abb, NAME=datasets::state.name))
us<-merge(us, us.dt, by="NAME")

dir.create(paste0(folder,sn))
dir.create(paste0(folder,sn,"/fia"))
path_db<-"~/Google Drive/Shared drives/Data/Original/FIA/"
library(DBI)
library(FIESTA)

con <- dbConnect(RSQLite::SQLite(), paste0(path_db, "SQLite_FIADB_ENTIRE.db"))
dbListTables(con)

fia.ref<-data.table::as.data.table(dbGetQuery(con, 'SELECT * FROM REF_SPECIES'))
fia.plot<-data.table::as.data.table(dbGetQuery(con, 'SELECT * FROM PLOT LIMIT 5000'))
fia.tree<-data.table::as.data.table(dbGetQuery(con, 'SELECT * FROM TREE LIMIT 5000'))
test<-data.table::as.data.table(dbGetQuery(con, 'SELECT * FROM SUBPLOT LIMIT 5000'))
ri<-UTdat <- DBgetPlots(states = "CT",allyrs = T,istree = T)
ri$tabs$
