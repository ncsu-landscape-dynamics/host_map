library(renv)
library(foreach)
library(rinat)

#GLOBAL: host data downloaders####
#GBIF####

#' @description function to get occurrences using a year or year and month query.
#' @param sn the scientific name of your target species
#' @param date the start date (yyyy-mm-dd) of the year you want to return
#' records for.
get_gbif_obs_ym <- function(sn, date, year, month, status,study.area,tax.key) {
  year_d <- lubridate::year(date)
  month_d <- lubridate::month(date)
  flds<-c(rgbif::occ_fields,"taxonomicStatus","acceptedScientificName")
    if(year==T & month==F){
      x <- rgbif::occ_search(taxonKey = tax.key,hasCoordinate = T, hasGeospatialIssue = F, 
                           geometry = study.area, year=year_d, limit=99999, occurrenceStatus = status,fields=flds)} 
    else if(month==T & year==F){
      x <- rgbif::occ_search(taxonKey = tax.key,hasCoordinate = T, hasGeospatialIssue = F, 
                           geometry = study.area, year=year_d, month=month_d,limit=99999, occurrenceStatus = status,fields=flds)
    }
}

#' @description this function takes the gbif outputs, cleans the data, and then
#' exports it to the folder as two csv files (1 for the raw data and 1 for the
#' cleaned data).
#' @param dt is a data.table containing the gbif output.

export_gbif_occs<-function(dt, tax.rank,tax.key){
  inat<-"28eb1a3f-1c15-4a95-931a-4af90ecb574d"
  dt<-dt[taxonomicStatus=="ACCEPTED" & publishingOrgKey!=inat,]
  dt$dateretrieved<-Sys.Date()
  if(tax.rank=="genus"){dt<-dt[genusKey==tax.key,]}else if(tax.rank=="species"){
    dt<-dt[speciesKey==tax.key,]}else{
      print("Warning: GBIF tax.rank NOT genus or species; this may affect host SDM data prep")}
  dt<-unique(dt)
  write.csv(dt,paste0(folder,"/", sn,"/gbif/",sn,".gbif.out.raw.csv"),row.names = F)
  dt.c<-unique(dt[,.(year=year,p_a=occurrenceStatus,
                                 lat=decimalLatitude, lon=decimalLongitude,
                                 db="gbif",dbname=acceptedScientificName,fkey=gbifID,
                                 retrieved=dateretrieved,species=species, below_sp=taxonRank)])
  dt.c$species<-gsub("× ", "×", dt.c$species,fixed=T)
  unique(dt.c$species)
  sp<-strsplit(dt.c$species," ")
  sp<-lapply(sp,function(x) cbind(genus=x[1],species=x[2]))
  sp<-data.table::as.data.table(do.call("rbind",sp))
  dt.c$species<-sp$species
  dt.c$genus<-sp$genus
  dt.c[below_sp!="GENUS"|below_sp!="SPECIES",]$below_sp<-"Y"
  dt.c[below_sp=="GENUS"|below_sp=="SPECIES",]$below_sp<-"N"
  write.csv(dt.c,paste0(folder,"/", sn,"/gbif/",sn,".gbif.out.clean.csv"),
            row.names = F)}

#' @description This function retrieves occurrences from GBIF using the rgbif package,
#' and saves the raw gbif data in a folder, and cleans the gbif data to prep it for
#' host mapping. Note: this function remove inaturalist records, and retains records
#' for the species-level or species and genus-level.
#' @param folder path to a folder to save the exported raw and cleaned host 
#' occurrence data. 
#' @param studyext a path to raster file that encompasses the area where you would
#' like to search for occurrence data. Enter "" if you do not want to filter
#' occurrences by spatial extent.
#' @param sn the scientific name of the focal host species or genus.
#' @param tax.rank the taxonomic rank of the *sn* provided. This can only be 
#' genus or species. If genus, then all species within the genus are retrieved.
#' @param date.range The date.range (e.g. c("2008-01-01","2022-12-01")) you'd like to request data for. Enter "" if you
#' do not want to filter occurrences by the date. For species already downloaded
#' to our GD data folder, the start date of the data.range will be updated to
#' reflect the last date represented in the metadata.

get_gbif_occs<-function(folder, studyext, sn, tax.rank,date.range){
  dir.create(paste0(folder,sn))
  dir.create(paste0(folder,sn,"/gbif"))
  if(nchar(studyext)>0){
    study.area<-terra::rast(studyext)
    study.area<-terra::project(study.area, 
                               "+proj=longlat +datum=WGS84 +no_defs +type=crs")
    study.area<-terra::ext(study.area)
    study.area<-rgbif::gbif_bbox2wkt(minx = study.area[1],miny =study.area[3],
                              maxx=study.area[2],maxy=study.area[4])}else{study.area<-gbif_bbox2wkt(minx = -179.9999999,miny=-89.9999999,
                                    maxx=179.9999999,maxy=89.9999999)}
  if(nchar(date.range[1])>0){date.range<-lubridate::ymd(date.range)}else{
    date.range<-lubridate::ymd("1600-01-01",Sys.Date())
  }
  tax.key<-rgbif::name_backbone(name=sn,rank=tax.rank)$usageKey
  n<-rgbif::occ_data(taxonKey = tax.key, hasCoordinate = T, hasGeospatialIssue = F,
                     geometry=study.area,year=paste0(lubridate::year(date.range[1]),",",
                                                     lubridate::year(date.range[2])),limit = 1)$meta$count
  flds<-c(rgbif::occ_fields,"taxonomicStatus","acceptedScientificName")
  
  if(n>0 & n<=99999){
    gbif.out<-rgbif::occ_search(taxonKey = tax.key, hasCoordinate = T, 
                                hasGeospatialIssue = F,geometry=study.area,
                                year=paste0(lubridate::year(date.range[1]),
                                            ",",lubridate::year(date.range[2])),
                                limit = 99999,fields=flds)
    gbif.out<-data.table::as.data.table(gbif.out$data)
    }else if (n>=99999){
      yr<-seq.Date(date.range[1],date.range[2], by="year")
      gbif.out<-foreach(i=1:(length(yr)),verbose=T) %do% {
        print(yr[i])
        x<-get_gbif_obs_ym(sn,date=yr[i],year=T, month=F,status="PRESENT", 
                           study.area,tax.key)
        n<-x$meta$count
        t1<-yr[i]
        Sys.sleep(3)
        if(n>99999){
          t1<-lubridate::ymd(paste0(lubridate::year(t1),"-01","-01"))
          t2<-lubridate::ymd(paste0(lubridate::year(t1)+1,"-01","-01"))
          mo<-seq.Date(t1,t2,by="month")
          mo<-mo[mo<=lubridate::ymd(date.range[2])]
          if(lubridate::year(tail(mo,n=1))==lubridate::year(date.range[2])){
            mo<-unique(c(mo,lubridate::ymd(date.range[2])))
            }
          x<-foreach(j=1:(length(mo)),verbose=T) %do% {
            print(mo[i])
            x<-get_gbif_obs_ym(sn,date=mo[j],month=T, year=F,status="PRESENT",
                               study.area,tax.key)
            data.table::as.data.table(x$data)
          }
          x<-data.table::rbindlist(x,fill=TRUE)
        }
        else if(n>0 & n<=99999){data.table::as.data.table(x$data)}}
    }else{}
  
  n.a<-rgbif::occ_data(taxonKey = tax.key, hasCoordinate = T, hasGeospatialIssue = F,
                     geometry=study.area,occurrenceStatus = "ABSENT",
                     year=paste0(lubridate::year(date.range[1]),",",
                                                     lubridate::year(date.range[2])),limit = 1)$meta$count
  if(n>0 & n.a>0){
    gbif.out.a<-rgbif::occ_search(taxonKey = tax.key, hasCoordinate = T, 
                                  hasGeospatialIssue = F,geometry=study.area,
                                  occurrenceStatus = "ABSENT",
                    year=paste0(lubridate::year(date.range[1]),",",
                                lubridate::year(date.range[2])),limit = 99999,
                    fields = flds)
    gbif.out.a<-data.table::as.data.table(gbif.out.a$data)
    if(n>=99999){gbif.out<-data.table::rbindlist(gbif.out,fill=TRUE)}else{gbif.out}
    gbif.out$occurrenceStatus<-1
    gbif.out.a$occurrenceStatus<-0
    gbif.out<-data.table::rbindlist(list(gbif.out,gbif.out.a),fill=TRUE)
    export_gbif_occs(dt=gbif.out,tax.rank,tax.key)
  }else if(n>0 & n.a==0){
    gbif.out<-data.table::rbindlist(gbif.out,fill=TRUE)
    gbif.out$occurrenceStatus<-1
    export_gbif_occs(dt=gbif.out,tax.rank,tax.key)
  }else if(n==0 & n.a>0){
    gbif.out.a<-rgbif::occ_search(taxonKey = tax.key, hasCoordinate = T, 
                                  hasGeospatialIssue = F,geometry=study.area,
                                  occurrenceStatus = "ABSENT",
                                  year=paste0(lubridate::year(date.range[1]),",",
                                              lubridate::year(date.range[2])),limit = 99999,
                                  fields = flds)
    gbif.out.a<-data.table::as.data.table(gbif.out.a$data)
    gbif.out.a$occurrenceStatus<-0
    export_gbif_occs(dt=gbif.out.a,tax.rank,tax.key)
  }
  else{print("No GBIF occurrences found")}}

#INAT####

#' @description function to get occurrences of species per year & month & date 
#' and ignore error when 0 results found.
#' @param sn the scientific name of your target species
#' @param date the start date (yyyy-mm-dd) of the year you want to return
#' records for.
#' @param max total number of results to return; must be between 1-10000
get_inat_obs_ymd <- function(sn, date, year, month, day, study.area) {
  year_d <- lubridate::year(date)
  month_d <- lubridate::month(date)
  day_d <- lubridate::day(date)
  tryCatch({
    if(year==T & month==F & day==F){
      x <- rinat::get_inat_obs(query = sn, maxresults = 9999, year = year_d, 
                               meta = TRUE, bounds=study.area, quality="research")} 
    else if(month==T & year==F & day==F){
      x <- rinat::get_inat_obs(query = sn, maxresults = 9999, year = year_d, 
                               month=month_d, meta = TRUE, bounds=study.area, quality="research")
    }
    else if(month==F & year==F & day==T){
      x <- rinat::get_inat_obs(query = sn, maxresults = 9999, year = year_d, 
                               month=month_d, day=day_d, meta = TRUE, bounds, 
                               quality="research")
    }
  }, warning = function(w) {
  }, error = function(e) {x$meta$found<-0
  })
}

#' @description This function retrieves occurrences from inaturalist using the rinat 
#' package, and saves the raw inat data in a folder, and cleans the inat data to 
#' prep it for host mapping.
#' @param folder path to a folder to save the exported raw and cleaned host 
#' occurrence data.
#' @param studyext a path to raster file that encompasses the area where you would
#' like to search for occurrence data. Enter "" if you do not want to filter
#' occurrences by spatial extent.
#' @param sn the scientific name of the focal host species or genus.
#' @param date.range The date.range (e.g. c("2008-01-01","2022-12-01")) you'd like 
#' to request data for. Enter "" if you do not want to filter occurrences by the date. 
#' For species already downloaded to our GD data folder, the start date of the 
#' data.range will be updated to reflect the last date represented in the metadata.

get_inat_occs<-function(folder, studyext, sn,date.range) {
  dir.create(paste0(folder,sn))
  dir.create(paste0(folder,sn,"/inat"))
  if(nchar(studyext)>0){
    study.area<-terra::rast(studyext)
    study.area<-terra::project(study.area, "+proj=longlat +datum=WGS84 +no_defs +type=crs")
    study.area<-terra::ext(study.area)
    study.area<-c(study.area[3],study.area[1],study.area[4],study.area[2])}else{study.area<-c(-90,-180,90,180)}
  if(nchar(date.range[1])>0){date.range<-lubridate::ymd(date.range)}else{
    date.range<-lubridate::ymd("2008-01-01",Sys.Date())
  }
  yr<-seq.Date(date.range[1], date.range[2], by="year")
  if(tail(lubridate::year(yr)[1])< lubridate::year(date.range[2])){yr<-c(yr[-1],date.range[2])}
  
  inat.out<-foreach(i=1:(length(yr)),.combine="rbind") %do% {
    x<-get_inat_obs_ymd(sn,date=yr[i],year=T, month=F,day=F,study.area)
    n<-x$meta$found
    t1<-yr[i]
    Sys.sleep(3)
    if(n>9999){
      t1<-lubridate::ymd(paste0(lubridate::year(t1),"-01","-01"))
      t2<-lubridate::ymd(paste0(lubridate::year(t1)+1,"-01","-01"))
      mo<-seq.Date(t1,t2,by="month")
      mo<-mo[mo<=lubridate::ymd(date.range[2])]
      if(lubridate::year(tail(mo,n=1))==lubridate::year(date.range[2])){
        mo<-unique(c(mo,lubridate::ymd(date.range[2])))}
      foreach(j=1:(length(mo)),.combine="rbind") %do% {
        x<-get_inat_obs_ymd(sn,date=mo[j],month=T, year=F,day=F,study.area)
        n<-x$meta$found
        t1<-mo[j]
        t2<-mo[j+1]
        Sys.sleep(3)
        if(n>9999){
          d<-seq.Date(t1,t2,by="day")
          foreach(k=1:(length(d)),.combine="rbind") %do% {
            x<-get_inat_obs_ymd(sn,date=d[k],day=T, year=F, month=F,study.area)
            n<-x$meta$found
            Sys.sleep(3)
            x$data}}else if(n<=9999){x$data}else{}}}else if(n<=9999){x$data}else{}}
  
  write.csv(inat.out,paste0(folder,sn,"/inat/",sn,".inat.out.raw.csv"), row.names = F)
  inat.out<-data.table::as.data.table(inat.out)
  inat.out<-inat.out[grep(sn, inat.out$scientific_name, ignore.case=T)]
  inat.out<-inat.out[geoprivacy!="obscured",]
  inat.out<-inat.out[taxon_geoprivacy!="obscured",]
  inat.out<-inat.out[quality_grade!="needs_id"]
  inat.out$datetime<-lubridate::ymd_hms(inat.out$datetime, tz="UTC")
  inat.out.c<-unique(inat.out[,.(year=lubridate::year(datetime),p_a=1,
                                 lat=latitude, lon=longitude,
                                 db="inat",dbname=scientific_name,fkey=id,
                                 retrieved=Sys.Date(),species=scientific_name)])
  inat.out.c$species<-gsub("× ", "×", inat.out.c$species,fixed=T)
  unique(inat.out.c$species)
  sp<-strsplit(inat.out.c$species," ")
  sp<-lapply(sp,function(x) cbind(genus=x[1],species=x[2], below_sp=x[3]))
  sp<-data.table::as.data.table(do.call("rbind",sp))
  inat.out.c$species<-sp$species
  inat.out.c$genus<-sp$genus
  inat.out.c$below_sp<-sp$below_sp
  inat.out.c[!is.na(below_sp),]$below_sp<-"Y"
  inat.out.c[is.na(below_sp),]$below_sp<-"N"
  write.csv(inat.out.c,paste0(folder,sn,"/inat/",sn,".inat.out.clean.csv"), row.names = F)}

#BIEN####
#' @description This function retrieves occurrences from BIEN using the BIEN 
#' package, and saves the raw BIEN data in a folder, and cleans the BIEN data to 
#' prep it for host mapping. FIA, inaturalist, and GBIF data will be removed prior
#' to exporting the raw and cleaned csvs. NOTE: BIEN does not allow the user to
#' filter records based on date. And therefore each call will download all
#' the records in the BIEN database for the spatial extent, sn, & 
#' @param folder path to a folder to save the exported raw and cleaned host 
#' occurrence data.
#' @param studyext a path to raster file that encompasses the area where you would
#' like to search for occurrence data. Enter "" if you do not want to filter
#' occurrences by spatial extent.
#' @param sn the scientific name of the focal host species or genus.
#' @param tax.rank the taxonomic rank of the *sn* provided. This can only be 
#' genus or species. If genus, then all species within the genus are retrieved.

get_bien_obs<-function(folder, studyext, sn, tax.rank){
  tax.rank<-tolower(tax.rank)
  dir.create(paste0(folder,sn))
  dir.create(paste0(folder,sn,"/bien"))
  sn<-strsplit(sn," ")
  if(length(sn[[1]])==1){sn<-stringr::str_to_title(sn[[1]][1])}else{sn<-paste0(stringr::str_to_title(sn[[1]][1]), " ",sn[[1]][-1])}
  if(nchar(studyext)>0){
    study.area<-terra::rast(studyext)
    study.area<-terra::project(study.area, "+proj=longlat +datum=WGS84 +no_defs +type=crs")
    study.area<-terra::ext(study.area)}else{study.area<-c(-180,180,-90,90)}
  if(tax.rank=='genus'){
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
  else{'tax.rank must be genus or species'}
  
  bien.out<-bien.out[datasource!='iNaturalist' & datasource!='GBIF' & datasource!= 'FIA',]
  bien.out<-bien.out[!is.na(latitude) | !is.na(longitude),]
  
  if(length(bien.out$scrubbed_species_binomial>1)){
    bien.out$retrieved<-Sys.Date()
    bien.out$uid<-paste0('bien', seq(1,length(bien.out$scrubbed_species_binomial),by=1))
    write.csv(bien.out,paste0(folder,sn,"/bien/",sn,".bien.out.raw.csv"), row.names = F)
    bien.out.c<-unique(bien.out[,.(year=lubridate::year(date_collected),p_a=1,
                                   lat=latitude,lon=longitude, db="bien",
                                   dbname=scrubbed_species_binomial,fkey=uid,
                                   retrieved=retrieved,species=scrubbed_species_binomial)])
    sp<-strsplit(bien.out.c$species," ")
    sp<-lapply(sp,function(x) cbind(genus=x[1],species=x[2], below_sp=x[3]))
    sp<-data.table::as.data.table(do.call("rbind",sp))
    bien.out.c$species<-sp$species
    bien.out.c$genus<-sp$genus
    bien.out.c$below_sp<-sp$below_sp
    bien.out.c[is.na(below_sp),]$below_sp<-"N"
    bien.out.c[!is.na(below_sp),]$below_sp<-"Y"
    write.csv(bien.out.c,paste0(folder,sn,"/bien/",sn,".bien.out.clean.csv"), row.names = F)}
  else if(length(bien.out$scrubbed_species_binomial==0)){print("No BIEN records found")}}
