## download iNaturalist Data for species of interest
library(spocc)
library(rinat)
species_of_interest <- "Ailanthus altissima"
extent <- c(24.396308,-124.848974, 49.384358, -66.885444) # bounding box for lower 48 states
inat_out <- get_inat_obs(taxon_name = species_of_interest, maxresults = 10000, geo = TRUE, bounds = extent, quality = "research")

## Download sentinenal 2 data for the study area
#devtools::install_github("16EAGLE/getSpatialData")
library(getSpatialData)
library(sf)
library(sp)

time_range <- c("2017-08-01","2018-01-01")
platform <- "Sentinel-2"
username <- "cmjone25"
password <- "I#3t8ert0t"
area_of_interest <- 

getSentinel_query(time_range = time_range, platform = platform, username = username,password = password, aoi = area_of_interest )