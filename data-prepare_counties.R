## -----------------------------------------------------------------------------
## 2020-09-01_prepare_counties.R
## Kyle Butts, CU Boulder Economics 
## 
## Downloads counties, simplifies geometries, create distance matrix and contiguous matrix, and saves it. This download will be used by all the other scripts
## -----------------------------------------------------------------------------


setwd("~/Documents/Projects/Spatial Spillover/")

library("tidyverse")
library("sf")
library("stars")
library("units")
library("tigris")
library("gstat")
library("exactextractr") 
library("fixest")
library("doParallel")  


## Prepare Spatial Data --------------------------------------------------------

# Load FIPS codes
data("fips_codes", package= "tigris")

# Run once
# sf::read_sf("data/2010_county/US_county_2010.shp") %>% rmapshaper::ms_simplify(keep = 0.02) %>% st_transform(st_crs(2163)) %>% sf::write_sf(., "data/2010_county.geojson")

# Load Counties Shapefile at low resolution
counties <- sf::read_sf("data/2010_county.geojson") %>% 
	sf::st_make_valid() %>% st_buffer(0) %>%
	select(state_code= STATEFP10, county_code= COUNTYFP10) %>%
	# Merge with fips_codes to get state and county names
	left_join(., fips_codes, by= c("state_code", "county_code")) %>%
	# Remove non-continental states
	filter(!state_name %in% c("Alaska", "Hawaii", "Puerto Rico", "U.S. Virgin Islands", "American Samoa", "Northern Mariana Islands", "Guam")) %>%
	mutate(state_county = paste(state_code, county_code, sep="-")) %>% 
	arrange(state_code, county_code)


# Load Counties Center of Popula"tion from NHGIS
counties_centpop <- sf::read_sf("data/2010_county_centpop/US_county_cenpop_2010.shp") %>%
	filter(!STNAME %in% c("Alaska", "Hawaii", "Puerto Rico", "U.S. Virgin Islands", "American Samoa", "Northern Mariana Islands", "Guam")) %>%
	select(state_code= STATEFP, county_code= COUNTYFP) %>%
	arrange(state_code, county_code) %>% 
	# Project to a prettier US project EPSG:2163
	st_transform(st_crs(2163))

## Create Distance Matrix ------------------------------------------------------


# st_crs(counties) will show you that units are in meters
dist <- st_distance(st_geometry(counties_centpop))

dist_mi <- dist %>%
	units::set_units(., "mi") %>% 
	units::drop_units()

contiguous <- st_touches(counties, counties, sparse= FALSE)


## Save Results ----------------------------------------------------------------

save(counties, dist_mi, contiguous, file = "data/counties_and_mat.RData")
