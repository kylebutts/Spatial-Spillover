## -----------------------------------------------------------------------------
## tva-analysis-replicate.R
## Kyle Butts, CU Boulder Economics 
## 
## Replicate Table 3 of Klein and Moretti (2014) 
## -----------------------------------------------------------------------------

library(tidyverse)
library(haven)
library(sf)
library(Rcpp)
library(RcppArmadillo)


setwd("~/Documents/Projects/Spatial Spillover/")

# Load theme_kyle()
source("https://raw.githubusercontent.com/kylebutts/templates/master/ggplot_theme/theme_kyle.R")

# Helper functions
source("tva-analysis-helpers.R")
source("helper-conley.R")

# Export
export <- FALSE
keep_slides <- TRUE

## Load Data -------------------------------------------------------------------

df <- read_dta("data/tva/build.dta") %>% 
	mutate(
		fips = sprintf("%05d", fipstat * 1000 + fipcnty)
	)
 
counties <- sf::read_sf("data/counties.geojson")

# Create TVA shape
tva_fips <- scan("data/tva/tvacounties.txt", character())

tva <- counties %>% 
	filter(fips %in% tva_fips) %>% 
	summarize()
	
dist_mat <- st_distance(counties %>% st_point_on_surface(), tva) 
counties$dist_to_tva <- as.vector(units::drop_units(units::set_units(dist_mat, "mi")))


# sanity check
ggplot() +
	geom_sf(data = counties, aes(fill = dist_to_tva)) + 
	geom_sf(data = tva, color = "red", fill = NA) + 
	coord_sf(datum = NULL) + 
	scale_fill_distiller()


# lat-long of county centroids
temp <- counties %>% st_transform(4326) %>% st_point_on_surface() %>% st_coordinates()
counties$long <- temp[,1]
counties$lat  <- temp[,2]

# Create spillover variable donuts
df_reg <- df %>% 
	left_join(., counties %>% st_drop_geometry(), by = "fips") %>% 
	mutate(
		tva_0_100 = !tva & dist_to_tva >= 0 & dist_to_tva < 100,
		tva_100_250 = !tva & dist_to_tva >= 100 & dist_to_tva < 250,
		tva_250_500 = !tva & dist_to_tva >= 250 & dist_to_tva < 500,
		tva_50_100 = !tva & dist_to_tva >= 50 & dist_to_tva < 100,
		year = "1"
	)




# Regression -------------------------------------------------------------------

df_reg <- df_reg %>% 
	mutate( 
		D_lnpop = winsorize_x((lnpop2000 - lnpop40)/6, 0.01),  
		D_lnwage = winsorize_x((lnwage2000 - lnwage40)/6, 0.01),
		D_lnagr = winsorize_x((lnagr2000 - lnagr40)/6, 0.01),
		D_lnmanuf = winsorize_x((lnmanuf2000 - lnmanuf40)/6, 0.01),
		D_lnvfprod = winsorize_x((lnvfprod2000 - lnvfprod40)/6, 0.01),
		D_lnmedfaminc = winsorize_x((lnmedfaminc2000 - lnmedfaminc50)/5, 0.01),
		D_lnfaval = winsorize_x((lnfaval2000 - lnfaval40)/6, 0.01)
	)

dep_names <- c(
		"Population" = "D_lnpop", 
		"Average manufacturing wage" = "D_lnwage", 
		"Agricultural employment" = "D_lnagr", 
		"Manufacturing employment" = "D_lnmanuf", 
		"Value of farm production" = "D_lnvfprod", 
		"Median family income" = "D_lnmedfaminc", 
		"Median housing value" = "D_lnfaval"
	)

controls <- c("lnelevmax", "lnelevrang", "lnarea", "lnpop20", "lnpop20sq", "lnpop30", "lnpop30sq", "popdifsq", "agrshr20", "agrshr20sq", "agrshr30", "agrshr30sq", "manufshr20", "manufshr30", "lnwage20", "lnwage30", "lntwage30", "lnemp20", "lnemp30", "urbshare20", "urbshare30", "lnfaval20", "lnfaval30", "lnmedhsval30", "lnmedrnt30", "white20", "white20sq", "white30", "white30sq", "pctil20", "pctil30", "urate30", "fbshr20", "fbshr30")

# missing: c("nowage20", "nowage30", "notwage30", "PRADIO")


# Wide Table Version -----------------------------------------------------------

table_tex <- ""
table_tex_slides <- ""

for(y in dep_names) {
	outcome_name <- names(dep_names[dep_names == y])
	cli::cli_h2("Starting on var: {outcome_name}")
	
	
	temp <- df_reg %>% 
		# Drop NAs from data
		drop_na(!!c(controls, y, "tva_0_100"))
	
	formula_oaxaca <- as.formula(paste0(y, " ~ ", paste(controls, collapse = " + ")))
	formula_controls <- as.formula(paste0(y, " ~ tva + ", paste(controls, collapse = " + ")))
	formula_controls_spill <- as.formula(paste0(y, " ~ tva + tva_0_100 + tva_100_250 + tva_250_500 + ", paste(controls, collapse = " + ")))
	
	# Oaxcac-Blinder Estimate
	reg_oaxaca <- oaxaca_estimate(temp %>% filter(is.na(border_county)), formula_oaxaca, "tva", cluster = "FIPSTAT1")
	
	# Diff-in-Diff with Controls
	reg_controls <- fixest::feols(formula_controls, data = temp, demeaned = TRUE)
	
	X <- reg_controls$X_demeaned
	e <- reg_controls$residuals
	coords <- as.matrix(temp[, c("lat", "long")])
	time <- rep(1, nrow(coords))
	# 1 mi = 1.60934 km
	dist_cutoff <- 500 * 1.60934
	lag_cutoff <- 1
	
	# cov_controls <- conley_ses(X, e, coords, dist_cutoff, lag_cutoff, cores = 4)$Spatial
	
	# Diff-in-Diff with Spillovers
	reg_spill <- fixest::feols(formula_controls_spill, data = temp, demeaned = TRUE)
	
	X <- reg_spill$X_demeaned
	e <- reg_spill$residuals
	coords <- as.matrix(temp[, c("lat", "long")])
	time <- rep(1, nrow(coords))
	# 1 mi = 1.60934 km
	dist_cutoff <- 500 * 1.60934
	lag_cutoff <- 1
	
	# cov_spill <- conley_ses(X, e, coords, dist_cutoff, lag_cutoff, cores = 4, time = time, id = time)$Spatial
	
	
	keep_slides <- y %in% c("lnagr", "lnmanuf", "lnmedfaminc")
	
	# Create row
	row    <- ""
	row_se <- ""
	row_slides    <- ""
	row_slides_se <- ""
	
	# Outcome Variable
	row    <- paste0(row, str_pad(outcome_name, 28, "right"), "& ")
	row_se <- paste0(row_se, str_pad("", 28, "right"), "& " )
	row_slides    <- paste0(row_slides, str_pad(outcome_name, 28, "right"), "& ")
	row_slides_se <- paste0(row_slides_se, str_pad("", 28, "right"), "& " )
	
	# Oaxaca-Binder 
	pt     <- reg_oaxaca[["te"]]
	se     <- reg_oaxaca[["se"]]
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	# Diff-in-Diff Controls TVA
	pt     <- coef(reg_controls)[["tva"]]
	se     <- sqrt(clubSandwich::vcovCR(reg_controls, temp$FIPSTAT1, "CR1S")[["tva", "tva"]])
	# se     <- cov_controls[["tva", "tva"]]
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	row_slides    <- paste0(row_slides, pt_str, "& ")
	row_slides_se <- paste0(row_slides_se, se_str, "& ")
	
	# Spillover TVA
	pt     <- coef(reg_spill)[["tva"]]
	se     <- sqrt(clubSandwich::vcovCR(reg_spill, temp$FIPSTAT1, "CR1S")[["tva", "tva"]])
	# se     <- cov_spill[["tva", "tva"]]
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	row_slides    <- paste0(row_slides, pt_str, "& ")
	row_slides_se <- paste0(row_slides_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_0_100TRUE"]]
	se     <- sqrt(clubSandwich::vcovCR(reg_spill, temp$FIPSTAT1, "CR1S")[["tva_0_100TRUE", "tva_0_100TRUE"]])
	# se     <- cov_spill[["tva_0_100TRUE", "tva_0_100TRUE"]]
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	row_slides    <- paste0(row_slides, pt_str, "& ")
	row_slides_se <- paste0(row_slides_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_100_250TRUE"]]
	se     <- sqrt(clubSandwich::vcovCR(reg_spill, temp$FIPSTAT1, "CR1S")[["tva_100_250TRUE", "tva_100_250TRUE"]])
	# se     <- cov_spill[["tva_100_250TRUE", "tva_100_250TRUE"]]
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	row_slides    <- paste0(row_slides, pt_str, "& ")
	row_slides_se <- paste0(row_slides_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_250_500TRUE"]]
	se     <- sqrt(clubSandwich::vcovCR(reg_spill, temp$FIPSTAT1, "CR1S")[["tva_250_500TRUE", "tva_250_500TRUE"]])
	# se     <- cov_spill[["tva_250_500TRUE", "tva_250_500TRUE"]]
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "\\\\\n")
	row_se <- paste0(row_se, se_str, "\\\\\n")
	row_slides    <- paste0(row_slides, pt_str, "\\\\\n")
	row_slides_se <- paste0(row_slides_se, se_str, "\\\\\n")
	
	cli::cat_line(row, row_se)
	
	table_tex <- paste(table_tex, row, row_se)
	if(keep_slides) table_tex_slides <- paste(table_tex_slides, row_slides, row_slides_se)
}

cat(table_tex)
if(export) cat(table_tex, file = "tables/tva_replication.tex")
cat(table_tex_slides)
if(export) cat(table_tex_slides, file = "tables/tva_replication_slides.tex")




# Long Table Version
# 
# for(y in dep_names) {
# 	outcome_name <- names(dep_names[dep_names == y])
# 	cli::cli_h2("Starting on var: {outcome_name}")
# 	
# 	formula <- as.formula(paste0("D_", y, " ~ tva"))
# 	formula_controls <- as.formula(paste0("D_", y, " ~ tva + ", paste(controls, collapse = " + ")))
# 	formula_controls_spill <- as.formula(paste0("D_", y, " ~ tva + tva_0_100 + tva_100_250 + tva_250_500 + ", paste(controls, collapse = " + ")))
# 	
# 	reg <- lm(formula, data = df_reg)
# 	reg_controls <- lm(formula_controls, data = df_reg)
# 	reg_spill <- lm(formula_controls_spill, data = df_reg)
# 	
# 	(tex <- texreg::texreg(
# 		l = list(reg, reg_controls, reg_spill), 
# 		custom.coef.map = list(
# 			"tva"             = "TVA", 
# 			"tva_0_100TRUE"   = "TVA between 0 and 100mi",
# 			"tva_100_250TRUE" = "TVA between 100 and 250mi",
# 			"tva_250_500TRUE" = "TVA between 250 and 500mi"
# 		)
# 	))
# 	
# 	cli::cli_text("{cat(tex)}")
# }





