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
export <- TRUE
slides <- FALSE

## Load Data -------------------------------------------------------------------

df <- read_dta("data/tva/build.dta") %>% 
	mutate(
		fips = sprintf("%05d", fipstat * 1000 + fipcnty)
	)
 
counties <- sf::read_sf("data/2010_county.geojson") %>% 
	mutate(fips = paste0(STATEFP10, COUNTYFP10)) %>% 
	select(state = STATEFP10, county = COUNTYFP10, fips) %>% 
	rmapshaper::ms_simplify(keep = 0.05)

# Create TVA shape
tva_fips <- scan("data/tva/tvacounties.txt", character())

tva <- counties %>% 
	filter(fips %in% tva_fips) %>% 
	summarize()
	
dist_mat <- st_distance(counties %>% st_point_on_surface(), tva) 
counties$dist_to_tva <- as.vector(units::drop_units(units::set_units(dist_mat, "mi")))


controls <- c("lnelevmax", "lnelevrang", "lnarea", "lnpop20", "lnpop20sq", "lnpop30", "lnpop30sq", "popdifsq", "agrshr20", "agrshr20sq", "agrshr30", "agrshr30sq", "manufshr20", "manufshr30", "lnwage20", "lnwage30", "lntwage30", "lnemp20", "lnemp30", "urbshare20", "urbshare30", "lnfaval20", "lnfaval30", "lnmedhsval30", "lnmedrnt30", "white20", "white20sq", "white30", "white30sq", "pctil20", "pctil30", "urate30", "fbshr20", "fbshr30", "PRADIO30", "nowage20dum", "nowage30dum", "notwage30dum")


## Create Spillover Variables --------------------------------------------------

# lat-long of county centroids
temp <- counties %>% st_transform(4326) %>% st_point_on_surface() %>% st_coordinates()
counties$long <- temp[,1]
counties$lat  <- temp[,2]

# Create spillover variable donuts
df_reg <- df %>% 
	left_join(., counties %>% st_drop_geometry(), by = "fips") %>% 
	mutate(
		tva_0_100 = !tva & dist_to_tva >= 0 & dist_to_tva < 100,
		tva_0_50 = !tva & dist_to_tva >= 0 & dist_to_tva < 50,
		tva_50_100 = !tva & dist_to_tva >= 50 & dist_to_tva < 100,
		tva_100_150 = !tva & dist_to_tva >= 100 & dist_to_tva < 150,
		tva_150_200 = !tva & dist_to_tva >= 150 & dist_to_tva < 200
	)


## Logit for subsample ---------------------------------------------------------

formula_logit <- as.formula(paste0("tva ~ ", paste(controls, collapse = " + ")))
tva_logit <- fixest::femlm(formula_logit, data = df_reg, family = "logit")
df_reg <- df_reg %>% 
	mutate(
		phat = predict(tva_logit, df_reg), 
		keep = phat > quantile(phat, probs = c(0.25), na.rm = TRUE)
	)


## Figure of counties and spillover variable -----------------------------------

fips_in_sample <- df_reg %>%
	drop_na(!!c(controls)) %>% 
	filter(keep == TRUE) %>% 
	pull(fips)
rings <- counties %>% 
	filter(fips %in% fips_in_sample) %>%
	mutate(
		spill = case_when(
			0 < dist_to_tva & dist_to_tva <= 50 ~ "0 to 50 miles",
			50 < dist_to_tva & dist_to_tva <= 100 ~ "50 to 100 miles",
			100 < dist_to_tva & dist_to_tva <= 150 ~ "100 to 150 miles",
			150 < dist_to_tva & dist_to_tva <= 200 ~ "150 to 200 miles",
		)
	) %>% 
	filter(!is.na(spill)) %>% 
	group_by(spill) %>% 
	summarize()

rings$spill <- factor(rings$spill, levels = c("0 to 50 miles", "50 to 100 miles", "100 to 150 miles", "150 to 200 miles"))

nord_palette <- c("#295080", "#BF616A", "#A3BE8C", "#D08770")
grey_palette <- c("grey30", "grey50", "grey70")

(spillover_map <- ggplot() +
		geom_sf(
			data = counties %>% filter(fips %in% fips_in_sample), 
			fill = NA, color = "grey40", size = 0.5
		) + 
		geom_sf(data = rings, aes(fill = spill), color = NA) + 
		# Outline of TVA
		geom_sf(data = tva, color = "Black", fill = NA, size = 1.1) + 
		coord_sf(datum = NA) + 
		theme_kyle(base_size = 14, slides = slides) +
		theme(
			legend.position = c(0.1, 0.15), 
			legend.background = element_rect(colour = "grey40")
		) +
		labs(fill = "Spillover Bin") + {
			if(slides) labs(title = "Effective Sample and Spillover Variables")
		} +
		scale_fill_manual(values = nord_palette, na.translate = FALSE)) 
		# scale_fill_manual(values = grey_palette, na.translate = FALSE)) 


if(export & !slides) ggsave("figures/figure-tva-sample.pdf", spillover_map, 
	   dpi= 300, width= 2400/300, height= 1350/300, bg= "white")

if(export & slides) ggsave("figures/figure-tva-sample_slides.pdf", spillover_map, 
						   dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")

if(export & !slides) ggsave("figures/figure-tva-sample.png", spillover_map, 
							dpi= 300, width= 2400/300, height= 1350/300, bg= "white")

if(export & slides) ggsave("figures/figure-tva-sample_slides.png", spillover_map, 
						   dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")




## Regression -------------------------------------------------------------------

df_reg <- df_reg %>% 
	mutate( 
		D_lnpop = winsorize_x((lnpop2000 - lnpop40)/6, 0.01),  
		D_lnwage = winsorize_x((lnwage2000 - lnwage40)/6, 0.01),
		D_lnagr = winsorize_x((lnagr2000 - lnagr40)/6, 0.01),
		D_lnmanuf = winsorize_x((lnmanuf2000 - lnmanuf40)/6, 0.01),
		D_lnvfprod = winsorize_x((lnvfprod2000 - lnvfprod40)/6, 0.01),
		D_lnmedfaminc = winsorize_x((lnmedfaminc2000 - lnmedfaminc50)/5, 0.01),
		D_lnfaval = winsorize_x((lnfaval2000 - lnfaval40)/6, 0.01),
		D_lnpop_short = winsorize_x((lnpop60 - lnpop40)/2, 0.01),  
		D_lnwage_short = winsorize_x((lnwage60 - lnwage40)/2, 0.01),
		D_lnagr_short = winsorize_x((lnagr60 - lnagr40)/2, 0.01),
		D_lnmanuf_short = winsorize_x((lnmanuf60 - lnmanuf40)/2, 0.01),
		D_lnvfprod_short = winsorize_x((lnvfprod60 - lnvfprod40)/2, 0.01),
		D_lnmedfaminc_short = winsorize_x((lnmedfaminc60 - lnmedfaminc50)/1, 0.01),
		D_lnfaval_short = winsorize_x((lnfaval60 - lnfaval40)/2, 0.01)
	)

controls <- c("lnelevmax", "lnelevrang", "lnarea", "lnpop20", "lnpop20sq", "lnpop30", "lnpop30sq", "popdifsq", "agrshr20", "agrshr20sq", "agrshr30", "agrshr30sq", "manufshr20", "manufshr30", "lnwage20", "lnwage30", "lntwage30", "lnemp20", "lnemp30", "urbshare20", "urbshare30", "lnfaval20", "lnfaval30", "lnmedhsval30", "lnmedrnt30", "white20", "white20sq", "white30", "white30sq", "pctil20", "pctil30", "urate30", "fbshr20", "fbshr30", "PRADIO30", "nowage20dum", "nowage30dum", "notwage30dum")


# 1940 - 2000 Version ----------------------------------------------------------


dep_names <- c(
	"Agricultural employment" = "D_lnagr", 
	"Manufacturing employment" = "D_lnmanuf"
)


table_tex <- ""

for(y in dep_names) {
	outcome_name <- names(dep_names[dep_names == y])
	cli::cli_h2("{outcome_name}")
	
	
	temp <- df_reg %>% 
		# Drop NAs from data
		drop_na(!!c(controls, y, "tva_0_100")) %>% 
		filter(keep == 1)
	
	# Diff-in-Diff with Controls -----------------------------------------------
	formula_controls <- as.formula(paste0(y, " ~ tva + ", paste(controls, collapse = " + ")))
	reg_controls <- fixest::feols(formula_controls, data = temp, demeaned = TRUE)
	
	X <- reg_controls$X_demeaned
	e <- reg_controls$residuals
	coords <- as.matrix(temp[, c("lat", "long")])
	# 1 mi = 1.60934 km
	dist_cutoff <- 200 * 1.60934
	
	cov_controls <- conley_ses(X, e, coords, dist_cutoff)$Spatial
	
	# Diff-in-Diff with Spillovers ---------------------------------------------
	formula_controls_spill <- as.formula(paste0(y, " ~ tva + tva_0_50 + tva_50_100 + tva_100_150 + tva_150_200 + ", paste(controls, collapse = " + ")))
	reg_spill <- fixest::feols(formula_controls_spill, data = temp, demeaned = TRUE)
	
	X <- reg_spill$X_demeaned
	e <- reg_spill$residuals
	coords <- as.matrix(temp[, c("lat", "long")])
	# 1 mi = 1.60934 km
	dist_cutoff <- 200 * 1.60934

	cov_spill <- conley_ses(X, e, coords, dist_cutoff)$Spatial
	
	
	# Create row
	row    <- ""
	row_se <- ""
	
	# Outcome Variable
	row    <- paste0(row, str_pad(outcome_name, 28, "right"), "& ")
	row_se <- paste0(row_se, str_pad("", 28, "right"), "& " )
	
	
	# Diff-in-Diff Controls TVA
	pt     <- coef(reg_controls)[["tva"]]
	se     <- sqrt(cov_controls[["tva", "tva"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	# Spillover TVA
	pt     <- coef(reg_spill)[["tva"]]
	se     <- sqrt(cov_spill[["tva", "tva"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_0_50TRUE"]]
	se     <- cov_spill[["tva_0_50TRUE", "tva_0_50TRUE"]]
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_50_100TRUE"]]
	se     <- sqrt(cov_spill[["tva_50_100TRUE", "tva_50_100TRUE"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_100_150TRUE"]]
	se     <- sqrt(cov_spill[["tva_100_150TRUE", "tva_100_150TRUE"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_150_200TRUE"]]
	se     <- sqrt(cov_spill[["tva_150_200TRUE", "tva_150_200TRUE"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "\\\\\n")
	if(y != dep_names[length(dep_names)]) {
		row_se <- paste0(row_se, se_str, "\\\\\n")
	} else {
		row_se <- paste0(row_se, se_str, "\n")
	}
	
	cli::cat_line(row, row_se)
	
	table_tex <- paste(table_tex, row, row_se)
}

cat(table_tex)
if(export) cat(table_tex, file = "tables/tva_replication.tex")




## 1940 - 1960 Version ---------------------------------------------------------

dep_names <- c(
	"Agricultural employment" = "D_lnagr_short", 
	"Manufacturing employment" = "D_lnmanuf_short"
)


table_tex <- ""

for(y in dep_names) {
	outcome_name <- names(dep_names[dep_names == y])
	cli::cli_h2("{outcome_name}")
	
	
	temp <- df_reg %>% 
		# Drop NAs from data
		drop_na(!!c(controls, y, "tva_0_100")) %>% 
		filter(keep == 1)

	# Diff-in-Diff with Controls -----------------------------------------------
	formula_controls <- as.formula(paste0(y, " ~ tva + ", paste(controls, collapse = " + ")))
	reg_controls <- fixest::feols(formula_controls, data = temp, demeaned = TRUE)
	
	X <- reg_controls$X_demeaned
	e <- reg_controls$residuals
	coords <- as.matrix(temp[, c("lat", "long")])
	# 1 mi = 1.60934 km
	dist_cutoff <- 200 * 1.60934
	cov_controls <- conley_ses(X, e, coords, dist_cutoff)$Spatial
	
	# Diff-in-Diff with Spillovers ---------------------------------------------
	formula_controls_spill <- as.formula(paste0(y, " ~ tva + tva_0_50 + tva_50_100 + tva_100_150 + tva_150_200 + ", paste(controls, collapse = " + ")))
	reg_spill <- fixest::feols(formula_controls_spill, data = temp, demeaned = TRUE)
	
	X <- reg_spill$X_demeaned
	e <- reg_spill$residuals
	coords <- as.matrix(temp[, c("lat", "long")])
	# 1 mi = 1.60934 km
	dist_cutoff <- 200 * 1.60934
	cov_spill <- conley_ses(X, e, coords, dist_cutoff)$Spatial
	
	
	# Create row
	row    <- ""
	row_se <- ""
	
	# Outcome Variable
	row    <- paste0(row, str_pad(outcome_name, 28, "right"), "& ")
	row_se <- paste0(row_se, str_pad("", 28, "right"), "& " )
	
	# Diff-in-Diff Controls TVA
	pt     <- coef(reg_controls)[["tva"]]
	se     <- sqrt(cov_controls[["tva", "tva"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	# Spillover TVA
	pt     <- coef(reg_spill)[["tva"]]
	se     <- sqrt(cov_spill[["tva", "tva"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_0_50TRUE"]]
	se     <- sqrt(cov_spill[["tva_0_50TRUE", "tva_0_50TRUE"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_50_100TRUE"]]
	se     <- sqrt(cov_spill[["tva_50_100TRUE", "tva_50_100TRUE"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_100_150TRUE"]]
	se     <- sqrt(cov_spill[["tva_100_150TRUE", "tva_100_150TRUE"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "& ")
	row_se <- paste0(row_se, se_str, "& ")
	
	pt     <- coef(reg_spill)[["tva_150_200TRUE"]]
	se     <- sqrt(cov_spill[["tva_150_200TRUE", "tva_150_200TRUE"]])
	pt_str <- str_pad(reg_format(pt, se), 15, "both")
	se_str <- str_pad(paste0("$(", sprintf("%0.4f", se), ")$"), 15, "both")
	row    <- paste0(row, pt_str, "\\\\\n")
	if(y != dep_names[length(dep_names)]) {
		row_se <- paste0(row_se, se_str, "\\\\\n")
	} else {
		row_se <- paste0(row_se, se_str, "\n")
	}
	
	cli::cat_line(row, row_se)
	
	table_tex <- paste(table_tex, row, row_se)
}

cat(table_tex)
if(export) cat(table_tex, file = "tables/tva_replication_short.tex")



