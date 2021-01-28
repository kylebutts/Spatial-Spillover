## -----------------------------------------------------------------------------
## sim-misspecification_of_spillover.R
## Kyle Butts, CU Boulder Economics 
## 
## This file estimates spillover functions using misspecified spillover variables and compares the biases found.
## -----------------------------------------------------------------------------

library("tidyverse")
library("glue")
library("sf")
library("stars")
library("units")
library("gstat")
library("exactextractr")
library("fixest")
library("doParallel")
library("gt")
library("magrittr")

# Macbook
setwd("~/Documents/Projects/Spatial Spillover/")
# Research Computing
# setwd("/projects/kybu6659/spatial_spillover/")

# Load theme_kyle()
source("https://gist.githubusercontent.com/kylebutts/7dc66a01ec7e499faa90b4f1fd46ef9f/raw/15196997e5aad41696b03c49be3bfaaca132fdf2/theme_kyle.R")

# Export
export <- TRUE
slides <- FALSE

# extract_body()
extract_body <- function(gt) {
	if(inherits(gt, "gt_tbl")) gt <- as.character(gt::as_latex(gt))
	stringr::str_match(gt, "(?s)\\\\midrule\\n(.*)\\\\bottomrule")[[2]]
}



# sim_data_misspecification()
source("helper-sim_function_misspecification.R")
 

## Load Spatial Data -----------------------------------------------------------
# From data_prepare_counties.R
load(file= "data/counties_and_mat.RData")



# Helper function estimate
estimate_treatment_effect <- function(data, formula, treat_var){
	feols(formula, data= data) %>%
		coefficients() %>% .[[treat_var]]
}



## Experimentation

sim <- sim_data_misspecification(drop_geometry=FALSE) %>% 
	st_as_sf()

# Check treatment effects with optional center of population dots
ggplot(data = sim %>% filter(year == 2019)) + 
	geom_sf(mapping = aes(fill = te_spill_decay_additive, color = as.factor(treat))) + 
	# geom_sf(mapping = aes(geometry = centroid), color="white", size= 0.25) +
	scale_color_manual(values = c("black", "red"))

# sim %>% 
# 	filter(year == 2019) %>% 
# 	select(starts_with("te_")) %>% 
# 	summary()
# 
# sim %>% 
# 	filter(year == 2019) %>% 
# 	select(starts_with("spill_")) %>% 
# 	summary()





## Simulation: Bias ------------------------------------------------------------
n_trials <- 10
doParallel::registerDoParallel(5)

dgp_types <- c(
	"Contiguous" = "spill_contig", 
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within (Additive)" = "spill_within_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive"
)

estimation_types <- c(
	"TWFE (No Spillovers)" = "twfe",
	"Contiguous" = "spill_contig",
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within 100mi." = "spill_within_100",
	"Within (Additive)" = "spill_within_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive",
	"Donuts (0-20, 20-40)" = "donuts_small",
	"Donuts (0-20, 20-40, 40-60, 60-80)" = "donuts",
	"Donuts (0-20, 20-40, 40-60, 60-80, 80-100)" = "donuts_large"
)

results <- foreach(i = 1:n_trials, .combine = 'rbind') %dopar% {
	df <- sim_data_misspecification()
	results_trial <- NULL
	
	for(y in dgp_types) {
		for(x in estimation_types) {
			if(x == "twfe") {
				formula <- as.formula(glue("y_{y} ~ treat_ind | state_county + year"))
			} else if(x == "donuts_small") {
				formula <- as.formula(glue("y_{y} ~ treat_ind + spill_0_10:post + spill_10_20:post + spill_20_40:post | state_county + year"))
			} else if(x == "donuts") {
				formula <- as.formula(glue("y_{y} ~ treat_ind + spill_0_10:post + spill_10_20:post + spill_20_40:post ", 
										   "+ spill_40_60:post + spill_60_80:post | state_county + year"))
			} else if(x == "donuts_large") { 
				formula <- as.formula(glue("y_{y} ~ treat_ind + spill_0_10:post + spill_10_20:post + spill_20_40:post ", 
										   "+ spill_40_60:post + spill_60_80:post + spill_80_100:post | state_county + year"))
			} else {
				formula <- as.formula(glue("y_{y} ~ treat_ind + ({x}:post) | state_county + year"))
			}
			
			treat_var <- as.character(glue("treat_ind"))
			
			te_hat <- estimate_treatment_effect(df, formula, treat_var)
			
			dgp <- names(dgp_types)[dgp_types == y]
			spec <- names(estimation_types)[estimation_types == x]
			
			
			results_trial <- bind_rows(results_trial, tibble(
				dgp = dgp, spec = spec, te_hat = te_hat
			))
		}
	}
	
	return(results_trial)
		
}

## Save results
# save(results, file="data/sim_misspecification.RData")
# load("data/sim_misspecification.RData")

results_tbl <- results %>% 
	group_by(spec, dgp) %>% 
	summarize(
		te_hat = mean(te_hat)
	) %>% 
	ungroup() %>% 
	mutate(bias = 2 - te_hat) %>%
	select(spec, dgp, bias) %>%
	pivot_wider(names_from = dgp, values_from = bias) %>%
	# Order columns
	select(spec, names(dgp_types)) %>%
	# Order rows
	arrange(match(spec, names(estimation_types)))


results_tbl %>% 
	gt::gt() %>%
	gt::fmt_number(
		columns = 2:7,
		decimals = 3
	) %>%
	extract_body() %T>% cat(., file="tables/misspecification.tex") %>% cat()



