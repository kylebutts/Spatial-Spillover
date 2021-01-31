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



## Experimentation

df <- sim_data_misspecification(drop_geometry=FALSE) %>% 
	st_as_sf()

# Check treatment effects with optional center of population dots
ggplot(data = df %>% filter(year == 2019)) + 
	geom_sf(mapping = aes(fill = te_spill_decay, color = as.factor(treat))) + 
	# geom_sf(mapping = aes(geometry = centroid), color="white", size= 0.25) +
	scale_color_manual(values = c("black", "red"))



## Exerimentation

for(i in 1:10) {
	
df <- sim_data_misspecification(drop_geometry=FALSE) %>% 
	st_as_sf()
	
y <- "spill_decay"
x <- "spill_within_large"

# Donuts 80
formula <- as.formula(glue("y_{y} ~ treat_ind + post:spill_0_20 + spill_20_30:post + spill_30_40:post ", 
						   "+ spill_40_60:post + spill_60_80:post | state_county + year"))

reg <- feols(formula, data= df %>% st_drop_geometry()) 

coef <- reg %>% coefficients()
b_0_20 <- coef[2]
b_20_30 <- coef[3]
b_30_40 <- coef[4]
b_40_60 <- coef[5]
b_60_80 <- coef[6]

var <- as.character(glue("te_{y}"))

# Within 80
formula <- as.formula(glue("y_{y} ~ treat_ind + ({x}:post) | state_county + year"))

reg <- feols(formula, data= df %>% st_drop_geometry()) 
coef <- reg %>% coefficients()
b <- coef[2]

var <- as.character(glue("te_{y}"))

df <- df %>%
	mutate(
		te_spill_hat_within_large = b * !!rlang::sym(x),
		te_spill_hat_donut = b_0_20 * spill_0_20 + b_20_30 * spill_20_30 + b_30_40 * spill_30_40 + b_40_60 * spill_40_60 + b_60_80 * spill_60_80
	)

temp1 <- df %>% filter(year == 2019 & treat == 0) %>% pull(te_spill_decay)
temp2 <- df %>% filter(year == 2019 & treat == 0) %>% pull(te_spill_hat_within_large)
temp3 <- df %>% filter(year == 2019 & treat == 0) %>% pull(te_spill_hat_donut)

#hist(temp1 - temp2)
#hist(temp1 - temp3)

print(glue("MSPE within_large: {sum((temp1-temp2)^2)}"))
print(glue("MSPE donut: {sum((temp1-temp3)^2)}"))
}




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
# Helper function estimate
estimate_treatment_effect <- function(df, formula, treat_var){
	feols(formula, data= df) %>%
		coefficients() %>% .[[treat_var]]
}


n_trials <- 250
doParallel::registerDoParallel(10)

dgp_types <- c(
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within (Additive)" = "spill_within_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive"
)

estimation_types <- c(
	"TWFE (No Spillovers)" = "twfe",
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within 100mi." = "spill_within_100",
	"Within (Additive)" = "spill_within_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive",
	"Donuts (0-20, 20-30, 30-40)" = "donuts_small",
	"Donuts (0-20, 20-30, 30-40, 40-60, 60-80)" = "donuts"
)

results <- foreach(i = 1:n_trials, .combine = 'rbind') %dopar% {
	df <- sim_data_misspecification()
	results_trial <- NULL
	
	for(y in dgp_types) {
		for(x in estimation_types) {
			if(x == "twfe") {
				formula <- as.formula(glue("y_{y} ~ treat_ind | state_county + year"))
			} else if(x == "donuts_small") {
				formula <- as.formula(glue("y_{y} ~ treat_ind + spill_0_20:post + spill_20_30:post + spill_30_40:post | state_county + year"))
			} else if(x == "donuts") {
				formula <- as.formula(glue("y_{y} ~ treat_ind + spill_0_20:post + spill_20_30:post + spill_30_40:post ", 
										   "+ spill_40_60:post + spill_60_80:post | state_county + year"))
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
		columns = 2:6,
		decimals = 3
	) %>%
	# extract_body() %T>% cat(., file="tables/misspecification.tex") %>% cat()
	{.}





## Simulation: MSPE of Spillovers ----------------------------------------------


n_trials <- 250
doParallel::registerDoParallel(10)

dgp_types <- c(
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within (Additive)" = "spill_within_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive"
)

estimation_types <- c(
	"TWFE (No Spillovers)" = "twfe",
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within 100mi." = "spill_within_100",
	"Within (Additive)" = "spill_within_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive",
	"Donuts (0-20, 20-30, 30-40)" = "donuts_small",
	"Donuts (0-20, 20-30, 30-40, 40-60, 60-80)" = "donuts"
)

results <- foreach(i = 1:n_trials, .combine = 'rbind') %dopar% {
	df <- sim_data_misspecification()
	results_trial <- NULL
	
	for(y in dgp_types) {
		for(x in estimation_types) {
			if(x == "twfe") {
				df <- df %>% 
					mutate(te_spill_hat = 0)
			} else if(x == "donuts_small") {
				formula <- as.formula(glue("y_{y} ~ treat_ind + post:spill_0_20 + spill_20_30:post + spill_30_40:post | state_county + year"))
				
				reg <- feols(formula, data= df) 
				
				coef <- reg %>% coefficients()
				b_0_20 <- coef[2]
				b_20_30 <- coef[3]
				b_30_40 <- coef[4]
				
				df <- df %>% 
					mutate(
						te_spill_hat = b_0_20 * spill_0_20 + b_20_30 * spill_20_30 + b_30_40 * spill_30_40
					) 
				
			} else if(x == "donuts") {
				formula <- as.formula(glue("y_{y} ~ treat_ind + post:spill_0_20 + spill_20_30:post + spill_30_40:post ", 
										   "+ spill_40_60:post + spill_60_80:post | state_county + year"))
				
				reg <- feols(formula, data= df) 
				
				coef <- reg %>% coefficients()
				b_0_20 <- coef[2]
				b_20_30 <- coef[3]
				b_30_40 <- coef[4]
				b_40_60 <- coef[5]
				b_60_80 <- coef[6]
				
				df <- df %>% 
					mutate(
						te_spill_hat = b_0_20 * spill_0_20 + b_20_30 * spill_20_30 + b_30_40 * spill_30_40 + b_40_60 * spill_40_60 + b_60_80 * spill_60_80
					)
				
			} else {
				formula <- as.formula(glue("y_{y} ~ treat_ind + ({x}:post) | state_county + year"))
				
				reg <- feols(formula, data= df) 
				coef <- reg %>% coefficients()
				b <- coef[2]
				
				df <- df %>% 
					mutate(
						te_spill_hat = b * !!rlang::sym(x)
					)
			}
			
			# Calculate MSPE and normalize
			
			true_spill <- df %>% filter(year == 2019 & treat == 0) %>% pull(glue("te_{y}"))
			estimated_spill <- df %>% filter(year == 2019 & treat == 0) %>% pull("te_spill_hat")
			
			mspe <- sum((true_spill - estimated_spill)^2)
			total_var <- sum((true_spill^2))
			
			normalized_mspe <- mspe/total_var
			
			
			dgp <- names(dgp_types)[dgp_types == y]
			spec <- names(estimation_types)[estimation_types == x]
			
			results_trial <- bind_rows(results_trial, tibble(
				dgp = dgp, spec = spec, normalized_mspe = normalized_mspe
			))
		}
	}
	
	return(results_trial)
	
}


results_tbl <- results %>% 
	group_by(spec, dgp) %>% 
	summarize(
		# Makes 1 the best, 0 the worst
		percent_explained = 1 - mean(normalized_mspe)
	) %>% 
	ungroup() %>% 
	select(spec, dgp, percent_explained) %>%
	pivot_wider(names_from = dgp, values_from = percent_explained) %>%
	# Order columns
	select(spec, names(dgp_types)) %>%
	# Order rows
	arrange(match(spec, names(estimation_types)))



results_tbl %>% 
	gt::gt() %>%
	gt::fmt_number(
		columns = 2:6,
		decimals = 3
	) %>%
	#extract_body() %T>% cat(., file="tables/misspecification_mspe.tex") %>% cat()
	{.}

results_tbl %>% 
	gt::gt() %>%
	gt::fmt_percent(
		columns = 2:6,
		decimals = 1
	) %>%
	# extract_body() %T>% cat(., file="tables/misspecification_mspe_percent.tex") %>% cat()
	{.}

