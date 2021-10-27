## -----------------------------------------------------------------------------
## sim-misspecification_of_spillover.R
## Kyle Butts, CU Boulder Economics 
## 
## This file performs a set of Monte Carlo simulations that estimate a set of 
## DGPs using various exposure mappings and estimates them using exposure 
## mappings and ring(s). The first simulation reports on average bias found in 
## the treatment effect estimate and the second simulation reports on the 
## ability to accurately predict spillover effects from fitted values.
## -----------------------------------------------------------------------------

library(tidyverse)
library(here)
library(glue)
library(sf)
library(stars)
library(fixest)
library(doParallel)
library(gt)
library(kfbmisc)

# Export
export <- TRUE
run <- FALSE

# extract_body()
extract_body <- function(gt) {
	if(inherits(gt, "gt_tbl")) gt <- as.character(gt::as_latex(gt))
	gt <- stringr::str_match(gt, "(?s)\\\\midrule\\n(.*)\\\\bottomrule")[[2]]
	stringr::str_remove(gt, " \\\\\\\\ \n$")
}


# sim_data_misspecification()
source(here::here("helper-sim_function_misspecification.R"))
 

# ---- Load Spatial Data -------------------------------------------------------
# From data_prepare_counties.R
load(file= here::here("data/counties_and_mat.RData"))


# ---- Simulation: Bias --------------------------------------------------------
# Simulates a set of DGPs with different kinds of exposure mappings
# For each DGP, estimates using many different parametric forms + ring(s) method
# Report on the bias of the treatment effect estimate

# Helper function estimate
estimate_treatment_effect <- function(df, formula, treat_var){
	feols(formula, data= df) %>%
		coefficients() %>% .[[treat_var]]
}


n_trials <- 1000
doParallel::registerDoParallel(10)

dgp_types <- c(
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Decay" = "spill_decay", 
	"Within 40mi. (Additive)" = "spill_within_additive", 
	"Within 80mi. (Additive)" = "spill_within_large_additive", 
	"Decay (Additive)" = "spill_decay_additive"
)

estimation_types <- c(
	"TWFE (No Spillovers)" = "twfe",
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within 40mi. (Additive)" = "spill_within_additive", 
	"Within 80mi. (Additive)" = "spill_within_large_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive",
	"Rings (0-20, 20-30, 30-40)" = "donuts_small",
	"Rings (0-20, 20-30, 30-40, 40-60, 60-80)" = "donuts",
	"Rings (0-20, 20-30, 30-40, 40-60, 60-80) (Additive)" = "donuts_additive"
)

if(run){
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
			} else if(x == "donuts_additive") {
				formula <- as.formula(glue("y_{y} ~ treat_ind + spill_0_20_additive:post + spill_20_30_additive:post + spill_30_40_additive:post ", 
										   "+ spill_40_60_additive:post + spill_60_80_additive:post | state_county + year"))
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
}

## Save results
if(run) save(results, file="data/sim_misspecification.RData")

## Load results if not running simulation
if(!run) load("data/sim_misspecification.RData")

results_tbl <- results %>% 
	group_by(spec, dgp) %>% 
	mutate(
		bias = 2 - te_hat
	) %>%
	summarize(
		# ci_lower = quantile(bias, probs = 0.05),
		# ci_upper = quantile(bias, probs = 0.95),
		mse = sum(bias^2)/n(),
		bias = mean(bias)
	) %>% 
	ungroup() %>% 
	select(spec, dgp, bias, mse) 


table_tex <- ""

for(spec in estimation_types) {
	
	# Create row
	row    <- ""
	row_se <- ""
	
	# Specification Name
	spec_name <- names(estimation_types)[estimation_types == spec]
	row    <- paste0(row, str_pad(spec_name, 55, "right"))
	row_se <- paste0(row_se, str_pad("", 55, "right"))
	
	
	for(dgp in dgp_types) {
		row <- paste0(row, "& ")
		row_se <- paste0(row_se, "& ")
		
		# DGP Names
		dgp_name <- names(dgp_types)[dgp_types == dgp]
		
		temp <- results_tbl %>% filter(spec == spec_name & dgp == dgp_name) 
		
		bias <- temp %>% pull(bias)
		bias <- round(bias, digits = 3) %>% format(., digits = 3, nsmall=3) %>% str_pad(., 8, side = "right")
		mse  <- temp %>% pull(mse)
		mse  <- paste0("[", format(round(mse, digits = 3), nsmall=3), "]") %>% str_pad(., 8, side = "right")
		
		row <- paste0(row, bias)
		row_se <- paste0(row_se, mse)
	}
	
	row <- paste0(row, "\\\\")
	if(spec != estimation_types[length(estimation_types)]){
		row_se <- paste0(row_se, "\\\\")
	}
	# cli::cat_line(row, "\n", row_se, "\n")
	
	table_tex <- paste(table_tex, row, row_se, sep="\n")
} 

cli::cli_h1("Misspecification Bias Results")
cat(table_tex)

if(export) cat(table_tex, file="tables/misspecification.tex")





# ---- Simulation: MSPE of Spillovers ------------------------------------------
# Simulates a set of DGPs with different kinds of exposure mappings
# For each DGP, estimates using many different parametric forms + ring(s) method
# Estimate spillover effect for control units and determine mspe of spillover

n_trials <- 1000
doParallel::registerDoParallel(10)

dgp_types <- c(
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Decay" = "spill_decay", 
	"Within 40mi. (Additive)" = "spill_within_additive", 
	"Within 80mi. (Additive)" = "spill_within_large_additive", 
	"Decay (Additive)" = "spill_decay_additive"
)

estimation_types <- c(
	"TWFE (No Spillovers)" = "twfe",
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within 40mi. (Additive)" = "spill_within_additive", 
	"Within 80mi. (Additive)" = "spill_within_large_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive",
	"Rings (0-20, 20-30, 30-40)" = "donuts_small",
	"Rings (0-20, 20-30, 30-40, 40-60, 60-80)" = "donuts",
	"Rings (0-20, 20-30, 30-40, 40-60, 60-80) (Additive)" = "donuts_additive"
)

if(run){
results_mspe <- foreach(i = 1:n_trials, .combine = 'rbind') %dopar% {
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
				
			} else if(x == "donuts_additive") {
				formula <- as.formula(glue("y_{y} ~ treat_ind + post:spill_0_20_additive + spill_20_30_additive:post + spill_30_40_additive:post ", 
										   "+ spill_40_60_additive:post + spill_60_80_additive:post | state_county + year"))
				
				reg <- feols(formula, data= df) 
				
				coef <- reg %>% coefficients()
				b_0_20 <- coef[2]
				b_20_30 <- coef[3]
				b_30_40 <- coef[4]
				b_40_60 <- coef[5]
				b_60_80 <- coef[6]
				
				df <- df %>% 
					mutate(
						te_spill_hat = b_0_20 * spill_0_20_additive + b_20_30 * spill_20_30_additive + b_30_40 * spill_30_40_additive + b_40_60 * spill_40_60_additive + b_60_80 * spill_60_80_additive
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
}

if(run) save(results_mspe, file=here::here("data/sim_misspecification_mspe.RData"))
if(!run) load(here::here("data/sim_misspecification_mspe.RData"))

results_tbl <- results_mspe %>% 
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
		columns = 2:7,
		decimals = 3
	) %>%
	extract_body() %T>% cat(., file=here::here("tables/misspecification_mspe.tex")) %>% cat()



results_tbl %>% 
	gt::gt() %>%
	gt::fmt_percent(
		columns = 2:7,
		decimals = 1
	) %>%
	extract_body() %T>% cat(., file=here::here("tables/misspecification_mspe_percent.tex")) %>% cat()

