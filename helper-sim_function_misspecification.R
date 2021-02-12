## -----------------------------------------------------------------------------
## helper-sim_function_misspecification.R
## Kyle Butts, CU Boulder Economics 
## 
## This file is going to create a function that will be used for monte carlo simulations. There will be many options to the function that can be used to test things. 
## -----------------------------------------------------------------------------

## Sim Function ----------------------------------------------------------------

#' Returns a dataframe generated with the given parameters to be used for diff-in-diff estimation
#'
#' @param treat_prob This is the unconditional probability of treatment
#' @param treat_spillover_distance The maximum distance in miles away a treated county can be to be considered treated 
#' @param treat_effect The true simulated treatment effect
#' @param treat_effect_spill The true simulated spillover effect.
#' @param treat_effect_spill_treat The true simulated spillover effect onto treated units.
#' @param spill_type Either "within", "contig", or "decay". Represents how spillover indicator is created
#' @param spatial_autocorr Should there be spatial autocorrelation of treament
#' @param zone_plus Value >= 0 and tells you how much to favor counties with values in the top 10% of the generated field
#' @param remove_spill Should control units with spillover effects be removed
#' 
#' @return a tibble of estimates of treatment effect under different scenarios. This will be used with dplyr::add_row
#'
sim_data_misspecification <- function(
	treat_prob = 0.05,
	treat_spillover_distance = 40,
	treat_spillover_distance_large = 80,
	treat_effect = 2,
	treat_effect_spill = 1,
	treat_effect_spill_treat = 0.5,
	spill_type = "contig",
	spatial_autocorr = FALSE,
	zone_plus = 0.4,
	drop_geometry = TRUE
){
	## Create Treatment and Spillover Variables --------------------------------
	
	# Logical matrix, =1 if within the maxmium distance for spillover
	within <- (dist_mi > 0 & dist_mi <= treat_spillover_distance)
	within_large <- (dist_mi > 0 & dist_mi <= treat_spillover_distance_large)
	within_100 <- (dist_mi > 0 & dist_mi <= 100)
	
	within_0_20 = (dist_mi > 0 & dist_mi <= 20)
	within_20_30 = (dist_mi > 20 & dist_mi <= 30)
	within_30_40 = (dist_mi > 30 & dist_mi <= 40)
	within_40_60 = (dist_mi > 40 & dist_mi <= 60)
	within_60_80 = (dist_mi > 60 & dist_mi <= 80)
	within_80_100 = (dist_mi > 60 & dist_mi <= 80)
	
	
	# Spatial Decay
	# Following https://fmwww.bc.edu/repec/bocode/s/spgen.pdf
	spatial_decay <- exp(-.02 * dist_mi) 
	
	# Cutoff by element-wise multiplication & no self-spillover
	spatial_decay <- spatial_decay * within_large
	

	# Treatment variable
	counties_treat <- counties %>%
		mutate(
			treat = as.numeric(runif(n()) <= treat_prob)
		)
	
	
	
	## Create Spillover Variables ----------------------------------------------
	
	counties_treat <- counties_treat %>%
		mutate(
			# Contiguous
			spill_contig = (treat == 0) * contiguous %*% treat,
			spill_contig = as.numeric(spill_contig > 0),
			spill_contig_treat = (treat == 1) * contiguous %*% treat,
			spill_contig_treat = as.numeric(spill_contig_treat > 0),
			# Within
			spill_within = (treat == 0) * within %*% treat,
			spill_within = as.numeric(spill_within > 0),
			spill_within_treat = (treat == 1) * within %*% treat,
			spill_within_treat = as.numeric(spill_within_treat > 0),
			# Within Large
			spill_within_large = (treat == 0) * within_large %*% treat,
			spill_within_large = as.numeric(spill_within_large > 0),
			spill_within_large_treat = (treat == 1) * within_large %*% treat,
			spill_within_large_treat = as.numeric(spill_within_large_treat > 0),
			# Within 100
			spill_within_100 = (treat == 0) * within_100 %*% treat,
			spill_within_100 = as.numeric(spill_within_100 > 0),
			spill_within_100_treat = (treat == 1) * within_100 %*% treat,
			spill_within_100_treat = as.numeric(spill_within_100_treat > 0),
			# Within Additive
			spill_within_additive = (treat == 0) * within %*% treat,
			spill_within_additive = as.numeric(spill_within_additive),
			spill_within_additive_treat = (treat == 1) * within %*% treat,
			spill_within_additive_treat = as.numeric(spill_within_additive_treat),
			# Within Large Additive
			spill_within_large_additive = (treat == 0) * within_large %*% treat,
			spill_within_large_additive = as.numeric(spill_within_large_additive),
			spill_within_large_additive_treat = (treat == 1) * within_large %*% treat,
			spill_within_large_additive_treat = as.numeric(spill_within_large_additive_treat),
			# Decay (Closest) w/ Cutoff Length 
			spill_decay = (treat == 0) * apply(spatial_decay * kronecker(matrix(1,n(),1), t(treat)), 1, max),
			spill_decay = as.numeric(spill_decay),
			spill_decay_treat = (treat == 1) * apply(spatial_decay * kronecker(matrix(1,n(),1), t(treat)), 1, max),
			spill_decay_treat = as.numeric(spill_decay_treat),
			# Decay Additive w/ Cutoff Length 
			spill_decay_additive = (treat == 0) * spatial_decay %*% treat,
			spill_decay_additive = as.numeric(spill_decay_additive),
			spill_decay_additive_treat = (treat == 1) * spatial_decay %*% treat,
			spill_decay_additive_treat = as.numeric(spill_decay_additive_treat),
			# Spillover Rings
			spill_0_20 = (treat == 0) * within_0_20 %*% treat,
			spill_0_20 = as.numeric(spill_0_20 > 0),
			spill_20_30 = (treat == 0) * within_20_30 %*% treat,
			spill_20_30 = as.numeric(spill_20_30 > 0),
			spill_30_40 = (treat == 0) * within_30_40 %*% treat,
			spill_30_40 = as.numeric(spill_30_40 > 0),
			spill_40_60 = (treat == 0) * within_40_60 %*% treat,
			spill_40_60 = as.numeric(spill_40_60 > 0),
			spill_60_80 = (treat == 0) * within_60_80 %*% treat,
			spill_60_80 = as.numeric(spill_60_80 > 0),
			spill_80_100 = (treat == 0) * within_80_100 %*% treat,
			spill_80_100 = as.numeric(spill_80_100 > 0),
			# Make rings only the closest
			spill_20_30 = as.numeric(spill_20_30 == 1 & (spill_0_20 == 0)), 
			spill_30_40 = as.numeric(spill_30_40 == 1 & (spill_20_30 == 0) & (spill_0_20 == 0)), 
			spill_40_60 = as.numeric(spill_40_60 == 1 & (spill_30_40 == 0) & (spill_20_30 == 0) & (spill_0_20 == 0)),
			spill_60_80 = as.numeric(spill_60_80 == 1 & (spill_40_60 == 0) & (spill_30_40 == 0) & (spill_20_30 == 0) & (spill_0_20 == 0)),
			spill_80_100 = as.numeric(spill_80_100 == 1 & (spill_60_80 == 0) & (spill_40_60 == 0) & (spill_30_40 == 0) & (spill_20_30 == 0) & (spill_0_20 == 0)),
			# Spillover Rings Additive
			spill_0_20_additive = (treat == 0) * within_0_20 %*% treat,
			spill_20_30_additive = (treat == 0) * within_20_30 %*% treat,
			spill_30_40_additive = (treat == 0) * within_30_40 %*% treat,
			spill_40_60_additive = (treat == 0) * within_40_60 %*% treat,
			spill_60_80_additive = (treat == 0) * within_60_80 %*% treat,
			spill_80_100_additive = (treat == 0) * within_80_100 %*% treat
		)
	
	## Spillover Variables Have Same Average Effect ----------------------------
	
	counties_treat <- counties_treat %>% 
		mutate(
			spill_contig = 1/4 * treat_effect_spill * spill_contig / mean(spill_contig),
			spill_within = 1/4 * treat_effect_spill * spill_within / mean(spill_within),
			spill_within_large = 1/4 * treat_effect_spill * spill_within_large / mean(spill_within_large),
			spill_within_additive = 1/4 * treat_effect_spill * spill_within_additive / mean(spill_within_additive),
			spill_within_large_additive = 1/4 * treat_effect_spill * spill_within_large_additive / mean(spill_within_large_additive),
			spill_decay = 1/4 * treat_effect_spill * spill_decay / mean(spill_decay),
			spill_decay_additive = 1/4 * treat_effect_spill * spill_decay_additive / mean(spill_decay_additive)
		)
	
	
	## Create Panel Data -------------------------------------------------------
	
	df <- counties_treat %>% 
		expand_grid(year = 2000:2019) %>%
		mutate(
			post = if_else(year >= 2010, 1, 0),
			treat_ind = treat * post,
			te = treat_effect * treat_ind,
			# Spill onto Controls
			te_spill_contig = treat_effect_spill * post * spill_contig,
			te_spill_within = treat_effect_spill * post * spill_within,
			te_spill_within_large = treat_effect_spill * post * spill_within_large,
			te_spill_within_additive = treat_effect_spill * post * spill_within_additive,
			te_spill_within_large_additive = treat_effect_spill * post * spill_within_large_additive,
			te_spill_decay = treat_effect_spill * post * spill_decay,
			te_spill_decay_additive = treat_effect_spill * post * spill_decay_additive,
			# Spill onto Also Treated
			# te_spill_contig_treat = treat_effect_spill_treat * post * spill_contig_treat,
			# te_spill_within_treat = treat_effect_spill_treat * post * spill_within_treat,
			# te_spill_within_large_treat = treat_effect_spill_treat * post * spill_within_large_treat,
			# te_spill_within_additive_treat = treat_effect_spill_treat * post * spill_within_additive_treat,
			# te_spill_decay_treat = treat_effect_spill_treat * post * spill_decay_treat
			# te_spill_decay_additive_treat = treat_effect_spill_treat * post * spill_decay_additive_treat,	
		) %>%
		# Year FE
		group_by(year) %>%
		mutate(
			year_fe= rnorm(1, mean= (year-2000) * 0.2, sd= 0.1^2)
		) %>%
		ungroup() %>%
		group_by(state_county) %>%
		mutate(
			county_fe= rnorm(1, mean= 6, sd= 2^2)
		) %>%
		ungroup() %>%
		mutate(
			epsilon = rnorm(n(), mean= 0, sd= 2^2),
			y_spill_contig = -2 + year_fe + county_fe + te + te_spill_contig + epsilon,
			y_spill_within = -2 + year_fe + county_fe + te + te_spill_within + epsilon,
			y_spill_within_large = -2 + year_fe + county_fe + te + te_spill_within_large + epsilon,
			y_spill_within_additive = -2 + year_fe + county_fe + te + te_spill_within_additive + epsilon,
			y_spill_within_large_additive = -2 + year_fe + county_fe + te + te_spill_within_large_additive + epsilon,
			y_spill_decay = -2 + year_fe + county_fe + te + te_spill_decay + epsilon,
			y_spill_decay_additive = -2 + year_fe + county_fe + te + te_spill_decay_additive + epsilon,
		)
	
	if(drop_geometry) df <- df %>% select(-centroid, -geometry)
	
	return(df)
}




