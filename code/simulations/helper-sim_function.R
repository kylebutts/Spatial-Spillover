## -----------------------------------------------------------------------------
## helper-sim_function.R
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
#' @param treat_effect_spill The true simulated spillover effect. This effect is non-additive
#' @param treat_effect_spill_treat The true simulated spillover effect onto treated units.
#' @param spill_type Either "within", "contig", or "decay". Represents how spillover indicator is created
#' @param spatial_autocorr Should there be spatial autocorrelation of treament
#' @param zone_plus Value >= 0 and tells you how much to favor counties with values in the top 10% of the generated field
#' @param remove_spill Should control units with spillover effects be removed
#' 
#' @return a tibble of estimates of treatment effect under different scenarios. This will be used with dplyr::add_row
#'
sim_data <- function(
	treat_prob = 0.1,
	treat_spillover_distance = 40,
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

	# Spatial Decay
	# Following https://fmwww.bc.edu/repec/bocode/s/spgen.pdf
	spatial_decay <- exp(-.01 * dist_mi)
	diag(spatial_decay) <- 0
	spatial_decay <- spatial_decay / rowSums(spatial_decay)
	
	# No Spatial Autocorrelation of Treatment Assignment
	if(!spatial_autocorr) {
		# Create Treatment and Spillover Variables
		counties_treat <- counties |>
			mutate(
				treat = sample(c(1L, 0L), size= n(), replace= TRUE, prob = c(treat_prob, 1-treat_prob)),
			)
	}
	
	# Spatial Autocorrelation of Treatment Assignment
	if(spatial_autocorr) {
		
		# Generate Spatially correlated treatment variable
		vgm_range <- 200000
		vgm_dummy <- gstat(
			formula= z ~ 1, 
			locations= ~ x + y, 
			dummy= TRUE, 
			beta= 1, 
			model= vgm(psill= 1, range= vgm_range, model= 'Exp'), 
			nmax= 20
		)
		
		# Prevent output
		capture.output(
			field <- st_as_stars(st_bbox(counties), nx = 750, ny = 250, values = 0) |>
				as_tibble() |>
				select(x, y) |>
				predict(vgm_dummy, newdata= ., nsim= 1) |>
				st_as_stars(crs= st_crs(counties)),
			file = "/dev/null"
		)
		
		
		# Counties Treat
		counties_treat <- 
			# Merge in simulation results
			cbind(counties, sim= exact_extract(as(field, "Raster"), counties, "mean")) |>
			mutate(
				# Treat in the top 1-treat_prob % of the field
				zone = sim > quantile(sim, 1-treat_prob),
				
				# `zone` has increased probability of treatment
				prob = 0.05 + zone_plus * zone,
				prob = prob * treat_prob/mean(prob),
				
				# Generate Treatment Variable
				treat = (runif(n()) <= prob),
				treat = as.numeric(treat)
			) 
	}
	
	
	
	## Create Spillover Variables ----------------------------------------------
	
	counties_treat <- counties_treat |>
		mutate(
			# Spill only on non-treated units for now
			spill_contig = (treat == 0) * contiguous %*% treat,
			spill_contig = as.numeric(spill_contig > 0),
			spill_contig_treat = (treat == 1) * contiguous %*% treat,
			spill_contig_treat = as.numeric(spill_contig_treat > 0),
			spill_within = (treat == 0) * within %*% treat,
			spill_within = as.numeric(spill_within > 0),
			spill_within_treat = (treat == 1) * within %*% treat,
			spill_within_treat = as.numeric(spill_within_treat > 0),
			spill_decay = (treat == 0) * spatial_decay %*% treat,
			spill_decay_treat = (treat == 1) * spatial_decay %*% treat,
		)

	if(spill_type == "contig") {
		counties_treat <- counties_treat |>
			mutate(
				spill = spill_contig,
				spill_treat = spill_contig_treat
			)
	}
	
	if(spill_type == "within") {
		counties_treat <- counties_treat |>
			mutate(
				spill = spill_within,
				spill_treat = spill_within_treat
			)
	}
	
	if(spill_type == "decay") {
		counties_treat <- counties_treat |>
			mutate(
				spill = spill_decay,
				spill_treat = spill_decay_treat
			)
	}

	
	
	## Create Panel Data -------------------------------------------------------

	df <- counties_treat |> 
			expand_grid(year = 2000:2019) |>
			mutate(
				post = if_else(year >= 2010, 1, 0),
				treat_ind = treat * post,
				spill_ind = spill * post, 
				spill_ind_treat = spill_treat * post, 
				te = treat_effect * treat_ind,
				te_spill = treat_effect_spill * spill_ind,
				te_spill_treat = treat_effect_spill_treat * spill_ind_treat
			) |>
			# Year FE
			group_by(year) |>
			mutate(
				year_fe= rnorm(1, mean= (year-2000) * 0.2, sd= 0.1^2)
			) |>
			ungroup() |>
			group_by(state_county) |>
			mutate(
				county_fe= rnorm(1, mean= 6, sd= 2^2)
			) |>
			ungroup() |>
			mutate(
				epsilon = rnorm(n(), mean= 0, sd= 2^2),
				y_spill = -2 + year_fe + county_fe + te + te_spill + epsilon,
				y_spill_treat = -2 + year_fe + county_fe + te + te_spill + te_spill_treat + epsilon,
			)
	
	if(drop_geometry) df <- df |> select(-centroid, -geometry)
	
	return(df)
}




