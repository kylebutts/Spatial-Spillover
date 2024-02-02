#' This file is going to create a function that will be used for monte carlo simulations. There will be many options to the function that can be used to test things.
#' 

#' Generate treatment and spillover effects
#' 
#' @param treat_prob This is the unconditional probability of treatment
#' @param treat_effect The true simulated treatment effect
#' @param treat_effect_spill The true simulated spillover onto control units
#' @param normalize Should the treatment effects be normalized across specifications
#' 
#' @return counties modified with treatment/spillover vars
gen_treat_and_spill <- function(
	counties, 
	treat_prob = 0.05,
	treat_effect = 2,
	treat_effect_spill = 1,
	normalize = TRUE
) {

	counties$treat = runif(nrow(counties)) <= treat_prob
	counties$control = !counties$treat
	counties$treat_effect = treat_effect

	# Closest
	counties$closest_treat = apply(dist_mi[, counties$treat], 1, min)

	# create treatment and spillover variables
  counties <- counties |>
    mutate(
      # Within
      spill_within = control & (closest_treat <= 40),
      # Within Large
      spill_within_large = control & (closest_treat <= 80),
      # Within 100
      spill_within_100 = control & (closest_treat <= 100),
      # Within Additive
      spill_within_additive = 
				control * rowSums(within[, treat]),
      # Within Large Additive
      spill_within_large_additive = 
				control * rowSums(within_large[, treat]),
      # Decay Additive w/ Cutoff Length
      spill_decay_additive = 
				control * rowSums(spatial_decay[, treat]),
      # Decay (Closest) w/ Cutoff Length
			# Following https://fmwww.bc.edu/repec/bocode/s/spgen.pdf
      spill_decay = 
				exp(-0.02 * closest_treat) * (closest_treat <= 80),
      # Spillover Rings
      spill_0_20 = control & between(closest_treat, 0, 20),
      spill_20_30 = control & between(closest_treat, 20, 30),
      spill_30_40 = control & between(closest_treat, 30, 40),
      spill_40_60 = control & between(closest_treat, 40, 60),
      spill_60_80 = control & between(closest_treat, 60, 80),
      spill_80_100 = control & between(closest_treat, 80, 100),
      # Spillover Rings Additive
      spill_0_20_additive = control * rowSums(within_0_20[, treat]),
      spill_20_30_additive = control * rowSums(within_20_30[, treat]),
      spill_30_40_additive = control * rowSums(within_30_40[, treat]),
      spill_40_60_additive = control * rowSums(within_40_60[, treat]),
      spill_60_80_additive = control * rowSums(within_60_80[, treat]),
      spill_80_100_additive = control * rowSums(within_80_100[, treat])
    )

  ## Spillover Variables Have Same Average Effect ----------------------------
  if (normalize) {
    counties <- counties |>
      mutate(
        spill_within = 1 / 4 * treat_effect_spill * spill_within / mean(spill_within),
        spill_within_large = 1 / 4 * treat_effect_spill * spill_within_large / mean(spill_within_large),
        spill_within_additive = 1 / 4 * treat_effect_spill * spill_within_additive / mean(spill_within_additive),
        spill_within_large_additive = 1 / 4 * treat_effect_spill * spill_within_large_additive / mean(spill_within_large_additive),
        spill_decay = 1 / 4 * treat_effect_spill * spill_decay / mean(spill_decay),
        spill_decay_additive = 1 / 4 * treat_effect_spill * spill_decay_additive / mean(spill_decay_additive)
      )
  }

	return(counties)
}

#' Returns a dataframe generated with the given parameters to be used for diff-in-diff estimation
#' @param treat_prob This is the unconditional probability of treatment
#' @param treat_effect The true simulated treatment effect
#' @param treat_effect_spill The true simulated spillover onto control units
#' @param normalize Should the treatment effects be normalized across specifications
#' 
#' @return a modified `df`
sim_data_misspecification <- function(
	counties,
	treat_prob = 0.05,
	treat_effect = 2,
	treat_effect_spill = 1,
	normalize = TRUE
) {
	
	# Assign treatment
	counties = gen_treat_and_spill(
		counties, 
		treat_prob = treat_prob,
		treat_effect = treat_effect,
		treat_effect_spill = treat_effect_spill,
		normalize = normalize
	)

	# Create panel
	df = bind_rows(
		lapply(2000:2019, function(y) { counties |> mutate(year = y) })
	)

	# Generate outcomes
  df <- df |>
  	# expand_grid(year = 2000:2019) |> 
    mutate(
      post = if_else(year >= 2010, 1, 0),
      treat_ind = treat * post,
      te = treat_effect * treat_ind,
      # Spill onto Controls
      te_spill_within = treat_effect_spill * post * spill_within,
      te_spill_within_large = treat_effect_spill * post * spill_within_large,
      te_spill_within_additive = treat_effect_spill * post * spill_within_additive,
      te_spill_within_large_additive = treat_effect_spill * post * spill_within_large_additive,
      te_spill_decay = treat_effect_spill * post * spill_decay,
      te_spill_decay_additive = treat_effect_spill * post * spill_decay_additive,
    ) |>
    mutate(
      year_fe = rnorm(1, mean = (year - 2000) * 0.2, sd = 0.1^2), 
			.by = year
    ) |>
    mutate(
      county_fe = rnorm(1, mean = 6, sd = 2^2), 
			.by = state_county
    ) |>
    mutate(
      epsilon = rnorm(n(), mean = 0, sd = 2^2),
      y_spill_within = -2 + year_fe + county_fe + te + te_spill_within + epsilon,
      y_spill_within_large = -2 + year_fe + county_fe + te + te_spill_within_large + epsilon,
      y_spill_within_additive = -2 + year_fe + county_fe + te + te_spill_within_additive + epsilon,
      y_spill_within_large_additive = -2 + year_fe + county_fe + te + te_spill_within_large_additive + epsilon,
      y_spill_decay = -2 + year_fe + county_fe + te + te_spill_decay + epsilon,
      y_spill_decay_additive = -2 + year_fe + county_fe + te + te_spill_decay_additive + epsilon,
    )

  return(df)
}
