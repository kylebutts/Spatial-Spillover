## -----------------------------------------------------------------------------
## 2020-08-22_bias_from_spatial_autocorr.R
## Kyle Butts, CU Boulder Economics 
## 
## Initial attempt at spatial autocorrelation of treatment variable with spillovers. Conditional on holding the number of units treated, does bias increase with the closeness of treated units. This answer will likely be different if treated units receive spillover effects or not.
## -----------------------------------------------------------------------------

library("tidyverse")
library("sf")
library("stars")
library("units")
library("gstat")
library("exactextractr")
library("fixest")
library("doParallel")

# Macbook
# setwd("~/Documents/Projects/Spatial Spillover/")
# Research Computing
setwd("/projects/kybu6659/spatial_spillover/")

# sim_data()
source("helper-sim_function.R")

# Load theme_kyle()
source("https://gist.githubusercontent.com/kylebutts/7dc66a01ec7e499faa90b4f1fd46ef9f/raw/15196997e5aad41696b03c49be3bfaaca132fdf2/theme_kyle.R")

# Export
export <- TRUE

## Load Spatial Data -----------------------------------------------------------
# From data-prepare_counties.R
load(file= "data/counties_and_mat.RData")



## Simulation Function ---------------------------------------------------------

#' Repeats n trials with specific parameters and returns a vector of te_hat
#'
#' @param n_trials the number of trials for each simulation
#' @param cl the number of cores to run this on
#' @return a tibble of estimates of treatment effect under different scenarios. This will be used with dplyr::add_row
treatment_effect_simulation <- function(
	# Simulation Parameters
	n_trials = 100, cl = 14,
	# sim_function parameters
	treat_prob = 0.1,
	treat_spillover_distance = 40,
	treat_effect = 2,
	treat_effect_spill = 1,
	treat_effect_spill_treat = 0.5,
	spill_type = "contig",
	spatial_autocorr = FALSE,
	zone_plus = 0.4,
	drop_geometry = TRUE
) {
	
	# Display what the function is doing
	usethis::ui_todo("Running simulation with parameters: (n_trials= {n_trials}, treat_prob= {treat_prob}, zone_plus = {zone_plus}).")
	cat("\n")
	
	doParallel::registerDoParallel(cl)
	results <- foreach(i = 1:n_trials, .combine = 'rbind') %dopar% {
		
		df <- sim_data(
			treat_prob = treat_prob,
			treat_spillover_distance = treat_spillover_distance,
			treat_effect = treat_effect,
			treat_effect_spill = treat_effect_spill,
			treat_effect_spill_treat = treat_effect_spill_treat,
			spill_type = spill_type,
			spatial_autocorr = spatial_autocorr,
			zone_plus = zone_plus,
			drop_geometry = drop_geometry
		)
		
		te_hat <- feols(y_spill_treat ~ treat_ind | state_county + year, 
						data= df) %>%
			coefficients() %>% .[["treat_ind"]]
		
		te_hat_drop_control <- feols(y_spill_treat ~ treat_ind | state_county + year, 
									 data= df %>% filter(spill != 1)) %>% 
			coefficients() %>% .[["treat_ind"]]
		
		te_hat_control_spill <- feols(y_spill_treat ~ treat_ind + spill_ind | state_county + year, 
									  data= df) %>%
			coefficients() %>% .[["treat_ind"]]
		
		te_hat_all_spill <- feols(y_spill_treat ~ treat_ind + spill_ind + spill_ind_treat | state_county + year, 
								  data= df) %>%
			coefficients() %>% .[["treat_ind"]]
		
		return(
			tibble(te_hat = te_hat, te_hat_drop_control = te_hat_drop_control, te_hat_control_spill = te_hat_control_spill, te_hat_all_spill = te_hat_all_spill)
		)
		
	}
	
	usethis::ui_done("Finished simulation.")
	return(results)
}		




## Run simulation with different paramter values --------------------------------

sims <- tribble(
	~n_trials, ~cl, ~treat_prob, ~treat_spillover_distance, ~treat_effect, ~treat_effect_spill, ~treat_effect_spill_treat, ~spill_type, ~spatial_autocorr, ~zone_plus, ~drop_geometry,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 0.0, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 0.2, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 0.4, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 0.6, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 0.8, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 1.0, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 1.2, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 1.4, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 1.6, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 1.8, TRUE,
	100, 14, 0.1, 40, 2, 1, -0.5, "within", TRUE, 2.0, TRUE
)

sim_results <- sims %>%
	rowwise() %>%
	mutate(
		ret= list(treatment_effect_simulation(
			# Simulation Parameters
			n_trials= n_trials, cl = cl, 
			# sim_data parameters
			treat_prob= treat_prob, 
			treat_spillover_distance= treat_spillover_distance,
			treat_effect= treat_effect,
			treat_effect_spill= treat_effect_spill,
			treat_effect_spill_treat= treat_effect_spill_treat,
			spill_type = spill_type, 
			spatial_autocorr = spatial_autocorr, zone_plus = zone_plus, 
			drop_geometry = drop_geometry
		))
	) %>% 
	unnest(cols = c("ret"))

if(export) save(sim_results, file = "data/sim-bias_from_spatial_autocorr.Rdata")


## Plot of Bias as a function of Treatment Probability -------------------------

# load("data/sim-bias_from_spatial_autocorr.Rdata")

bias <- sim_results %>% 
	group_by(zone_plus) %>% 
	# Bias from each simulation and standard error
	mutate(
		bias= mean(te_hat - 2), 
		bias_se= sd(te_hat - 2),
		bias_upper= bias + 1.96*bias_se,
		bias_lower= bias - 1.96*bias_se
	) %>%
	ungroup()

(bias_plot <- ggplot(bias, aes(x= zone_plus, y= bias)) +
 	geom_point(color = "#B48EAD") + geom_line(color = "#B48EAD") +
 	geom_ribbon(aes(ymin= bias_lower, ymax= bias_upper), fill = "#B48EAD", color = "white", alpha= 0.2) + 
 	xlim(0, 2.0) + 
	ylim(NA, 0) +
 	labs(
 		x= "Spatial Autocorrelation Measure", 
 		y= "Bias of Tau hat", 
 		title= "Bias as spatial correlation changes", 
 		subtitle= "Monte Carlo simulations with 100 trials for each simulation",
 		color= "Spillover on to Treated Units"
	) +
 	theme_kyle(slides = TRUE, title_pos = "left", has_subtitle = TRUE) + 
	theme(plot.title = ggplot2::element_text(size = 14))
)


if(export) ggsave("figures/figure-bias_from_spatial_autocorr.png", plot= bias_plot, 
	   dpi= 300, width= 2400/300, height= 1800/300, bg= "#ECECEC")


## With and Without removing Control units
bias_by_method <- sim_results %>% 
	rename(te_hat_keep_control = te_hat) %>% 
	pivot_longer(
		cols = starts_with("te_hat_"),
		names_to = "method",
		values_to = "te"
	) %>% 
	mutate(
		method = case_when(
			method == "te_hat_keep_control" ~ "Keep Contiguous Controls",
			method == "te_hat_drop_control" ~ "Drop Contiguous Controls",
			TRUE ~ "other"
		)
	) %>% 
	filter(method != "other") %>%
	group_by(zone_plus, method) %>% 
	summarize(
		bias = mean(te - treat_effect),
		se = sd(te), 
		bias_lower = quantile(te - treat_effect, 0.025),
		bias_upper = quantile(te - treat_effect, 0.975)
	)

(bias_fix_plot <- ggplot(bias_by_method, aes(x = zone_plus, y = bias)) +
		geom_point(aes(color = method)) + geom_line(aes(color = method)) +
		geom_ribbon(aes(ymin = bias_lower, ymax = bias_upper, fill = method), color = "white", alpha = 0.2) + 
		xlim(0, 2.0) + 
		ylim(NA, 0) +
		labs(
			x = "Spatial Autocorrelation Measure", 
			y = "Bias of Tau hat", 
			title = "Dropping control units no longer effectively removes all bias", 
			subtitle = "Monte Carlo simulations with 100 trials per simulation",
			color = "Estimation Strategy", fill = "Estimation Strategy"
		) +
		scale_fill_manual(values = c("#5E81AC", "#B48EAD")) +
		scale_color_manual(values = c("#5E81AC", "#B48EAD")) +
		theme_kyle(slides = TRUE, title_pos = "left", has_subtitle = TRUE) +
		# Put Legend on Bottom
		guides(col = guide_legend(title.position = "top", label.position = "bottom", nrow = 1)) +
		theme(
			plot.title = ggplot2::element_text(size = 14),
			legend.position = "bottom",
			legend.spacing.x = unit(5, "points")
		)
)

if(export) ggsave("figures/figure-bias_fix_spatial_autocorr.png", plot= bias_fix_plot, 
	   dpi= 300, width= 2400/300, height= 1800/300, bg= "#ECECEC")
