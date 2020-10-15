## -----------------------------------------------------------------------------
## sim-bias_as_treat_increases.R
## Kyle Butts, CU Boulder Economics 
## 
## This file will run a simulation of the bias arising from Diff-in-diffs in spatial situations where spillover effects occur. The end result is a plot of the magnitude of the bias occurring as the number of treated units, and hence the "corrupted" control units, increase. 
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
setwd("~/Documents/Projects/Spatial Spillover/")
# Research Computing
# setwd("/projects/kybu6659/spatial_spillover/")

# sim_data()
source("helper-sim_function.R")

# Load theme_kyle()
source("https://gist.githubusercontent.com/kylebutts/7dc66a01ec7e499faa90b4f1fd46ef9f/raw/15196997e5aad41696b03c49be3bfaaca132fdf2/theme_kyle.R")

# Export
export <- TRUE
slides <- FALSE


## Load Spatial Data -----------------------------------------------------------
# From data_prepare_counties.R
load(file= "data/counties_and_mat.RData")



## Simulation Function ---------------------------------------------------------

#' Repeats n trials with specific parameters and returns a vector of te_hat
#'
#' @param n_trials the number of trials for each simulation
#' @param treat_prob This is the unconditional probability of treatment
#' @param treat_spillover_distance The maximum distance in meters away a treated county can be to be considered treated 
#' @param treat_effect The true simulated treatment effect
#' @param treat_effect_spill The true simulated spillover effect. This effect is as of 2020-07-15 is non-additive so it won't be added multiple times if multiple treated units are nearby.
#' @return a tibble of estimates of treatment effect under different scenarios. This will be used with dplyr::add_row
#' 
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
	usethis::ui_todo("Running simulation with parameters: (n_trials= {n_trials}, treat_prob= {treat_prob}).")
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
		
		te_hat <- feols(y_spill ~ treat_ind | state_county + year, data= df) %>%
			coefficients() %>% .[["treat_ind"]]
		
		te_hat_drop_control <- feols(y_spill ~ treat_ind | state_county + year, 
										   data= df %>% filter(spill_ind != 1)) %>% 
			coefficients() %>% .[["treat_ind"]]
		
		te_hat_control_spill <- feols(y_spill ~ treat_ind + spill_ind | state_county + year, 
									  data= df) %>%
			coefficients() %>% .[["treat_ind"]]
		
		te_hat_all_spill <- feols(y_spill ~ treat_ind + spill_ind + spill_ind_treat | state_county + year, 
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
	100, 14,  0.03, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.06, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.09, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.12, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.15, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.18, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.21, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.24, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.27, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.30, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.33, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.36, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.39, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.42, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.45, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.48, 40, 2, 1, 0, "within", FALSE, 0, TRUE,
	100, 14,  0.50, 40, 2, 1, 0, "within", FALSE, 0, TRUE
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
	unnest(cols = c(ret))

# Export Simulations
# if(export) save(sim_results, file = "data/sim-bias_as_treat_increase.RData")


## Plot of Bias as a function of Treatment Probability -------------------------

# load("data/sim-bias_as_treat_increase.RData")

bias <- sim_results %>% 
	group_by(treat_prob) %>%
	# Bias from each simulation and standard error
	summarize(
		bias = mean(te_hat - 2), 
		bias_se = sd(te_hat - 2),
		bias_lower = quantile(te_hat - treat_effect, 0.025),
		bias_upper = quantile(te_hat - treat_effect, 0.975),
	) %>%
	ungroup()

bias_plot <- ggplot(bias, aes(x= treat_prob, y= bias)) +
	geom_point(color = "#B48EAD") + geom_line(color = "#B48EAD") +
	geom_ribbon(aes(ymin = bias_lower, ymax = bias_upper), fill = "#B48EAD", color = "white", alpha = 0.2) + 
	xlim(0, 0.5) +
	ylim(NA, 0) +
	labs(
		x = "Treatment Probability", 
		y = "Bias of Tau hat"
	) +
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 1)) 


if(slides) {
	bias_plot <- bias_plot + 
		labs(	
			title = "As more units are treated, the bias increases", 
			subtitle = "Monte Carlo simulations with 100 trials per simulation"
		) + 
		theme_kyle(slides = TRUE, title_pos = "left", has_subtitle = TRUE)
}
if(!slides) {
	bias_plot <- bias_plot +
		theme_kyle()
}

bias_plot

if(export & slides) ggsave("figures/figure-bias_from_prob_treat_slides.png", plot= bias_plot, 
	   dpi= 300, width= 2400/300, height= 1800/300, bg= "#ECECEC")

if(export & !slides) ggsave("figures/figure-bias_from_prob_treat.png", plot= bias_plot, 
						   dpi= 300, width= 2400/300, height= 1800/300, bg= "transparent")


## Plot Bias with and without removing spillover -------------------------------

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
	group_by(treat_prob, method) %>% 
	summarize(
		bias = mean(te - treat_effect),
		se = sd(te), 
		bias_lower = quantile(te - treat_effect, 0.025),
		bias_upper = quantile(te - treat_effect, 0.975)
	) %>% 
	ungroup()

bias_fix_plot <- ggplot(bias_by_method, aes(x = treat_prob, y = bias)) +
 	geom_point(aes(color = method)) + geom_line(aes(color = method)) +
 	geom_ribbon(aes(ymin = bias_lower, ymax = bias_upper, fill = method), color = "white", alpha = 0.2) + 
 	xlim(0, 0.5) + 
 	labs(
 		x = "Treatment Probability", 
 		y = "Bias of Tau hat",
 		color = "Estimation Strategy", fill = "Estimation Strategy"
 		# caption = "spillovers: contiguous, tau = 2, tau_spill,control = 1, tau_spill,treat = 0.5"
 	) +
	scale_fill_manual(values = c("#5E81AC", "#B48EAD")) +
	scale_color_manual(values = c("#5E81AC", "#B48EAD")) +
	# Put Legend on Bottom
	guides(col = guide_legend(title.position = "top", label.position = "bottom", nrow = 1))


if(slides) {
	bias_fix_plot <- bias_fix_plot + 
		labs(	
			title = "Dropping control units removes bias, but increases variance", 
			subtitle = "Monte Carlo simulations with 100 trials per simulation"
		) + 
		theme_kyle(slides = TRUE, title_pos = "left", has_subtitle = TRUE) + 
		theme(
			plot.title = ggplot2::element_text(size = 14),
			legend.position = "bottom",
			legend.spacing.x = unit(5, "points")
		)
}
if(!slides) {
	bias_fix_plot <- bias_fix_plot +
		theme_kyle() +
		theme(
			plot.title = ggplot2::element_text(size = 14),
			legend.position = "bottom",
			legend.spacing.x = unit(5, "points")
		)
}

bias_fix_plot



if(export & slides) ggsave("figures/figure-bias_fix_slides.png", plot= bias_fix_plot, 
	   dpi= 300, width= 2400/300, height= 1800/300, bg= "#ECECEC")

if(export & !slides) ggsave("figures/figure-bias_fix.png", plot= bias_fix_plot, 
		dpi= 300, width= 2400/300, height= 1800/300, bg= "transparent")
