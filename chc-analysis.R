## -----------------------------------------------------------------------------
## chc_replicate.R
## Kyle Butts, CU Boulder Economics 
## 
## Replication of Bailey and Goodman-Bacon (2015) with extension for event study 
## framework.
## -----------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(cli)
library(did)
library(fixest)

source("https://raw.githubusercontent.com/kylebutts/templates/master/ggplot_theme/theme_kyle.R")

setwd("~/Documents/Projects/Spatial Spillover/")

slides <- FALSE
h_w_ratio <- 9/16


# Set-up Data ------------------------------------------------------------------

df <- 
	haven::read_dta("data/chc/aer_data.dta") %>% 
	# Remove LA, New York, Chicago
	filter(
		!(stfips==36 & cofips==61) & !(stfips==6  & cofips==37) & !(stfips==17 & cofips==31)
	) %>%
	# Create urban dummy
	group_by(fips) %>% 
	mutate(
		`_urb` = sum(`_60pcturban` * (year == 1960), na.rm = TRUE)
	) %>% 
	ungroup() %>%
	# cut every 25% with indicator from 0
	mutate(
		Durb = Hmisc::cut2(`_urb`, cuts = c(0, 1, 24.999, 49.999, 74.999, 110))
	) %>% 
	replace_na(list(chc_year_exp = 0))

ctr_pop <- read_sf("data/2010_county_centpop/") %>% 
	mutate(fips = as.numeric(paste0(STATEFP, COUNTYFP))) %>% 
	select(fips) %>% 
	arrange(fips)

df <- df %>% filter(fips != 12025) %>% left_join(., ctr_pop, by = "fips") %>% st_as_sf()


# Create within indicator

## helper function: for each year, calculate within
within_x_treat <- function(dist, year, chc_year_exp, geometry) {
	y <- min(year)
	treated <- (chc_year_exp <= y) & (chc_year_exp > 0)
	
	
	dist_mat <- st_distance(geometry, geometry) %>%
		units::set_units("mi") %>% units::drop_units()
	
	within <- (dist_mat > 0) & (dist_mat < dist)
	
	return(as.numeric( (within %*% treated) > 0))
} 

## helper function: find minimum year that within == 1, requires sorted!
min_year <- function(year, within) {
	if(any(within == 1)){
		idx <- which(within == 1)[1]
		y <- year[idx]
	}else{
		y <- 0
	}
	
	return(y)
}

df_spill <- df %>% 
	group_by(year) %>% 
	filter(year <= 1988) %>% 
	# Create within indicator
	mutate(
		treat = as.numeric((chc_year_exp < year) & (chc_year_exp != 0)),
		within_25 = within_x_treat(25, year, chc_year_exp, geometry),
		within_25 = within_25 * (1 - treat)
	) %>% 
	ungroup() 

df_spill <- df_spill %>% 
	# Create year when within_25 first turns on
	arrange(fips, year) %>% 
	group_by(fips) %>% 
	mutate(year_within_25 = min_year(year, within_25)) 


# Replicate Event Study --------------------------------------------------------

es_pts <- df %>% 
	filter(year <= 1988) %>%
	st_drop_geometry() %>% 
	fixest::feols(
		amr ~ i(exp1, ref = -1) + D_tot_act_md_t + D_60pctnonwhit_t + D_60pctrurf_t + D_60pcturban_t + D_pct59inclt3k_t + R_tranpcret + R_tranpcpa1 + H_bpc + H_hpc | fips + year + year^Durb + year^stfips, 
		data = .,
		weights = .$popwt
	) %>%
	broom::tidy() %>%
	filter(stringr::str_detect(term, "exp1::")) %>% 
	mutate(
		time = as.numeric(stringr::str_remove(term, "exp1::"))
	) %>% 
	select(time, estimate, std.error) %>%
	add_row(tibble(time = -1, estimate = 0, std.error = 0)) %>%
	mutate(
		ub = estimate + 1.96 * std.error,
		lb = estimate - 1.96 * std.error
	) %>%
	filter(time <= 14)

ggplot(es_pts) +
	geom_point(aes(x = time, y = estimate)) +
	geom_errorbar(aes(x = time, ymin = lb, ymax = ub), alpha = 0.8) + 
	geom_vline(xintercept = -0.5, color = "grey50") + 
	geom_hline(yintercept = 0, color = "black") +
	ylim(-70, 40) +
	theme_kyle() + 
	labs(y = NULL, x = "Years since CHC establishment")
	



# Modern Event Study -----------------------------------------------------------

df <- df %>% 
	mutate(
		chc_open = if_else(chc_year_exp <= 1974, chc_year_exp, 0)
	)

es_1 <- att_gt(
	yname = "amr",
	tname = "year",
	idname = "fips",
	weightsname = "popwt",
	est_method = "reg",
	xformla = ~ paste(year, Durb, sep="_") + paste(year, stfips, sep="_") + D_tot_act_md_t + D_60pctnonwhit_t + D_60pctrurf_t + D_60pcturban_t + D_pct59inclt3k_t + R_tranpcret + R_tranpcpa1 + H_bpc + H_hpc,
	gname = "chc_open",
	data = df %>% filter(year <= 1988) %>% drop_na(D_tot_act_md_t, D_60pctnonwhit_t, D_60pctrurf_t, D_60pcturban_t, D_pct59inclt3k_t, R_tranpcret, R_tranpcpa1, H_bpc, H_hpc)
)

agg_es <- aggte(es_1, type = "dynamic", na.rm = TRUE)

es_pts_did <- tibble(time = unlist(agg_es["egt"]), estimate = unlist(agg_es["att.egt"]), std.error = unlist(agg_es["se.egt"])) %>% 
	mutate(
		ub = estimate + 1.96 * std.error,
		lb = estimate - 1.96 * std.error
	) %>%
	filter(time >= -7 & time <= 14) 

ggplot(es_pts_did) +
	geom_point(aes(x = time, y = estimate)) +
	geom_errorbar(aes(x = time, ymin = lb, ymax = ub), alpha = 0.8) + 
	geom_vline(xintercept = -0.5, color = "grey50") + 
	geom_hline(yintercept = 0, color = "black") +
	ylim(-70, 40) +
	theme_kyle() + 
	labs(y = NULL, x = "Years since CHC establishment")



# Modern Event Study w/ Spillover Controls -------------------------------------

df_spill <- df_spill %>% 
	mutate(
		chc_open = if_else(chc_year_exp <= 1974, chc_year_exp, 0)
	)

es_spillover_control <- att_gt(
	yname = "amr",
	tname = "year",
	idname = "fips",
	weightsname = "popwt",
	est_method = "reg",
	xformla = ~ paste(year, Durb, sep="_") + paste(year, stfips, sep="_") + within_25 + D_tot_act_md_t + D_60pctnonwhit_t + D_60pctrurf_t + D_60pcturban_t + D_pct59inclt3k_t + R_tranpcret + R_tranpcpa1 + H_bpc + H_hpc,
	gname = "chc_year_exp",
	data = df_spill %>% filter(year <= 1988) %>% st_drop_geometry() %>% drop_na(D_tot_act_md_t, D_60pctnonwhit_t, D_60pctrurf_t, D_60pcturban_t, D_pct59inclt3k_t, R_tranpcret, R_tranpcpa1, H_bpc, H_hpc)
)


agg_es <- aggte(es_spillover_control, type = "dynamic", na.rm = TRUE)

es_pts_spillover_control <- tibble(time = unlist(agg_es["egt"]), estimate = unlist(agg_es["att.egt"]), std.error = unlist(agg_es["se.egt"])) %>% 
	mutate(
		ub = estimate + 1.96 * std.error,
		lb = estimate - 1.96 * std.error
	) %>%
	filter(time >= -7 & time <= 14) 

ggplot(es_pts_spillover_control) +
	geom_point(aes(x = time, y = estimate)) +
	geom_errorbar(aes(x = time, ymin = lb, ymax = ub), alpha = 0.8) + 
	geom_vline(xintercept = -0.5, color = "grey50") + 
	geom_hline(yintercept = 0, color = "black") +
	theme_kyle() + 
	labs(y = NULL, x = "Years since CHC establishment")


# Spillover Effects ------------------------------------------------------------


es_spillover_effect <- att_gt(
	yname = "amr",
	tname = "year",
	idname = "fips",
	weightsname = "popwt",
	est_method = "reg",
	xformla = ~ paste(year, Durb, sep="_") + paste(year, stfips, sep="_") + D_tot_act_md_t + D_60pctnonwhit_t + D_60pctrurf_t + D_60pcturban_t + D_pct59inclt3k_t + R_tranpcret + R_tranpcpa1 + H_bpc + H_hpc,
	gname = "year_within_25",
	data = df_spill %>% filter(year <= 1988) %>% st_drop_geometry() %>% filter(!(chc_year_exp <= 1974 & chc_year_exp > 0)) %>% drop_na(D_tot_act_md_t, D_60pctnonwhit_t, D_60pctrurf_t, D_60pcturban_t, D_pct59inclt3k_t, R_tranpcret, R_tranpcpa1, H_bpc, H_hpc)
)


agg_es <- aggte(es_spillover_effect, type = "dynamic", na.rm = TRUE)

es_pts_spillover_effect <- tibble(time = unlist(agg_es["egt"]), estimate = unlist(agg_es["att.egt"]), std.error = unlist(agg_es["se.egt"])) %>% 
	mutate(
		ub = estimate + 1.96 * std.error,
		lb = estimate - 1.96 * std.error
	) %>%
	filter(time >= -7 & time <= 14) 

ggplot(es_pts_spillover_effect) +
	geom_point(aes(x = time, y = estimate)) +
	geom_errorbar(aes(x = time, ymin = lb, ymax = ub), alpha = 0.8) + 
	geom_vline(xintercept = -0.5, color = "grey50") + 
	geom_hline(yintercept = 0, color = "black") +
	ylim(-75, 25) +
	theme_kyle() + 
	labs(y = NULL, x = "Years since CHC establishment")



# Save att_gt objects
# save(es_pts, es_1, es_spillover_control, es_spillover_effect, file = "data/chc-results.RData")



# Export figures ---------------------------------------------------------------

# treatment effect original
(es_plot_original <- ggplot(es_pts_did) +
 	geom_vline(xintercept = -0.5, color = "grey50") + 
 	geom_hline(yintercept = 0, color = "black") +
 	geom_point(aes(x = time, y = estimate), color = "#bf616a") +
 	geom_errorbar(aes(x = time, ymin = lb, ymax = ub), color = "#bf616a", alpha = 0.8) + 
 	theme_kyle(base_size = 16, slides = slides) + 
 	theme(title = element_text(size = 12, margin = margin(b = 0, unit = "pt"))) + 
 	scale_x_continuous(minor_breaks = seq(-7, 14, 1)) + 
 	scale_y_continuous(minor_breaks = seq(-75, 40, 5), limits = c(-75, 40)) + 
 	labs(y = "Deaths per 100,000 persons", x = "Years since CHC establishment", color = NULL) +
 	{ if(slides) labs(title = "Effect of CHC Establishment", subtitle = "(Bailey and Goodman-Bacon, 2015)") })


if(slides) ggsave("figures/figure-chc-es_original_slides.pdf", es_plot_original, width = 8, height = 8 * h_w_ratio, dpi = 300)
if(!slides) ggsave("figures/figure-chc-es_original.pdf", es_plot_original, width = 8, height = 5, dpi = 300)

if(slides) ggsave("figures/figure-chc-es_original_slides.png", es_plot_original, width = 8, height = 8 * h_w_ratio, dpi = 300)
if(!slides) ggsave("figures/figure-chc-es_original.png", es_plot_original, width = 8, height = 5, dpi = 300)


# treatment effect with controls
(es_plot_spillover_control <- ggplot(es_pts_spillover_control) +
 	geom_vline(xintercept = -0.5, color = "grey50") + 
 	geom_hline(yintercept = 0, color = "black") +
 	geom_point(aes(x = time, y = estimate), color = "#bf616a") +
 	geom_errorbar(aes(x = time, ymin = lb, ymax = ub), color = "#bf616a", alpha = 0.8) + 
 	theme_kyle(base_size = 16, slides = slides) + 
 	scale_x_continuous(minor_breaks = seq(-7, 14, 1)) + 
 	scale_y_continuous(minor_breaks = seq(-75, 25, 5), limits = c(-75, 25)) + 
 	labs(y = "Deaths per 100,000 persons", x = "Years since CHC establishment", color = NULL) +
 	{ if(slides) labs(title = "Effect of CHC Establishment", subtitle = "(Controlling for Spillovers)") })


if(slides) ggsave("figures/figure-chc-es_spillover_control_slides.pdf", es_plot_spillover_control, width = 8, height = 8 * h_w_ratio, dpi = 300)
if(!slides) ggsave("figures/figure-chc-es_spillover_control.pdf", es_plot_spillover_control, width = 8, height = 5, dpi = 300)

if(slides) ggsave("figures/figure-chc-es_spillover_control_slides.png", es_plot_spillover_control, width = 8, height = 8 * h_w_ratio, dpi = 300)
if(!slides) ggsave("figures/figure-chc-es_spillover_control.png", es_plot_spillover_control, width = 8, height = 5, dpi = 300)



# spillover effect
(es_plot_spillover_effect <- ggplot(es_pts_spillover_effect) +
		geom_vline(xintercept = -0.5, color = "grey50") + 
		geom_hline(yintercept = 0, color = "black") +
		geom_point(aes(x = time, y = estimate), color = "#5e81ac") +
		geom_errorbar(aes(x = time, ymin = lb, ymax = ub), color = "#5e81ac", alpha = 0.8) + 
		theme_kyle(base_size = 16, slides = slides) + 
		scale_x_continuous(minor_breaks = seq(-7, 14, 1)) + 
		scale_y_continuous(minor_breaks = seq(-75, 25, 5), limits = c(-75, 25)) + 
		labs(y = "Deaths per 100,000 persons", x = "Years since CHC establishment", color = NULL) +
		{ if(slides) labs(title = "Spillover onto Control Units within 25 miles") })


if(slides) ggsave("figures/figure-chc-es_spillover_effect_slides.pdf", es_plot_spillover_effect, width = 8, height = 8 * h_w_ratio, dpi = 300)
if(!slides) ggsave("figures/figure-chc-es_spillover_effect.pdf", es_plot_spillover_effect, width = 8, height = 5, dpi = 300)

if(slides) ggsave("figures/figure-chc-es_spillover_effect_slides.png", es_plot_spillover_effect, width = 8, height = 8 * h_w_ratio, dpi = 300)
if(!slides) ggsave("figures/figure-chc-es_spillover_effect.png", es_plot_spillover_effect, width = 8, height = 5, dpi = 300)



# combined 
es_pts_combined <- bind_rows(
	es_pts_spillover_control %>% mutate(group = "Treatment Effect", time = time - 0.1),
	es_pts_spillover_effect %>% mutate(group = "Spillover Onto Control", time = time + 0.1)
)

(es_plot_combined <- ggplot(es_pts_combined) +
	geom_vline(xintercept = -0.5, color = "grey50") + 
	geom_hline(yintercept = 0, color = "black") +
	geom_point(aes(x = time, y = estimate, color = group)) +
	geom_errorbar(aes(x = time, ymin = lb, ymax = ub, color = group), alpha = 0.8) + 
	theme_kyle(base_size = 16, slides = slides) + 
	theme(
		legend.position = c(0.15, 0.25),
		legend.spacing.x = unit(0, "pt"),
		legend.spacing.y = unit(0, "pt"),
		legend.background = element_rect(fill = if_else(slides, "#ECECEC", "white"))
	) +
	scale_shape_manual(values = c(16, 18)) + 
	scale_color_manual(values = c("#5e81ac", "#bf616a")) +
	scale_x_continuous(minor_breaks = seq(-7, 14, 1)) + 
	scale_y_continuous(minor_breaks = seq(-75, 25, 5), limits = c(-75, 25)) + 
	labs(y = "Deaths per 100,000 persons", x = "Years since CHC establishment", color = NULL) +
	{ if(slides) labs(title = "Effect of CHC Establishment") })


if(slides) ggsave("figures/figure-chc-es_combined_slides.pdf", es_plot_combined, width = 8, height = 8 * h_w_ratio, dpi = 300)
if(!slides) ggsave("figures/figure-chc-es_combined.pdf", es_plot_combined, width = 8, height = 5, dpi = 300)

if(slides) ggsave("figures/figure-chc-es_combined_slides.png", es_plot_combined, width = 8, height = 8 * h_w_ratio, dpi = 300)
if(!slides) ggsave("figures/figure-chc-es_combined.png", es_plot_combined, width = 8, height = 5, dpi = 300)


# Preview
# ggpreview <- function(..., device = "png", cairo = TRUE) {
# 	fname <- tempfile(fileext = paste0(".", device))
# 	
# 	if (cairo & device == "pdf") {
# 		ggplot2::ggsave(filename = fname, device = cairo_pdf, ...)
# 	} else if (cairo & device == "png") {
# 		ggplot2::ggsave(filename = fname, device = device, type = "cairo", ...)
# 	} else {
# 		ggplot2::ggsave(filename = fname, device = device, ...)
# 	}
# 	
# 	system2("open", fname)
# 	invisible(NULL)
# }
# 
# ggpreview(es_plot_combined, device = "pdf", cairo = TRUE, width = 8, height = 5, dpi = 300)














