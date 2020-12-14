## -----------------------------------------------------------------------------
## figure-gen_map_examples_with_bias.R
## Kyle Butts, CU Boulder Economics 
## 
## Generate maps of Spatial Spillover, Spatial Autocorrelation of Treatment, etc. and their bias
## -----------------------------------------------------------------------------

setwd("~/Documents/Projects/Spatial Spillover/")

library("tidyverse")
library("sf")
library("stars")
library("units")
library("tigris")
library("gstat")
library("exactextractr") 
library("fixest")
library("doParallel")  
library("RColorBrewer")
library("patchwork")

# sim_data(i)
source("helper-sim_function.R")

# Load theme_kyle()
source("https://raw.githubusercontent.com/kylebutts/templates/master/ggplot_theme/theme_kyle.R")

# Is for slides or paper?
slides <- TRUE
# h_w_ratio <- 9/16
h_w_ratio <- 3/4


## Load Parameters -------------------------------------------------------------
te <- 2
te_spill <- 1
te_spill_treat <- -0.5
low_zone_plus <- 0.4
high_zone_plus <- 1.4


##------------------------------------------------------------------------------
##  Prepare Spatial Data                                                    ----
##------------------------------------------------------------------------------

## Load Spatial Data -----------------------------------------------------------
# From data_prepare_counties.R
load(file= "data/counties_and_mat.RData")

us <- counties %>% st_make_valid() %>% st_buffer(0) %>% st_set_precision(1e4) %>% summarize() 


##------------------------------------------------------------------------------
##  Generate Maps with No Spatial Autocorrelation                           ----
##------------------------------------------------------------------------------

set.seed(70)

treat_prob <- 0.05
treat_prob_percent <- paste0(round(treat_prob * 100, 1), "%")

counties_treat <- sim_data(
		treat_prob = treat_prob, treat_spillover_distance = 30, treat_effect = te, 
		treat_effect_spill = te_spill, treat_effect_spill_treat = te_spill_treat, spill_type = "contig", 
		spatial_autocorr = FALSE, zone_plus = 0, drop_geometry = FALSE
	) %>% 
	filter(year == 2019) %>% select(-centroid) %>% st_as_sf()

bias_control <- round(te_spill * sum(counties_treat$spill_within) / sum(counties_treat$treat == 0), 2)
bias_treat <- round(te_spill_treat * sum(counties_treat$spill_within_treat) / sum(counties_treat$treat == 1), 2)
estimate <- te
estimate_control <- te - bias_control
estimate_both <- te + bias_treat - bias_control

# Treatment Effect Maps
(plot_te <- ggplot() +
	geom_sf(data= counties_treat, aes(fill= factor(te)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat %>% filter(treat == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title= "Direct Effect of Treatment",
		subtitle= glue::glue("Treatment Effect Estimate: {estimate} \n Treated units outlined in black"),
		fill= "Effect Size"
	) + 
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	# theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) + 
	theme(title = element_text(lineheight = 1.2)))


plot_te_spill_control <- ggplot(counties_treat) +
	geom_sf(data= counties_treat, aes(fill= factor(te_spill + te)), color = "grey60", size = 0.1) + 
	geom_sf(data= counties_treat %>% filter(treat == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title= "Direct Effect + Spillover on Control", 
		subtitle= glue::glue("Treatment Effect Estimate: {estimate_control} \n Treated units outlined in black"),
		fill= "Effect Size"
	) + 
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	#theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) + 
	theme(title = element_text(lineheight = 1.2))

plot_te_spill_all <- ggplot(counties_treat) +
	geom_sf(data= counties_treat, aes(fill= factor(te + te_spill + te_spill_treat)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat %>% filter(treat == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title= "Direct Effect + Spillover on Control and Treated", 
		subtitle= glue::glue("Treatment Effect Estimate: {estimate_both} \n Treated units outlined in black"),
		fill= "Effect Size"
	) +
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	#theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) + 
	theme(title = element_text(lineheight = 1.2))


# Export Maps for Beamer Slides
ggsave("figures/figure-map_te.png", plot_te, 
	   dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
ggsave("figures/figure-map_te_spill_control.png", plot_te_spill_control, 
	   dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
ggsave("figures/figure-map_te_spill_all.png", plot_te_spill_all, 
	   dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")

cli::cli_alert_success("Finished exporting no sp_corr maps")




##------------------------------------------------------------------------------
##  Generate Maps with Spatial Autocorrelation                              ----
##------------------------------------------------------------------------------

## Low Zone Plus ---------------------------------------------------------------

# Generate Treatment Variables
counties_treat_low_zp <- sim_data(
	treat_prob = treat_prob, treat_spillover_distance = 40, treat_effect = 2, 
	treat_effect_spill = 1, treat_effect_spill_treat = -0.5, spill_type = "contig", 
	spatial_autocorr = TRUE, zone_plus = low_zone_plus, drop_geometry = FALSE
) %>% 
	filter(year == 2019) %>% select(-centroid) %>% st_as_sf()
	

bias_control <- round(te_spill * sum(counties_treat_low_zp$spill_within) / sum(counties_treat_low_zp$treat_ind == 0), 2)
bias_treat <- round(te_spill_treat * sum(counties_treat_low_zp$spill_within_treat) / sum(counties_treat_low_zp$treat_ind == 1), 2)


# Treatment Effect Maps, spcorr_low_zp
(
plot_te_low_zp <- ggplot() +
	geom_sf(data= counties_treat_low_zp, aes(fill= as.factor(te)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_low_zp %>% filter(treat_ind == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title= "Small Spatial Autocorrelation",
		subtitle= glue::glue("Treatment Effect Estimate: {estimate} \n Spatial Autocorrelation Measure = {low_zone_plus}"),
		fill= "Effect Size"
	) + 
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	# theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) 
)


(
plot_te_spill_control_low_zp <- ggplot() +
	geom_sf(data= counties_treat_low_zp, aes(fill= as.factor(te_spill + te)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_low_zp %>% filter(treat_ind == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title= "Small Spatial Autocorrelation", 
		subtitle= glue::glue("Treatment Effect Estimate: {estimate_control} \n Spatial Autocorrelation Measure = {low_zone_plus}"),
		fill= "Effect Size"
	) + 
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	#theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A"))
)

(
plot_te_spill_all_low_zp <- ggplot() +
	geom_sf(data= counties_treat_low_zp, aes(fill= as.factor(te + te_spill + te_spill_treat)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_low_zp %>% filter(treat == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title= "Small Spatial Autocorrelation", 
		subtitle= glue::glue("Treatment Effect Estimate: {estimate_both} \n Spatial Autocorrelation Measure = {low_zone_plus}"),
		fill= "Effect Size"
	) +
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	# theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) 
)


ggsave("figures/figure-spcorr_low_map_te.png", plot_te_low_zp, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
ggsave("figures/figure-spcorr_low_map_te_spill_control.png", plot_te_spill_control_low_zp, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
ggsave("figures/figure-spcorr_low_map_te_spill_all.png", plot_te_spill_all_low_zp, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
usethis::ui_done("Finished exporting low_zp maps")


## High Zone Plus --------------------------------------------------------------

# Generate Treatment Variables
counties_treat_high_zp <- sim_data(
	treat_prob = treat_prob, treat_spillover_distance = 40, treat_effect = 2, 
	treat_effect_spill = 1, treat_effect_spill_treat = -0.5, spill_type = "contig", 
	spatial_autocorr = TRUE, zone_plus = high_zone_plus, drop_geometry = FALSE
) %>% 
	filter(year == 2019) %>% select(-centroid) %>% st_as_sf()

bias_control <- round(te_spill * sum(counties_treat_high_zp$spill_within) / sum(counties_treat_high_zp$treat == 0), 2)
bias_treat <- round(te_spill_treat * sum(counties_treat_high_zp$spill_within_treat) / sum(counties_treat_high_zp$treat == 1), 2)

# Treatment Effect Maps, spcorr_high_zp
(
plot_te_high_zp <- ggplot() +
	geom_sf(data= counties_treat_high_zp, aes(fill= as.factor(te)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_high_zp %>% filter(treat_ind == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title= "Large Spatial Autocorrelation",
		subtitle= glue::glue("Treatment Effect Estimate: {estimate} \n Spatial Autocorrelation Measure = {high_zone_plus}"),
		fill= "Effect Size"
	) + 
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	# theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) 
)


(
plot_te_spill_control_high_zp <- ggplot() +
	geom_sf(data= counties_treat_high_zp, aes(fill= as.factor(te_spill + te)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_high_zp %>% filter(treat_ind == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title= "Large Spatial Autocorrelation", 
		subtitle= glue::glue("Treatment Effect Estimate: {estimate_both} \n Spatial Autocorrelation Measure = {high_zone_plus}"),
		fill= "Effect Size"
	) + 
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	# theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) 
)

(
plot_te_spill_all_high_zp <- ggplot() +
	geom_sf(data= counties_treat_high_zp, aes(fill= as.factor(te + te_spill + te_spill_treat)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_high_zp %>% filter(treat == 1), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		# title= "Direct Effect + Spillover on Control and Treated", 
		title = "Large Spatial Autocorrelation",
		subtitle= glue::glue("Treatment Effect Estimate: {estimate_both} \n Spatial Autocorrelation Measure = {high_zone_plus}"),
		fill= "Effect Size"
	) +
	theme_kyle(slides = TRUE) + 
	# Put Legend on Bottom
	guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	# theme(legend.position = "bottom") +
	# Fill Scale, colors from R color Brewer
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) 
)


ggsave("figures/figure-spcorr_high_map_te.png", plot_te_high_zp, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
ggsave("figures/figure-spcorr_high_map_te_spill_control.png", plot_te_spill_control_high_zp, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
ggsave("figures/figure-spcorr_high_map_te_spill_all.png", plot_te_spill_all_high_zp, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
usethis::ui_done("Finished exporting high_zp maps")





## Kriging Example -------------------------------------------------------------
plot_krig <- ggplot() +
	geom_sf(data= counties_treat_high_zp, aes(fill= sim), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_high_zp %>% filter(zone == TRUE), fill = NA, color = "black") +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	theme_kyle(base_size = 18, slides = slides) + {
		if(slides) labs(
			title = "Kriging Example"
		)
	} +
	labs(fill= NULL) +
	# Put Legend on Bottom
	# guides(fill = guide_legend(title.position = "top", nrow = 4)) +
	# Fill Scale, colors from R color Brewer
	scale_fill_distiller(palette = "Spectral")

plot_krig

plot_krig_highzp <- ggplot() +
	geom_sf(data= counties_treat_high_zp, aes(fill= as.factor(treat)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_high_zp %>% filter(zone == TRUE), fill = NA, color = "black", size = 0.4) +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title = "Zone Plus = 1.4",
		fill= "Treated Unit"
	) +
	theme_kyle(slides = slides) +
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) +
	guides(fill = guide_legend(title.position = "top", nrow = 4))

counties_treat_high_zp <- counties_treat_high_zp %>% mutate(
	# `zone` has increased probability of treatment
	prob_low = 0.1 + low_zone_plus * zone,
	prob_low = prob_low * treat_prob/mean(prob_low),
	
	# Generate Treatment Variable
	treat_low = (runif(n()) <= prob_low),
	treat_low = as.numeric(treat_low)
) 

plot_krig_lowzp <- ggplot() +
	geom_sf(data= counties_treat_high_zp, aes(fill= as.factor(treat_low)), color= "grey60", size= 0.2) + 
	geom_sf(data= counties_treat_high_zp %>% filter(zone == TRUE), fill = NA, color = "black", size = 0.4) +
	geom_sf(data= us, fill = NA, color= "grey30", size= 0.2) + 
	# Remove Coordinates, leaving just the map
	coord_sf(datum = NA) +
	labs(
		title = "Zone Plus = 0.4",
		fill= "Treated Unit"
	) +
	theme_kyle(slides = slides) +
	scale_fill_manual(values= c("0" = "#ffffff", "1" = "#FBB4B9", "1.5" = "#F768A1", "2" = "#C51B8A")) + 
	guides(fill = guide_legend(title.position = "top", nrow = 4))



if(slides) ggsave("figures/figure-krig_slides.png", plot_krig, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "#ECECEC")
	
if(!slides) ggsave("figures/figure-krig.png", plot_krig, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "transparent")

if(!slides) ggsave("figures/figure-krig_lowzp.png", plot_krig_lowzp, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "transparent")

if(!slides) ggsave("figures/figure-krig_highzp.png", plot_krig_highzp, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "transparent")
