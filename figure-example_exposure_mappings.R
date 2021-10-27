## -----------------------------------------------------------------------------
## figure-example_exposure_mappings.R
## Kyle Butts, CU Boulder Economics 
## 
## Create examples of each exposure mapping using US counties
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

# Export Options
slides <- TRUE
h_w_ratio <- 9/16

# sim_data_misspecification()
source("helper-sim_function_misspecification.R")


## Load Spatial Data -----------------------------------------------------------
# From data_prepare_counties.R
load(file= here::here("data/counties_and_mat.RData"))

us <- counties %>% summarize()

set.seed(20210215)

df <- sim_data_misspecification(drop_geometry = FALSE, normalize = TRUE) %>% 
	st_as_sf()


## Treatment Effect ------------------------------------------------------------

types <- c(
	"Within 40mi." = "spill_within", 
	"Within 80mi." = "spill_within_large", 
	"Within 40mi. (Additive)" = "spill_within_additive", 
	"Within 80mi. (Additive)" = "spill_within_large_additive", 
	"Decay" = "spill_decay", 
	"Decay (Additive)" = "spill_decay_additive"
)



for(type in types){
	title <- names(types)[types == type]
	cli::cli_rule(glue("{title}: {type}"))

	temp <- df %>% 
		filter(year == 2019) %>% 
		mutate(!!rlang::sym(glue("te_{type}")) := if_else(treat == 1, NA_real_, !!rlang::sym(glue("te_{type}"))))
	

	
	(
		plot <- ggplot(data = temp) + 
			geom_sf(mapping = aes(fill = !!rlang::sym(glue("te_{type}")), color = as.factor(treat)), size= 0.2) + 
			geom_sf(data = us, color = "black", fill = NA) + 
			scale_color_manual(values = c("white", "black"), labels = c("Not Treated", "Treated")) +
			coord_sf(datum = NA) +
			labs(
				title = title,
				fill = "Spillover", color = "Treated"
			) + 
			scale_fill_gradient(low = "white", high = "grey10", na.value = NA, limits = c(0, 2)) +
			theme_kyle() + 
			guides(
				fill = guide_legend(title.position = "top"),
				color = guide_legend(title.position = "top")
			) +
			theme(title = element_text(lineheight = 1.2))
	)

	ggsave(glue("figures/figure-{type}.png"), plot, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "transparent")	
}







## Rings -----------------------------------------------------------------------

df <- df %>% 
	mutate(
		spill = case_when(
			spill_0_20 == 1 ~ "0 to 20 miles",
			spill_20_30 == 1 ~ "20 to 30 miles",
			spill_30_40 == 1 ~ "30 to 40 miles",
			spill_40_60 == 1 ~ "40 to 60 miles",
			spill_60_80 == 1 ~ "60 to 80 miles",
		)
	)

df$spill <- factor(df$spill, levels = c("0 to 20 miles", "20 to 30 miles", "30 to 40 miles", "40 to 60 miles", "60 to 80 miles"))

grey_palette <- glue("grey{seq(10, 80, by=14)}")

(plot <- ggplot(data = df) + 
		geom_sf(mapping = aes(fill = spill, color = as.factor(treat)), size= 0.2) + 
		geom_sf(data = us, color = "black", fill = NA) + 
		scale_color_manual(values = c("white", "black"), labels = c("Not Treated", "Treated")) +
		coord_sf(datum = NA) +
		labs(
			title = "Rings (0-20, 20-30, 30-40, 40-60, 60-80)",
			fill = "Exposure", color = "Treated"
		) + 
		scale_fill_manual(values = grey_palette, na.value="white") +
		# scale_fill_gradient(low = "white", high = "#E34A33", na.value = NA) +
		theme_kyle() + 
		guides(
			fill = guide_legend(title.position = "top"),
			color = guide_legend(title.position = "top")
		) +
		theme(title = element_text(lineheight = 1.2))
)



ggsave(glue("figures/figure-spill_ring.png"), plot, dpi= 300, width= 2400/300, height= 2400/300 * h_w_ratio, bg= "transparent")	












