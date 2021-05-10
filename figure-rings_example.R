## -----------------------------------------------------------------------------
## figure-binscatter_example.R
## Kyle Butts, CU Boulder Economics 
## 
## On Binscatter (2020)
## -----------------------------------------------------------------------------

library(tidyverse)
library(glue)
library(binsreg)

logname <- Sys.getenv("LOGNAME")
cli::cli_alert_success("User: {logname}")

if (logname == "kylebutts") setwd("~/Documents/Projects/Spatial Spillover/")
rm(logname)

# Load Custom Theme 
source("https://raw.githubusercontent.com/kylebutts/templates/master/ggplot_theme/theme_kyle.R")

# ggpreview
ggpreview <- function(..., device = "png", cairo = TRUE) {
	fname <- tempfile(fileext = paste0(".", device))
	
	if (cairo & device == "pdf") {
		ggplot2::ggsave(filename = fname, device = cairo_pdf, ...)
	} else if (cairo & device == "png") {
		ggplot2::ggsave(filename = fname, device = device, type = "cairo", ...)
	} else {
		ggplot2::ggsave(filename = fname, device = device, ...)
	}
	
	system2("open", fname)
	invisible(NULL)
}


# Generate Data ----------------------------------------------------------------

set.seed(100)

df <- tibble(id = 1:500) %>%
	# Generate Treatment and Distance
	mutate(
		treat = runif(n()) > 0.8,
		treat = as.numeric(treat),
		control = 1 - treat,
		dist = runif(n(), 0, 100),
		dist = dist * (1-treat)
	) %>% 
	# Make pre-post
	expand_grid(t = c(0,1)) %>%
	# Treatment Effects
	mutate(
		te = treat * t * 2,
		spill1 = 1.5 * t * exp(-0.03 * dist) * (1-treat) * (dist <= 40),
		eps = rnorm(n(), mean = 0, sd = 0.05),
		y1 = 2*t + te + spill1 + eps
	)


# Take First Differences
first_diff <- df %>% 
	arrange(id, t) %>%
	group_by(id) %>% 
	mutate(
		d_y1 = y1 - lag(y1)
	) %>% 
	filter(t == 1)

ggplot() + 
	geom_point(data=first_diff, aes(x=dist, y=d_y1), shape=16, col="lightgrey")

# Ad-Hoc Approach --------------------------------------------------------------

# Bins for 0-15, 15-30, 30-45
first_diff <- first_diff %>% 
	mutate(
		within_0_60  = dist > 0  & dist <= 60,
		within_0_15  = dist > 0  & dist <= 15,
		within_15_30 = dist > 15 & dist <= 30,
		within_30_45 = dist > 30 & dist <= 45,
		within_45_60 = dist > 45 & dist <= 60
	)

reg <- fixest::feols(d_y1 ~ treat + within_0_15 + within_15_30 + within_30_45 + within_45_60, data = first_diff) %>% coef()

dots <- tribble(
	~x, ~fit,
	0, reg[["treat"]] + reg[["(Intercept)"]],
	7.5, reg[["within_0_15TRUE"]] + reg[["(Intercept)"]],
	22.5, reg[["within_15_30TRUE"]] + reg[["(Intercept)"]],
	37.5, reg[["within_30_45TRUE"]] + reg[["(Intercept)"]],
	52.5, reg[["within_45_60TRUE"]] + reg[["(Intercept)"]],
)

line <- tribble(
	~x, ~fit,
	0, reg[["within_0_15TRUE"]] + reg[["(Intercept)"]],
	15, reg[["within_0_15TRUE"]] + reg[["(Intercept)"]], 
	15, NA,
	15, reg[["within_15_30TRUE"]] + reg[["(Intercept)"]],
	30, reg[["within_15_30TRUE"]] + reg[["(Intercept)"]],
	30, NA,
	30, reg[["within_30_45TRUE"]] + reg[["(Intercept)"]],
	45, reg[["within_30_45TRUE"]] + reg[["(Intercept)"]],
	45, NA,
	45, reg[["within_45_60TRUE"]] + reg[["(Intercept)"]],
	60, reg[["within_45_60TRUE"]] + reg[["(Intercept)"]]
)

reg_within <- fixest::feols(d_y1 ~ treat + within_0_60, data = first_diff) %>% coef()


dots_within <- tribble(
	~x, ~fit, 
	0, reg_within[["treat"]] + reg_within[["(Intercept)"]],
	30, reg_within[["within_0_60TRUE"]] + reg_within[["(Intercept)"]]
)

line_within <- tribble(
	~x, ~fit, 
	0, reg_within[["within_0_60TRUE"]] + reg_within[["(Intercept)"]],
	60, reg_within[["within_0_60TRUE"]] + reg_within[["(Intercept)"]], 
	60, NA,
)


(plot <- ggplot() + 
	geom_point(data=first_diff, aes(x=dist, y=d_y1), shape=16, col="lightgrey") + 
	geom_hline(yintercept = reg[[("(Intercept)")]], linetype = 'longdash') + 
	geom_point(data=dots, aes(x=x, y=fit), col="#5e81ac", size=2, shape=19) +
	geom_line(data=line, aes(x=x, y=fit), col="#5e81ac", size=0.5) + 
	theme_kyle(base_size = 16) + 
	labs(x = "Distance from Nearest Treated", y = "Change in Outcome") +
	scale_x_continuous(breaks=c(0, 15, 30, 45, 60, 100), labels=c("Treated", 15, 30, 45, 60, 100)) +
	# Annotation: Counterfactual Trend
	geom_curve(data = data.frame(x = 35, y = 1.85, xend = 42.5, yend = reg[[("(Intercept)")]] - 0.03),
			   mapping = aes(x = x, y = y, xend = xend, yend = yend),
			   arrow = arrow(30L, unit(0.1, "inches"), "last", "closed"),
			   inherit.aes = FALSE) + 
	geom_text(data = data.frame(x = -1, y = 1.82, label = "Counterfactual Trend Estimate"),
			  mapping = aes(x = x, y = y, label = label), size = 4, 
			  hjust = 0L, vjust = 0L, family = "fira_sans", fontface = 3, inherit.aes = FALSE))

(plot_within <- ggplot() + 
		geom_point(data=first_diff, aes(x=dist, y=d_y1), shape=16, col="lightgrey") + 
		geom_hline(yintercept = reg_within[[("(Intercept)")]], linetype = 'longdash') + 
		geom_point(data=dots_within, aes(x=x, y=fit), col="#bf616a", size=2, shape=19) +
		geom_line(data=line_within, aes(x=x, y=fit), col="#bf616a", size=0.5) + 
		theme_kyle(base_size = 16) + 
		labs(x = "Distance from Nearest Treated", y = "Change in Outcome") +
		scale_x_continuous(breaks=c(0, 15, 30, 45, 60, 100), labels=c("Treated", 15, 30, 45, 60, 100)) +
		# Annotation: Counterfactual Trend
		geom_curve(data = data.frame(x = 35, y = 1.85, xend = 42.5, yend = reg_within[[("(Intercept)")]] - 0.03),
				   mapping = aes(x = x, y = y, xend = xend, yend = yend),
				   arrow = arrow(30L, unit(0.1, "inches"), "last", "closed"),
				   inherit.aes = FALSE) + 
		geom_text(data = data.frame(x = -1, y = 1.82, label = "Counterfactual Trend Estimate"),
				  mapping = aes(x = x, y = y, label = label), size = 4, 
				  hjust = 0L, vjust = 0L, family = "fira_sans", fontface = 3, inherit.aes = FALSE))

# ggpreview(plot, width = 2400/300, height = 2400/300 * 9/16, dpi = 300, device="pdf")
# ggpreview(plot_within, width = 2400/300, height = 2400/300 * 9/16, dpi = 300, device="pdf")

ggsave("figures/figure-rings_example.pdf", plot, width = 2400/300, height = 2400/300 * 9/16, dpi = 300)
ggsave("figures/figure-within_example.pdf", plot_within, width = 2400/300, height = 2400/300 * 9/16, dpi = 300)


ggsave("figures/figure-rings_example.png", plot, width = 2400/300, height = 2400/300 * 9/16, dpi = 300)
ggsave("figures/figure-within_example.png", plot_within, width = 2400/300, height = 2400/300 * 9/16, dpi = 300)










