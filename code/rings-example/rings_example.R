#' ---
#' title: "Example single indicator vs. multiple rings"
#' ---

# %%
#| message: false
#| warning: false
library(tidyverse, warn.conflicts = FALSE)
library(glue)
library(here)
library(binsreg)
library(kfbmisc)

# Counterfactual Trends Everywhere ---------------------------------------------
#' First we will generate a setting where parallel trends is constant across distance from nearest treated (e.g. small neighborhood around treatment where treatment location is random). This makes rings very effective. 
 
# %%
set.seed(100)
df <- tibble(id = 1:500) |>
  # Generate Treatment and Distance
  mutate(
    treat = runif(n()) > 0.8,
    treat = as.numeric(treat),
    control = 1 - treat,
    dist = runif(n(), 0, 90),
    dist = dist * (1 - treat)
  ) |>
  # Make pre-post
  expand_grid(t = c(0, 1)) |>
  # Treatment Effects
  mutate(
    te = treat * t * 2,
    spill1 = 1.5 * t * exp(-0.03 * dist) * (1 - treat) * (dist <= 40),
    eps = rnorm(n(), mean = 0, sd = 0.05),
    y1 = 2 * t + te + spill1 + eps
  )

# %%
first_diff <- df |>
  arrange(id, t) |>
  group_by(id) |>
  mutate(
    d_y1 = y1 - lag(y1)
  ) |>
  filter(t == 1)

ggplot() +
  geom_point(
    data = first_diff, aes(x = dist, y = d_y1), 
    shape = 16, col = "lightgrey"
  ) +
  labs(
    x = "Distance to nearest treated", 
    y = "Change in y"
  ) +
  kfbmisc::theme_kyle(base_size = 12)

## Ad-Hoc Approach -------------------------------------------------------------

# %%
# Bins for 0-15, 15-30, 30-45
first_diff <- first_diff |>
  mutate(
    within_0_60  = dist > 0 & dist <= 60,
    within_0_15  = dist > 0 & dist <= 15,
    within_15_30 = dist > 15 & dist <= 30,
    within_30_45 = dist > 30 & dist <= 45,
    within_45_60 = dist > 45 & dist <= 60
  )

reg <- fixest::feols(
  d_y1 ~ treat + within_0_15 + within_15_30 + within_30_45 + within_45_60,
  data = first_diff
) |> coef()

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

reg_within <- fixest::feols(
  d_y1 ~ treat + within_0_60,
  data = first_diff
) |> coef()

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

# %%
(plot_rings <- ggplot() +
  geom_point(data = first_diff, aes(x = dist, y = d_y1), shape = 16, col = "lightgrey") +
  geom_hline(yintercept = reg[[("(Intercept)")]], linetype = "longdash") +
  geom_point(data = dots, aes(x = x, y = fit), col = "#5e81ac", size = 3, shape = 19) +
  geom_line(data = line, aes(x = x, y = fit), col = "#5e81ac", linewidth = 1.5) +
  theme_kyle(base_size = 12) +
  labs(
    x = "Distance from Nearest Treated", 
    y = "Change in Outcome ($\\Delta Y$)"
  ) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90), labels = c("Treated", 15, 30, 45, 60, 75, 90)) +
  scale_y_continuous(limits = c(1, 4.25)) +
  # Annotation: Counterfactual Trend
  geom_text(
    data = data.frame(x = 1, y = 1.78, label = "Counterfactual Trend Estimate"),
    mapping = aes(x = x, y = y, label = label), size = 4,
    hjust = 0L, vjust = 0L, fontface = 3, inherit.aes = FALSE
  ))

# %%
(plot_within <- ggplot() +
  geom_point(data = first_diff, aes(x = dist, y = d_y1), shape = 16, col = "lightgrey") +
  geom_hline(yintercept = reg_within[[("(Intercept)")]], linetype = "longdash") +
  geom_point(data = dots_within, aes(x = x, y = fit), col = "#bf616a", size = 3, shape = 19) +
  geom_line(data = line_within, aes(x = x, y = fit), col = "#bf616a", linewidth = 1.5) +
  theme_kyle(base_size = 12) +
  labs(
    x = "Distance from Nearest Treated", 
    y = "Change in Outcome ($\\Delta Y$)"
  ) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90), labels = c("Treated", 15, 30, 45, 60, 75, 90)) +
  scale_y_continuous(limits = c(1, 4.25)) +
  # Annotation: Counterfactual Trend
  geom_text(
    data = data.frame(x = 1, y = 1.78, label = "Counterfactual Trend Estimate"),
    mapping = aes(x = x, y = y, label = label), size = 4,
    hjust = 0L, vjust = 0L, fontface = 3, inherit.aes = FALSE
  ))

# %%
kfbmisc::tikzsave(
  here("figures/rings-example/rings_pt_everywhere.pdf"),
  plot_rings, width = 8, height = 8 * 9 / 16
)
kfbmisc::tikzsave(
  here("figures/rings-example/within_pt_everywhere.pdf"),
  plot_within, width = 8, height = 8 * 9 / 16
)

# %%
(plot_both <- ggplot() +
  geom_point(
    data = first_diff |> filter(dist <= 90),
    aes(x = dist, y = d_y1), shape = 16, col = "lightgrey"
  ) +
  geom_hline(yintercept = reg[[("(Intercept)")]], linetype = "dotted") +
  # Within
  geom_point(data = dots_within, aes(x = x, y = fit, color = "#bf616a"), size = 3, shape = 19) +
  geom_line(data = line_within, aes(x = x, y = fit, color = "#bf616a", linetype = 1), linewidth = 1.5) +
  # Rings
  geom_point(data = dots, aes(x = x, y = fit, color = "#5e81ac"), size = 3, shape = 19) +
  geom_line(data = line, aes(x = x, y = fit, color = "#5e81ac", linetype = 1), linewidth = 1.5) +
  theme_kyle(base_size = 12) +
  labs(
    x = "Distance from Nearest Treated", y = "Change in Outcome ($\\Delta Y$)",
    color = "Estimation Method"
  ) +
  scale_x_continuous(
    breaks = c(0, 15, 30, 45, 60, 75, 90),
    labels = c("Treated", 15, 30, 45, 60, 75, 90)
  ) +
  scale_y_continuous(limits = c(1, 4.25)) +
  scale_color_identity(
    breaks = c("#bf616a", "#5e81ac"),
    labels = c("One Ring", "Multiple Rings"),
    guide = guide_legend(
      override.aes = list(linetype = c(1, 1), shape = c(NA, NA), size = c(2, 2))
    )
  ) +
  scale_linetype_identity(guide = "none") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white")
    # legend.key.width = unit(4,"cm")
  ) +
  # Annotation: Counterfactual Trend
  geom_text(
    data = data.frame(x = 1, y = 1.78, label = "Counterfactual Trend Estimate"),
    mapping = aes(x = x, y = y, label = label), size = 4,
    hjust = 0L, vjust = 0L, fontface = 3, inherit.aes = FALSE
  ))

# %%
kfbmisc::tikzsave(
  here::here("figures/rings-example/rings_v_within_pt_everywhere.pdf"), 
  plot_both, width = 8, height = 8 * 9 / 16
)



# Average Counterfactual Trends Only -------------------------------------------
#' No we will generate a setting where parallel trends holds on average for the treated and the control units (for a given $\bar{d}$), but does not hold at every distance. In this setting, rings can not be trusted since they are based on the outer-most rings' counterfactual trend.

# %% 
set.seed(100)

df <- tibble(id = 1:500) |>
  # Generate Treatment and Distance
  mutate(
    treat = runif(n()) > 0.8,
    treat = as.numeric(treat),
    control = 1 - treat,
    dist = runif(n(), 0, 90),
    dist = dist * (1 - treat)
  ) |>
  # Make pre-post
  expand_grid(t = c(0, 1)) |>
  # Treatment Effects
  mutate(
    te = treat * t * 2,
    trend = 0.6 * sin(pi / 20 * dist),
    spill1 = 1.5 * t * exp(-0.03 * dist) * (1 - treat) * (dist <= 40),
    eps = rnorm(n(), mean = 0, sd = 0.05),
    y0 = trend * t + eps,
    y1 = trend * t + 2 * t + te + spill1 + eps
  )

# Take First Differences
first_diff <- df |>
  arrange(id, t) |>
  group_by(id) |>
  mutate(
    d_y1 = y1 - lag(y1),
    d_y0 = y0 - lag(y0)
  ) |>
  filter(t == 1)

ggplot() +
  geom_point(data = first_diff, aes(x = dist, y = d_y1), shape = 16, col = "lightgrey")

## Ad-Hoc Approach -------------------------------------------------------------
# %%
# Bins for 0-15, 15-30, 30-45
first_diff <- first_diff |>
  mutate(
    within_0_60  = dist > 0 & dist <= 60,
    within_0_15  = dist > 0 & dist <= 15,
    within_15_30 = dist > 15 & dist <= 30,
    within_30_45 = dist > 30 & dist <= 45,
    within_45_60 = dist > 45 & dist <= 60
  )

reg <- fixest::feols(
  d_y1 ~ treat + within_0_15 + within_15_30 + within_30_45 + within_45_60, 
  data = first_diff
) |> coef()

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

reg_within <- fixest::feols(d_y1 ~ treat + within_0_60, data = first_diff) |> coef()

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

# %%
(plot_rings <- ggplot() +
  geom_point(data = first_diff, aes(x = dist, y = d_y1), shape = 16, col = "lightgrey") +
  geom_hline(yintercept = reg[[("(Intercept)")]], linetype = "longdash") +
  geom_point(data = dots, aes(x = x, y = fit), col = "#5e81ac", size = 3, shape = 19) +
  geom_line(data = line, aes(x = x, y = fit), col = "#5e81ac", linewidth = 1.5) +
  theme_kyle(base_size = 12) +
  labs(x = "Distance from Nearest Treated", y = "Change in Outcome ($\\Delta Y$)") +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90), labels = c("Treated", 15, 30, 45, 60, 75, 90)) +
  scale_y_continuous(limits = c(1, 4.25)) +
  # Annotation: Counterfactual Trend
  geom_text(
    data = data.frame(x = 1, y = 1.78, label = "Counterfactual Trend Estimate"),
    mapping = aes(x = x, y = y, label = label), size = 4,
    hjust = 0L, vjust = 0L, fontface = 3, inherit.aes = FALSE
  ))

# %% 
(plot_within <- ggplot() +
  geom_point(data = first_diff, aes(x = dist, y = d_y1), shape = 16, col = "lightgrey") +
  geom_hline(yintercept = reg_within[[("(Intercept)")]], linetype = "longdash") +
  geom_point(data = dots_within, aes(x = x, y = fit), col = "#bf616a", size = 3, shape = 19) +
  geom_line(data = line_within, aes(x = x, y = fit), col = "#bf616a", linewidth = 1.5) +
  theme_kyle(base_size = 12) +
  labs(x = "Distance from Nearest Treated", y = "Change in Outcome ($\\Delta Y$)") +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90), labels = c("Treated", 15, 30, 45, 60, 75, 90)) +
  scale_y_continuous(limits = c(1, 4.25)) +
  # Annotation: Counterfactual Trend
  geom_text(
    data = data.frame(x = 1, y = 1.78, label = "Counterfactual Trend Estimate"),
    mapping = aes(x = x, y = y, label = label), size = 4,
    hjust = 0L, vjust = 0L, fontface = 3, inherit.aes = FALSE
  ))

# %% 
kfbmisc::tikzsave(
  here::here("figures/rings-example/rings_pt_on_avg.pdf"), 
  plot_rings, width = 8, height = 8 * 9 / 16
)
kfbmisc::tikzsave(
  here::here("figures/rings-example/within_pt_on_avg.pdf"), 
  plot_within, width = 8, height = 8 * 9 / 16
)

# %% 
(plot_both <- ggplot() +
  geom_point(
    data = first_diff |> filter(dist <= 90),
    aes(x = dist, y = d_y1), shape = 16, col = "lightgrey"
  ) +
  geom_hline(yintercept = reg[[("(Intercept)")]], linetype = "dotted") +
  # Within
  geom_point(data = dots_within, aes(x = x, y = fit, color = "#bf616a"), size = 3, shape = 19) +
  geom_line(data = line_within, aes(x = x, y = fit, color = "#bf616a", linetype = 1), linewidth = 1.5) +
  # Rings
  geom_point(data = dots, aes(x = x, y = fit, color = "#5e81ac"), size = 3, shape = 19) +
  geom_line(data = line, aes(x = x, y = fit, color = "#5e81ac", linetype = 1), linewidth = 1.5) +
  theme_kyle(base_size = 12) +
  labs(
    x = "Distance from Nearest Treated", y = "Change in Outcome ($\\Delta Y$)",
    color = "Estimation Method"
  ) +
  scale_x_continuous(
    breaks = c(0, 15, 30, 45, 60, 75, 90),
    labels = c("Treated", 15, 30, 45, 60, 75, 90)
  ) +
  scale_y_continuous(limits = c(1, 4.25)) +
  scale_color_identity(
    breaks = c("#bf616a", "#5e81ac"),
    labels = c("One Ring", "Multiple Rings"),
    guide = guide_legend(
      override.aes = list(linetype = c(1, 1), shape = c(NA, NA), size = c(2, 2))
    )
  ) +
  scale_linetype_identity(guide = "none") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = "white")
    # legend.key.width = unit(4,"cm")
  ) +
  # Annotation: Counterfactual Trend
  geom_curve(
    data = data.frame(x = 35, y = 1.85, xend = 42.5, yend = reg[[("(Intercept)")]] - 0.03),
    mapping = aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(30L, unit(0.1, "inches"), "last", "closed"),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = data.frame(x = 1, y = 1.78, label = "Counterfactual Trend Estimate"),
    mapping = aes(x = x, y = y, label = label), size = 4,
    hjust = 0L, vjust = 0L, fontface = 3, inherit.aes = FALSE
  ))

# %% 
kfbmisc::tikzsave(
  here::here("figures/rings-example/rings_v_within_pt_on_avg.pdf"), plot_both,
  width = 8, height = 8 * 9 / 16, dpi = 300
)
