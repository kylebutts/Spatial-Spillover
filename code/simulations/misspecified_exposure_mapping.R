#' ---
#' title: "Misspecification of Exposure Mapping"
#' ---
#' This file performs a set of Monte Carlo simulations that estimate a set of
#' DGPs using various exposure mappings and estimates them using exposure
#' mappings and ring(s). The first simulation reports on average bias found in
#' the treatment effect estimate and the second simulation reports on the
#' ability to accurately predict spillover effects from fitted values.

# %%
#| message: false
#| warning: false
library(tidyverse, warn.conflicts = FALSE)
library(data.table)
library(here)
library(glue)
library(sf)
library(stars)
library(fixest)
library(doParallel)
library(gt)
library(kfbmisc)
library(Matrix)

# flags
EXPORT <- FALSE
RUN <- TRUE

# %%
# extract_body()
extract_body <- function(gt) {
  if (inherits(gt, "gt_tbl")) gt <- as.character(gt::as_latex(gt))
	gt = stringr::str_split_1(gt, "\\n")
	
	idx = c(
		which(str_detect(gt, "\\\\midrule")) + 1,
		which(str_detect(gt, "\\\\bottomrule")) - 1
	)
	gt = paste0(gt[idx[1]:idx[2]], collapse = "\n")
  stringr::str_remove(gt, " \\\\\\\\ \n$")
}

# sim_data_misspecification()
source(here::here("code/simulations/helper-sim_function_misspecification.R"))

# Helper function estimate
estimate_treatment_effect <- function(df, formula, treat_var) {
  coef(feols(formula, data = df, warn = FALSE, notes = FALSE))[[treat_var]]
}


#' From data_prepare_counties.R
# %%
# counties, dist_mi
load(file = here::here("data/counties_and_mat.RData"))

counties = counties |> st_drop_geometry() |> as.data.table()

# Logical matrix, =1 if within the maxmium distance for spillover
within <- (dist_mi > 0 & dist_mi <= 40)
within_large <- (dist_mi > 0 & dist_mi <= 80)
within_100 <- (dist_mi > 0 & dist_mi <= 100)
within_0_20 <- (dist_mi > 0 & dist_mi <= 20)
within_20_30 <- (dist_mi > 20 & dist_mi <= 30)
within_30_40 <- (dist_mi > 30 & dist_mi <= 40)
within_40_60 <- (dist_mi > 40 & dist_mi <= 60)
within_60_80 <- (dist_mi > 60 & dist_mi <= 80)
within_80_100 <- (dist_mi > 60 & dist_mi <= 80)

# Spatial Decay
# Cutoff by element-wise multiplication & no self-spillover
# Following https://fmwww.bc.edu/repec/bocode/s/spgen.pdf
spatial_decay <- exp(-.02 * dist_mi)
spatial_decay <- spatial_decay * within_large

# Convert to sparse matrices
within <- as(within, "dgCMatrix")
within_large <- as(within_large, "dgCMatrix")
within_100 <- as(within_100, "dgCMatrix")
within_0_20 <- as(within_0_20, "dgCMatrix")
within_20_30 <- as(within_20_30, "dgCMatrix")
within_30_40 <- as(within_30_40, "dgCMatrix")
within_40_60 <- as(within_40_60, "dgCMatrix")
within_60_80 <- as(within_60_80, "dgCMatrix")
within_80_100 <- as(within_80_100, "dgCMatrix")
spatial_decay <- as(spatial_decay, "dgCMatrix")


# Simulation: Bias -------------------------------------------------------------
#' Simulates a set of DGPs with different kinds of exposure mappings
#' For each DGP, estimates using many different parametric forms + ring(s) method
#' Report on the bias of the treatment effect estimate

# %%
n_trials <- 100
doParallel::registerDoParallel(10)

dgp_types <- c(
  "Within 40mi." = "spill_within",
  "Within 80mi." = "spill_within_large",
  "Decay" = "spill_decay",
  "Within 40mi. (Additive)" = "spill_within_additive",
  "Within 80mi. (Additive)" = "spill_within_large_additive",
  "Decay (Additive)" = "spill_decay_additive"
)

estimation_types <- c(
  "TWFE (No Spillovers)" = "twfe",
  "Within 40mi." = "spill_within",
  "Within 80mi." = "spill_within_large",
  "Within 40mi. (Additive)" = "spill_within_additive",
  "Within 80mi. (Additive)" = "spill_within_large_additive",
  "Decay" = "spill_decay",
  "Decay (Additive)" = "spill_decay_additive",
  "Rings (0-20, 20-30, 30-40)" = "donuts_small",
  "Rings (0-20, 20-30, 30-40, 40-60, 60-80)" = "donuts",
  "Rings (0-20, 20-30, 30-40, 40-60, 60-80) (Additive)" = "donuts_additive"
)

if (RUN) {
  results <- foreach(i = 1:n_trials, .combine = "rbind") %dopar% {
    df <- sim_data_misspecification(counties)
    results_trial <- NULL

    for (y in dgp_types) {
      for (x in estimation_types) {
        if (x == "twfe") {
          formula <- as.formula(glue("y_{y} ~ treat_ind | state_county + year"))
        } else if (x == "donuts_small") {
          formula <- as.formula(glue("y_{y} ~ treat_ind + spill_0_20:post + spill_20_30:post + spill_30_40:post | state_county + year"))
        } else if (x == "donuts") {
          formula <- as.formula(glue(
            "y_{y} ~ treat_ind + spill_0_20:post + spill_20_30:post + spill_30_40:post ",
            "+ spill_40_60:post + spill_60_80:post | state_county + year"
          ))
        } else if (x == "donuts_additive") {
          formula <- as.formula(glue(
            "y_{y} ~ treat_ind + spill_0_20_additive:post + spill_20_30_additive:post + spill_30_40_additive:post ",
            "+ spill_40_60_additive:post + spill_60_80_additive:post | state_county + year"
          ))
        } else {
          formula <- as.formula(glue("y_{y} ~ treat_ind + ({x}:post) | state_county + year"))
        }

        treat_var <- as.character(glue("treat_ind"))

        te_hat <- estimate_treatment_effect(df, formula, treat_var)

        dgp <- names(dgp_types)[dgp_types == y]
        spec <- names(estimation_types)[estimation_types == x]


        results_trial <- bind_rows(results_trial, tibble(
          dgp = dgp, spec = spec, te_hat = te_hat
        ))
      }
    }

    return(results_trial)
  }
}

## Save results
if (RUN) save(results, file = "data/sim_misspecification.RData")

## Load results if not running simulation
if (!RUN) load("data/sim_misspecification.RData")

results_tbl <- results |>
  mutate(
    bias = 2 - te_hat
  ) |>
  summarize(
    # ci_lower = quantile(bias, probs = 0.05),
    # ci_upper = quantile(bias, probs = 0.95),
    mse = sum(bias^2) / n(),
    bias = mean(bias), 
		.by = c(spec, dgp)
  ) |>
  select(spec, dgp, bias, mse)


# %%
# Format into tex table
table_tex <- ""

for (spec in estimation_types) {
  # Create row
  row <- ""
  row_se <- ""

  # Specification Name
  spec_name <- names(estimation_types)[estimation_types == spec]
  row <- paste0(row, str_pad(spec_name, 55, "right"))
  row_se <- paste0(row_se, str_pad("", 55, "right"))

  for (dgp in dgp_types) {
    row <- paste0(row, "& ")
    row_se <- paste0(row_se, "& ")

    # DGP Names
    dgp_name <- names(dgp_types)[dgp_types == dgp]
    temp <- results_tbl |> filter(spec == spec_name & dgp == dgp_name)

    bias <- round(temp$bias, digits = 3) |>
      format(digits = 3, nsmall = 3) |>
      str_pad(8, side = "right")
    mse <- 
			paste0("[", format(round(temp$mse, digits = 3), nsmall = 3), "]") |> 
			str_pad(8, side = "right")

    row <- paste0(row, bias)
    row_se <- paste0(row_se, mse)
  }

  row <- paste0(row, "\\\\")
  if (spec != estimation_types[length(estimation_types)]) {
    row_se <- paste0(row_se, "\\\\")
  }
  # cli::cat_line(row, "\n", row_se, "\n")

  table_tex <- paste(table_tex, row, row_se, sep = "\n")
}

cli::cli_h1("Misspecification Bias Results")
cat(table_tex)

if (EXPORT) cat(table_tex, file = "tables/misspecification.tex")


# Simulation: MSPE of Spillovers -----------------------------------------------
#' Simulates a set of DGPs with different kinds of exposure mappings
#' For each DGP, estimates using many different parametric forms + ring(s) method
#' Estimate spillover effect for control units and determine mspe of spillover

# %%
n_trials <- 1000
doParallel::registerDoParallel(10)

dgp_types <- c(
  "Within 40mi." = "spill_within",
  "Within 80mi." = "spill_within_large",
  "Decay" = "spill_decay",
  "Within 40mi. (Additive)" = "spill_within_additive",
  "Within 80mi. (Additive)" = "spill_within_large_additive",
  "Decay (Additive)" = "spill_decay_additive"
)

estimation_types <- c(
  "TWFE (No Spillovers)" = "twfe",
  "Within 40mi." = "spill_within",
  "Within 80mi." = "spill_within_large",
  "Within 40mi. (Additive)" = "spill_within_additive",
  "Within 80mi. (Additive)" = "spill_within_large_additive",
  "Decay" = "spill_decay",
  "Decay (Additive)" = "spill_decay_additive",
  "Rings (0-20, 20-30, 30-40)" = "donuts_small",
  "Rings (0-20, 20-30, 30-40, 40-60, 60-80)" = "donuts",
  "Rings (0-20, 20-30, 30-40, 40-60, 60-80) (Additive)" = "donuts_additive"
)

if (RUN) {
  results_mspe <- foreach(i = 1:n_trials, .combine = "rbind") %dopar% {
    df <- sim_data_misspecification(counties)
    results_trial <- NULL

    for (y in dgp_types) {
      for (x in estimation_types) {
        if (x == "twfe") {
          df <- df |>
            mutate(te_spill_hat = 0)
        } else if (x == "donuts_small") {
          formula <- as.formula(glue("y_{y} ~ treat_ind + post:spill_0_20 + spill_20_30:post + spill_30_40:post | state_county + year"))

          reg <- feols(formula, data = df, warn = FALSE, notes = FALSE)

          coef <- reg |> coefficients()
          b_0_20 <- coef[2]
          b_20_30 <- coef[3]
          b_30_40 <- coef[4]

          df <- df |>
            mutate(
              te_spill_hat = b_0_20 * spill_0_20 + b_20_30 * spill_20_30 + b_30_40 * spill_30_40
            )
        } else if (x == "donuts") {
          formula <- as.formula(glue(
            "y_{y} ~ treat_ind + post:spill_0_20 + spill_20_30:post + spill_30_40:post ",
            "+ spill_40_60:post + spill_60_80:post | state_county + year"
          ))

          reg <- feols(formula, data = df, warn = FALSE, notes = FALSE)

          coef <- reg |> coefficients()
          b_0_20 <- coef[2]
          b_20_30 <- coef[3]
          b_30_40 <- coef[4]
          b_40_60 <- coef[5]
          b_60_80 <- coef[6]

          df <- df |>
            mutate(
              te_spill_hat = b_0_20 * spill_0_20 + b_20_30 * spill_20_30 + b_30_40 * spill_30_40 + b_40_60 * spill_40_60 + b_60_80 * spill_60_80
            )
        } else if (x == "donuts_additive") {
          formula <- as.formula(glue(
            "y_{y} ~ treat_ind + post:spill_0_20_additive + spill_20_30_additive:post + spill_30_40_additive:post ",
            "+ spill_40_60_additive:post + spill_60_80_additive:post | state_county + year"
          ))

          reg <- feols(formula, data = df, warn = FALSE, notes = FALSE)

          coef <- reg |> coefficients()
          b_0_20 <- coef[2]
          b_20_30 <- coef[3]
          b_30_40 <- coef[4]
          b_40_60 <- coef[5]
          b_60_80 <- coef[6]

          df <- df |>
            mutate(
              te_spill_hat = b_0_20 * spill_0_20_additive + b_20_30 * spill_20_30_additive + b_30_40 * spill_30_40_additive + b_40_60 * spill_40_60_additive + b_60_80 * spill_60_80_additive
            )
        } else {
          formula <- as.formula(glue("y_{y} ~ treat_ind + ({x}:post) | state_county + year"))

          reg <- feols(formula, data = df, warn = FALSE, notes = FALSE)
          coef <- reg |> coefficients()
          b <- coef[2]

          df <- df |>
            mutate(
              te_spill_hat = b * !!rlang::sym(x)
            )
        }

        # Calculate MSPE and normalize

        true_spill <- df |>
          filter(year == 2019 & treat == 0) |>
          pull(glue("te_{y}"))
        estimated_spill <- df |>
          filter(year == 2019 & treat == 0) |>
          pull("te_spill_hat")

        mspe <- sum((true_spill - estimated_spill)^2)
        total_var <- sum((true_spill^2))

        normalized_mspe <- mspe / total_var


        dgp <- names(dgp_types)[dgp_types == y]
        spec <- names(estimation_types)[estimation_types == x]

        results_trial <- bind_rows(results_trial, tibble(
          dgp = dgp, spec = spec, normalized_mspe = normalized_mspe
        ))
      }
    }

    return(results_trial)
  }
}

if (RUN) save(results_mspe, file = here::here("data/sim_misspecification_mspe.RData"))
if (!RUN) load(here::here("data/sim_misspecification_mspe.RData"))

# %%
results_tbl <- results_mspe |>
  summarize(
    # Makes 1 the best, 0 the worst
    percent_explained = 1 - mean(normalized_mspe),
		.by = c(spec, dgp)
  ) |>
  select(spec, dgp, percent_explained) |>
  pivot_wider(names_from = dgp, values_from = percent_explained) |>
  # Order columns
  select(spec, names(dgp_types)) |>
  # Order rows
  arrange(match(spec, names(estimation_types)))

# %% 
tab_mspe = results_tbl |>
  gt::gt() |>
  gt::fmt_number(
    columns = 2:7,
    decimals = 3
  ) |>
  extract_body()
	
cat(tab_mspe)
cat("\n\n")
if (EXPORT) {
	cat(tab_mspe, file = here::here("tables/misspecification_mspe.tex"))
}

# %%
tab_mspe_pct = results_tbl |>
  gt::gt() |>
  gt::fmt_percent(
    columns = 2:7,
    decimals = 1
  ) |>
  extract_body() 

cat(tab_mspe_pct)
cat("\n\n")
if (EXPORT) {
	cat(tab_mspe_pct, file = here::here("tables/misspecification_mspe_percent.tex"))
}
