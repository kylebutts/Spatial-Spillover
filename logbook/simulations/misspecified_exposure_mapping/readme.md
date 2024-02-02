# Misspecification of Exposure Mapping


This file performs a set of Monte Carlo simulations that estimate a set
of DGPs using various exposure mappings and estimates them using
exposure mappings and ring(s). The first simulation reports on average
bias found in the treatment effect estimate and the second simulation
reports on the ability to accurately predict spillover effects from
fitted values.

``` r
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
```

``` r
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
```

From data_prepare_counties.R

``` r
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
```

## Simulation: Bias

Simulates a set of DGPs with different kinds of exposure mappings For
each DGP, estimates using many different parametric forms + ring(s)
method Report on the bias of the treatment effect estimate

``` r
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
```

``` r
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
```

    ── Misspecification Bias Results ───────────────────────────────────────────────

``` r
cat(table_tex)
```


    TWFE (No Spillovers)                                   & 0.237   & 0.237   & -0.547  & 0.237   & 0.237   & 0.237   \\
                                                           & [0.076] & [0.076] & [0.321] & [0.076] & [0.076] & [0.076] \\
    Within 40mi.                                           & -0.021  & 0.198   & -0.626  & -0.021  & 0.163   & 0.129   \\
                                                           & [0.020] & [0.059] & [0.414] & [0.020] & [0.046] & [0.037] \\
    Within 80mi.                                           & -0.029  & -0.029  & -0.775  & -0.029  & -0.029  & -0.029  \\
                                                           & [0.022] & [0.022] & [0.622] & [0.022] & [0.022] & [0.022] \\
    Within 40mi. (Additive)                                & 0.026   & 0.206   & -0.614  & -0.022  & 0.164   & 0.129   \\
                                                           & [0.020] & [0.062] & [0.399] & [0.020] & [0.047] & [0.037] \\
    Within 80mi. (Additive)                                & 0.015   & 0.113   & -0.683  & -0.030  & -0.028  & -0.029  \\
                                                           & [0.021] & [0.034] & [0.490] & [0.022] & [0.022] & [0.022] \\
    Decay                                                  & 1.196   & 0.659   & -0.032  & 1.231   & 0.782   & 0.894   \\
                                                           & [1.468] & [0.468] & [0.034] & [1.553] & [0.647] & [0.835] \\
    Decay (Additive)                                       & -0.044  & 0.129   & -0.691  & -0.104  & 0.003   & -0.025  \\
                                                           & [0.023] & [0.038] & [0.500] & [0.032] & [0.021] & [0.021] \\
    Rings (0-20, 20-30, 30-40)                             & -0.021  & 0.198   & -0.626  & -0.021  & 0.163   & 0.129   \\
                                                           & [0.020] & [0.059] & [0.414] & [0.020] & [0.046] & [0.037] \\
    Rings (0-20, 20-30, 30-40, 40-60, 60-80)               & -0.029  & -0.029  & -0.775  & -0.029  & -0.029  & -0.029  \\
                                                           & [0.022] & [0.022] & [0.622] & [0.022] & [0.022] & [0.022] \\
    Rings (0-20, 20-30, 30-40, 40-60, 60-80) (Additive)    & 0.016   & 0.113   & -0.682  & -0.027  & -0.027  & -0.028  \\
                                                           & [0.021] & [0.034] & [0.489] & [0.022] & [0.022] & [0.022] 

``` r
if (EXPORT) cat(table_tex, file = "tables/misspecification.tex")
```

## Simulation: MSPE of Spillovers

Simulates a set of DGPs with different kinds of exposure mappings For
each DGP, estimates using many different parametric forms + ring(s)
method Estimate spillover effect for control units and determine mspe of
spillover

``` r
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
```

``` r
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
```

``` r
tab_mspe = results_tbl |>
  gt::gt() |>
  gt::fmt_number(
    columns = 2:7,
    decimals = 3
  ) |>
  extract_body()
    
cat(tab_mspe)
```

    TWFE (No Spillovers) & $0.000$ & $0.000$ & $0.000$ & $0.000$ & $0.000$ & $0.000$ \\ 
    Within 40mi. & $0.993$ & $0.254$ & $0.586$ & $0.854$ & $0.386$ & $0.558$ \\ 
    Within 80mi. & $0.394$ & $0.963$ & $0.843$ & $0.339$ & $0.715$ & $0.676$ \\ 
    Within 40mi. (Additive) & $0.851$ & $0.208$ & $0.512$ & $0.995$ & $0.406$ & $0.607$ \\ 
    Within 80mi. (Additive) & $0.456$ & $0.616$ & $0.702$ & $0.471$ & $0.983$ & $0.934$ \\ 
    Decay & $0.599$ & $0.821$ & $0.963$ & $0.524$ & $0.755$ & $0.818$ \\ 
    Decay (Additive) & $0.604$ & $0.563$ & $0.777$ & $0.636$ & $0.932$ & $0.984$ \\ 
    Rings (0-20, 20-30, 30-40) & $0.223$ & $0.070$ & $-0.081$ & $-0.039$ & $0.048$ & $-0.054$ \\ 
    Rings (0-20, 20-30, 30-40, 40-60, 60-80) & $0.144$ & $0.545$ & $0.019$ & $-0.094$ & $0.228$ & $-0.005$ \\ 
    Rings (0-20, 20-30, 30-40, 40-60, 60-80) (Additive) & $0.834$ & $0.574$ & $0.712$ & $0.976$ & $0.950$ & $0.950$ \\ 

``` r
cat("\n\n")
```

``` r
if (EXPORT) {
    cat(tab_mspe, file = here::here("tables/misspecification_mspe.tex"))
}
```

``` r
tab_mspe_pct = results_tbl |>
  gt::gt() |>
  gt::fmt_percent(
    columns = 2:7,
    decimals = 1
  ) |>
  extract_body() 

cat(tab_mspe_pct)
```

    TWFE (No Spillovers) & $0.0\%$ & $0.0\%$ & $0.0\%$ & $0.0\%$ & $0.0\%$ & $0.0\%$ \\ 
    Within 40mi. & $99.3\%$ & $25.4\%$ & $58.6\%$ & $85.4\%$ & $38.6\%$ & $55.8\%$ \\ 
    Within 80mi. & $39.4\%$ & $96.3\%$ & $84.3\%$ & $33.9\%$ & $71.5\%$ & $67.6\%$ \\ 
    Within 40mi. (Additive) & $85.1\%$ & $20.8\%$ & $51.2\%$ & $99.5\%$ & $40.6\%$ & $60.7\%$ \\ 
    Within 80mi. (Additive) & $45.6\%$ & $61.6\%$ & $70.2\%$ & $47.1\%$ & $98.3\%$ & $93.4\%$ \\ 
    Decay & $59.9\%$ & $82.1\%$ & $96.3\%$ & $52.4\%$ & $75.5\%$ & $81.8\%$ \\ 
    Decay (Additive) & $60.4\%$ & $56.3\%$ & $77.7\%$ & $63.6\%$ & $93.2\%$ & $98.4\%$ \\ 
    Rings (0-20, 20-30, 30-40) & $22.3\%$ & $7.0\%$ & $-8.1\%$ & $-3.9\%$ & $4.8\%$ & $-5.4\%$ \\ 
    Rings (0-20, 20-30, 30-40, 40-60, 60-80) & $14.4\%$ & $54.5\%$ & $1.9\%$ & $-9.4\%$ & $22.8\%$ & $-0.5\%$ \\ 
    Rings (0-20, 20-30, 30-40, 40-60, 60-80) (Additive) & $83.4\%$ & $57.4\%$ & $71.2\%$ & $97.6\%$ & $95.0\%$ & $95.0\%$ \\ 

``` r
cat("\n\n")
```

``` r
if (EXPORT) {
    cat(tab_mspe_pct, file = here::here("tables/misspecification_mspe_percent.tex"))
}
```
