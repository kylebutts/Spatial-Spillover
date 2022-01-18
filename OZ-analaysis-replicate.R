library(tidyverse)
library(arrow)
library(here)
library(fixest)
library(broom)

# ---- Load Data ---------------------------------------------------------------

zillow_data <- read_feather(here("data/tracts_data.feather"))
zillow_data <- zillow_data |> filter(year != 2020)
zillow_data$tract <- factor(zillow_data$tract)
zillow_data$treat_x_post <- zillow_data$treatment * zillow_data$post

# ------------------------------------------------------------------------------
# Difference-in-Differences
# ------------------------------------------------------------------------------

# ---- Helper Function ---------------------------------------------------------

fit_did_fixest <- function(data, fmla, pretest_fmla, pretest_cols) {
  model_pretest <- fixest::feols(
    fml = pretest_fmla,
    data = data, cluster = ~state_abbr
  )

  pretest <- wald(model_pretest, keep = pretest_cols)

  model <- fixest::feols(
    fml = fmla, 
    data = data, cluster = ~state_abbr
  )

	return(list(
    model_pretest = model_pretest,
    lh_pretest = pretest,
    model = model
  ))
}

# ---- No Covariates -----------------------------------------------------------

zillow_model <- fit_did_fixest(
  data = zillow_data, 
  fmla = annual_change ~ treatment_and_post | year + tract,
  pretest_fmla = annual_change ~  i(treatment, i.year, ref = F, ref2 = 2017) | year + tract,
  pretest_cols = c("treatmentTRUE:year2014", "treatmentTRUE:year2015", "treatmentTRUE:year2016")
)
zillow_model


model_results <- tidy(zillow_model$model)
tau <- zillow_model$model$estimate 
se <- zillow_model$model$std.error
tstat <- zillow_model$model$statistics
pval <- zillow_model$model$p.value


# ---- Sparse Covariates -------------------------------------------------------

zillow_model_cov_sp <- fit_did_fixest(
  data = zillow_data,
  fmla = annual_change ~ 1 + treatment * post + i(year, log_median_household_income, ref = 2017) + i(year, pct_white, ref = 2017) | year + tract, 
  pretest_fmla = annual_change ~ 1 + i(treatment, i.year, ref = F) + i(year, log_median_household_income, ref = 2017) + i(year, pct_white, ref = 2017) | year + tract, 
  pretest_cols = c("treatmentTRUE:year2014", "treatmentTRUE:year2015", "treatmentTRUE:year2016")
)


# ---- Covariates --------------------------------------------------------------

zillow_model_cov <- fit_did_fixest(
  data = zillow_data,
  fmla = annual_change ~ 1 + treatment * post + i(year, log_median_household_income, ref = 2017) + i(year, total_housing, ref = 2017) + i(year, pct_white, ref = 2017) + i(year, pct_higher_ed, ref = 2017) + i(year, pct_rent, ref = 2017) + i(year, pct_native_hc_covered, ref = 2017) + i(year, pct_poverty, ref = 2017) + i(year, pct_supplemental_income, ref = 2017) + i(year, pct_employed, ref = 2017)  | year + tract, 
  pretest_fmla = annual_change ~ 1 + treatment * year + i(year, log_median_household_income, ref = 2017) + i(year, total_housing, ref = 2017) + i(year, pct_white, ref = 2017) + i(year, pct_higher_ed, ref = 2017) + i(year, pct_rent, ref = 2017) + i(year, pct_native_hc_covered, ref = 2017) + i(year, pct_poverty, ref = 2017) + i(year, pct_supplemental_income, ref = 2017) + i(year, pct_employed, ref = 2017) | year + tract, 
  pretest_cols = c("treatmentTRUE:year2014", "treatmentTRUE:year2015", "treatmentTRUE:year2016")
)

# ------------------------------------------------------------------------------
# Weighting Estimators
# ------------------------------------------------------------------------------

# ---- Callway and Sant'Anna ---------------------------------------------------

library(did)

zillow_data$tract <- as.numeric(zillow_data$tract)

out <- att_gt(
    yname = "annual_change",
    tname = "year",
    idname = "tract",
    gname = "first_treat",
    xformla = ~log_median_household_income + total_housing +
    pct_white + pct_higher_ed + pct_rent +
    pct_native_hc_covered + pct_poverty +
    pct_supplemental_income + pct_employed,
    data = zillow_data,
    bstrap = FALSE,
    cband = FALSE
)

did::aggte(out, type="simple")

# ---- Doubly Robust -----------------------------------------------------------

library(DRDID)
two_periods <- read_feather(here("data/two_periods.feather"))

two_periods$treatment <- as.numeric(two_periods$treatment)
two_periods$tract <- as.numeric(two_periods$tract)

drdid_out <- drdid(yname = "annual_change", tname = "year", idname = "tract", dname = "treatment",
             xformla= ~log_median_household_income + total_housing +
    pct_white + pct_higher_ed + pct_rent +
    pct_native_hc_covered + pct_poverty +
    pct_supplemental_income + pct_employed,
             data = two_periods, panel = TRUE)

# ------------------------------------------------------------------------------
# Matched Pair Design
# ------------------------------------------------------------------------------


pair_diff_panel <- read_feather(here("data/pair_diff_panel.feather"))

pair_diff_panel <- pair_diff_panel |> filter(year != 2020)
pair_diff_panel$post <- pair_diff_panel$year == "2018" 
pair_diff_panel$year <- relevel(factor(pair_diff_panel$year), ref = "2017")

model_diff <- fixest::feols(
  treated_minus_untreated ~ i(post, ref = F) | pair_id,
  data = pair_diff_panel
)


