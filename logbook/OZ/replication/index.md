---
title: "Replication of Chen, Glaeser, and Wessel (2023)"
---


::: {.cell}

```{.r .cell-code}
library(tidyverse, warn.conflicts = FALSE)
library(arrow)
library(here)
library(fixest)
library(broom)
```
:::


## Load Data 


::: {.cell}

```{.r .cell-code}
zillow_data <- read_feather(here("data/OZ/tracts_data.feather"))
zillow_data <- zillow_data |> filter(year != 2020)
zillow_data$tract <- factor(zillow_data$tract)
zillow_data$treat_x_post <- zillow_data$treatment * zillow_data$post
```
:::


## Classic Difference-in-Differences 


::: {.cell}

```{.r .cell-code}
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
```
:::


### No Covariates


::: {.cell}

```{.r .cell-code}
zillow_model <- fit_did_fixest(
  data = zillow_data,
  fmla = annual_change ~ treatment_and_post | year + tract,
  pretest_fmla = annual_change ~ 
    i(treatment, i.year, ref = FALSE, ref2 = 2017) | 
    year + tract,
  pretest_cols = c("treatment::TRUE:year::2014", "treatment::TRUE:year::2015", "treatment::TRUE:year::2016")
)
```

::: {.cell-output .cell-output-stdout}

```
Wald test, H0: joint nullity of treatment::TRUE:year::2014, treatment::TRUE:year::2015 and treatment::TRUE:year::2016
 stat = 0.101862, p-value = 0.958974, on 3 and 101,340 DoF, VCOV: Clustered (state_abbr).
```


:::

```{.r .cell-code}
zillow_model
```

::: {.cell-output .cell-output-stdout}

```
$model_pretest
OLS estimation, Dep. Var.: annual_change
Observations: 121,620 
Fixed-effects: year: 6,  tract: 20,270
Standard-errors: Clustered (state_abbr) 
                           Estimate Std. Error  t value Pr(>|t|) 
treatment::TRUE:year::2014 0.010515   0.253057 0.041551  0.96702 
treatment::TRUE:year::2015 0.024419   0.187918 0.129947  0.89713 
treatment::TRUE:year::2016 0.118437   0.235842 0.502190  0.61774 
treatment::TRUE:year::2018 0.395984   0.333319 1.188002  0.24045 
treatment::TRUE:year::2019 0.300182   0.198673 1.510937  0.13710 
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 7.09563     Adj. R2: 0.050177
                Within R2: 5.614e-5

$lh_pretest
$lh_pretest$stat
[1] 0.1018617

$lh_pretest$p
[1] 0.9589745

$lh_pretest$df1
[1] 3

$lh_pretest$df2
[1] 101340

$lh_pretest$vcov
[1] "Clustered (state_abbr)"


$model
OLS estimation, Dep. Var.: annual_change
Observations: 121,620 
Fixed-effects: year: 6,  tract: 20,270
Standard-errors: Clustered (state_abbr) 
                   Estimate Std. Error t value Pr(>|t|) 
treatment_and_post  0.30974   0.202635 1.52856  0.13268 
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RMSE: 7.09565     Adj. R2: 0.050209
                Within R2: 5.08e-5 
```


:::

```{.r .cell-code}
model_results <- tidy(zillow_model$model)
tau <- zillow_model$model$estimate
se <- zillow_model$model$std.error
tstat <- zillow_model$model$statistics
pval <- zillow_model$model$p.value
```
:::


### Sparse covariates


::: {.cell}

```{.r .cell-code}
zillow_model_cov_sp <- fit_did_fixest(
  data = zillow_data,
  fmla = annual_change ~ 1 + treatment * post + i(year, log_median_household_income, ref = 2017) + i(year, pct_white, ref = 2017) | year + tract,
  pretest_fmla = annual_change ~ 1 + i(treatment, i.year, ref = FALSE) + i(year, log_median_household_income, ref = 2017) + i(year, pct_white, ref = 2017) | year + tract,
  pretest_cols = c("treatment::TRUE:year::2014", "treatment::TRUE:year::2015", "treatment::TRUE:year::2016")
)
```

::: {.cell-output .cell-output-stdout}

```
Wald test, H0: joint nullity of treatment::TRUE:year::2014, treatment::TRUE:year::2015 and treatment::TRUE:year::2016
 stat = 0.227689, p-value = 0.877179, on 3 and 101,329 DoF, VCOV: Clustered (state_abbr).
```


:::

::: {.cell-output .cell-output-stderr}

```
The variables 'treatmentTRUE' and 'post' have been removed because of collinearity (see $collin.var).
```


:::
:::


### Full Covariates


::: {.cell}

```{.r .cell-code}
zillow_model_cov <- fit_did_fixest(
  data = zillow_data,
  fmla = annual_change ~ 1 + treatment * post + i(year, log_median_household_income, ref = 2017) + i(year, total_housing, ref = 2017) + i(year, pct_white, ref = 2017) + i(year, pct_higher_ed, ref = 2017) + i(year, pct_rent, ref = 2017) + i(year, pct_native_hc_covered, ref = 2017) + i(year, pct_poverty, ref = 2017) + i(year, pct_supplemental_income, ref = 2017) + i(year, pct_employed, ref = 2017) | year + tract,
  pretest_fmla = annual_change ~ 1 + treatment * year + i(year, log_median_household_income, ref = 2017) + i(year, total_housing, ref = 2017) + i(year, pct_white, ref = 2017) + i(year, pct_higher_ed, ref = 2017) + i(year, pct_rent, ref = 2017) + i(year, pct_native_hc_covered, ref = 2017) + i(year, pct_poverty, ref = 2017) + i(year, pct_supplemental_income, ref = 2017) + i(year, pct_employed, ref = 2017) | year + tract,
  pretest_cols = c("treatment::TRUE:year::2014", "treatment::TRUE:year::2015", "treatment::TRUE:year::2016")
)
```

::: {.cell-output .cell-output-stderr}

```
The variables 'treatmentTRUE' and 'year' have been removed because of collinearity (see $collin.var).
```


:::

::: {.cell-output .cell-output-stderr}

```
The variables 'treatmentTRUE' and 'post' have been removed because of collinearity (see $collin.var).
```


:::
:::


## Callway and Sant'Anna 


::: {.cell}

```{.r .cell-code}
library(did)
zillow_data$tract <- as.numeric(zillow_data$tract)

out <- att_gt(
  yname = "annual_change",
  tname = "year",
  idname = "tract",
  gname = "first_treat",
  xformla = ~ log_median_household_income + total_housing +
    pct_white + pct_higher_ed + pct_rent +
    pct_native_hc_covered + pct_poverty +
    pct_supplemental_income + pct_employed,
  data = zillow_data,
  bstrap = FALSE,
  cband = FALSE
)
did::aggte(out, type = "simple")
```

::: {.cell-output .cell-output-stdout}

```

Call:
did::aggte(MP = out, type = "simple")

Reference: Callaway, Brantly and Pedro H.C. Sant'Anna.  "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics, Vol. 225, No. 2, pp. 200-230, 2021. <https://doi.org/10.1016/j.jeconom.2020.12.001>, <https://arxiv.org/abs/1803.09015> 

    ATT    Std. Error     [ 95%  Conf. Int.] 
 0.3575        0.2314    -0.0961      0.8111 


---
Signif. codes: `*' confidence band does not cover 0

Control Group:  Never Treated,  Anticipation Periods:  0
Estimation Method:  Doubly Robust
```


:::
:::


## Doubly Robust 


::: {.cell}

```{.r .cell-code}
# library(DRDID)
# two_periods <- read_feather(here("data/OZ/two_periods.feather"))
#
# two_periods$treatment <- as.numeric(two_periods$treatment)
# two_periods$tract <- as.numeric(two_periods$tract)
#
# drdid_out <- drdid(yname = "annual_change", tname = "year", idname = "tract", dname = "treatment",
#              xformla= ~log_median_household_income + total_housing +
#     pct_white + pct_higher_ed + pct_rent +
#     pct_native_hc_covered + pct_poverty +
#     pct_supplemental_income + pct_employed,
#              data = two_periods, panel = TRUE)
```
:::


## Matched Pair Design 


::: {.cell}

```{.r .cell-code}
pair_diff_panel <- read_feather(here("data/OZ/pair_diff_panel.feather"))

pair_diff_panel <- pair_diff_panel |> filter(year != 2020)
pair_diff_panel$post <- pair_diff_panel$year == "2018"
pair_diff_panel$year <- relevel(factor(pair_diff_panel$year), ref = "2017")

model_diff <- fixest::feols(
  treated_minus_untreated ~ i(post, ref = FALSE) | pair_id,
  data = pair_diff_panel
)
```
:::
