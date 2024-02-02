# Adding spillovers to Opportunity Zone analysis


Opportunity Zone analysis where I run Selected vs. Eligible but Not
Selected and have spillovers on adjacent census tracts. See if
differences in effect is due to identification strategy.

``` r
library(tidyverse, warn.conflicts = FALSE)
library(sf)
```

    Warning: package 'sf' was built under R version 4.3.1

``` r
library(units)
```

    Warning: package 'units' was built under R version 4.3.1

``` r
library(arrow)
library(here)
library(fixest)
library(broom)
```

Find distance to nearest treated:

``` r
geos <- read_sf(here("data/OZ/tracts/all_tracts_with_status.shp"))
geos <- geos |>
    st_point_on_surface() |>
    mutate(
        status = ifelse(is.na(status), "Ineligible", status),
        intptlat = as.numeric(intptlat),
        intptlon = as.numeric(intptlon),
    )
```

    Warning: st_point_on_surface assumes attributes are constant over geometries

    Warning in st_point_on_surface.sfc(st_geometry(x)): st_point_on_surface may not
    give correct results for longitude/latitude data

Code is inefficient because of memory concerns

``` r
selected <- geos[geos$status == "Selected", ]
geos$distance_to_nearest <- kfbmisc::st_nearest_distance_rcpp(geos, selected)[, 2]
```

    [1] "Distance in miles. To use kilometers use option `unit == 'km'`"

``` r
geos <- geos |>
    st_as_sf() |>
    st_drop_geometry() |>
    mutate(tract = as.character(geoid)) |>
    select(tract, status, distance_to_nearest)
```

Loading the data and merging with spatial data

``` r
tracts <- read_feather(here("data/OZ/tracts_data.feather")) |> 
    filter(year != 2020) |> 
    left_join(geos, by = "tract")
```

## Regression Analysis

### Selected vs. Eligible

``` r
(eligible <- feols(
    annual_change ~ treatment_and_post + 
        i(year, log_median_household_income, ref = 2017) + 
        i(year, total_housing, ref = 2017) + 
        i(year, pct_white, ref = 2017) + 
        i(year, pct_higher_ed, ref = 2017) + 
        i(year, pct_rent, ref = 2017) + 
        i(year, pct_native_hc_covered, ref = 2017) + 
        i(year, pct_poverty, ref = 2017) + 
        i(year, pct_supplemental_income, ref = 2017) + 
        i(year, pct_employed, ref = 2017) | 
        year + tract,
    data = tracts, cluster = ~state_abbr
))
```

    OLS estimation, Dep. Var.: annual_change
    Observations: 121,620 
    Fixed-effects: year: 6,  tract: 20,270
    Standard-errors: Clustered (state_abbr) 
                                            Estimate Std. Error   t value
    treatment_and_post                      0.303266   0.166103  1.825776
    year::2014:log_median_household_income 10.311156   3.307763  3.117259
    year::2015:log_median_household_income  0.352815   0.632039  0.558217
    year::2016:log_median_household_income  1.241923   0.898994  1.381459
    year::2018:log_median_household_income -1.802364   0.899350 -2.004073
    year::2019:log_median_household_income -4.599137   1.183967 -3.884515
    year::2014:total_housing               -0.000483   0.000273 -1.767482
    year::2015:total_housing               -0.000045   0.000131 -0.341487
                                             Pr(>|t|)    
    treatment_and_post                     0.07385888 .  
    year::2014:log_median_household_income 0.00302404 ** 
    year::2015:log_median_household_income 0.57918742    
    year::2016:log_median_household_income 0.17328039    
    year::2018:log_median_household_income 0.05049443 .  
    year::2019:log_median_household_income 0.00030188 ***
    year::2014:total_housing               0.08324988 .  
    year::2015:total_housing               0.73416776    
    ... 38 coefficients remaining (display them with summary() or use argument n)
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    RMSE: 7.00077     Adj. R2: 0.07503 
                    Within R2: 0.026614

### Selected vs. Nearby

``` r
pair_diff_panel <- read_feather(here("data/OZ/pair_diff_panel.feather"))

pair_diff_panel <- pair_diff_panel |> 
    filter(year != 2020) |> 
    mutate(
        treatment_and_post = as.numeric(year == "2018"),
        year = relevel(factor(year), ref = "2017")
    )
```

``` r
(nearby <- fixest::feols(
    treated_minus_untreated ~ treatment_and_post | pair_id,
    data = pair_diff_panel
))
```

    OLS estimation, Dep. Var.: treated_minus_untreated
    Observations: 15,882 
    Fixed-effects: pair_id: 2,647
    Standard-errors: Clustered (pair_id) 
                       Estimate Std. Error t value  Pr(>|t|)    
    treatment_and_post 0.647793   0.245663 2.63691 0.0084153 ** 
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    RMSE: 9.44519     Adj. R2: -0.14431 
                    Within R2:  6.529e-4

## Spillover Analysis

``` r
tracts$within_half <- tracts$distance_to_nearest <= 1 / 2 & tracts$distance_to_nearest > 0
tracts$within_1 <- tracts$distance_to_nearest <= 1 & tracts$distance_to_nearest > 1 / 2
tracts$within_half_and_post <- as.numeric(tracts$within_half & tracts$year >= 2018)
tracts$within_1_and_post <- as.numeric(tracts$within_1 & tracts$year >= 2018)
```

``` r
(spill <- feols(
    annual_change ~ 
        treatment_and_post + within_half_and_post + within_1_and_post + 
        i(year, log_median_household_income, ref = 2017) + 
        i(year, total_housing, ref = 2017) + 
        i(year, pct_white, ref = 2017) + 
        i(year, pct_higher_ed, ref = 2017) + 
        i(year, pct_rent, ref = 2017) + 
        i(year, pct_native_hc_covered, ref = 2017) + 
        i(year, pct_poverty, ref = 2017) + 
        i(year, pct_supplemental_income, ref = 2017) + 
        i(year, pct_employed, ref = 2017) | 
        year + tract,
    data = tracts, cluster = ~state_abbr
))
```

    OLS estimation, Dep. Var.: annual_change
    Observations: 121,620 
    Fixed-effects: year: 6,  tract: 20,270
    Standard-errors: Clustered (state_abbr) 
                                            Estimate Std. Error   t value
    treatment_and_post                      0.178783   0.169177  1.056782
    within_half_and_post                   -1.056679   0.361845 -2.920253
    within_1_and_post                      -0.743020   0.192220 -3.865464
    year::2014:log_median_household_income 10.311156   3.307791  3.117234
    year::2015:log_median_household_income  0.352815   0.632044  0.558212
    year::2016:log_median_household_income  1.241923   0.899001  1.381447
    year::2018:log_median_household_income -1.747633   0.888325 -1.967336
    year::2019:log_median_household_income -4.544406   1.184462 -3.836684
                                             Pr(>|t|)    
    treatment_and_post                     0.29569047    
    within_half_and_post                   0.00523418 ** 
    within_1_and_post                      0.00032057 ***
    year::2014:log_median_household_income 0.00302426 ** 
    year::2015:log_median_household_income 0.57919053    
    year::2016:log_median_household_income 0.17328386    
    year::2018:log_median_household_income 0.05470491 .  
    year::2019:log_median_household_income 0.00035093 ***
    ... 40 coefficients remaining (display them with summary() or use argument n)
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    RMSE: 6.99994     Adj. R2: 0.07523 
                    Within R2: 0.026844

## Results Table

``` r
(table <- fixest::esttex(eligible, nearby, spill,
    dict = c(treatment_and_post = "Treat X Post", within_half_and_post = "< 1/2mi. X Post", within_1_and_post = "< 1mi. X Post"),
    keep = c("Treat X Post", "< 1/2mi. X Post", "< 1mi. X Post")
))
```

    \begin{table}[htbp]
       \caption{no title}
       \centering
       \begin{tabular}{lccc}
          \tabularnewline \midrule \midrule
          Dependent Variables: & annual\_change  & treated\_minus\_untreated   & annual\_change\\   
          Model:               & (1)             & (2)                         & (3)\\  
          \midrule
          \emph{Variables}\\
          Treat X Post         & 0.3033$^{*}$    & 0.6478$^{***}$              & 0.1788\\   
                               & (0.1661)        & (0.2457)                    & (0.1692)\\   
          < 1/2mi. X Post      &                 &                             & -1.057$^{***}$\\   
                               &                 &                             & (0.3618)\\   
          < 1mi. X Post        &                 &                             & -0.7430$^{***}$\\   
                               &                 &                             & (0.1922)\\   
          \midrule
          \emph{Fixed-effects}\\
          year                 & Yes             &                             & Yes\\  
          tract                & Yes             &                             & Yes\\  
          pair\_id             &                 & Yes                         & \\  
          \midrule
          \emph{Fit statistics}\\
          Observations         & 121,620         & 15,882                      & 121,620\\  
          R$^2$                & 0.22957         & 0.04642                     & 0.22975\\  
          Within R$^2$         & 0.02661         & 0.00065                     & 0.02684\\  
          \midrule \midrule
          \multicolumn{4}{l}{\emph{Signif. Codes: ***: 0.01, **: 0.05, *: 0.1}}\\
       \end{tabular}
    \end{table}

``` r
table = as.character(table)
table[10:15] |>
    stringr::str_replace("X", "$\\\\times$") |>
    paste0("\n") |>
    cat(file = here("tables/OZ_replication.tex"))

cat(paste0(table[10:15], "\n"), file = here("tables/OZ_replication.tex"))
```
