#' ---
#' title: "Adding spillovers to Opportunity Zone analysis"
#' ---
#' Opportunity Zone analysis where I run Selected vs. Eligible but Not Selected and have spillovers on adjacent census tracts. See if differences in effect is due to identification strategy.

# %%
#| message: false
#| warning: false
library(tidyverse, warn.conflicts = FALSE)
library(sf)
library(units)
library(arrow)
library(here)
library(fixest)
library(broom)

#' Find distance to nearest treated:
# %%
geos <- read_sf(here("data/OZ/tracts/all_tracts_with_status.shp"))
geos <- geos |>
    st_point_on_surface() |>
    mutate(
        status = ifelse(is.na(status), "Ineligible", status),
        intptlat = as.numeric(intptlat),
        intptlon = as.numeric(intptlon),
    )

#' Code is inefficient because of memory concerns
# %%
selected <- geos[geos$status == "Selected", ]
geos$distance_to_nearest <- kfbmisc::st_nearest_distance_rcpp(geos, selected)[, 2]

geos <- geos |>
    st_as_sf() |>
    st_drop_geometry() |>
    mutate(tract = as.character(geoid)) |>
    select(tract, status, distance_to_nearest)

#' Loading the data and merging with spatial data
# %% 
tracts <- read_feather(here("data/OZ/tracts_data.feather")) |> 
    filter(year != 2020) |> 
    left_join(geos, by = "tract")



# Regression Analysis ----------------------------------------------------------
## Selected vs. Eligible -------------------------------------------------------
# %% 
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

## Selected vs. Nearby ---------------------------------------------------------
# %%
pair_diff_panel <- read_feather(here("data/OZ/pair_diff_panel.feather"))

pair_diff_panel <- pair_diff_panel |> 
    filter(year != 2020) |> 
    mutate(
        treatment_and_post = as.numeric(year == "2018"),
        year = relevel(factor(year), ref = "2017")
    )

# %% 
(nearby <- fixest::feols(
    treated_minus_untreated ~ treatment_and_post | pair_id,
    data = pair_diff_panel
))


# Spillover Analysis -----------------------------------------------------------
# %% 
tracts$within_half <- tracts$distance_to_nearest <= 1 / 2 & tracts$distance_to_nearest > 0
tracts$within_1 <- tracts$distance_to_nearest <= 1 & tracts$distance_to_nearest > 1 / 2
tracts$within_half_and_post <- as.numeric(tracts$within_half & tracts$year >= 2018)
tracts$within_1_and_post <- as.numeric(tracts$within_1 & tracts$year >= 2018)

# %%
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


# Results Table ----------------------------------------------------------------
# %%
(table <- fixest::esttex(eligible, nearby, spill,
    dict = c(treatment_and_post = "Treat X Post", within_half_and_post = "< 1/2mi. X Post", within_1_and_post = "< 1mi. X Post"),
    keep = c("Treat X Post", "< 1/2mi. X Post", "< 1mi. X Post")
))

# %% 
table = as.character(table)
table[10:15] |>
    stringr::str_replace("X", "$\\\\times$") |>
    paste0("\n") |>
    cat(file = here("tables/OZ_replication.tex"))

cat(paste0(table[10:15], "\n"), file = here("tables/OZ_replication.tex"))
