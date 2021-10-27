# **Differences-in-Differences with Spatial Spillovers**

## Abstract 

Empirical work often uses treatment assigned following geographic boundaries. When the effects of treatment cross over borders, classical difference-in-differences estimation produces biased estimates for the average treatment effect. In this paper, I introduce a potential outcomes framework to model spillover effects and decompose the estimate's bias in two parts: (1) the control group no longer identifies the counterfactual trend because their outcomes are affected by treatment and (2) changes in treated units' outcomes reflect the effect of their own treatment status and the effect from the treatment status of "close" units. I propose estimation strategies that can remove both sources of bias and semi-parametrically estimate the spillover effects themselves. I extend Callaway and Sant'anna (2020) to allow for event-study estimates that control for spillovers. To highlight the importance of spillover effects, I revisit analyses of three place-based interventions. 

## Example of Bias from Spillovers

![Example](https://raw.githubusercontent.com/kylebutts/Spatial-Spillover/master/figures/figure-map_te_spill_all.png?token=AEYJMX2LRDHTKLXC5G3SEW27ZLJSG)



## To Recreate Exhibits

Figure 1 -- Comparison of Single vs. Multiple Rings Estimation of Spillover Effects

- figure-rings_example.R

Figure 2 -- TVA Effective Sample and Spillover Variables

- tva-data-build.do
- tva-analysis-replicate.R

Table 1 -- Effects of Tennessee Valley Authority on Decadel Growth

- tva-data-build.do
- tva-analysis-replicate.R

Figure 3 -- Direct and Spillover Effects of Community Health Centers

- chc-analysis.R
