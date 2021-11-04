# **Differences-in-Differences with Spatial Spillovers**

<a href="https://raw.githubusercontent.com/kylebutts/Spatial-Spillover/master/latex/paper/auxiliary/Spillover.pdf">![Paper](https://img.shields.io/badge/Paper-%23323330.svg?style=for-the-badge&logo=adobe&logoColor=white)</a>
<a href="https://raw.githubusercontent.com/kylebutts/Spatial-Spillover/master/latex/slides/auxiliary/Spillover_slides.pdf">![Slides](https://img.shields.io/badge/Slides-%23323330.svg?style=for-the-badge&logo=adobe&logoColor=white)</a>


## Abstract 

Empirical work often uses treatment assigned following geographic boundaries. When the effects of treatment cross over borders, classical difference-in-differences estimation produces biased estimates for the average treatment effect. In this paper, I introduce a potential outcomes framework to model spillover effects and decompose the estimate's bias in two parts: (1) the control group no longer identifies the counterfactual trend because their outcomes are affected by treatment and (2) changes in treated units' outcomes reflect the effect of their own treatment status and the effect from the treatment status of ``close'' units. I propose estimation strategies that can remove both sources of bias and semi-parametrically estimate the spillover effects themselves including in settings with staggered treatment timing. To highlight the importance of spillover effects, I revisit analyses of three place-based interventions.

## Figures from Paper

Figure 1 -- Comparison of Single vs. Multiple Rings Estimation of Spillover Effects

![Within vs. Rings](https://raw.githubusercontent.com/kylebutts/Spatial-Spillover/master/figures/figure-rings_v_within.png)

<caption>An example figure that describes the advantages of multiple rings when estimating spillover effects.</caption>

Figure 3 -- Direct and Spillover Effects of Community Health Centers

![Revisiting Bailey and Goodman-Bacon (2016)](https://raw.githubusercontent.com/kylebutts/Spatial-Spillover/master/figures/figure-chc-es_combined.png)

<caption>Revisiting Bailey and Goodman-Bacon (2016) to see if results are robust to spillovers (they are) and to see if the effects of community health centers extended to other counties (they do not).</caption>



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
