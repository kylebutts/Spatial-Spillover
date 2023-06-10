# Difference-in-Differences with Spatial Spillovers

[Kyle Butts](https://www.kylebutts.com/)<sup>1</sup>
<br>
<sup>1</sup>University of Colorado: Boulder

#### [Paper](https://arxiv.org/abs/2105.03737) | [Five-minute Summary](https://www.kylebutts.com/papers/spatial-spillovers/)


## Abstract

Empirical work often uses treatment assigned following geographic boundaries. When the effects of treatment cross over borders, classical difference-in-differences estimation produces biased estimates for the average treatment effect. In this paper, I introduce a potential outcomes framework to model spillover effects and decompose the estimate's bias in two parts: (1) the control group no longer identifies the counterfactual trend because their outcomes are affected by treatment and (2) changes in treated units' outcomes reflect the effect of their own treatment status and the effect from the treatment status of ``close'' units. I propose estimation strategies that can remove both sources of bias and semi-parametrically estimate the spillover effects themselves including in settings with staggered treatment timing. To highlight the importance of spillover effects, I revisit analyses of three place-based interventions.


## Citation

```
@article{butts2023difference,
  title={Difference-in-Differences Estimation with Spatial Spillovers},
  author={Butts, Kyle},
  journal={arXiv preprint arXiv:2105.03737},
  year={2023}
}
```


## Replication

**Figure 1:** Comparison of Single vs. Multiple Rings Estimation of Spillover Effects

- `figure-rings_example.R`

**Figure 2:** TVA Effective Sample and Spillover Variables

- `tva-data-build.do`
- `tva-analysis-replicate.R`

**Table 1:** Effects of Tennessee Valley Authority on Decadel Growth

- `tva-data-build.do`
- `tva-analysis-replicate.R`

**Table B1:** Effects of Opportunity Zones on Annual Home Price Growth

- `OZ-analaysis-replicate.R`
- `OZ-analaysis-spillovers.R`

**Figure C1:** Total and Spillover Effects of Community Health Centers

- `chc-analysis.R`


