[![DOI](https://zenodo.org/badge/931844817.svg)](https://doi.org/10.5281/zenodo.17128110)

This repository contains code associated with **Hernández-Carrasco et al.** manuscript presenting *Generalised Dissimilarity Uniqueness Models (GDUM)*.

NOTE: An installation of the **gdmmTMB R package** is required to run the scripts in this repository. `gdmmTMB` can be installed from GitHub using:  

```r
devtools::install_github("dhercar/gdmmTMB")
```

### Contents
- `Simulation_study/Simulation_param_recover.R` — performs the simulation study.  
- `Case-study/` — contains scripts:  
  - `conceptual_fig.R` — recreates Fig. 1 
  - `CaseStudyMicrobial.R` — Microbial communities along a pH gradient (case study 1)
  - `CaseStudySW.R` — Plant community change along environmental gradients (case study 2)
