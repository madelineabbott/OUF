# Dynamic factor model for ILD

This repository contains example code for fitting the Ornstein-Uhlenbeck factor (OUF) model along with code for simulating datasets. For more details, see https://arxiv.org/abs/2307.15681.

* Measurements of 4 longitudinal outcomes can be simulated from OU factor models with either 1, 2, or 3 latent factors using **generate_data/ouf_sim_dat.R**.  As an example, **generate_data/data/sim_dat_v1_g1.csv** contains data generated from an OU factor model with 2 latent factors.

* The OUF model can be fit using the code in **fit_model/fit_ouf.R**.  This code uses functions that are provided in **fit_model/ouf_functions.R**, along with some code written in C++ (see **fit_model/ou_precision_matrix.cpp**).  Additional C++ functions can be installed in R using 
```
devtools::install_github("madelineabbott/OUFgrad")
```
* Initial parameter estimates can be supplied by the user or can be estimated empirically using **fit_model/init_ouf.R**.

* AIC and BIC can be calculated using **fit_model/calc_aic_bic.R**

To fit the OUF model, users should begin with the **fit_ouf.R** file.  All other files that contain code used in the estimation algorithm are called within this main file.  For questions, please contact mrabbott@umich.edu.
