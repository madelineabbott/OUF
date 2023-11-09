# Dynamic factor model for ILD

This repository contains example code for fitting the Ornstein-Uhlenbeck factor (OUF) model along with a simulated dataset for illustrative purposes. For more details, see https://arxiv.org/abs/2307.15681.

* Data simulated from a bivariate factor model with 4 longitudinal outcomes are contained in **sim_dat.csv**.

* The OUF model can be fit using the code in **fit_ouf.R**.  This code uses functions that are provided in **ouf_functions.R**, along with some additional functions written in C++.  The C++ functions can be installed in R using 
```
devtools::install_github("madelineabbott/FABOUP_fast")
```
* Initial parameter estimates can be supplied by the user or can be estimated empirically using **init_ouf.R**.

To fit the OUF model, users should begin with the **fit_ouf.R** file.  All other files are called within this main file.  For questions, please contact mrabbott@umich.edu.
