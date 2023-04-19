# Dynamic factor model for ILD

This repository constains code for fitting the Ornstein-Uhlenbeck factor (OUF) model along with a simulated dataset for illustrative purposes.

* Data simulated from a bivariate factor model with 4 longitudinal outcomes is contained in **sim_dat.csv**.

* The OUF model can be fit using the code in **fit_ouf.R**.  This code uses functions that are provided in **ouf_functions.R**, along with some additional functions written in C++.  The C++ functions can be installed in R using 
```
devtools::install_github("https://github.com/madelineabbott/FABOUP_fast.git", ref="main")
```
Initial parameter estimates can be supplied by the user or can be estimated empirically using **init_ouf.R**.

For questions, please contact mrabbott@umich.edu.
