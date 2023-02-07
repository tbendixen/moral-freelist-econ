# **Moral and cultural beliefs in behavioral economics**

`MoralizingGods_data set_V1.0_wIDs.csv`: Wave 1 and 2 game and demographic data, incl. original IDs for merging across waves (see: https://osf.io/epkbw/).

`FreeList_CERC_V0.1_FIN.csv`: Wave 1 free-list data (see: https://github.com/bgpurzycki/Evolution-of-Religion-and-Morality/tree/master/Wave%20I and https://github.com/tbendixen/Cross-cultural-Free-list-Project).

`FreeList_Wave2_BGD.csv`: Wave 2 free-list data (BGD only).

`data_prepped.csv`: Data prepped for analysis.

`data_prep.R`: Script for prepping data for analysis (outputs `data_prepped.csv`).

`codebook.txt`: Codebook for `data_prepped.csv`.

`model_fits.Rmd`: Notebook for fitting all models. Reports coefficient plots, trace rank plots, and posterior predictive checks for all models as well as some simple diagnostics (e.g., number of divergent transitions, if any) and model comparison metrics. Outputs `model_fits.pdf`.

`rag_result_plots.Rmd`: Reports plots for all supplementary models of the Random Allocation Game. Outputs `rag_result_plots.pdf`.

`dg_result_plots.Rmd`: Reports plots for all supplementary models of the Dictator Game. Outputs `dg_result_plots.pdf`.

`imp_plots.Rmd`: Reports plots from missing data imputations of all variables and across all main models. Outputs `imp_plots.pdf`.

`pre-registration.pdf`: Pre-registration document.

The **pre-registration** folder contains the following files:

`pre-reg_dg-sim.R` and `pre-reg_rag-sim.R`: Simulation and analysis script for the Dictator Game and Random Allocation Game models, respectively.

`pre-reg_ppc.R`: Compute and plot prior and posterior predictive checks and coefficient plots for focal parameters.

`pre-reg_dg-prob.R` and `pre-reg_rag-prob.R`: Compute and plot posterior predicted probabilities from the fitted Dictator Game and Random Allocation Game models, respectively.

-

The `.Rmd` files lists all `R` packages, their dependencies, and version number used.

The `.Rmd` files require the following software and their dependencies:
 - `R`: https://cran.r-project.org/bin/windows/base/
 - `RStudio`: https://www.rstudio.com/products/rstudio/
 - `R Markdown`: https://rmarkdown.rstudio.com/
 - `renv`: https://cran.r-project.org/web/packages/renv/index.html

Additional system details:
 - `version  R version 4.1.2 (2021-11-01)`
 - `os       Windows 10 x64 (build 19043)`
 - `system   x86_64, mingw32`
 - `rstudio  2022.07.2+576 Spotted Wakerobin (desktop)`
 - `pandoc   2.19.2 @ C:/PROGRA~1/Pandoc/ (via rmarkdown)`
