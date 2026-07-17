# Little evidence for meaningful additive intervention effect heterogeneity of resistance training upon strength and hypertrophy: Evidence from causal inference and meta-analytic variance comparisons in a large dataset of randomised controlled trials

## Abstract
Resistance training studies consistently demonstrate substantial variation in observed outcomes after interventions between participants. This variation is commonly interpreted as evidence that individuals differ meaningfully in their responses to resistance training interventions. However, observed variation in outcomes does not necessarily imply heterogeneity in the underlying causal effects of an intervention. Rather, observed outcomes reflect multiple sources of variation, including between-participant differences, within-participant variability, measurement error, and potentially true intervention effect heterogeneity. Within the potential outcomes framework, heterogeneity of intervention effects is fundamentally a question about the variance of causal effects. Although directly estimating this variance generally requires specialised study designs such as randomised replicated crossover trials, traditional randomised controlled trials (RCTs) comparing an intervention against an inert control can provide indirect evidence regarding intervention effect heterogeneity through comparisons of outcome variances. Here, I outline the causal framework underpinning variance comparisons under additive intervention effects, discuss the assumptions required for causal interpretation, and demonstrate how meta-analytic models can synthesise evidence across studies to compare variances while accounting for the ubiquitous relationship between means and variances. Applying these methods to a large (102 studies, ) meta-analytic dataset of RCTs of resistance training interventions reveals little evidence that interventions meaningfully increase outcome variability beyond that expected from the increase in mean outcomes themselves. Although the conclusions drawn necessarily depend upon assumptions regarding additivity of intervention effects, that controls are inert, and the correlation between potential outcomes, these assumptions appear reasonable for resistance training studies in healthy participants. Taken together, current evidence provides little support for the existence of substantial intervention effect heterogeneity in resistance training adaptations. As such, recommendations, policies, and guidelines for resistance training interventions can be simplified and it can be justifiably assumed that the average intervention effects estimated in adequately powered and precise RCTs are constant and experienced by all undergoing the intervention.

# Supplementary materials
Any supplementary materials including analyses and plots are available [here](https://jamessteeleii.github.io/ma_variation_irv_rt/manuscript/supplementary).

## Reproducibility
This repository contains the necessary files and code to reproduce the analyses, figures, and the manuscript. 

## Usage
To reproduce the analyses, you will need to have R (https://cran.r-project.org/) and RStudio (https://www.rstudio.com/products/rstudio/download/#download) installed on your computer.

To help with reproducibility, this project uses the `renv` R package (see https://rstudio.github.io/renv/articles/renv.html). With `renv`, the state of this R project can be easily loaded as `renv` keeps track of the required R packages (including version), and (if known) the external source from which packages were retrieved (e.g., CRAN, Github). With `renv`, packages are installed to a project specific library rather than your user or system library. The `renv` package must be installed on your machine before being able to benefit from its features. The package can be installed using the following command:

``` r
install.packages("renv")
```

Once you have `renv` installed, you can get a copy of this repository on your machine by clicking the green Code button then choose Download zip. Save to your machine and extract. After extraction, double click the `ma_variation_irv_rt.Rproj` file in the root directory. This will automatically open RStudio. This will ensure all paths work on your system as the working directory will be set to the location of the `.Rproj` file. Upon opening, RStudio will recognize the `renv` files and you will be informed that the project library is out of sync with the lockfile. At shown in the console pane of RStudio, running `renv::restore()` will install the packages recorded in the lockfile. This could take some time depending on your machine and internet connection.

## Targets analysis pipeline

This project also uses a function based analysis pipeline using
[`targets`](https://books.ropensci.org/targets/). Instead of script based pipelines the `targets` package makes use of functions applied to targets specified within the pipeline. The targets can be viewed in the `_targets.R` file, and any user defined functions are available in `R/functions.r`.

You can view the existing targets pipeline by clicking [here](https://jamessteeleii.github.io/ma_variation_irv_rt/targets_pipeline.html).

Useful console functions:

- `tar_edit()` opens the make file
- `tar_make()` to run targets
- `tar_visnetwork()` to view pipeline

## Software and packages used

The [`grateful`](https://pakillo.github.io/grateful/index.html) package was used to create citations to all software and packages used in the analysis. The `grateful` report can be viewed by downloading by clicking [here](https://jamessteeleii.github.io/ma_variation_irv_rt/grateful-report.html).

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
  [cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

