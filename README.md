# ADMUR
Ancient Demographic Modelling Using Radiocarbon

Tools to directly model underlying population dynamics using date datasets (radiocarbon and other) with a Continuous Piecewise Linear (CPL) model framework, and model comparison framework using BIC. Package also calibrates 14C samples, and generates Summed Probability Distributions (SPD).  CPL modelling directly estimates the most likely population trajectory given a dataset, using SPD simulation analysis to generate a Goodness-of-fit test for the best selected model.

Please contact a.timpson@ucl.ac.uk  in the first instance to make suggestions, report bugs or request help.

## Installation

``` r
# Installing and loading 'devtools' give you the ability to install direct from GitHub:
install.packages('devtools')
library(devtools)

# The ADMUR package can then be installed and loaded:
install_github('UCL/ADMUR')
library(ADMUR)
```

## Guide

Please refer to the vignette 'guide' for detailed support and examples.

``` r
help(AMDUR)
```

## References

<div id="ref-timpson-rstb.2020">

Timpson A., Barberena R., Thomas M. G., Mendez C., Manning K. 2020. "Directly modelling population dynamics in the South American "Arid Diagonal using 14C dates",Philosophical Transactions B. <https://doi.org/10.1098/rstb.2019.0723>.

</div>