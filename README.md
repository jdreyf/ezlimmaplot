# ezlimmaplot
R package for plotting bioinformatics results, especially those from the `ezlimma` package.

[![Travis-CI Build Status](https://travis-ci.org/jdreyf/ezlimmaplot.svg?branch=master)](https://travis-ci.org/jdreyf/ezlimmaplot)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/ezlimmaplot/master.svg)](https://codecov.io/github/jdreyf/ezlimmaplot?branch=master)

## Install
Install `ezlimmaplot` from GitHub using `devtools`  within R. You must install `devtools` if you haven't before. `ezlimmaplot` is intended for use with `ezlimma`, which depends on `limma`, so you should install these if you haven't before.
```
source("http://bioconductor.org/biocLite.R")
biocLite("limma") #if haven't already installed limma
install.packages("devtools") #if haven't already installed devtools
library(devtools)
install_github(repo="jdreyf/ezlimma", build_vignettes = TRUE)
devtools::install_github(repo="jdreyf/ezlimmaplot", build_vignettes = TRUE)
```

## Usage
The vignette presents a tutorial. To see the vignette:
```
library(limma)
library(ezlimma)
library(ezlimmaplot)
browseVignettes(package="ezlimmaplot")
```
and click on HTML.