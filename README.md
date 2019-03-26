# ezlimmaplot
R package for plotting bioinformatics results, especially those from the `ezlimma` package.

[![Travis-CI Build Status](https://travis-ci.org/jdreyf/ezlimmaplot.svg?branch=master)](https://travis-ci.org/jdreyf/ezlimmaplot)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/ezlimmaplot/master.svg)](https://codecov.io/github/jdreyf/ezlimmaplot?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

## Install
On Windows, you need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

Install `ezlimmaplot` from GitHub using `remotes`  within R. You must install `remotes`, e.g. with install.packages("remotes", repos="https://cloud.r-project.org"), if you haven't before. `ezlimmaplot` is intended for use with `ezlimma`, which depends on `limma`, so you should install these with instructions below if you haven't before.
```
#if haven't already installed limma
install.packages("BiocManager", repos="https://cloud.r-project.org") #if haven't already installed BiocManager
library(BiocManager)
BiocManager::install("limma")

library(remotes)
remotes::install_github(repo="jdreyf/ezlimma", build_opts = c("--no-resave-data", "--no-manual"))
remotes::install_github(repo="jdreyf/ezlimmaplot", build_opts = c("--no-resave-data", "--no-manual"))
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
