# ezlimmaplot
R package for plotting bioinformatics results, especially those from the `ezlimma` package.

[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/ezlimmaplot/master.svg)](https://codecov.io/github/jdreyf/ezlimmaplot?branch=master)

## Install
On Windows, you should have [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

Install `ezlimmaplot` from GitHub using `remotes`  within R. You must install `remotes`, e.g. with install.packages("remotes"), if you haven't before. `ezlimmaplot` is intended for use with `ezlimma`, which depends on `limma`, so you should install these with instructions below if you haven't before.
```
#if haven't already installed limma
install.packages("BiocManager") #if haven't already installed BiocManager
library(BiocManager)
BiocManager::install("limma")

library(remotes)
remotes::install_github(repo="jdreyf/ezlimma", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
remotes::install_github(repo="jdreyf/Hitman", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
remotes::install_github(repo="jdreyf/ezlimmaplot", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
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

## Code of Conduct
This project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.
