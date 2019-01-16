# ezlimmaplot
R package for plotting bioinformatics results, especially those from the `ezlimma` package.

[![Travis-CI Build Status](https://travis-ci.org/jdreyf/ezlimmaplot.svg?branch=master)](https://travis-ci.org/jdreyf/ezlimmaplot)
[![Coverage Status](https://img.shields.io/codecov/c/github/jdreyf/ezlimmaplot/master.svg)](https://codecov.io/github/jdreyf/ezlimmaplot?branch=master)

## Install
On Windows, you need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

Install `ezlimmaplot` from GitHub using `remotes`  within R. You must install `remotes` if you haven't before. `ezlimmaplot` is intended for use with `ezlimma`, which depends on `limma`, so you should install these if you haven't before.
```
#if haven't already installed limma
source("http://bioconductor.org/biocLite.R")
biocLite("limma")

install.packages("remotes") #if haven't already installed remotes
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