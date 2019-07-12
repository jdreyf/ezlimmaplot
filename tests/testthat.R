library(testthat)
library(ezlimmaplot)

Sys.setenv(VDIFFR_RUN_TESTS=FALSE)

test_check("ezlimmaplot")
