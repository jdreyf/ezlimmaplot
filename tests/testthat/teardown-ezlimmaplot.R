# cannot remove Rplots.pdf; maybe vdiffr makes it after testthat?
pdfs <- grep("\\.pdf$", dir(recursive = TRUE), value=TRUE)
unlink(pdfs)

rplots.dir <- test_path("Rplots.pdf")
unlink(rplots.dir, force=TRUE)
