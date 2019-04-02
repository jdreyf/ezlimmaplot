# cannot remove Rplots.pdf; maybe vdiffr makes it after testthat?
pdfs <- grep("\\.pdf$", dir(recursive = TRUE), value=TRUE)
unlink(pdfs)
