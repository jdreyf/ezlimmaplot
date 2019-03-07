context("plot_by_grp")

test_that("pbg", {
  pbg1 <- plot_by_grp(object=M["gene1",], grp=pheno$grp, name=NA, ylab="Log2 abundance", add.se=TRUE, main="gene 1", dotsize = 0.7)
  vdiffr::expect_doppelganger("gene1", pbg1)

  pbg1b <- plot_by_grp(object=M["gene1",], grp=pheno$grp, name=NA, ylab="Log2 abundance", main="gene 1", dotsize = 1, bins=100)
  vdiffr::expect_doppelganger("gene1-nbins100", pbg1b)

  pbg2 <- plot_by_grp(object=M["gene2",], grp=pheno$grp, name=NA, type="box", main="gene 2", xlab="grp")
  vdiffr::expect_doppelganger("gene2", pbg2)

  #pbg returns last row as ggp object
  pbg3 <- plot_by_grp(object=M[1:3,], grp=pheno$grp, name=NA, x.angle=15, manual.color = c("black", "orange"))
  vdiffr::expect_doppelganger("gene3", pbg3)

  #plot_by_group can error when object is a df, maybe b/c it doesn't drop columns in object[i,]
  pbg2df <- plot_by_grp(object=data.frame(M)["gene2",], grp=pheno$grp, name=NA, type="box", main="gene 2", xlab="grp")
  vdiffr::expect_doppelganger("gene2", pbg2df)
})

teardown(unlink("Rplots.pdf"))
