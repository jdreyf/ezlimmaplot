context("plot_by_grp")

# need to name plot to avoid vdiffr error "svglite only supports one page"
test_that("pbg", {
  pbg1 <- plot_by_grp(object=M["gene1",], grp=pheno$grp, ylab="Log2 abundance", add.se=TRUE, main="gene 1",
                      dotsize = 0.7, name="tmp")
  vdiffr::expect_doppelganger("pbg_gene1", pbg1)

  pbg1b <- plot_by_grp(object=M["gene1",], grp=pheno$grp, name="tmp", ylab="Log2 abundance", main="gene 1", dotsize = 1,
                       bins=100)
  vdiffr::expect_doppelganger("pbg_gene1-nbins100", pbg1b)

  pbg2 <- plot_by_grp(object=M["gene2",], grp=pheno$grp, name="tmp", type="box", main="gene 2", xlab="grp")
  vdiffr::expect_doppelganger("pbg_gene2", pbg2)

  #pbg returns last row as ggp object
  pbg3 <- plot_by_grp(object=M[1:3,], grp=pheno$grp, name="tmp", x.angle=15, manual.color = c("black", "orange"))
  vdiffr::expect_doppelganger("pbg_gene3", pbg3)

  #plot_by_group can error when object is a df, maybe b/c it doesn't drop columns in object[i,]
  pbg2df <- plot_by_grp(object=data.frame(M)["gene2",], grp=pheno$grp, name="tmp", type="box", main="gene 2", xlab="grp")
  vdiffr::expect_doppelganger("pbg_gene2_df", pbg2df)
})

test_that("pbg: non-visual tests",{
  pbg1 <- plot_by_grp(object=M["gene1",], grp=pheno$grp, ylab="Log2 abundance", add.se=TRUE, main="gene 1",
                      dotsize = 0.7, name="tmp")
  # Validating the location of the points
  expect_equal(pbg1$data["sample4", "Exprs"], M["gene1", "sample4"])
  # Validating the type of the point
  expect_equal(as.character(pbg1$data["sample4", "Group"]), pheno["sample4", "grp"])
  pbg2 <- plot_by_grp(object=M["gene2",], grp=pheno$grp, name="tmp", type="box", main="gene 2", xlab="grp")
  # Validating the location of the points
  expect_equal(pbg2$data["sample4", "Exprs"], M["gene2", "sample4"])
  # Validating the type of the point
  expect_equal(as.character(pbg2$data["sample4", "Group"]), pheno["sample4", "grp"])

  pbg1b <- plot_by_grp(object=M["gene1",], grp=pheno$grp, name="tmp", ylab="Log2 abundance", main="gene 1", dotsize = 1,
                       bins=100)
  expect_equal(pbg1b$data["sample4", "Exprs"], M["gene1", "sample4"])
  expect_equal(as.character(pbg1b$data["sample4", "Group"]), pheno["sample4", "grp"])
  pbg3 <- plot_by_grp(object=M[1:3,], grp=pheno$grp, name="tmp", x.angle=15, manual.color = c("black", "orange"))
  expect_equal(pbg3$data["sample4", "Exprs"], M[3, "sample4"])
  expect_equal(as.character(pbg3$data["sample4", "Group"]), pheno["sample4", "grp"])
  pbg2df <- plot_by_grp(object=data.frame(M)["gene2",], grp=pheno$grp, name="tmp", type="box", main="gene 2", xlab="grp")
  expect_equal(pbg2df$data["sample4", "Exprs"], M['gene2', "sample4"])
  expect_equal(as.character(pbg2df$data["sample4", "Group"]), pheno["sample4", "grp"])

  pbg3 <- plot_by_grp(object=M["gene1",], grp=pheno$grp, name="tmp", ylab="Log2 abundance", main="gene 1",  type="violin")
  expect_equal(pbg3$data["sample4", "Exprs"], M["gene1", "sample4"])

  pbg4 <- plot_by_grp(object=M["gene1",], grp=pheno$grp, name="tmp", ylab="Log2 abundance", main="gene 1", type="bar", add.dot=TRUE)
  expect_equal(pbg4$data[pbg4$data$Group=="First3",] |> dplyr::pull(Mean), mean(M["gene1", pheno$grp=="First3"]))
  expect_equal(pbg4$data[pbg4$data$Group=="First3",] |> dplyr::pull(SE), sd(M["gene1", pheno$grp=="First3"])/sqrt(sum( pheno$grp=="First3")))
})

