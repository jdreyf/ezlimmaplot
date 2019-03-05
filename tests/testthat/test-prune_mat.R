context("prune mat")

test_that("does not corrupt", {
 expect_equal(M, ezlimmaplot:::prune_mat(M, verbose=FALSE))
 expect_equal(M, ezlimmaplot:::prune_mat(data.frame(M), verbose=FALSE))
})

test_that("labrows", {
  pm <- ezlimmaplot:::prune_mat(M, labrows = sym.v, only.labrows = TRUE, verbose=FALSE)
  expect_equal(rownames(pm), c(letters[1:9], "a", "~"))

  pmu <- ezlimmaplot:::prune_mat(M, labrows = sym.v, only.labrows = TRUE, unique.rows = TRUE, verbose=FALSE)
  expect_equal(rownames(pmu), c(letters[1:9], "~"))

  pmu <- ezlimmaplot:::prune_mat(M, labrows = sym.v, only.labrows = TRUE, unique.rows = TRUE, na.lab=c("~", "---", ""),
                                 verbose=FALSE)
  expect_equal(rownames(pmu), letters[1:9])
})

test_that("ntop", {
  pmnt <- ezlimmaplot:::prune_mat(M, labrows = sym.v, only.labrows = TRUE, ntop=5, verbose=FALSE)
  expect_equal(rownames(pmnt), letters[1:5])

  pmnt2 <- ezlimmaplot:::prune_mat(M, labrows = sym.v, only.labrows = TRUE, ntop=15, verbose=FALSE)
  expect_equal(rownames(pmnt2), c(letters[1:9], "a", "~"))
})
