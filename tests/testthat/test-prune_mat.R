context("prune mat")

test_that("does not corrupt", {
  expect_equal(setNames(rownames(M), nm=rownames(M)), prune_mat(M, verbose=FALSE))
  expect_equal(setNames(rownames(M), nm=rownames(M)), prune_mat(data.frame(M), verbose=FALSE))
})

test_that("labrows", {
  pm <- prune_mat(M, labrows = sym.v, only.labrows = TRUE, verbose=FALSE)
  expect_equal(pm, setNames(c(letters[1:9], "a", "~"), nm=names(pm)))

  # pm2 <- prune_mat(M, labrows = sym.v, verbose=FALSE)
  # expect_equal(setNames(pm2, nm=NULL), c(letters[1:9], "gene10", "a", "~", rownames(M)[-(1:12)]))

  pm3 <- prune_mat(M, labrows = rep("", nrow(M)), verbose=FALSE)
  expect_equal(setNames(pm3, nm=NULL), rep("", nrow(M)))

  pmu <- prune_mat(object=M, labrows = sym.v, only.labrows = TRUE, unique.rows = TRUE, verbose=FALSE)
  expect_equal(pmu, setNames(c(letters[1:9], "~"), nm=names(pmu)))

  pmu2 <- prune_mat(M, labrows = sym.v, only.labrows = TRUE, unique.rows = TRUE, na.lab=c("~", "---", ""),
                                 verbose=FALSE)
  expect_equal(setNames(pmu2, nm=NULL), letters[1:9])
})

test_that("ntop", {
  pmnt <- prune_mat(M, labrows = sym.v, only.labrows = TRUE, ntop=5, verbose=FALSE)
  expect_equal(setNames(pmnt, nm=NULL), letters[1:5])

  expect_warning(pmnt2 <- prune_mat(M, labrows = sym.v, only.labrows = TRUE, ntop=15, verbose=FALSE))
  expect_equal(setNames(pmnt2, nm=NULL), c(letters[1:9], "a", "~"))
})
