context("prune mat")

test_that("does not corrupt", {
 expect_equal(M, prune_mat(M, verbose=FALSE))
 expect_equal(M, prune_mat(data.frame(M), verbose=FALSE))
})

test_that("symbols", {
  pm <- prune_mat(M, symbols = sym.v, only.symbols = TRUE, verbose=FALSE)
  expect_equal(rownames(pm), c(letters[1:9], "a", "~"))

  pmu <- prune_mat(M, symbols = sym.v, only.symbols = TRUE, unique.rows = TRUE, verbose=FALSE)
  expect_equal(rownames(pmu), c(letters[1:9], "~"))

  pmu <- prune_mat(M, symbols = sym.v, only.symbols = TRUE, unique.rows = TRUE, na.lab=c("~", "---", ""), verbose=FALSE)
  expect_equal(rownames(pmu), letters[1:9])
})

test_that("ntop", {
  pmnt <- prune_mat(M, symbols = sym.v, only.symbols = TRUE, ntop=5, verbose=FALSE)
  expect_equal(rownames(pmnt), letters[1:5])

  pmnt2 <- prune_mat(M, symbols = sym.v, only.symbols = TRUE, ntop=15, verbose=FALSE)
  expect_equal(rownames(pmnt2), c(letters[1:9], "a", "~"))
})
