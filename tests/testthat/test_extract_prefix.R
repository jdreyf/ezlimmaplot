context("extract prefix")


test_that("extract_prefix", {
  expect_error(extract_prefix(nm="p", suffix="p"))
  expect_true(is.na(extract_prefix(nm="p", suffix="p", sep=NA)))

  expect_equal(extract_prefix(nm=".p", suffix="p"), "")
  expect_equal(extract_prefix(nm=c("x.p", "x.fdr"), suffix="p"), "x")
  expect_error(extract_prefix(nm=c("x.p", "x.fdr"), suffix="p", sep=NA))

  expect_error(extract_prefix(nm=c("x_p", "x_fdr"), suffix="p", sep="."))
  expect_equal(extract_prefix(nm=c("x_p", "x_fdr"), suffix="p", sep="_"), "x")

  expect_equal(extract_prefix(nm=c("x.y.p", "x.y.p", "x.p.fdr"), suffix="p"), c("x.y", "x.y"))
})
