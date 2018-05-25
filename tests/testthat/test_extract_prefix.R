context("extract prefix")

test_that("nm", {
  expect_equal(extract_prefix(nm=c("x.y.p", "x.y.p", "x.p.fdr"), suffix="p"), c("x.y", "x.y"))
  expect_error(extract_prefix(nm=""))
  expect_error(extract_prefix(nm=NA))
  expect_error(extract_prefix(nm=NULL))
  expect_equal(extract_prefix(nm=".p", suffix="p"), "")
  expect_error(extract_prefix(nm="_p", suffix="p"))
})

test_that("suffix", {
  expect_error(extract_prefix(nm=c("x.y.p", "x.y.p", "x.p.fdr"), suffix=""))
  expect_error(extract_prefix(nm=c("x.y.p", "x.y.p", "x.p.fdr"), suffix=NA))
  expect_error(extract_prefix(nm=c("x.y.p", "x.y.p", "x.p.fdr"), suffix=NULL))
})

test_that("sep", {
  expect_error(extract_prefix(nm=c("x.y.p", "x.y.p", "x.p.fdr"), suffix="p", sep=".."))
  expect_error(extract_prefix(nm=c("x.y.p", "x.y.p", "x.p.fdr"), suffix="p", sep=NA))
  expect_true(is.na(extract_prefix(nm=c("x.y.p", "x.y.p", "p"), suffix="p", sep=NA)))
  expect_true(is.na(extract_prefix(nm="p", suffix="p", sep=NA)))
  expect_error(extract_prefix(nm=c("x.y.p", "x.y.p", "p"), suffix="p", sep="NA"))
  expect_error(extract_prefix(nm=c("x.y.p", "x.y.p", "x.p.fdr"), suffix="p", sep=NULL))
  expect_error(extract_prefix(nm=c("x.y..p", "x.y..p", "x.p..fdr"), suffix="p", sep=".."))

  expect_equal(extract_prefix(nm=c("x_p", "x_fdr"), suffix="p", sep="_"), "x")
  expect_error(extract_prefix(nm=c("x_p", "x_fdr"), suffix="p", sep="."))
})

test_that("integration", {
  expect_error(extract_prefix(nm="p", suffix="p"))
  expect_equal(extract_prefix(nm=c("x.p", "x.fdr"), suffix="p"), "x")
})
