context("sub_labels")

test_that("Labels replace names", {
  month.days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  x <- 1:12
  y <- letters[x]
  expect_equal(sub_labels(month.days, x ,'-'), x)
  expect_equal(sub_labels(month.days, y, '-'), y)
})

test_that("When x includes negative numbers", {
  n <- 22
  month.days <- rnorm(n=n)
  x <- c(-1,NA,2,NA,3,NA,4,NA,NA,5,6,7,8,9,10,NA,NA,NA,NA,NA,11,12)
  y <- letters[1:n]
  # Test for checking length when integers is available to be replaced.
  #expect_length (length(sub_labels(month.days,x,'-')), 1)
  expect_error(sub_labels(x, month.days, '-'))
  #Test for known output
  expect_equal(sub_labels(month.days, y, '-'), y)
})

test_that("x should have some value except NA", {
  month.days <- c(31, 28, 31, 30, 31, 30,31,31,30,31,30,31)
  x <- rep(NA, length(month.days))

  #Test for any error
  expect_error(sub_labels(month.days, x, na.lab = c("--", 31)))
  # Test that all the arguments are defined
  expect_error(sub_labels(month.days, x, na.lab))
})

test_that("x should have some value except NA", {
  month.days <- c(31, 28, 31, 30, 31, 30,31,31,30,31,30,31)
  x <- rep(NA, length(month.days))
  x[1] <- 1
  t <- c(31, 28, 31, 30, 1, 1, 1, 1, 1, 1, 1, 1)
  names(t) <- month.days
  #Test to check that we receive the named vector
  expect_named(sub_labels(month.days, x, na.lab = c("---","")))
  #Test to check that we receive the known value
  expect_equal(sub_labels(month.days, t, na.lab = c("---","")), t)
  #Test to check the replacement length
  expect_error(sub_labels(month.days, x, na.lab = c(31, "")))
  #all NA
  expect_equal(sub_labels(month.days, rep(NA, length(month.days))), setNames(month.days, nm=month.days))
})
