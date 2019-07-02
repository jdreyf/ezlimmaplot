context(sub_labels)
test_that("Labels replace names", {

  month.days<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  x <- c(1:10)
  y <- letters
  # Test for checking length when integers is available to be replaced.
  #expect_length (length(sub_labels(month.days,x,'-')), 1)
  #Test for known output
  expect_equal(sub_labels(month.days,x,'-'), (c(1 :10)))
  #Test for checking the length when characters are present to be replaced.
  #expect_length ((sub_labels(month.days,y,'-')), 26)
})

  test_that("When x includes negative numbers", {

    month.days <- c(31, 28, 31, 30, 31, 30,31,31,30,31,30,31)
    x <- c(-1,NA,2,NA,3,NA,4,NA,NA,5,6,7,8,9,10,NA,NA,NA,NA,NA,11,12)
    y <- letters
    # Test for checking length when integers is available to be replaced.
    #expect_length (length(sub_labels(month.days,x,'-')), 1)
    expect_warning (sub_labels(month.days,x,'-'))


    #Test for known output
    expect_equal(sub_labels(month.days,y,'-'), y)
    #Test for checking the length when characters are present to be replaced.
    #expect_length ((sub_labels (month.days, y,'-')), 26)

})
  test_that("x should have some value except NA", {

    month.days <- c(31, 28, 31, 30, 31, 30,31,31,30,31,30,31)
    x <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

    #Test for any error
    expect_that(sub_labels(month.days, x, na.lab = c("---","")),throws_error())
    # Test that all the arguments are being used.
    expect_that(sub_labels(month.days, x), throws_error())
    # Test that all the arguments are defined
    expect_that(sub_labels(month.days, x, na.lab), throws_error())
  })

  test_that("x should have some value except NA", {

    month.days <- c(31, 28, 31, 30, 31, 30,31,31,30,31,30,31)
    x <- c(NA, NA, NA, NA, 1)
    t <- c(31, 28, 31, 30, 1, 1, 1, 1, 1, 1, 1, 1)
    names(t) <- month.days
    #Test to check that we receive the named vector
    expect_named(sub_labels(month.days, x, na.lab = c("---","")))
    #Test to check that we receive the known value
    expect_equal(sub_labels(month.days, x, na.lab = c("---","")), t)
    #Test to check the replacement length
    expect_that(sub_labels(month.days, x, na.lab = c("1","")), throws_error)

  })

  test_that("x should have some value except NA", {

    month.days <- c(31, 28, 31, 30, 31, 30,31,31,30,31,30,31)
    x <- c(NA,-1, -2)
    t <- c(31, -1, -2, -1, -2, -1, -2, -1, -2, -1, -2, -1)
    names(t) <- month.days
    r <- c(31, 28, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2)
    names(r) <- month.days
    s <- c(31, -1, 31, -1, -1, -1, -1, -1, -1, -1, -1, -1)
    names(s) <- month.days

    #Test to check that we receive the known value
    expect_equal(sub_labels(month.days, x, na.lab = c("1","")), t)

    expect_equal(sub_labels(month.days, x, na.lab = c("-1","")), r)

    expect_equal(sub_labels(month.days, x, na.lab = c("-2","")), s)
  })

