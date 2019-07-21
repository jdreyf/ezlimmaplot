context("sif2edgelist")

test_that("clean pwycommons SIF", {
  sif <- cbind(c(letters[c(2,1,2,2,3)], ""), "interacts-with", letters[c(1,3,1,3,1,1)])
  el.o <- t(apply(sif[,c(1,3)], 1, FUN=sort))
  el2 <- sif2edgelist(sif)
  expect_equal(el2, el.o[c(1,2,4,6),])
  el3 <- sif2edgelist(sif[c(1,2,4,6),])
  expect_equal(el3, el.o[c(1,2,4,6),])
  
  el4 <- sif2edgelist(sif, rm.ids = c("", "c", "water"))
  expect_equal(el4, el2[1,,drop=FALSE])
})