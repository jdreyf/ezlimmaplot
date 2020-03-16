context("neighbor nms")

test_that("nn", {
  expect_equal(neighbor_nms(gr, "a"), c("b", "d"))
  expect_equal(neighbor_nms(gr, "d"), c("a", "b", "c"))
})
