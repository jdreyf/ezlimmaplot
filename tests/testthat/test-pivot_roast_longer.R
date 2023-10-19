test_that("from helper", {
  expect_equal(nrow(pivot_roast_longer(rc2, prefix.v = "First3", direction.v = "Down")), 1)
  expect_equal(nrow(pivot_roast_longer(rc2, prefix.v = "First3", direction.v = "Mixed")), nrow(rc2))
  expect_equal(nrow(pivot_roast_longer(rc2[-3,], prefix.v = "First3", direction.v = "Down")), 0)
  expect_equal(pivot_roast_longer(rc2, prefix.v = "Last3", direction.v = "Down")$Count, 0)
  expect_equal(signif(pivot_roast_longer(rc2, prefix.v = "Last3", direction.v = "Down")$FDR, 2), 0.6)
  expect_equal(pivot_roast_longer(rc2, prefix.v = "Last3vsFirst3", direction.v = "Down")$Count[2], 2)
})
