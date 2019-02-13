context("multi_venn")

# verified by eye that these match ezvenn figs; diff is main

test_that("multi_venn", {
  pref <- list(alone=names(contr.v)[1:2], all=names(contr.v))

  #One of p.cutoff or fdr.cutoff must be given.
  expect_error(multi_venn(res.df, prefix.lst = pref, prefix.lst = pref))

  mv <- multi_venn(tab=res.df, prefix.lst = pref, p.cutoff = 0.05, name=NA)
  expect_equal(ncol(mv), 3)
  expect_equal(as.numeric(mv["gene1",]), c(1,0,-1))

  mv <- multi_venn(res.df, prefix.lst = pref, p.cutoff = 0.05, logfc.cutoff = 1)
  expect_equal(sum(abs(mv)), 2)
  expect_equal(as.numeric(mv["gene1",]), c(1,0,-1))

  mv <- multi_venn(res.df, prefix.lst = pref, fdr.cutoff = 0.01)
  expect_equal(sum(abs(mv)), 2)
  expect_equal(as.numeric(mv["gene1",]), c(1,0,-1))

  mv.fn <- function() ret <- multi_venn(res.df, prefix.lst = pref[2], p.cutoff = 0.05, name=NA)
  vdiffr::expect_doppelganger("mvenn", mv.fn)

  mv.fn2 <- function() ret <- multi_venn(res.df, prefix.lst = pref[2], p.cutoff = 0.05, logfc.cutoff = 1, name=NA)
  vdiffr::expect_doppelganger("mvenn2", mv.fn2)

  res2 <- res.df[,-grep("logFC$", colnames(res.df))]
  mv.fn3 <- function() ret <- multi_venn(res2, prefix.lst = pref[2], p.cutoff = 0.05, name=NA)
  vdiffr::expect_doppelganger("mvenn_noLFC", mv.fn3)

  mv.fn4 <- function() ret <- multi_venn(res2, prefix.lst = pref[2], fdr.cutoff = 0.01, name=NA)
  vdiffr::expect_doppelganger("mvenn_noLFC2", mv.fn4)
})
