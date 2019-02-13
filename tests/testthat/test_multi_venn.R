context("multi_venn")

# verified by eye that these match ezvenn figs; diff is main title

test_that("multi_venn", {
  pref <- list(alone=names(contr.v)[1:2], all=names(contr.v))

  #One of p.cutoff or fdr.cutoff must be given
  expect_error(multi_venn(res.df, prefix.lst = pref, prefix.lst = pref))

  mv <- multi_venn(tab=res.df, prefix.lst = pref, p.cutoff = 0.05, name=NA)
  ezv <- ezvenn(tab=res.df, p.cutoff = 0.05, name=NA)
  #row names differ in ties
  expect_equal(as.numeric(rowSums(abs(mv))), as.numeric(rowSums(abs(ezv))))
  expect_equal(mv[rownames(ezv),], ezv)

  mv <- multi_venn(res.df, prefix.lst = pref, p.cutoff = 0.05, logfc.cutoff = 1, name=NA)
  ezv <- ezvenn(tab=res.df, p.cutoff = 0.05, logfc.cutoff = 1, name=NA)
  expect_equal(as.numeric(rowSums(abs(mv))), as.numeric(rowSums(abs(ezv))))
  expect_equal(mv[rownames(ezv),], ezv)

  mv <- multi_venn(res.df, prefix.lst = pref, p.cutoff = 0.05, logfc.cutoff = 1, name=NA)
  ezv <- ezvenn(tab=res.df, p.cutoff = 0.05, logfc.cutoff = 1, name=NA)
  expect_equal(as.numeric(rowSums(abs(mv))), as.numeric(rowSums(abs(ezv))))
  expect_equal(mv[rownames(ezv),], ezv)

  #vdiffr
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
