context("plot_pwy")

# set plot=T to test more of the fcn, since teardown rm PDFs

test_that("returned object", {
  # tbl_graph is subclass of `igraph`
  # vertex a is most significant or tied
  expect_gte(V(pp)$EMY.chisq[match("a", V(pp)$name)], max(V(pp)$EMY.chisq[-match("a", V(pp)$name)]))
})

test_that("ntop & seed", {
  # ntop must be >= 2
  expect_error(plot_pwy(feat.tab = hm, G.pwy = gmt[[1]], stat.colnm = "EMY.chisq", annot.col = "symbol",
           gr=gr, name = NA, colorbar.nm = "chisq", ntop = 1, seed = 0, plot = F, alternative="greater"))
  pp2 <- plot_pwy(feat.tab = hm, G.pwy = gmt[[1]], stat.colnm = "EMY.chisq", annot.col = "symbol",
                  gr=gr, name = NA, colorbar.nm = "chisq", ntop = 2, seed = 0, plot = FALSE, alternative="greater")$gg.pwy
  expect_equal(names(V(pp2)), c("a", "b"))
  expect_equal(V(pp2)$symbol, c("A", "B"))
})

test_that("annot", {
  abc <- c("a", "b", "c")
  feat.tab[abc, "symbol"] <- c("abijabee3.4", "abijabee3-_4", "oh-not-so-Much")
  pp.ann <- plot_pwy(feat.tab = feat.tab, G.pwy = G.pwy, stat.colnm = "EMY.chisq", annot.col = "symbol",
                     gr=gr, name = NA, colorbar.nm = "chisq", ntop = 7, seed = 1, plot = FALSE, alternative="greater")$gg.pwy
  # only returns genes in pwy
  expect_equal(sort(V(pp.ann)$symbol), c("abijabee3-_4", "abijabee3.4", "oh-not-so-Much"))
})

test_that("non-NA annot overrides feature name", {
  feat.tab["a", "symbol"] <- "abijabee3.4"
  pp.ann2 <- plot_pwy(feat.tab = feat.tab, G.pwy = G.pwy, stat.colnm = "EMY.chisq", annot.col = "symbol",
                      gr=gr, name = NA, colorbar.nm = "chisq", ntop = 7, seed = 1, plot = FALSE, alternative="greater")$gg.pwy
  expect_equal(V(pp.ann2)$symbol[1], feat.tab["a", "symbol"])
})

test_that("analyte in G & not feat.tab & connected to top.nodes in plot, but nodes already connected", {
  ft <- feat.tab
  ft["b", "EMY.chisq"] <- NA
  pp.na <- plot_pwy(feat.tab = ft, G.pwy = G.pwy, stat.colnm = "EMY.chisq",
                    annot.col = "symbol", gr=gr, name = NA, colorbar.nm = "chisq", ntop = 2, seed = 1,
                    plot = FALSE, alternative="greater")$gg.pwy
  expect_true(is.na(V(pp.na)$EMY.chisq[V(pp.na)$name == "b"]))
})
