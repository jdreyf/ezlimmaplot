context("plot_pwy")

test_that("returned object", {
  # tbl_graph is subclass of `igraph`
  # vertex a is most significant or tied
  expect_gte(V(pp)$EMY.z[match("A", V(pp)$name)], max(V(pp)$EMY.z[-match("A", V(pp)$name)]))
})

test_that("ntop & seed", {
  pp2 <- plot_pwy(feat.tab = hm, G.pwy = gmt[[1]], stat.colnm = "EMY.z", annot.col = "symbol",
                 gr=gr, name = NA, colorbar.nm = "z", ntop = 1, seed = 0, plot = T, alternative="greater")
  expect_equal(names(V(pp2)), c("A", "B"))
})

test_that("annot vdiffr", {
  abc <- c("a", "b", "c")
  V(gr)$name[match(abc, V(gr)$name)] <- rownames(feat.tab)[match(abc, rownames(feat.tab))] <-
    G.pwy$genes[match(abc, G.pwy$genes)] <-
    c("abijabee3.4", "abijabee3-_4", "oh-not-so-Much")

  pp.ann <- plot_pwy(feat.tab = feat.tab, G.pwy = G.pwy, stat.colnm = "EMY.z",
                     gr=gr, name = NA, colorbar.nm = "z", ntop = 7, seed = 1, plot = F, alternative="greater")
  # only returns genes in pwy
  expect_equal(sort(V(pp.ann)$name), c("abijabee3-_4", "abijabee3.4", "oh-not-so-Much"))
})

test_that("non-NA annot does not override feature name", {
  feat.tab["a", "symbol"] <- "abijabee3.4"
  pp.ann2 <- plot_pwy(feat.tab = feat.tab, G.pwy = G.pwy, stat.colnm = "EMY.z",
                      gr=gr, name = NA, colorbar.nm = "z", ntop = 7, seed = 1, plot = FALSE, alternative="greater")
  expect_equal(names(pp.ann2[[1]]), "a")
})

test_that("analyte in G & not feat.tab & connected to top.nodes is in plot", {
  pp.na <- plot_pwy(feat.tab = feat.tab[-1,], G.pwy = G.pwy, stat.colnm = "EMY.z",
                    gr=gr, name = NA, colorbar.nm = "z", ntop = 3, seed = 1, plot = T, alternative="greater")
  expect_true("a" %in% names(pp.na[[1]]))
})
