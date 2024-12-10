context("multi_heat")

test_that("mat", {
  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  # warning prune_mat n rows < 50
  expect_warning(mat0 <- multi_heat(tab=res.ss, object=M[topgenes,], labrows = lr, pheno.df=pheno.df, name="tmp", clip=1, plot=FALSE)[[1]]$mat)
  gene.d <- rownames(res.df)[4]
  expect_gt(M[gene.d, "sample5"], max(M[gene.d, setdiff(colnames(M), "sample5")]))
  expect_lte(max(abs(mat0)), 1)

  # error due to labrows not having names
  expect_error(mat1 <- multi_heat(tab=res.ss, object=M[topgenes,], labrows = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name="tmp", clip=1, only.labrows = TRUE,
                                    plot=FALSE, verbose=TRUE, ntop=10)[[1]]$mat)

  # msg is about removing rows with is.na(labrow)
  lr.tmp <- setNames(res.df[topgenes, "Gene.Symbol"], nm=topgenes)
  expect_message(mat1 <- multi_heat(tab=res.ss, object=M[topgenes,], labrows = lr.tmp, pheno.df=pheno.df, name="tmp", clip=1, only.labrows = TRUE,
                                    plot=FALSE, verbose=TRUE, ntop=10)[[1]]$mat)
})

test_that("non vdiffr single heatmap", {
  #pheatmap creates two pages, but vdiffr supports only one, so using returned gtable
  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  p.cols <- grep("\\.p$", colnames(res.ss), value=TRUE)
  contr.names <- sub("\\.(p)$", "", p.cols)

  expect_warning(ret <- multi_heat(tab=res.ss, M[topgenes,], labrows = lr, pheno.df=pheno.df, name="tmp"))
  expect_equal(sort(rownames(ret[[1]]$mat)), sort(rownames(as.data.frame(lr))))
  expect_equal(sort(colnames(ret[[1]]$mat)), sort(row.names(pheno.df)))
  expect_equal(sort(rownames(ret[[1]]$mat)), sort(topgenes))
  expect_equal(ret[[1]]$mat["gene1",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene1",])
  expect_equal(ret[[1]]$mat["gene60",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene60",])
  expect_equal(ret[[1]]$mat["gene74",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene74",])
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression", contr.names))

  expect_warning(ret <- multi_heat(tab=res.ss, M[topgenes,], labrows = lr, pheno.df=pheno2, name="tmp", sc="z"))
  expect_equal(ret[[1]]$mat["gene1",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene1",])
  expect_equal(ret[[1]]$mat["gene60",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene60",])
  expect_equal(ret[[1]]$mat["gene74",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene74",])
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Z-scored","Log2 Expression", contr.names))

  expect_warning(ret <- multi_heat(tab=res.ss, object=res.df[topgenes, 1:2]+10, labrows = lr, name="tmp", sc="none",
                  stat.tab=res.df[topgenes, c("First3.p", "Last3.p")]))
  mat <- res.df[topgenes, 1:2]+10
  expect_equal(as.data.frame(ret[[1]]$mat)["gene1",], as.data.frame(mat)["gene1",] )
  expect_equal(as.data.frame(ret[[1]]$mat)["gene74",], as.data.frame(mat)["gene74",] )
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Log2 Expression", contr.names))

  expect_warning(ret <- multi_heat(tab=res.ss, M[topgenes,], labrows = lr, pheno.df=pheno2, name="tmp",
                  reorder_rows = TRUE, reorder_cols = TRUE))
  expect_setequal(rownames(ret[[1]]$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret[[1]]$mat), row.names(pheno.df))
  expect_setequal(rownames(ret[[1]]$mat), topgenes)
  expect_setequal(ret[[1]]$mat, t(scale(x=t(M[topgenes,]), scale=FALSE)))
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression", contr.names))

  expect_warning(ret <- multi_heat(tab=res.ss, M[topgenes,], labrows = lr, pheno.df=pheno.df, name="tmp", sc="none",
                  reorder_rows = TRUE, reorder_cols = TRUE, labcols=letters[1:6]))
  expect_setequal(rownames(ret[[1]]$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret[[1]]$mat), row.names(pheno.df))
  expect_setequal(rownames(ret[[1]]$mat), topgenes)
  expect_setequal(ret[[1]]$mat, t(M[topgenes,]))
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Log2 Expression", contr.names))

  #verify only unique rows
  expect_warning(ret <- multi_heat(tab=res.ss, rbind(M[topgenes,], M[topgenes,]), name="tmp", unique.rows=TRUE, sc="none")[[1]]$mat)
  expect_setequal(ret, unique(rbind(M[topgenes,], M[topgenes,])))

  #verify rows are reordering
  expect_warning(ret <- multi_heat(tab=res.ss, M[topgenes,], name="tmp", sc="none", reorder_rows=TRUE)[[1]]$mat)
  expect_setequal(row.names(ret), row.names(M[topgenes,]))

  #verify columns are reordering
  expect_warning(ret <- multi_heat(tab=res.ss, M[topgenes,], name="tmp", sc="none", reorder_cols=FALSE, reorder_rows=TRUE)[[1]]$mat)
  expect_setequal(colnames(ret), colnames(M[topgenes,]))

  #verify rows and columns are reordering
  expect_warning(ret <- multi_heat(tab=res.ss, M[topgenes,], name="tmp", sc="none", reorder_cols=FALSE, reorder_rows=TRUE)[[1]]$mat)
  expect_setequal(row.names(ret), row.names(M[topgenes,]))
  expect_setequal(colnames(ret), colnames(M[topgenes,]))

  #verify any value >=1 in M is clipped to 1
  expect_warning(ret <- multi_heat(tab=res.ss, M[topgenes,], name="tmp", sc="none", clip=1)[[1]]$mat)
  expect_equal(max(ret), 1)
})

test_that("non vdiffr multiple heatmap", {
  #pheatmap creates two pages, but vdiffr supports only one, so using returned gtable
  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  p.cols <- grep(paste0("\\.p$"), colnames(res), value=TRUE)
  contr.names <- sub(paste0("\\.(p)$"), "", p.cols)

  expect_warning(ret <- multi_heat(tab=res, M[topgenes,], labrows = lr, pheno.df=pheno.df, name="tmp"))
  expect_equal(sort(rownames(ret[[1]]$mat)), sort(rownames(as.data.frame(lr))))
  expect_equal(sort(colnames(ret[[1]]$mat)), sort(row.names(pheno.df)))
  expect_equal(sort(rownames(ret[[1]]$mat)), sort(topgenes))
  expect_equal(ret[[1]]$mat["gene1",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene1",])
  expect_equal(ret[[1]]$mat["gene60",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene60",])
  expect_equal(ret[[1]]$mat["gene74",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene74",])
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression", contr.names[1]))
  expect_equal(sort(rownames(ret[[2]]$mat)), sort(rownames(as.data.frame(lr))))
  expect_equal(sort(colnames(ret[[2]]$mat)), sort(row.names(pheno.df)))
  expect_equal(sort(rownames(ret[[2]]$mat)), sort(topgenes))
  expect_equal(ret[[2]]$mat["gene1",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene1",])
  expect_equal(ret[[2]]$mat["gene60",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene60",])
  expect_equal(ret[[2]]$mat["gene74",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene74",])
  expect_equal(ret[[2]]$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression", contr.names[2]))
  expect_equal(sort(rownames(ret[[3]]$mat)), sort(rownames(as.data.frame(lr))))
  expect_equal(sort(colnames(ret[[3]]$mat)), sort(row.names(pheno.df)))
  expect_equal(sort(rownames(ret[[3]]$mat)), sort(topgenes))
  expect_equal(ret[[3]]$mat["gene1",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene1",])
  expect_equal(ret[[3]]$mat["gene60",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene60",])
  expect_equal(ret[[3]]$mat["gene74",], t(scale(x=t(M[topgenes,]), scale=FALSE))["gene74",])
  expect_equal(ret[[3]]$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression", contr.names[3]))

  expect_warning(ret <- multi_heat(tab=res, M[topgenes,], labrows = lr, pheno.df=pheno2, name="tmp", sc="z"))
  expect_equal(ret[[1]]$mat["gene1",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene1",])
  expect_equal(ret[[1]]$mat["gene60",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene60",])
  expect_equal(ret[[1]]$mat["gene74",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene74",])
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Z-scored","Log2 Expression", contr.names[1]))
  expect_equal(ret[[2]]$mat["gene1",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene1",])
  expect_equal(ret[[2]]$mat["gene60",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene60",])
  expect_equal(ret[[2]]$mat["gene74",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene74",])
  expect_equal(ret[[2]]$gtable$grobs[[1]]$label, paste("Z-scored","Log2 Expression", contr.names[2]))
  expect_equal(ret[[3]]$mat["gene1",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene1",])
  expect_equal(ret[[3]]$mat["gene60",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene60",])
  expect_equal(ret[[3]]$mat["gene74",], t(scale(x=t(M[topgenes,]), scale=TRUE))["gene74",])
  expect_equal(ret[[3]]$gtable$grobs[[1]]$label, paste("Z-scored","Log2 Expression", contr.names[3]))

  expect_warning(ret <- multi_heat(tab=res, object=res.df[topgenes, 1:2]+10, labrows = lr, name="tmp", sc="none",
                  stat.tab=res.df[topgenes, c("First3.p", "Last3.p")]))
  mat <- res.df[topgenes, 1:2]+10
  expect_equal(as.data.frame(ret[[1]]$mat)["gene1",], as.data.frame(mat)["gene1",] )
  expect_equal(as.data.frame(ret[[1]]$mat)["gene74",], as.data.frame(mat)["gene74",] )
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Log2 Expression", contr.names[1]))
  expect_equal(as.data.frame(ret[[2]]$mat)["gene1",], as.data.frame(mat)["gene1",] )
  expect_equal(as.data.frame(ret[[2]]$mat)["gene74",], as.data.frame(mat)["gene74",] )
  expect_equal(ret[[2]]$gtable$grobs[[1]]$label, paste("Log2 Expression", contr.names[2]))
  expect_equal(as.data.frame(ret[[3]]$mat)["gene1",], as.data.frame(mat)["gene1",] )
  expect_equal(as.data.frame(ret[[3]]$mat)["gene74",], as.data.frame(mat)["gene74",] )
  expect_equal(ret[[3]]$gtable$grobs[[1]]$label, paste("Log2 Expression", contr.names[3]))

  expect_warning(ret <- multi_heat(tab=res, M[topgenes,], labrows = lr, pheno.df=pheno2, name="tmp",
                  reorder_rows = TRUE, reorder_cols = TRUE))
  expect_setequal(rownames(ret[[3]]$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret[[3]]$mat), row.names(pheno.df))
  expect_setequal(rownames(ret[[3]]$mat), topgenes)
  expect_setequal(ret[[3]]$mat, t(scale(x=t(M[topgenes,]), scale=FALSE)))
  expect_equal(ret[[3]]$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression", contr.names[3]))
  expect_setequal(rownames(ret[[2]]$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret[[2]]$mat), row.names(pheno.df))
  expect_setequal(rownames(ret[[2]]$mat), topgenes)
  expect_setequal(ret[[2]]$mat, t(scale(x=t(M[topgenes,]), scale=FALSE)))
  expect_equal(ret[[2]]$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression", contr.names[2]))
  expect_setequal(rownames(ret[[1]]$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret[[1]]$mat), row.names(pheno.df))
  expect_setequal(rownames(ret[[1]]$mat), topgenes)
  expect_setequal(ret[[1]]$mat, t(scale(x=t(M[topgenes,]), scale=FALSE)))
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression", contr.names[1]))

  expect_warning(ret <- multi_heat(tab=res, M[topgenes,], labrows = lr, pheno.df=pheno.df, name="tmp", sc="none",
                  reorder_rows = TRUE, reorder_cols = TRUE, labcols=letters[1:6]))
  expect_setequal(rownames(ret[[1]]$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret[[1]]$mat), row.names(pheno.df))
  expect_setequal(rownames(ret[[1]]$mat), topgenes)
  expect_setequal(ret[[1]]$mat, t(M[topgenes,]))
  expect_equal(ret[[1]]$gtable$grobs[[1]]$label, paste("Log2 Expression", contr.names[1]))
  expect_setequal(rownames(ret[[2]]$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret[[2]]$mat), row.names(pheno.df))
  expect_setequal(rownames(ret[[2]]$mat), topgenes)
  expect_setequal(ret[[2]]$mat, t(M[topgenes,]))
  expect_equal(ret[[2]]$gtable$grobs[[1]]$label, paste("Log2 Expression", contr.names[2]))
  expect_setequal(rownames(ret[[3]]$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret[[3]]$mat), row.names(pheno.df))
  expect_setequal(rownames(ret[[3]]$mat), topgenes)
  expect_setequal(ret[[3]]$mat, t(M[topgenes,]))
  expect_equal(ret[[3]]$gtable$grobs[[1]]$label, paste("Log2 Expression", contr.names[3]))

  #verify only unique rows
  expect_warning(ret <- multi_heat(tab=res, rbind(M[topgenes,], M[topgenes,]), name="tmp", unique.rows=TRUE, sc="none"))
  expect_setequal(ret[[1]]$mat, unique(rbind(M[topgenes,], M[topgenes,])))
  expect_setequal(ret[[2]]$mat, unique(rbind(M[topgenes,], M[topgenes,])))
  expect_setequal(ret[[3]]$mat, unique(rbind(M[topgenes,], M[topgenes,])))

  #verify rows are reordering
  expect_warning(ret <- multi_heat(tab=res, M[topgenes,], name="tmp", sc="none", reorder_rows=TRUE))
  expect_setequal(row.names(ret[[1]]$mat), row.names(M[topgenes,]))
  expect_setequal(row.names(ret[[2]]$mat), row.names(M[topgenes,]))
  expect_setequal(row.names(ret[[3]]$mat), row.names(M[topgenes,]))

  #verify columns are reordering
  expect_warning(ret <- multi_heat(tab=res, M[topgenes,], name="tmp", sc="none", reorder_cols=FALSE, reorder_rows=TRUE))
  expect_setequal(colnames(ret[[1]]$mat), colnames(M[topgenes,]))
  expect_setequal(colnames(ret[[2]]$mat), colnames(M[topgenes,]))
  expect_setequal(colnames(ret[[3]]$mat), colnames(M[topgenes,]))

  #verify rows and columns are reordering
  expect_warning(ret <- multi_heat(tab=res, M[topgenes,], name="tmp", sc="none", reorder_cols=FALSE, reorder_rows=TRUE))
  expect_setequal(row.names(ret[[1]]$mat), row.names(M[topgenes,]))
  expect_setequal(colnames(ret[[1]]$mat), colnames(M[topgenes,]))
  expect_setequal(row.names(ret[[2]]$mat), row.names(M[topgenes,]))
  expect_setequal(colnames(ret[[2]]$mat), colnames(M[topgenes,]))
  expect_setequal(row.names(ret[[3]]$mat), row.names(M[topgenes,]))
  expect_setequal(colnames(ret[[3]]$mat), colnames(M[topgenes,]))

  #verify any value >=1 in M is clipped to 1
  expect_warning(ret <- multi_heat(tab=res, M[topgenes,], name="tmp", sc="none", clip=1))
  expect_equal(max(ret[[1]]$mat), 1)
  expect_equal(max(ret[[2]]$mat), 1)
  expect_equal(max(ret[[3]]$mat), 1)
})
