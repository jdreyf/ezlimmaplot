context("ezheat")

test_that("mat", {
  mat0 <- ezheat(M[topgenes,], labrows = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA, clip=1,
                 plot=FALSE)$mat
  expect_gt(max(mat0[1,1:2]), 0.75)
  expect_lt(min(mat0[1,4:6]), -0.5)
  gene.d <- rownames(res.df)[4]
  expect_gt(M[gene.d, "sample5"], max(M[gene.d, setdiff(colnames(M), "sample5")]))
  expect_equal(rownames(mat0)[11],  topgenes[11])
  expect_lte(max(abs(mat0)), 1)

  expect_message(mat1 <- ezheat(M[topgenes,], labrows = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA,
                                clip=1, only.labrows = TRUE, plot=FALSE, verbose=TRUE)$mat)
  expect_equal(nrow(mat1), sum(!is.na(res.df$Gene.Symbol[1:20])))

  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno2, name=NA, reorder_rows = TRUE, reorder_cols = TRUE, plot=FALSE)
  ret.mat=ret[[1]] # since ret is a list
  expect_false(all(colnames(ret.mat) == colnames(M)))
  expect_false(all(rownames(ret.mat) == rownames(M)))
})

test_that("sc", {
  ezh.sc <- ezheat(M, pheno.df=pheno.df, name=NA, sc="z", plot=FALSE)$mat
  #same as pheatmap:::scale_mat(M, scale = "row")
  base.sc <- t(scale(t(M), center=TRUE, scale = TRUE))
  expect_equal(ezh.sc, base.sc)
})

test_that("non vdiffr", {
  #pheatmap creates two pages, but vdiffr supports only one, so using returned gtable
  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno.df, name=NA)
  expect_equal(rownames(ret$mat), rownames(as.data.frame(lr)))
  expect_equal(colnames(ret$mat), row.names(pheno.df))
  expect_equal(rownames(ret$mat), topgenes)
  expect_equal(ret$mat, t(scale(x=t(M[topgenes,]), scale=FALSE)))
  expect_equal(ret$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression"))

  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno2, name=NA, sc="z")
  expect_equal(ret$mat, t(scale(x=t(M[topgenes,]), center=TRUE, scale=TRUE)))
  expect_equal(ret$gtable$grobs[[1]]$label, paste("Z-scored","Log2 Expression"))

  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  ret <- ezheat(object=res.df[topgenes, 1:2]+10, labrows = lr, name=NA, sc="none",
                  stat.tab=res.df[topgenes, c("First3.p", "Last3.p")])
  mat <- res.df[topgenes, 1:2]+10
  expect_equal(as.data.frame(ret$mat), as.data.frame(mat) )
  expect_equal(ret$gtable$grobs[[1]]$label, "Log2 Expression")

  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno2, name=NA,
                  reorder_rows = TRUE, reorder_cols = TRUE)
  expect_setequal(rownames(ret$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret$mat), row.names(pheno.df))
  expect_setequal(rownames(ret$mat), topgenes)
  expect_setequal(ret$mat, t(scale(x=t(M[topgenes,]), scale=FALSE)))
  expect_equal(ret$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression"))

  #verify by eye that columns & rows are reordered; title is only "Log2 Expression"
  lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
  ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno.df, name=NA, sc="none",
                  reorder_rows = TRUE, reorder_cols = TRUE, labcols=letters[1:6])
  expect_setequal(rownames(ret$mat), rownames(as.data.frame(lr)))
  expect_setequal(colnames(ret$mat), row.names(pheno.df))
  expect_setequal(rownames(ret$mat), topgenes)
  expect_setequal(ret$mat, t(M[topgenes,]))
  expect_equal(ret$gtable$grobs[[1]]$label, "Log2 Expression")

  #verify only unique rows
  ret <- ezheat(rbind(M[topgenes,], M[topgenes,]), name=NA, unique.rows=TRUE, sc="none")$mat
  expect_equal(ret, unique(rbind(M[topgenes,], M[topgenes,])))

  #verify only one row is included
  ret <- ezheat(M[topgenes[1],, drop=FALSE], name=NA, sc="none")$mat
  expect_equal(ret, head(M[topgenes[1],, drop=FALSE]))

  #verify only top row is included
  ret <- ezheat(M[topgenes,], name=NA, sc="none", ntop=1)$mat
  expect_equal(ret, head(M[topgenes,, drop=FALSE], n=1))

  #verify only top 10 rows are included
  ret <- ezheat(M[topgenes,], name=NA, sc="none", ntop=10)$mat
  expect_equal(ret, head(M[topgenes,], n=10))

  #verify only top 10 rows are included after taking unique set of rows
  ret <- ezheat(rbind(M[topgenes,], M[topgenes,]), name=NA, unique.rows=TRUE, sc="none", ntop=10)$mat
  expect_equal(ret, head(unique(rbind(M[topgenes,], M[topgenes,])), n=10))

  #verify rows are reordering
  ret <- ezheat(M[topgenes,], name=NA, sc="none", reorder_rows=TRUE)$mat
  expect_setequal(row.names(ret), row.names(M[topgenes,]))

  #verify columns are reordering
  ret <- ezheat(M[topgenes,], name=NA, sc="none", reorder_cols=FALSE, reorder_rows=TRUE)$mat
  expect_setequal(colnames(ret), colnames(M[topgenes,]))

  #verify rows and columns are reordering
  ret <- ezheat(M[topgenes,], name=NA, sc="none", reorder_cols=FALSE, reorder_rows=TRUE)$mat
  expect_setequal(row.names(ret), row.names(M[topgenes,]))
  expect_setequal(colnames(ret), colnames(M[topgenes,]))

  #verify any value >=1 in M is clipped to 1
  ret <- ezheat(M[topgenes,], name=NA, sc="none", clip=1)$mat
  expect_equal(max(ret), 1)
 })

test_that("non vdiffr low res", {
  # low res heat w/ clusters & w/o rowlabs
  cl.df <- data.frame(cl=LETTERS[as.numeric(sub("gene", "", rownames(M))) %% 5 + 1])
  rownames(cl.df) <- rownames(M)
  ret <- ezheat(M, pheno.df=pheno.df, labrows = "", name=NA, annotation_row = cl.df)
  expect_equal(ret$mat, t(scale(x=t(M), scale=FALSE)))
  expect_equal(colnames(ret$mat), row.names(pheno.df))
  expect_equal(ret$gtable$grobs[[1]]$label, paste("Centered","Log2 Expression"))
})

test_that("vdiffr", {
  #pheatmap creates two pages, but vdiffr supports only one, so using returned gtable
  ezh <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno.df, name=NA)$gtable
  }
  vdiffr::expect_doppelganger(title="heat", fig=ezh)

  ezh2 <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno2, name=NA, sc="z")$gtable
  }
  vdiffr::expect_doppelganger(title="heat2", fig=ezh2)

  #all avg+10 > 0
  ezh3 <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(object=res.df[topgenes, 1:2]+10, labrows = lr, name=NA, sc="none",
                  stat.tab=res.df[topgenes, c("First3.p", "Last3.p")])$gtable
  }
  vdiffr::expect_doppelganger(title="heat3", fig=ezh3)

  ezh4 <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno2, name=NA,
                  reorder_rows = TRUE, reorder_cols = TRUE)$gtable
  }
  #verify by eye that columns & rows are reordered; title is only "Log2 Expression"
  vdiffr::expect_doppelganger(title="heat4", fig=ezh4)

  ezh5 <- function(){
    lr <- sub_labels(topgenes, res.df[topgenes, "Gene.Symbol"])
    ret <- ezheat(M[topgenes,], labrows = lr, pheno.df=pheno.df, name=NA,
                  reorder_rows = TRUE, reorder_cols = TRUE, labcols=letters[1:6])$gtable
  }
  #verify by eye that columns & rows are reordered; title is only "Log2 Expression"
  vdiffr::expect_doppelganger(title="heat5", fig=ezh5)
})

test_that("vdiffr low res", {
  # low res heat w/ clusters & w/o rowlabs
  cl.df <- data.frame(cl=LETTERS[as.numeric(sub("gene", "", rownames(M))) %% 5 + 1])
  rownames(cl.df) <- rownames(M)

  ezh6 <- function(){
    ret <- ezheat(M, pheno.df=pheno.df, labrows = "", name=NA, annotation_row = cl.df)
  }
  #verify by eye that columns & rows are reordered; title is only "Log2 Expression"
  vdiffr::expect_doppelganger(title="heat6", fig=ezh6)
})
