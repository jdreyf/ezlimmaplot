library(covr)
library(testthat)
library(vdiffr)

set.seed(42)
M <- matrix(rnorm(100*6, sd=0.3), nrow=100, ncol=6)
dimnames(M) <- list(paste0("gene", 1:nrow(M)), paste0("sample", 1:ncol(M)))
grp <- rep(c("First3", "Last3"), each=3)
tissue <- rep(c("muscle", "liver"), times=3)
pheno <- data.frame(row.names = colnames(M), sample=colnames(M), grp=grp, tissue=tissue, stringsAsFactors = FALSE)
M[1,1:3] <- M[1,1:3] + 2

contr.v <- c(First3="First3", Last3="Last3", Last3vsFirst3="Last3-First3")
res <- ezlimma::limma_contrasts(M, grp=grp, contrast.v = contr.v)
res.df <- data.frame(signif(res, 3), Gene.Symbol=NA)
res.df$Gene.Symbol[1:10] <- LETTERS[1:10]

sym.v <- rep(NA, nrow(M))
sym.v[1:9] <- letters[1:9]
sym.v[10:13] <- c("---", "a", "~", "")

pheno.df <- pheno[,"grp",drop=F]
pheno2 <- pheno[,c("grp", "tissue"),drop=F]
topgenes <- rownames(res.df)[1:20]
