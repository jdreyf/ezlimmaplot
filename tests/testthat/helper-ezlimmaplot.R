Sys.setenv(VDIFFR_RUN_TESTS=FALSE)

library(covr)
library(ezlimma)
library(ggplot2)
library(igraph)
library(testthat)
library(vdiffr)

set.seed(42)
M <- matrix(rnorm(100*6, sd=0.3), nrow=100, ncol=6)
dimnames(M) <- list(paste0("gene", 1:nrow(M)), paste0("sample", 1:ncol(M)))
grp <- rep(c("First3", "Last3"), each=3)
tissue <- rep(c("muscle", "liver"), times=3)
covar_num <- rnorm(n = ncol(M))
pheno <- data.frame(row.names = colnames(M), sample=colnames(M), grp, tissue, covar_num, stringsAsFactors = FALSE)
M[1,1:3] <- M[1,1:3] + 2

contr.v <- c(First3="First3", Last3="Last3", Last3vsFirst3="Last3-First3")
res <- ezlimma::limma_contrasts(M, grp=grp, contrast.v = contr.v)
res.df <- data.frame(signif(res, 3), Gene.Symbol=NA, stringsAsFactors = FALSE)
res.df$Gene.Symbol[1:10] <- LETTERS[1:10]

sym.v <- rep(NA, nrow(M))
sym.v[1:9] <- letters[1:9]
sym.v[10:13] <- c("---", "a", "~", "")

pheno.df <- pheno[,"grp",drop=F]
pheno2 <- pheno[,c("grp", "tissue"),drop=F]
topgenes <- rownames(res.df)[1:20]

res.ss <- res.df[topgenes, grep("avg$|^Last3vsFirst|^Gene", colnames(res.df))]

## plot pwy
#sq w/ all nodes connected except a <-> c; a <-> d repeated, as in sif
el <- rbind(t(combn(letters[1:4], 2))[-2,], c("a", "d"))
gr <- igraph::graph_from_edgelist(el, directed = FALSE)
gr2 <- igraph::add_edges(gr, edges=c("a", "b"))
gr <- edgelist2graph(el)

set.seed(0)
object <- matrix(rnorm(n=90), ncol=9, dimnames=list(letters[1:10], paste0("s", 1:9)))
object["a", 1:3] <- object["a", 1:3]+5
E <- rep(1:0, times=c(3, 6))
Y <- object["a",]
pheno.num <- Y
pheno.chr <- rep(c("trt1", "trt2", "ctrl"), times=3)
hm <- Hitman::hitman(M=object, Y=Y, E=E)
hm$symbol <- toupper(rownames(hm))

# plot_pwy
gmt <- list(pwy1=list(name="pwy1", description="pwy1", genes=c("a", "b", "c")),
            pwy2=list(name="pwy2", description="pwy2", genes=c("b", "c", "d")))

feat.tab <- hm
G.pwy = gmt[[1]]
pp <- plot_pwy(feat.tab = feat.tab, G.pwy = gmt[[1]], stat.colnm = "EMY.chisq", annot.colnm = "symbol", repel=TRUE,
               gr=gr, name = NA, colorbar.nm = "chisq", ntop = 7, seed = 0, plot = FALSE, alternative="greater")$gg.pwy

# barplot pwy
G <- list(pwy1=list(name="pwy1", description="pwy1", genes=paste0("gene", 1:3)),
          pwy2=list(name="pwy2", description="pwy2", genes=paste0("gene", 2:4)))

G2 <- G
G2$pwy3 <- list(name="pwy3", description="pwy3", genes=paste0("gene", 5:9))
G2$mother_of_all_pathways <- list(name="mother_of_all_pathways", description="mother_of_all_pathways", genes=rownames(res)[1:3])

rc <- ezlimma::roast_contrasts(object=M, G=G, feat.tab=res, grp=grp, contrast.v = contr.v, fun="fry")
rc2 <- ezlimma::roast_contrasts(object=M, G=G2, feat.tab=res, grp=grp, contrast.v = contr.v, fun="fry")
