context("ezheat")

sym <- c("A", "B", "C", "D", "E")
pheno.df <- data.frame(row.names = rownames(pheno), grp=pheno$grp)
topgenes <- rownames(res.df)[1:5]
ezh <- ezheat(M[topgenes,], symbols = res.df[topgenes, "Gene.Symbol"], pheno.df=pheno.df, name=NA)
