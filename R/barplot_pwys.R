#' Make barplots for the p-Values of Pathway Analysis
#'
#' Make barplots for the p-Values of pathway analysis using \pkg{ggplot2}.
#' @param ntop Number of top pathways.
#' @param direction Direction of change.
#' @param pwys_nm_size The maximum number of characters for pathway names. Longer names will be truncated.
#' @inheritParams ezheat
#' @inheritParams ezvenn
#' @return NULL
#' @export
#'

barplot_pwys <- function(tab, prefix.v=NULL, name = NA, width = 10, height = 4, ntop = 20, direction = c("Up", "Down"), pwys_nm_size = 100){

  if(!is.na(name)) {
    name <- paste0(name, "_barplots.pdf")
    pdf(name, width, height)
    on.exit(dev.off())
  }

  if(is.null(prefix.v)){
    p.cols <- grep("\\.p$", colnames(tab), value = TRUE)
    p.cols <- p.cols[-grep("Mixed", p.cols)]
    prefix.v <- gsub("\\.p", "", p.cols)
  }

  for(prefix in prefix.v) {
    tab.sub <- tab[, c("NGenes", paste(prefix, c("Direction", "p"), sep = "."))]
    colnames(tab.sub) <- gsub(paste0(prefix, "."), "", colnames(tab.sub))
    tab.sub <- tab.sub[order(tab.sub$p), ]

    for(d in direction) {
      dat2p <- tab.sub[tab.sub$Direction == d, ]
      dat2p <- dat2p[1:min(ntop, nrow(dat2p)), ]
      dat2p$Pathway <- substr(rownames(dat2p), 1, pwys_nm_size)
      dat2p$Pathway <- gsub("_", " ", dat2p$Pathway)
      dat2p$Pathway <- factor(dat2p$Pathway, levels = rev(dat2p$Pathway))
      dat2p$NGenes <- paste0("(", dat2p$NGenes, ")")
      dat2p$neglog10p <- - log10(dat2p$p)

      ggp <- ggplot(dat2p, aes(Pathway, neglog10p)) + theme_bw() +  coord_flip()
      ggp <- ggp + labs(x = "", y = expression("-"*log[10]~p*"-"*value), title = paste0(prefix, ", ", d, "-regulated"))

      if (d == "Up") {
        ggp <- ggp +geom_bar(stat = "identity", width = 0.7, fill = "red")
        ggp <- ggp + geom_text(aes(label = NGenes), hjust = -0.1) + ylim(0, 1.05*dat2p$neglog10p[1])
      } else {
        ggp <- ggp +geom_bar(stat = "identity", width = 0.7, fill = "blue") + scale_y_reverse()
        # msg: Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.
        ggp <- ggp + geom_text(aes(label = NGenes), hjust = 1.1) + ylim(1.05*dat2p$neglog10p[1], 0)
      }
      plot(ggp)
    }
  }
}
