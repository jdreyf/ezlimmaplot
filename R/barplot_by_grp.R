#' Make bar plots for a feature per group
#'
#' Make bar plots for a feature per group using \pkg{ggplot2}.
#' @param bar.width Width of the bars
#' @param errorbar.width Width of the errorbars
#' @inheritParams plot_by_grp
#' @inheritParams ezheat
#' @return NULL
#' @export

barplot_by_grp <- function(object, grp, name="top_genes", width=7, height=7,  main.v="", xlab="Group",  ylab="Log2 Expression",
                           manual.color=NULL, x.angle=0, bar.width=0.7, errorbar.width=0.3){

  if (is.vector(object)){ object <- t(as.matrix(object)) }
  if (all(main.v=="") & !is.null(rownames(object))){ main.v <- rownames(object) }
  grp <- factor(grp)

  stopifnot(ncol(object)==length(grp), nrow(object)==length(main.v), colnames(object)==names(grp))

  if(!is.na(name)) {
    pdf(paste0(name, "_barplots.pdf"), width=width, height=height)
    on.exit(dev.off())
  }

  for (i in 1:nrow(object)){
    dat2p <- data.frame(Exprs=object[i, ], Group=grp)
    dat2p <- dat2p[complete.cases(dat2p), ]
    dat2p <- plyr::ddply(dat2p, .(Group), summarise, N=length(Exprs), Mean=mean(Exprs), SE=sd(Exprs)/sqrt(N))

    ggp <- ggplot(data=dat2p, mapping=aes(x=Group, y=Mean, fill=Group)) + theme_bw()
    ggp <- ggp + labs(title=main.v[i], x=xlab, y=ylab) + theme(legend.position="none")
    ggp <- ggp + geom_bar(position=position_dodge(), stat="identity", width=bar.width, color="black")
    ggp <- ggp + geom_errorbar(mapping=aes(ymin=Mean, ymax=Mean+SE), position=position_dodge(), width=errorbar.width)

    if (!is.null( manual.color)) { ggp <- ggp + scale_fill_manual(values=manual.color)}
    if (x.angle !=0){ ggp <- ggp + theme(axis.text.x=element_text(angle=x.angle, hjust=1)) }
    plot(ggp)
  }
}
