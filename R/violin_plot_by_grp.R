#' Make violin plots for a feature per group
#'
#' Make violin plots for a feature per group using \pkg{ggplot2}.
#' @param violin.alpha Transparency of violins
#' @param dot.alpha Transparency of dots
#' @inheritParams plot_by_grp
#' @inheritParams ezheat
#' @return NULL
#' @export

violin_plot_by_grp <- function(object, grp, name="top_genes", width=7, height=7,  main.v="", xlab="Group",  ylab="Log2 Expression",
                               manual.color=NULL, x.angle=0, violin.alpha=0.3, dot.alpha=0.3){

  if (is.vector(object)){ object <- t(as.matrix(object)) }
  if (all(main.v=="") & !is.null(rownames(object))){ main.v <- rownames(object) }
  grp <- factor(grp)

  stopifnot(ncol(object)==length(grp), nrow(object)==length(main.v), colnames(object)==names(grp))

  if(!is.na(name)){
    pdf(paste0(name, "_violin_plots.pdf"), width=width, height=height)
    on.exit(dev.off())
  }

  for (i in 1:nrow(object)){
    dat2p <- data.frame(Exprs=object[i, ], Group=grp)
    dat2p <- dat2p[complete.cases(dat2p), ]

    ggp <- ggplot(data=dat2p, mapping=aes(x=Group, y=Exprs, color=Group, fill=Group)) + theme_bw()
    ggp <- ggp + labs(title=main.v[i], x=xlab, y=ylab) + theme(legend.position="none")
    ggp <- ggp + geom_violin(trim=FALSE, alpha=violin.alpha)
    ggp <- ggp + geom_jitter(shape=21, position=position_jitter(0.2), alpha=dot.alpha)

    if (!is.null(manual.color)) { ggp <- ggp + scale_color_manual(values=manual.color) + scale_fill_manual(values=manual.color)}
    if (x.angle !=0){ ggp <- ggp + theme(axis.text.x=element_text(angle=x.angle, hjust=1)) }
    plot(ggp)
  }
}
