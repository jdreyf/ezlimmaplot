#' Make Venn diagrams
#' 
#' Make Venn diagrams using output from \code{ezlimma} package.

ezvenn <- function(dat, contrast.v, fdr.cutoff = NULL, p.cutoff = NULL, logfc.cutoff = NULL,
                    circle.name = NULL, main = '', name = NA, cex = c(1, 1, 1), plot = TRUE) {
  
  contrasts <- names(contrast.v)
  if(!is.null(fdr.cutoff)){
    fdr.col <- paste0(contrasts, '.FDR')
    stopifnot(fdr.col %in% colnames(dat))
    dat.sig <- dat[, fdr.col]
    dat.sig <- dat.sig < fdr.cutoff
  } else if (!is.null(p.cutoff)){
    p.col <- paste0(contrasts, '.p')
    stopifnot(p.col %in% colnames(dat))
    dat.sig <- dat[, p.col]
    dat.sig <- dat.sig < p.cutoff
  } else {
    cat('Please set FDR/p-value cutoff')
    stop
  }
  
  if(is.null(circle.name)) circle.name <- contrasts
  
  if(!is.na(name) & plot) {
    name <- paste0(name, "_venn.pdf")
    pdf(name)
  }
  
  logfc.col <- paste0(contrasts, '.logFC')
  if (all(logfc.col %in% colnames(dat))){
    dat.logfc <- dat[, logfc.col]
    if(!is.null(logfc.cutoff)){
      dat.sig <- dat.sig * (abs(dat.logfc) > logfc.cutoff)
    }
    
    dat.sig <- dat.sig * sign(dat.logfc)
    dat.sig[is.na(dat.sig)] <- 0
    
    if (plot) limma::vennDiagram(dat.sig, include = c('up', 'down'), names = circle.name,
                         circle.col = rainbow(length(contrasts)), counts.col = c('red', 'blue'),
                         main = main, cex = cex)
  } else {
    dat.sig[is.na(dat.sig)] <- 0
    if (plot) limma::vennDiagram(dat.sig, include = 'both', names = circle.name,
                         circle.col = rainbow(length(contrasts)), counts.col = 'blue',
                         main = main, cex = cex)
  }
  
  if (!is.na(name) & plot) dev.off()
  
  colnames(dat.sig) <- gsub('\\.(logFC|p|FDR)$', '.sig', colnames(dat.sig))
  dat.sig <- dat.sig[order(-abs(rowSums(dat.sig))), ]
  return(dat.sig)
}