#' Plot mediation exposure, mediator, and outcome per mediator
#'
#' Plot mediation exposure, mediator, and outcome per mediator in \code{M}, as used in \code{\link[Hitman]{lotman}} &
#' \code{\link[Hitman]{hitman}}, using \pkg{ggplot2}.
#'
#' @inheritParams ezheat
#' @inheritParams plot_by_grp
#' @inheritParams Hitman::lotman
#' @inheritParams Hitman::hitman
#' @return Invisibly, a \pkg{ggplot2} object from the last row of \code{M} that was plotted.
#' @export

plot_man <- function(E, M, Y, name="top_mediators", main.v="", xlab = "Log2 abundance",  ylab="Outcome",
                     check.names=TRUE, manual.color = NULL){
  stopifnot(limma::isNumeric(M), is.numeric(Y), !is.na(E), !is.na(Y), is.null(dim(E)), is.null(dim(Y)),
            stats::var(Y) > 0, nrow(M) >= 1, length(E)==ncol(M), length(Y)==ncol(M))
  if (is.data.frame(M)){ M <- data.matrix(M) }
  if (all(main.v == "") & !is.null(rownames(M))){ main.v <- rownames(M) }
  stopifnot(length(main.v)==nrow(M))
  if (check.names){
    stopifnot(names(E)==colnames(M), colnames(M)==names(Y))
  }

  if (!is.na(name)){ grDevices::pdf(paste0(name, ".pdf")) }
  for (ind in 1:nrow(M)){
    M2p <- data.frame(Exposure = E, Exprs = M[ind, ], Outcome=Y)
    ggp <- ggplot2::ggplot(data = M2p, mapping = ggplot2::aes(x = Exprs, y = Outcome, color=Exposure)) +
      ggplot2::labs(title = main.v[ind], x = xlab, y = ylab) + ggplot2::theme_bw() +
      ggplot2::geom_point() + ggplot2::geom_smooth(method="lm", se=FALSE)
    if(!is.null(manual.color)){ ggp <- ggp + ggplot2::scale_fill_manual(values = manual.color) +
      ggplot2::scale_color_manual(values = manual.color) }
    graphics::plot(ggp)
  }
  if (!is.na(name)){ grDevices::dev.off() }
  return(invisible(ggp))
}
