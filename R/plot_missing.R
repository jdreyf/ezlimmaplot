#' Plot missingness probability against average log abundance
#'
#' Plot missingness probability against average log abundance.
#'
#' @param name File name of the plot.
#' @param alpha Transparency of color for data points. The regression line uses half this transparency.
#' @param link Binomial link to use. Possibilities are logit or probit.
#' @inheritParams ezheat
#' @return Invisibly, a ggplot object.
#' @export

plot_missing <- function(object, name=NA, alpha=0.5, link=c("logit", "probit")){
  link <- match.arg(link)
  stopifnot(sum(is.na(object)) > 0, rowMeans(is.na(object)) < 1)

  na.prop <- rowSums(is.na(object))/ncol(object)
  # na.prop <- na.prop[na.prop != 1]
  # na.logit <- log2(na.prop/(1 - na.prop))
  mu <- rowMeans(object, na.rm=TRUE)
  fm <- suppressWarnings( stats::glm(na.prop ~ mu, family=stats::binomial(link = link)) )
  coeff <- stats::coefficients(fm)
  pv.mu <- summary(fm)$coefficients["mu", "Pr(>|z|)"]
  pv.ttle <- ifelse(pv.mu < 10**-15, "p < 10^(-15)", paste("p =", signif(pv.mu, 3)))

  xx <- seq(from=min(mu), to=max(mu), length.out=10**3)
  y.exp <- exp(coeff[1] + coeff[2]*xx)/(1 + exp(coeff[1] + coeff[2]*xx))

  gg.df1 <- data.frame(mu, na.prop)
  gg.df2 <- data.frame(xx, y.exp)

  ggp <- ggplot2::ggplot() + ggplot2::geom_point(mapping = ggplot2::aes(mu, na.prop), data = gg.df1, alpha=alpha)+
    ggplot2::geom_point(mapping = ggplot2::aes(xx, y.exp), data=gg.df2, color="red", alpha=alpha/5)+
    ggplot2::labs(x="Average abundance of observed psite values", y="Proportion of samples psite is missing",
         title = paste("Abundance-dependent missingness", pv.ttle),
         subtitle = "Logistic fit in red")

  if (!is.na(name)){
    grDevices::pdf(paste0(name, ".pdf"))
    on.exit(grDevices::dev.off())
  }
  base::plot(ggp)
  return(invisible(ggp))
}
