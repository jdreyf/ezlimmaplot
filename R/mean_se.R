#' Mean and standard error
#'
#' Gives mean, mean - standard error, & mean + standard error
#'
#' @param x A numeric vector

mean_se <- function(x){
  stopifnot(sum(!is.na(x)) > 0, is.numeric(x), is.vector(x), stats::sd(x, na.rm = TRUE) > 0)
  x <- x[!is.na(x)]
  m <- mean(x)
  se <- stats::sd(x)/sqrt(length(x))
  return(c(y = m, ymin = m - se, ymax = m + se))
}
