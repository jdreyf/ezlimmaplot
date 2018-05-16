#' Mean and standard error
#'
#' Gives mean, mean - standard error, & mean + standard error
#'
#' @param x A numeric vector
#' @import stats

mean_se <- function(x){
  x <- x[!is.na(x)]
  m <- mean(x)
  se <- stats::sd(x)/sqrt(length(x))
  return(c(y = m, ymin = m - se, ymax = m + se))
}
