% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_missing.R
\name{plot_missing}
\alias{plot_missing}
\title{Plot missingness probability against average log abundance}
\usage{
plot_missing(object, name = NA, alpha = 0.5, link = c("logit", "probit"))
}
\arguments{
\item{object}{Matrix-like object with features (e.g. genes) as rows and samples as columns.}

\item{name}{File name of the plot.}

\item{alpha}{Transparency of color for data points. The regression line uses half this transparency.}

\item{link}{Binomial link to use. Possibilities are logit or probit.}
}
\value{
Invisibly, a ggplot object.
}
\description{
Plot missingness probability against average log abundance.
}
