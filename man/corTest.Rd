% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MatCor.R
\name{corTest}
\alias{corTest}
\title{Calculate correlations between two vectors}
\usage{
corTest(x, y, method, ...)
}
\arguments{
\item{x}{numeric vectors of data values.}

\item{y}{numeric vectors of data values. Must have the same length as x.}

\item{method}{character value indicating which correlation method to be used
One of "peason", "kendall", "spearman", or "lm". If "lm" is used, a design
argument must be specified.}

\item{...}{other arguments. If either of "person", "kendal", or "spearman"
is specified, any other arguments in the \code{\link{cor.test}} function can
be used. If "lm" is specified, a design matrix \strong{must} be passed. The
design matrix can be built using the \code{\link{model.matrix}} function.}
}
\value{
A named vector with length of 3 is returned

\item{\strong{stat}}{the value of the test statistic.}
\item{\strong{estimate}}{the estimated measure of association. Corresponds to the "r"
for pearson, "rho" for spearman, "tau" for kendal, and coefficient for lm}
\item{\strong{pval}}{the p-value of the test.}
}
\description{
This function takes 2 numeric vectors and calculate the
correlation between them. The \code{\link{cor.test}} function from the
\code{\link{stats-package}} is used to calculate pearson, spearman, or
kendal correlation. A more complex model is supported by setting the method
argument to "lm", and assign a model matrix to the design argument.
}
\author{
Chenghao Zhu
}
