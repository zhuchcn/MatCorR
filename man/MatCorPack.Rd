% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MatCor.R
\name{MatCorPack}
\alias{MatCorPack}
\title{Compute correlation in mutiple methods}
\usage{
MatCorPack(X, Y, methods = c("pearson", "spearman", "kendall", "lm"),
  adjust.method = "BH", design)
}
\arguments{
\item{X}{numaric matrix of data values.}

\item{Y}{numeric matrix of data values. Must have same number of columns as X.}

\item{adjust.method}{character value indicating the method to use for
p value adjustement. See \code{\link{p.adjust}}.}

\item{design}{a design matrix \strong{must} if "lm" is in the methods. The
design matrix can be built using the \code{\link{model.matrix}} function.}

\item{method}{character value indicating the method to use. Must be one of
"pearson", "spearman", "kendall", or "lm".}
}
\value{
A nested list, with the methods being the first level, X variables
being the second, and Y variables the bottom level.
}
\description{
Compute correlation tests between two matrices using multiple
methods and return a nested list. The length of the first level list equals
to the numebr of methods input. Each member of the first level list is a
output of the \code{\link{MatCor}} function. The methods supported are
pearson, spearman, kendall, and lm.
}
\author{
Chenghao Zhu
}
