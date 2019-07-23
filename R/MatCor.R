################################################################################
#' @title Calculate correlations between two vectors
#'
#' @description This function takes 2 numeric vectors and calculate the
#' correlation between them. The \code{\link{cor.test}} function from the
#' \code{\link{stats-package}} is used to calculate pearson, spearman, or
#' kendal correlation. A more complex model is supported by setting the method
#' argument to "lm", and assign a model matrix to the design argument.
#'
#' @param x numeric vectors of data values.
#' @param y numeric vectors of data values. Must have the same length as x.
#' @param method character value indicating which correlation method to be used
#' One of "peason", "kendall", "spearman", or "lm". If "lm" is used, a design
#' argument must be specified.
#' @param ... other arguments. If either of "person", "kendal", or "spearman"
#' is specified, any other arguments in the \code{\link{cor.test}} function can
#' be used. If "lm" is specified, a design matrix \strong{must} be passed. The
#' design matrix can be built using the \code{\link{model.matrix}} function.
#'
#'
#' @return A named vector with length of 3 is returned
#'
#' \item{\strong{stat}}{the value of the test statistic.}
#' \item{\strong{estimate}}{the estimated measure of association. Corresponds to the "r"
#' for pearson, "rho" for spearman, "tau" for kendal, and coefficient for lm}
#' \item{\strong{pval}}{the p-value of the test.}
#'
#' @export
#' @import dplyr
#' @author Chenghao Zhu
corTest = function(x, y, method, ...){
    if(method %in% c("pearson", "spearman", "kendall")){
        cor = tryCatch(
            cor.test(x=x, y=y, method = method, ...),
            error = function(e){
                return(list(statistic = 0, estimate = 0, p.value = 1))
            }
        )
        result = c(cor$statistic, cor$estimate, cor$p.value)
        names(result) = c("stat", "estimate", "pval")
    }else if(method == "lm"){
        args = list(...)
        if(!"design" %in% names(args))
            stop("[ MatCorR::corTest ] the design is missing", call. = FALSE)
        fit = tryCatch(
            summary(lm(y ~ x + design + 0))$coefficients,
            error = function(e){
                out = matrix(c(0, 0, 1), nrow = 1)
                rownames(out) = "x"
                colnames(out) = c("t value", "Estimate", "Pr(>|t|)")
                return(out)
            }
        )
        if("x" %in% rownames(fit)){
            result = c(fit["x", c("t value", "Estimate", "Pr(>|t|)")])
        }else{
            result = c(0, 0, 1)
        }
        names(result) = c("stat", "estimate", "pval")
        if(is.nan(result["stat"])) result["stat"] = 0
        if(is.nan(result["estimate"])) result["estimate"] = 0
        if(is.nan(result["pval"])) result["pval"] = 1
    }else stop("[ MatCorR::corTest ] Method not found", call. = FALSE)
    return(result)
}

################################################################################
#' @title Calculate correlation between two matrices
#'
#' @description This funciton calculate correlation between 2 matrices. The 2
#' matrices must have samples (observations) on columns, and variables to
#' compute correltaion on rows. The 2 matrices must have same number of
#' columns. Methods supported are pearson, spearman, kendall, and linear model
#' with more complex design.
#'
#' @param X numaric matrix of data values.
#' @param Y numeric matrix of data values. Must have same number of columns as X.
#' @param method character value indicating the method to use. Must be one of
#' "pearson", "spearman", "kendall", or "lm".
#' @param adjust.method character value indicating the method to use for
#' p value adjustement. See \code{\link{p.adjust}}.
#' @param ... other arguments. If either of "person", "kendal", or "spearman"
#' is specified, any other arguments in the \code{\link{cor.test}} function can
#' be used. If "lm" is specified, a design matrix \strong{must} be passed. The
#' design matrix can be built using the \code{\link{model.matrix}} function.
#'
#' @return A list of data frames. The length of the list equals to the number
#' of rows of X, and the number of rows of each data fram equals to the number
#' of rows of Y. Each data frame has 4 columns.
#'
#' \item{\strong{stat}}{the value of the test statistic.}
#' \item{\strong{estimate}}{the estimated measure of association. Corresponds to the "r"
#' for pearson, "rho" for spearman, "tau" for kendal, and coefficient for lm}
#' \item{\strong{pval}}{the p-value of the test.}
#' \item{\strong{padj}}{the p-value of the test after multiple test adjustment.}
#'
#' @export
#' @importFrom plyr "adply"
#' @importFrom tibble "column_to_rownames"
#' @author Chenghao Zhu
MatCor = function(X, Y, method, adjust.method = "BH", ...){
    apply(X, 1, function(x){
        dat = adply(Y, 1, function(y){
            corTest(x=x, y=y, method = method, ...)
        }) %>%
            mutate(padj = p.adjust(pval, method = adjust.method)) %>%
            tibble::column_to_rownames("X1")
    })
}

################################################################################
#' @title Compute correlation in mutiple methods
#' @description Compute correlation tests between two matrices using multiple
#' methods and return a nested list. The length of the first level list equals
#' to the numebr of methods input. Each member of the first level list is a
#' output of the \code{\link{MatCor}} function. The methods supported are
#' pearson, spearman, kendall, and lm.
#'
#' @param X numaric matrix of data values.
#' @param Y numeric matrix of data values. Must have same number of columns as X.
#' @param method character value indicating the method to use. Must be one of
#' "pearson", "spearman", "kendall", or "lm".
#' @param adjust.method character value indicating the method to use for
#' p value adjustement. See \code{\link{p.adjust}}.
#' @param design a design matrix \strong{must} if "lm" is in the methods. The
#' design matrix can be built using the \code{\link{model.matrix}} function.
#'
#' @return A nested list, with the methods being the first level, X variables
#' being the second, and Y variables the bottom level.
#' @export
#' @author Chenghao Zhu
MatCorPack = function(X, Y, methods = c("pearson", "spearman", "kendall", "lm"), adjust.method = "BH", design){
    if(any(!methods %in% c("pearson", "spearman", "kendall", "lm")))
        stop("[ MatCorR::MatCorPack ] only support pearson, spearman, kendall, and lm",
             .call = FALSE)
    if("lm" %in% methods & missing(design))
        stop("[ MatCorR::MatCorPack ] the design matrix must be given to use lm method",
             .call = FALSE)
    li = lapply(methods, function(method){
        args = list(X=X, Y=Y, adjust.method = adjust.method, method = method)
        if(method == "lm") args$design = design
        do.call(MatCor, args)
    })
    names(li) = methods
    return(li)
}
