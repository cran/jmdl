#' @import methods
NULL

#' Class "jmdlMod" of Fitted Joint Mean-Correlation Models.
#'
#' @slot call the matched call
#' @slot opt the optimization result returned by optimizeJmdl
#' @slot args arguments m, Y, X, W, time
#' @slot q degree of polynomial of the time lag to model the lower triangular matrix
#' @slot std standard error for parameters
#' @slot tval t statistic components list
#' @slot p p.value components list
#' @slot devcomp the deviance components list
#' @slot offset a vecter to be added to a linear predictor
#'
#' @exportClass jmdlMod
setClass("jmdlMod", representation( call = "call", opt = "list", args = "list", std = "list",
                                    tval = "list",  p = "list", q = "numeric", devcomp = "list",
                                    offset = "numeric"))
