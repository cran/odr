#' Optimal sample allocation calculation for two-level GRTs
#'     under fixed budget and cost structure
#'
#' @description The optimal design of two-level
#'     group randomized trials (GRTs) is to choose
#'     the average level-1 sample size per level-2 unit (\code{n})
#'     and the proportion of groups assigned to treatment (\code{p})
#'     such that the variances of treatment effects
#'     are minimized under fixed budget and cost structure.
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without constraints.
#'
#' @inheritParams power.2
#' @param m total budget, default is the total costs of sampling 60
#'     level-2 units across treatment conditions.
#' @param plot logical, print plot if TRUE, otherwise not, default is TRUE.
#' @param plot.by specify plot by n and/or p; default is both.
#' @param nlim the plot range for n, default is c(2, 50)
#' @param plim the plot range for p, default is c(0, 1)
#' @param varlim the plot range of variance, default is c(0, 0.2)
#' @param nlab the plot label for n,
#'     default is "Level-1 Sample Size: n".
#' @param plab the plot label for p,
#'     default is "Proportion Level-2 Units in Treatment: p"
#' @param varlab the plot label for variance,
#'     default is "Variance".
#' @param vartitle the title of variance plots, default is NULL.
#' @param verbose logical; print the values of n and p if TRUE,
#'    otherwise not; default is TRUE.
#' @return
#'     \code{n} and \code{p}. The function also returns
#'     the function name, design type,
#'     and the list of parameters used in the calculation.
#'
#' @export od.2
#'
#' @references Shen, Z., & Kelcey, B. (under review). Optimal design of cluster
#'     randomized trials under condition- and unit-specific cost structures.
#'     2018 American Educational Research Association (AERA) annual conference.
#'
#' @examples
#' # unconstrained optimal design
#' myod1 <- od.2(ICC = 0.2, R12 = 0.5, R22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               varlim = c(0.01, 0.02))
#' myod1$out # output
#'
#' # constrained optimal design with n = 20
#' myod2 <- od.2(ICC = 0.2, R12 = 0.5, R22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               n = 20, varlim = c(0.005, 0.025))
#' myod2$out
#'
#' # constrained optimal design with p = 0.5
#' myod3 <- od.2(ICC = 0.2, R12 = 0.5, R22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'              p = 0.5, varlim = c(0.005, 0.025))
#' myod3$out
#'
#' # constrained n and p, no calculation performed
#' myod4 <- od.2(ICC = 0.2, R12 = 0.5, R22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               n = 20, p = 0.5, varlim = c(0.005, 0.025))
#' myod4$out
#'
od.2 <- function(n = NULL, p = NULL, ICC = NULL, R12 = NULL, R22 = NULL,
                 c1 = NULL, c2 = NULL, c1t = NULL, c2t = NULL, m = NULL,
                 plot.by = NULL, nlim = NULL, plim = NULL, varlim = NULL,
                 nlab = NULL, plab = NULL, varlab = NULL, plot = TRUE,
                 vartitle = NULL,verbose = TRUE) {
  funName <- "od.2"
  designType <- "CRT.2"
  if (sum(sapply(list(ICC, R12, R22, c1, c2, c1t, c2t),
                 function(x) is.null(x))) >= 1)
    stop("all of 'ICC', 'R12', 'R22', 'c1', 'c2',
         'c1t', 'c2t' must be specified")
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (NumberCheck(ICC) || any(0 >= ICC | ICC >= 1))
    stop("'ICC' must be numeric in (0, 1)")
  if (sum(sapply(list(R12, R22), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1)
    stop("'R12', 'R22' must be numeric in [0, 1)")
  if (sum(sapply(list(c1, c2, c1t, c2t), function(x) {
    NumberCheck(x) || x < 0})) >= 1)
    stop("'c1', 'c2', 'c1t', 'c2t' must be numeric in [0, inf)")
  if (!is.null(plot.by) && !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n'))")
  derivative <- quote({
    p <- sqrt((c1 * n + c2) / (c1t * n + c2t)) /
      (1 + sqrt((c1 * n + c2) / (c1t * n + c2t)))
    n - sqrt((1 - ICC) * (1 - R12) / (ICC * (1 - R22) )) *
      sqrt(((1 - p) * c2 + p * c2t)/(( 1 - p) * c1 + p * c1t))
  })
  par <- list(ICC = ICC, R12 = R12, R22 = R22, c1 = c1, c2 = c2,
              c1t =c1t, c2t = c2t, n = n, p = p)
  if (is.null(p) && is.null(n)) {
    n <- stats::uniroot(function(n) eval(derivative), interval = c(1, 10e10))$root
    p <- sqrt((c1 * n + c2) / (c1t * n + c2t)) /
      (1 + sqrt((c1 * n + c2) / (c1t * n + c2t)))
  } else if (!is.null(n) && is.null(p)) {
    if (!is.numeric(n) || n < 2)
      stop("constrained 'n' must be nu meric in [2, inf)")
    p <- sqrt((c1 * n + c2) / (c1t * n + c2t)) /
      (1 + sqrt((c1 * n + c2)/(c1t * n + c2t)))
  } else if (!is.null(p) && is.null(n)) {
    if (!is.numeric(p) || any(p <=0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    n <- sqrt(( 1 - ICC) * (1 - R12) / (ICC * (1 - R22))) *
         sqrt(((1 - p) * c2 + p * c2t) / ((1 - p) * c1 + p * c1t))
  } else if (!is.null(p) && !is.null(n)) {
    if (!is.numeric(n) || n < 2)
      stop("constrained 'n' must be nu meric in [2, inf)")
    if (!is.numeric(p) || any(p <=0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
      cat("===============================\n",
        "both p and n are constrained, there is no calculation from other parameters",
        ".\n===============================\n", sep = "")
    }
  if (verbose) {
    if (!is.null(par$n)) {
      cat("The constrained level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    } else {
      cat("The optimal level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    }
    if (!is.null(par$p)) {
      cat("The constrained proportion of level-2 units in treatment (p) is ", p, ".\n", "\n", sep = "")
    } else {
      cat("The optimal proportion of level-2 units in treatment (p) is ", p, ".\n", "\n" ,sep = "")
    }
  }
  m <- ifelse(!is.null(m), m, 60 * (p * (c1t * n + c2t) + (1 - p) * (c1 * n + c2)))
  if (NumberCheck(m) || m < 0) stop("'m' must be numeric in (0, inf)")
  var.expr <- quote({
    J <- m / ((1 - p) * (c1 * n + c2)
              + p * (c1t * n + c2t))
    (ICC * (1 - R22) + (1 - ICC) * (1 - R12) / n )/ (p * (1 - p) * J)
  })
  var <- eval(var.expr)
  par <- c(par, list(m = m, var = eval(var.expr)), var.expr = var.expr)
  out <- list(n = n, p = p)
  od.out <- list(funName = funName, designType = designType,
                 par = par, out = out)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(2, 50))
  plim <- limFun(x = plim, y = c(0, 1))
  varlim <- limFun(x = varlim, y = c(0, 0.2))
  labFun <- function(x, y) {
    if (!is.null(x) && length(x) == 1 && is.character(x)) {x} else {y}
  }
  nlab <- labFun(x = nlab, y = "Level-1 Sample Size: n")
  plab <- labFun(x = plab, y = "Proportion Level-2 Units in Treatment: p")
  varlab <- labFun(x = varlab, y = "Variance")
  vartitle <- labFun(x = vartitle, y = "")
  plotbyFun <- function(x, y) {
    if (!is.null(x) && is.list(x)) {x} else {y}
  }
  plot.by <- plotbyFun(x = plot.by, y = list(n = "n", p = "p"))
  nrange <- seq(nlim[1], nlim[2], by = 1)
  prange <- seq(plim[1] + 0.05, plim[2] - 0.05, by = 0.01)
  if (length(plot.by) == 2) figure <- par(mfrow = c (1, 2))
  if (length(plot.by) == 1) figure <- par(mfrow = c (1, 1))
  if (plot) {
    if (!is.null(plot.by$n)) {
      plot.y <- NULL
      for (n in nrange)
        plot.y <- c(plot.y, eval(var.expr))
      plot(nrange, plot.y,
           type = "l", lty = 1,
           xlim = nlim, ylim = varlim,
           xlab = nlab, ylab = varlab,
           main = vartitle, col = "black")
      n <- out$n
      graphics::abline(v = n, lty = 2, col = "Blue")
    }
    if (!is.null(plot.by$p)) {
      plot.y <- NULL
      for (p in prange)
        plot.y <- c(plot.y, eval(var.expr))
      plot(prange, plot.y,
           type = "l", lty = 1,
           xlim = plim, ylim = varlim,
           xlab = plab, ylab = varlab,
           main = vartitle, col = "black")
      p <- out$p
      graphics::abline(v = p, lty = 2, col = "Blue")
    }
  }
  par(figure)
  class(od.out) <- c("pars")
  return(od.out)
}
