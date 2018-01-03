#' Required budget (and/or sample size), power, MDES calculation for four-level CRTs
#'
#' @description This function can calculate required budget for desired power,
#'     power or minimum detectable effect size (MDES) under fixed budget
#'     for four-level cluster randomized trials (CRTs).
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, and MDES calculation).
#'
#' @param cost.model logical, perform power analyses accommodating costs and budget
#'     (e.g., required budget for desired power, power/MDES under fixed budget)
#'     if TRUE, perform conventional power analyses
#'     (e.g., required sample size, power, or MDES calculation) if FALSE; default is TRUE.
#' @param expr returned objects from function \code{\link{od.4}}; default is NULL;
#'     if \code{expr} is specified, parameter values of \code{icc2}, \code{icc3}, \code{icc4},
#'     \code{r12}, \code{r22}, \code{r32}, \code{r42},
#'     \code{c1}, \code{c2}, \code{c3}, \code{c4},
#'     \code{c1t}, \code{c2t}, \code{c3t}, \code{c4t}, \code{p}, \code{n}, \code{J}, and \code{K}
#'     used or solved in function \code{\link{od.4}} will
#'     be passed to current function;
#'     only the values of \code{p}, \code{n}, \code{J}, and/or \code{K} that specified or solved in
#'     function \code{\link{od.4}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint specify the constrained values of \code{p}, \code{n}, \code{J},
#'     and/or \code{K} in list format to overwrite
#'     those from \code{expr}; default is NULL.
#' @param sig.level significance level or type I error rate, default is 0.05.
#' @param two.tailed logical; perform two-tailed tests if TRUE,
#'     perform one-tailed tests if FALSE; default is TRUE.
#' @param es effect size.
#' @param power statistical power.
#' @param m total budget.
#' @param icc2 the unconditional intraclass correlation coefficient (icc) at level 2.
#' @param icc3 the unconditional intraclass correlation coefficient (icc) at level 3.
#' @param icc4 the unconditional intraclass correlation coefficient (icc) at level 4.
#' @param r12 the proportion of level-1 variance explained by covariates.
#' @param r22 the proportion of level-2 variance explained by covariates.
#' @param r32 the proportion of level-3 variance explained by covariates.
#' @param r42 the proportion of level-4 variance explained by covariates.
#' @param c1 the cost of sampling one level-1 unit in control condition.
#' @param c2 the cost of sampling one level-2 unit in control condition.
#' @param c3 the cost of sampling one level-3 unit in control condition.
#' @param c4 the cost of sampling one level-4 unit in control condition.
#' @param c1t the cost of sampling one level-1 unit in treatment condition.
#' @param c2t the cost of sampling one level-2 unit in treatment condition.
#' @param c3t the cost of sampling one level-3 unit in treatment condition.
#' @param c4t the cost of sampling one level-4 unit in treatment condition.
#' @param n the average level-1 sample size per level-2 unit.
#' @param J the average level-2 sample size per level-3 unit.
#' @param K the average level-3 sample size per level-4 unit.
#' @param L the level-4 sample size across treatment conditions.
#' @param p the proportion of level-4 clusters/units assigned to treatment.
#' @param q the number of covariates at level 4.
#' @param eslim the range for searching the root of effect size (\code{es}) numerically,
#'     default is c(0, 5).
#' @param powerlim the range for searching the root of power (\code{power}) numerically,
#'     default is c(1e-10, 1 - 1e-10).
#' @param Llim the range for searching the root of level-4 sample size (\code{L}) numerically,
#'     default is c(4, 10e10).
#' @param mlim the range for searching the root of budget (\code{m}) numerically,
#'     default is the costs sampling \code{Llim} level-4 units across treatment conditions
#'     or c(4 * Lcost, 10e10 * Lcost), with Lcost =
#'     ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4) +
#'     p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t)).
#' @param rounded logical; round the values of \code{p}, \code{n}/\code{J}/\code{K}
#'     that are from functions \code{\link{od.4}}
#'     to two decimal places and integer, respectively if TRUE,
#'     no rounding if FALSE; default is TRUE.
#'
#' @return Required budget (and/or required level-4 sample size), statistical power, or MDES
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.4
#'
#' @references
#'   Shen, Z.(in progress). Using optimal sample allocation to
#'   improve design efficiency for multilevel randomized trials.
#'   (Unpublished doctoral dissertation). University of Cincinnati, Cincinnati, OH.
#'
#' @examples
#' # unconstrained optimal design
#' myod1 <- od.4(icc2 = 0.2, icc3 = 0.1, icc4 = 0.05,
#'               r12 = 0.5, r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'               c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'               c1t = 1, c2t = 50, c3t = 250, c4t = 2500)
#' myod1$out # output # n = 7.1, J = 3.2, K = 4.2, p = 0.23
#'
#' # ------- power analyses by default considering costs and budget -------
#' # required budget and sample size
#' mym.1 <- power.4(expr = myod1, es = 0.3, q = 1, power = 0.8)
#' mym.1$out  # m = 33188, L = 26.6
#' mym.1$par  # parameters and their values used for the function
#' # or equivalently, specify every argument in the function
#' mym.1 <- power.4(es = 0.3, power = 0.8, q = 1,
#'                  icc2 = 0.2, icc3 = 0.1, icc4 = 0.05,
#'                  r12 = 0.5, r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'                  c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'                  c1t = 1, c2t = 50, c3t = 250, c4t = 2500,
#'                  n = 7, J = 3, K = 4, p = 0.23)
#' # required budget and sample size with constrained p
#' mym.2 <- power.4(expr = myod1, es = 0.3, q = 1, power = 0.8,
#'                  constraint = list(p = 0.5))
#' mym.2$out  # m = 44584, L = 19.6
#' # required budget and sample size with constrained p and K
#' mym.3 <- power.4(expr = myod1, es = 0.3, q = 1, power = 0.8,
#'                  constraint = list(p = 0.5, K = 20))
#' mym.3$out  # m = 79108, L = 12.9
#'
#' # Power calculation
#' mypower <- power.4(expr = myod1, q = 1, es = 0.3, m = 33188)
#' mypower$out  # power = 0.80
#' # Power calculation under constrained p (p = 0.5)
#' mypower.1 <- power.4(expr = myod1, q = 1, es = 0.3, m = 33188,
#'                  constraint = list(p = 0.5))
#' mypower.1$out  # power = 0.65
#'
#' # MDES calculation
#' mymdes <- power.4(expr = myod1, q = 1, power = 0.80, m = 33188)
#' mymdes$out  # MDES = 0.30
#'
#'
#' # ------- conventional power analyses with cost.model = FALSE-------
#' # Required sample size
#' myL <- power.4(cost.model = FALSE, expr = myod1, es = 0.3, q = 1, power = 0.8)
#' myL$out  # L = 26.6
#' myL$par  # parameters and their values used for the function
#' # or equivalently, specify every argument in the function
#' myL <- power.4(cost.model = FALSE, es = 0.3, power = 0.8, q = 1,
#'                   icc2 = 0.2, icc3 = 0.1, icc4 = 0.05,
#'                   r12 = 0.5, r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'                   n = 7, J = 3, K = 4, p = 0.23)
#'
#' # Power calculation
#'  mypower1 <- power.4(cost.model = FALSE, expr = myod1, L = 27, es = 0.3, q = 1)
#'  mypower1$out  # power = 0.806
#'
#' # MDES calculation
#' mymdes1 <- power.4(cost.model = FALSE, expr = myod1, L = 27, power = 0.8, q = 1)
#' mymdes1$out  # es = 0.30
#'
#'
power.4 <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                    sig.level = 0.05, two.tailed = TRUE,
                    es = NULL, power = NULL, m = NULL,
                    n = NULL, J = NULL, K = NULL, L = NULL, p = NULL,
                    icc2 = NULL, icc3 = NULL, icc4 = NULL,
                    r12 = NULL, r22 = NULL, r32 = NULL, r42 = NULL, q = NULL,
                    c1 = NULL, c2 = NULL, c3 = NULL, c4 = NULL,
                    c1t = NULL, c2t = NULL, c3t = NULL, c4t = NULL,
                    eslim = NULL, powerlim = NULL, Llim = NULL, mlim = NULL,
                    rounded = TRUE) {
  funName <- "power.4"
  designType <- "four-level CRTs"
  if (cost.model) {
    if (sum(sapply(list(m, es, power), is.null)) != 1)
      stop("exactly one of m, es, and power must be NULL
           when cost.model is TRUE")
    if (!is.null(L))
      stop("'L' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(L, es, power), is.null)) != 1)
      stop("exactly one of L, es, and power must be NULL
           when cost.model is FALSE")
    if (!is.null(m))
      stop("'m' must be NULL when cost.model is FALSE")
  }
  if (!is.null(expr)) {
    if (expr$funName != "od.4") {
      stop("'expr' can only be NULL or
           the return from the function of 'od.4'")
    } else {
      if (sum(sapply(list(icc2, icc3, icc4, r12, r22, r32, r42, c1, c2, c3, c4,
                          c1t, c2t, c3t, c4t, n, J, K, p), function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'icc2', 'icc3', 'icc4', 'r12', 'r22', 'r32', 'r42',
             'c1', 'c2', 'c3', 'c4', 'c1t', 'c2t', 'c3t', 'c4t', 'n', 'J', 'K', and 'p'
             have been specified in expr of 'od.4'")
      icc2 <- expr$par$icc2
      icc3 <- expr$par$icc3
      icc4 <- expr$par$icc4
      r12 <- expr$par$r12
      r22 <- expr$par$r22
      r32 <- expr$par$r32
      r42 <- expr$par$r42
      c1 <- expr$par$c1
      c2 <- expr$par$c2
      c3 <- expr$par$c3
      c4 <- expr$par$c4
      c1t <- expr$par$c1t
      c2t <- expr$par$c2t
      c3t <- expr$par$c3t
      c4t <- expr$par$c4t
      if (rounded) {
        n <- round(expr$out$n, 0)
        J <- round(expr$out$J, 0)
        K <- round(expr$out$K, 0)
        p <- round(expr$out$p, 2)
      } else {
        n <- expr$out$n
        J <- expr$out$J
        K <- expr$out$K
        p <- expr$out$p
      }
    }
  } else {
    if (!is.null(constraint))
      stop("'constraint' must be NULL when 'expr' is NULL")
  }
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (!is.null(constraint) && !is.list(constraint))
    stop("'constraint' must be in list format
         (e.g., constraint = list(p = 0.5))")
  if (length(constraint) > 4)
    stop("'constraint' must be limited to 'n', 'J', 'K', and/or 'p'")
  if (!is.null(constraint$n)) {
    if(NumberCheck(constraint$n) || constraint$n <= 0)
      stop("constrained 'n' must be numeric in (0, 10e10)")
    n <- constraint$n
  }
  if (!is.null(constraint$J)) {
    if(NumberCheck(constraint$J) || constraint$J <= 0)
      stop("constrained 'J' must be numeric in (0, 10e10)")
    J <- constraint$J
  }
  if (!is.null(constraint$K)) {
    if(NumberCheck(constraint$K) || constraint$K <= 0)
      stop("constrained 'K' must be numeric in (0, 10e10)")
    K <- constraint$K
  }
  if (!is.null(constraint$p)) {
    if(NumberCheck(constraint$p) ||
       any (0 >= constraint$p | constraint$p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p <- constraint$p
  }
  if (sum(sapply(list(icc2, icc3, icc4, p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'icc2', 'icc3', 'icc4', 'p', 'power', and 'sig.level'
                 must be numeric in (0, 1)")
  if (sum(sapply(list(r12, r22, r32, r42), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'r12', 'r22', 'r32', and 'r42' must be numeric in [0, 1)")
  if (cost.model){
    if (sum(sapply(list(c1, c2, c3, c4, c1t, c2t, c3t, c4t), function(x) {
      NumberCheck(x) || x < 0})) >= 1)
      stop("'c1', 'c2', 'c3', 'c4', 'c1t', 'c2t', 'c3t', 'c4t' must be numeric in [0, Inf)")
    if (NumberCheck(m))
      stop("'m' must be numeric in [0, Inf)")
  }
  if (NumberCheck(q) | q < 0)
    stop("'q' must be numeric in [0, 10e3]")
  if (NumberCheck(n) || n <= 0)
    stop("'n' must be numeric in (0, 10e10)")
  if (NumberCheck(J) || J <= 0)
    stop("'J' must be numeric in (0, 10e10)")
  if (NumberCheck(K) || K <= 0)
    stop("'K' must be numeric in (0, 10e10)")
  if (NumberCheck(es) || any(0 > es | es > 5))
    stop("'es' must be numeric in [0, 5],
         please transfer negative effect size to positive one if needed")
  if (r42 > 0 && q == 0)
    stop("'q' must be q >= 1 when r42 != 0")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              es = es, icc2 = icc2, icc3 = icc3, icc4 = icc4,
              r12 = r12, r22 = r22, r32 = r32, r42 = r42,
              c1 = c1, c2 = c2, c3 = c3, c4 = c4,
              c1t = c1t, c2t = c2t, c3t = c3t, c4t = c4t,
              n = n, J = J, K = K, L = L, p = p,
              q = q, m = m, power = power)
  tside <- ifelse(two.tailed, 2, 1)
  if (cost.model) {
    Lcost <- ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4) +
                p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t))
    if (two.tailed) {
      p.expr <- quote({
        L <- m / ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4) +
                    p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t));
        lambda <- es * sqrt(p * (1 - p) * m /
                              ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4) +
                                 p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t))) /
          sqrt(icc4 * (1 - r42) + icc3 * (1 - r32) / K + icc2 * (1 - r22) / (J * K ) +
                 (1 - icc2 - icc3 - icc4) * (1 - r12) / (n * J * K));
        1 - pt(qt(1 - sig.level / tside, df = L - q - 2) ,
               df = L - q - 2, lambda) +
          pt(qt(sig.level / tside, df = L - q - 2),
             df = L - q - 2, lambda)
      })
    } else {
      p.expr <- quote({
        L <- m / ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4) +
                    p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t));
        lambda <- es * sqrt(p * (1 - p) * m /
                              ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4) +
                                 p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t))) /
          sqrt(icc4 * (1 - r42) + icc3 * (1 - r32) / K + icc2 * (1 - r22) / (J * K ) +
                 (1 - icc2 - icc3 - icc4) * (1 - r12) / (n * J * K));
        1 - pt(qt(1 - sig.level / tside, df = L - q - 2),
               df = L - q - 2, lambda)
      })
    }
  } else {
    if (two.tailed) {
      p.expr <- quote({
        lambda <- es * sqrt(p * (1 - p) * L) /
          sqrt(icc4 * (1 - r42) + icc3 * (1 - r32) / K + icc2 * (1 - r22) / (J * K ) +
                 (1 - icc2 - icc3 - icc4) * (1 - r12) / (n * J * K));
        1 - pt(qt(1 - sig.level / tside, df = L - q - 2),
               df = L - q - 2, lambda) +
          pt(qt(sig.level / tside, df = L - q - 2),
             df = L - q - 2, lambda)
      })
    } else {
      p.expr <- quote({
        lambda <- es * sqrt(p * (1 - p) * L) /
          sqrt(icc4 * (1 - r42) + icc3 * (1 - r32) / K + icc2 * (1 - r22) / (J * K ) +
                 (1 - icc2 - icc3 - icc4) * (1 - r12) / (n * J * K));
        1 - pt(qt(1 - sig.level / tside, df = L - q - 2),
               df = L - q - 2, lambda)
      })
    }
  }
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  Llim <- limFun(x = Llim, y = c(4, 10e10))
  powerlim <- limFun(x = powerlim, y = c(1e-10, 1 - 1e-10))
  eslim <- limFun(x = eslim, y = c(0, 5))
  if(cost.model) {
    mlim <- limFun(x = mlim, y = c(Llim[1] * Lcost, Llim[2] * Lcost))
    if (is.null(power)) {
      out <- list(power = eval(p.expr))
    } else if (is.null(m)) {
      out <- list(m = stats::uniroot(function(m)
        eval(p.expr) - power, mlim)$root)
      out <- c(out, list(L = out$m / ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4) +
                                        p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t))))
    } else if (is.null(es)) {
      out <- list(es = stats::uniroot(function(es)
        eval(p.expr) - power, eslim)$root)
    }
  } else {
    if (is.null(power)) {
      out <- list(power = eval(p.expr))
    } else if (is.null(L)) {
      out <- list(L = stats::uniroot(function(L)
        eval(p.expr) - power, Llim)$root)
    } else if (is.null(es)) {
      out <- list(es = stats::uniroot(function(es)
        eval(p.expr) - power, eslim)$root)
    }
  }
  power.out <- list(funName = funName,
                    designType = designType,
                    par = par, out = out)
  return(power.out)
}
