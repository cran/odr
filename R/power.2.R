#' Required budget, power, MDES, required sample size calculation for two-level GRTs
#'
#' @description This function can calculate required budget for desired power,
#'     power or the minimum detectable effect size (MDES) under fixed budget
#'     for two-level group randomized trials (GRTs).
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, and MDES calculation).
#'
#' @param expr returns from the function of \code{\link{od.2}}; default is NULL;
#'     if \code{expr} is specified, parameter values of \code{ICC},
#'     \code{R12}, \code{R22},
#'     \code{c1}, \code{c2}, \code{c1t}, \code{c2t}, \code{n}, and \code{p}
#'     used or solved in the function of \code{\link{od.2}} will
#'     be passed to the current function;
#'      only the values of \code{n} and \code{p} that specified or solved in
#'      the function of \code{\link{od.2}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint specify the values of n and/or p to overwrite those
#'     from \code{expr}; default is NULL.
#' @param es effect size.
#' @param power statistical power.
#' @param m total budget.
#' @param ICC the unconditional intraclass correlation coefficient in population or in
#'     each treatment condition.
#' @param R12 the proportion level-1 variance explained by covariates in population or in
#'     each treatment condition.
#' @param R22 the proportion level-2 variance explained by covariates in population or in
#'     each treatment condition.
#' @param c1 the cost of sampling one level-1 unit in control condition.
#' @param c2 the cost of sampling one level-2 unit in control condition.
#' @param c1t the cost of sampling one level-1 unit in treatment condition.
#' @param c2t the cost of sampling one level-2 unit in treatment condition.
#' @param n the average level-1 sample size per level-2 unit.
#' @param J the level-2 sample size across treatment conditions.
#' @param p the proportion of level-2 units assigned to treatment.
#' @param q the number of level-2 covariates, default is 0.
#' @param eslim the range to search the root of effect size (es) numerically,
#'     default is c(0, 5)
#' @param powerlim the range for searching the root of power (power) numerically,
#'     default is c(1e-10, 1 - 1e-10)
#' @param Jlim the range for searching the root of level-2 sample size (J) numerically,
#'     default is c(4, 10e10)
#' @param mlim the range for searching the root of budget (m) numerically,
#'     default is the costs for sampling same number of level-2 units across treatment conditions
#'     or c(4 * Jcost
#'     , 10e10 * Jcost), with Jcost = ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))
#' @param sig.level significance level or type I error rate, default is 0.05.
#' @param two.tailed logical; TRUE for two-tailed tests,
#'     FALSE for one-tailed tests; default is TRUE.
#' @param cost.model logical, perform power analyses accommodating costs and budget
#'     (e.g., required budget for desired power, power/MDES under fixed budget)
#'     if TRUE, perform conventional power analyses
#'     (e.g., required sample size, power, or MDES) if FALSE; default is TRUE.
#'
#' @return Required budget / sample size, statistical power, or MDES
#'     depending on the parameter specification.
#'     The function also returns the function name, design type,
#'     and the list of parameters used in the calculation.
#'
#' @export power.2
#'
#' @references Shen, Z., & Kelcey, B. (under review). Optimal design of cluster
#'     randomized trials under condition- and unit-specific cost structures.
#'     2018 American Educational Research Association (AERA) annual conference.
#'
#' @examples
#' # unconstrained optimal design
#' myod1 <- od.2(ICC = 0.2, R12 = 0.5, R22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50)
#' myod1$out   # n = 8.9, p = 0.33
#'
#' # ------- power analyses by default considering costs and budget -------
#' mym <- power.2(expr = myod1, es = 0.3, q = 1, power = 0.8)
#' mym$out  # m =1702
#'
#' mymdes <- power.2(expr = myod1, q = 1, power = 0.80, m = 1702)
#' mymdes$out  # MDES = 0.30
#'
#' mypower <- power.2(expr = myod1, q = 1, es = 0.3, m = 1702)
#' mypower$out  # power = 0.80
#'
#' # ------- conventional power analyses with cost.model = FALSE-------
#' # Required J
#' myJ1 <- power.2(cost.model = FALSE, expr = myod1, es = 0.3, q = 1, power = 0.8)
#' myJ1$out  # J = 59
#' myJ1$par  # parameters and their values used for the function
#'
#' # or equivalently, specify every arguments in the function
#' myJ1.1 <- power.2(cost.model = FALSE, es = 0.3, power = 0.8, ICC = 0.2,
#'                   R12 = 0.5, R22 = 0.5, n = 9, p = 0.33, q = 1)
#' myJ1.1$out
#'
#' # power
#'  mypower1 <- power.2(cost.model = FALSE, expr = myod1, J = 59, es = 0.3, q = 1)
#' mypower1$out  # power = 0.80
#'
#' # MDES
#' mymdes1 <- power.2(cost.model = FALSE, expr = myod1, J = 59, power = 0.8, q = 1)
#' mymdes1$out  # es = 0.30
#'
#'
power.2 <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                    sig.level = 0.05, two.tailed = TRUE,
                    es = NULL, n = NULL, J = NULL, p = NULL,
                    m = NULL, power = NULL, q = NULL,
                    ICC = NULL, R12 = NULL, R22 = NULL,
                    c1 = NULL, c2 = NULL, c1t = NULL, c2t = NULL,
                    eslim = NULL, mlim = NULL, powerlim = NULL,
                    Jlim = NULL) {
  funName <- "power.2"
  designType <- "CRT.2"
  if (cost.model) {
    if (sum(sapply(list(m, es, power), is.null)) != 1)
      stop("exactly one of m, es, and power must be NULL
           when cost.model is TRUE")
    if (!is.null(J))
      stop("'J' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(J, es, power), is.null)) != 1)
      stop("exactly one of J, es, and power must be NULL
           when cost.model is FALSE")
    if (!is.null(m))
      stop("'m' must be NULL when cost.model is FALSE")
  }
  if (!is.null(expr)) {
    if (expr$funName != "od.2") {
      stop("'expr' can only be NULL or
            the return from the function of 'od.2'")
    } else {
      if (sum(sapply(list(ICC, R12, R22, c1, c2, c1t, c2t, n, p),
                     function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'ICC', 'R12', 'R22', 'c1', 'c2',
             'c1t', 'c2t', 'n', 'p'
             have been specified in expr of 'od.2'")
      ICC <- expr$par$ICC
      R12 <- expr$par$R12
      R22 <- expr$par$R22
      c1 <- expr$par$c1
      c2 <- expr$par$c2
      c1t <- expr$par$c1t
      c2t <- expr$par$c2t
      n <- round(expr$out$n, 0)
      p <- round(expr$out$p, 2)
    }
  } else {
    if (!is.null(constraint))
      stop("'constraint' must be NULL when 'expr' is NULL")
    }
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (!is.null(constraint) && !is.list(constraint))
    stop("'constraint' must be in list format
         (e.g., constraint = list(p = 0.5))")
  if (length(constraint) > 2)
    stop("'constraint' must be limited to 'n' and/or 'p'")
  if (!is.null(constraint$n)) {
    if(NumberCheck(constraint$n) || constraint$n < 1)
      stop("constrained 'n' must be numeric in [1, inf)")
    n <- constraint$n
  }
  if (!is.null(constraint$p)) {
    if(NumberCheck(constraint$p) ||
       any (0 >= constraint$p | constraint$p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p <- constraint$p
  }
  if (sum(sapply(list(ICC, p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'ICC', 'p', 'power', and 'sig.level'
                 must be numeric in (0, 1)")
  if (sum(sapply(list(R12, R22), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'R12', 'R22' must be numeric in [0, 1)")
  if (cost.model){
   if (sum(sapply(list(c1, c2, c1t, c2t), function(x) {
    NumberCheck(x) || x < 0})) >= 1)
    stop("'c1', 'c2', 'c1t', 'c2t' must be numeric in [0, inf)")
   if (NumberCheck(m))
    stop("'m' must be numeric in [0, inf)")
   }
  if (NumberCheck(q) | q < 0)
   stop("'q' must be numeric in [0, inf)")
  if (NumberCheck(n) || n <= 1)
    stop("'n' must be numeric in (1, inf)")
  if (NumberCheck(es) || any(0 > es | es > 3))
    stop("'es' must be numeric in [0, 3],
         please transfer negative effect size to positive one if needed")
  if (R22 > 0 && q == 0)
    stop("'q' must be q >= 1 when R22 != 0")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              es = es, ICC = ICC, R12 = R12, R22 = R22,
              c1 = c1, c2 = c2, c1t = c1t, c2t = c2t,
              n = n, J = J, p = p,
              q = q, m = m, power = power)
  if (!(cost.model)) {
    par$c1 <- par$c2 <- par$c1t <- par$c2t <- NULL
  }
  tside <- ifelse(two.tailed, 2, 1)
  if (cost.model) {
    Jcost <- ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))
    if (two.tailed) {
      p.expr <- quote({
        J <- m / ((1 - p) * (c1 * n + c2)
                  + p * (c1t * n + c2t));
        lamda <- es * sqrt(p * (1 - p) * J) /
          sqrt(ICC * (1 - R22) + (1 - ICC) * (1 - R12) / n);
        1 - pt(qt(1 - sig.level / tside, df = J - q - 2) ,
               df = J - q - 2, lamda) +
          pt(qt(sig.level / tside, df = J - q - 2),
             df = J - q - 2, lamda)
      })
    } else {
      p.expr <- quote({
        J <- m / ((1 - p) * (c1 * n + c2)
                  + p * (c1t * n + c2t));
        lamda <- es * sqrt(p * (1 - p) * J) /
          sqrt(ICC * (1 - R22) + (1 - ICC) * (1 - R12) / n);
        1 - pt(qt(1 - sig.level / tside, df = J - q - 2),
               df = J - q - 2, lamda)
      })
    }
  } else {
    if (two.tailed) {
      p.expr <- quote({
        lamda <- es * sqrt(p * (1 - p) * J) /
          sqrt(ICC * (1 - R22) + (1 - ICC) * (1 - R12) / n);
        1 - pt(qt(1 - sig.level / tside, df = J - q - 2),
               df = J - q - 2, lamda) +
          pt(qt(sig.level / tside, df = J - q - 2),
             df = J - q - 2, lamda)
      })
    } else {
      p.expr <- quote({
        lamda <- es * sqrt(p * (1 - p) * J) /
          sqrt(ICC * (1 - R22) + (1 - ICC) * (1 - R12) / n);
        1 - pt(qt(1 - sig.level / tside, df = J - q - 2),
               df = J - q - 2, lamda)
      })
    }
  }
  par$p.expr <- p.expr
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  Jlim <- limFun(x = Jlim, y = c(4, 10e10))
  powerlim <- limFun(x = powerlim, y = c(1e-10, 1 - 1e-10))
  eslim <- limFun(x = eslim, y = c(0, 5))
  if (cost.model)
    mlim <- limFun(x = mlim, y = c(Jlim[1] * Jcost, Jlim[2] * Jcost))
  if(cost.model){
    if (is.null(power)) {
      out <- list(power = eval(p.expr))
    } else if (is.null(m)) {
      out <- list(m = stats::uniroot(function(m)
        eval(p.expr) - power, mlim)$root)
    } else if (is.null(es)) {
      out <- list(es = stats::uniroot(function(es)
        eval(p.expr) - power, eslim)$root)
    }
  } else {
    if (is.null(power)) {
      out <- list(power = eval(p.expr))
    } else if (is.null(J)) {
      out <- list(J = stats::uniroot(function(J)
        eval(p.expr) - power, Jlim)$root)
    } else if (is.null(es)) {
      out <- list(es = stats::uniroot(function(es)
        eval(p.expr) - power, eslim)$root)
    }
  }
  power.out <- list(funName = funName,
                    designType = designType,
                    par = par, out = out)
  class(power.out) <- c("pars")
  return(power.out)
}

