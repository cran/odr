#' Budget and/or sample size, power, MDES calculation for CRTs probing
#'     mediation effects with cluster-level mediators
#' @description This function can calculate required budget for desired power and
#'     power under a fixed budget
#'     for experimental studies with group mediators probing mediation effects.
#'     It also can perform conventional power analyses
#'     (e.g., required sample size and power calculation).
#' @inheritParams power.2
#' @inheritParams od.2.221
#' @param expr returned object from function \code{\link{od.2.221}}; default value is NULL;
#'     if \code{expr} is specified, parameter values of \code{a}, \code{b},
#'     \code{c1}, \code{c1t}, and \code{p}
#'     used or solved in function \code{\link{od.2.221}} will
#'     be passed to the current function;
#'     only the values of \code{p} and \code{n} that specified or solved in
#'     function \code{\link{od.2.221}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint specify the constrained value of
#'     \code{p} and/or \code{n} in a list format to overwrite that/those
#'     from \code{expr}; default value is NULL.
#' @param mlim the range for searching the root of budget (\code{m}) numerically,
#'     default value is the costs sampling \code{nlim} units across treatment conditions
#'     or c(4 * ncost, 10e10 * ncost) with ncost = ((1 - p) * c1 + p * c1t)
#'
#' @return Required budget (or required sample size), statistical power, or MDES
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.2.221
#'


power.2.221 <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                        sig.level = 0.05, two.tailed = TRUE,
                        a = NULL, b = NULL,
                        power = NULL, m = NULL, test = "joint",
                        n = NULL, p = NULL,
                        c1 = NULL, c1t = NULL,
                        c2 = NULL, c2t = NULL,
                        r12 = 0, r22 = 0, r2m = 0,
                        icc = NULL, J = NULL, q = 0,
                        q.a = 0, q.b = 0,
                        powerlim = NULL, Jlim = NULL, mlim = NULL,
                        rounded = TRUE) {
  funName <- "power.2.221"
  designType <- "2-2-1 mediation in two-level CRTs"
  if (cost.model == TRUE) {
    if (sum(sapply(list(m, power), is.null)) != 1)
      stop("exactly one of 'm' and 'power' must be NULL
           when cost.model is TRUE")
    if (!is.null(J))
      stop("'J' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(J, power), is.null)) != 1)
      stop("exactly one of 'J' and 'power' must be NULL
           when cost.model is FALSE")
    if (!is.null(m))
      stop("'m' must be NULL when cost.model is FALSE")
  }
  if (!is.null(expr)) {
    if (expr$funName != "od.2.221") {
      stop("'expr' can only be NULL or
            the return from the function of 'od.2.221'")
    } else {
      if (sum(sapply(list(a, b, c1, c1t, c2, c2t, n, p),
                     function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'a', 'b', 'c1', 'c1t', 'c2', 'c2t', 'n', 'p'
             have been specified in expr of 'od.2.221'")
      a <- expr$par$a
      b <- expr$par$b
      c1 <- expr$par$c1
      c1t <- expr$par$c1t
      c2 <- expr$par$c2
      c2t <- expr$par$c2t
      icc <- expr$par$icc
      r12 <- expr$par$r12
      r22 <- expr$par$r22
      r2m <- expr$par$r2m
      if(!is.null(test)){
        if (test != expr$par$test)
        {cat('Tests are different in power and the od function ', '(', test, ' and ', expr$par$test,
            '). \n', 'The power analysis is for the ', test, ' test.', sep = "")}
      } else {
        test <- expr$par$test
      }
      n <- expr$out$n
      p <- expr$out$p
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
    stop("'constraint' must be limited to 'p' and 'n'")
  if (!is.null(constraint$p)) {
    if(NumberCheck(constraint$p) ||
       any (0 >= constraint$p | constraint$p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p <- constraint$p
  }
  if (!is.null(constraint$n)) {
    if(NumberCheck(constraint$n) ||
       (0 >= constraint$n))
      stop("constrained 'n' must be numeric in (0.05, 1e10)")
    n <- constraint$n
  }
  if (sum(sapply(list(p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'p', 'power', and 'sig.level' must be numeric in (0, 1]")

  if (sum(sapply(list(b), function(x) {
    NumberCheck(x) || any(-1 >= x | x >= 1)
  })) >= 1) stop("'b' must be numeric in (-1, 1)")
  if (cost.model == TRUE){
    if (sum(sapply(list(c1, c1t, c2, c2t), function(x) {
      NumberCheck(x) || x < 0})) >= 1)
      stop("'c1', 'c1t', 'c2', 'c2t' must be numeric in [0, Inf)")
    if (NumberCheck(m))
      stop("'m' must be numeric in [0, Inf)")
  }
  if (NumberCheck(a) || any(-5 > a | a > 5))
    stop("'a' must be numeric in [-5, 5]")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              a = a, b = b,
              r12 = r12, r22 = r22, r2m = r2m,
              c1 = c1, c1t = c1t, c2 = c2, c2t = c2t,
              n = n, p = p, J = J, funName = funName,
              q = q, m = m, power = power)
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  Jlim <- limFun(x = Jlim, y = c(6, 1e6))
  powerlim <- limFun(x = powerlim, y = c(5e-10, 1 - 5e-10))
  mlim <- limFun(x = mlim, y = c(4 * ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t)),
                                 1e10 * ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))))
  if (test == "sobel" | test == "Sobel" | test == "SOBEL"){
    if (cost.model){
      # when cost.model is true for the Sobel test
      if (two.tailed) {
        pwr.sobel <- quote({
          1 - pnorm(qnorm(1 - sig.level / tside),
                    a * b / sqrt(a^2 * (icc * (1 - r22) + (1 - icc) * (1 - r12)/n)/
                                   (m * (1 - r2m) / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))) +
                                   b^2 * (1 -r2m) / (p * (1 - p) *
                                                       (m / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t)))))) +
            pnorm(qnorm(sig.level / tside),  a * b / sqrt(a^2 * (icc * (1 - r22) + (1 - icc) * (1 - r12)/n)/
                                                            (m * (1 - r2m) / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))) +
                                                            b^2 * (1 -r2m) / (p * (1 - p) * (m / ((1 - p) *
                                                                                                    (c1 * n + c2) + p * (c1t * n + c2t))))))
        })} else {
          pwr.sobel <- quote({
            1 - pnorm(qnorm(1 - sig.level / tside),
                      a * b / sqrt(a^2 * (icc * (1 - r22) + (1 - icc) * (1 - r12)/n)/
                                     (m * (1 - r2m) / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))) +
                                     b^2 * (1 -r2m) / (p * (1 - p) *
                                                         (m / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))))))
          })}
    } else {
      # when cost.model is not true for the Sobel test
      if (two.tailed) {
        pwr.sobel <- quote({
          1 - pnorm(qnorm(1 - sig.level / tside),
                    a * b / sqrt(a^2 * (icc * (1 - r22) + (1 - icc) * (1 - r12)/n)/
                                   (J* (1 - r2m)) +
                                   b^2 * (1 -r2m) / (p * (1 - p) * J))) +
            pnorm(qnorm(sig.level / tside),   a * b / sqrt(a^2 * (icc * (1 - r22) + (1 - icc) * (1 - r12)/n)/
                                                             (J* (1 - r2m)) +
                                                             b^2 * (1 -r2m) / (p * (1 - p) * J)))
        })} else {
          pwr.sobel <- quote({
            1 - pnorm(qnorm(1 - sig.level / tside),
                      a * b / sqrt(a^2 * (icc * (1 - r22) + (1 - icc) * (1 - r12)/n)/
                                     (J* (1 - r2m)) +
                                     b^2 * (1 -r2m) / (p * (1 - p) * J)))
          })}
    }

    if (is.null(power)) {
      sobel.out <- list(power = eval(pwr.sobel))
    } else if (is.null(m) & is.null(J)) {
      if(cost.model){
        sobel.out <- list(m = stats::uniroot(function(m)
          eval(pwr.sobel) - power, mlim)$root);
        sobel.out <- c(sobel.out, list(J = sobel.out$m / (((1 - p) * (c1 * n + c2)
                                                           + p * (c1t * n + c2t)))))
      } else {
        sobel.out <- list(J = stats::uniroot(function(J)
          eval(pwr.sobel) - power, Jlim)$root)
        }
    }
    power.out <- list(funName = funName,
                      designType = designType,
                      par = par, sobel.out = sobel.out)
    return(power.out)
  }

  if (test == "joint" | test == "Joint" | test == "JOINT"){
    if (cost.model){
      # when cost.model is true for the Joint test
      if (two.tailed) {
        pwr.joint <- quote({
          J <- m / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t));
          lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
          lambda.b <- b * sqrt(J * (1 - r2m) /
                                 (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

          (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                  df = J - q.a - 2, lambda.a) +
             pt(qt(sig.level/tside, df = J - q.a - 2),
                df = J - q.a - 2, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b) +
               pt(qt(sig.level/tside, df = J - q.b - 3),
                  df = J - q.b - 3, lambda.b))
        })
      } else {
        pwr.joint <- quote({
          J <- m / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t));
          lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
          lambda.b <- b * sqrt(J * (1 - r2m) /
                                 (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

          (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                  df = J - q.a - 2, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b))
        })
      }
    } else {
      # when cost.model is not true for the Joint test
      if (two.tailed) {
        pwr.joint <- quote({
          lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
          lambda.b <- b * sqrt(J * (1 - r2m) /
                                 (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

          (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                  df = J - q.a - 2, lambda.a) +
             pt(qt(sig.level/tside, df = J - q.a - 2),
                df = J - q.a - 2, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b) +
               pt(qt(sig.level/tside, df = J - q.b - 3),
                  df = J - q.b - 3, lambda.b))
        })
      } else {
        pwr.joint <- quote({
          lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
          lambda.b <- b * sqrt(J * (1 - r2m) /
                                 (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

          (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                  df = J - q.a - 2, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b))
        })
      }
    }


    if (is.null(power)) {
      joint.out <- list(power = eval(pwr.joint))
    } else if (is.null(m) & is.null(J)) {
      if(cost.model){
        joint.out <- list(m = stats::uniroot(function(m)
          eval(pwr.joint) - power, mlim)$root);
        joint.out <- c(joint.out, list(J = joint.out$m / (((1 - p) * (c1 * n + c2)
                                                           + p * (c1t * n + c2t)))))
      } else {
        joint.out <- list(J = stats::uniroot(function(J)
          eval(pwr.joint) - power, Jlim)$root);
      }

    }
    power.out <- list(funName = funName,
                      designType = designType,
                      par = par, joint.out = joint.out)
    return(power.out)
  }

}




