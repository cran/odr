% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.2m.R
\name{power.2m}
\alias{power.2m}
\title{Budget and/or sample size, power, MDES calculation for two-level
MRTs detecting main effects}
\usage{
power.2m(
  cost.model = TRUE,
  expr = NULL,
  constraint = NULL,
  sig.level = 0.05,
  two.tailed = TRUE,
  d = NULL,
  power = NULL,
  m = NULL,
  n = NULL,
  J = NULL,
  p = NULL,
  icc = NULL,
  r12 = NULL,
  r22m = NULL,
  q = NULL,
  c1 = NULL,
  c2 = NULL,
  c1t = NULL,
  omega = NULL,
  dlim = NULL,
  powerlim = NULL,
  Jlim = NULL,
  mlim = NULL,
  rounded = TRUE
)
}
\arguments{
\item{cost.model}{Logical; power analyses accommodating costs and budget
(e.g., required budget for desired power, power/MDES under fixed budget)
if TRUE, otherwise conventional power analyses
(e.g., required sample size, power, or MDES calculation); default value is TRUE.}

\item{expr}{Returned objects from function \code{\link{od.2m}} or
\code{\link{od.2m.111m}}; default is NULL;
    if \code{expr} is specified, parameter values of \code{icc},
    \code{r12}, \code{r22m},
    \code{c1}, \code{c2},
    \code{c1t}, \code{p}, and \code{n}
    used or solved in function \code{\link{od.2m}} will
    be passed to current function;
    only the values of \code{p} and \code{n} that specified or solved in
    function \code{\link{od.2m}} can be overwritten
    if \code{constraint} is specified.}

\item{constraint}{Specify the constrained values of \code{p} and/or \code{n}
in list format to overwrite those from \code{expr}; default value is NULL.}

\item{sig.level}{Significance level or type I error rate, default value is 0.05.}

\item{two.tailed}{Logical; two-tailed tests if TRUE,
otherwise one-tailed tests; default value is TRUE.}

\item{d}{Effect size.}

\item{power}{Statistical power.}

\item{m}{Total budget.}

\item{n}{The level-1 sample size per level-2 unit.}

\item{J}{The number of sites.}

\item{p}{The proportion of level-1 units to be assigned to treatment.}

\item{icc}{The unconditional intraclass correlation coefficient (ICC) in population or in
each treatment condition.}

\item{r12}{The proportion of level-1 variance explained by covariates.}

\item{r22m}{The proportion of variance of site-specific treatment
effect explained by covariates.}

\item{q}{The number of covariates at level 2.}

\item{c1}{The cost of sampling one level-1 unit in control condition.}

\item{c2}{The cost of sampling one level-2 unit (site).}

\item{c1t}{The cost of sampling one level-1 unit in treatment condition.}

\item{omega}{The standardized variance of site-specific treatment effect.}

\item{dlim}{The range for solving the root of effect size (\code{d}) numerically,
default value is c(0, 5).}

\item{powerlim}{The range for solving the root of power (\code{power}) numerically,
default value is c(1e-10, 1 - 1e-10).}

\item{Jlim}{The range for searching the root of level-2 sample size (\code{J}) numerically,
default is c(4, 10e10).}

\item{mlim}{The range for searching the root of budget (\code{m}) numerically,
default is the costs sampling \code{Jlim} level-2 units
or c(4 * Jcost, 1e+10 * Jcost) with Jcost =
(1 - p) * c1 * n + p * c1t * n + c2.}

\item{rounded}{Logical; round the values of \code{p}, \code{n}/\code{J}/\code{K}
that are from functions \code{\link{od.4}}
to two decimal places and integer, respectively if TRUE,
otherwise no rounding; default value is TRUE.}
}
\value{
Required budget (and/or required level-2 sample size), s
    tatistical power, or MDES  depending on the specification of parameters.
    The function also returns the function name, design type,
    and parameters used in the calculation.
}
\description{
This function can calculate required budget for desired power,
    power or minimum detectable effect size (MDES) under fixed budget
    for two-level multisite randomized trials (MRTs).
    It also can perform conventional power analyses
    (e.g., required sample size, power, and MDES calculation).
}
\examples{
# Unconstrained optimal design #---------
  myod1 <- od.2m(icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
              c1 = 1, c2 = 10, c1t = 10,
              varlim = c(0, 0.005))
  myod1$out # n = 23.5, p = 0.24

# ------- Power analyses by default considering costs and budget -------
# Required budget and sample size
  mym.1 <- power.2m(expr = myod1, d = 0.2, q = 1, power = 0.8)
  mym.1$out  # m = 1882, J = 22.8
  # mym.1$par  # parameters and their values used for the function
# Or, equivalently, specify every argument in the function
  mym.1 <- power.2m(d = 0.2, power = 0.8, q = 1,
                 icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
                 c1 = 1, c2 = 10, c1t = 10,
                 n = 23.5, p = 0.24)
# Required budget and sample size with constrained p
  mym.2 <- power.2m(expr = myod1, d = 0.2, q = 1, power = 0.8,
                 constraint = list(p = 0.5))
  mym.2$out  # m = 2424, J = 18
# Required budget and sample size with constrained p and n
  mym.3 <- power.2m(expr = myod1, d = 0.2, q = 1, power = 0.8,
                 constraint = list(p = 0.5, n = 5))
  mym.3$out  # m = 2502, J = 66.7

# Power calculation
  mypower <- power.2m(expr = myod1, q = 1, d = 0.2, m = 1882)
  mypower$out  # power = 0.80
# Power calculation under constrained p (p = 0.5)
  mypower.1 <- power.2m(expr = myod1, q = 1, d = 0.2, m = 1882,
                 constraint = list(p = 0.5))
  mypower.1$out  # power = 0.68

# MDES calculation
  mymdes <- power.2m(expr = myod1, q = 1, power = 0.80, m = 1882)
  mymdes$out  # d = 0.20


# ------- Conventional power analyses with cost.model = FALSE-------
# Required sample size
  myJ <- power.2m(cost.model = FALSE, expr = myod1, d = 0.2,
  q = 1, power = 0.8)
  myJ$out  # J = 22.8
  # myL$par  # parameters and their values used for the function
# Or, equivalently, specify every argument in the function
  myJ <- power.2m(cost.model = FALSE, d = 0.2, power = 0.8, q = 1,
                 icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
                 c1 = 1, c2 = 10, c1t = 10,
                 n = 23.5, p = 0.24)

# Power calculation
  mypower1 <- power.2m(cost.model = FALSE, expr = myod1, J = 22.8,
  d = 0.2, q = 1)
  mypower1$out  # power = 0.80

# MDES calculation
  mymdes1 <- power.2m(cost.model = FALSE, expr = myod1, J = 22.8,
  power = 0.8, q = 1)
  mymdes1$out  # d = 0.20

}
\references{
Shen, Z., & Kelcey, B. (2022). Optimal sample
   allocation in multisite randomized trials.
   The Journal of Experimental Education, 90(3), 693-711.
   <https://doi.org/10.1080/00220973.2020.1830361
}
