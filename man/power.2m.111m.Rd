% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.2m.111m.R
\name{power.2m.111m}
\alias{power.2m.111m}
\title{Budget and/or sample size, power, MDES calculation for two-level
MRTs detecting moderation effects with moderators at level 1}
\usage{
power.2m.111m(
  cost.model = TRUE,
  expr = NULL,
  constraint = NULL,
  sig.level = 0.05,
  two.tailed = TRUE,
  omega = NULL,
  gamma = NULL,
  power = NULL,
  m = NULL,
  n = NULL,
  J = NULL,
  p = NULL,
  icc = NULL,
  r12 = NULL,
  q = 1,
  c1 = NULL,
  c2 = NULL,
  c1t = NULL,
  gammalim = c(0, 5),
  powerlim = c(1e-10, 1 - 1e-10),
  Jlim = c(6, 1e+10),
  mlim = NULL,
  rounded = TRUE,
  binary = TRUE,
  Q = 0.5,
  random.slope = TRUE
)
}
\arguments{
\item{cost.model}{Logical; power analyses accommodating costs and budget
(e.g., required budget for desired power, power/MDES under fixed budget)
if TRUE, otherwise conventional power analyses
(e.g., required sample size, power, or MDES calculation); default value is TRUE.}

\item{expr}{Returned objects from function \code{\link{od.2m}}; default is NULL;
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

\item{omega}{The treatment-by-site variance of the outcome.}

\item{gamma}{The standardized moderated treatment effect.}

\item{power}{Statistical power.}

\item{m}{Total budget, default is the total costs of sampling 60
sites.}

\item{n}{The level-1 sample size per level-2 unit.}

\item{J}{The number of sites.}

\item{p}{The proportion of level-1 units to be assigned to treatment.}

\item{icc}{The unconditional intraclass correlation coefficient (ICC) in population or in
each treatment condition.}

\item{r12}{The proportion of level-1 variance explained by covariates.}

\item{q}{The number of additional covariates at level 1 beyond the treatment
indicator, covariate, and the interaction between the moderator and
the treatment. Default is 1.}

\item{c1}{The cost of sampling one level-1 unit in control condition.}

\item{c2}{The cost of sampling one level-2 unit (site).}

\item{c1t}{The cost of sampling one level-1 unit in treatment condition.}

\item{gammalim}{The range for searching the root of standardized
moderation effect (gamma). Default is c(0, 5).}

\item{powerlim}{The range for solving the root of power (\code{power}) numerically,
default value is c(1e-10, 1 - 1e-10).}

\item{Jlim}{The range for solving the root of level-2 sample size
(\code{J}) numerically, default value is c(6, 1e+10).}

\item{mlim}{The range for searching the root of budget (\code{m}) numerically,
default is the costs sampling \code{Jlim} level-2 units
or c(4 * Jcost, 1e+10 * Jcost) with Jcost =
(1 - p) * c1 * n + p * c1t * n + c2.}

\item{rounded}{Logical; round the values of \code{p}, \code{n}/\code{J}/\code{K}
that are from functions \code{\link{od.4}}
to two decimal places and integer, respectively if TRUE,
otherwise no rounding; default value is TRUE.}

\item{binary}{Logical; binary moderator if TURE and continuous moderator if
FALSE. Default is TRUE.}

\item{Q}{The proportion of units in one group for the binary moderator.
Default is 0.5.}

\item{random.slope}{Logical, the model is a random slope one if TURE. Default
is TRUE.}
}
\value{
Required budget (and/or required level-2 sample size), statistical power, or MDES
    depending on the specification of parameters.
    The function also returns the function name, design type,
    and parameters used in the calculation.
}
\description{
This function can calculate required budget for desired power,
    power or minimum detectable effect size (MDES) under fixed budget
    for two-level multisite randomized trials (MRTs)
    detecting moderation effects with moderators at level 1.
    It also can perform conventional power analyses
    (e.g., required sample size, power, and MDES calculation).
}
\examples{
myod <- od.2m.111m(icc = .2, r12 = .5, r22m = .5,
                   c1 = 10, c1t = 100, c2 = 50, omega = .01, gamma = 0.1)
mypowwer <- power.2m.111m(expr = myod, gamma = .1, power = .8)
}
