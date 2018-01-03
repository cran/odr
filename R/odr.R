#' Optimal Design and Statistical Power of
#'     Multilevel Randomized Trials
#'
#' This package is to help researchers design cost-efficient and
#'     well-powered multilevel randomized trials by solving optimal sample allocation
#'     and explicitly accommodating costs and budget in power analyses.
#'
#' This package has two categorical functions for each type of
#'     multilevel randomized trials.
#'     The two categorical functions are
#'     'od' and 'power'. The 'od' function can
#'     calculate the optimal sample allocation with and without constraints for multilevel
#'     randomized trials under fixed budget and cost structure that minimizes
#'     the variance of a treatment effect. The optimal design parameters
#'     include the optimal sample sizes at each level and
#'     proportion of units assigned to treatment.
#'     The 'power' function by default can calculate required budget for desired
#'     power, minimum detectable effect size (MDES) or power under fixed budget.
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, MDES calculation).
#'     There is an additional function 're' to compare
#'     the relative efficiency to two studies.
#'
#' @author Zuchao "William" Shen
#'
#' Maintainer: Zuchao "William" Shen \href{mailto: shenzo@mail.uc.edu}{shenzo@mail.uc.edu}
#'
#' @docType package
#' @name odr
NULL
