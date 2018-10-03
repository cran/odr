#' Optimal Design and Statistical Power of Multilevel Randomized Trials
#'
#' This package is to help researchers design cost-efficient multilevel randomized trials
#'     with adequate statistical precision by (1) solving optimal sample allocations,
#'     (2) comparing design precision and efficiency between different sample allocations,
#'     and (3) explicitly accommodating costs and budget in power analyses.
#'
#' The package covers seven types of trials with continuous outcomes and these trials are
#'     individual randomized controlled trials (RCTs), two-,
#'     three-, and four-level cluster-randomized trials (CRTs), and
#'     two-, three-, and four-level multisite randomized trials (MRTs).
#'     There are two categorical functions for each type of
#'     trial and a uniform function for all types of trials.
#'     The two categorical functions are
#'     'od' and 'power'. The 'od' function can
#'     calculate the optimal sample allocation with and without a constraint(s) for
#'     each type of trial. The optimal design parameters in 'od' function
#'     include the sample sizes at each level and
#'     proportion of units to be assigned to the treatment condition.
#'     The 'power' function by default can calculate required budget
#'     (and required sample size) for desired
#'     power, minimum detectable effect size (MDES) under a fixed budget,
#'     statistical power under a fixed budget.
#'     The 'power' function also can perform conventional power analyses
#'     (e.g., required sample size, power, MDES calculation).
#'     The uniform function 're' (or 'rpe') is to compare
#'     the relative (precision and) efficiency between two designs with different sample allocations.
#'
#' @author Zuchao Shen, Ben Kelcey
#'
#' Maintainer: Zuchao Shen \href{mailto: shenzo@mail.uc.edu}{shenzo@mail.uc.edu}
#'
#' @docType package
"_PACKAGE"
