#' Relative efficiency (RE) calculation
#'
#' @description Calculate the relative efficiency (RE) of two designs
#'
#' @param od returned object of first design (e.g., unconstrained optimal design)
#'     from function \code{\link{od.2}},
#'     \code{\link{od.3}}, or \code{\link{od.4}}
#' @param subod returned object of second design (e.g., constrained optimal design)
#'     from function \code{\link{od.2}},
#'     \code{\link{od.3}}, or \code{\link{od.4}}
#' @param verbose logical; print the values of relative efficiency if TRUE,
#'    otherwise not; default is TRUE.
#' @param rounded logical; round the values of \code{p}, \code{n}/\code{J}/\code{K}
#'     that are from functions to two decimal places and integer, respectively if TRUE,
#'     no rounding if FALSE; default is TRUE.
#' @return
#'     Relative efficiency value
#'
#' @export re
#'
#' @references
#'   (1) Shen, Z., & Kelcey, B. (2018, April). Optimal design of cluster
#'   randomized trials under condition- and unit-specific cost structures. Roundtable
#'   discussion to be presented at American Educational Research Association (AERA)
#'   annual conference; (2) Shen, Z., & Kelcey, B. (under review).
#'   Optimal sample allocation under unequal costs in cluster-randomized trials.
#'   Journal of Educational and Behavioral Statistics. (3) Shen, Z.(in progress).
#'   Using optimal sample allocation to
#'   improve design efficiency for multilevel randomized trials.
#'   (Unpublished doctoral dissertation). University of Cincinnati, Cincinnati, OH.
#'
#' @examples
#' # unconstrained optimal design of 2-level CRT #----------
#' myod1 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               varlim = c(0.01, 0.02))
#' # constrained optimal design with n = 20
#' myod2 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               n = 20, varlim = c(0.005, 0.025))
#' # relative efficiency (RE)
#' myre <- re(od = myod1, subod= myod2)
#' myre$out # RE = 0.88
#' # constrained optimal design with p = 0.5
#' myod2 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'              p = 0.5, varlim = c(0.005, 0.025))
#' # relative efficiency (RE)
#' myre <- re(od = myod1, subod= myod2)
#' myre$out # RE = 0.90
#'
#' # unconstrained optimal design of 3-level CRT #----------
#' myod1 <- od.3(icc2 = 0.2, icc3 = 0.1, r12 = 0.5, r22 = 0.5, r32 = 0.5,
#'              c1 = 1, c2 = 5, c3 = 25, c1t = 1, c2t = 50, c3t = 250,
#'              varlim = c(0.005, 0.025))
#' # constrained optimal design with J = 20
#' myod2 <- od.3(icc2 = 0.2, icc3 = 0.1, r12 = 0.5, r22 = 0.5, r32 = 0.5, J = 20,
#'              c1 = 1, c2 = 5, c3 = 25, c1t = 1, c2t = 50, c3t = 250,
#'              varlim = c(0, 0.025))
#' # relative efficiency (RE)
#' myre <- re(od = myod1, subod= myod2)
#' myre$out # RE = 0.53
#'
#' # unconstrained optimal design of 4-level CRT #---------
#' myod1 <- od.4(icc2 = 0.2, icc3 = 0.1, icc4 = 0.05, r12 = 0.5,
#'               r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'               c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'               c1t = 1, c2t = 50, c3t = 250, c4t = 2500,
#'               varlim = c(0, 0.01))
#' # constrained optimal design with p = 0.5
#' myod2 <- od.4(icc2 = 0.2, icc3 = 0.1, icc4 = 0.05, r12 = 0.5, p = 0.5,
#'               r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'               c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'               c1t = 1, c2t = 50, c3t = 250, c4t = 2500,
#'               varlim = c(0, 0.01))
#' # relative efficiency (RE)
#' myre <- re(od = myod1, subod= myod2)
#' myre$out # RE = 0.78
#'
re <- function(od, subod, rounded = TRUE, verbose = TRUE) {
  funName <- "re"
  if (od$funName==subod$funName) {
    designType <- od$designType
  } else {
    stop("re function can only compare relative efficiency (RE) of
         two studies with same design type")
  }
  if (designType == "two-level CRTs") {
    icc2 <- od$par$icc
    r12 <- od$par$r12
    r22 <- od$par$r22
    c1 <- od$par$c1
    c2 <- od$par$c2
    c1t <- od$par$c1t
    c2t <- od$par$c2t
    if (rounded) {
      no <- round(od$out$n, 0)
      po <- round(od$out$p, 2)
    } else {
      no <- od$out$n
      po <- od$out$p
    }
    if (
      od$par$icc != subod$par$icc ||
      od$par$r12 != subod$par$r12 ||
      od$par$r22 != subod$par$r22 ||
      od$par$c1 != subod$par$c1 ||
      od$par$c2 != subod$par$c2 ||
      od$par$c1t != subod$par$c1t ||
      od$par$c2t != subod$par$c2t
    )
    stop("Each of 'icc',  'r12', 'r22',
         'c1', 'c2', 'c1t', and 'c2t'
         must be equal in two designs")
    if (rounded) {
      n <- round(subod$out$n, 0)
      p <- round(subod$out$p, 2)
    } else {
      n <- subod$out$n
      p <- subod$out$p
    }
    RE <- ( no * icc2 * (1 - r22) + (1 - icc2) * (1 - r12)) /
      (n * icc2 * (1 - r22) + (1 - icc2) * (1 - r12)) *
      ((1 - po) * (c1 * no + c2)
       + po * (c1t * no + c2t)) /
      ((1 - p) * (c1 * n + c2)
       + p * (c1t * n + c2t)) *
      (p * (1 - p) * n) / (po * (1 - po) * no)
  } else if (designType == "three-level CRTs") {
    icc2 <- od$par$icc2
    icc3 <- od$par$icc3
    r12 <- od$par$r12
    r22 <- od$par$r22
    r32 <- od$par$r32
    c1 <- od$par$c1
    c2 <- od$par$c2
    c3 <- od$par$c3
    c1t <- od$par$c1t
    c2t <- od$par$c2t
    c3t <- od$par$c3t
    if (rounded) {
      no <- round(od$out$n, 0)
      Jo <- round(od$out$J, 0)
      po <- round(od$out$p, 2)
    } else {
      no <- od$out$n
      Jo <- od$out$J
      po <- od$out$p
    }
    if (
      od$par$icc2 != subod$par$icc2 ||
      od$par$icc3 != subod$par$icc3 ||
      od$par$r12 != subod$par$r12 ||
      od$par$r22 != subod$par$r22 ||
      od$par$r32 != subod$par$r32 ||
      od$par$c1 != subod$par$c1 ||
      od$par$c2 != subod$par$c2 ||
      od$par$c3 != subod$par$c3 ||
      od$par$c1t != subod$par$c1t ||
      od$par$c2t != subod$par$c2t ||
      od$par$c3t != subod$par$c3t
    )
    stop("Each of 'icc2', 'icc3', 'r12', 'r22', 'r32',
         'c1', 'c2', 'c3', 'c1t', 'c2t', and 'c3t'
         must be equal in two designs")
    if (rounded) {
      n <- round(subod$out$n, 0)
      J <- round(subod$out$J, 0)
      p <- round(subod$out$p, 2)
    } else {
      n <- subod$out$n
      J <- subod$out$J
      p <- subod$out$p
    }
    RE <- ( no * Jo * icc3 * (1 - r32)
           + no * icc2 * (1 - r22) + (1 - icc2 - icc3 ) * (1 - r12)) /
      ( n * J * icc3 * (1 - r32)
       + n * icc2 * (1 - r22) + (1 - icc2 - icc3) * (1 - r12)) *
      ((1 - po) * (c1 * no * Jo  + c2 * Jo  + c3)
       + po * (c1t * no * Jo + c2t * Jo + c3t )) /
      ((1 - p) * (c1 * n * J + c2 * J + c3)
       + p * (c1t * n * J + c2t * J + c3t)) *
      (p * (1 - p) * n * J) / (po * (1 - po) * no * Jo)
  } else if (designType == "four-level CRTs") {
    icc2 <- od$par$icc2
    icc3 <- od$par$icc3
    icc4 <- od$par$icc4
    r12 <- od$par$r12
    r22 <- od$par$r22
    r32 <- od$par$r32
    r42 <- od$par$r42
    c1 <- od$par$c1
    c2 <- od$par$c2
    c3 <- od$par$c3
    c4 <- od$par$c4
    c1t <- od$par$c1t
    c2t <- od$par$c2t
    c3t <- od$par$c3t
    c4t <- od$par$c4t
    if (rounded) {
      no <- round(od$out$n, 0)
      Jo <- round(od$out$J, 0)
      Ko <- round(od$out$K, 0)
      po <- round(od$out$p, 2)
    } else {
      no <- od$out$n
      Jo <- od$out$J
      Ko <- od$out$K
      po <- od$out$p
    }
    if (
      od$par$icc2 != subod$par$icc2 ||
      od$par$icc3 != subod$par$icc3 ||
      od$par$icc4 != subod$par$icc4 ||
      od$par$r12 != subod$par$r12 ||
      od$par$r22 != subod$par$r22 ||
      od$par$r32 != subod$par$r32 ||
      od$par$r42 != subod$par$r42 ||
      od$par$c1 != subod$par$c1 ||
      od$par$c2 != subod$par$c2 ||
      od$par$c3 != subod$par$c3 ||
      od$par$c4 != subod$par$c4 ||
      od$par$c1t != subod$par$c1t ||
      od$par$c2t != subod$par$c2t ||
      od$par$c3t != subod$par$c3t ||
      od$par$c4t != subod$par$c4t
    )
    stop("Each of 'icc2', 'icc3', 'icc4', 'r12', 'r22', 'r32', 'r42',
         'c1', 'c2', 'c3', 'c4', 'c1t', 'c2t', 'c3t', and 'c4t'
         must be equal in two designs")
    if (rounded) {
      n <- round(subod$out$n, 0)
      J <- round(subod$out$J, 0)
      K <- round(subod$out$K, 0)
      p <- round(subod$out$p, 2)
    } else {
      n <- subod$out$n
      J <- subod$out$J
      K <- subod$out$K
      p <- subod$out$p
    }
    RE <- (no * Jo * Ko * icc4 * (1 - r42) + no * Jo * icc3 * (1 - r32)
           + no * icc2 * (1 - r22) + (1 - icc2 - icc3 - icc4) * (1 - r12)) /
          (n * J * K * icc4 * (1 - r42) + n * J * icc3 * (1 - r32)
          + n * icc2 * (1 - r22) + (1 - icc2 - icc3 - icc4) * (1 - r12)) *
          ((1 - po) * (c1 * no * Jo * Ko + c2 * Jo * Ko + c3 * Ko + c4)
           + po * (c1t * no * Jo * Ko + c2t * Jo * Ko + c3t * Ko + c4t)) /
          ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4)
           + p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t)) *
          (p * (1 - p) * n * J * K) / (po * (1 - po) * no * Jo * Ko)
  }
    out <- RE
    if (verbose)  cat("The relative efficiency (RE) of the two ",
                      designType, " is ", RE, ".\n", sep = "")
    if (out > 1) {
      cat("===============================\n",
          "Please switch the objuects for 'od' and 'subod' to have 0 < RE <= 1",
          ".\n===============================\n", sep = "")
    }
    re.out <- list(funName = funName,
                        designType = designType,
                        par = list(od, subod), out = out)
    return(re.out)
  }




