


#' @title Simulate defective-count contributions
#'
#' @description
#' Simulates defective-count contributions from a fitted SABER object. The
#' returned Monte Carlo matrix contains one row for each unit, lot, or process
#' element in the fitted object and one column for each Monte Carlo replicate.
#' Cumulative summaries can be computed from the returned matrix using
#' [sim_ndefect_summary()].
#'
#' @param saber_fit A fitted SABER object, such as an object of class
#'   `saber_cps` or `saber_sls`.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' This is an S3 generic. Methods are provided for SABER fit subclasses such as
#' `saber_cps` and `saber_sls`.
#'
#' The returned matrix uses the common `Xsim` convention: rows correspond to
#' sequential units or lots, and columns correspond to Monte Carlo replicates.
#' Each matrix entry is the simulated defective-count contribution for that row
#' under the selected frame of reference.
#'
#' @return
#' A numeric matrix of simulated defective-count contributions. Rows correspond
#' to units, lots, or process elements in `saber_fit`; columns correspond to
#' Monte Carlo replicates.
#' 
#' @examples
#' \dontrun{
#' 
#' # Continuous production sampling (cps)
#' control <- saber_control(
#'   prior = c(1, 1),
#'   lambda = 0.95,
#'   forget = "exponential",
#'   frame = "total"
#' )
#'
#' x <- c(0, 1, 0, NA, NA)
#' tidx <- seq_along(x)
#'
#' fit <- cps_bound(
#'   x = x,
#'   tidx = tidx,
#'   alpha = 0.05,
#'   control = control
#' )
#'
#' Xsim <- sim_ndefect(fit, nSim = 1000)
#' dim(Xsim)
#'
#' sim_ndefect_summary(Xsim, alpha = 0.05)
#' 
#' # Series-of-lots (SLS)
#' fit_sls <- sls_bound(
#'   x = c(0, 1, 0, 2),
#'   n = c(10, 10, 10, 10),
#'   N = c(100, 100, 100, 100),
#'   tidx = seq_len(4),
#'   alpha = 0.05,
#'   control = control
#' )
#'
#' Xsim_sls <- sim_ndefect(fit_sls, nSim = 1000)
#' dim(Xsim_sls)
#'
#' sim_ndefect_summary(Xsim_sls, alpha = 0.05)
#' 
#' }
#' 
#' @seealso [sim_ndefect_summary()], [sim_ndefect_quantiles()],
#'   [cps_bound()], [sls_bound()]
#'
#' 
#' @export
#' 
sim_ndefect <- function(saber_fit, ...) {
  UseMethod("sim_ndefect")
}


#' @describeIn sim_ndefect Default method for unsupported object classes.
#'
#' @export
#'  
sim_ndefect.default <- function(saber_fit, ...) {
  stop(
    "No 'sim_ndefect' method is available for objects of class: ",
    paste(class(saber_fit), collapse = ", "),
    call. = FALSE
  )
}



#' @describeIn sim_ndefect Method for continuous production sampling fits.
#'
#' @param nSim Number of Monte Carlo replicates.
#' @param frame Frame of reference for simulation. If `NULL`, taken from
#'   `saber_fit`.
#' @param Xsim_prefix Optional matrix containing already-simulated prefix rows.
#'   Used when extending simulations sequentially.
#' @param seed Optional random-number seed.
#' 
#' 
#' 
#' @export
#' 
sim_ndefect.saber_cps <- function(saber_fit, nSim = 1000,
                                  frame = NULL,
                                  Xsim_prefix = NULL,
                                  seed = NULL, ...) {
  
  validate_saber_fit(saber_fit, subclass = "saber_cps")
  
  if (is.null(frame)) {
    frame <- get_saber_frame(saber_fit)
    if (is.null(frame)) {
      frame <- "total"
    }
  }
  frame <- match_saber_frame(frame)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  nUnits <- nrow(saber_fit)
  Xsim <- matrix(0, nrow = nUnits, ncol = nSim)
  
  x <- saber_fit[["x"]]
  sampled <- !is.na(x)
  avec <- saber_fit[["a"]]
  bvec <- saber_fit[["b"]]
  
  if (!is.null(Xsim_prefix)) {
    if (!is.matrix(Xsim_prefix)) {
      stop("Argument 'Xsim_prefix' must be a matrix.", call. = FALSE)
    }
    
    if (ncol(Xsim_prefix) != nSim) {
      stop("Argument 'Xsim_prefix' must have nSim columns.", call. = FALSE)
    }
    
    if (nrow(Xsim_prefix) > nUnits) {
      stop("Argument 'Xsim_prefix' cannot have more rows than 'saber_fit'.", call. = FALSE)
    }
    
    Xsim[seq_len(nrow(Xsim_prefix)), ] <- Xsim_prefix
    ii_start <- nrow(Xsim_prefix) + 1
  } else {
    ii_start <- 1
  }
  
  if (ii_start <= nUnits) {
    for (ii in ii_start:nUnits) {
      if (sampled[ii]) {
        Xsim[ii, ] <- rep(0, nSim)
      } else {
        Xsim[ii, ] <- rsaber::rbetabinomial(
          n = nSim,
          Nn = 1,
          a = avec[ii],
          b = bvec[ii]
        )
      }
    }
  }
  
  if (frame == "total") {
    idx_sampled <- which(sampled)
    idx_sampled <- idx_sampled[idx_sampled >= ii_start]
    Xsim[idx_sampled, ] <- Xsim[idx_sampled, ] + x[idx_sampled]
  }
  
  return(Xsim)
}


#' @describeIn sim_ndefect Method for series-of-lots fits.
#'
#' @export
sim_ndefect.saber_sls <- function(saber_fit, nSim = 1000,
                                  frame = NULL,
                                  seed = NULL, ...) {
  
  validate_saber_fit(saber_fit, subclass = "saber_sls")
  
  if (is.null(frame)) {
    frame <- get_saber_frame(saber_fit)
    if (is.null(frame)) {
      frame <- "total"
    }
  }
  frame <- match_saber_frame(frame)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  Nvec <- saber_fit[["N"]]
  nvec <- saber_fit[["n"]]
  xvec <- saber_fit[["x"]]
  avec <- saber_fit[["a"]]
  bvec <- saber_fit[["b"]]
  
  nLots <- nrow(saber_fit)
  nu <- Nvec - nvec
  
  Xsim <- matrix(0, nrow = nLots, ncol = nSim)
  
  for (ii in seq_len(nLots)) {
    Xsim[ii, ] <- rsaber::rbetabinomial(
      n = nSim,
      Nn = nu[ii],
      a = avec[ii],
      b = bvec[ii]
    )
  }
  
  if (frame == "total") {
    Xsim <- Xsim + xvec
  }
  
  return(Xsim)
}









# Helper functions for sim_ndefect function(s)


#' @title Quantiles for Monte Carlo simulations of defectives
#'
#' @description
#' Computes row-wise quantiles of the cumulative defective-count simulations
#' represented by an `Xsim` matrix. Rows of `Xsim` are interpreted as sequential
#' units, lots, or process elements, and columns are interpreted as Monte Carlo
#' replicates.
#'
#' @param Xsim Numeric matrix of simulated defective-count contributions, such
#'   as the result of [sim_ndefect()].
#' @param probs Numeric vector of probabilities for the desired quantiles.
#'
#' @return A numeric matrix with one row for each row of `Xsim` and one column
#'   for each value in `probs`.
#' 
#' @seealso [sim_ndefect()], [sim_ndefect_summary()]
#' 
#' @examples
#' \dontrun{
#' control <- saber_control(
#'   prior = c(1, 1),
#'   lambda = 0.95,
#'   forget = "exponential",
#'   frame = "total"
#' )
#'
#' x <- c(0, 1, 0, NA, NA)
#' tidx <- seq_along(x)
#'
#' fit <- cps_bound(
#'   x = x,
#'   tidx = tidx,
#'   alpha = 0.05,
#'   control = control
#' )
#'
#' Xsim <- sim_ndefect(fit, nSim = 1000)
#'
#' sim_ndefect_quantiles(
#'   Xsim = Xsim,
#'   probs = c(0.05, 0.50, 0.95)
#' )
#' } 
#' 
#'
#' @export
#' 
sim_ndefect_quantiles <- function(Xsim, probs) {
  Xsim_csum <- apply(Xsim, 2, FUN=cumsum)
  Xsim_qnt  <- t(apply(Xsim_csum, 1, FUN=quantile, probs=probs))
  return(Xsim_qnt)
}








#' @title Summarize Monte Carlo simulations of defectives
#'
#' @description
#' Computes median and interval summaries from cumulative defective-count
#' simulations represented by an `Xsim` matrix. Rows of `Xsim` are interpreted
#' as sequential units, lots, or process elements, and columns are interpreted
#' as Monte Carlo replicates.
#'
#' @param Xsim Numeric matrix of simulated defective-count contributions, such
#'   as the result of [sim_ndefect()].
#' @param alpha Numeric tail probability used for lower and upper interval
#'   summaries. The returned interval uses quantiles `alpha` and `1 - alpha`.
#'
#' @return A numeric matrix with one row for each row of `Xsim` and columns
#'   `Xest`, `Xlcl`, and `Xucl`.
#'
#' @seealso [sim_ndefect()], [sim_ndefect_quantiles()]
#' 
#' 
#' @examples
#' \dontrun{
#' control <- saber_control(
#'   prior = c(1, 1),
#'   lambda = 0.95,
#'   forget = "exponential",
#'   frame = "total"
#' )
#'
#' x <- c(0, 1, 0, NA, NA)
#' tidx <- seq_along(x)
#'
#' fit <- cps_bound(
#'   x = x,
#'   tidx = tidx,
#'   alpha = 0.05,
#'   control = control
#' )
#'
#' Xsim <- sim_ndefect(fit, nSim = 1000)
#'
#' sim_ndefect_quantiles(
#'   Xsim = Xsim,
#'   probs = c(0.05, 0.50, 0.95)
#' )
#' } 
#' 
#' 
#' 
#' @export
sim_ndefect_summary <- function(Xsim, alpha) {
  probs <- c(0.5, alpha[1], 1 - alpha[1])
  Xsummary <- sim_ndefect_quantiles(
    Xsim = Xsim,
    probs = probs
  )
  colnames(Xsummary) <- c("Xest", "Xlcl", "Xucl")
  return(Xsummary)
}
