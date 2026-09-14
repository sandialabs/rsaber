#' Sequential lot sampling (SLS) evaluation functions
#'
#' These functions implement lot-by-lot evaluation rules used by \pkg{rsaber}'s
#' sequential lot sampling (SLS) simulation utilities. Each evaluator uses a
#' common argument set so that it can be selected programmatically (e.g., via a
#' control list such as `sim_control`).
#'
#'
#' @param Lot (numeric) Vector of 0's and 1's to denote good and defective units.
#' @param preq (numeric) Requirement on defect rate.
#' @param alpha (numeric) Tail area for credibility level (`alpha = 0.05` means 95% credibility).
#' @param kvec (numeric) vector with the sequence of number defective in each lot.
#' @param nvec (numeric) vector with the sequence of number sampled from each lot.
#' @param idx (numeric) The index of the current lot.
#' @param ab (numeric) vector of length 2 with the most recent posterior.
#' @param control An object of class `saber_control` containing SABER model
#'   settings.
#' @param sls_control An object of class `sls_control` containing SLS
#'   sampling/evaluation settings. 
#'
#' @return Returns a numeric vector with:
#' 
#' - Number of defects sampled from the current lot.
#' - Number of units sampled from the current lot.
#' - First parameter of the posterior.
#' - Second parameter of the posterior.
#' - Number of defective units accepted.
#' - Number of units accepted.
#' - Whether the lot was accepted (1) or rejected (0).
#'
#' @details
#' An *evaluation function* is called once per lot during a simulation. It
#' updates any required posterior state and returns a numeric vector of summary
#' outputs for the current lot (e.g., sample size, defects observed, posterior
#' parameters, and accept/reject indicators).
#'
#' To use a custom evaluator, supply a function with the same argument list and
#' expected return structure, and pass it to [sls_control()] through the
#' `eval_function` argument.
#' 
#' @note
#' Not all of the arguments are necessarily used. They are in place to ensure a 
#' consistent set of arguments for evaluation functions. 
#' If a user writes their own evaluation function, it may be passed by adding it
#' to `sim_control` with the name "eval_method".
#'
#' @seealso
#' [sls_control()] for specifying the SLS evaluation method and policy settings;
#' [sls_sim_case()] for running a single SLS simulation case;
#' [saber_control()] for SABER model settings.
#'
#' 
#' 
#' @examples
#' control <- saber_control(
#'   prior = c(1, 1),
#'   lambda = 0.95,
#'   forget = "exponential",
#'   frame = "total"
#' )
#'
#' sls_ctrl <- sls_control(
#'   x0 = 0,
#'   nn_min = 5,
#'   eval_function = "saber_fixed"
#' )
#'
#' Lot <- c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
#'
#' eval_sls_saber_fixed(
#'   Lot = Lot,
#'   preq = 0.10,
#'   alpha = 0.05,
#'   kvec = numeric(0),
#'   nvec = numeric(0),
#'   idx = 1,
#'   ab = control[["prior"]],
#'   control = control,
#'   sls_control = sls_ctrl
#' )
#' 
#' 
#' 
#' @name sls_evaluation
#' @aliases sls_evaluation
#' @keywords internal
NULL


#' @describeIn sls_evaluation Applies the adaptive SABER SLS sampling rule.
#' 
#' @export
#'  
eval_sls_saber_adapt <- function( Lot, preq, alpha, 
                                  kvec, nvec, idx, ab, 
                                  control, sls_control ){
  
  # Extract values from sim_control to improve readability
  x0 <- sls_control[["x0"]]
  Ni <- length(Lot)
  
  # Account for unspecified values
  kappa <- sls_control[["kappa"]]
  if(is.null(kappa) ){
    kappa <- 0.5 
  }
  
  nn_min <- sls_control[["nn_min"]]
  if (is.null(nn_min)) {
    if( alpha == 0.5 ){
      stop_message <- c(
        "Cannot automatically determine a rational minimum sample size nn_min.\n",
        "Options to resolve this include:\n",
        "- Set nn_min in sim_control\n",
        "- Set kappa  in sim_control\n",
        "- Set alpha < 0.5\n"
      )
    }
    nn_min <- rsaber::prop_n_bayes(
      x     = x0,
      N     = Ni,
      preq  = preq,
      alpha = kappa,
      prior = control[["prior"]],
      frame = control[["frame"]], 
      verbose = FALSE
    )[["n"]]
  }
  
  # Get sample size
  nn_saber <- rsaber::prop_n_bayes(
    x = x0, 
    N = Ni, preq = preq, alpha = alpha, 
    prior = ab,
    frame = control[["frame"]], 
    verbose = FALSE
  )[["n"]]
  # Rescue for if no sample size is sufficient
  if( is.na(nn_saber) ){ nn_saber <- Ni }
  
  # Sample size and number of defectives sampled
  nn <- max(nn_min, nn_saber)
  ki <- sum( Lot[1:nn] )
  
  # Evaluate lot acceptances
  if (ki > x0) {
    dAcc <- 0
    nAcc <- 0
    lAcc <- 0
  } else {
    dAcc <- sum(Lot[-(1:nn)])
    nAcc <- length(Lot[-(1:nn)])
    lAcc <- 1
  }
  
  # Update kvec and nvec, then get the posterior
  kvec <- c( kvec, ki )
  nvec <- c( nvec, nn )
  sab <- saber_posterior(
    x = kvec, 
    n = nvec,
    tidx = seq_len(idx),
    prior  = control[["prior"]],
    lambda = control[["lambda"]],
    forget = control[["forget"]]
  )
  
  # Ship the results
  return(
    c(  
      ki, nn, sab[1], sab[2],
      dAcc, nAcc, lAcc
    )
  )
}



#' @describeIn sls_evaluation Applies SABER SLS model, but with fixed sample 
#' size rather than adaptive.
#'
#' @export
#' 
eval_sls_saber_fixed <- function( Lot, preq, alpha, 
                                  kvec, nvec, idx, ab, 
                                  control, sls_control ){
  
  # Extract values from sim_control to improve readability
  x0 <- sls_control[["x0"]]
  
  # Sample size and number of defectives sampled
  nn <- sls_control[["nn_min"]]
  if (is.null(nn)) {
    stop("For fixed sample size method, sim_control must contain nn_min.")
  }
  ki <- sum( Lot[1:nn] )
  
  # Evaluate lot acceptances
  if (ki > x0) {
    dAcc <- 0
    nAcc <- 0
    lAcc <- 0
  } else {
    dAcc <- sum(Lot[-(1:nn)])
    nAcc <- length(Lot[-(1:nn)])
    lAcc <- 1
  }
  
  # Update kvec and nvec, then get the posterior
  kvec <- c( kvec, ki )
  nvec <- c( nvec, nn )
  sab <- saber_posterior(
    x = kvec, n = nvec,
    tidx = seq_len(idx),
    prior  = control[["prior"]],
    lambda = control[["lambda"]],
    forget = control[["forget"]]
  )
  
  # Ship the results
  return(
    c(  
      ki, nn, sab[1], sab[2],
      dAcc, nAcc, lAcc
    )
  )
}


#' @describeIn sls_evaluation Applies an independent Bayesian sampling rule
#'   without accumulating information across lots.
#'
#' @export
#'  
eval_sls_bayesian_indp <- function( Lot, preq, alpha, 
                                    kvec, nvec, idx, ab, 
                                    control, sls_control ){
  
  if( is.null(sls_control[["nn_min"]]) | is.null(control[["prior"]]) ) {
    stop("For independent Bayesian sampling, sim_control must contain: nn_min, prior")
  }
  
  # Extract values from sim_control to improve readability
  x0     <- sls_control[["x0"]]
  prior  <- control[["prior"]]
  
  # Sample size and number of defectives sampled
  nn <- sls_control[["nn_min"]]
  ki <- sum( Lot[1:nn] )
  
  if (ki > x0) {
    dAcc <- 0
    nAcc <- 0
    lAcc <- 0
  } else {
    dAcc <- sum(Lot[-(1:nn)])
    nAcc <- length(Lot[-(1:nn)])
    lAcc <- 1
  }
  
  # Get the posterior
  sab <- c( prior[1] + ki, prior[2] + nn - ki )
  
  # Ship the results
  return(
    c(  
      ki, nn, sab[1], sab[2],
      dAcc, nAcc, lAcc
    )
  )
}


#' @describeIn sls_evaluation Applies an independent Bayesian skip-lot sampling
#'   rule.
#'
#' @export
#'  
eval_sls_skiplot_indp <- function( Lot, preq, alpha, 
                                   kvec, nvec, idx, ab, 
                                   control, sls_control ){
  # Extract values from sim_control to improve readability
  x0        <- sls_control[["x0"]]
  nn_min    <- sls_control[["nn_min"]]
  prior     <- control[["prior"]]
  max_skip  <- sls_control[["skip_max"]]
  skip_prob <- sls_control[["skip_prob"]]
  
  if( is.null(nn_min) | is.null(max_skip) | is.null(skip_prob) ) {
    stop("For skip-lot sampling, sim_control must contain: nn_min, max_skip, skip_prob.")
  }
  
  # Evaluate skipping lot & find the sample size (nn = 0 if skipped)
  nn <- nn_min
  if( idx > max_skip ){ # Cannot skip if less than the largest number of skips permitted
    idx_skip <- seq( idx-1, idx-max_skip, -1 )
    k_last   <- sum( kvec[ idx_skip ] )
    n_last   <-      nvec[ idx_skip ]
    skip_lot <- runif( 1, 0, 1 ) <= skip_prob
    # Logic should be:
    # - There WAS a sample in the last [max_skip]
    # - There were NOT defects in the last [max_skip]
    # - The randomizer DID choose to skip a lot
    if( any(n_last>0) & (k_last==0) & (skip_lot==1) ){
      nn <- 0
    }
  }
  
  # Evaluate lot acceptances
  if( nn > 0 ){
    # Number of defectives sampled
    ki <- sum( Lot[1:nn] )
    if (ki > x0) {
      dAcc <- 0
      nAcc <- 0
      lAcc <- 0
    } else {
      dAcc <- sum(Lot[-(1:nn)])
      nAcc <- length(Lot[-(1:nn)])
      lAcc <- 1
    }
    sab <- c( prior[1] + ki, prior[2] + nn - ki )
  } else{
    # Skipped lot, so no defectives sampled
    ki <- 0
    dAcc <- sum(Lot)
    nAcc <- length(Lot)
    lAcc <- 1
    sab <- prior
  }
  
  # Ship the results
  return(
    c(  
      ki, nn, sab[1], sab[2],
      dAcc, nAcc, lAcc
    )
  )
}



#' @describeIn sls_evaluation Applies a SABER-based skip-lot sampling rule.
#'
#' @export
#'  
eval_sls_skiplot_saber <- function( Lot, preq, alpha, 
                                    kvec, nvec, idx, ab, 
                                    control, sls_control ){
  
  # Extract values from sim_control to improve readability
  x0        <- sls_control[["x0"]]
  nn_min    <- sls_control[["nn_min"]]
  prior     <- control[["prior"]]
  max_skip  <- sls_control[["skip_max"]]
  skip_prob <- sls_control[["skip_prob"]]
  
  if( is.null(nn_min) | is.null(max_skip) | is.null(skip_prob) ) {
    stop("For skip-lot sampling, sim_control must contain: nn_min, max_skip, skip_prob.")
  }
  
  # Evaluate skipping lot & find the sample size (nn = 0 if skipped)
  nn <- nn_min
  if( idx > max_skip ){ # Cannot skip if less than the largest number of skips permitted
    idx_skip <- seq( idx-1, idx-max_skip, -1 )
    k_last   <- sum( kvec[ idx_skip ] )
    n_last   <-      nvec[ idx_skip ]
    skip_lot <- runif( 1, 0, 1 ) <= skip_prob
    # Logic should be:
    # - There WAS a sample in the last [max_skip]
    # - There were NOT defects in the last [max_skip]
    # - The randomizer DID choose to skip a lot
    if( any(n_last>0) & (k_last==0) & (skip_lot==1) ){
      nn <- 0
    }
  }
  
  # Evaluate lot acceptances
  if( nn > 0 ){
    # Number of defectives sampled
    ki <- sum( Lot[1:nn] )
    if (ki > x0) {
      dAcc <- 0
      nAcc <- 0
      lAcc <- 0
    } else {
      dAcc <- sum(Lot[-(1:nn)])
      nAcc <- length(Lot[-(1:nn)])
      lAcc <- 1
    }
  } else{
    # Skipped lot, so no defectives sampled
    ki <- 0
    dAcc <- sum(Lot)
    nAcc <- length(Lot)
    lAcc <- 1
  }
  
  # Update kvec and nvec, then get the posterior
  kvec <- c( kvec, ki )
  nvec <- c( nvec, nn )
  
  sab <- saber_posterior(
    x = kvec, n = nvec,
    tidx = seq_len(idx),
    prior  = control[["prior"]],
    lambda = control[["lambda"]],
    forget = control[["forget"]]
  )
  
  # Ship the results
  return(
    c(  
      ki, nn, sab[1], sab[2],
      dAcc, nAcc, lAcc
    )
  )
}

