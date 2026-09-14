



#' @title SABER for Sequential Lot Sampling
#' 
#' @description Perform estimation and inference using the SABER method. Applies 
#' to continuous production (one-at-a-time) sampling.
#' 
#' @param x Numeric; vector denoting whether unit is defect (x=1) or not (x=0)
#' @param tidx Numeric; the vector containing the time point
#' @param alpha (numeric) Setting the credibility level.
#' @param control object of class [saber_control()]
#' @param saber_fit Object of class `saber_cps`
#' 
#' @details
#' Argument `control` should contain the following elements:
#' 
#' - `prior`: A vector of length 2 denoting the starting prior.
#' - `lambda`: The parameter for the forgetting function.
#' - `forget`: Character string or function. 
#' - `frame`: The frame of reference; "total", "unsampled", or "accepted".
#' 
#' 
#' @return Returns an object of class `saber_cps` 
#' 
#' @examples
#' 
#' control <- saber_control(
#'   prior = c(1, 1),
#'   lambda = 0.95,
#'   forget = "exponential",
#'   frame = "total"
#' )
#'
#' x <- c(0, 1, 0, NA, NA, 0, NA)
#' tidx <- seq_along(x)
#'
#' fit <- cps_bound(
#'   x = x,
#'   tidx = tidx,
#'   alpha = 0.05,
#'   control = control
#' )
#'
#' fit
#' 
#' @export
#' 
cps_bound <- function( x, tidx, alpha=NULL, control = NULL, saber_fit=NULL ){
  
  if (is.null(control)) {
    if (!is.null(saber_fit)) {
      validate_saber_fit(saber_fit, subclass = "saber_cps")
      control <- get_saber_control(saber_fit)
    }
  }
  control <- validate_saber_control(control)
  
  if (is.null(alpha)) {
    if (!is.null(saber_fit)) {
      alpha <- get_saber_alpha(saber_fit)
    }
    if (is.null(alpha)) {
      stop(
        "Argument 'alpha' must be supplied or recoverable from 'saber_fit'.",
        call. = FALSE
      )
    }
  }
  
  prior  <- control[["prior"]]
  lambda <- control[["lambda"]]
  forget <- control[["forget"]]
  frame  <- control[["frame"]]
  
  nUnits <- length(x)
  s  <- !is.na(x)
  nc <- cumsum(s)  # number SAMPLED to-date
  Nc <- cumsum(rep(1,nUnits)) # number UNITS to-date
  Nu <- Nc - nc # number UNSAMPLED to-date
  
  dat_process <- matrix( 0, nrow=nUnits, ncol=9 )
  
  # cols 1,2,3,4  : Time, x, n, N
  # cols 5,6      : Posterior parameters
  # cols 7,8,9    : Posterior estimates on p
  col_names <- c("tidx", "x", "n", "N", 
                 "a", "b", 
                 "pest", "plcl", "pucl" )
  colnames(dat_process) <- col_names
  dat_process[,1] <- tidx
  dat_process[,2] <- x
  dat_process[,3] <- nc
  dat_process[,4] <- Nc
  
  ii_start <- 1
  if( !is.null(saber_fit) ){
    dat_process[ 1:nrow(saber_fit), ] <- as.matrix(saber_fit)
    ii_start <- nrow(saber_fit)+1
  }
  
  # Inference about the PROCESS
  for( ii in ii_start:nUnits ){
    # Compute the posterior
    posterior <- saber_posterior(
      x = x[1:ii], 
      n = rep(1,ii),
      tidx = tidx[1:ii], 
      prior=prior, lambda=lambda, forget=forget 
    )
    
    # Estimate for SAMPLED results
    pmed <- qbeta(     0.5, posterior[1], posterior[2] )
    plcl <- qbeta(   alpha, posterior[1], posterior[2] )
    pucl <- qbeta( 1-alpha, posterior[1], posterior[2] )
    
    dat_process[ii,5:9] <- c(
      posterior[1], posterior[2],
      pmed, plcl, pucl
    )
  }
  
  ## Create output object
  out <- as.data.frame(dat_process)
  
  out <- new_saber_fit(
    x = out,
    subclass = "saber_cps",
    control = control,
    alpha = alpha,
    call = match.call()
  )
  
  return(out)
  
}




#' @title Probability for sampling next unit in CPS
#' 
#' @description 
#' Computes a probability for sampling subsequent units in a continuous 
#' production sampling framework.
#' 
#' 
#' @param Xsim_mc (matrix) The MCMC for the current time step.
#' @param tidx (numeric) Vector of time points associated with each unit.
#' @param preq (numeric) The requirement on the defect rate
#' @param forecast_type (string) Type of forecast, see details.
#' @param control object of class [saber_control()]
#' 
#' 
#' @details
#' Argument `forecast_type` can take values:
#' 
#' - "optimistic" : assumes the next unit is not a defect
#' - "pessimistic": assumes the next unit is a defect
#' - "estimated": uses a weighted average of the optimistic and pessimistic probabilities
#' 
#' 
#' @return 
#' A numeric scalar giving the probability of sampling the next unit.
#' 
#' @seealso [cps_bound()]
#' 
#' @examples
#' 
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
#' cps_prob_sample(
#'   saber_fit = fit,
#'   pest = fit[["pest"]][nrow(fit)],
#'   preq = 0.10,
#'   tidx_next = max(tidx) + 1,
#'   forecast_type = "estimated"
#' )
#' 
#' 
#' @export
#' 
cps_prob_sample <- function(Xsim_mc, tidx, preq, forecast_type, control ){
  
  indx <- length(tidx)
  prior <- control[["prior"]]
  lambda <- control[["lambda"]]
  forget <- control[["forget"]]
  
  # Current estimate of defect rate for each MC sequence
  post_pest <- t(apply( Xsim_mc[1:(indx-1), ], MARGIN=2, FUN=function(XX, tidx_ii){
    out1 <- saber_posterior( XX, rep(1,indx-1), tidx_ii, 
                             prior=prior,  lambda=lambda, forget=forget)
    out1
  },
  tidx_ii = tidx[ 1:(indx-1) ]
  ))
  
  iseq <- seq_len( indx-1 )
  
  # Probability to tail requirement if next x=0
  Xsim_mc[ indx, ] <- 0
  post0 <- t(apply( Xsim_mc[1:indx, ], MARGIN=2, FUN=function(XX, tidx_ii){
    out1 <- saber_posterior( XX, rep(1,indx), tidx_ii, 
                             prior=prior,  lambda=lambda, forget=forget)
    out1
  },
  tidx_ii = tidx
  ))
  
  # Probability to tail requirement if next x=1
  Xsim_mc[ indx, ] <- 1
  post1 <- t(apply( Xsim_mc[1:indx, ], MARGIN=2, FUN=function(XX, tidx_ii){
    out1 <- saber_posterior( XX, rep(1,indx), tidx_ii, 
                             prior=prior,  lambda=lambda, forget=forget)
    out1
  },
  tidx_ii = tidx
  ))
  
  pests <- apply( post_pest, MARGIN=1, FUN=function(ab){ qbeta(0.5,ab[1],ab[2]) } )
  prob_sample_x0 <- apply( 
    post0, MARGIN=1, 
    FUN=function(ab){ min(1, 1 - pbeta(preq, ab[1], ab[2])) } 
  )
  prob_sample_x1 <- apply( 
    post1, MARGIN=1, 
    FUN=function(ab){ min(1, 1 - pbeta(preq, ab[1], ab[2])) }
  )
  
  prob_sample_wt <- median( prob_sample_x0*(1-pests) + prob_sample_x1*(pests) )
  prob_sample <- switch(
    forecast_type,
    optimistic  = prob_sample_x0,
    pessimistic = prob_sample_x1,
    estimated   = prob_sample_wt
  )
  
  return(prob_sample)
  
}



#' @title Simulate a case for continuous production sampling
#' 
#' @description 
#' Simulates a single run for continuous production sampling and applies 
#' SABER-CPS. Used for assessing performance and parameterizing the forgetting 
#' function.
#' 
#' @param y (numeric) Vector of the true responses.
#' @param n (numeric) Vector denoting whether a unit was samples (n=1) or not (n=0), or not yet produced (n=NA).
#' @param p (numeric) Scalar or vector of the true defect rate.
#' @param tidx (numeric) Scalar or vector of the true defect rate.
#' @param preq (numeric) Requirement on defect rate.
#' @param alpha (numeric) Tail area for credibility level (`alpha = 0.05` means 95% credibility).
#' @param control object of class `saber_control`.
#' @param cps_control object of class `cps_control`.
#' @param seed (numeric) Seed for the random number generator.
#' 
#' 
#' @details
#' Either `y` or `p` must be provided. The total number of units is determined by the length of `n`.
#' If provided, `y` is used instead of `p`. When providing `p`, it may be vectorized (e.g., to 
#' simulate drift or step-changes).
#' 
#' Argument `control` should be the result of [saber_control()], setting the 
#' parameterization for the SABER model. It must contain the following elements:
#' 
#' - `prior`: A vector of length 2 denoting the starting prior.
#' - `lambda`: The parameter for the forgetting function.
#' - `forget`: Character string or function. 
#' - `frame`: The frame of reference; "total", "unsampled", or "accepted".
#' 
#' Argument `cps_control` should be the result of [cps_control()], setting the 
#' parameterization for the simulation. It must contain the following elements:
#' 
#' - `nSim`  : Number of Monte Carlo simulations for defect-count estimation.
#' - `forecast_type`  : Forecasting mode.
#' 
#' 
#' @return 
#' A list with the following elements:
#' 
#' - `saber_fit`: An object of class `saber_cps` containing the final fitted
#'   SABER-CPS trajectory.
#' - `y`: The true response for each unit (may have been supplied, may have been generated).
#' - `x`: The observed response for each unit.
#' - `Xsim_mc`: A matrix of Monte Carlo simulated defective-count contributions
#'   for the CPS sequence.
#' - `prob_sampled`: A numeric vector giving the sampling probability assigned
#'   to each unit.
#' - `control`: The `saber_control` object used for the SABER model.
#' - `cps_control`: The `cps_control` object used for the CPS simulation.
#' 
#'  
#' @examples
#' 
#' \dontrun{
#' control <- saber_control(
#'   prior = c(1, 1),
#'   lambda = 0.95,
#'   forget = "exponential",
#'   frame = "total"
#' )
#' 
#' cps_ctrl <- cps_control(
#'   nSim = 100,
#'   forecast_type = "estimated"
#' )
#' 
#' n_vector   <- rep( c(1,NA), c(20,80) )
#' unit_times <- cumsum(rexp(length(n_vector)))
#' 
#' sim <- cps_sim_case(
#'   p = 0.05,
#'   n = n_vector,
#'   tidx = unit_times,
#'   preq = 0.10,
#'   alpha = 0.05,
#'   control = control,
#'   cps_control = cps_ctrl,
#'   seed = 123
#' )
#'
#' names(sim)
#' class(sim$saber_fit)
#' }
#' 
#' @export
#' 
cps_sim_case <- function( y, n, p, tidx, preq, alpha, control, cps_control, seed=NULL, verbose=FALSE ){
  
  # y = true data
  # n = vector of observed (n=1) / not (n=0). Should be NA if unit not yet "seen"
  
  control     <- rsaber:::validate_saber_control(control)
  cps_control <- rsaber:::validate_cps_control(cps_control)
  
  frame <- control[["frame"]]
  prior <- control[["prior"]]
  lambda <- control[["lambda"]]
  forget <- control[["forget"]]
  
  nSim          <- cps_control[["nSim"]]
  forecast_type <- cps_control[["forecast_type"]]
  
  if( !is.null(seed) ){ set.seed(seed) }
  if( !missing(y) ){
    N     <- length(n)
    y_pop <- y
  } else if( missing(y) & !missing(p) ){
    N     <- length(n)
    y_pop <- rbinom( N, 1, p ) # True fail/pass
  } else{
    stop( "Must provide either 'y' or both 'p' ")
  }
  
  n0 <- sum( !is.na(n) )
  s_pop <- rep( c(1,0), c(n0,N-n0) )
  x_pop <- rep( NA, N )
  
  indx_init <- which( !is.na(n[1:n0]) )
  x_pop[ indx_init ] <- y_pop[ indx_init ] # observed fail/pass
  
  # initial fit and monte carlo
  saber_out <- cps_bound(
    x       = x_pop[ 1:n0 ],
    tidx    = tidx[ 1:n0 ],
    alpha   = alpha,
    control = control
  )
  
  # Initialize the posterior matrix.
  posterior <- saber_posterior( 
    x_pop[1:n0], rep(1,n0), tidx[1:n0], 
    prior=prior, lambda=lambda, forget=forget
  )
  
  # Initialize the Defect simulation
  Xsim_mc <- matrix(NA, nrow = N, ncol = nSim )
  Xsim_mc[ seq_len(n0), ] <- sim_ndefect( saber_out, nSim=nSim )
  
  ## results objects
  prob_sampled <- c( rep(1,n0), rep(NA,N-n0) )
  
  for( indx in (n0+1):N ){
    # Update the indexes
    if( verbose & (indx %% 10 == 0) ){ cat( "n =", indx, "\n") }
    
    
    prob_sample <- cps_prob_sample(
      Xsim_mc, tidx=tidx[1:indx], preq=preq,
      forecast_type=forecast_type, control=control
    )
    prob_sampled[indx] <- prob_sample
    
    if( rbinom(1,1,prob_sample)==1 ){
      x_pop[ indx ]  <- y_pop[ indx ]
      Xsim_mc[indx,] <- y_pop[ indx ]
    } else{
      Xsim_mc[indx,] <- rsaber::rbetabinomial( nSim, Nn=1, 
                                               posterior[1], 
                                               posterior[2] )
    }
    
    # Update posterior
    posterior <-  saber_posterior(
      x = x_pop[1:indx],
      n = ifelse( is.na(x_pop[1:indx]), NA, 1),
      tidx = tidx[ 1:indx ],
      prior=control[["prior"]],
      lambda=control[["lambda"]],
      forget=control[["forget"]]
    )
    
  } # end-while
  
  saber_out <- cps_bound(
    x       = x_pop,
    tidx    = tidx,
    alpha   = alpha,
    control = control
  )
  
  
  return(
    list(
      saber_fit = saber_out,
      y = y_pop,
      x = x_pop,
      Xsim_mc  = Xsim_mc ,
      prob_sampled = prob_sampled,
      control = control,
      cps_control = cps_control
    )
  )
  
}








