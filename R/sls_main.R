





#' @title SABER for Sequential Lot Sampling
#' 
#' @description Perform estimation and inference using the SABER method. Applies 
#' to a sequence of lots.
#' 
#' @param x (numeric vector) Sequence containing the number of defects.
#' @param n (numeric vector) Sequence containing the number of units samples.
#' @param N (numeric vector) Sequence containing the lot size.
#' @param tidx (numeric vector) Sequence of the time the lot was produced/tested.
#' @param alpha (numeric) Setting the credibility level.
#' @param control An object of class `saber_control`, typically created by
#'   [saber_control()].
#'  
#' @return 
#' An object of class `saber_sls`, inheriting from `saber_fit` and
#' `data.frame`. The returned data frame contains one row per lot and columns:
#' `tidx`, `x`, `n`, `N`, `a`, `b`, `pest`, `plcl`, `pucl`, `Xest`, `Xlcl`,
#' and `Xucl`.
#' 
#' @details
#' The SABER posterior is updated sequentially over lots using the forgetting
#' settings in `control`. The fitted posterior parameters are stored in columns
#' `a` and `b`.
#'
#' The frame of reference is set by `control[["frame"]]`. Available frames are:
#'
#' - `"total"`: considers the full lot, including sampled and unsampled units.
#' - `"unsampled"`: considers only the unsampled portion of the lot.
#' - `"accepted"`: considers tested-good units plus unsampled units, excluding
#'   observed defective units.
#'
#' For details on SABER controls and forgetting functions, see
#' [saber_control()] and [saber_weights()].
#'
#' @seealso [saber_control()], [saber_weights()], [saber_posterior()],
#'   [sim_ndefect()] 
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
#' fit <- sls_bound(
#'   x = c(0, 1, 0, 2),
#'   n = c(10, 10, 10, 10),
#'   N = c(100, 100, 100, 100),
#'   tidx = seq_len(4),
#'   alpha = 0.05,
#'   control = control
#' )
#'
#' fit
#' 
#' 
#' @export
#' 
sls_bound <- function( x, n, N, tidx, alpha, control ){
  
  control <- validate_saber_control(control)
  
  prior  <- control[["prior"]]
  lambda <- control[["lambda"]]
  forget <- control[["forget"]]
  frame  <- control[["frame"]]
  
  nLots <- length(x)
  dat_process <- matrix( 0, nrow=nLots, ncol=12 )
  
  # cols 1,2,3,4  : Time, x, n, N
  # cols 5,6      : Posterior parameters
  # cols 7,8,9    : Posterior estimates on p
  # cols 10,11,12 : Posterior estimates on number of defects (U or K)
  col_names <- c("tidx", "x", "n", "N", 
                 "a", "b", 
                 "pest", "plcl", "pucl", 
                 "Xest", "Xlcl", "Xucl")
  dat_process[,1] <- tidx
  dat_process[,2] <- x
  dat_process[,3] <- n
  dat_process[,4] <- N
  
  # Inference about the PROCESS, not about the LOT
  for( ii in 1:nLots ){
    
    # Get some current values for simplicity
    nui <- N[ii] - n[ii]
    
    # Compute the posterior
    posterior <- saber_posterior(
      x = x[ seq_len(ii) ], 
      n = n[ seq_len(ii) ], 
      tidx = tidx[ seq_len(ii) ], 
      prior=prior, lambda=lambda, forget=forget 
    )
    
    # Compute the posterior
    if( is.infinite(dat_process[ii,4]) ){
      # Infinite population case
      pmed <- qbeta(     0.5, posterior[1], posterior[2] )
      plcl <- qbeta(   alpha, posterior[1], posterior[2] )
      pucl <- qbeta( 1-alpha, posterior[1], posterior[2] )
      xmed <- NA
      xlcl <- NA
      xucl <- NA
    } else{
      # Finite population case
      
      # Inference on unsampled portion, NUMBER of defective
      umed <- rsaber::qbetabinomial(     0.5, nui, posterior[1], posterior[2] )
      ulcl <- rsaber::qbetabinomial(   alpha, nui, posterior[1], posterior[2] )
      uucl <- rsaber::qbetabinomial( 1-alpha, nui, posterior[1], posterior[2] )
      
      if( frame=="total"){
        # Consider everything, the full population
        xmed <- (x[ii] + umed) 
        xlcl <- (x[ii] + ulcl) 
        xucl <- (x[ii] + uucl) 
        
        pmed <- (x[ii] + umed) / N[ii]
        plcl <- (x[ii] + ulcl) / N[ii] 
        pucl <- (x[ii] + uucl) / N[ii]
      }
      if( frame=="accepted"){
        # Consider only accepted units: Tested & good + untested 
        # Context: Non-destructive tests, non-repairable units
        xmed <- umed
        xlcl <- ulcl
        xucl <- uucl
        
        pmed <- umed / (N[ii]-x[ii])
        plcl <- ulcl / (N[ii]-x[ii])
        pucl <- uucl / (N[ii]-x[ii])
      }
      if( frame=="unsampled"){
        # Consider only untested units
        # Context: Destructive tests.
        # Use with cauton, can give highly unintuitive results
        xmed <- umed
        xlcl <- ulcl
        xucl <- uucl
        
        pmed <- umed / nui
        plcl <- ulcl / nui
        pucl <- uucl / nui
      }
    }
    
    # Collect the results
    dat_process[ii,5:12] <- c(
      posterior[1], posterior[2],
      pmed, plcl, pucl,
      xmed, xlcl, xucl
    )
    
  }
  
  # create the return object
  out <- as.data.frame(dat_process)
  colnames(out) <- col_names
  out <- new_saber_fit(
    x = out,
    subclass = "saber_sls",
    control = control,
    alpha = alpha,
    call = match.call()
  )
  return(out)
}

