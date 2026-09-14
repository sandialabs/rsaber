




#' @title Compute SABER forgetting weights
#' 
#' @description Evaluates the phi function in SABER for weighting data.
#' 
#' 
#' @param n Numeric vector giving the sample size or binomial denominator
#'   associated with each observation.
#' @param tidx Numeric vector giving the time/order index for each observation.
#'   Used to determine recency in the forgetting function.
#' @param lambda Numeric scalar or vector giving the parameter(s) for the forgetting function.
#' @param forget Character string or function specifying the forgetting method.
#'   Character options are `"exponential"`, `"fixed_window"`,
#'   `"total_samples"`, and `"none"`. 
#'  
#' @return Returns a numeric vector
#' 
#' @note
#' The meaning of `n` and `lambda` changes depending on the forgetting function 
#' and type of SABER model. For SABER-SLS (when running [sls_bound()]), `n` 
#' denotes the sample size from each lot. For SABER-CPS (when running 
#' [cps_bound()]), `n` is a vector with all 1's.
#' 
#' @details
#' Custom functions for `forget` should accept arguments `n`, `tidx`, and 
#' `lambda`, and return a numeric weight vector with the same length as `n`.
#' 
#' 
#' @examples
#' n <- rep(10, 8)
#' tidx <- seq_along(n)
#'
#' # Exponential forgetting: older observations receive smaller weights.
#' saber_weights(
#'   n = n,
#'   tidx = tidx,
#'   lambda = 0.90,
#'   forget = "exponential"
#' )
#'
#' # Fixed-window forgetting: only observations within the window receive weight.
#' saber_weights(
#'   n = n,
#'   tidx = tidx,
#'   lambda = 3,
#'   forget = "fixed_window"
#' )
#'
#' # No forgetting: all observations receive weight one.
#' # Can also be accomplished withforget= "exponential" and lambda=1
#' saber_weights(
#'   n = n,
#'   tidx = tidx,
#'   lambda = 1,
#'   forget = "none"
#' ) 
#' 
#' 
#' @export
#' 
saber_weights <- function( n, tidx, lambda, forget="exponential" ){
  
  tdiff <- abs(tidx - max(tidx))
  
  if( is.function(forget) ){
    iLambda <- forget( n, tidx, lambda )
    # Some checks on the result of a custom forgetting function
    if( !is.numeric(iLambda) ){
      stop("Output of custom forgetting functon must be numeric.")
    }
    if( length(iLambda) != length(n) ){
      stop("Output of custom forgetting functon must have same length as n.")
    }
    if( max(iLambda) > 1 ){ stop("Output of custom forgetting function should not exceed 1.") }
    if( min(iLambda) < 0 ){ stop("Output of custom forgetting function should not fall below 0.") }
    if( iLambda[length(n)] < 1 ){ 
      ## See about only issuing this warning once per R session.
      ## Way to do this with rlang package ... can we do it manually instead?
      warning("Output of forgetting function does not weight most recent lot at phi=1. This is suspect.") 
    }
  } else if( is.character(forget) ){
    # Gets just the weights, this facilitates plotting the weights
    forget <- tolower(forget)
    
    if( forget =="exponential" ){
      # Exponential decay of weight ... may be "too fast"
      iLambda <- lambda^tdiff
      #iLambda <- iLambda / iLambda[length(iLambda)]
    } else if( forget == "fixed_window" ){
      # Fixed number of batches
      iLambda <- ifelse( tdiff <= lambda, 1, 0 )
    } else if( forget == "total_samples" ){
      # Maintain approximate total sample size
      rev_sum_sum <- rev(cumsum(rev(n)))
      idx <- max( c(1, which.max( !(rev_sum_sum >= lambda) )-1 ) )
      iLambda <- rep( 0, length(tidx) )
      iLambda[(idx):length(tidx)] <- 1
    } else if( forget == "none" ){
      iLambda <- rep( 1, length(tidx) )
    }
  } else{
    stop_message <- "Argument 'forget' must be either function or character 
    matching one of the included forgetting methods."
    stop(stop_message)
  }
  return( iLambda )
}



#' @title Posterior parameters for SABER
#' 
#' @description Obtain the posterior parameters for SABER
#' 
#' @param x Numeric vector giving the observed number of defectives for each
#'   observation.
#' @param n Numeric vector giving the sample size or binomial denominator
#'   associated with each observation.
#' @param tidx Numeric vector giving the time/order index for each observation.
#'   Used to determine recency in the forgetting function.
#' @param prior Numeric vector of length 2 giving the beta prior parameters.
#' @param lambda Numeric scalar or vector giving the parameter for the forgetting function.
#' @param forget String; the type of forgetting function.
#' 
#' @return A numeric vector of length 2 containing the posterior 
#' beta parameters `c(a, b)`.
#'  
#' @seealso [saber_weights()]
#' 
#' 
#' @examples
#' x <- c(0, 1, 0, 2)
#' n <- c(10, 10, 10, 10)
#' tidx <- seq_along(x)
#'
#' saber_posterior(
#'   x = x,
#'   n = n,
#'   tidx = tidx,
#'   prior = c(1, 1),
#'   lambda = 0.95,
#'   forget = "exponential"
#' )
#'
#' saber_posterior(
#'   x = x,
#'   n = n,
#'   tidx = tidx,
#'   prior = c(1, 1),
#'   lambda = 2,
#'   forget = "fixed_window"
#' ) 
#' 
#' @export
#' 
saber_posterior <- function( x, n, tidx, prior=c(1,1), lambda, forget="exponential" ){
  
  ## Argument validation
  arg_lens <- c(length(x), length(n), length(tidx))
  if( diff(range(arg_lens)) > 0 ){
    stop("Vector of x, n, and tidx have different lengths")
  }
  
  iLambda <- saber_weights( n, tidx, lambda, forget )
  aa <- prior[1] + sum( iLambda *  x   , na.rm=TRUE )
  bb <- prior[2] + sum( iLambda * (n-x), na.rm=TRUE )
  
  return( c("a"=aa, "b"=bb) )
}





# saber_control , possibly split into separate file eventually ------

#' Construct SABER model control settings
#'
#' @param prior Numeric vector of length 2 giving beta prior parameters.
#' @param lambda Forgetting parameter.
#' @param forget Forgetting method or custom forgetting function.
#' @param frame Frame of reference for defect estimation.
#' 
#' 
#' @examples
#' # Exponential forgetting with the total-population frame.
#' control <- saber_control(
#'   prior = c(1, 1),
#'   lambda = 0.95,
#'   forget = "exponential",
#'   frame = "total"
#' )
#'
#' control
#'
#' # No forgetting.
#' saber_control(
#'   prior = c(1, 1),
#'   forget = "none",
#'   frame = "total"
#' ) 
#' 
#' @export
#' 
saber_control <- function(prior = c(1, 1),
                          lambda = NULL,
                          forget = c("exponential", "fixed_window", "total_samples", "none"),
                          frame = c("total", "unsampled", "accepted")) {
  
  frame <- match_saber_frame(frame)
  
  if (is.character(forget)) {
    forget <- match.arg(forget)
  }
  
  if (!is.numeric(prior) || length(prior) != 2 || any(!is.finite(prior)) || any(prior <= 0)) {
    stop("Argument 'prior' must be a positive numeric vector of length 2.", call. = FALSE)
  }
  
  if (!is.character(forget) && !is.function(forget)) {
    stop("Argument 'forget' must be a character string or function.", call. = FALSE)
  }
  
  if (is.null(lambda)) {
    if (identical(forget, "none")) {
      lambda <- 1
    } else {
      stop("Argument 'lambda' must be supplied unless forget = 'none'.", call. = FALSE)
    }
  }
  
  if (!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda)) {
    stop("Argument 'lambda' must be a finite numeric scalar.", call. = FALSE)
  }
  
  structure(
    list(
      prior = prior,
      lambda = lambda,
      forget = forget,
      frame = frame
    ),
    class = "saber_control"
  )
}





#' Validate SABER model control settings
#'
#' @keywords internal
validate_saber_control <- function(control) {
  
  if (is.null(control)) {
    stop("Argument 'control' must be supplied.", call. = FALSE)
  }
  
  if (!is.list(control)) {
    stop("Argument 'control' must be a list, preferably created by saber_control().", call. = FALSE)
  }
  
  required_names <- c("prior", "lambda", "forget", "frame")
  missing_names <- setdiff(required_names, names(control))
  
  if (length(missing_names) > 0) {
    stop(
      "Argument 'control' is missing required elements: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  prior <- control[["prior"]]
  lambda <- control[["lambda"]]
  forget <- control[["forget"]]
  frame <- control[["frame"]]
  
  if (!is.numeric(prior) || length(prior) != 2 || any(!is.finite(prior)) || any(prior <= 0)) {
    stop("control[['prior']] must be a positive numeric vector of length 2.", call. = FALSE)
  }
  
  if (!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda)) {
    stop("control[['lambda']] must be a finite numeric scalar.", call. = FALSE)
  }
  
  if (!is.character(forget) && !is.function(forget)) {
    stop("control[['forget']] must be a character string or function.", call. = FALSE)
  }
  
  if (is.character(forget)) {
    control[["forget"]] <- match.arg(
      forget,
      choices = c("exponential", "fixed_window", "total_samples", "none")
    )
  }
  
  control[["frame"]] <- match_saber_frame(frame)
  
  if (!inherits(control, "saber_control")) {
    class(control) <- c("saber_control", class(control))
  }
  
  control
}



















