


#' @title Bayesian point and interval estimation for proportions
#' 
#' @description 
#' Uses Bayesian approach to estimate the proportion for binary event data. 
#' For the finite-population case, a superpopulation approach is used, leading to a 
#' Beta-Binomial posterior distribution. The prior distribution defaults to a uniform 
#' (non-informative) distribution.
#' 
#' 
#' @param x Numeric (integer); the number of events observed
#' @param n Numeric (integer); the number of samples taken
#' @param N Numeric (integer); Size of the population (can be `Inf` for infinite-population case)
#' @param alpha numeric, complement of the confidence level
#' @param prior numeric vector of the Beta prior, defautls to `c(1,1)`.
#' @param mesh character string identifying how to mesh the parameters, see `parameter_mesh`.
#' @param frame character string identifying the frame of reference for finite populations, 
#' currently accepted values are "total", "accepted", or "unsampled."
#' 
#' @details Computes a "margin of error" as the larger of the differences 
#'          $\eqn{X_{1-\alpha}-X_{0.5}}$ and $\eqn{X_{0.5}-X_{\alpha}}$
#' 
#' For finite populations, the posterior distribution is Beta-Binomial with 
#' population size N-n (that is, there are N-n untested units in the population). 
#' There are several frames of reference which may be of interest here.
#' 
#' - `frame="total"` corresponds to the entire population, including the 
#' units that were tested, and those which are defective.
#' - `frame="accepted"` corresponds to only accepted units. This will *include* tested units 
#' but will *exclude* tested units that are defective.
#' - `frame="unsampled"` corresponds to only untested units. This will exclude 
#' all tested units, regardless of whether they are defective or not.
#' 
#' The "accepted" case relates to non-destructive (and non-degrading) testing, 
#' while "unsampled" relates to destructive testing.
#' 
#' 
#' @return 
#' Returns a list, or a tibble with columns for point and interval estimates. 
#' For the finite-population case, this includes estimates of both the proportion 
#' of events as well as the number of events in the population (*including* the observed events).
#' 
#' @note 
#' The Hypergeometric interval on p for finite population uses Bayesian estimation.
#' The default is a Beta(1,1) prior on p. An example the BetaBinomial is noninformative: 
#' hist(  rsaber::qbetabinomial( runif(10000,0,1) , Nn=100, a=1, b=1)  )
#' 
#' 
#' @examples
#' prop_interval_bayes(
#'   x = 0,
#'   n = c(10,20,30,40),
#'   N = 100,
#'   alpha = c(0.1, 0.05),
#'   prior = c(1,1)
#' )
#' 
#' 
#' 
#' @export
#' 
prop_interval_bayes <- function( x, n, N, alpha, prior=c(1,1), mesh="crossing", frame="total" ){
  
  # Mesh the input parameters
  frame = match_saber_frame(frame)
  params <- list( x, n, N, alpha )
  names(params) <- c("x", "n", "N", "alpha")
  dat00 <- parameter_mesh( mesh, params ) |> 
    subset( x <= n ) |> 
    subset( n <= N )
  if( nrow(dat00)==0 ){ 
    stop("Removing cases where x>n and n>N left zero cases remaining. Check inputs.")
  }
  
  xx <- dat00[["x"]]
  nn <- dat00[["n"]]
  NN <- dat00[["N"]]
  
  aa <- xx + prior[1]
  bb <- nn - xx + prior[2]
  Nn <- NN - nn
  cc <- dat00[["alpha"]]
  
  ## Infinite-population estimation
  ## xmed, xlcl, xucl do not apply, but are needed to fill the dataframe correctly
  pmed <- qbeta( 0.5    , aa, bb ) 
  plcl <- qbeta( cc  , aa, bb ) 
  pucl <- qbeta( 1-cc, aa, bb ) 
  xmed <- rep(NA, nrow(dat00) )
  xlcl <- rep(NA, nrow(dat00) )
  xucl <- rep(NA, nrow(dat00) )
  
  ## Finite-population estimation
  ## Needs to be in and if() since it's vectorized
  ## Otherwise rmutil::qbetabinom fails for 
  finite_N <- is.finite( dat00[["N"]] )
  if( any(finite_N) ){
    
    xxf <- xx[ finite_N ]
    NNf <- NN[ finite_N ]
    Nnf <- Nn[ finite_N ]
    aaf <- aa[ finite_N ] 
    bbf <- bb[ finite_N ] 
    ccf <- cc[ finite_N ]
    
    umed <- rsaber::qbetabinomial( 0.5 , Nnf, aaf, bbf )
    ulcl <- rsaber::qbetabinomial( ccf  , Nnf, aaf, bbf )
    uucl <- rsaber::qbetabinomial( 1-ccf, Nnf, aaf, bbf )
    
    if( frame=="total" ){
      xmed[ finite_N ] <- umed + xxf
      xlcl[ finite_N ] <- ulcl + xxf
      xucl[ finite_N ] <- uucl + xxf

      pmed[ finite_N ] <- (umed + xxf) / NNf
      plcl[ finite_N ] <- (ulcl + xxf) / NNf
      pucl[ finite_N ] <- (uucl + xxf) / NNf
    }
    if( frame=="accepted" ){
      xmed[ finite_N ] <- umed
      xlcl[ finite_N ] <- ulcl
      xucl[ finite_N ] <- uucl

      pmed[ finite_N ] <- (umed) / (NNf-xxf)
      plcl[ finite_N ] <- (ulcl) / (NNf-xxf)
      pucl[ finite_N ] <- (uucl) / (NNf-xxf)
    }
    if( frame=="unsampled" ){
      xmed[ finite_N ] <- umed
      xlcl[ finite_N ] <- ulcl
      xucl[ finite_N ] <- uucl

      pmed[ finite_N ] <- (umed) / (Nnf)
      plcl[ finite_N ] <- (ulcl) / (Nnf)
      pucl[ finite_N ] <- (uucl) / (Nnf)
    }
    
  }
  
  dat00 <- dat00 |>
    transform(
      pest = pmed,
      plcl = plcl,
      pucl = pucl,
      xest = xmed,
      xlcl = xlcl,
      xucl = xucl
    )
  
  # Return the results
  return( dat00 )
  
}




#' @title Bayesian method to determine sample size for reliability requirement
#' 
#' @description 
#' Bayesian approach of obtaining minimum sample size needed to assert
#' a specific reliability requirement (defined in terms of the defect rate).
#' For the finite-population case, a superpopulation approach is used, leading to a 
#' Beta-Binomial posterior distribution. The prior distribution defaults to a uniform 
#' (non-informative) distribution.
#' 
#' 
#' @param x Numeric (integer); the number of events observed
#' @param N Numeric (integer); Size of the population (can be `Inf` for infinite-population case)
#' @param preq Numeric (decimal); the requirement on the defect rate
#' @param alpha Numeric (decimal); complement of the confidence level
#' @param prior Numeric (float); value of \eqn{\alpha} and \eqn{\beta} in the prior distribution for \eqn{\pi}. Defautls to `c(1,1)`.
#' @param mesh character string identifying how to mesh the parameters, see `parameter_mesh`.
#' @param frame Character; the frame of reference for finite populations, currently accepted values are "total", "accepted", or "unsampled."
#' @param maxiter Numeric (integer); maximum number of iterations to permit. See notes.
#' @param verbose Logical; whether or not to print non-error messages.
#' 
#' @return 
#' Returns a tibble with columns for N, preq, alpha, and smallest sample size.
#' 
#' @note
#' The function will estimate an upper bound based on a zero-failure scheme. If 
#' x>0 defects are allowed, this may not be a large enough upper bound. The 
#' function uses a while loop to slowly increase the sample size, the `maxiter` 
#' parameter is there to prevent infinite loops.
#' 
#' @examples
#' prop_n_bayes(
#'   x = 0,
#'   N = 100,
#'   preq = 0.05,
#'   alpha = c(0.1, 0.05),
#'   prior = c(1,1)
#' )
#' 
#' @export
#' 
prop_n_bayes <- function( x, N, preq, alpha, prior=c(1,1), mesh="crossing", frame="total", maxiter=100, verbose=TRUE ){
  
  frame = match_saber_frame(frame)
  # Mesh the input parameters
  params <- list( x, N, preq, alpha )
  names(params) <- c("x", "N", "preq", "alpha")
  dat00 <- parameter_mesh( mesh, params ) |>
    transform(
      n = NA
    )
  
  for( ii in 1:nrow(dat00) ){
    xi <- dat00[["x"]][ii]
    Ni <- dat00[["N"]][ii]
    pi <- dat00[["preq"]][ii]
    ai <- dat00[["alpha"]][ii]
    
    ## Accommodate the infinite-population case
    if( is.infinite(Ni) ){
      
      flag_while <- TRUE
      iter <- 0
      nstep <- ceiling( log(ai)/log(1-pi) )*2

      while( flag_while ){
        iter <- iter + 1
        nlb <- max( xi  , 1 + nstep*(iter-1))
        nub <- max( xi+1, nstep*iter )
        n <- seq( nlb, nub, 1 )
        pib_out <- prop_interval_bayes( xi, n, Ni, ai, prior=prior, frame=frame )

        if(  any( pib_out[["pucl"]] <= pi ) ){
          flag_while <- FALSE
        }
        if( iter == maxiter ){
          stop("Reached maximum iterations. Increase maxiter or verify inputs are sensible.")
        }
      }

    } else{
      n <- seq( max(c(xi,1)), Ni-1, 1)
    }

    dtmp00 <- prop_interval_bayes( xi, n, Ni, ai, prior=prior, frame=frame )
    dtmp01 <- dtmp00[ dtmp00[["pucl"]] <= pi , ]
    if (length(dtmp01[["n"]]) > 0) {
      dat00[["n"]][ii] <- min( dtmp01[["n"]] )
    }
  }
  
  # Return the results
  if ( (verbose==TRUE) & any(is.na(dat00[["n"]])) ) {
    message("Note: preq cannot be achieved in some cases due to allowable x")
  }
  return( dat00 )
}







#' @title Frequentist point and interval estimation proportions
#' 
#' @description 
#' Implements Frequentist approaches to estimate the proportion for binary 
#' event data. Methods available are Wilson, Agresti-Coull, Clopper-Pearson, Jeffreys, and Wald. 
#' Applies a finite-population correction factor (FPC) when argument `N` is non-infinite.
#' 
#' 
#' @param x Numeric (integer); the number of events observed
#' @param n Numeric (integer); the number of samples taken
#' @param N Numeric (integer); Size of the population (can be `Inf` for infinite-population case)
#' @param alpha Numeric (decimal); complement of the confidence level
#' @param mesh character string identifying how to mesh the parameters, see `parameter_mesh`.
#' @param type Character; the type of bound, currently supported options are 
#' "wilson" (default), "agresti-coull", "clopper-pearson", and "wald"
#' 
#' @details
#' The argument `type` is converted to lowercase, and several abbreviations are permitted 
#' (wils for wilson, agresti, ac, and a-c for agresti-coull). 
#' 
#' The Wald interval tends to demonstrate poor coverage, particularly when p is close 
#' to 0 or 1 and/or n is small. It is included for completeness, but should rarely be 
#' the preferred method. The default method is the Wilson (score) interval.
#' 
#' For the Clopper-Pearson interval is obtained using the Beta representation. This implies 
#' two different "priors". To resolve this for the estimate, we use the result from the 
#' Jeffreys interval.
#' 
#' @return 
#' Returns a list or tibble with columns for point and interval estimates.
#' 
#' @examples
#' prop_interval_freq(
#'   x = 0,
#'   n = c(10,20,30,40),
#'   N = 100,
#'   alpha = c(0.1, 0.05),
#' )
#' 
#' @export
#' 
prop_interval_freq <- function( x, n, N, alpha, mesh="crossing", type="wilson"  ){
  
  accepted_types <- c(
    "wald", "wilson", "agresti-coull", "clopper-pearson"
  )
  type <- tolower( type )
  if( !(type %in% accepted_types ) ){
    warning("Argument `type` not in known list, attempting common matches")
    type <- switch( type , 
                    "wils" = "wilson",
                    "agresti"="agresti-coull",
                    "ac"="agresti-coull",
                    "a-c"="agresti-coull",
                    "cp" = "clopper-pearson",
                    "c-p" = "clopper-pearson"
    )
  }
  
  # Mesh the input parameters
  params <- list( x, n, N, alpha )
  names(params) <- c("x", "n", "N", "alpha")
  dat00 <- parameter_mesh( mesh, params )
  
  xx <- dat00[["x"]]
  nn <- dat00[["n"]]
  NN <- dat00[["N"]]
  alpha <- dat00[["alpha"]]
  kk       <- qnorm( 1 - alpha/2 )
  phat0    <- xx/nn
  pwils    <- (xx + (kk^2)/2) / (nn + kk^2)
  
  fpc <- ifelse( is.infinite(NN),
                 1,
                 sqrt( (NN - nn)/(NN-1) )
                 )
  
  
  if( type=="wilson" ){
    phat <- pwils
    se1_wils <- ( sqrt(nn) / (nn + kk^2) )
    se2_wils <- sqrt( phat0*(1-phat0) + (kk^2 / (4*nn)) )  ## Wilson form of SE
    se_prop <- se1_wils * se2_wils * fpc
    plcl <- phat - kk*se_prop
    pucl <- phat + kk*se_prop
  } else if( type=="agresti-coull" ){
    phat <- pwils
    se_prop <- sqrt( pwils*(1-pwils)/(nn+kk^2) ) * fpc
    plcl <- phat - kk*se_prop
    pucl <- phat + kk*se_prop
  } else if( type == "wald" ){
    phat <- phat0
    se_prop <- sqrt( phat0*(1-phat0)/nn ) * fpc
    plcl <- phat - kk*se_prop
    pucl <- phat + kk*se_prop
  } else if( type %in% c("clopper-pearson") ){
    
    # Beta equivalent for infinite populations
    a1 <- 0 + sqrt(.Machine$double.eps)
    a2 <- 1
    b1 <- 1
    b2 <- 0 + sqrt(.Machine$double.eps)
    
    phat <- qbeta(       0.5, xx+0.5, nn-xx+0.5 )  # Jeffreys-style point estimate
    plcl <- qbeta(   alpha/2, xx+a1, nn-xx+b1 )
    pucl <- qbeta( 1-alpha/2, xx+a2, nn-xx+b2 )
    se_prop <- rep( NA, length(phat) )
    
    finite_N <- is.finite(NN)
    if (any(finite_N)) {
      
      xxf    <- xx[finite_N]
      nnf    <- nn[finite_N]
      NNf    <- NN[finite_N]
      alphaf <- alpha[finite_N]
      
      estf <- lclf <- uclf <- numeric(length(xxf))
      
      for(ii in seq_along(xxf) ){
        xx <- xxf[ii]
        nn <- nnf[ii]
        NN <- NNf[ii]
        
        uu <- 0:(NN - nn)
        KK <- uu + xx
        pp <- KK / NN
        
        prbs_lo <- 1 - phyper(q = xx - 1, m = KK, n = NN - KK, k = nn)
        prbs_hi <- 0 + phyper(q = xx - 0, m = KK, n = NN - KK, k = nn)
        indx_lo <- which( prbs_lo <= alphaf[ii]/2 )
        indx_hi <- which( prbs_hi <= alphaf[ii]/2 )
        
        lclf[ii] <- ifelse( 
          length(indx_lo) == 0, 
          pp[1],
          pp[max(indx_lo)]
        )
        uclf[ii] <- ifelse( 
          length(indx_hi) == 0, 
          pp[length(pp)],
          pp[min(indx_hi)]
        )
        
        indx_m1 <- which( prbs_lo <= 0.5 )
        indx_m2 <- which( prbs_hi <= 0.5 )
        pmed1 <- ifelse( 
          length(indx_m1) == 0, 
          pp[1],
          pp[max(indx_m1)]
        )
        pmed2 <- ifelse( 
          length(indx_m2) == 0, 
          pp[length(pp)],
          pp[min(indx_m2)]
        )
        estf <- 0.5*(pmed1 + pmed2)
        
      }
      
      phat[finite_N] <- estf
      plcl[finite_N] <- lclf
      pucl[finite_N] <- uclf
    }
    
  }
  
  dat00 <- dat00  |>
    transform(
      est = phat,
      se_prop = se_prop,
      plcl = plcl,
      pucl = pucl,
      type = type
    )
  
  return( dat00 )
  
}






#' @title Sample size estimation for reliability
#' 
#' @description 
#' Implementation of several a method of obtaining minimum sample size needed 
#' to assert a specific reliability requirement (defined in terms of the defect 
#' rate). Leverages Clopper-Pearson type upper bound on defect rate.
#' 
#' @rdname prop_freq
#' 
#' @param x Numeric (integer); the number of events observed, defaults to 0
#' @param N Numeric (integer); Size of the population (can be `Inf` for infinite-population case)
#' @param preq Numeric (decimal); the requirement on the defect rate
#' @param alpha Numeric (decimal); complement of the confidence level
#' @param mesh character string identifying how to mesh the parameters, see `parameter_mesh`.
#' @param method Character; Use original Darby ("darby") formulation, 
#' or Lauren Wilson's modification ("wilson").
#' @param maxiter Numeric (integer); maximum number of iterations to permit. See notes.
#' @param preqdigits Numeric (integer); number of digits to round preq, defaults to 8
#' @return 
#' Returns a tibble with columns for N, preq, alpha, and smallest sample size.
#' 
#' @details
#' The method implemented assumes a Clopper-Pearson type (or "exact") one-sided 
#' upper confidence bound on the defect rate, and sets the sample size such that 
#' this upper bound will be less than or equal to the requirement.
#' 
#' @note
#' The function will estimate an upper bound based on a zero-failure scheme. If 
#' x>0 defects are allowed, this may not be a large enough upper bound. The 
#' function uses a while loop to slowly increase the sample size, the `maxiter` 
#' parameter is there to prevent infinite loops. 
#' The `preqdigits` argument was added due to machine precision issues. For 
#' example, a value of 0.05 may produce a different results to 1 - 0.95. This 
#' value defaults to 8, which should be suitable for most cases.
#' 
#' 
#' prop_n_freq(
#'   x = 0,
#'   N = 100,
#'   preq = 0.05,
#'   alpha = c(0.1, 0.05),
#' )
#' 
#' 
#' @export
#' 
prop_n_freq <- function( x=0, N=Inf, preq, alpha, mesh="crossing", method="wilson", maxiter=100, preqdigits=8 ){
  
  # Address small numerical precision issues
  # e.g., setting preq=0.05 vs preq = 1-0.95
  preq <- round(preq, preqdigits)
  
  # Mesh the input parameters
  params <- list( x, N, preq, alpha )
  names(params) <- c("x", "N", "preq", "alpha")
  dat00 <- parameter_mesh( mesh, params ) |> transform(n = NA)
  
  for( ii in 1:nrow(dat00) ){
    
    xi <- dat00[["x"]][ii]
    Ni <- dat00[["N"]][ii]
    pii <- dat00[["preq"]][ii]
    ai <- dat00[["alpha"]][ii]
    
    if( is.infinite(Ni) ){
      ## Infinite-population / Binomial case (Darby/Wilson)
      flag_while <- TRUE
      iter <- 0
      Nstep <- ceiling( log(ai)/log(1-pii) )*2
      while( flag_while ){
        iter <- iter + 1
        Nlb <- 1 + Nstep*(iter-1)
        Nub <- Nstep*iter
        nn <- seq( Nlb, Nub, 1 )
        pp <- pbinom( xi, nn, pii )
        
        if( any(pp <= ai ) ){
          flag_while <- FALSE
        }
        if( iter == maxiter ){
          stop("Reached maximum iterations. Increase maxiter or verify inputs are sensible.")
        }
      }
      
    } else{
      ## Finite-population / Hypergeometric case (Darby/Wilson)
      if( method == "darby" ){
        kval <- round( pii * Ni )   ## Original Darby method
      }
      if( method == "wilson" ){
        kval <- ceiling( pii * Ni )   ## Lauren's implementation
      }
      nn <- seq(1,Ni-1,1)
      pp <- phyper( xi, kval, Ni-kval, nn )
    }
    
    dat01 <- data.frame(
      n     = nn,
      probs = pp
    ) # |> subset( probs <= ai )
    dat01 <- dat01[ dat01[["probs"]] <= ai , ]
    dat00[["n"]][ii] <- min( dat01[["n"]] )
  }
  
  return( dat00 )
}













