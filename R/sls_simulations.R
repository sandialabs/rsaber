




#' @title Simulate lots for SABER SLS
#' 
#' @description 
#' Generates a sequence of lots for use in SABER-SLS simulation studies. For
#' each lot, the total number of defective units is simulated from a binomial
#' distribution with lot size `N` and defect probability `p`.
#' 
#' @param p Numeric scalar or vector of length `nLots` giving the true defect
#'   probability for each lot.
#' @param N Numeric scalar or vector of length `nLots` giving the lot size.
#' @param nLots Positive integer giving the number of lots to generate.
#' @param seed Optional numeric seed for the random number generator.
#' 
#' @return 
#' A data frame with one row per lot and columns:
#'
#' - `Lot`: Lot index.
#' - `N`: Lot size.
#' - `K`: Simulated number of defective units in the lot.
#' 
#' @examples
#' 
#' # Consistent defect rate
#' sls_gen_lots( p=0.05, N=100, nLots = 20, seed=42 )
#' 
#' # Changes in defect rate
#' sls_gen_lots( 
#'   p=rep(c(0.05,0.15,0.05), c(30,30,30) ), 
#'   N=100, 
#'   nLots = 90, 
#'   seed=42 
#' )
#' 
#' # Variable lot size
#' sls_gen_lots( 
#'   p=0.25, 
#'   N=rep( c(100,50), c(10,10) ), 
#'   nLots = 20, 
#'   seed=42 
#' )
#' 
#' @export
#' 
sls_gen_lots <- function( p, N, nLots, seed ){
  if( length(p)==1 ){ p <- rep(p, nLots) }
  if( length(N)==1 ){ N <- rep(N, nLots) }
  if( !missing(seed) ){ set.seed(seed) }
  Kseq <- rbinom( nLots, size=N, prob=p )
  return_dat <- data.frame(
    Lot = seq_len(nLots),
    N = N, 
    K = Kseq
  )
  return( return_dat )
}




#' @title Bayesian estimates for sequential lots
#' 
#' @description 
#' Computes posterior summaries for the defect probability and posterior
#' predictive summaries for the number of defectives in the unsampled portion
#' of the current lot. This is the default estimation function used by
#' [sls_sim_case()], [sls_sim_oc()], and [sls_sim_timeline()]. 
#' 
#' @param kvec Numeric vector giving the sequence of sampled defect counts from
#'   previous and current lots. Included for compatibility with the SLS
#'   estimation-function interface.
#' @param nvec Numeric vector giving the sequence of sample sizes from previous
#'   and current lots.
#' @param Nvec Numeric vector giving the sequence of lot sizes from previous and
#'   current lots.
#' @param posterior Numeric vector of length 2 giving the current posterior beta
#'   parameters.
#' @param alpha Numeric posterior tail probability used for interval summaries.
#' @param frame String What frame of reference is relevant for number of defectives.
#'   
#' @return 
#' Returns a numeric vector of length 6 which contains:
#' 
#' - `pest`: Posterior median of the defect probability.
#' - `plcl`: Posterior lower credible bound of the defect probability.
#' - `pucl`: Posterior upper credible bound of the defect probability.
#' - `Xest`: Posterior predictive median of the number of defectives in the
#'   unsampled portion of the current lot.
#' - `Xlcl`: Posterior predictive lower credible bound for the number of
#'   defectives in the unsampled portion of the current lot.
#' - `Xucl`: Posterior predictive upper credible bound for the number of
#'   defectives in the unsampled portion of the current lot.
#'   
#' @note 
#' Custom estimation functions supplied through [sls_control()] should accept
#' the same arguments and return the same output structure.
#' 
#' @seealso [sls_control()], [sls_sim_case()], [sls_sim_oc()],
#'   [sls_sim_timeline()]
#'   
#' @examples
#' sls_sim_estimate(
#'   kvec = c(0, 1, 0),
#'   nvec = c(10, 10, 10),
#'   Nvec = c(100, 100, 100),
#'   posterior = c(2, 30),
#'   alpha = 0.05,
#'   frame = "total"
#' )
#'  
#' @export
#'    
sls_sim_estimate <- function( kvec, nvec, Nvec, posterior, alpha, frame=c("total", "unsampled", "accepted") ){
  
  frame <- match.arg(frame)
  NN <- Nvec[ length(Nvec) ]
  nn <- nvec[ length(Nvec) ]
  xx <- kvec[ length(Nvec) ]
  
  if( is.infinite(NN) ){
    pest <- qbeta(     0.5, posterior[1], posterior[2] ) 
    plcl <- qbeta(   alpha, posterior[1], posterior[2] ) 
    pucl <- qbeta( 1-alpha, posterior[1], posterior[2] ) 
    Xest <- NA
    Xlcl <- NA
    Xucl <- NA
  } else{
    
    if( nn < NN ){
      umed <- rsaber::qbetabinomial(    0.5 , NN - nn, posterior[1], posterior[2] )
      ulcl <- rsaber::qbetabinomial(   alpha, NN - nn, posterior[1], posterior[2] )
      uucl <- rsaber::qbetabinomial( 1-alpha, NN - nn, posterior[1], posterior[2] )
      if( frame=="total" ){
        Xest <- umed + xx
        Xlcl <- ulcl + xx
        Xucl <- uucl + xx
        
        pest <- (umed + xx) / NN
        plcl <- (ulcl + xx) / NN
        pucl <- (uucl + xx) / NN
      }
      if( frame=="accepted" ){
        Xest <- umed
        Xlcl <- ulcl
        Xucl <- uucl
        
        pest <- (umed) / (NN-xx)
        plcl <- (ulcl) / (NN-xx)
        pucl <- (uucl) / (NN-xx)
      }
      if( frame=="unsampled" ){
        Xest <- umed
        Xlcl <- ulcl
        Xucl <- uucl
        
        pest <- (umed) / (NN - nn)
        plcl <- (ulcl) / (NN - nn)
        pucl <- (uucl) / (NN - nn)
      }
    } else{
      # When NN == nn 
      if( frame=="total" ){
        Xest <- xx
        Xlcl <- xx
        Xucl <- xx
        
        pest <- (xx) / NN
        plcl <- (xx) / NN
        pucl <- (xx) / NN
      }
      if( frame=="accepted" ){
        Xest <- 0
        Xlcl <- 0
        Xucl <- 0
        
        pest <- 0
        plcl <- 0
        pucl <- 0
      }
      if( frame=="unsampled" ){
        Xest <- 0
        Xlcl <- 0
        Xucl <- 0
        
        pest <- 0
        plcl <- 0
        pucl <- 0
      }
    }
    
  }
  
  
  return(
    # c(
    #   pest = qbeta(     0.50, posterior[1], posterior[2]),
    #   plcl = qbeta(    alpha, posterior[1], posterior[2]),
    #   pucl = qbeta(1 - alpha, posterior[1], posterior[2]),
    #   Xest = rsaber::qbetabinomial(     0.50, NN - nn, posterior[1], posterior[2]),
    #   Xlcl = rsaber::qbetabinomial(    alpha, NN - nn, posterior[1], posterior[2]),
    #   Xucl = rsaber::qbetabinomial(1 - alpha, NN - nn, posterior[1], posterior[2])
    # )
    c( pest, plcl, pucl,
       Xest, Xlcl, Xucl )
  )
}





#' @title Simulate one case for evaluating a sequence of lots 
#' 
#' @description 
#' Simulates a sequence of lots and applies a selected sequential lot sampling
#' evaluation rule to each lot. This function is used to assess SLS operating
#' behavior for a single simulated sequence.
#'  
#' @param p Numeric scalar or vector of length `nLots` giving the true defect
#'   probability for each lot.
#' @param N Numeric scalar or vector of length `nLots` giving the lot size.
#' @param preq Numeric requirement on the defect probability.
#' @param alpha Numeric posterior tail probability used for lower and upper
#'   interval summaries.
#' @param nLots Positive integer giving the number of lots to simulate.
#' @param control An object of class `saber_control`, typically created by
#'   [saber_control()].
#' @param sls_control An object of class `sls_control`, typically created by
#'   [sls_control()]. This controls the SLS evaluation method and policy
#'   settings.
#' @param seed Optional numeric seed for the random number generator.
#'  
#' @return 
#' A numeric matrix with one row per lot and columns:
#'
#' - `Lot`: Lot index.
#' - `N`: Lot size.
#' - `K`: True simulated number of defective units in the lot.
#' - `k`: Number of defective units sampled from the lot.
#' - `n`: Number of units sampled from the lot.
#' - `a`: First posterior beta parameter after evaluating the lot.
#' - `b`: Second posterior beta parameter after evaluating the lot.
#' - `dAcc`: Number of defective units accepted from the lot.
#' - `nAcc`: Number of units accepted from the lot.
#' - `lAcc`: Lot acceptance indicator, where 1 indicates accepted and 0
#'   indicates rejected.
#' - `pest`: Posterior median of the defect probability.
#' - `plcl`: Posterior lower credible bound of the defect probability.
#' - `pucl`: Posterior upper credible bound of the defect probability.
#' - `Xest`: Posterior predictive median of the number of defectives in the
#'   unsampled portion of the lot.
#' - `Xlcl`: Posterior predictive lower credible bound for the number of
#'   defectives in the unsampled portion of the lot.
#' - `Xucl`: Posterior predictive upper credible bound for the number of
#'   defectives in the unsampled portion of the lot.
#' 
#' @details
#' The evaluation rule is selected through `sls_control`. Built-in character
#' aliases include:
#'
#' - `"saber"`: adaptive SABER SLS sampling.
#' - `"saber_fixed"`: SABER SLS sampling with a fixed sample size.
#' - `"bayesian_indp"`: independent Bayesian sampling without accumulating
#'   information across lots.
#' - `"skiplot_indp"`: independent Bayesian skip-lot sampling.
#' - `"skiplot_saber"`: SABER-based skip-lot sampling.
#'
#' A custom evaluation function may also be supplied through
#' [sls_control()]. Custom evaluators should use the argument structure
#' documented in [sls_evaluation] and return the same seven-element numeric
#' vector.
#'
#' @seealso [sls_control()], [saber_control()], [sls_evaluation],
#'   [sls_sim_oc()], [sls_sim_timeline()]
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
#' sls_ctrl <- sls_control(
#'   x0 = 0,
#'   nn_min = 5,
#'   eval_function = "saber_fixed"
#' )
#'
#' sls_sim_case(
#'   p = 0.05,
#'   N = 100,
#'   preq = 0.10,
#'   alpha = 0.05,
#'   nLots = 5,
#'   control = control,
#'   sls_control = sls_ctrl,
#'   seed = 123
#' ) 
#' 
#' 
#' @export
#' 
sls_sim_case <- function( p, N, preq, alpha, nLots, 
                          control, sls_control, seed=NULL ) {
  
  control <- validate_saber_control(control)
  sls_control <- validate_sls_control(sls_control)
  
  # Lot, N, K, k, n, a, b, dAcc, nAcc, lAcc, pest, plcl, pucl, Xest, Xlcl, Xucl
  case_results <- sls_gen_lots(p = p, N = N, nLots = nLots, seed=seed ) |> 
    transform(
      k = NA, n = NA, a = NA, b = NA,
      dAcc = NA, nAcc = NA, lAcc = NA,
      pest = NA, plcl = NA, pucl = NA,
      Xest = NA, Xlcl = NA, Xucl = NA
    ) |> as.matrix() # matrix is easier/faster for indexing
  
  eval_function <- sls_control[["eval_function"]]
  estimation_function <- sls_control[["estimation_function"]]
  
  if( is.character(eval_function) ){
    eval_function <- switch( tolower(eval_function),
                             "saber"         = eval_sls_saber_adapt,
                             "saber_fixed"   = eval_sls_saber_fixed,
                             "bayesian_indp" = eval_sls_bayesian_indp,
                             "skiplot_indp"  = eval_sls_skiplot_indp,
                             "skiplot_saber" = eval_sls_skiplot_saber
    )
  }
  
  ab <- control[["prior"]] ## initialize the prior
  
  for (ii in 1:nLots) {
    # Universal, make the lot and scramble it
    Ni <- case_results[ii,2]
    Ki <- case_results[ii,3]
    Lot <- rep(c(1, 0), c(Ki, Ni - Ki))
    Lot <- Lot[sample(1:Ni, Ni, replace = FALSE)]
    
    case_results[ii,4:10] <- eval_function(
      Lot=Lot, preq=preq, alpha=alpha, 
      kvec=case_results[0:(ii-1),4], 
      nvec=case_results[0:(ii-1),5],
      idx=ii, ab=ab, 
      control = control,
      sls_control=sls_control
    )
    
    Nveci <- case_results[ seq_len(ii), 2 ]
    kveci <- case_results[ seq_len(ii), 4 ]
    nveci <- case_results[ seq_len(ii), 5 ]
    ab <- case_results[ii,6:7]
    case_results[ii,11:16] <- estimation_function( 
      kveci, nveci, Nveci, ab, alpha, frame=control[["frame"]]
    )
    
  }
  return( case_results )
}





#' Convert an SLS simulation case to operating-characteristic summaries
#'
#' Internal helper used by [sls_sim_oc()] to summarize the output from
#' [sls_sim_case()] for one simulated sequence of lots.
#'
#' @param case_results Matrix or data frame returned by [sls_sim_case()].
#'
#' @return A named numeric vector containing operating-characteristic summaries.
#'
#' @noRd 
#' 
sls_sim_case_to_oc <- function( case_results ) {
  
  case_results <- case_results |> as.data.frame()
  
  # r1 <- c( sum(  case_results[["nAcc"]] ) / sum( case_results[["N"]]    ) ,
  #          sum(  case_results[["dAcc"]] ) / sum( case_results[["nAcc"]] ),
  #          sum(  case_results[["k"]]),
  #          sum(  case_results[["n"]]),
  #          sum(  case_results[["n"]])/sum(case_results[["N"]]),
  #          sum(  case_results[["nAcc"]]),
  #          mean( case_results[["lAcc"]])
  # )
  r1 <- c(
    prop_accepted           = sum(case_results[["nAcc"]]) / sum(case_results[["N"]]),
    prop_accepted_defective = sum(case_results[["dAcc"]]) / sum(case_results[["nAcc"]]),
    number_defective_seen   = sum(case_results[["k"]]),
    number_sampled          = sum(case_results[["n"]]),
    fraction_sampled        = sum(case_results[["n"]]) / sum(case_results[["N"]]),
    number_accepted         = sum(case_results[["nAcc"]]),
    prop_lots_accepted      = mean(case_results[["lAcc"]])
  )
  return( r1 )
}





#' @title Operating characteristics for sequential lot sampling
#'
#' @description
#' Generates Monte Carlo estimates of operating-characteristic summaries for
#' sequential lot sampling over one or more true defect probabilities.
#'
#' @param p Numeric vector of true defect probabilities at which to evaluate
#'   operating characteristics.
#' @param N Numeric scalar or vector giving the lot size.
#' @param preq Numeric requirement on the defect probability.
#' @param alpha Numeric posterior tail probability used for lower and upper
#'   interval summaries.
#' @param nLots Positive integer giving the number of lots per simulated case.
#' @param nRep Positive integer giving the number of Monte Carlo replications
#'   for each value of `p`.
#' @param control An object of class `saber_control`, typically created by
#'   [saber_control()].
#' @param sls_control An object of class `sls_control`, typically created by
#'   [sls_control()]. This controls the SLS evaluation method and policy
#'   settings.
#' @param seed Optional numeric seed for the random number generator. If
#'   omitted, a seed is generated from the local time.
#'
#' @return
#' A data frame with one row per replication for each value of `p` and columns:
#'
#' - `p`: True defect probability used for the simulated cases.
#' - `prop_accepted`: Proportion of units accepted.
#' - `prop_accepted_defective`: Proportion of accepted units that are defective.
#' - `number_defective_seen`: Number of defective units observed in sampling.
#' - `number_sampled`: Number of units sampled.
#' - `fraction_sampled`: Fraction of units sampled.
#' - `number_accepted`: Number of units accepted.
#' - `prop_lots_accepted`: Proportion of lots accepted.
#'
#' @seealso [sls_sim_case()], [sls_sim_timeline()], [sls_control()],
#'   [saber_control()], [rsaber_future_plan()]
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture %dofuture%
#' 
#' @examples
#' \dontrun{
#' rsaber_future_plan("sequential")
#'
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
#' oc <- sls_sim_oc(
#'   p = c(0.02, 0.05, 0.10),
#'   N = 100,
#'   preq = 0.10,
#'   alpha = 0.05,
#'   nLots = 20,
#'   nRep = 100,
#'   control = control,
#'   sls_control = sls_ctrl,
#'   seed = 123
#' )
#'
#' head(oc)
#' } 
#' 
#' @export
#' 
#'  
sls_sim_oc <- function( p, N, preq, alpha, nLots, nRep, control, sls_control, seed ){
  
  # sim_control <- validate_sim_control( sim_control )
  control <- validate_saber_control(control)
  sls_control <- validate_sls_control(sls_control)
  
  # Initial setup/calculations
  dfoc <- matrix( NA, nrow=length(p)*nRep, ncol=8 )
  
  if( missing(seed) ){
    seed <- as.integer(as.numeric(format(Sys.time(), "%S%M%H%d%m%Y")) %% .Machine$integer.max)
  }
  # set.seed(seed)
  seed_from_index <- function(jj, ii, base_seed = 1L) {
    as.integer((base_seed + 100000L * jj + ii) %% (.Machine$integer.max - 1L) + 1L)
  }
  
  for( jj in 1:length(p) ){
    pj <- p[jj]
    results <- foreach( ii = seq_len(nRep), .combine=rbind, 
                        .options.future = list(packages = sls_control[["packages"]], seed=TRUE) 
    ) %dofuture% {
      sls_sim_case_to_oc(
        sls_sim_case( 
          p=pj, N=N, preq=preq, alpha=alpha, nLots=nLots, 
          control = control,
          sls_control=sls_control, 
          seed = seed_from_index(jj,ii,seed)
        )
      )
    }
    
    idx1 <- (jj-1)*nRep + 1
    idx2 <- jj*nRep
    dfoc[ idx1:idx2, ] <- cbind( rep(pj, nrow(results)), results )
  }
  
  dfoc_comb <- as.data.frame(dfoc)
  colnames(dfoc_comb) <- c("p",
                           "prop_accepted",
                           "prop_accepted_defective", 
                           "number_defective_seen",
                           "number_sampled",
                           "fraction_sampled",
                           "number_accepted",
                           "prop_lots_accepted")
  return(dfoc_comb)
}





#' @title Timeline simulation for sequential lot sampling
#'
#' @description
#' Generates Monte Carlo simulation results for visualizing sequential
#' lot-by-lot behavior under a selected SLS evaluation rule.
#'
#' @param p Numeric scalar or vector giving the true defect probability used to
#'   simulate lots. If a vector is supplied, it should be compatible with
#'   `nLots`.
#' @param N Numeric scalar or vector giving the lot size.
#' @param preq Numeric requirement on the defect probability.
#' @param alpha Numeric posterior tail probability used for lower and upper
#'   interval summaries.
#' @param nLots Positive integer giving the number of lots per simulated case.
#' @param nRep Positive integer giving the number of Monte Carlo replications.
#' @param control An object of class `saber_control`, typically created by
#'   [saber_control()].
#' @param sls_control An object of class `sls_control`, typically created by
#'   [sls_control()]. This controls the SLS evaluation method and policy
#'   settings.
#' @param seed Optional numeric seed for the random number generator. If
#'   omitted, a seed is generated from the local time.
#'
#' @return
#' A data frame with one row per lot per replication and columns:
#'
#' - `rpl`: Replication index.
#' - `lot`: Lot index.
#' - `N`: Lot size.
#' - `K`: True simulated number of defective units in the lot.
#' - `k`: Number of defective units sampled from the lot.
#' - `n`: Number of units sampled from the lot.
#' - `a`: First posterior beta parameter after evaluating the lot.
#' - `b`: Second posterior beta parameter after evaluating the lot.
#' - `dAcc`: Number of defective units accepted from the lot.
#' - `nAcc`: Number of units accepted from the lot.
#' - `lAcc`: Lot acceptance indicator, where 1 indicates accepted and 0
#'   indicates rejected.
#' - `pest`: Posterior median of the defect probability.
#' - `plcl`: Posterior lower credible bound of the defect probability.
#' - `pucl`: Posterior upper credible bound of the defect probability.
#' - `Xest`: Posterior predictive median of the number of defectives in the
#'   unsampled portion of the lot.
#' - `Xlcl`: Posterior predictive lower credible bound for the number of
#'   defectives in the unsampled portion of the lot.
#' - `Xucl`: Posterior predictive upper credible bound for the number of
#'   defectives in the unsampled portion of the lot.
#'
#' @seealso [sls_sim_case()], [sls_sim_oc()], [sls_control()],
#'   [saber_control()], [rsaber_future_plan()]
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doFuture %dofuture%
#' 
#' 
#' @examples
#' \dontrun{
#' rsaber_future_plan("sequential")
#'
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
#' timeline <- sls_sim_timeline(
#'   p = 0.05,
#'   N = 100,
#'   preq = 0.10,
#'   alpha = 0.05,
#'   nLots = 20,
#'   nRep = 100,
#'   control = control,
#'   sls_control = sls_ctrl,
#'   seed = 123
#' )
#'
#' head(timeline)
#' } 
#' 
#' @export
#' 
sls_sim_timeline <- function( p, N, preq, alpha, nLots, nRep, control, sls_control, seed ){
  
  # sim_control <- validate_sim_control( sim_control )
  control <- validate_saber_control(control)
  sls_control <- validate_sls_control(sls_control)
  
  if( missing(seed) ){
    seed <- as.integer(as.numeric(format(Sys.time(), "%S%M%H%d%m%Y")) %% .Machine$integer.max)
  }
  # set.seed(seed)
  seed_from_index <- function(ii, base_seed = 1L) {
    as.integer((base_seed + ii) %% (.Machine$integer.max - 1L) + 1L)
  }
  
  results <- foreach( ii = seq_len(nRep), .combine=rbind, 
                      .options.future = list(packages = sls_control[["packages"]], seed=TRUE)
  ) %dofuture% {
    case_out <- sls_sim_case( 
      p=p, N=N, preq=preq, alpha=alpha, nLots=nLots, 
      control = control,
      sls_control=sls_control, 
      seed = seed_from_index( ii, seed )
    )
    
    est_mat <- case_out 
    rep_vec <- rep( ii, nLots )
    cbind( rep_vec, est_mat )
  }
  
  results <- data.frame(results)
  colnames( results ) <- c( "rpl", "lot", 
                            "N", "K", "k", "n", "a", "b", 
                            "dAcc", "nAcc", "lAcc" , 
                            "pest", "plcl", "pucl", 
                            "Xest", "Xlcl", "Xucl")
  return( results )
}
















