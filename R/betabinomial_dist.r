



#' @title The Beta-Binomial Distribution
#' 
#' @description Density, distribution function, quantile function and random 
#' generation for the Beta-Binomial distribution.
#' 
#' @param y Vector of quantiles
#' @param Nn number of observations
#' @param a alpha parameter
#' @param b beta parameter
#' 
#' @details The Beta-Binomial distribution has the mass function
#' 
#' \deqn{
#' p(x) = \binom{Nn}{y}\dfrac{B(y+a, Nn-y+b)}{B(a,b)}
#' }
#' 
#' where \eqn{B()} is the Beta function. The syntax `Nn` is used so that
#' `rbetabinomial` can use `n` for the number of random variates to draw.
#' 
#' @return \code{dbetabinomial} gives the mass function (PMF), 
#' \code{pbetabinomial} gives the distribution function (CDF), 
#' \code{qbetabinomial} gives the quantile function, and 
#' \code{rbetabinomial} generates random deviates.
#' 
#' @examples
#' 
#' dbetabinomial( 1, 10, 1, 5 )
#' dbetabinomial( 0:10, 10, 1, 5 )
#' 
#' 
#' @export
#' 
dbetabinomial <- function( y, Nn, a, b ) {
  
  # Do some input validation
  arg_lengths <- c( length(y), length(Nn), length(a), length(b) )
  
  if(!(  (sd(arg_lengths)==0) | (sum(arg_lengths>1)==1)  )){
    stop("Either all arguments must have equal length, or only one argument may exceed length 1")
  }
  
  # Equalize the lengths of arguments
  max_len <- max( arg_lengths )
  if( length(y)==1 ){  y  <- rep( y,  max_len ) }
  if( length(Nn)==1 ){ Nn <- rep( Nn, max_len ) }
  if( length(a)==1 ){  a  <- rep( a,  max_len ) }
  if( length(b)==1 ){  b  <- rep( b,  max_len ) }
  
  if( any(y<0) | any(y>Nn) ){ stop("y must be in [0, Nn]") }
  
  combnt <- lchoose( Nn, y )
  bfun_num <- lbeta( y + a, Nn-y+b )
  bfun_den <- lbeta( a, b )
  
  pdf_val <- combnt + bfun_num - bfun_den
  
  return( exp(pdf_val) )
  
}



#' @title The Beta-Binomial Distribution
#' 
#' @rdname dbetabinomial
#' 
#' @examples
#' 
#' dbetabinomial( 1, 10, 1, 5 )
#' dbetabinomial( 0:10, 10, 1, 5 )
#' 
#' 
#' 
#' @export
#' 
pbetabinomial <- function( y, Nn, a, b ) {
  
  # Do some input validation
  arg_lengths <- c( length(y), length(Nn), length(a), length(b) )
  
  if(!(  (sd(arg_lengths)==0) | (sum(arg_lengths>1)==1)  )){
    stop("Either all arguments must have equal length, or only one argument may exceed length 1")
  }
  
  # Equalize the lengths of arguments
  max_len <- max( arg_lengths )
  if( length(y)==1 ){  y  <- rep( y,  max_len ) }
  if( length(Nn)==1 ){ Nn <- rep( Nn, max_len ) }
  if( length(a)==1 ){  a  <- rep( a,  max_len ) }
  if( length(b)==1 ){  b  <- rep( b,  max_len ) }
  
  if( any(y<0) | any(y>Nn) ){ stop("y must be in [0, Nn]") }
  
  cdf_val <- numeric( length=max_len )
  for( ii in 1:max_len ){
    # pdf_vals <- dbetabinomial( y=seq(0,y[ii],1), 
    #                            Nn = Nn[ii], a=a[ii], b=b[ii] )
    yseq <- seq(0,y[ii],1)
    combnt <- lchoose( Nn[ii], yseq )
    bfun_num <- lbeta( yseq + a[ii], Nn[ii]-yseq+b[ii] )
    bfun_den <- lbeta( a[ii], b[ii] )
    pdf_vals <- combnt + bfun_num - bfun_den
    cdf_val[ii] <- sum( exp(pdf_vals) )
  }
  return(cdf_val)
}




#' @title The Beta-Binomial Distribution
#' 
#' @rdname dbetabinomial
#' 
#' @param p Vector of probabilities
#' `round` (default), `ceiling`, or `floor`. If 
#' 
#' @details
#' The quantile function `qbetabinomial` uses an optimization which does not 
#' necessarily return an integer. The final result is rounded to the nearest 
#' integer. Investigations show this to be a very small rounding (on the order 
#' of 1e-4 or less).
#' 
#' 
#' @examples
#' 
#' qbetabinomial( 0.35, 10, 1, 5 )
#' qbetabinomial( seq(0,1,0.10), 10, 1, 5 )
#' 
#' 
#' @export
#' 
qbetabinomial <- function (p, Nn, a, b) {
  
  qbb_uniroot <- function(y, Nn, a, b, p ) {
    yseq <- seq(0, y, 1)
    sum(exp(
      lchoose(Nn, yseq) +
        lbeta(yseq + a , Nn - yseq + b) -
        lbeta(a, b)
    )) - p
  }
  
  # Do some input validation here?
  arg_lengths <- c( length(p), length(Nn), length(a), length(b) )
  
  if( any(p<0) | any(p>1) ){ stop("p must be in [0, 1]") }
  
  # Equalize the lengths of arguments
  max_len <- max( arg_lengths )
  if( length(p)==1 ){  p  <- rep( p,  max_len ) }
  if( length(Nn)==1 ){ Nn <- rep( Nn, max_len ) }
  if( length(a)==1 ){  a  <- rep( a,  max_len ) }
  if( length(b)==1 ){  b  <- rep( b,  max_len ) }
  
  qnt_val <- numeric( length=max_len )
  for (ii in 1:max_len) {
    if( p[ii]==0 ){
      qnt_val[ii] <- 0
    } else if( p[ii]==1 ){
      qnt_val[ii] <- Nn[ii]
    } else if( qbb_uniroot(0, Nn[ii], a[ii], b[ii], p[ii]) > 0 ){
      qnt_val[ii] <- 0
    } else{
      qnt_val[ii] <- uniroot( qbb_uniroot , interval=c(0, Nn[ii]), 
                              Nn[ii], a[ii], b[ii], p[ii])$root
    }
  }
  round(qnt_val)
}


#' @title The Beta-Binomial Distribution
#' 
#' @rdname dbetabinomial
#' 
#' @param n Number of random variates to draw
#' 
#' @details
#' `rbetabinomial` does not use the usual inverse of the CDF. For numerical 
#' simplicity and vectorization it uses the fact that the BetaBinomial is a 
#' Binomial random variable with unknown probability. Variates are obtained by 
#' sampling the Beta distribution with the appropriate prior, then sampling the 
#' Binomial distribution with that Beta variate.
#' 
#' @examples
#' 
#' rbetabinomial( 20, 10, 1, 5 )
#' 
#' 
#' @export
#' 
rbetabinomial <- function( n, Nn, a, b ){
  
  arg_lengths <- c( length(Nn), length(a), length(b) )
  arg_sensible <- all( arg_lengths %in% c(1,n) ) 
  if( !all(arg_sensible) ){
    stop( "Arguments Nn, a, and b should be length 1 or have length equal to n")
  }
  
  # Equalize the lengths of arguments
  if( length(Nn)==1 ){ Nn <- rep( Nn, n ) }
  if( length(a)==1 ){  a  <- rep( a,  n ) }
  if( length(b)==1 ){  b  <- rep( b,  n ) }
  
  # Generate random probability
  # Then binomial variate from the random probability
  random_probs <- rbeta( n=n, shape1=a, shape2=b )
  y_val <- rbinom(n = n, size = Nn, prob = random_probs)
  return( y_val )
}

# Old version
# rbetabinomial <- function( n, Nn, a, b ){
#   # For this function, assume that the parameters are all scalers
#   # Unusual to want to generate random variates with different parameters
#   
#   if( (length(Nn)>1) | (length(a)>1) | (length(b)>1) ){
#     stop("Parameters Nn, a, and b must be scaler for rbetabinomial")
#   }
#   
#   ppi <- runif( n, 0, 1 )
#   y_val <- rsaber::qbetabinomial( p=ppi, Nn=Nn, a=a, b=b )
#   return( y_val )
# }














