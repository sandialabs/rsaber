##
## This is the PACKAGE documentation
##

#' R package for the SABER method
#' 
#' @docType package
#' @name rsaber-package
#' @rdname rsaber-package
#' 
#' @description
#' Sequential Adaptive Bayesian Estimation for Reliability (SABER) is a 
#' Bayesian model for binary responses over time. SABER is intended primarily 
#' for two goals: Estimate defect rates over time; and adaptively set sampling 
#' rates (e.g., for acceptance sampling).
#' 
#' 
#' 
#' @details
#' The model assumes a Binomial likelihood and conjugate Beta prior. For finite 
#' populations (such as lots/batches), a superpopulation framework is used. A 
#' weighting function permits the model to "forget" older data and prioritize 
#' recent information. Functions to conduct simulation studies to parameterize 
#' the forgetting function.
#' 
#' Two main models are provided:
#' 
#' - SABER for Sequential Lot Sampling: For when data are considered in batches 
#' or lots over time.
#' - SABER for Continuous Production Sampling: For when units are produced or 
#' inspected one-at-a-time.
#' 
#' The fundamental idea of each is that the posterior is defined by 
#' Beta(\eqn{\alpha_{t}}{at}, \eqn{\beta_{t}}{bt}), where the parameters of the 
#' posterior are:
#' 
#' \deqn{ \alpha_{t} = \alpha_{0} + \sum_{i=1}^{t} \phi_{i}x_{i} }{ at = a0 + SUM phi_i xi}
#' 
#' \deqn{ \beta_{t} = \beta_{0} + \sum_{i=1}^{t} \phi_{i}(n_{i} - x_{i}) }{ bt = b0 + SUM phi_i(ni - xi)}
#' 
#' Defining \eqn{K} as the total number of defective units in a batch (for the 
#' case of finite populations), the posterior predictive distribution of 
#' \eqn{U=K-x} is also used. Estimation of the defect rate for finite populations 
#' can leverage this quantity in several different frames of references based on 
#' the nature of testing.
#' 
#' @references
#' Under review
#' 
#' 
#' @import methods
#' @import stats
#' 
"_PACKAGE"


