



#' @title Set the number of decimal places for a value
#' 
#' @description Forces a particular decimal representation of a numeric value. This
#' converts the result to character.
#' 
#' 
#' @param x (Numeric) the value to be set
#' @param k (numeric integer) the number of decimal places
#' 
#' @return Returns a string value.
#' 
#' 
specify_decimal <- function(x, k) trimws( format(round(x, k), nsmall=k) )





#' @title Create a mesh of parameters
#' 
#' @description Takes a list of input parameters and creates a mesh in one of 
#' several manners.
#' 
#' 
#' @param mesh (string) a character string denoting the type of meshing to be done.
#' @param params (list) a named list of parameters to be meshed.
#' 
#' @return Returns a data frame.
#' 
#' @details Argument `mesh` has two options:
#' - "crossing" will create all combinations (using `expand.grid`)
#' - "recycle" will cycle through the elements, repeating as needed to fill the grid.
#'   Users should carefully inspect the resulting mesh if this option is used.
#' 
parameter_mesh <- function( mesh, params ){
  
  # Input validation
  if( !(mesh %in% c("crossing", "recycle")) ){
    stop( "Argument `mesh` should be 'crossing' or 'recycle'" )
  }
  if( is.null(names(params)) ){
    stop( "Argument `params` must be a named list" )
  }
  
  # Make the dataframe
  if( tolower(mesh) == "crossing" ){
    dat00 <- expand.grid( params )

  } else if( tolower(mesh) == "recycle" ){
    # If we can get this to display only once per session, that would be great.
    # Maybe a global variable that gets initialized on load of package, and then overwritten?
    # There's a way with rlang, but that's a dependency, and I'd prefer to avoid those as much as possible.
    # warn_message <- paste0(
    #   "Recycling inputs (x, n, N, alpha), inspect output to ensure cases are correct."
    # )
    # rlang::warn( paste0(warn_message, collapse="\n"), .frequency="once" )
    
    arg_lengths <- unname(unlist( lapply( params, FUN=length ) ))
    max_length  <- max(arg_lengths)
    arg_reps    <- ceiling(max_length/arg_lengths)
    for( ii in 1:length(params) ){
      params[[ii]] <- rep( params[[ii]], arg_reps[ii] )[1:max_length]
    }
    dat00 <- data.frame(params)
  } else{
    stop( "Unknown error attempting to create the parameter mesh" )
  }
  
  # Return the mesh
  return( dat00 )
}









#' Configure the future plan used by rsaber simulations
#'
#' Sets up the \pkg{future} backend used by \pkg{doFuture}/\pkg{foreach} when
#' running rsaber simulation functions that use \code{\%dofuture\%}.
#'
#' This function is intended to be called once per R session (or whenever you
#' want to change parallel settings). It modifies the global future plan; i.e.,
#' it affects other code in the same R session that uses \pkg{future}.
#'
#' @param strategy (character) Future strategy to use. Options are
#'   \code{"sequential"} (no parallelism), \code{"multisession"} (separate R
#'   sessions; recommended on Windows), and \code{"multicore"} (forking;
#'   available on Linux/macOS only).
#' @param workers (integer) Number of parallel workers to use for strategies
#'   that support it. Defaults to \code{future::availableCores()}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects of
#'   registering \pkg{doFuture} and setting the global future plan.
#'
#' @details
#' Internally, this function calls \code{doFuture::registerDoFuture()} and then
#' \code{future::plan()} using the selected strategy.
#'
#' @seealso
#' \code{\link[future]{plan}}, \code{\link[future]{multisession}},
#' \code{\link[future]{multicore}}, \code{\link[future]{sequential}},
#' \code{\link[doFuture]{registerDoFuture}}
#'
#' @examples
#' 
#' \dontrun{
#' # Run sequentially (no parallelism)
#' rsaber_future_plan("sequential")
#'
#' # Run in parallel using separate R sessions
#' rsaber_future_plan("multisession", workers = 4)
#' }
#' 
#' @importFrom future availableCores plan sequential multisession multicore
#'
#' @export
#' 
rsaber_future_plan <- function(strategy = c("sequential", "multisession", "multicore"),
                               workers = future::availableCores()) {
  strategy <- match.arg(strategy)
  doFuture::registerDoFuture()
  if (strategy == "sequential") {
    future::plan(future::sequential)
  } else if (strategy == "multisession") {
    future::plan(future::multisession, workers = workers)
  } else if (strategy == "multicore") {
    if (.Platform$OS.type == "windows") {
      stop("The 'multicore' future strategy is not supported on Windows. Use 'multisession' instead.")
    }
    future::plan(future::multicore, workers = workers)
  }
  invisible(NULL)
}











