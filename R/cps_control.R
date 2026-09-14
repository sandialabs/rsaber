



#' Construct CPS simulation control settings
#'
#' @param n_init Initial number of units sampled.
#' @param nSim Number of Monte Carlo simulations for defect-count estimation.
#' @param n_forecast Number of future units forecast at a time.
#' @param nSim_forecast Number of forecast simulation paths.
#' @param forecast_type Forecasting mode.
#' 
#' 
#' @examples
#' 
#' cps_control( 
#'   n_init=100,
#'   nSim = 1000,
#'   n_forecast = 10,
#'   nSim_forecast = 100,
#'   forecast_type = "estimated"
#' )
#' 
#' 
#' @export
cps_control <- function(n_init,
                        nSim = 1000,
                        n_forecast,
                        nSim_forecast = 100,
                        forecast_type = c("estimated", "optimistic", "pessimistic")) {
  
  forecast_type <- match.arg(forecast_type)
  
  if (!is.numeric(n_init) || length(n_init) != 1 || n_init < 1 || n_init != floor(n_init)) {
    stop("Argument 'n_init' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.numeric(nSim) || length(nSim) != 1 || nSim < 1 || nSim != floor(nSim)) {
    stop("Argument 'nSim' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.numeric(n_forecast) || length(n_forecast) != 1 || n_forecast < 1 || n_forecast != floor(n_forecast)) {
    stop("Argument 'n_forecast' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.numeric(nSim_forecast) || length(nSim_forecast) != 1 || nSim_forecast < 1 || nSim_forecast != floor(nSim_forecast)) {
    stop("Argument 'nSim_forecast' must be a positive integer.", call. = FALSE)
  }
  
  structure(
    list(
      n_init = as.integer(n_init),
      nSim = as.integer(nSim),
      n_forecast = as.integer(n_forecast),
      nSim_forecast = as.integer(nSim_forecast),
      forecast_type = forecast_type
    ),
    class = "cps_control"
  )
}


#' Validate CPS simulation control settings
#'
#' @keywords internal
validate_cps_control <- function(control) {
  
  if (is.null(control)) {
    stop("Argument 'cps_control' must be supplied.", call. = FALSE)
  }
  
  if (!is.list(control)) {
    stop("Argument 'cps_control' must be a list, preferably created by cps_control().", call. = FALSE)
  }
  
  required_names <- c(
    "n_init",
    "nSim",
    "n_forecast",
    "nSim_forecast",
    "forecast_type"
  )
  
  missing_names <- setdiff(required_names, names(control))
  
  if (length(missing_names) > 0) {
    stop(
      "Argument 'cps_control' is missing required elements: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  forecast_type <- control[["forecast_type"]]
  control[["forecast_type"]] <- match.arg(
    forecast_type,
    choices = c("estimated", "optimistic", "pessimistic")
  )
  
  int_names <- c("n_init", "nSim", "n_forecast", "nSim_forecast")
  
  for (nm in int_names) {
    val <- control[[nm]]
    
    if (!is.numeric(val) || length(val) != 1 || val < 1 || val != floor(val)) {
      stop(
        "cps_control[['", nm, "']] must be a positive integer.",
        call. = FALSE
      )
    }
    
    control[[nm]] <- as.integer(val)
  }
  
  if (!inherits(control, "cps_control")) {
    class(control) <- c("cps_control", class(control))
  }
  
  control
}

