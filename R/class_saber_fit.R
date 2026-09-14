


#' Construct a SABER fit object
#'
#' Internal constructor for data-frame subclasses used to represent fitted
#' SABER objects.
#'
#' @param x A data frame containing the fitted trajectory.
#' @param subclass Character scalar giving the SABER subclass, such as
#'   `"saber_cps"` or `"saber_sls"`.
#' @param control A `saber_control` object or compatible list.
#' @param alpha Posterior tail probability used for interval summaries.
#' @param call Original matched call, if available.
#'
#' @return A data frame with classes `subclass`, `"saber_fit"`, and
#'   `"data.frame"`, with SABER metadata stored as attributes.
#'
#' @keywords internal
#' 
new_saber_fit <- function(x,
                          subclass,
                          control,
                          alpha = NULL,
                          call = NULL) {
  stopifnot(is.data.frame(x))
  stopifnot(is.character(subclass), length(subclass) == 1)
  
  control <- validate_saber_control(control)
  attr(x, "saber_control") <- control
  
  attr(x, "saber_alpha") <- alpha
  # attr(x, "saber_frame") <- frame
  attr(x, "saber_frame") <- control[["frame"]]
  attr(x, "saber_call") <- call
  
  class(x) <- c(subclass, "saber_fit", "data.frame")
  
  x
}




#' Check whether an object is a SABER fit
#'
#' @keywords internal
is_saber_fit <- function(x) {
  inherits(x, "saber_fit")
}


#' Extract SABER control metadata
#'
#' @keywords internal
get_saber_control <- function(x) {
  attr(x, "saber_control", exact = TRUE)
}


#' Extract SABER alpha metadata
#'
#' @keywords internal
get_saber_alpha <- function(x) {
  attr(x, "saber_alpha", exact = TRUE)
}


#' Extract SABER frame metadata
#'
#' @keywords internal
get_saber_frame <- function(x) {
  attr(x, "saber_frame", exact = TRUE)
}


#' Extract original SABER fit call
#'
#' @keywords internal
get_saber_call <- function(x) {
  attr(x, "saber_call", exact = TRUE)
}


#' Validate a SABER fit object
#'
#' @keywords internal
validate_saber_fit <- function(x,
                               subclass = NULL,
                               required_cols = NULL,
                               require_control = TRUE) {
  
  if (!inherits(x, "saber_fit")) {
    stop("Object must inherit from class 'saber_fit'.", call. = FALSE)
  }
  
  if (!is.null(subclass) && !inherits(x, subclass)) {
    stop(
      "Object must inherit from class '", subclass, "'.",
      call. = FALSE
    )
  }
  
  if (is.null(required_cols)) {
    required_cols <- c("tidx", "x", "n", "N", "a", "b", "pest", "plcl", "pucl")
  }
  
  missing_cols <- setdiff(required_cols, names(x))
  
  if (length(missing_cols) > 0) {
    stop(
      "SABER fit object is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (require_control) {
    control <- get_saber_control(x)
    
    if (!is.list(control)) {
      stop(
        "SABER fit object is missing attribute 'saber_control'.",
        call. = FALSE
      )
    }
    
    required_control <- c("prior", "lambda", "forget")
    missing_control <- setdiff(required_control, names(control))
    
    if (length(missing_control) > 0) {
      stop(
        "SABER fit object has incomplete 'saber_control' metadata. Missing: ",
        paste(missing_control, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  invisible(TRUE)
}



#' Validate CPS frame
#'
#' @keywords internal
match_saber_frame <- function(frame) {
  match.arg(frame, choices = c("total", "unsampled", "accepted"))
}









