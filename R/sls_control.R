#' @title Construct SLS simulation/control settings
#'
#' @description
#' Constructs a control object for series-of-lots simulation and evaluation
#' functions. These settings define the sampling/evaluation policy used by
#' SLS simulation functions.
#'
#' @param x0 Nonnegative integer giving the maximum allowed sampled defectives
#'   for lot acceptance.
#' @param nn_min Nonnegative integer or `NULL`. Minimum or fixed sample size,
#'   depending on the evaluation method.
#' @param kappa Numeric value in `(0, 1)` or `NULL`. Tail probability used for
#'   automatic minimum sample-size selection when applicable.
#' @param skip_max Positive integer or `NULL`. Maximum number of lots to skip
#'   for skip-lot methods.
#' @param skip_prob Numeric value in `[0, 1]` or `NULL`. Probability of skipping
#'   an eligible lot for skip-lot methods.
#' @param eval_function Evaluation function or character alias. Character
#'   options are `"saber"`, `"saber_fixed"`, `"bayesian_indp"`,
#'   `"skiplot_indp"`, and `"skiplot_saber"`.
#' @param estimation_function Function used to compute posterior estimates from
#'   simulated SLS results.
#' @param packages Character vector of package names required for parallel
#'   execution, especially when custom evaluation or estimation functions use
#'   functions from other packages.
#'
#' @return
#' An object of class `sls_control`, represented as a list of SLS
#' simulation/control settings.
#'
#' @seealso [saber_control()], [sls_sim_case()], [sls_sim_oc()],
#'   [sls_sim_timeline()]
#'
#' @examples
#' sls_ctrl <- sls_control(
#'   x0 = 0,
#'   nn_min = 5,
#'   eval_function = "saber"
#' )
#'
#' sls_ctrl
#'
#' @export
sls_control <- function(x0 = 0,
                        nn_min = NULL,
                        kappa = NULL,
                        skip_max = NULL,
                        skip_prob = NULL,
                        eval_function = c("saber", "saber_fixed", "bayesian_indp",
                                          "skiplot_indp", "skiplot_saber"),
                        estimation_function = sls_sim_estimate,
                        packages = "rsaber") {
  
  if (is.character(eval_function)) {
    eval_function <- match.arg(eval_function)
  }
  
  if (!is.numeric(x0) || length(x0) != 1 || x0 < 0 || x0 != floor(x0)) {
    stop("Argument 'x0' must be a nonnegative integer.", call. = FALSE)
  }
  
  if (!is.null(nn_min)) {
    if (!is.numeric(nn_min) || length(nn_min) != 1 || nn_min < 0 || nn_min != floor(nn_min)) {
      stop("Argument 'nn_min' must be a nonnegative integer or NULL.", call. = FALSE)
    }
    nn_min <- as.integer(nn_min)
  }
  
  if (!is.null(kappa)) {
    if (!is.numeric(kappa) || length(kappa) != 1 || kappa <= 0 || kappa >= 1) {
      stop("Argument 'kappa' must be in (0, 1), or NULL.", call. = FALSE)
    }
  }
  
  if (!is.null(skip_max)) {
    if (!is.numeric(skip_max) || length(skip_max) != 1 || skip_max < 1 || skip_max != floor(skip_max)) {
      stop("Argument 'skip_max' must be a positive integer or NULL.", call. = FALSE)
    }
    skip_max <- as.integer(skip_max)
  }
  
  if (!is.null(skip_prob)) {
    if (!is.numeric(skip_prob) || length(skip_prob) != 1 || skip_prob < 0 || skip_prob > 1) {
      stop("Argument 'skip_prob' must be in [0, 1], or NULL.", call. = FALSE)
    }
  }
  
  if (!is.function(estimation_function)) {
    stop("Argument 'estimation_function' must be a function.", call. = FALSE)
  }
  
  packages <- validate_sls_packages(packages)
  
  structure(
    list(
      x0 = as.integer(x0),
      nn_min = nn_min,
      kappa = kappa,
      skip_max = skip_max,
      skip_prob = skip_prob,
      eval_function = eval_function,
      estimation_function = estimation_function,
      packages = packages
    ),
    class = "sls_control"
  )
}





#' Validate SLS simulation/control settings
#'
#' @keywords internal
validate_sls_control <- function(control) {
  
  if (is.null(control)) {
    stop("Argument 'sls_control' must be supplied.", call. = FALSE)
  }
  
  if (!is.list(control)) {
    stop("Argument 'sls_control' must be a list, preferably created by sls_control().", call. = FALSE)
  }
  
  required_names <- c(
    "x0", "nn_min", "kappa",
    "skip_max", "skip_prob",
    "eval_function", "estimation_function",
    "packages"
  )
  
  missing_names <- setdiff(required_names, names(control))
  
  if (length(missing_names) > 0) {
    stop(
      "Argument 'sls_control' is missing required elements: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  ## Validate package vector first, because function checks need it.
  control[["packages"]] <- validate_sls_packages(control[["packages"]])
  
  ## Validate eval_function.
  if (is.character(control[["eval_function"]])) {
    control[["eval_function"]] <- match.arg(
      control[["eval_function"]],
      choices = c(
        "saber",
        "saber_fixed",
        "bayesian_indp",
        "skiplot_indp",
        "skiplot_saber"
      )
    )
  } else if (is.function(control[["eval_function"]])) {
    check_sls_control_function(
      f = control[["eval_function"]],
      packages = control[["packages"]],
      name = "eval_function"
    )
  } else {
    stop(
      "sls_control[['eval_function']] must be a character alias or function.",
      call. = FALSE
    )
  }
  
  
  ## Validate estimation_function.
  if (!is.function(control[["estimation_function"]])) {
    stop("sls_control[['estimation_function']] must be a function.", call. = FALSE)
  }
  
  check_sls_control_function(
    f = control[["estimation_function"]],
    packages = control[["packages"]],
    name = "estimation_function"
  )  
  
  
  if (!inherits(control, "sls_control")) {
    class(control) <- c("sls_control", class(control))
  }
  
  return(control)
}







#' Validate package list for SLS parallel execution
#'
#' @keywords internal
validate_sls_packages <- function(packages) {
  
  if (is.null(packages)) {
    packages <- "rsaber"
  }
  
  if (!is.character(packages) || length(packages) < 1) {
    stop("sls_control[['packages']] must be a character vector of package names.", call. = FALSE)
  }
  
  if (any(is.na(packages)) || any(packages == "")) {
    stop("sls_control[['packages']] contains missing or empty package names.", call. = FALSE)
  }
  
  packages <- unique(packages)
  
  return(packages)
}


#' Check package dependencies inside a control function
#'
#' @keywords internal
check_sls_control_function <- function(f, packages = "rsaber", name = NULL) {
  
  if (!is.function(f)) {
    return(invisible(TRUE))
  }
  
  packages <- validate_sls_packages(packages)
  
  txt <- paste0(deparse(body(f), width.cutoff = 500L), collapse = "\n")
  
  ## 1) Disallow explicit package dependencies via pkg::fun or pkg:::fun
  m <- gregexpr(
    "\\b([[:alnum:].]+)\\s*:::{0,1}\\s*[[:alnum:]_.]+\\b",
    txt,
    perl = TRUE
  )
  
  hits <- regmatches(txt, m)[[1]]
  
  if (length(hits) > 0 && !identical(hits, character(0))) {
    pkgs <- sub("^([[:alnum:].]+)\\s*:::{0,1}.*$", "\\1", hits)
    pkgs <- unique(pkgs)
    
    bad_pkgs <- setdiff(pkgs, packages)
    
    if (length(bad_pkgs) > 0) {
      fn_label <- if (is.null(name)) {
        "A function in 'sls_control'"
      } else {
        paste0("sls_control[['", name, "']]")
      }
      
      stop(
        fn_label,
        " references packages via '::' or ':::' that are not listed in ",
        "sls_control[['packages']]: ",
        paste(bad_pkgs, collapse = ", "),
        ". Add these package names to sls_control[['packages']].",
        call. = FALSE
      )
    }
  }
  
  ## 2) Flag common tidyverse verbs used without qualification
  tidyverse_verbs <- c(
    "mutate", "group_by", "summarize", "summarise",
    "tibble", "filter", "select", "arrange"
  )
  
  used_verbs <- tidyverse_verbs[
    vapply(
      tidyverse_verbs,
      function(v) grepl(paste0("\\b", v, "\\s*\\("), txt, perl = TRUE),
      logical(1)
    )
  ]
  
  if (length(used_verbs) > 0) {
    fn_label <- if (is.null(name)) {
      "A function in 'sls_control'"
    } else {
      paste0("sls_control[['", name, "']]")
    }
    
    warning(
      fn_label,
      " appears to call common tidyverse functions without package qualification: ",
      paste(used_verbs, collapse = ", "),
      ". If this is intentional, use explicit qualification, e.g. dplyr::mutate(), ",
      "and add the package to sls_control[['packages']].",
      call. = FALSE
    )
  }
  
  ## 3) Discourage loading packages inside the function
  if (grepl("\\b(library|require)\\s*\\(", txt, perl = TRUE)) {
    fn_label <- if (is.null(name)) {
      "A function in 'sls_control'"
    } else {
      paste0("sls_control[['", name, "']]")
    }
    
    warning(
      fn_label,
      " appears to call library() or require(). This is discouraged for ",
      "parallel execution. Prefer explicit pkg::fun() calls and list required ",
      "packages in sls_control[['packages']].",
      call. = FALSE
    )
  }
  
  return(invisible(TRUE))
}

