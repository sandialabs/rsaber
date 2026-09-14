

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "rsaber uses future and doFuture for parallel processing in simulations.\n",
    "Please configure future if you plan to use the simulation features.\n",
    "The function rsaber::rsaber_future_plan() can perform this configuration.\n",
    "Example:\n    rsaber_future_plan('multisession', workers = 4)"
  )
}


# Prevent following Note in the CRAN check: 
#     no visible binding for global variable 'ii'
# This gets raised at least in two functions using foreach()
#   - sls_sim_oc
#   - sls_sim_timeline
utils::globalVariables(c("ii", "jj"))