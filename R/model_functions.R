####https://cran.r-project.org/web/packages/rstantools/vignettes/minimal-rstan-package.html

#' #' Run Model
#' #'
#' stan_surv_model = function(inputlist, ...){
#'
#'   # mod1 = stan(
#'   #   file = "BayesSurv STANcode v15.stan",
#'   #   data = mod1.list,
#'   #   chains = 4,
#'   #   cores = 4,
#'   #   iter = 2000,
#'   #   init = list(list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#'   #               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#'   #               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#'   #               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)))
#'   #
#'   # )
#'
#'   mod = rstan::stan(
#'     file = stanmodels$BayesSurv_STANcode_v15,
#'     data = inputlist
#'   )
#'   return(mod)
#' }
#'
#'
#'
#'
#'
#'
