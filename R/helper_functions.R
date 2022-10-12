#' Gompertez Survival Function
#'
#'@export
gomp_surv = function(x, a = NULL, b, M = NULL, by.M = FALSE){
  if(by.M){ # by modal age
    # https://www.monicaalexander.com/posts/2018-02-15-gompertz/
    Sx = exp(-exp(-b*M)*(exp(b*x)-1))
  } else { # by baseline mortality
    Sx = exp(-1*(a/b)*(exp(b*x)-1))
  }
  return(Sx)
}

#' Gompertez Prob Distribution Function
#'
#'@export
gomp.pdf = function(x, x.range = seq(0:150) ,a, b){
  pd = numeric(length(x))
  for(i in 1:length(x)){
    lx = gomp_surv(x = x[i], a = a, b = b)
    sigma_lx = sum(gomp_surv(x = x.range, a = a, b = b))
    pd[i] = lx/sigma_lx
  }
  return(pd)
}

#' Get expected lifespan at a given age
#'
#'@export
get_ex = function(ageseq, sc.l, ex.at){
  if(dim(sc.l) == 1){
    i = which(abs(0.5-round(sc.l, 2)) == min(abs(0.5-round(sc.l, 2))))
    age.ex = ageseq[i]
  } else {
    i.1 = which(abs(0.5-round(sc.l[1,], 2)) == min(abs(0.5-round(sc.l[1,], 2))))
    i.2 = which(abs(0.5-round(sc.l[2,], 2)) == min(abs(0.5-round(sc.l[2,], 2))))
    age.ex = c(ageseq[i.1], ageseq[i.2])
  }
  return(age.ex)
}

#' Get age at which X prop of survival remains
#'
#'@export
analytical_age_of_lx = function(a, b, X = 0.5){
  # see surv at age alegbra v1.jpg
  x = (log(log(X)*(b/-a)+1))/b
  return(x)
}

#' Get Prior Shapes
#'
#'@export
get_gomp_prior_shapes = function(species = "odontocete"){

  ## Gets gompterez priors
  #import data as odontocete.gomp (withod changing wd) and then run this
  # usethis::use_data(odotocete.gomp)

  if(species == "odontocete"){
    real.gomp = odotocete.gomp
  }

  shape1_a = mean(real.gomp$initital) /
    ((stats::sd(real.gomp$initital)*2) * (stats::sd(real.gomp$initital)*2)) # 2
  shape2_a = (1-mean(real.gomp$initital)) /
    ((stats::sd(real.gomp$initital)*2) * (stats::sd(real.gomp$initital)*2)) #2
  shape1_b = mean(real.gomp$sense, na.rm = TRUE) /
    ((stats::sd(real.gomp$sense, na.rm = TRUE)*1.5) * (stats::sd(real.gomp$sense, na.rm = TRUE)*1.5)) #1.5
  shape2_b = (1- mean(real.gomp$sense, na.rm = TRUE)) /
    ((stats::sd(real.gomp$sense, na.rm = TRUE)*1.5) * (stats::sd(real.gomp$sense, na.rm = TRUE)*1.5))# 1.5
  comb = list(shape1_a = shape1_a,
              shape2_a = shape2_a,
              shape1_b = shape1_b,
              shape2_b = shape2_b
  )

  return(comb)
}

#' Generate intitial values
#'
#' @export
generate_inits = function(chains, input.list){

  if(length(unique(input.list$species_sex_vector))==1){
    inits = vector(chains, mode = "list")
    for(i in 1:chains){
      inits[[i]] =
        list(b = stats::runif(input.list$Nspecies, max = 0.05),
             a = stats::runif(input.list$Nspecies, max = 0.2),
             bbar = stats::runif(1, max = 0.03),
             abar = stats::runif(1, max = 0.2))
    }

  } else {
    inits = vector(chains, mode = "list")
    for(i in 1:chains){
      inits[[i]] =
        list(b = stats::runif(input.list$Nspecies, max = 0.05),
             a = stats::runif(input.list$Nspecies, max = 0.2),
             bbar = rep(stats::runif(1, max = 0.03),2),
             abar = rep(stats::runif(1, max = 0.2),2)
        )
    }
  }

  return(inits)
}


#' Get list of parameters actually active in the model
#'
#' @export
get_activepars = function(out, mod, input.list){
  activepars = rownames(out)[!stringr::str_detect(rownames(out), "true_age")]
  for(i in 1:input.list$Ndatasets){
    if(input.list$include_samplebias_error[i] == 0){
      activepars <- activepars[activepars != paste("s", "[", i, "]", sep = "")]
    }
    if(input.list$include_popchange_error[i] == 0){
      activepars <- activepars[activepars != paste("r", "[", i, "]", sep = "")]
    }
  }
  activepars = activepars[!(activepars %in% c("abar_theta", "bbar_theta"))]
  activepars = activepars[!stringr::str_detect(activepars, "rho")]
  return(activepars)
}

