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
  if(dims(sc.l) == 1){
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

