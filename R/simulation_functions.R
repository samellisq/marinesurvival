#' Simulate static survival sample
#'
#'@export
run_simulation = function(
  sim.in = NULL,
  include_r = TRUE,
  include_e = TRUE,
  include_s = TRUE,
  sample_size = 200,
  age.maturity = 0,
  show.plot = TRUE,
  max.age = 150
  ){

  real.gmps = get_gomp_prior_shapes()

  if(is.null(sim.in)){
    #TODO: need a function to check sim.in input here

    sim.in =
      data.frame(
        # a = stats::runif(1,0.01, 0.04), # baseline mortality taken from real examples
        # b = stats::runif(1,0.04, 0.1), # ageing
        a = rbeta(1, real.gmps$shape1_a, real.gmps$shape2_a), # baseline mortality taken from real examples
        b = rbeta(1, real.gmps$shape1_b, real.gmps$shape2_b), # ageing
        r = stats::runif(1, -0.5, 0.5), # 0 = no change
        e = 0.05, # standard deviation in of age estimation error
        s =  stats::runif(1, 0.5, 2) # change in probability of being sampled above age S
        # ageS = 25
      )

    all.gomp.surv = gomp_surv(x = 0:1000, a = sim.in$a, b = sim.in$b)
    real.max.age = min(which(all.gomp.surv <0.0001))
    sim.in$ageS = sample(floor(real.max.age/2):real.max.age, 1)


    # print(sim.in)
  }
  else{
    heads.in = names(sim.in)

    if(!("a" %in% heads.in)){
      sim.in$a =  rbeta(1, real.gmps$shape1_a, real.gmps$shape2_a)
    }
    if(!("b" %in% heads.in)){
      sim.in$b =  rbeta(1, real.gmps$shape1_b, real.gmps$shape2_b)
    }
    if(!("r" %in% heads.in)){
      sim.in$r = stats::runif(1, -0.5, 0.5)
    }
    if(!("e" %in% heads.in)){
      sim.in$e = 0.05
    }
    if(!("s" %in% heads.in)){
      sim.in$s = stats::runif(1, 0.5, 2)
    }
    if(!("ageS" %in% heads.in)){
      all.gomp.surv = gomp_surv(x = 0:1000, a = sim.in$a, b = sim.in$b)
      # print(sim.in)
      real.max.age = min(which(all.gomp.surv <0.0001))
      # print(real.max.age)
      sim.in$ageS = sample(floor(real.max.age/2):real.max.age, 1) # choose an age in the second half of life
    }

  }

  maxage = max.age
  sampling.probs =
    data.frame(
      age = seq(age.maturity, maxage, 1),
      lx.pdf = gomp.pdf(x = seq(age.maturity, maxage, 1),
                        x.range = seq(age.maturity, maxage*1.5,1),
                        a = sim.in$a, b = sim.in$b
                        ),
      r.weight = 1,
      s.weight = 1
      )
  # sampling.probs$obs.p = sampling.probs$lx.pdf

  #used for rho below but also to get xlim for plots
  all.gomp.surv = gomp_surv(x = 0:1000, a = sim.in$a, b = sim.in$b)
  real.max.age = min(which(all.gomp.surv <0.0001)) # becasue =0 is just limited by rs float point
  if(include_r){

    rho = sim.in$r / (0.5 * real.max.age)

    sampling.probs$r.weight = (1-rho)^seq(age.maturity, maxage, 1)
    # print(sampling.probs)
  }

  if(include_s){
    sampling.probs$s.weight = ifelse(sampling.probs$age >=sim.in$ageS,
                                     sim.in$s,
                                     1
    )



  }


  # if(include_r == TRUE & include_s == FALSE){
    # sampling.probs$obs.p = sampling.probs$lx.pdf*sampling.probs$r.weight
  # }
  # if(include_r == FALSE & include_s == TRUE){
    # print(sampling.probs)
    # sampling.probs$obs.p = sampling.probs$lx.pdf*sampling.probs$s.weight
  # }
  # if(include_r == TRUE & include_s == TRUE){
    sampling.probs$obs.p = sampling.probs$lx.pdf*sampling.probs$r.weight*sampling.probs$s.weight
  # }

  # print(sampling.probs)


  sampling.probs = dplyr::relocate(sampling.probs,obs.p,.after = last_col())
  gplot =
    ggplot2::ggplot(tidyr::pivot_longer(filter(sampling.probs, age <=real.max.age), cols = lx.pdf:obs.p, names_to = "p"),
                    ggplot2::aes(x = age, y = value, colour = p))+
    ggplot2::geom_line()+
    # ggplot2::xlim(c(age.maturity, real.max.age))+
    ggplot2::facet_wrap(.~p, ncol = 2, scales = "free_y")


  SAMPLE_SIZE = sample_size
  sample.cohort = sample(sampling.probs$age,
                         SAMPLE_SIZE,
                         prob = sampling.probs$obs.p,
                         replace = TRUE)
  age = sample.cohort
  sim = data.frame(age , surv = 1)

  gplot2 =  ggplot2::ggplot(sim, ggplot2::aes(age))+ggplot2::geom_histogram(binwidth = 1)+ggplot2::xlim(age.maturity, real.max.age)

  # sim.in$r = sim.in$rho*(max(sim)/2)

  if(include_e){
    sim$age.noerror = sim$age
    sim$e.sd = sim$age.noerror*sim.in$e
    for(i in 1:length(sim$age)){
      sim$age[i] = floor(stats::rnorm(1, sim$age[i]+0.5, sim$e.sd[i]))
    }
    # }

    sim = dplyr::filter(sim, age >=0)
  }

  if(show.plot){
    gridExtra::grid.arrange(gplot, gplot2, ncol = 1)
  }
  return(list(sim.in = sim.in, simulation = sim))

}
