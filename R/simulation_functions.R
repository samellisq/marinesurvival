#' Simulate static survival sample
#'
#'@export
run_simulation = function(
  sim.in = NULL,
  include_r = TRUE,
  include_e = FALSE,
  include_s = FALSE,
  sample_size = 200,
  age.maturity = 0){

  if(is.null(sim.in)){
    #TODO: need a function to check sim.in input here
    sim.in =
      data.frame(
        a = stats::runif(1,0.01, 0.04), # baseline mortality taken from real examples
        b = stats::runif(1,0.04, 0.1), # ageing
        rho = stats::runif(1, -0.01, 0.01), # 0 = no change
        e = 0.05, # standard deviation in of age estimation error
        s =  stats::runif(1, 0.8, 1.2), # change in probability of being sampled above age S
        ageS = 25
      )
    #print(sim.in)
  }

  maxage = 75 #TODO make dynamics
  sampling.probs =
    data.frame(
      age = seq(age.maturity, maxage, 1),
      lx.pdf = gomp.pdf(x = seq(age.maturity, maxage, 1),
                        x.range = seq(age.maturity, maxage*1.5,1),
                        a = sim.in$a, b = sim.in$b
                        )
      )
  sampling.probs$obs.p = sampling.probs$lx.pdf


  if(include_r){
    sampling.probs$r.weight = (1-sim.in$rho)^seq(age.maturity, maxage, 1)
  }

  if(include_s){
    sampling.probs$s.weight = ifelse(sampling.probs$age >=sim.in$ageS,
                                     sim.in$s,
                                     1
    )
  }


  if(include_r == TRUE & include_s == FALSE){
    sampling.probs$obs.p = sampling.probs$lx.pdf*sampling.probs$r.weight
  }
  if(include_r == FALSE & include_s == TRUE){
    # print(sampling.probs)
    sampling.probs$obs.p = sampling.probs$lx.pdf*sampling.probs$s.weight
  }
  if(include_r == TRUE & include_s == TRUE){
    sampling.probs$obs.p = sampling.probs$lx.pdf*sampling.probs$r.weight*sampling.probs$s.weight
  }


  sampling.probs = dplyr::relocate(sampling.probs,obs.p,.after = last_col())
  gplot =
    ggplot2::ggplot(tidyr::pivot_longer(sampling.probs, cols = lx.pdf:obs.p, names_to = "p"),
                    ggplot2::aes(x = age, y = value, colour = p))+
    ggplot2::geom_line()+
    ggplot2::facet_wrap(.~p, ncol = 2, scales = "free_y")


  SAMPLE_SIZE = sample_size
  sample.cohort = sample(sampling.probs$age,
                         SAMPLE_SIZE,
                         prob = sampling.probs$obs.p,
                         replace = TRUE)
  age = sample.cohort
  sim = data.frame(age , surv = 1)

  gplot2 =  ggplot2::ggplot(sim, ggplot2::aes(age))+ggplot2::geom_histogram(binwidth = 1)

  sim.in$r = sim.in$rho*(max(sim)/2)

  if(include_e){
    sim$age.noerror = sim$age
    sim$e.sd = sim$age.noerror*sim.in$e
    for(i in 1:length(sim$age)){
      sim$age[i] = floor(stats::rnorm(1, sim$age[i]+0.5, sim$e.sd[i]))
    }
    # }

    sim = dplyr::filter(sim, age >=0)
  }

  #TODO: optional plotting
  gridExtra::grid.arrange(gplot, gplot2, ncol = 1)
  return(list(sim.in = sim.in, simulation = sim))

}
