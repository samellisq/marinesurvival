require(tidyverse)
require(gridExtra)


###############################################################
# Simulations
###############################################################

run_simulation = function(sim.in = NULL, include_r = TRUE, include_e = FALSE, include_s = FALSE, sample_size = 200, age.maturity = 0){
  if(is.null(sim.in)){
    sim.in = 
      data.frame(
        a = runif(1,0.01, 0.04), # baseline mortality taken from real examples
        b = runif(1,0.04, 0.1), # ageing
        rho = runif(1, -0.01, 0.01), # 0 = no change
        e = 0.05, # standard deviation in of age estimation error
        s =  runif(1, 0.8, 1.2), # change in probability of being sampled above age S
        ageS = 25
      )
    print(sim.in)
  }
  
  maxage = 75
  sampling.probs = data.frame(age = seq(age.maturity, maxage, 1),
                              lx.pdf = gomp.pdf(x = seq(age.maturity, maxage, 1), x.range = seq(age.maturity, maxage*1.5,1), a = sim.in$a, b = sim.in$b)
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
  
  
  
  sampling.probs = relocate(sampling.probs,obs.p,.after = last_col()) 
  gplot = ggplot(pivot_longer(sampling.probs, cols = lx.pdf:obs.p, names_to = "p"),
                 aes(x = age, y = value, colour = p))+
    geom_line()+
    facet_wrap(.~p, ncol = 2, scales = "free_y")
  
  
  SAMPLE_SIZE = sample_size
  sample.cohort = sample(sampling.probs$age, 
                         SAMPLE_SIZE, 
                         prob = sampling.probs$obs.p,
                         replace = TRUE)
  age = sample.cohort
  sim = data.frame(age , surv = 1)
  
  gplot2 =  ggplot(sim, aes(age))+geom_histogram(binwidth = 1)
  
  sim.in$r = sim.in$rho*(max(sim)/2)
  
  if(include_e){
    sim$age.noerror = sim$age
    sim$e.sd = sim$age.noerror*sim.in$e
    for(i in 1:length(sim$age)){
      sim$age[i] = floor(rnorm(1, sim$age[i]+0.5, sim$e.sd[i]))
    }
    # }
    
    sim = filter(sim, age >=0)
  }
  
  grid.arrange(gplot, gplot2, ncol = 1)
  return(list(sim.in = sim.in, simulation = sim))
  
}

#######################################################
# ACCESORY FUNCTIONS
##########################################################

### Functions
gomp.pdf = function(x, x.range = seq(0:150) ,a, b){
  pd = numeric(length(x))
  for(i in 1:length(x)){
    lx = gomp_surv(x = x[i], a = a, b = b)
    sigma_lx = sum(gomp_surv(x = x.range, a = a, b = b))
    pd[i] = lx/sigma_lx
  }
  return(pd)
}


gomp_surv = function(x, a = NULL, b, M = NULL, by.M = FALSE){
  if(by.M){ # by modal age
    # https://www.monicaalexander.com/posts/2018-02-15-gompertz/
    Sx = exp(-exp(-b*M)*(exp(b*x)-1))
  } else { # by baseline mortality
    Sx = exp(-1*(a/b)*(exp(b*x)-1))
  }
  return(Sx)
} 

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

analytical_age_of_lx = function(a, b, X = 0.5){
  # see surv at age alegbra v1.jpg
  x = (log(log(X)*(b/-a)+1))/b
  return(x)
}


##############################################################
# GOMP priors from real data
##############################################################
## Gets gompterez priors
real.gomp = read_csv("gomp.surv.egs.csv")
shape1_a = mean(real.gomp$initital) / ((sd(real.gomp$initital)*2) * (sd(real.gomp$initital)*2)) # 2
shape2_a = (1-mean(real.gomp$initital)) / ((sd(real.gomp$initital)*2) * (sd(real.gomp$initital)*2)) #2 
shape1_b = mean(real.gomp$sense, na.rm = TRUE) / ((sd(real.gomp$sense, na.rm = TRUE)*1.5) * (sd(real.gomp$sense, na.rm = TRUE)*1.5)) #1.5
shape2_b = (1- mean(real.gomp$sense, na.rm = TRUE)) / ((sd(real.gomp$sense, na.rm = TRUE)*1.5) * (sd(real.gomp$sense, na.rm = TRUE)*1.5))# 1.5



###############################################################
# Exploring the data
################################################################
# p_link = function(ageseq, post, input.list){
#   # R = length(post[[1]])
#   R = length(ageseq)
#   lam = list()
#   for(i in 1:R){
#     lams = numeric(length(post[[1]]))
#     for(j in 1:length(post[[1]])){
#       a = post$a[j]
#       b = post$b[j]
#       if(input.list$include_popchange_error > 0 ){
#         r = post$r[j]
#         rho = post$r[j]/(max(input.list$sample_ages)/2)
#         R = (1-rho)^ageseq[i]
#       } else {
#         R = 1 
#       }
#       if(input.list$include_sampling_error > 0){
#         S = input.list$s_vector[i]*post$s[j] +1 
#       } else {
#         S = 1
#       }
#       
#       lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
#       sigma_lx = sum(exp(-1 * (a/b) * (exp(b * ageseq) - 1)))
#       
#       num = lx*R*S
#       den = sigma_lx*R*S
#       #lams[j] = (lx / sigma_lx)*input.list$sum_dead*R*S
#       lams[j] = (num/den)*input.list$sum_dead # this just turns the theta into a prediction about how many will be drawn. This is not elegant there is a btter way
#     }
#     lam[[i]] = lams
#   }
#   return(lam)
# }

p_link = function(ageseq, post, input.list){
  # R = length(post[[1]])
  R = length(ageseq)
  lam = list()
  for(i in 1:R){
    lams = numeric(length(post[[1]]))
    for(j in 1:length(post[[1]])){
      a = post$a[j]
      b = post$b[j]
      r = post$r[j]
      rho = post$r[j]/(max(input.list$sample_ages)/2)
      if(input.list$include_sampling_error > 0){
        S = input.list$s_vector[i]*post$s[j] +1 
      } else {
        S = 1
      }
      
      lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
      sigma_lx = sum(exp(-1 * (a/b) * (exp(b * ageseq) - 1)))
      R = (1-rho)^ageseq[i]
      lams[j] = (lx / sigma_lx)*input.list$sum_dead*R*S
    }
    lam[[i]] = lams
  }
  return(lam)
}

modplot_samplescheck = function(age.seq, post, input.list, staninputdf, simdata, minage = 0){
  lambda = p_link(ageseq = age.seq, post = post, input.list = input.list)#*input.list$sum_dead 
  lmu = unlist(lapply(lambda,  mean))#*input.list$sum_dead
  lci = bind_rows(lapply(lambda,  PI))#*input.list$sum_dead
  gplot = 
  ggplot(data = simdata, aes())+
    geom_line(data = data.frame(age = age.seq + minage, mu = lmu), aes(age, mu), colour = "red", size = 2)+
    geom_ribbon(data = data.frame(age = age.seq+ minage, mu = lmu, lCI = lci[,1], uCI = lci[,2]), aes(x = age, ymin = X5., ymax = X94.), fill = "red", alpha = 0.5)+
    geom_bar(data = staninputdf, aes(x = age.adj+ minage, y = n.dead), stat = "identity", alpha = 0.5)
  
  return(gplot)
  
}


modplot_curvecomp = function(sim.in = NULL, post, age.seq){
    
    R = length(post$a)
    post.curves = list()
    for(i in 1:R){
      a = post$a[i]
      b = post$b[i]
      post.curves[[i]] = data.frame(age = age.seq, 
                                    lx = exp(-1 * (a/b) * (exp(b * age.seq) - 1)),
                                    i,
                                    type = "predicted")
    }
    plot.curves = bind_rows(post.curves[sample(seq(1:length(post.curves)), 150)])
    post.curves$i = as.numeric(as.character(post.curves$i))
    
    if(!is.null(sim.in)){
      inputted = function(x)exp(-1*(sim.in$a/sim.in$b)*(exp(sim.in$b*x)-1))
      plot.curves = bind_rows(
        plot.curves,
        data.frame(
          age = age.seq,
          lx = inputted(age.seq),
          i = max(plot.curves$i)+1,
          type = "real"
        )
      )
      
    } 
    
    plot.curves$i = as.factor(plot.curves$i)
    
    gplot = 
    ggplot(plot.curves, aes(x = age, y = lx, group = i, colour = type))+
      geom_line(alpha = ifelse(plot.curves$type == "real", 1, 0.25), 
                size = ifelse(plot.curves$type == "real", 1.5, 0.25))+
      scale_colour_manual(values = c("black", "red"))
    
    
  return(gplot)
}


modplot_ageXaccuracy = function(X, post, age.seq, sim.in){
  
  R = length(post$a)
  post.curves = list()
  for(i in 1:R){
    a = post$a[i]
    b = post$b[i]
    post.curves[[i]] = data.frame(age = age.seq, 
                                  lx = exp(-1 * (a/b) * (exp(b * age.seq) - 1)),
                                  i,
                                  type = "predicted")
  }

  #plot.curves = bind_rows(post.curves[sample(seq(1:length(post.curves)), 150)])
  #post.curves$i = as.numeric(as.character(post.curves$i))
 
  ageX= numeric(length(post.curves))
  for(i in 1:(length(post.curves))){
    x = post.curves[[i]]$age[which(abs(post.curves[[i]]$lx - X) == min(abs(post.curves[[i]]$lx - X)))]
    ageX[i] = x
  }
  
  
  plot.df = data.frame(
    x = " Modelled Age X(50% ex)",
    mean = mean(ageX), 
    lCI = quantile(ageX, 0.05), 
    uCI =quantile(ageX, 0.95)
  )
  
  gplot =
  ggplot(plot.df, aes(x = x, y = mean))+
    geom_point()+
    geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0)+
    geom_hline(yintercept = analytical_age_of_lx(a = sim.in$a, b = sim.in$b, X= X), linetype = "dashed", colour = "red")+
    ylim(c(min(ageX)-5,max(ageX)+5))+
    ylab("Age (mean +/- 95% CI)")+
    xlab("")
  
  return(gplot)
  
}





