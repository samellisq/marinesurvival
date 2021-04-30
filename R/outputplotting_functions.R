#' Link function
#'
#' @export
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
      if(input.list$include_samplebias_error > 0){
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


#' Model to Sample Comparison plot
#'
#'@export
modplot_samplescheck = function(age.seq, post, input.list, staninputdf, simdata, minage = 0){
  lambda = p_link(ageseq = age.seq, post = post, input.list = input.list)#*input.list$sum_dead
  lmu = unlist(lapply(lambda,  mean))#*input.list$sum_dead
  lci = dplyr::bind_rows(lapply(lambda,  PI))#*input.list$sum_dead
  gplot =
    ggplot2::ggplot(data = simdata, ggplot2::aes())+
    ggplot2::geom_line(data = data.frame(age = age.seq + minage, mu = lmu), ggplot2::aes(age, mu), colour = "red", size = 2)+
    ggplot2::geom_ribbon(data = data.frame(age = age.seq+ minage, mu = lmu, lCI = lci[,1], uCI = lci[,2]), ggplot2::aes(x = age, ymin = X5., ymax = X94.), fill = "red", alpha = 0.5)+
    ggplot2::geom_bar(data = staninputdf, ggplot2::aes(x = age.adj+ minage, y = n.dead), stat = "identity", alpha = 0.5)

  return(gplot)

}

#' Model predicted survival curves
#'
#'@export
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
  plot.curves = dplyr::bind_rows(post.curves[sample(seq(1:length(post.curves)), 150)])
  plot.curves$i = as.numeric(as.character(plot.curves$i))

  if(!is.null(sim.in)){
    inputted = function(x)exp(-1*(sim.in$a/sim.in$b)*(exp(sim.in$b*x)-1))
    plot.curves = dplyr::bind_rows(
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
    ggplot2::ggplot(plot.curves, ggplot2::aes(x = age, y = lx, group = i, colour = type))+
    ggplot2::geom_line(alpha = ifelse(plot.curves$type == "real", 1, 0.25),
              size = ifelse(plot.curves$type == "real", 1.5, 0.25))+
    scale_colour_manual(values = c("black", "red"))


  return(gplot)
}


#' Prediction Survial to Age plot
#'
#'@export
modplot_ageXaccuracy = function(X, post, age.seq, sim.in = NULL){

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
    lCI = stats::quantile(ageX, 0.05),
    uCI =stats::quantile(ageX, 0.95)
  )

  if(is.null(sim.in)){
    gplot =
      ggplot2::ggplot(plot.df, ggplot2::aes(x = x, y = mean))+
      ggplot2::geom_point()+
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lCI, ymax = uCI), width = 0)+
      ylim(c(min(ageX)-5,max(ageX)+5))+
      ylab("Age (mean +/- 95% CI)")+
      xlab("")
  } else {
    gplot =
      ggplot2::ggplot(plot.df, ggplot2::aes(x = x, y = mean))+
      ggplot2::geom_point()+
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lCI, ymax = uCI), width = 0)+
      ggplot2::geom_hline(yintercept = analytical_age_of_lx(a = sim.in$a, b = sim.in$b, X= X), linetype = "dashed", colour = "red")+
      ylim(c(min(ageX)-5,max(ageX)+5))+
      ylab("Age (mean +/- 95% CI)")+
      xlab("")
  }


  return(gplot)

}


#####################################################################
# LATENT MODEL FUNCTIONS
#############################################

#' Latent Link function
#'
#' @export
latent_p_link = function(ageseq, post, input.list, species.i){
  R = length(ageseq)
  lam = list()
  speceis.sample.size = length(input.list$sample_ages[input.list$sample_species == species.i])
  for(i in 1:R){
    S = max(dim(post[[1]]))
    lams = numeric(S)
    for(j in 1:S){
      a = post$a[j,species.i]
      b = post$b[j, species.i]
      r = post$r[j, species.i]
      if(input.list$include_popchange_error[species.i] > 0){
        rho = post$r[j]/(input.list$species_maxages[species.i]/2)
      } else{
        rho =0
      }
      if(input.list$include_samplebias_error[species.i] > 0){
        S = input.list$BiasMat[i, species.i]*post$s[j, species.i] +1
      } else {
        S = 1
      }

      lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
      sigma_lx = sum(exp(-1 * (a/b) * (exp(b * ageseq) - 1)))
      R = (1-rho)^ageseq[i]
      lams[j] = (lx / sigma_lx)*R*S*speceis.sample.size
    }
    lam[[i]] = lams
  }
  return(lam)
}

#' Latent Model to Sample Comparison plot
#'
#'@export
latentmodplot_samplescheck = function(age.seq, post, input.list, staninputdf, species = "all" , minages = NULL, names.key = NULL){

  if(species == "all"){
    species.is = seq(1,input.list$Nspecies,1)
  } else{
    species.is = species
  }

  if(is.null(minages)){
    minages = rep.int(0, length(species.is))
  }

  species.plots = list()
  for(k in 1:length(species.is)){

    lambda = latent_p_link(ageseq = age.seq, post = post, input.list = mod.list, species.i = species.is[k])
    lmu = unlist(lapply(lambda,  mean))
    lci = dplyr::bind_rows(lapply(lambda,  PI))
    minage = minages[k]

    sp.title = ifelse(is.null(names.key), k, names.key$species[k])
    plotdata = filter(staninputdf, species.num == k)
    plotdata$real.age = plotdata$age.adj + minages[k]
    gplot =
      ggplot(data = plotdata, aes())+
      geom_line(data = data.frame(age = age.seq + minage, mu = lmu), aes(age, mu), colour = "red", size = 2)+
      geom_ribbon(data = data.frame(age = age.seq+ minage, mu = lmu, lCI = lci[,1], uCI = lci[,2]), aes(x = age, ymin = X5., ymax = X94.), fill = "red", alpha = 0.5)+
      geom_bar(data = plotdata, aes(x = real.age, y = n.dead), stat = "identity", alpha = 0.5)+
      ggtitle(sp.title)+
      scale_x_continuous(limits = c(minages[k]-1, input.list$species_maxages[k]*1.25+minages[k] ))

    species.plots[[k]] = gplot
  }

  # colchooser = ifelse(length(species.plots) >3, 3, min(c(2,length(species.plots))))
  colchooser = floor(sqrt(length(species.plots)))
  print(grid.arrange(grobs = species.plots, ncol = colchooser))
  return(species.plots)
}

#' Latent Model predicted survival curves
#'
#'@export
latentmodplot_curvecomp = function(post, age.seq, input.list, names.key = NULL){
  allspec.plot.curves = list()

  for(j in 1:input.list$Nspecies){

    R = max(dim(post[[1]]))
    post.curves = list()
    for(i in 1:R){
      a = post$a[i,j]
      b = post$b[i,j]
      post.curves[[i]] = data.frame(age = age.seq,
                                    lx = exp(-1 * (a/b) * (exp(b * age.seq) - 1)),
                                    i,
                                    type = "predicted")
    }
    plot.curves = dplyr::bind_rows(post.curves[sample(seq(1:length(post.curves)), 150)])
    plot.curves$i = as.numeric(as.character(plot.curves$i))

    plot.curves$i = as.factor(plot.curves$i)
    plot.curves$species = j
    allspec.plot.curves[[j]] = plot.curves
  }
  plot.curves = bind_rows(allspec.plot.curves)

  colchooser = floor(sqrt(input.list$Nspecies))
  if(!is.null(names.key)){
    plot.curves = left_join(plot.curves, select(names.key, species = species.num, species.name = species))
    plot.curves$species = plot.curves$species.name
  }
  gplot =
    ggplot(plot.curves, aes(x = age, y = lx, group = i, colour = type))+
    geom_line(alpha =  0.25, size = 0.25)+
    scale_colour_manual(values = c("black"))+
    facet_wrap(.~species, ncol = colchooser)+
    theme(legend.position = "none")
  print(gplot)
  return(gplot)


}


#' Prediction Survial to Age plot latent model version
#'
#'@export
latentmodplot_ageXaccuracy = function(X, post, age.seq, species.is, minages = NULL, names.key = NULL, return.data = FALSE){

  if(is.null(minages)){
    minages = rep.int(0, length(species.is))
  }

  ageX.plotdata = list()
  for(j in 1:length(species.is)){

    R = max(dim(post[[1]]))
    post.curves = list()
    for(i in 1:R){
      a = post$a[i, j]
      b = post$b[i, j]
      post.curves[[i]] = data.frame(age = age.seq,
                                    lx = exp(-1 * (a/b) * (exp(b * age.seq) - 1)),
                                    i
      )
    }

    ageX= numeric(length(post.curves))
    for(i in 1:(length(post.curves))){
      x = post.curves[[i]]$age[which(abs(post.curves[[i]]$lx - X) == min(abs(post.curves[[i]]$lx - X)))]
      ageX[i] = x
    }

    ageX.plotdata[[j]] = data.frame(
      species = j,
      mean = mean(ageX) + minages[j],
      lCI = stats::quantile(ageX, 0.05) +minages[j],
      uCI =stats::quantile(ageX, 0.95)+ minages[j],
      sd = sd(ageX)
    )

  }

  ageX.plotdata = bind_rows(ageX.plotdata)
  if(!is.null(names.key)){
    ageX.plotdata = left_join(ageX.plotdata, select(names.key, species = species.num, species.name = species))
    ageX.plotdata$species.num = ageX.plotdata$species
    ageX.plotdata$species = ageX.plotdata$species.name
    ageX.plotdata = ageX.plotdata[,names(ageX.plotdata) != "species.name"]
  }
  ageX.plotdata$species = as.factor(ageX.plotdata$species)
  gplot =
    ggplot(ageX.plotdata, aes(x = species, y = mean))+
    geom_point()+
    geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0)+
    ylab(paste("Age at",X ," lx remaining (mean +/- 95% CI)", sep = ""))+
    xlab("Species")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


  if(return.data){
    print(gplot)
    return(ageX.plotdata)
  } else {
    print(gplot)
    return(gplot)
  }



}




######################################################
# Combined
########################################################
#' Latent Link function
#'
#' @export
comb_p_link = function(ageseq, post, input.list, species.i){
  if(input.list$Nspecies >1){
    R = length(ageseq)
    lam = list()
    speceis.sample.size = length(input.list$sample_ages[input.list$sample_species == species.i])
    for(i in 1:R){
      S = max(dim(post[[1]]))
      lams = numeric(S)
      for(j in 1:S){
        a = post$a[j,species.i]
        b = post$b[j, species.i]
        r = post$r[j, species.i]
        if(input.list$include_popchange_error[species.i] > 0){
          rho = post$r[j]/(input.list$species_maxages[species.i]/2)
        } else{
          rho =0
        }
        # if(input.list$include_samplebias_error[species.i] > 0){
          S = input.list$BiasMat[i, species.i]*post$s[j, species.i] +1
        # } else {
        #   S = 1
        # }

        lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
        sigma_lx = sum(exp(-1 * (a/b) * (exp(b * ageseq) - 1)))
        R = (1-rho)^ageseq[i]
        lams[j] = (lx / sigma_lx)*R*S*speceis.sample.size
      }
      lam[[i]] = lams
    }

  } else { # if only 1 species
    R = length(ageseq)
    lam = list()
    speceis.sample.size = length(input.list$sample_ages)
    for(i in 1:R){
      S = max(dim(post[[1]]))
      lams = numeric(S)
      for(j in 1:S){
        a = post$a[j]
        b = post$b[j]
        r = post$r[j]
        if(input.list$include_popchange_error[1] > 0){
          rho = post$r[j]/(input.list$species_maxages[1]/2)
        } else{
          rho =0
        }
        # if(input.list$include_samplebias_error[species.i] > 0){
          S = input.list$BiasMat[i, 1]*post$s[j] +1
        # } else {
        #   S = 1
        # }

        lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
        sigma_lx = sum(exp(-1 * (a/b) * (exp(b * ageseq) - 1)))
        R = (1-rho)^ageseq[i]
        lams[j] = (lx / sigma_lx)*R*S*speceis.sample.size
      }
      lam[[i]] = lams
    }

  }

  return(lam)

}


#' Latent Model to Sample Comparison plot
#'
#'@export
combmodplot_samplescheck = function(age.seq, post, input.list, staninputdf, species = "all" , minages = NULL, names.key = NULL){

  if(species == "all"){
    species.is = seq(1,input.list$Nspecies,1)
  } else{
    species.is = species
  }

  if(is.null(minages)){
    minages = rep.int(0, length(species.is))
  }

  species.plots = list()
  plot.titles = list()
  for(k in 1:length(species.is)){

    lambda = comb_p_link(ageseq = age.seq, post = post, input.list = mod.list, species.i = species.is[k])
    lmu = unlist(lapply(lambda,  mean))
    lci = dplyr::bind_rows(lapply(lambda,  PI))
    minage = minages[k]

    sp.title = ifelse(is.null(names.key), k, names.key$species[k])
    plotdata = filter(staninputdf, species.num == k)
    plotdata$real.age = plotdata$age.adj + minages[k]
    gplot =
      ggplot(data = plotdata, aes())+
      geom_line(data = data.frame(age = age.seq + minage, mu = lmu), aes(age, mu), colour = "red", size = 2)+
      geom_ribbon(data = data.frame(age = age.seq+ minage, mu = lmu, lCI = lci[,1], uCI = lci[,2]), aes(x = age, ymin = X5., ymax = X94.), fill = "red", alpha = 0.5)+
      geom_bar(data = plotdata, aes(x = real.age, y = n.dead), stat = "identity", alpha = 0.5)+
      ggtitle(sp.title)+
      scale_x_continuous(limits = c(minages[k]-1, input.list$species_maxages[k]*1.25+minages[k] ))

    species.plots[[k]] = gplot
  }

  # colchooser = ifelse(length(species.plots) >3, 3, min(c(2,length(species.plots))))
  colchooser = floor(sqrt(length(species.plots)))
  # print()
  # sp.grid = grid.arrange(grobs = species.plots, ncol = colchooser)
  sp.grid = cowplot::plot_grid(plotlist = species.plots, ncol = colchooser)
  return(sp.grid)
}

#' Latent Model predicted survival curves
#'
#'@export
combmodplot_curvecomp = function(post, age.seq, input.list, names.key = NULL, sim.in = NULL, N =150){
  allspec.plot.curves = list()

  if(input.list$Nspecies > 1) {
    for(j in 1:input.list$Nspecies){

      R = max(dim(post[[1]]))
      post.curves = list()
      for(i in 1:R){
        a = post$a[i,j]
        b = post$b[i,j]
        post.curves[[i]] = data.frame(age = age.seq,
                                      lx = exp(-1 * (a/b) * (exp(b * age.seq) - 1)),
                                      i,
                                      type = "predicted")
      }
      plot.curves = dplyr::bind_rows(post.curves[sample(seq(1:length(post.curves)), N)])
      plot.curves$i = as.numeric(as.character(plot.curves$i))

      plot.curves$i = as.factor(plot.curves$i)
      plot.curves$species = j
      allspec.plot.curves[[j]] = plot.curves
    }
    plot.curves = bind_rows(allspec.plot.curves)
  } else { # if n species >1
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
    plot.curves = dplyr::bind_rows(post.curves[sample(seq(1:length(post.curves)), N)])
    plot.curves$i = as.numeric(as.character(plot.curves$i))
    plot.curves$species = 1

  }


  colchooser = floor(sqrt(input.list$Nspecies))
  if(!is.null(names.key)){
    plot.curves = left_join(plot.curves, select(names.key, species = species.num, species.name = species))
    plot.curves$species = plot.curves$species.name
  }
  gplot =
    ggplot(plot.curves, aes(x = age, y = lx, group = i, colour = type))+
    geom_line(alpha =  0.25, size = 0.25)+
    scale_colour_manual(values = c("black"))+
    facet_wrap(.~species, ncol = colchooser)+
    theme(legend.position = "none")
  print(gplot)
  return(gplot)


}

#' Combined age X plot
#'
#' @export
combmodplot_ageXaccuracy = function(X, post, age.seq, species.is, minages = NULL, names.key = NULL, return.data = FALSE){

  if(is.null(minages)){
    minages = rep.int(0, length(species.is))
  }

  if(length(species.is) > 1){ # for multiple species
    ageX.plotdata = list()
    for(j in 1:length(species.is)){

      R = max(dim(post[[1]]))
      post.curves = list()
      for(i in 1:R){
        a = post$a[i, j]
        b = post$b[i, j]
        post.curves[[i]] = data.frame(age = age.seq,
                                      lx = exp(-1 * (a/b) * (exp(b * age.seq) - 1)),
                                      i
        )
      }
      ageX= numeric(length(post.curves))
      for(i in 1:(length(post.curves))){
        x = post.curves[[i]]$age[which(abs(post.curves[[i]]$lx - X) == min(abs(post.curves[[i]]$lx - X)))]
        ageX[i] = x
      }

      ageX.plotdata[[j]] = data.frame(
        species = j,
        mean = mean(ageX) + minages[j],
        lCI = stats::quantile(ageX, 0.05) +minages[j],
        uCI =stats::quantile(ageX, 0.95)+ minages[j],
        sd = sd(ageX)
      )

    }

    ageX.plotdata = bind_rows(ageX.plotdata)
  } else { # if there is only one species
    ageX.plotdata = list()
      R = length(post$a)
      post.curves = list()
      for(i in 1:R){
        a = post$a[i]
        b = post$b[i]
        post.curves[[i]] = data.frame(age = age.seq,
                                      lx = exp(-1 * (a/b) * (exp(b * age.seq) - 1)),
                                      i
        )
      }
      ageX= numeric(length(post.curves))
      for(i in 1:(length(post.curves))){
        x = post.curves[[i]]$age[which(abs(post.curves[[i]]$lx - X) == min(abs(post.curves[[i]]$lx - X)))]
        ageX[i] = x
      }

      ageX.plotdata = data.frame(
        species = 1,
        mean = mean(ageX) + minages,
        lCI = stats::quantile(ageX, 0.05) +minages,
        uCI =stats::quantile(ageX, 0.95)+ minages,
        sd = sd(ageX)
      )

    }


  if(!is.null(names.key)){
    ageX.plotdata = left_join(ageX.plotdata, select(names.key, species = species.num, species.name = species))
    ageX.plotdata$species.num = ageX.plotdata$species
    ageX.plotdata$species = ageX.plotdata$species.name
    ageX.plotdata = ageX.plotdata[,names(ageX.plotdata) != "species.name"]
  }
  ageX.plotdata$species = as.factor(ageX.plotdata$species)
  gplot =
    ggplot(ageX.plotdata, aes(x = species, y = mean))+
    geom_point()+
    geom_errorbar(aes(ymin = lCI, ymax = uCI), width = 0)+
    ylab(paste("Age at",X ," lx remaining (mean +/- 95% CI)", sep = ""))+
    xlab("Species")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


  if(return.data){
    print(gplot)
    return(ageX.plotdata)
  } else {
    print(gplot)
    return(gplot)
  }



}



