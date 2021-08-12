#' Link Posterior estimates into output
#'
#' Works for both 1 and many species
#'
#' @export
p_link = function(ageseq, post, input.list, species.i, dataset.size = NULL){
  if(input.list$Nspecies >1){
    R = length(ageseq)
    lam = list()
    if(is.null(dataset.size)){
      #so if there is no sample size provided calcualte it for all examples of the species (which could be one)
      speceis.sample.size = length(input.list$sample_ages[input.list$sample_species == species.i])
    } else {
      speceis.sample.size = dataset.size
    }

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


#' Model to Sample Comparison plot
#'
#'@export
plot_modtosample = function(age.seq, post, input.list, staninputdf, datasets = "all" , minages = NULL, names.key = NULL, return.data = FALSE){

  if(datasets == "all"){
    datasets.is = seq(1,input.list$Ndatasets,1)
  } else{
    datasets.is = datasets
  }

  species.is = input.list$species_vector[datasets.is]

  if(is.null(minages)){
    minages = rep.int(0, length(datasets.is))
  }

  if(!("n.dead" %in% names(staninputdf))){
    staninputdf$n.dead = staninputdf$n
  }

  datasets.plots = list()
  datasets.plots.data = list()




  for(k in 1:length(datasets.is)){

    lambda = p_link(ageseq = age.seq,
                    post = post,
                    input.list = mod.list,
                    species.i = species.is[k],
                    dataset.size = sum(byage$n.dead[byage$dataset.num == k]))
    lmu = unlist(lapply(lambda,  mean))
    lci = dplyr::bind_rows(lapply(lambda,  PI))
    minage = minages[k]

    if(!is.null(names.key)){
      sp.title = names.key$species[species.is][k]
    } else {
      sp.title = paste("dataset", k, sep = " ")
    }

    #sp.title = ifelse(is.null(datasets), k, names.key$species[k])
    plotdata = dplyr::filter(staninputdf, dataset.num == k)
    plotdata$real.age = plotdata$age.adj + minages[k]
    datasets.plots.data[[k]] = plotdata

    if(!return.data){
      gplot =
        ggplot2::ggplot(data = plotdata, ggplot2::aes())+
        ggplot2::geom_line(data = data.frame(age = age.seq + minage, mu = lmu), ggplot2::aes(age, mu), colour = "red", size = 2)+
        ggplot2::geom_ribbon(data = data.frame(age = age.seq+ minage, mu = lmu, lCI = lci[,1], uCI = lci[,2]), ggplot2::aes(x = age, ymin = X5., ymax = X94.), fill = "red", alpha = 0.5)+
        ggplot2::geom_bar(data = plotdata, ggplot2::aes(x = real.age, y = n.dead), stat = "identity", alpha = 0.5)+
        ggplot2::ggtitle(sp.title)+
        ggplot2::scale_x_continuous(limits = c(minages[k]-1, input.list$species_maxages[k]*1.25+minages[k] ))

      datasets.plots[[k]] = gplot
    }

  }

  if(!return.data){
    # colchooser = ifelse(length(species.plots) >3, 3, min(c(2,length(datasets.plots))))
    colchooser = floor(sqrt(length(datasets.plots)))
    # print()
    # sp.grid = grid.arrange(grobs = species.plots, ncol = colchooser)
    sp.grid = cowplot::plot_grid(plotlist = datasets.plots, ncol = colchooser)
    return(sp.grid)
  } else{
    return(datasets.plots.data)
  }


}

#' Model predicted survival curves
#'
#'@export
plot_posteriorsurvival = function(post, age.seq, input.list, names.key = NULL, sim.in = NULL, N =150){
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
    plot.curves = dplyr::bind_rows(allspec.plot.curves)
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
    plot.curves = dplyr::left_join(plot.curves, dplyr::select(names.key, species = species.num, species.name = species))
    plot.curves$species = plot.curves$species.name
  }
  gplot =
    ggplot2::ggplot(plot.curves, ggplot2::aes(x = age, y = lx, group = i, colour = type))+
    ggplot2::geom_line(alpha =  0.25, size = 0.25)+
    ggplot2::scale_colour_manual(values = c("black"))+
    ggplot2::facet_wrap(.~species, ncol = colchooser)+
    ggplot2::theme(legend.position = "none")
  print(gplot)
  return(gplot)


}

#' Age X plot
#'
#' @export
plot_ageX = function(X, post, age.seq, species.is = NULL, minages = NULL, names.key = NULL, return.data = FALSE){

  if(is.null(species.is)){
    species.is = 1
  }

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
        sd = stats::sd(ageX)
      )

    }

    ageX.plotdata = dplyr::bind_rows(ageX.plotdata)
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
        sd = stats::sd(ageX)
      )

    }


  if(!is.null(names.key)){
    ageX.plotdata = dplyr::left_join(ageX.plotdata, dplyr::select(names.key, species = species.num, species.name = species))
    ageX.plotdata$species.num = ageX.plotdata$species
    ageX.plotdata$species = ageX.plotdata$species.name
    ageX.plotdata = ageX.plotdata[,names(ageX.plotdata) != "species.name"]
  }
  ageX.plotdata$species = as.factor(ageX.plotdata$species)
  gplot =
    ggplot2::ggplot(ageX.plotdata, ggplot2::aes(x = species, y = mean))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lCI, ymax = uCI), width = 0)+
    ggplot2::ylab(paste("Age at",X ," lx remaining (mean +/- 95% CI)", sep = ""))+
    ggplot2::xlab("Species")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1))


  if(return.data){
    print(gplot)
    return(ageX.plotdata)
  } else {
    print(gplot)
    return(gplot)
  }



}



