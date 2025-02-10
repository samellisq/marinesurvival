#'
#'@export
p_link = function(){
  return(1)
}

#' Link Posterior estimates for a dataset into output
#'
#' Works for both 1 and many species
#'
#' @export
p_link_dataset = function(ageseq, post, input.list, dataset.i, dataset.size = NULL){
  if(input.list$Nspecies >1){

    dataset.sample.size = length(input.list$sample_ages[input.list$dataset_vector == dataset.i])
    dataset.species.i = input.list$species_vector[dataset.i]
    dataset.pop.i = input.list$population_vector[dataset.i]

    R = length(ageseq)
    lam = list()
    for(i in 1:R){
      S = max(dim(post[[1]]))
      lams = numeric(S)
      for(j in 1:S){
        a = post$a[j,dataset.species.i]
        b = post$b[j, dataset.species.i]
        r = post$r[j, dataset.pop.i]
        if(input.list$include_popchange_error[dataset.i] > 0){
          rho = post$r[j]/(input.list$species_maxages[dataset.i]/2)
        } else{
          rho =0
        }
        if(input.list$include_samplebias_error[dataset.i] > 0){
          S = input.list$BiasMat[i, dataset.i]*post$s[j, dataset.i] +1
        } else {
          S = 1
        }

        lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
        sigma_lx = sum(exp(-1 * (a/b) * (exp(b * ageseq) - 1)))
        R = (1-rho)^ageseq[i]
        lams[j] = (lx / sigma_lx)*R*S*dataset.sample.size
      }
      lam[[i]] = lams
    }

  } else { # if only 1 species
    dataset.sample.size = length(input.list$sample_ages[input.list$dataset_vector == dataset.i])
    dataset.species.i = input.list$species_vector[dataset.i]
    dataset.pop.i = input.list$population_vector[dataset.i]

    R = length(ageseq)
    lam = list()
    for(i in 1:R){
      S = max(dim(post[[1]]))
      lams = numeric(S)
      for(j in 1:S){
        a = post$a[j]
        b = post$b[j]
        r = post$r[dataset.pop.i]
        if(input.list$include_popchange_error[1] > 0){
          rho = post$r[dataset.pop.i]/(input.list$species_maxages[1]/2)
        } else{
          rho =0
        }
        if(input.list$include_samplebias_error[dataset.i] > 0){
          S = input.list$BiasMat[i, 1]*post$s[j] +1
        } else {
          S = 1
        }

        lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
        sigma_lx = sum(exp(-1 * (a/b) * (exp(b * ageseq) - 1)))
        R = (1-rho)^ageseq[i]
        lams[j] = (lx / sigma_lx)*R*S*dataset.sample.size
      }
      lam[[i]] = lams
    }

  }

  return(lam)

}

#' Link Posterior estimates for a species into output
#'
#' Works for both 1 and many species
#'
#' @export
p_link_species = function(ageseq, post, input.list, species.i, dataset.size = NULL){
  if(input.list$Nspecies >1){

    R = length(ageseq)
    lam = list()
    for(i in 1:R){
      S = max(dim(post[[1]]))
      lams = numeric(S)
      for(j in 1:S){
        a = post$a[j,species.i]
        b = post$b[j, species.i]
        lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
        lams[j] = (lx )
      }
      lam[[i]] = lams
    }

  } else { # if only 1 species

    R = length(ageseq)
    lam = list()
    for(i in 1:R){
      S = max(dim(post[[1]]))
      lams = numeric(S)
      for(j in 1:S){
        a = post$a[j]
        b = post$b[j]
        lx = exp(-1 * (a/b) * (exp(b * ageseq[i]) - 1))
        lams[j] = (lx)
      }
      lam[[i]] = lams
    }

  }

  return(lam)

}



#' Model to Sample Comparison plot
#'
#'@export
plot_modtosample = function(age.seq, post, input.list, dataset.is = NULL , minages = NULL, names.key = NULL, return.data = FALSE, datasetskey = NULL){

  if(is.null(dataset.is)){
    dataset.is = seq(1,input.list$Ndatasets,1)
  } else{
    dataset.is = dataset.is
  }


  species.is = input.list$species_vector[dataset.is]

  if(is.null(minages)){
    minages = rep.int(0, length(dataset.is))
  }


  datasets.plots = list()
  datasets.plots.data = list()


  for(k in 1:length(dataset.is)){

    # dataset.n = dataset.is[k]

    lambda = p_link_dataset(ageseq = age.seq,
                    post = post,
                    input.list = input.list,
                    dataset.i = dataset.is[k])
    lmu = unlist(lapply(lambda,  mean))
    lci = dplyr::bind_rows(lapply(lambda,  rethinking::PI))
    minage = minages[dataset.is[k]]

    if(!is.null(names.key)){

      if(!is.null(datasetskey)){
        dataset.n = filter(datasetskey, dataset.num == dataset.is[k])$dataset
        datasetinfo = paste(": dataset #", dataset.n,sep = "")
      } else {
        dataset.n = dataset.is[k]
        datasetinfo = paste(": dataset.num #", dataset.n,sep = "")
      }

      sp.deets = filter(names.key, species.num == species.is[k])
      sp.title = paste(sp.deets$species, "_", sp.deets$sex, ": dataset.num #", dataset.n, sep = "")
    } else {
      sp.title = paste("dataset", k, sep = " ")
    }

    #sp.title = ifelse(is.null(datasets), k, names.key$species[k])
    sample = input.list$sample_ages[input.list$dataset_vector == dataset.is[k]]
    plotdata = data.frame(age.adj = age.seq, n.dead = 0)
    for(a in 1:nrow(plotdata)){
      plotdata$n.dead[a] = sum(sample == plotdata$age.adj[a])
    }
    # plotdata = data.frame(age.adj = input.list$sample_ages[input.list$dataset_vector == k])
    # plotdata = dplyr::filter(staninputdf, dataset.num == dataset.is[k])

    plotdata$real.age = plotdata$age.adj + minages[species.is[k]]
    datasets.plots.data[[k]] = plotdata

    modplotdata = data.frame(
      age = age.seq +  minages[species.is[k]],
      mu = lmu,
      lCI = lci[,1],
      uCI = lci[,2]
    )
    names(modplotdata)[3:4] = c("lCI", "uCI")

    modplotdata$lCI = ifelse(modplotdata$lCI<0, 0, modplotdata$lCI)
    modplotdata$uCI = ifelse(modplotdata$uCI> max(plotdata$n.dead)*5,max(plotdata$n.dead)*5, modplotdata$uCI )
    modplotdata$mu = ifelse(modplotdata$mu> max(plotdata$n.dead)*5,max(plotdata$n.dead)*5, modplotdata$mu )
    modplotdata$lCI = ifelse(modplotdata$lCI> max(plotdata$n.dead)*5,max(plotdata$n.dead)*5, modplotdata$lCI )

    if(!return.data){
      print(k)
      gplot =
        ggplot2::ggplot(data = plotdata, ggplot2::aes())+
        ggplot2::geom_line(data = modplotdata, ggplot2::aes(age, mu), colour = "red", size = 2)+
        ggplot2::geom_ribbon(data = modplotdata, ggplot2::aes(x = age, ymin = lCI, ymax = uCI), fill = "red", alpha = 0.5)+
        ggplot2::geom_bar(data = plotdata, ggplot2::aes(x = real.age, y = n.dead), stat = "identity", alpha = 0.5)+
        ggplot2::ggtitle(sp.title)#+
      # ggplot2::scale_x_continuous(limits = c(minages[dataset.is[k]]-1, input.list$species_maxages[dataset.is[k]]*1.25+minages[dataset.is[k]] ))

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
plot_posteriorsurvival = function(post, age.seq, input.list, names.key = NULL, N =150, species.is = NULL){

  if(is.null(species.is)){
    species.is = seq(1:input.list$Nspecies)
  }


  allspec.plot.curves = list()



  for(i in 1:length(species.is)){
    p_lam = p_link_species(ageseq = age.seq, post = post, input.list = input.list, species.i = i)
    n.samples = length(p_lam[[1]])
    iter = rep(seq(1:n.samples), times = length(age.seq))
    age = rep(age.seq, each = n.samples)
    curves.df = lapply(p_lam, as.data.frame)
    curves.df = bind_rows(curves.df)
    curves.df = data.frame(iter, age, lx = curves.df[,1])
    samples = sample(seq(1:n.samples), N)
    curves.df = filter(curves.df, iter %in% samples)
    curves.df$species = species.is[i]
    curves.df$curve = seq(1:N)
    allspec.plot.curves[[i]] = curves.df
  }
  plot.curves = dplyr::bind_rows(allspec.plot.curves)


  colchooser = floor(sqrt(input.list$Nspecies))
  if(!is.null(names.key)){
    plot.curves = dplyr::left_join(plot.curves, dplyr::select(names.key, species = species.num, species.name = species, sex = sex))
    plot.curves$species = paste(plot.curves$species.name, plot.curves$sex, sep = "_")
  }
  gplot =
    ggplot2::ggplot(plot.curves, ggplot2::aes(x = age, y = lx, group = as.factor(curve)))+
    ggplot2::geom_line(alpha =  0.25, size = 0.25)+
    # ggplot2::scale_colour_manual(values = c("black"))+
    ggplot2::facet_wrap(.~species, ncol = colchooser)+
    ggplot2::theme(legend.position = "none")
  print(gplot)
  return(gplot)


}

#' Age X plot
#'
#' @export
plot_ageX = function(X, post, input.list, age.seq, species.is = NULL, minages = NULL, THIN = NULL, names.key = NULL, return.data = FALSE, return.summary = TRUE){

  ##THIS DOESNT QUITE WORK
  if(is.null(species.is)){
    species.is = seq(1:input.list$Nspecies)

  }

  if(is.null(minages)){
    minages = rep.int(0, length(species.is))
  }

  # if(length(species.is) > 1){ # for multiple species
  ageX.plotdata = list()
  ageX.all = list()

  for(i in 1:length(species.is)){
    p_lam = p_link_species(ageseq = age.seq, post = post, input.list = input.list, species.i = i)
    if(!is.null(THIN)){
      sample = sample(seq(1, length(p_lam[[1]]), 1),THIN)
      for(t in 1:length(p_lam))
        p_lam[[t]] = p_lam[[t]][sample]
    }
    n.samples = length(p_lam[[1]])
    iter = rep(seq(1:n.samples), times = length(age.seq))
    age = rep(age.seq, each = n.samples)
    curves.df = lapply(p_lam, as.data.frame)
    curves.df = bind_rows(curves.df)
    curves.df = data.frame(iter, age, lx = curves.df[,1])
    curves.df$species = species.is[i]

    S = max(curves.df$iter)
    ageX= numeric(S) # this is too many. its all the iters * all the ages
    for(j in 1:S){
      curve.df = dplyr::filter(curves.df, iter == j)
      x = curve.df$age[which(abs(curve.df$lx - X) == min(abs(curve.df$lx - X)))] # in here need to indicate wchi iter
      ageX[j] = x
    }

    ageX.plotdata[[i]] = data.frame(
      species = species.is[i],
      mean = mean(ageX) + minages[i],
      lCI = stats::quantile(ageX, 0.05) +minages[i],
      uCI =stats::quantile(ageX, 0.95)+ minages[i],
      sd = stats::sd(ageX)
    )
    ageX.all[[i]] = ageX

  }


  ageX.plotdata = dplyr::bind_rows(ageX.plotdata)

  if(!is.null(names.key)){
    ageX.plotdata = dplyr::left_join(ageX.plotdata, dplyr::select(names.key, species = species.num, species.name = species, sex = sex) %>% distinct())
    ageX.plotdata$species.num = ageX.plotdata$species
    ageX.plotdata$species = ageX.plotdata$species.name
    ageX.plotdata = ageX.plotdata[,names(ageX.plotdata) != "species.name"]
  } else{
    ageX.plotdata$sex = "U"
  }
  ageX.plotdata$species = as.factor(ageX.plotdata$species)
  gplot =
    ggplot2::ggplot(ageX.plotdata, ggplot2::aes(x = species, y = mean))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lCI, ymax = uCI), width = 0)+
    ggplot2::ylab(paste("Age at",X ," lx remaining (mean +/- 95% CI)", sep = ""))+
    ggplot2::xlab("Species")+
    ggplot2::facet_grid(rows = vars(sex))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust=1))


  if(return.data){
    if(return.summary == TRUE){
      return(ageX.plotdata)
    } else {
      sexnames = ifelse(input.list$species_sex_vector == 1, "F", "M")
      names(ageX.all) = sexnames

      return(ageX.all)
    }

  } else {
    print(gplot)
    return(gplot)
  }



}


#' Wrapper for plots
#'
#' @export
plot_wrapper = function(species, stanmodel, plot.type, thin = 100){

  if(!plot.type %in% c("a", "b", "c", "d", "mod.to.sample", "post.survival", "ageX.plot", "ageX.data")){
    warning("Error: plot.type must be either: a, b, c, d, mod.to.sample, post.survival, ageX.plot or ageX.data")
    stop()
  }

  input.dat = marinelifehistdata::get_lifehist_data(data.type = "age-structure", species = species, sex = NULL)
  mod.list = marinelifehistdata::create_marinesurvival_modinput(input.dat)

  age.seq = seq(0, 100,1)
  post = rstan::extract(stanmodel)
  SP = species
  ages.at.mat =
    marinelifehistdata::marine.lifehist.speciesdata$species_age.maturity %>%
    filter(species == SP) %>%
    arrange(sex)
  ages.at.mat = ages.at.mat$age.mat

  datasets.key = get_lifehist_data(data.type = "age-structure", species = species, sex = NULL, return.key = TRUE)[[1]]
  datasets.key$species.num = mod.list$species_vector
  datasets.key$dataset.num = seq(1,nrow(datasets.key), 1)

  if(plot.type == "a" | plot.type == "mod.to.sample"){
    plot = marinesurvival::plot_modtosample(age.seq = age.seq, post = post, input.list = mod.list, minages = ages.at.mat, names.key = datasets.key, datasetskey = datasets.key)
  }

  if(plot.type == "b" | plot.type == "post.survival"){
    plot = marinesurvival::plot_posteriorsurvival(post = post, age.seq = age.seq, input.list = mod.list, names.key = datasets.key, N = thin)
  }

  if(plot.type == "c" | plot.type == "ageX.plot"){
    plot = marinesurvival::plot_ageX(X = 0.1, post = post, age.seq = age.seq, input.list = mod.list, minages = ages.at.mat, names.key = datasets.key, THIN = thin)
    # plot = plot_ageX2(X = 0.1, post = post, age.seq = age.seq, input.list = mod.list, minages = ages.at.mat, names.key = datasets.key, THIN = thin)
  }

  if(plot.type == "d" | plot.type == "ageX.data"){
    plot = marinesurvival::plot_ageX(X = 0.1, post = post, age.seq = age.seq, input.list = mod.list, minages = ages.at.mat, names.key = datasets.key, THIN = thin, return.data = TRUE)
    # plot = plot_ageX2(X = 0.1, post = post, age.seq = age.seq, input.list = mod.list, minages = ages.at.mat, names.key = datasets.key, THIN = thin, return.data = TRUE)
  }

  return(plot)
}

