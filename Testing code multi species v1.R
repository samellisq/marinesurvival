rm(list = ls())
require(marinesurvival)
require(tidyverse)
require(rethinking)

# The model only works on adults. Age at maturity is the age of adulthood. So the age of sexual maturity in practice


##################################
# Simulate some data
####################################
Nspecies = 2
AGE.AT.MATURITY = c(3,3)
Ndatasets.per.species = c(2,2)
simmed.datasets =list()
counter = 1
ageing.rates = c(0.1, 0.04)
for(i in 1:Nspecies){

  sim.in =
    #Every species has the same baseline rates
    data.frame(
      a = 0.02, # baseline mortality
      b = ageing.rates[i], # ageing
      rho = -0.05, # 0 = no change
      e = 0.05, # standard deviation in of age estimation error
      s =  2, # change in probability of being sampled above age S
      ageS = 35
    )

  for(j in 1:Ndatasets.per.species[i]){
    sim = run_simulation(sim.in = sim.in,
                         include_r = FALSE,
                         include_e = FALSE,
                         include_s = FALSE,
                         sample_size = c(40, 150)[i], # 100
                         age.maturity = AGE.AT.MATURITY[i]
    )
    simmed.datasets[[counter]] = sim$simulation$age
    counter = counter +1
  }

}
simmed.datasets
simmed.species = rep(seq(1, Nspecies),Ndatasets.per.species)
simmed.species
AGE.AT.MATURITY = rep(AGE.AT.MATURITY,Ndatasets.per.species)
AGE.AT.MATURITY


############################################
# Prepare the data for the model
#############################################

# Make a new data frame of n individuals of any given age
byage = list()
for(i in 1:length(simmed.datasets)){
  byage[[i]] =
    data.frame(
      age = seq(min(simmed.datasets[[i]]), 100+AGE.AT.MATURITY[i]),
      n.dead = 0,
      age.adj = 0,
      dataset.num = i,
      species.num = simmed.species[i],
      s.vector = 0
    )
  for(j in 1:length(simmed.datasets[[i]])){
    byage[[i]]$n.dead[byage[[i]]$age == simmed.datasets[[i]][j]] = byage[[i]]$n.dead[byage[[i]]$age == simmed.datasets[[i]][j]] +1
  }
  byage[[i]]$age.adj = byage[[i]]$age - AGE.AT.MATURITY[i]
  byage[[i]] = filter(byage[[i]], age.adj >=0)
}
byage



# Adjusted age is age reformated so 0 is AGE AT MATURITY.
sample.df = list()
for(i in 1:length(simmed.datasets)){
  sample.df[[i]] = data.frame(
    age = simmed.datasets[[i]],
    age.adj = simmed.datasets[[i]] - AGE.AT.MATURITY[i],
    dataset.num = i,
    species.num = simmed.species[i]
  )
  sample.df[[i]] = filter(sample.df[[i]], age.adj >=0)
}
sample.df

byage = bind_rows(byage)
sample.df = bind_rows(sample.df)

biasmat = list()
for(i in 1:max(byage$dataset.num)){
  biasmat[[i]] = filter(byage, dataset.num == i)$s.vector ##

}
biasmat = do.call(cbind, biasmat)
biasmat

maxages = numeric(max(byage$dataset.num))
for(i in 1:max(byage$dataset.num)){
  (maxages[i] = max(sample.df$age.adj[sample.df$dataset.num == i]))
}
maxages


popchange.direction = rep.int(0, length(simmed.datasets))

#################################
# Generate the model input data
######################################
# THe data has to be input as this horrible complex list.
# Some notes throughout but none are actually too complex

#Get some starting priors. These are based on some real data and don't matter too much but just give the model somewhere to start.
priors = get_gomp_prior_shapes()

mod.list = list(
  Nages = 101, # Number of ages in the data (e.g. 0-10, Nages = 11)
  Nsamples = nrow(sample.df), # Number of samples
  age = seq(0, max(byage$age.adj), 1), # All ages
  sample_ages = sample.df$age.adj, # Ages of sampled animals

  Nspecies = max(byage$species.num), # how many species in the data
  Ndatasets = max(byage$dataset.num),
  species_vector = as.array(simmed.species),
  sample_species = sample.df$dataset.num, # species id of each sample. For our purposes they are all the
  species_maxages = as.array(maxages), # Maximum possible age. Needed for some technical reasons.
  prior_a_s1 = priors$shape1_a, # Just the priors from before
  prior_a_s2 = priors$shape2_a,
  prior_b_s1 = priors$shape1_b,
  prior_b_s2 = priors$shape2_b,


  # all the following need to be 'as.array' to make it work. But it doesn;t actually change the data. just its 'type'. A necessarry quirk.
  include_age_est_error = as.array(rep.int(0, max(byage$dataset.num))), # 0 for do not include 1 for include
  age_error_sd = ifelse(sample.df$age == 0, 0.01, sample.df$age*0.05), # error based on real rather than adjusted age

  include_samplebias_error = as.array(rep.int(0, max(byage$dataset.num))), # 0 for do not include 1 for include
  BiasMat = biasmat, # A vector showing if a particualr age is subject toa biased sample (1) or not (0).
  direction_samplebias = as.array(rep.int(0, max(byage$dataset.num))), # +1/-1/0. 0 for UNK/Both

  include_popchange_error = as.array(rep.int(0, max(byage$dataset.num))), # 0 for do not include 1 for include
  direction_popchange = as.array(popchange.direction) # +1/-1/0. 0 for UNK/Both
)


########################################
# Run the model
##########################################

# This generates the code to run the model

stancode = get_marinesurvival_stancode(mod.list)

setwd("C:/Users/samel/Google Drive/Work/Comparative/Projects/Menopause Evolution/Analysis/Menopause Evolution II")

###
mod = stan(
  # file = "onespecies stan code v1.stan",
  model_code = stancode,
  data = mod.list,
  chains = 4,
  cores = 4,
  iter = 500,
  init = generate_inits(4, mod.list)
)


out = precis(mod, digits = 4, depth = 2)
activepars = get_activepars(out = out, mod = mod, input.list = mod.list)
out[rownames(out) %in% activepars,]
traceplot(mod, pars = activepars)
par(mfrow= c(1,1))

age.seq = seq(min(byage$age.adj), max(byage$age.adj),1)
post = extract.samples(mod)
# byage.data$age.adj = byage.data$age
minages = AGE.AT.MATURITY

names.key = data.frame(species.num = c(1,2,3), species = c("Aadvark","Bilby","Coati"))

plot_modtosample(age.seq = age.seq, post = post, input.list = mod.list, staninputdf = byage, minages = minages, names.key = names.key)
plot_posteriorsurvival(post = post, age.seq = age.seq, input.list = mod.list, N = 150, names.key = names.key)
plot_ageX(X = 0.1, post = post, age.seq = age.seq, species.is = c(1,2), minages = minages, return.data = FALSE, names.key = names.key)
