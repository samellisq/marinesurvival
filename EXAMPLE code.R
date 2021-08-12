require(marinesurvival)
require(tidyverse)
require(rethinking)

# The model only works on adults. Age at maturity is the age of adulthood. So the age of sexual maturity in practice
AGE.AT.MATURITY = 3

##################################
# Simulate some data
####################################
sim1.in =
  data.frame(
    a = 0.01, # baseline mortality
    b = 0.07, # ageing
    rho = -0.05, # 0 = no change
    e = 0.05, # standard deviation in of age estimation error
    s =  2, # change in probability of being sampled above age S
    ageS = 35
  )

sim1.list = run_simulation(sim.in = sim1.in,
                           include_r = FALSE,
                           include_e = TRUE,
                           include_s = FALSE,
                           sample_size = 50,
                           age.maturity = AGE.AT.MATURITY # if this isn't zero something in the mod1el intepretation
)
sim1.in = sim1.list$sim.in
sample = sim1.list$simulation$age # This is what any real data should look like
sample

############################################
# Prepare the data for the model
#############################################

# Make a new data frame of n individuals of any given age
byage =
  data.frame(
    age = seq(min(sample), max(sample)*1.5),
    n.dead = 0
    )
for(i in 1:length(sample)){
  byage$n.dead[byage$age == sample[i]] = byage$n.dead[byage$age == sample[i]] +1
}
byage

# Adjusted age is age reformated so 0 is AGE AT MATURITY.
byage$age.adj = byage$age - AGE.AT.MATURITY

# Sample data made into a dataframe with the adjusted ages calculated
sample.df = data.frame(
  age = sample,
  age.adj = sample - AGE.AT.MATURITY
)

# Filter the dataframes so that they only include adults
sample.df = filter(sample.df, age.adj >=0)
byage = filter(byage, age.adj >=0)


##Sampling bias
byage$s.bias = c(rep.int(0, 35), rep.int(1, nrow(byage)-35))
byage$s.bias

# An annoying quirk but the need a column called 'species.num'. Comes in useful later but can just be 1.
sample.df$species.num = 1
byage$species.num = 1



#################################
# Generate the model input data
######################################
# THe data has to be input as this horrible complex list.
# Some notes throughout but none are actually too complex

#Get some starting priors. These are based on some real data and don't matter too much but just give the model somewhere to start.
priors = get_gomp_prior_shapes()

mod.list = list(
  Nages = nrow(byage), # Number of ages in the data (e.g. 0-10, Nages = 11)
  Nsamples = nrow(sample.df), # Number of samples
  age = byage$age.adj, # All ages
  sample_ages = sample.df$age.adj, # Ages of sampled animals

  Nspecies = 1, # how many species in the data
  sample_species = sample.df$species.num, # species id of each sample. For our purposes they are all the
  species_maxages = as.array(max(byage$age.adj)), # Maximum possible age. Needed for some technical reasons.
  prior_a_s1 = priors$shape1_a, # Just the priors from before
  prior_a_s2 = priors$shape2_a,
  prior_b_s1 = priors$shape1_b,
  prior_b_s2 = priors$shape2_b,


  # all the following need to be 'as.array' to make it work. But it doesn;t actually change the data. just its 'type'. A necessarry quirk.
  include_age_est_error = as.array(1), # 0 for do not include 1 for include
  age_error_sd = ifelse(sample.df$age == 0, 0.01, sample.df$age*0.05), # error based on real rather than adjusted age

  include_samplebias_error = as.array(0), # 0 for do not include 1 for include
  BiasMat = cbind(byage$s.bias), # A vector showing if a particualr age is subject toa biased sample (1) or not (0).
  direction_samplebias = as.array(1), # +1/-1/0. 0 for UNK/Both

  include_popchange_error = as.array(0), # 0 for do not include 1 for include
  direction_popchange = as.array(0) # +1/-1/0. 0 for UNK/Both
)


########################################
# Run the model
##########################################

# This generates the code to run the model
stancode = get_marinesurvival_stancode(mod.list)

###
mod = stan(
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

plot_modtosample(age.seq = age.seq, post = post, input.list = mod.list, staninputdf = byage, minages = minages)
plot_posteriorsurvival(post = post, age.seq = age.seq, input.list = mod.list, N = 150, sim.in = sim1.in)
plot_ageX(X = 0.1, post = post, age.seq = age.seq, species.is = 1, minages = minages, return.data = TRUE)
