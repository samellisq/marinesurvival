# rm(list = ls())
#
# library(marinesurvival)
# require(rethinking)
# require(gridExtra)
#
# sim1.in =
#   data.frame(
#     a = 0.01, # baseline mortality
#     b = 0.07, # ageing
#     rho = -0.007, # 0 = no change
#     e = 0.05, # standard deviation in of age estimation error
#     s =  1, # change in probability of being sampled above age S
#     ageS = 35
#   )
#
# sim1.list = run_simulation(sim.in = sim1.in,
#                            include_r = FALSE,
#                            include_e = FALSE,
#                            include_s = FALSE,
#                            sample_size = 250,
#                            age.maturity = 0 # if this isn't zero something in the mod1el intepretation
# )
# sim1.in = sim1.list$sim.in
# sim1.data = sim1.list$simulation
#
#
# mod1.df =
#   data.frame(age = seq(min(sim1.data$age), max(sim1.data$age)*1.5,1))
# mod1.df$n.dead = numeric(nrow(mod1.df))
# for(i in 1:nrow(mod1.df)){mod1.df$n.dead[i] = sum(mod1.df$age[i] == sim1.data$age)}
# mod1.df$age.adj = mod1.df$age - min(mod1.df$age)
# mod1.df$s.vector = rep.int(0, nrow(mod1.df))
#
# priors = get_gomp_prior_shapes()
#
# mod1.list = list(
#   sample_ages = sim1.data$age,
#   Ndead = sum(mod1.df$n.dead),
#   age = mod1.df$age.adj,
#   prior_a_s1 = priors$shape1_a,
#   prior_a_s2 = priors$shape2_a,
#   prior_b_s1 = priors$shape1_b,
#   prior_b_s2 = priors$shape2_b,
#   N = nrow(mod1.df),
#   sum_dead = sum(mod1.df$n.dead),
#   age_error_sd = ifelse(sim1.data$age ==0, 0.01, (sim1.data$age+min(mod1.df$age))*0.05), # error based on real rather than adjusted age
#   s_vector = mod1.df$s.vector,
#   include_age_est_error = 0, # 0 for do not include 1 for include
#   include_sampling_error = 0, # 0 for do not include 1 for include
#   include_popchange_error = 0
# )
#
# stancode = get_marinesurvival_stancode()
# mod1 = stan(
#   model_code = stancode,
#   data = mod1.list,
#   chains = 4,
#   cores = 4,
#   iter = 2000,
#   init = list(list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)))
#
# )
#
# out = precis(mod1, digits = 4, depth = 2, pars = c("b", "a"))
# out
# sim1.in
#
# traceplot(mod1, pars = c("b", "a"))
# par(mfrow= c(1,1))
#
# age.seq = seq(min(mod1.df$age.adj), max(mod1.df$age.adj),1)
# post = extract.samples(mod1)
# modplot_samplescheck(age.seq = age.seq,
#                      post = post,
#                      input.list = mod1.list,
#                      staninputdf = mod1.df,
#                      simdata = sim1.data)
# modplot_curvecomp(sim.in = sim1.in, post = post, age.seq = age.seq)
#
#
# #########################################################
# # 4 Adding Age Estimation Error
# #######################################################
#
# sim4.in =
#   data.frame(
#     a = 0.01, # baseline mortality
#     b = 0.07, # ageing
#     rho = -0.01, # 0 = no change
#     e = 0.05, # standard deviation in of age estimation error
#     s =  1, # change in probability of being sampled above age S
#     ageS = 30
#   )
#
# sim4.list = run_simulation(sim.in = sim4.in,
#                            include_r = FALSE,
#                            include_e = TRUE,
#                            include_s = FALSE,
#                            sample_size = 200,
#                            age.maturity = 0 # if this isn't zero something in the mod1el intepretation
# )
# sim4.in = sim4.list$sim.in
# sim4.data = sim4.list$simulation
#
#
# mod4.df =
#   data.frame(age = seq(min(sim4.data$age), max(sim4.data$age)*1.5,1))
# mod4.df$n.dead = numeric(nrow(mod4.df))
# for(i in 1:nrow(mod4.df)){mod4.df$n.dead[i] = sum(mod4.df$age[i] == sim4.data$age)}
# mod4.df$age.adj = mod4.df$age - min(mod4.df$age)
# mod4.df$s.vector = ifelse(mod4.df$age.adj <sim4.in$ageS,0, 1 )
#
# # Without e parameter
# mod4a.list = list(
#   sample_ages = sim4.data$age,
#   Ndead = sum(mod4.df$n.dead),
#   age = mod4.df$age.adj,
#   prior_a_s1 = priors$shape1_a,
#   prior_a_s2 = priors$shape2_a,
#   prior_b_s1 = priors$shape1_b,
#   prior_b_s2 = priors$shape2_b,
#   N = nrow(mod4.df),
#   sum_dead = sum(mod4.df$n.dead),
#   age_error_sd = ifelse(sim4.data$age ==0, 0.01, (sim4.data$age+min(mod4.df$age))*0.05), # error based on real rather than adjusted age
#   s_vector = mod4.df$s.vector,
#   include_age_est_error = 0, # 0 for do not include 1 for include
#   include_sampling_error = 0, # 0 for do not include 1 for include
#   include_popchange_error = 0
# )
#
# stancode = get_marinesurvival_stancode()
# mod4a = stan(
#   model_code =stancode,
#   data = mod4a.list,
#   chains = 4,
#   cores = 4,
#   iter = 500,
#   init = list(list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)))
#
# )
#
#
# # # With e parameter
# mod4b.list = list(
#   sample_ages = sim4.data$age,
#   Ndead = sum(mod4.df$n.dead),
#   age = mod4.df$age.adj,
#   prior_a_s1 = priors$shape1_a,
#   prior_a_s2 = priors$shape2_a,
#   prior_b_s1 = priors$shape1_b,
#   prior_b_s2 = priors$shape2_b,
#   N = nrow(mod4.df),
#   sum_dead = sum(mod4.df$n.dead),
#   age_error_sd = ifelse(sim4.data$age ==0, 0.01, (sim4.data$age+min(mod4.df$age))*0.05), # error based on real rather than adjusted age
#   s_vector = mod4.df$s.vector,
#   include_age_est_error = 1, # 0 for do not include 1 for include
#   include_sampling_error = 0, # 0 for do not include 1 for include
#   include_popchange_error = 0
# )
#
# mod4b = stan(
#   model_code =stancode,
#   data = mod4b.list,
#   chains = 4,
#   cores = 4,
#   iter = 500,
#   init = list(list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)),
#               list(b = runif(1, max = 0.05), a = runif(1, max = 0.2)))
#
# )
#
#
# precis(mod4a, digits = 4,  depth = 1)
# precis(mod4b, digits = 4, depth = 1)
#
# age.seq = seq(min(mod4.df$age.adj), max(mod4.df$age.adj),1)
#
# post4a = extract.samples(mod4a)
# post4b = extract.samples(mod4b)
#
# grid.arrange(
#   modplot_samplescheck(age.seq = age.seq,
#                        post = post4a,
#                        input.list = mod4a.list,
#                        staninputdf = mod4.df,
#                        simdata = sim4.data),
#   modplot_samplescheck(age.seq = age.seq,
#                        post = post4b,
#                        input.list = mod4b.list,
#                        staninputdf = mod4.df,
#                        simdata = sim4.data),
#   ncol = 1
#
# )
#
# grid.arrange(
#   modplot_curvecomp(sim.in = sim4.in, post = post4a, age.seq = age.seq),
#   modplot_curvecomp(sim.in = sim4.in, post = post4b, age.seq = age.seq),
#   ncol = 1
# )
#
# grid.arrange(
#   modplot_ageXaccuracy(X = 0.5, post = post4a, age.seq = age.seq, sim.in = sim4.in),
#   modplot_ageXaccuracy(X = 0.5, post = post4b, age.seq = age.seq, sim.in = sim4.in),
#   ncol = 1
# )
#
#
