#' Choose model and generate appropriate stan code
#'
#'From the same modelling framework different modelling structures are needed for
#'different combinations are required depending on if (a) one or many datasets
#'are being analysed; (b) if those datasets belong to one or many species and (c)
#'if for a given species there are data present for one or both sexes.
#'
#'The framework for each model is identical, but each requires different indexing
#'and latent variable structure. See https://doi.org/10.1101/2023.02.22.529527 for
#'details.
#'
#'This function simply chooses the right code for a given dataset
#'
#' @export

get_marinesurvival_stancode = function(input){


  #if the input is a model input list
  if(is.list(input)){

    mod.data.list = input
    if(length(unique(mod.data.list$species_sex_vector)) == 1){


      if(mod.data.list$Nspecies == 1 & mod.data.list$Ndatasets ==1){
        code = get_marinesurvival_onespecies.onedataset_stancode()
        print("One Species, One Dataset, One Sex")
      } else {
        if(mod.data.list$Nspecies ==1 & mod.data.list$Ndatasets >1){
          code = get_marinesurvival_onespecies.manydatasets_stancode()
          print("One Species, Many Datasets, One Sex")
        } else {
          code =  get_marinesurvival_manyspecies.oneormanydatasets_stancode()
          print("Many Species, Many Datasets, One Sex")
        }
      }

    } else {

      if(mod.data.list$Nspecies == 2){
        code =  get_marinesurvival_onespecies.oneormanydatasets_twosexes_stancode()
        print("One Species, Two Sexes")
      } else {
        code =  get_marinesurvival_manyspecies.oneormanydatasets_twosexes_stancode()
        print("Many Species, Two Sexes")
      }


    }

  } else{

    valid = 0
    if(input == "One Species, One Dataset, One Sex"){
      code = get_marinesurvival_onespecies.onedataset_stancode()
      print("One Species, One Dataset, One Sex")
      valid = 1
    }
    if(input == "One Species, Many Datasets, One Sex"){
      code = get_marinesurvival_onespecies.manydatasets_stancode()
      print("One Species, Many Datasets, One Sex")
      valid = 1
    }
    if(input == "Many Species, Many Datasets, One Sex"){
      code = get_marinesurvival_manyspecies.oneormanydatasets_stancode()
      print("Many Species, Many Datasets, One Sex")
      valid = 1
    }
    if(input == "One Species, Two Sexes"){
      code = get_marinesurvival_onespecies.oneormanydatasets_twosexes_stancode()
      print("One Species, Two Sexes")
      valid = 1
    }
    if(input == "Many Species, Two Sexes"){
      code = get_marinesurvival_manyspecies.oneormanydatasets_twosexes_stancode()
      print("Many Species, Two Sexes")
      valid = 1
    }
    if(valid == 0){
      warning("Model text input not correct. Must be of the form:
              One species, Many Datasets, One Sex | One Species, Two Sexes")
    }



  }

  return(code)
}


#' Generate a text file of the stan mode to be used
#'
#'@export
get_marinesurvival_onespecies.onedataset_stancode = function(){

  code =
    "
    functions{
  // This function coverts a real value to an integer.
  // from https://discourse.mc-stan.org/t/real-to-integer-conversion/5622/8
  int bin_search(real x, int min_val, int max_val){
    int range = (max_val - min_val+1)/2;
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2;
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range;
      }
    }
    return out;
  }

}

data{
  int Nages;
  int Nsamples;
  vector[Nages] age;
  vector[Nsamples] sample_ages;
  vector[Nsamples] age_error_sd;
  int Ndatasets;
  int Nspecies;
  int Npopulations;
  int species_vector[Ndatasets];
  int population_vector[Ndatasets];
  int dataset_vector[Nsamples];
  int species_sex_vector[Nspecies];
  real prior_b_s2;
  real prior_b_s1;
  real prior_a_s2;
  real prior_a_s1;
  int include_age_est_error [Ndatasets]; // development option to turn off age estimation error part of model to speed up or test the model.
  int species_maxages [Ndatasets];
  int include_samplebias_error [Ndatasets];
  int BiasMat [Nages, Ndatasets];
  int direction_samplebias [Ndatasets];
  int include_popchange_error [Ndatasets];
  int direction_popchange [Npopulations];
}

parameters{
  real<lower=0,upper=1> b;
  real<lower=0,upper=1> a;
  real <lower=-0.5,upper=0.5> r; // so this is r change over half the max lifespan.
  real <lower=0,upper=0.5> s;
  vector[Nsamples] true_age;

}

transformed parameters{
  real rho;
  rho = r / species_maxages[1]/2;
}

model{
  matrix[Nages, Nspecies] theta;
  real num;
  real den;
  real lx;
  real sigma_lx;
  real R;
  real S;
  int DeadMat [Nages, Nspecies];
  vector[Nsamples] ages_true;
  int done_indicator;
  int counter;

  // PHASE 1: ESTIMATE TRUE AGES
  //empty new count variable
  for( j in 1:Nspecies){
      for(i in 1:Nages){
      DeadMat[i,j] =0;
    }
  }

  //Do the fitting
  for(i in 1:Nsamples){
    if(include_age_est_error[dataset_vector[i]] >0){
      true_age ~ normal(sample_ages, age_error_sd);
      ages_true[i] = bin_search(round(true_age[i]), 0, Nages);
    } else {
      true_age ~ normal(sample_ages, 1);
      ages_true[i] = sample_ages[i];
    }
    //Recount age cohorts
    if(ages_true[i] >= 0){ // Becasue error means that ages can drop out
      done_indicator = 0;
      counter = 1; // note counter is one ahead of age number becasue of age 0
      while(done_indicator == 0){ // while to stop it going through each age every time
      if(counter > Nages){ // stop going up the counter if we are beyond the maximum possible age. Deals with initalisation issue (initalising beyond the range)
      done_indicator = 1;
      } else {
        if(age[counter] == ages_true[i]){
          DeadMat[counter, dataset_vector[i]] = DeadMat[counter, dataset_vector[i]] +1;
          done_indicator = 1;
        } else
          counter = counter +1;
      }
    }
  }


  }


  // PHASE 2: FIT MODEL TO AGE COHORTS
  a ~ beta( prior_a_s1 , prior_a_s2 );
  b ~ beta( prior_b_s1 , prior_b_s2 );

  if(direction_popchange[1] <0){
      r ~ normal(-0.25,0.2);
    } else {
      if(direction_popchange[1] >0){
        r ~ normal(0.25, 0.2);
      } else {
        r ~ normal(0,1);
      }
    }
    if(direction_samplebias[1] <0){
      s ~ normal(-1,0.5);
    } else {
      if(direction_samplebias[1] >0){
        s ~ normal(1, 0.5);
      } else {
        s ~ normal(0,2);
      }
    }

      for ( i in 1:Nages ) {
      lx = exp(-1 * (a/b) * (exp(b * age[i]) - 1));
      sigma_lx = sum(exp(-1 * (a/b) * (exp(b * age) - 1)));

      if(include_popchange_error[1] > 0){
        R = (1-rho)^age[i];
      } else{
        R = 1;
      }

      S = (BiasMat[i,1]*s)+1; // Bias Mat will be 0 if not in use so will have no effect

      num = lx*R*S;
      den = sigma_lx*R*S;
      theta[i, 1] = (num / den);
    }
  DeadMat[,1] ~ multinomial(theta[,1]);
}
"

return(code)
}


#' Generate a text file of the stan mode to be used for the latent model adapted for if there is only one species
#'
#'@export
get_marinesurvival_onespecies.manydatasets_stancode= function(){

  code =
    "
  functions{
  // This function coverts a real value to an integer.
  // from https://discourse.mc-stan.org/t/real-to-integer-conversion/5622/8
  int bin_search(real x, int min_val, int max_val){
    int range = (max_val - min_val+1)/2;
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2;
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range;
      }
    }
    return out;
  }

}

data{
  int Nages;
  int Nsamples;
  vector[Nages] age;
  vector[Nsamples] sample_ages;
  vector[Nsamples] age_error_sd;
  int Ndatasets;
  int Nspecies;
  int Npopulations;
  int species_vector[Ndatasets];
  int population_vector[Ndatasets];
  int dataset_vector[Nsamples];
  int species_sex_vector[Nspecies];
  real prior_b_s2;
  real prior_b_s1;
  real prior_a_s2;
  real prior_a_s1;
  int include_age_est_error [Ndatasets]; // development option to turn off age estimation error part of model to speed up or test the model.
  int species_maxages [Ndatasets];
  int include_samplebias_error [Ndatasets];
  int BiasMat [Nages, Ndatasets];
  int direction_samplebias [Ndatasets];
  int include_popchange_error [Ndatasets];
  int direction_popchange [Npopulations];
}

parameters{
  real<lower=0,upper=1> b;
  real<lower=0,upper=1> a;
  vector <lower=-0.5,upper=0.5>[Npopulations] r; // so this is r change over half the max lifespan.
  vector<lower=-2,upper=2>[Ndatasets] s;
  vector[Nsamples] true_age;

}

transformed parameters{
  vector[Ndatasets] rho;
  for(i in 1:Ndatasets){
    rho[population_vector[i]] = (r[population_vector[i]] / (species_maxages[i]))/2;
  }
}

model{
  matrix[Nages, Ndatasets] theta;
  real num;
  real den;
  real lx;
  real sigma_lx;
  real R;
  real S;
  int DeadMat [Nages, Ndatasets];
  vector[Nsamples] ages_true;
  int done_indicator;
  int counter;

  // PHASE 1: ESTIMATE TRUE AGES
  //empty new count variable
  for( j in 1:Ndatasets){
      for(i in 1:Nages){
      DeadMat[i,j] =0;
    }
  }

  //Do the fitting
  for(i in 1:Nsamples){
    if(include_age_est_error[dataset_vector[i]] >0){
      true_age ~ normal(sample_ages, age_error_sd);
      ages_true[i] = bin_search(round(true_age[i]), 0, Nages);
    } else {
      true_age ~ normal(sample_ages, 1);
      ages_true[i] = sample_ages[i];
    }
    //Recount age cohorts
    if(ages_true[i] >= 0){ // Becasue error means that ages can drop out
      done_indicator = 0;
      counter = 1; // note counter is one ahead of age number becasue of age 0
      while(done_indicator == 0){ // while to stop it going through each age every time
      if(counter > Nages){ // stop going up the counter if we are beyond the maximum possible age. Deals with initalisation issue (initalising beyond the range)
      done_indicator = 1;
      } else {
        if(age[counter] == ages_true[i]){
          DeadMat[counter, dataset_vector[i]] = DeadMat[counter, dataset_vector[i]] +1;
          done_indicator = 1;
        } else
          counter = counter +1;
      }
    }
  }

  }


  // PHASE 2: FIT MODEL TO AGE COHORTS
  a ~ beta(prior_a_s1, prior_a_s2);

  //as above
  b ~ beta(prior_b_s1, prior_b_s2);

  for(j in 1:Ndatasets){ //datasets

    if(direction_popchange[population_vector[j]] <0){
      r[population_vector[j]] ~ normal(-0.25,0.2);
    } else {
      if(direction_popchange[population_vector[j]] >0){
        r[population_vector[population_vector[j]]] ~ normal(0.25, 0.2);
      } else {
        r[population_vector[j]] ~ normal(0,1);
      }
    }
    if(direction_samplebias[j] <0){
      s[j] ~ normal(-1,0.5);
    } else {
      if(direction_samplebias[j] >0){
        s[j] ~ normal(1, 0.5);
      } else {
        s[j] ~ normal(0,2);
      }
    }


    for ( i in 1:Nages ) {
      lx = exp(-1 * (a/b) * (exp(b * age[i]) - 1));
      sigma_lx = sum(exp(-1 * (a/b) * (exp(b * age) - 1)));

      if(include_popchange_error[j] > 0){
        R = (1-rho[population_vector[j]])^age[i];
      } else{
        R = 1;
      }
      S = (BiasMat[i,j]*s[j])+1; // Bias Mat will be 0 if not in use so will have no effect

      num = lx*R*S;
      den = sigma_lx*R*S;
      theta[i, j] = (num / den);
    }
    DeadMat[,j] ~ multinomial(theta[,j]);
  }

}


"
return(code)
}



#' Generate a text file of the stan mode to be used for the latent model
#'
#'@export
get_marinesurvival_manyspecies.oneormanydatasets_stancode= function(){

  code =
    "
  functions{
  // This function coverts a real value to an integer.
  // from https://discourse.mc-stan.org/t/real-to-integer-conversion/5622/8
  int bin_search(real x, int min_val, int max_val){
    int range = (max_val - min_val+1)/2;
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2;
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range;
      }
    }
    return out;
  }

}

data{
  int Nages;
  int Nsamples;
  vector[Nages] age;
  vector[Nsamples] sample_ages;
  vector[Nsamples] age_error_sd;
  int Ndatasets;
  int Nspecies;
  int Npopulations;
  int species_vector[Ndatasets];
  int population_vector[Ndatasets];
  int dataset_vector[Nsamples];
  int species_sex_vector[Nspecies];
  real prior_b_s2;
  real prior_b_s1;
  real prior_a_s2;
  real prior_a_s1;
  int include_age_est_error [Ndatasets]; // development option to turn off age estimation error part of model to speed up or test the model.
  int species_maxages [Ndatasets];
  int include_samplebias_error [Ndatasets];
  int BiasMat [Nages, Ndatasets];
  int direction_samplebias [Ndatasets];
  int include_popchange_error [Ndatasets];
  int direction_popchange [Npopulations];
}

parameters{
  vector<lower=0,upper=1>[Nspecies] a; // this needs to go to datasets
  real<lower=0,upper=1> abar; // For all whales
  real<lower = 0> abar_phi; // for all whales

  //same as above
  vector<lower=0,upper=1>[Nspecies] b;
  real<lower=0,upper=1> bbar;
  real<lower = 0> bbar_phi;

  vector <lower=-0.5,upper=0.5>[Ndatasets] r; // so this is r change over half the max lifespan.
  vector<lower=-2,upper=2>[Ndatasets] s;
  vector[Nsamples] true_age;

}

transformed parameters{
  vector[Ndatasets] rho;
  real abar_theta;
  real bbar_theta;
  for(i in 1:Ndatasets){
    rho[population_vector[i]] = (r[population_vector[i]] / (species_maxages[i]))/2;
  }
  abar_theta = abar_phi +2;
  bbar_theta = bbar_phi +2;

}

model{
  matrix[Nages, Ndatasets] theta;// datasets
  real num;
  real den;
  real lx;
  real sigma_lx;
  real R;
  real S;
  int DeadMat [Nages, Ndatasets]; // datasets
  vector[Nsamples] ages_true;
  int done_indicator;
  int counter;
  int k; //species index

  // PHASE 1: ESTIMATE TRUE AGES
  //empty new count variable
  for( j in 1:Ndatasets){ // datsets
      for(i in 1:Nages){
      DeadMat[i,j] =0;
    }
  }

  //Do the fitting
  for(i in 1:Nsamples){
    if(include_age_est_error[dataset_vector[i]] >0){
      true_age ~ normal(sample_ages, age_error_sd);
      ages_true[i] = bin_search(round(true_age[i]), 0, Nages);
    } else {
      true_age ~ normal(sample_ages, 1);
      ages_true[i] = sample_ages[i];
    }
    //Recount age cohorts
    if(ages_true[i] >= 0){ // Becasue error means that ages can drop out
      done_indicator = 0;
      counter = 1; // note counter is one ahead of age number becasue of age 0
      while(done_indicator == 0){ // while to stop it going through each age every time
      if(counter > Nages){ // stop going up the counter if we are beyond the maximum possible age. Deals with initalisation issue (initalising beyond the range)
      done_indicator = 1;
      } else {
        if(age[counter] == ages_true[i]){
          DeadMat[counter, dataset_vector[i]] = DeadMat[counter, dataset_vector[i]] +1;
          done_indicator = 1;
        } else
          counter = counter +1;
      }
    }
  }


  }


  // PHASE 2: FIT MODEL TO AGE COHORTS
  abar ~ beta(prior_a_s1, prior_a_s2); // so this the mean whale is taken from a distribution given by the real data
  abar_phi ~ exponential(1); // keep these two about whales in general
  //new layer here for asp and asp_phi
  a ~ beta( abar*abar_theta , (1-abar)*abar_theta );

  //as above
  bbar ~ beta (prior_b_s1, prior_b_s2);
  bbar_phi ~ exponential(1);
  b ~ beta( bbar*bbar_theta , (1-bbar)*bbar_theta );

  for(j in 1:Ndatasets){ //datasets

    if(direction_popchange[population_vector[j]] <0){
      r[population_vector[j]] ~ normal(-0.25,0.2);
    } else {
      if(direction_popchange[population_vector[j]] >0){
        r[population_vector[j]] ~ normal(0.25, 0.2);
      } else {
        r[population_vector[j]] ~ normal(0,1);
      }
    }
    if(direction_samplebias[j] <0){
      s[j] ~ normal(-1,0.5);
    } else {
      if(direction_samplebias[j] >0){
        s[j] ~ normal(1, 0.5);
      } else {
        s[j] ~ normal(0,2);
      }
    }

    k = species_vector[j];

    for ( i in 1:Nages ) {
      lx = exp(-1 * (a[k]/b[k]) * (exp(b[k] * age[i]) - 1));
      sigma_lx = sum(exp(-1 * (a[k]/b[k]) * (exp(b[k] * age) - 1)));

      if(include_popchange_error[j] > 0){
        R = (1-rho[population_vector[j]])^age[i];
      } else{
        R = 1;
      }
      S = (BiasMat[i,j]*s[j])+1; // Bias Mat will be 0 if not in use so will have no effect

      num = lx*R*S;
      den = sigma_lx*R*S;
      theta[i, j] = (num / den);
    }
    DeadMat[,j] ~ multinomial(theta[,j]);
  }

}


"
return(code)
}



get_marinesurvival_onespecies.oneormanydatasets_twosexes_stancode = function(){
  "
  functions{
  // This function coverts a real value to an integer.
  // from https://discourse.mc-stan.org/t/real-to-integer-conversion/5622/8
  int bin_search(real x, int min_val, int max_val){
    int range = (max_val - min_val+1)/2;
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2;
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range;
      }
    }
    return out;
  }

}

data{
  int Nages;
  int Nsamples;
  vector[Nages] age;
  vector[Nsamples] sample_ages;
  vector[Nsamples] age_error_sd;
  int Ndatasets;
  int Nspecies;
  int Npopulations;
  int species_vector[Ndatasets];
  int population_vector[Ndatasets];
  int dataset_vector[Nsamples];
  int species_sex_vector[Nspecies];
  real prior_b_s2;
  real prior_b_s1;
  real prior_a_s2;
  real prior_a_s1;
  int include_age_est_error [Ndatasets]; // development option to turn off age estimation error part of model to speed up or test the model.
  int species_maxages [Ndatasets];
  int include_samplebias_error [Ndatasets];
  int BiasMat [Nages, Ndatasets];
  int direction_samplebias [Ndatasets];
  int include_popchange_error [Ndatasets];
  int direction_popchange [Npopulations];
}

parameters{
  vector<lower=0,upper=1>[2] b;
  vector<lower=0,upper=1>[2] a;
  vector <lower=-0.5,upper=0.5>[Npopulations] r; // so this is r change over half the max lifespan.
  vector<lower=-2,upper=2>[Ndatasets] s;
  vector[Nsamples] true_age;

}

transformed parameters{
  vector[Npopulations] rho;
  for(i in 1:Ndatasets){
    rho[population_vector[i]] = (r[population_vector[i]] / (species_maxages[i]))/2;
  }
}

model{
  matrix[Nages, Ndatasets] theta;
  real num;
  real den;
  real lx;
  real sigma_lx;
  real R;
  real S;
  int DeadMat [Nages, Ndatasets];
  vector[Nsamples] ages_true;
  int done_indicator;
  int counter;
  int k; //sex index

  // PHASE 1: ESTIMATE TRUE AGES
  //empty new count variable
  for( j in 1:Ndatasets){
      for(i in 1:Nages){
    //newdead_obs[i] = 0;
      DeadMat[i,j] =0;
    }
  }

  //Do the fitting
  for(i in 1:Nsamples){
    if(include_age_est_error[dataset_vector[i]] >0){
      true_age ~ normal(sample_ages, age_error_sd);
      ages_true[i] = bin_search(round(true_age[i]), 0, Nages);
    } else {
      true_age ~ normal(sample_ages, 1);
      ages_true[i] = sample_ages[i];
    }
    //Recount age cohorts
    if(ages_true[i] >= 0){ // Becasue error means that ages can drop out
      done_indicator = 0;
      counter = 1; // note counter is one ahead of age number becasue of age 0
      while(done_indicator == 0){ // while to stop it going through each age every time
      if(counter > Nages){ // stop going up the counter if we are beyond the maximum possible age. Deals with initalisation issue (initalising beyond the range)
      done_indicator = 1;
      } else {
        if(age[counter] == ages_true[i]){
          DeadMat[counter, dataset_vector[i]] = DeadMat[counter, dataset_vector[i]] +1;
          done_indicator = 1;
        } else
          counter = counter +1;
      }
    }
  }


  }


  // PHASE 2: FIT MODEL TO AGE COHORTS
  a ~ beta(prior_a_s1, prior_a_s2);

  //as above
  b ~ beta(prior_b_s1, prior_b_s2);


  for(j in 1:Ndatasets){ //datasets

    if(direction_popchange[population_vector[j]] <0){
      r[population_vector[j]] ~ normal(-0.25,0.2);
    } else {
      if(direction_popchange[population_vector[j]] >0){
        r[population_vector[population_vector[j]]] ~ normal(0.25, 0.2);
      } else {
        r[population_vector[j]] ~ normal(0,1);
      }
    }
    if(direction_samplebias[j] <0){
      s[j] ~ normal(-1,0.5);
    } else {
      if(direction_samplebias[j] >0){
        s[j] ~ normal(1, 0.5);
      } else {
        s[j] ~ normal(0,2);
      }
    }

	k = species_vector[j];

    for ( i in 1:Nages ) {
      lx = exp(-1 * (a[k]/b[k]) * (exp(b[k] * age[i]) - 1));
      sigma_lx = sum(exp(-1 * (a[k]/b[k]) * (exp(b[k] * age) - 1)));

      if(include_popchange_error[j] > 0){
        R = (1-rho[population_vector[j]])^age[i];
      } else{
        R = 1;
      }
      S = (BiasMat[i,j]*s[j])+1; // Bias Mat will be 0 if not in use so will have no effect

      num = lx*R*S;
      den = sigma_lx*R*S;
      theta[i, j] = (num / den);
    }
    DeadMat[,j] ~ multinomial(theta[,j]);
  }

}



  "
}


get_marinesurvival_manyspecies.oneormanydatasets_twosexes_stancode = function(){
  code =
  "
  functions{
  // This function coverts a real value to an integer.
  // from https://discourse.mc-stan.org/t/real-to-integer-conversion/5622/8
  int bin_search(real x, int min_val, int max_val){
    int range = (max_val - min_val+1)/2;
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2;
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range;
      }
    }
    return out;
  }

}

data{
  int Nages;
  int Nsamples;
  vector[Nages] age;
  vector[Nsamples] sample_ages;
  vector[Nsamples] age_error_sd;
  int Ndatasets;
  int Nspecies;
  int Npopulations;
  int species_vector[Ndatasets];
  int population_vector[Ndatasets];
  int dataset_vector[Nsamples];
  int species_sex_vector[Nspecies];
  real prior_b_s2;
  real prior_b_s1;
  real prior_a_s2;
  real prior_a_s1;
  int include_age_est_error [Ndatasets]; // development option to turn off age estimation error part of model to speed up or test the model.
  int species_maxages [Ndatasets];
  int include_samplebias_error [Ndatasets];
  int BiasMat [Nages, Ndatasets];
  int direction_samplebias [Ndatasets];
  int include_popchange_error [Ndatasets];
  int direction_popchange [Npopulations];
}

parameters{
  vector<lower=0,upper=1>[Nspecies] a; // this needs to go to datasets
  vector<lower=0,upper=1>[2] abar; // For all whales
  vector<lower = 0>[2] abar_phi; // for all whales

  //same as above
  vector<lower=0,upper=1>[Nspecies] b;
  vector<lower=0,upper=1>[2] bbar;
  vector<lower = 0>[2] bbar_phi;

  vector <lower=-0.5,upper=0.5>[Npopulations] r; // so this is r change over half the max lifespan.
  vector<lower=-2,upper=2>[Ndatasets] s;
  vector[Nsamples] true_age;

}

transformed parameters{
  vector[Npopulations] rho;
  real abar_theta[2];
  real bbar_theta[2];
  for(i in 1:Ndatasets){
    rho[population_vector[i]] = (r[population_vector[i]] / (species_maxages[i]))/2;
  }

	abar_theta[1] = abar_phi[1] +2;
	bbar_theta[1] = bbar_phi[1] +2;
	abar_theta[2] = abar_phi[2] +2;
	bbar_theta[2] = bbar_phi[2] +2;


}

model{
  matrix[Nages, Ndatasets] theta;// datasets
  real num;
  real den;
  real lx;
  real sigma_lx;
  real R;
  real S;
  int DeadMat [Nages, Ndatasets]; // datasets
  vector[Nsamples] ages_true;
  int done_indicator;
  int counter;
  int k; //species index

  // PHASE 1: ESTIMATE TRUE AGES
  //empty new count variable
  for( j in 1:Ndatasets){ // datsets
      for(i in 1:Nages){
      DeadMat[i,j] =0;
    }
  }

  //Do the fitting
  for(i in 1:Nsamples){
    if(include_age_est_error[dataset_vector[i]] >0){
      true_age ~ normal(sample_ages, age_error_sd);
      ages_true[i] = bin_search(round(true_age[i]), 0, Nages);
    } else {
      true_age ~ normal(sample_ages, 1);
      ages_true[i] = sample_ages[i];
    }
    //Recount age cohorts
    if(ages_true[i] >= 0){ // Becasue error means that ages can drop out
      done_indicator = 0;
      counter = 1; // note counter is one ahead of age number becasue of age 0
      while(done_indicator == 0){ // while to stop it going through each age every time
      if(counter > Nages){ // stop going up the counter if we are beyond the maximum possible age. Deals with initalisation issue (initalising beyond the range)
      done_indicator = 1;
      } else {
        if(age[counter] == ages_true[i]){
          DeadMat[counter, dataset_vector[i]] = DeadMat[counter, dataset_vector[i]] +1;
          done_indicator = 1;
        } else
          counter = counter +1;
      }
    }
  }


  }


  // PHASE 2: FIT MODEL TO AGE COHORTS
  abar ~ beta(prior_a_s1, prior_a_s2); // so this the mean whale is taken from a distribution given by the real data
  abar_phi ~ exponential(1); // keep these two about whales in general
  bbar ~ beta (prior_b_s1, prior_b_s2);
  bbar_phi ~ exponential(1);
  //new layer here for asp and asp_phi
  for(i in 1:Nspecies){
	a ~ beta( abar[species_sex_vector[i]]*abar_theta[species_sex_vector[i]] , (1-abar[species_sex_vector[i]])*abar_theta[species_sex_vector[i]] );
	 b ~ beta( bbar[species_sex_vector[i]]*bbar_theta[species_sex_vector[i]] , (1-bbar[species_sex_vector[i]])*bbar_theta[species_sex_vector[i]] );
  }


  //as above

  for(j in 1:Ndatasets){ //datasets

    if(direction_popchange[population_vector[j]] <0){
      r[population_vector[j]] ~ normal(-0.25,0.2);
    } else {
      if(direction_popchange[population_vector[j]] >0){
        r[population_vector[j]] ~ normal(0.25, 0.2);
      } else {
        r[population_vector[j]] ~ normal(0,1);
      }
    }
    if(direction_samplebias[j] <0){
      s[j] ~ normal(-1,0.5);
    } else {
      if(direction_samplebias[j] >0){
        s[j] ~ normal(1, 0.5);
      } else {
        s[j] ~ normal(0,2);
      }
    }

    k = species_vector[j];

    for ( i in 1:Nages ) {

      lx = exp(-1 * (a[k]/b[k]) * (exp(b[k] * age[i]) - 1));
      sigma_lx = sum(exp(-1 * (a[k]/b[k]) * (exp(b[k] * age) - 1)));

      if(include_popchange_error[j] > 0){
        R = (1-rho[population_vector[j]])^age[i];
      } else{
        R = 1;
      }
      S = (BiasMat[i,j]*s[j])+1; // Bias Mat will be 0 if not in use so will have no effect

      num = lx*R*S;
      den = sigma_lx*R*S;
      theta[i, j] = (num / den);
    }
    DeadMat[,j] ~ multinomial(theta[,j]);
  }

}
"
return(code)
}
