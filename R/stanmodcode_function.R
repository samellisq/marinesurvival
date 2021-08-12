#' Choose model and generate approiate code
#'
#'Becasue the code is an indexing nightmare need three versions
#'1. One species, one datasets
#'2. One species, many datasets
#'3. Many species, one or many datasets
#'
#'All are very similar but with slightly differetn indexing structure (and no
#'latent variable if only one species).
#'
#'This code chooses the right code for a given dataset
#'
#' @export

get_marinesurvival_stancode = function(mod.data.list){
  if(mod.data.list$Nspecies == 1 & mod.data.list$Ndatasets ==1){
    code = get_marinesurvival_onespecies.onedataset_stancode()
  } else {
    if(mod.data.list$Nspecies ==1 & mod.data.list$Ndatasets >1){
      code = get_marinesurvival_onespecies.manydatasets_stancode()
    } else {
      code =  get_marinesurvival_manyspecies.oneormanydatasets_stancode()
    }
  }
  return(code)
}


#' Generate a text file of the stan mode to be used
#'
#'@export
get_marinesurvival_onespecies.onedataset_stancode = function(){
  ##copied from Menopause Evolution II/onespecies stan code v1.stan on 10/08/2021
  code =
    "functions{
  // This function coverts a real value to an integer. A pain but necessarry
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
  int Nspecies;
  int sample_species[Nsamples];
  real prior_b_s2;
  real prior_b_s1;
  real prior_a_s2;
  real prior_a_s1;
  int include_age_est_error [Nspecies];
  int species_maxages [Nspecies];
  int include_samplebias_error [Nspecies];
  int BiasMat [Nages, Nspecies];
  int direction_samplebias [Nspecies];
  int include_popchange_error [Nspecies];
  int direction_popchange [Nspecies];
}

parameters{
  real<lower=0,upper=1> b;
  real<lower=0,upper=1> a;
  real <lower=-0.5,upper=0.5> r; // so this is r change over half the max lifespan. Currently both +ve and negative change but might have more info in some species.
  real <lower=0,upper=0.5> s; // at the moment it assumes that indiivuals above age S are more likely to be sampled but that can be changed
  vector[Nsamples] true_age; // there is no way to exclude this factor (or others) if not being used. Bit messy but makes no difference

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
  //int newdead_obs[N];
  int DeadMat [Nages, Nspecies];
  vector[Nsamples] ages_true;
  int done_indicator;
  int counter;

  // PHASE 1: ESTIMATE TRUE AGES
  //empty new count variable
  for( j in 1:Nspecies){
      for(i in 1:Nages){
    //newdead_obs[i] = 0;
      DeadMat[i,j] =0;
    }
  }

  //Do the fitting
  for(i in 1:Nsamples){
    if(include_age_est_error[sample_species[i]] >0){
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
        if(age[counter] == ages_true[i]){
          //newdead_obs[counter] = newdead_obs[counter]+1;
          DeadMat[counter, sample_species[i]] = DeadMat[counter, sample_species[i]] +1;
          done_indicator = 1;
        } else
          counter = counter +1;
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

      S = (BiasMat[i,1]*s)+1;

      num = lx*R*S;
      den = sigma_lx*R*S;
      theta[i, 1] = (num / den);
    }
  DeadMat[,1] ~ multinomial(theta[,1]);
}"

return(code)
}


#' Generate a text file of the stan mode to be used for the latent model adapted for if there is only one species
#'
#'@export
get_marinesurvival_onespecies.manydatasets_stancode= function(){
  # copied from Menopause Evolution II/ latent variable stan code v5c.stan on 10/08/2021
  code =
    "
  functions{
  // This function coverts a real value to an integer. A pain but necessarry
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
  int species_vector[Ndatasets];
  int sample_species[Nsamples];
  real prior_b_s2;
  real prior_b_s1;
  real prior_a_s2;
  real prior_a_s1;
  int include_age_est_error [Ndatasets];
  int species_maxages [Ndatasets];
  int include_samplebias_error [Ndatasets];
  int BiasMat [Nages, Ndatasets];
  int direction_samplebias [Ndatasets];
  int include_popchange_error [Ndatasets];
  int direction_popchange [Ndatasets];
}

parameters{
  real<lower=0,upper=1> b;
  real<lower=0,upper=1> a;
  vector <lower=-0.5,upper=0.5>[Ndatasets] r; // so this is r change over half the max lifespan. Currently both +ve and negative change but might have more info in some species.
  vector<lower=-2,upper=2>[Ndatasets] s;  // at the moment it assumes that indiivuals above age S are more likely to be sampled but that can be changed
  vector[Nsamples] true_age; // there is no way to exclude this factor (or others) if not being used. Bit messy but makes no difference

}

transformed parameters{
  vector[Ndatasets] rho;
  for(i in 1:Ndatasets){
    rho[i] = (r[i] / (species_maxages[i]))/2;
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
  //int newdead_obs[N];
  int DeadMat [Nages, Ndatasets];
  vector[Nsamples] ages_true;
  int done_indicator;
  int counter;

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
    if(include_age_est_error[sample_species[i]] >0){
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
        if(age[counter] == ages_true[i]){
          //newdead_obs[counter] = newdead_obs[counter]+1;
          DeadMat[counter, sample_species[i]] = DeadMat[counter, sample_species[i]] +1;
          done_indicator = 1;
        } else
          counter = counter +1;
      }
    }


  }


  // PHASE 2: FIT MODEL TO AGE COHORTS
  a ~ beta(prior_a_s1, prior_a_s2);

  //as above
  b ~ beta(prior_b_s1, prior_b_s2);

  // r ~ normal(0,1);
  //s ~ normal(0, 0.5);

  for(j in 1:Ndatasets){ //datasets

    if(direction_popchange[j] <0){
      r[j] ~ normal(-0.25,0.2);
    } else {
      if(direction_popchange[j] >0){
        r[j] ~ normal(0.25, 0.2);
      } else {
        r[j] ~ normal(0,1);
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

      //lx = exp(-1 * (a[j]/b[j]) * (exp(b[j] * AgeMat[i,j]) - 1));
      lx = exp(-1 * (a/b) * (exp(b * age[i]) - 1));
      //sigma_lx = sum(exp(-1 * (a[j]/b[j]) * (exp(b[j] * col(AgeMat,j)) - 1)));
      sigma_lx = sum(exp(-1 * (a/b) * (exp(b * age) - 1)));

      if(include_popchange_error[j] > 0){
        //R = (1-rho[j])^AgeMat[i,j];
        R = (1-rho[j])^age[i];
      } else{
        R = 1;
      }
      // if(include_samplebias_error[j] > 0){
      //   S = (BiasMat[i,j]*s[j])+1;
      // } else {
      //   S = 1;
      // } // unecessarry because Bia Mat will be 0 when needed
      S = (BiasMat[i,j]*s[j])+1;

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
  # copied from Menopause Evolution II/ latent variable stan code v5b.stan on 10/08/2021
  code =
    "
  functions{
  // This function coverts a real value to an integer. A pain but necessarry
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
  int species_vector[Ndatasets];
  int sample_species[Nsamples];
  real prior_b_s2;
  real prior_b_s1;
  real prior_a_s2;
  real prior_a_s1;
  int include_age_est_error [Ndatasets];
  int species_maxages [Ndatasets];
  int include_samplebias_error [Ndatasets];
  int BiasMat [Nages, Ndatasets];
  int direction_samplebias [Ndatasets];
  int include_popchange_error [Ndatasets];
  int direction_popchange [Ndatasets];
}

parameters{
  vector<lower=0,upper=1>[Nspecies] a; // this needs to go to datasets
  real<lower=0,upper=1> abar; // For all whales
  real<lower = 0> abar_phi; // for all whales

  //same as above
  vector<lower=0,upper=1>[Nspecies] b;
  real<lower=0,upper=1> bbar;
  real<lower = 0> bbar_phi;

  vector <lower=-0.5,upper=0.5>[Ndatasets] r; // so this is r change over half the max lifespan. Currently both +ve and negative change but might have more info in some species.
  vector<lower=-2,upper=2>[Ndatasets] s;
  vector[Nsamples] true_age; // there is no way to exclude this factor (or others) if not being used. Bit messy but makes no difference

}

transformed parameters{
  vector[Ndatasets] rho;
  real abar_theta;
  real bbar_theta;
  for(i in 1:Ndatasets){
    rho[i] = (r[i] / (species_maxages[i]))/2;
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
  //int newdead_obs[N];
  int DeadMat [Nages, Ndatasets]; // datasets
  vector[Nsamples] ages_true;
  int done_indicator;
  int counter;
  int k; //species index

  // PHASE 1: ESTIMATE TRUE AGES
  //empty new count variable
  for( j in 1:Ndatasets){ // datsets
      for(i in 1:Nages){
    //newdead_obs[i] = 0;
      DeadMat[i,j] =0;
    }
  }

  //Do the fitting
  for(i in 1:Nsamples){
    if(include_age_est_error[sample_species[i]] >0){
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
        if(age[counter] == ages_true[i]){
          //newdead_obs[counter] = newdead_obs[counter]+1;
          DeadMat[counter, sample_species[i]] = DeadMat[counter, sample_species[i]] +1; //datasets
          done_indicator = 1;
        } else
          counter = counter +1;
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

  // r ~ normal(0,1);
  //s ~ normal(0, 0.5);

  for(j in 1:Ndatasets){ //datasets

    if(direction_popchange[j] <0){
      r[j] ~ normal(-0.25,0.2);
    } else {
      if(direction_popchange[j] >0){
        r[j] ~ normal(0.25, 0.2);
      } else {
        r[j] ~ normal(0,1);
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

      //lx = exp(-1 * (a[j]/b[j]) * (exp(b[j] * AgeMat[i,j]) - 1));
      lx = exp(-1 * (a[k]/b[k]) * (exp(b[k] * age[i]) - 1));
      //sigma_lx = sum(exp(-1 * (a[j]/b[j]) * (exp(b[j] * col(AgeMat,j)) - 1)));
      sigma_lx = sum(exp(-1 * (a[k]/b[k]) * (exp(b[k] * age) - 1)));

      if(include_popchange_error[j] > 0){
        //R = (1-rho[j])^AgeMat[i,j];
        R = (1-rho[j])^age[i];
      } else{
        R = 1;
      }
      // if(include_samplebias_error[j] > 0){
      //   S = (BiasMat[i,j]*s[j])+1;
      // } else {
      //   S = 1;
      // } // unecessarry because Bia Mat will be 0 when needed
      S = (BiasMat[i,j]*s[j])+1;

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
