
data {
  int<lower=1> Nobs;  // number of observations
  int pres[Nobs];  // response variable
  int<lower=1> NX;  // number of initial population-level effects in X
  int<lower=1> Ntrait;  // number of trait-productivity interactions
  matrix[Nobs, NX] X;  // population-level design matrix
  
  
  // data for group-level effects of species
  int<lower=1> Nplot;  // number of grouping levels
  int<lower=1> Mplot;  // number of coefficients per level. makes is possible to vectorize sd, thus much faster
  int<lower=1> plot[Nobs];  // grouping indicator per observation


  // data for group-level effects of species
  int<lower=1> Nspecies;  // number of grouping levels
  int<lower=1> Mspecies;  // number of coefficients per level. makes is possible to vectorize sd, thus much faster
  int<lower=1> species[Nobs];  // grouping indicator per observation


  // latent productivity
  int<lower=1> Nlatbm; // number of biomass observations
  int<lower=1> plotbm[Nlatbm]; // latent variables ID variable
  vector[Nlatbm] latbm; // measured biomass
  vector[Nplot] meansbm;  // biomass means as informative priors
  vector[Nplot] sdbm;  // biomass SDs as informative priors
  
}

transformed data {
  
  matrix[Nobs, NX] c_X;  // centered version of X without an intercept
  vector[NX] means_X;  // column means of X before centering
 
  for (i in 1:NX) {
    means_X[i] = mean(X[, i]);
    c_X[, i] = X[, i] - means_X[i];
  }
}

parameters {
  vector[NX + 1 + Ntrait] b;  // population-level effects, plus productivity and 
                              // trait-productivity interactions

  real c_Int;  // temporary intercept for centered predictors
  vector<lower=0>[Mplot] sd_plot;  // group-level standard deviations
  vector[Nplot] z_plot[Mplot];  // standardized group-level effects
  vector<lower=0>[Mspecies] sd_species;  // group-level standard deviations
  vector[Nspecies] z_species[Mspecies];  // standardized group-level effects
  
  
  // latent productivity
  vector<lower=0>[Nplot] realbm; // define the 'real' latent productivity, constrained to be positive
  real<lower=0> sigma_bm; // define sigmas of latents to be non-negative

}

transformed parameters {
  
  // set matrix and means to merge observed and latent predictor variable for delta15 response
  // transform latent to log and log-centered variable


  // create new productivity vector with length of predictor matrix
  vector[Nobs] log_realbmX;
  real sd_log_realbmX;
  vector[Nobs] log_realbmX_scl;

  
  // create new variables: centred productivity and interactions with productivity
  vector[Nobs] c_log_realbmX_scl;
  real mean_log_realbm_scl;

  vector[NX + 1] means_Xtmp;
  vector[NX + 1 + Ntrait] means_Xall;
  
  matrix[Nobs, NX + 1] c_Xtmp;
  matrix[Nobs, NX + 1 + Ntrait] c_Xall;


  matrix[Nobs, Ntrait] c_X_int;
  vector[Ntrait] means_X_interact;
  vector[Ntrait] sd_X_interact;
  
  
  // transform group-level effects
  vector[Nplot] r_plot;  // actual group-level effects
  vector[Nspecies] r_species;  // actual group-level effects
  
  r_plot = (sd_plot[1] * (z_plot[1]));
  r_species = (sd_species[1] * (z_species[1]));




  // map latent variable of group level plot to length of observations
  for (i in 1:Nobs) {
    log_realbmX[i] = log(realbm[ plot[i] ] ) ;
    };

  // get log productivity standard deviation to scale this parameter
  sd_log_realbmX = sd(log_realbmX);
  log_realbmX_scl = log_realbmX/sd_log_realbmX;

  // get mean productivity and center
  mean_log_realbm_scl = mean(log_realbmX_scl);
  c_log_realbmX_scl = log_realbmX_scl - mean_log_realbm_scl;

  

  // get interaction of traits with latent productivity and center those
  for(i in 1: Ntrait) {
    sd_X_interact[i] = sd(X[ , i] .* log_realbmX_scl);
    // vector .* vector multiplies element-wise 
    means_X_interact[i] = mean( (X[ , i] .* log_realbmX_scl) / sd_X_interact[i] ) ;
    c_X_int[ , i] =  ( (X[ , i] .* log_realbmX_scl) / sd_X_interact[i] ) - means_X_interact[i] ;
  };


  // append to design matrix 
  // append means to means vector of predictor matrix 
  means_Xtmp = append_row(means_X, mean_log_realbm_scl);
  means_Xall = append_row(means_Xtmp, means_X_interact);
  
  
  // add centered latent variables to the predictor matrix
  c_Xtmp = append_col(c_X, c_log_realbmX_scl);
  c_Xall = append_col(c_Xtmp, c_X_int);

 
}

model {

  // initialize linear predictor term
  vector[Nobs] mu;


  // initialize linear predictor term
  mu = c_Int + c_Xall * b;

  for (n in 1:Nobs) {
    // add more terms to the linear predictor
    mu[n] += r_plot[plot[n]] + r_species[species[n]] ;
  }

  // priors including all constants
  target += normal_lpdf(b | 0, 5);
  target += student_t_lpdf(c_Int | 3, 0, 10);
  target += student_t_lpdf(sd_plot | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_plot[1] | 0, 1);
  target += student_t_lpdf(sd_species | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_species[1] | 0, 1);
  // likelihood including all constants
  target += bernoulli_logit_lpmf(pres | mu);
    
      
  // model the 'real' productivity 
  // iterate over the Nlat biomass observations
    for (i in 1:Nlatbm) {
      latbm[i] ~ normal(realbm[ plotbm[i] ], sigma_bm);
      };

  // informative priors for latent productivity intercept, vector of grassland arithmetic means and SDs
  target += normal_lpdf(realbm | meansbm, sdbm);
  // prior for latent sigmas
  target += cauchy_lpdf(sigma_bm | 0, 150);

}

generated quantities {
  // actual population-level intercept for uncentered variables
  real Intercept = c_Int - dot_product(means_Xall, b);
  
  // get log-likelihood
  vector[Nobs] log_lik;

  // generate y_pred values
  int post[Nobs] ;  // response variable
  

  // get log-likelihood
  for (n in 1:Nobs) {
    log_lik[n] = bernoulli_logit_lpmf(pres[n] | c_Int + c_Xall[n, ] * b + r_plot[plot[n]] + r_species[species[n]]);
  }

  // generate ypred
  for (n in 1:Nobs) {
    post[n] = bernoulli_logit_rng(c_Int + c_Xall[n, ] * b + r_plot[plot[n]] + r_species[species[n]]);
    };


  
}
