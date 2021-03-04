

data {
  // Number of observations for responses, equals the number of plot IDs for latent variables
  int<lower=1> Nplot;


  // delta15 response
  vector[Nplot] delta15;
  int<lower=1> NXdelta15; // number of design matrix columns
  matrix[Nplot, NXdelta15] Xdelta15; // design matrix


  // delta19 response
  vector[Nplot] delta19;
  int<lower=1> NXdelta19; // number of design matrix columns
  matrix[Nplot, NXdelta19] Xdelta19; // design matrix


  // productivity response - design matrix is only one variable
  int<lower=1> NXbm; // number of design matrix columns
  matrix[Nplot, NXbm] Xbm; // design matrix


  // latent productivity
  int<lower=1> Nlatbm; // number of biomass observations
  int<lower=1> plotbm[Nlatbm]; // latent variables ID variable
  vector[Nlatbm]   latbm; // measured biomass
  vector[Nplot] meansbm;  // biomass means as informative priors
  vector[Nplot] sdbm;  // biomass SDs as informative priors


}


transformed data {

  // delta 15 response
  // center the design matrix variables to speed up sampling
  matrix[Nplot, NXdelta15] c_Xdelta15; // define centered design matrix 
  // create vector of design matrix means for response delta15
  vector[NXdelta15] meansdelta15tmp;



  // delta 19 response
  // center the design matrix variables to speed up sampling
  matrix[Nplot, NXdelta19] c_Xdelta19; // define centered design matrix 
  // create vector of design matrix means for response delta19
  vector[NXdelta19] meansdelta19tmp;



  // productivity response
  // center the design matrix variables to speed up sampling
  matrix[Nplot, NXbm] c_Xbm; // define centered design matrix 
  // create vector of design matrix means for productivity response
  vector[NXbm] meansrealbm;



  // delta 15 response
  //get means to center design  matrix to speed up sampling
  for (i in 1:NXdelta15) {
    meansdelta15tmp[i] = mean(Xdelta15[, i]);
    c_Xdelta15[, i] = Xdelta15[, i] - meansdelta15tmp[i];
    };



  // delta 19 response
  //get means to center design  matrix to speed up sampling
  for (i in 1:NXdelta19) {
    meansdelta19tmp[i] = mean(Xdelta19[, i]);
    c_Xdelta19[, i] = Xdelta19[, i] - meansdelta19tmp[i];
    };


  // productivity response
  //get means to center design  matrix to speed up sampling
  for (i in 1:NXbm) {
    meansrealbm[i] = mean(Xbm[, i]);
    c_Xbm[, i] = Xbm[, i] - meansrealbm[i];
    };


}


parameters {

  // delta15 response
  real c_Int_delta15; // define temporary intercept with centered predictor matrix 
  vector[NXdelta15 + 1] b_delta15; // define betas for predictor matrix with appended latent variables
  real<lower=0> sigma_delta15; // define sigmas of delta15 response variable


  // delta19 response
  real c_Int_delta19; // define temporary intercept with centered predictor matrix 
  vector[NXdelta19 + 1] b_delta19; // define betas for predictor matrix with appended latent variables
  real<lower=0> sigma_delta19; // define sigmas of delta15 response variable


  // productivity response
  real c_Int_bm;  // temporary intercept for centered predictors
  vector[NXbm] b_bm;  // population-level effects
  real<lower=0> phi_bm;  // variance parameter


  // latent productivity
  vector<lower=0>[Nplot] realbm; // define the 'real' latent productivity, constrained to be positive
  real<lower=0> sigma_bm; // define sigmas of latents to be non-negative



}

transformed parameters {

  // set matrix and means to merge observed and latent predictor variable for delta15 response
  // transform latent to log and log-centered variable

  vector[Nplot] log_realbm = log(realbm);
  real sd_log_realbm = sd(log_realbm);
  
  real mean_log_realbm = mean( log_realbm / sd_log_realbm );
  vector[Nplot] c_log_realbm = ( log_realbm / sd_log_realbm ) - mean_log_realbm;



  // delta15 response
  // append means to means vector of predictor matrix 
  vector[NXdelta15 + 1] meansdelta15 = append_row(meansdelta15tmp, mean_log_realbm);
  // add centered latent variables to the predictor matrix
  matrix[Nplot, NXdelta15 + 1] c_Xdelta15all = append_col(c_Xdelta15, c_log_realbm);



  // delta19 response
  // append means to means vector of predictor matrix 
  vector[NXdelta19 + 1] meansdelta19 = append_row(meansdelta19tmp, mean_log_realbm);
    // add centered latent variables to the predictor matrix
  matrix[Nplot, NXdelta19 + 1] c_Xdelta19all = append_col(c_Xdelta19, c_log_realbm);




  // transformations to estimate shape and rate parameters in Gamma distribution models
  // see https://datascienceplus.com/bayesian-regression-with-stan-beyond-normality/ for parametrization
  // initialize linear predictor term and shape and rate parameters
 
  // productivity response
  vector[Nplot] mu_bm = exp(c_Int_bm + c_Xbm * b_bm); // inverse link is exp
  vector[Nplot] alpha_bm = mu_bm .* mu_bm / phi_bm;   // shape parameter, mu²*phi
  vector[Nplot] beta_bm = mu_bm / phi_bm;             // rate parameter


  
}



model {


  // priors for delta15 response model
  target += student_t_lpdf(c_Int_delta15 | 3, 0, 10);
  target += normal_lpdf(b_delta15 | 0, 10);
  target += cauchy_lpdf(sigma_delta15 | 0, 10);



  // priors for delta19 response model
  target += student_t_lpdf(c_Int_delta19 | 3, 0, 10);
  target += normal_lpdf(b_delta19 | 0, 10);
  target += cauchy_lpdf(sigma_delta19 | 0, 10);



  // priors for real productivity response model
  target += student_t_lpdf(c_Int_bm | 3, 0, 10);
  target += normal_lpdf(b_bm | 0, 2);
  target += gamma_lpdf( (1 / phi_bm ) | 0.01, 0.01);



  // delta15 response variable follows a normal distribution
  target += normal_id_glm_lpdf(delta15 | c_Xdelta15all, c_Int_delta15, b_delta15, sigma_delta15);


  // delta19 response variable follows a normal distribution
  target += normal_id_glm_lpdf(delta19 | c_Xdelta19all, c_Int_delta19, b_delta19, sigma_delta19);


  // productivity response variable follows a Gamma distribution
  target += gamma_lpdf(realbm | alpha_bm, beta_bm);



// model the 'real' latent productivity 
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

  // population-level intercepts for uncentered design matrix
  real Int_delta15 = c_Int_delta15 - dot_product(meansdelta15, b_delta15)    ;
  real Int_delta19 = c_Int_delta19 - dot_product(meansdelta19, b_delta19)    ;
  real Int_bm      = c_Int_bm      - dot_product(meansrealbm,  b_bm)         ;
  
  
  // fitted means and initializing posterior predictions for delta15
  vector[Nplot] mu_delta15 = c_Int_delta15 + c_Xdelta15all * b_delta15 ;
  vector[Nplot] post_delta15;


  // fitted means and initializing posterior predictions for delta19
  vector[Nplot] mu_delta19 = c_Int_delta19 + c_Xdelta19all * b_delta19 ;
  vector[Nplot] post_delta19;


  // posterior predictions for real productivity. yhat values equal mu_bm
  vector[Nplot] post_bm;

  
  // get log-likelihood for delta 19 response
  vector[Nplot] log_lik;



  // predict delta 19 responses with residual sigma
  for (i in 1:Nplot) {
    post_delta15[i] = normal_rng(mu_delta15[i], sigma_delta15);
    };


  // predict delta 19 responses with residual sigma
  for (i in 1:Nplot) {
    post_delta19[i] = normal_rng(mu_delta19[i], sigma_delta19);
    };


  // predict real productivity responses with shape and rate parameters
  for (i in 1:Nplot) {
    post_bm[i] = gamma_rng(alpha_bm[i], beta_bm[i]);
    };


  // get log-likelihood for delta 19 response
  for (n in 1:Nplot) {
    log_lik[n] = normal_lpdf(delta19[n] | mu_delta19[n], sigma_delta19);
  }



}

