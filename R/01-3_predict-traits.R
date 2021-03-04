


# predict establishment probabilities for trait values (mean +- SD)
# along the productivity gradient

# at the end of the script, I also calculate a Bayesian R². This is 
# done here in the script, although I could also have added a 'y_pred'
# vector in the .stan model file. But the Stan model objects were already 
# quite large and consumed large parts of the computer RAM, so 
# I had to take this rather lengthy option


predict_traits <- function(model = model,
         iter = iter, chains = chains,
         X = X, newly.sown.data = newly.sown.data) {
  
  ### I will make predictions for length = 200 
  require(rstan)
  require(brms)
  
  # provide latent productivity newdata to predict for
  realbm <- as.data.frame(model, pars = "realbm") # get latent productivity estimates
  median.realbmX <- summarise_draws(realbm, "median")  #get median of observations (N = 73) to define the range of newdata
  
  # create length = 200 newdata vector. Predictions are going to be plotted on the log scale.
  newdat.bm.seq <- seq(from = min(log(median.realbmX$median)), 
                       to = max(log(median.realbmX$median)), 
                       length.out = 200)
  newdat.orig.bm <- exp(newdat.bm.seq) # bring back to original scale to add to prediction data frame for plotting
  
  # repeat this vector n = iterations*chains times to scale the log_productivity by the log_productivity SD (stored in the model)
  newdat.bm <- matrix(rep(newdat.bm.seq, iter/2 * chains), 
                      ncol = iter/2 * chains, byrow = FALSE)
  sd.log.realbmX <- as.vector(extract(model, pars = "sd_log_realbmX", 
                                      permuted = FALSE)) # get log.realbm standard deviation
  newdat.bm <- t( t(newdat.bm) / sd.log.realbmX)
  
  
  
  ### Now, get matrices for all parameters. It's going to be ugly.
  
  
  # get Intercept for uncentered predictors
  traits.Int <- as.matrix(extract(model, pars = "Intercept", permuted = FALSE))
  # this call here and later repeats the posterior 200 times to meet the length of the newdata vector
  traits.Int <- matrix(rep(traits.Int, 200), nrow = 200, byrow = TRUE) 
  
  
  
  # get the betas
  traits.b <- as.data.frame(model, pars = "b")
  # append variable names from design matrix and add productivity and interaction names
  names(traits.b) <- c(colnames(X), "realbm", paste(colnames(X)[1:3], ".realbm", sep = "")) 
  
  
  # get the single-parameter betas. 
  # for traits, these refer to predictors scaled to standard deviation
  
  traits.seeddens <- traits.b["live.seeddensity.scl"]
  traits.seeddens <- as.matrix(traits.seeddens)
  traits.seeddens <- matrix(rep(traits.seeddens, 200), nrow = 200, byrow = TRUE)

  traits.logheight <- traits.b["logheight.scl"]
  traits.logheight <- as.matrix(traits.logheight)
  traits.logheight <- matrix(rep(traits.logheight, 200), nrow = 200, byrow = TRUE)
  
  traits.logseedMass <- traits.b["logseedMass.scl"]
  traits.logseedMass <- as.matrix(traits.logseedMass)
  traits.logseedMass <- matrix(rep(traits.logseedMass, 200), nrow = 200, byrow = TRUE)
  
  traits.SLA <- traits.b["SLA.scl"]
  traits.SLA <- as.matrix(traits.SLA)
  traits.SLA <- matrix(rep(traits.SLA, 200), nrow = 200, byrow = TRUE)
  
  traits.bm <- traits.b["realbm"]
  traits.bm <- as.matrix(traits.bm)
  traits.bm <- matrix(rep(traits.bm, 200), nrow = 200, byrow = TRUE)
  
  
  # get standard deviation of interactions, to bring interaction terms to the right scale
  traits.interactions_sd <- as.data.frame(model, pars = "sd_X_interact")
  # append variable names from design matrix and add productivity and interaction names
  names(traits.interactions_sd) <- c(paste(colnames(X)[1:3], ".realbm", sep = "")) 
  
  
  # get interaction sd for each interaction term
  traits.logheight.realbm.sd   <- as.matrix(traits.interactions_sd["logheight.scl.realbm"])
  traits.logseedMass.realbm.sd <- as.matrix(traits.interactions_sd["logseedMass.scl.realbm"])
  traits.SLA.realbm.sd         <- as.matrix(traits.interactions_sd["SLA.scl.realbm"])
  
  
  
  # get the interaction-parameter betas. these refer to interactions of predictors 
  # that were scaled to unit standard deviation, 
  # and interactions were scaled to interaction sd as well. 
  # This is a workaround to make predictions work:
  # instead of scaling the trait*realbiomass parameters to sd, 
  # I scale the beta of the interaction to sd
  
  traits.logheight.realbm <- traits.b["logheight.scl.realbm"]
  traits.logheight.realbm <- as.matrix(traits.logheight.realbm)
  traits.logheight.realbm <- traits.logheight.realbm / traits.logheight.realbm.sd
  traits.logheight.realbm <- matrix(rep(traits.logheight.realbm, 200), nrow = 200, byrow = TRUE)
  
  traits.logseedMass.realbm <- traits.b["logseedMass.scl.realbm"]
  traits.logseedMass.realbm <- as.matrix(traits.logseedMass.realbm)
  traits.logseedMass.realbm <- traits.logseedMass.realbm / traits.logseedMass.realbm.sd
  traits.logseedMass.realbm <- matrix(rep(traits.logseedMass.realbm, 200), nrow = 200, byrow = TRUE)
  
  traits.SLA.realbm <- traits.b["SLA.scl.realbm"]
  traits.SLA.realbm <- as.matrix(traits.SLA.realbm)
  traits.SLA.realbm <- traits.SLA.realbm / traits.SLA.realbm.sd
  traits.SLA.realbm <- matrix(rep(traits.SLA.realbm, 200), nrow = 200, byrow = TRUE)
  
  
  
  
  
  ### predict for each low, mean and high trait values along the latent productivity gradient
  ### keep trait-environment interactions in mind
  
  logheight.high     <- mean(X$logheight.scl)   + sd(X$logheight.scl)  
  logseedMass.high   <- mean(X$logseedMass.scl) + sd(X$logseedMass)
  SLA.high           <- mean(X$SLA.scl)         + sd(X$SLA.scl)        
  
  logheight.low      <- mean(X$logheight.scl)   - sd(X$logheight.scl)  
  logseedMass.low    <- mean(X$logseedMass.scl) - sd(X$logseedMass.scl)
  SLA.low            <- mean(X$SLA.scl)         - sd(X$SLA.scl)        
  
  logheight.mean     <- mean(X$logheight.scl)   
  logseedMass.mean   <- mean(X$logseedMass.scl) 
  SLA.mean           <- mean(X$SLA.scl)         
  
  
  
  
  ### Now calculate the posterior predictions for trait values along the productivity gradient, including interactions 
  ### directly apply the inverse logit link function from the brms package, scaling parameters default to range [0,1]: binomial
  ### Here, I calculate the posterior intervals directly: Median, 5% and 95% interval
  
  # predict for mean trait values and condition on mean live seeding density, 
  # i.e. live seeding density = 0 (hence not included below), as the predictor is scaled to zero mean)

  ### Log. height
  
  pred.logheight.high <- brms::inv_logit_scaled(
    traits.Int + 

      traits.logheight   * logheight.high   +
      traits.logseedMass * logseedMass.mean +
      traits.SLA         * SLA.mean         +
      
      traits.bm * newdat.bm + 
      
      traits.logheight.realbm   * (logheight.high   * newdat.bm) +
      traits.logseedMass.realbm * (logseedMass.mean * newdat.bm) +
      traits.SLA.realbm         * (SLA.mean         * newdat.bm) 
  )
  
  pred.logheight.high <- t(apply(pred.logheight.high, MARGIN = 1, function(x) quantile(x, probs = c(0.05, 0.50, 0.95)) ))
  pred.logheight.high <- data.frame(trait = "(a) Height (log.)", value = rep("Mean + SD", 200), pred.logheight.high)
  
  
  
  pred.logheight.low <- brms::inv_logit_scaled(
    traits.Int + 
      
      traits.logheight   * logheight.low   +
      traits.logseedMass * logseedMass.mean +
      traits.SLA         * SLA.mean         +
      
      traits.bm * newdat.bm + 
      
      traits.logheight.realbm   * (logheight.low    * newdat.bm) +
      traits.logseedMass.realbm * (logseedMass.mean * newdat.bm) +
      traits.SLA.realbm         * (SLA.mean         * newdat.bm)
  )
  
  pred.logheight.low <- t(apply(pred.logheight.low, MARGIN = 1, function(x) quantile(x, probs = c(0.05, 0.50, 0.95)) ))
  pred.logheight.low <- data.frame(trait = "(a) Height (log.)", value = rep("Mean - SD", 200), pred.logheight.low)
  
  
  # rbind posterior predictions for log. height values into one data frame
  pred.logheight <- rbind(pred.logheight.high, pred.logheight.low)
  
  
  
  
  
  ### Log. seed mass
  
  pred.logseedMass.high <- brms::inv_logit_scaled(
    traits.Int + 
      
      traits.logheight   * logheight.mean   +
      traits.logseedMass * logseedMass.high +
      traits.SLA         * SLA.mean         +
      
      traits.bm * newdat.bm                 + 
      
      traits.logheight.realbm   * (logheight.mean   * newdat.bm)  +
      traits.logseedMass.realbm * (logseedMass.high * newdat.bm)  +
      traits.SLA.realbm         * (SLA.mean         * newdat.bm)
  )
  
  pred.logseedMass.high <- t(apply(pred.logseedMass.high, MARGIN = 1, function(x) quantile(x, probs = c(0.05, 0.50, 0.95)) ))
  pred.logseedMass.high <- data.frame(trait = "(b) Seed mass (log.)", value = rep("Mean + SD", 200), pred.logseedMass.high)
  
  
  pred.logseedMass.low <- brms::inv_logit_scaled(
    traits.Int + 
      
      traits.logheight   * logheight.mean   +
      traits.logseedMass * logseedMass.low  +
      traits.SLA         * SLA.mean         +
      
      traits.bm * newdat.bm                 + 
      
      traits.logheight.realbm   * (logheight.mean   * newdat.bm)  +
      traits.logseedMass.realbm * (logseedMass.low  * newdat.bm)  +
      traits.SLA.realbm         * (SLA.mean         * newdat.bm) 
  )
  pred.logseedMass.low <- t(apply(pred.logseedMass.low, MARGIN = 1, function(x) quantile(x, probs = c(0.05, 0.50, 0.95)) ))
  pred.logseedMass.low <- data.frame(trait = "(b) Seed mass (log.)", value = rep("Mean - SD", 200), pred.logseedMass.low)
  
  # rbind into one data frame
  pred.logseedMass <- rbind(pred.logseedMass.high, pred.logseedMass.low)
  
  
  
  
  
  ### Specific leaf area SLA
  
  pred.SLA.high <- brms::inv_logit_scaled(
    traits.Int + 
      
      traits.logheight   * logheight.mean   +
      traits.logseedMass * logseedMass.mean +
      traits.SLA         * SLA.high         +
      
      traits.bm * newdat.bm                 + 
      
      traits.logheight.realbm   * (logheight.mean   * newdat.bm)  +
      traits.logseedMass.realbm * (logseedMass.mean * newdat.bm)  +
      traits.SLA.realbm         * (SLA.high         * newdat.bm)
  )
  
  pred.SLA.high <- t(apply(pred.SLA.high, MARGIN = 1, function(x) quantile(x, probs = c(0.05, 0.50, 0.95)) ))
  pred.SLA.high <- data.frame(trait = "(c) Specific leaf area", value = rep("Mean + SD", 200), pred.SLA.high)
  
  
  
  pred.SLA.low <- brms::inv_logit_scaled(
    traits.Int + 
      
      traits.logheight   * logheight.mean   +
      traits.logseedMass * logseedMass.mean +
      traits.SLA         * SLA.low          +
      
      traits.bm * newdat.bm                 + 
      
      traits.logheight.realbm   * (logheight.mean   * newdat.bm)  +
      traits.logseedMass.realbm * (logseedMass.mean * newdat.bm)  +
      traits.SLA.realbm         * (SLA.low          * newdat.bm)
  )
  
  pred.SLA.low <- t(apply(pred.SLA.low, MARGIN = 1, function(x) quantile(x, probs = c(0.05, 0.50, 0.95)) ))
  pred.SLA.low <- data.frame(trait = "(c) Specific leaf area", value = rep("Mean - SD", 200), pred.SLA.low)
  
  # rbind into one data frame
  pred.SLA <- rbind(pred.SLA.high, pred.SLA.low)
  
  
  
  
  ### and now combine all posterior predictions into one data frame
  pred.traits <- rbind( pred.logheight, 
                        pred.logseedMass,
                        pred.SLA)

  pred.traits$trait <- as.factor(pred.traits$trait)
  pred.traits$value <- as.factor(pred.traits$value)
  
  # add original-scale productivity values to the data frame
  # repeat 6 times (each high and low SD for three traits)
  pred.traits$productivity <- rep(newdat.orig.bm, 6) 
  
  
  
  
  
  
  
  ################################################################################
  ### 5 calculate fitted y_hats of observations for Bayes R2            ##########
  ################################################################################
  
  
  # get the posteriors of plot random intercepts
  r.plot <- t(as.matrix(model, pars = "r_plot"))

  # get the posteriors of species random intercepts
  r.species <- t(as.matrix(model, pars = "r_species"))

  # posteriors of latent productivity (for N = 73 plots) already stored as 'realbm'. Log as in the model
  t.realbm <- t(realbm)
  
  
  log.realbmX <- list()
  plotID <- as.integer(newly.sown.data$plot)
  
  for (i in 1:nrow(newly.sown.data)) {
    log.realbmX[[i]] <- log(t.realbm[ plotID[i] , ] ) ;
  };
  
  log.realbmX <- do.call(rbind, log.realbmX)
  
  
  
  
  
  
  ### start with the conditional R2 including random intercepts --------
  
  # create empty list to loop over
  pred.obs.cond <- list()
  
  
  # all parameters were repeated n = 200 rows for the previous predictions, so take only one row of those for predictions
  for(i in 1:nrow(newly.sown.data)){
    # predict  
    pred.obs.i  <- brms::inv_logit_scaled(
        traits.Int[1,] + 
          
        traits.seeddens[1,]    * newly.sown.data$live.seeddensity.scl[i]  +  
          
        traits.logheight[1,]   * newly.sown.data$logheight.scl[i]   +
        traits.logseedMass[1,] * newly.sown.data$logseedMass.scl[i] +
        traits.SLA[1,]         * newly.sown.data$SLA.scl[i]         +

        traits.bm[1,] * log.realbmX[i,]                              + 
        
        traits.logheight.realbm[1,]   * (newly.sown.data$logheight.scl[i]   * log.realbmX[i,] )  +
        traits.logseedMass.realbm[1,] * (newly.sown.data$logseedMass.scl[i] * log.realbmX[i,] )  +
        traits.SLA.realbm[1,]         * (newly.sown.data$SLA.scl[i]         * log.realbmX[i,] )  + 
        
        r.species[as.integer(newly.sown.data$species[i])] + 
        r.plot[as.integer(newly.sown.data$plot[i])]  
    )
    
    pred.obs.cond[[i]] <- pred.obs.i
    
  }
  
  pred.obs.cond <- t(do.call(rbind, pred.obs.cond)) # transpose to meet format for rstantools::bayes_R2
  
  R2.cond <- rstantools::bayes_R2(object = pred.obs.cond, 
                                            y = newly.sown.data$presence)

  
  
  ### now the same for the marginal R2 without random intercepts --------
  
  # create empty list to loop over
  pred.obs.marg <- list()
  
  for(i in 1:nrow(newly.sown.data)){
    # predict  
    pred.obs.i  <- brms::inv_logit_scaled(
        traits.Int[1,] + 
        
        traits.seeddens[1,]    * newly.sown.data$live.seeddensity.scl[i]  +  
          
        traits.logheight[1,]   * newly.sown.data$logheight.scl[i]   +
        traits.logseedMass[1,] * newly.sown.data$logseedMass.scl[i] +
        traits.SLA[1,]         * newly.sown.data$SLA.scl[i]         +

        traits.bm[1,] * log.realbmX[i,]                              + 
        
        traits.logheight.realbm[1,]   * (newly.sown.data$logheight.scl[i]   * log.realbmX[i,] )  +
        traits.logseedMass.realbm[1,] * (newly.sown.data$logseedMass.scl[i] * log.realbmX[i,] )  +
        traits.SLA.realbm[1,]         * (newly.sown.data$SLA.scl[i]         * log.realbmX[i,] )  
      
    )
    
    pred.obs.marg[[i]] <- pred.obs.i
    
  }
  
  pred.obs.marg <- t(do.call(rbind, pred.obs.marg)) # transpose to meet format for rstantools::bayes_R2
  
  R2.marg <- rstantools::bayes_R2(object = pred.obs.marg, 
                                            y = newly.sown.data$presence)
  
  
  results <- list(pred.traits = pred.traits, 
                  R2.cond = R2.cond, 
                  R2.marg = R2.marg)
  
  return(results)
  
}
