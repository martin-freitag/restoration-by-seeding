

predict_established_year5 <- function(object = object, 
                                      response.year5 = response.year5, 
                                      iter = iter, 
                                      chains = chains) {

### I will make predictions for length=200 

# provide latent productivity newdata to predict for
median.realbm <- get_ci(object, pars = "realbm", probs = 0.5)[ , 1] # get latent productivity estimates

# create length=200 newdata vector. Predictions are going to be plotted on the log scale.
newdat.bm.seq <- seq(from = min(log(median.realbm)), to = max(log(median.realbm)), length.out = 200)
newdat.orig.bm <- exp(newdat.bm.seq) # bring back to original scale to add to prediction data frame for plotting

# repeat this vector n=25000 iterations * 6 chains times to scale the log_productivity by the log_productivity SD (stored in the model)
newdat.bm <- matrix(rep(newdat.bm.seq, iter/2 * chains), ncol = iter/2 * chains, byrow=FALSE)
sd.log.realbm <- as.vector(extract(object, pars = "sd_log_realbm", permuted = TRUE)[[1]] ) # get log_realbm standard deviation
newdat.bm <- t( t(newdat.bm) / sd.log.realbm)



### Now, get matrices for all parameters. It's going to be ugly.


# get Intercept for centered predictors
Int.delta19 <- as.matrix(extract(object, pars = "Int_delta19", permuted = TRUE)[[1]] )
# this call here and later repeats the posterior 200 times to meet the length of the newdata vector
Int.delta19 <- matrix(rep(Int.delta19, 200), nrow = 200, byrow = TRUE) 



# get the beta estimate for productivity
b.delta19 <- as.data.frame(extract(object, pars = "b_delta19", permuted = TRUE)[[1]])
# append variable names from design matrix and add productivity and interaction names
names(b.delta19) <- names.b.delta19


# get the single-parameter betas. these refer to predictors scaled to standard deviation

delta19.HAI <- b.delta19[ , "delta19_HAI"]
delta19.HAI <- as.matrix(delta19.HAI)
delta19.HAI <- matrix(rep(delta19.HAI, 200), nrow = 200, byrow = TRUE)


delta19.SCH <- b.delta19[ , "delta19_SCH"]
delta19.SCH <- as.matrix( delta19.SCH)
delta19.SCH <- matrix(rep(delta19.SCH, 200), nrow = 200, byrow = TRUE)


delta19.logrealbm <- b.delta19["delta19_logrealbm"]
delta19.logrealbm <- as.matrix( delta19.logrealbm)
delta19.logrealbm <- matrix(rep(delta19.logrealbm, 200), nrow = 200, byrow = TRUE)




### Now calculate the posterior predictions for the three regions along the productivity gradient
#### Here, I calculate the posterior intervals directly: Median, 5% and 95% interval


# because I use the intercept for centred predictors, the predictions
# are already conditioned on mean values of mowing and grazing intensities
# and mean delta richness in the first year

# predict

pred.mean <- Int.delta19 + 0.5 * delta19.HAI   + 0.5 * delta19.SCH +
  newdat.bm * delta19.logrealbm

pred.mean     <- t(apply(pred.mean, MARGIN = 1, function(x) quantile(x, probs = c(0.05, 0.50, 0.95)) ))

### add productivity values used for prediction (original, not log-transformed scale)
pred.mean  <- data.frame(productivity = newdat.orig.bm, 
                         pred.mean)



# now predict for the three regions separately

pred.ALB.mean <- Int.delta19 + 0 * delta19.HAI   + 0 * delta19.SCH +
  newdat.bm * delta19.logrealbm

pred.HAI.mean <- Int.delta19 + 1 * delta19.HAI   + 0 * delta19.SCH +
  newdat.bm * delta19.logrealbm

pred.SCH.mean <- Int.delta19 + 0 * delta19.HAI   + 1 * delta19.SCH +
  newdat.bm * delta19.logrealbm


pred.regions.mean     <- rbind(pred.ALB.mean, pred.HAI.mean, pred.SCH.mean)
pred.regions.mean     <- t(apply(pred.regions.mean, MARGIN = 1, function(x) quantile(x, probs = c(0.05, 0.50, 0.95)) ))



pred.regions.mean  <- data.frame(productivity = rep(newdat.orig.bm, 3), 
                                 region = factor(c(rep("ALB", 200), rep("HAI", 200), rep("SCH", 200)) ),
                                 pred.regions.mean)



# To plot not only the predictions but also the observed data on plots (73 obs.),
# add the log latent productivity median to the response data in the pred.mean data table
# observed data is ordered by region: 25 obs from ALB, 23 obs from HAI, 25 obs from SCH

# median.realbm object contains the latent productivity medians for the 73 grassland plots. 

# and add observed number of established species in fifth year
obs.delta <- data.frame(productivity = median.realbm, 
                             delta   = response.year5, 
                             region  = factor(c(rep("ALB", 25), 
                                                rep("HAI", 23), 
                                                rep("SCH", 25) )  )
)


# because the regions have different ranges of productivity, 
# set predictions for productivity values outside of the range to NA

minrealbm.ALB <- min(obs.delta[ obs.delta$region == "ALB" , "productivity" ] )
maxrealbm.ALB <- max(obs.delta[ obs.delta$region == "ALB" , "productivity" ] )

minrealbm.HAI <- min(obs.delta[ obs.delta$region == "HAI" , "productivity" ] )
maxrealbm.HAI <- max(obs.delta[ obs.delta$region == "HAI" , "productivity" ] )

minrealbm.SCH <- min(obs.delta[ obs.delta$region == "SCH" , "productivity" ] )
maxrealbm.SCH <- max(obs.delta[ obs.delta$region == "SCH" , "productivity" ] )




pred.regions.mean[ which(pred.regions.mean$region == "ALB" &
                           (pred.regions.mean$productivity > maxrealbm.ALB |
                              pred.regions.mean$productivity < minrealbm.ALB 
                           )
) , 3:5 ]  <- NA

pred.regions.mean[ which(pred.regions.mean$region == "HAI" &
                           (pred.regions.mean$productivity > maxrealbm.HAI |
                              pred.regions.mean$productivity < minrealbm.HAI 
                           )
) , 3:5 ]  <- NA

pred.regions.mean[ which(pred.regions.mean$region == "SCH" &
                           (pred.regions.mean$productivity > maxrealbm.SCH |
                              pred.regions.mean$productivity < minrealbm.SCH 
                           ) 
), 3:5 ]  <- NA


out <- list(pred.mean = pred.mean, 
            obs.delta = obs.delta, 
            pred.regions.mean = pred.regions.mean)

return(out)

}
