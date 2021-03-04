



# load some packages. Add version numbers as reported in the manuscript

library(data.table)
library(rstan)
library(brms)
library(bayesplot)
library(posterior)
library(ggpubr)
library(ggplot2)

# functions for Bayes R2 and panel.cor correlation in pairs plots
source("R/helpers.R")





################################################################################
# 0.1 Load and prepare the data                                        #########
################################################################################




# we are going to use the 'metadata' data, which contains info about  
# the experimental design and metadata like grazing and mowing intensities, 
# but also productivity (aboveground biomass) sampling data
metadata <- read.table("data/metadata.txt", header = TRUE,
                          sep = "\t", dec = ".", stringsAsFactors = TRUE)



# here, we source a script to prepare the 'sown.long' data frame, 
# which contains species cover data of only sown species. 
# Species that were already present on control and 
# disturbance treatments were removed. The full 'species.long' 
# data frame is sourced within the script as well.

source("R/00-1_sown-and-established-species-data-wrangling.R") 



# also load the species cover data from our vegetation surveys, long format
species.long <- read.table("data/species-long.txt", header = TRUE,
                           sep = "\t", dec = ".", stringsAsFactors = TRUE)





# Convert long format species data to wide format to calculate species richness


# first for all species, afterwards only sown species
# convert species data to wide format to calculate diversity indices
species.wide <- tidyr::pivot_wider(species.long, 
                                   names_from = species, 
                                   values_from = cover,
                                   values_fill = 0)

sown.wide <- tidyr::pivot_wider(sown.long, 
                                names_from = species, 
                                values_from = cover, 
                                values_fill = 0) 
# Note that NA's get replaced by 0. Necessary for species sown 
# in one region but not in the other region



# first make sure data frames in the same order
identical(metadata$plotIdent, species.wide$plotIdent, sown.wide$plotIdent)



# Three grasslands were not accessible in 2018 (3 grasslands * 4 treatments) 
# This gives 73 grasslands * 4 treatments * 5 years = 1,460 obs. - 12 missing obs.

nrow(metadata) #1,448 observations
nrow(species.wide)
nrow(sown.wide)
# same number of observations


# it is safe now to remove the metadata columns from the species data frame
species.wide <- species.wide[ , -1] # only 'plotIdent' identifier
sown.wide <- sown.wide[ , -c(1:5)]  # some more metadate on year, treatment etc.

any(is.na(species.wide)) # no NA's in the data
any(is.na(sown.wide))




################################################################################
### 0.2 Calculate diversity metrics                                      #######
################################################################################


# we'll go with overall species richness and richness of sown species


# species richness
metadata$richness <- apply(species.wide, MARGIN = 1, function(x) sum(x > 0)) 

metadata$sown.richness <- apply(sown.wide, MARGIN = 1, function(x) sum(x > 0)) 



# calculate differences ('Delta') in richness

setDT(metadata) #set to data.table format to use the diff-like function
metadata$plot.year <- interaction(metadata$plot, metadata$year) # identifier


# calculate differences between treatments within plot/year with control as reference
metadata <- metadata[ , ":=" (delta.rich =  richness - richness[treatment == "control"]  ,
                                    delta.sown =  sown.richness - sown.richness[treatment == "control"] 
                                    ) ,
                                    by = plot.year ] # differences within plot and year



###   # subset the data here to test whether the three ALB grasslands with
###   # very low productivity have a large influence on parameter estimates
###
###   metadata <- droplevels(metadata[metadata$plot != "AEG07" &
###                                     metadata$plot != "AEG09" &
###                                     metadata$plot != "AEG26"
###                                   , ])






################################################################################
### 1 set up the design matrix and data                                  #######
################################################################################



# LUI variable means are the same over all years (and treatments of course)

# We will sample four multivariate models: 
#
# - seeding only treatment
#     + difference in      species richness (2015 and 2019)
#     + difference in sown species richness (2015 and 2019)
#
# - seeding and disturbance treatment
#     + difference in      species richness (2015 and 2019)
#     + difference in sown species richness (2015 and 2019)



# first prepare the responses

# seeding only responses
delta.rich.seed.only15      = subset(metadata, treatment == "seed" 
                                & year == 2015)$delta.rich # response in first year
delta.rich.seed.only19      = subset(metadata, treatment == "seed" 
                                & year == 2019)$delta.rich # response in fifth year
delta.sown.seed.only15      = subset(metadata, treatment == "seed" 
                                & year == 2015)$delta.sown # response in first year
delta.sown.seed.only19      = subset(metadata, treatment == "seed" 
                                & year == 2019)$delta.sown # response in fifth year



# seeding and disturbance responses
delta.rich.seed.disturb15      = subset(metadata, treatment == "seeddisturb" 
                                      & year == 2015)$delta.rich # response in first year
delta.rich.seed.disturb19      = subset(metadata, treatment == "seeddisturb" 
                                      & year == 2019)$delta.rich # response in fifth year
delta.sown.seed.disturb15      = subset(metadata, treatment == "seeddisturb" 
                                      & year == 2015)$delta.sown # response in first year
delta.sown.seed.disturb19      = subset(metadata, treatment == "seeddisturb" 
                                      & year == 2019)$delta.sown # response in fifth year




# The predictor matrix for the 2015 data is the same for all models, 
# it is mainly the response and the 2019 predictor matrix which differ
# between models. Region and mowing and grazing intensities 
# do not differ between years or treatments, subset here to control in 2015 
# to get the right format


Xdelta15 <- with(subset(metadata, treatment == "control" & year == 2015),
                    cbind(
                      HAI      = ifelse(region == "HAI", 1, 0),
                      SCH      = ifelse(region == "SCH", 1, 0),
                      Mow      = Mow / sd(Mow),
                      logGraz  = log(Graz + 1) / sd(log(Graz + 1) )
                    ) )

# predictor matrix for 2019 response will be defined later on



# now the predictors for productivity, the same for all models

Xbm <- with(subset(metadata, treatment == "control" & year == 2015),
                 cbind(bm_HAI = ifelse(region == "HAI", 1, 0),
                       bm_SCH = ifelse(region == "SCH", 1, 0),
                       bm_logFert  = log(Fert + 1) / sd(log(Fert + 1) )
                 ) )



# Productivity will be modelled as latent using the biomass 
# measurements from 2015 to 2019 from the control treatment. These 
# five measurements will be modelled as emerging from a 'true' 
# productivity with mean mu and measurement error sigma
# (to avoid shrinkage towards the global mean, we specify priors 
# for the mean for each grassland separately)


latbiomass  <- droplevels(subset(metadata, 
                                 treatment == "control", 
                                 select = c("plot", "biomass.g") ) )

latbiomass <- latbiomass[complete.cases(latbiomass) , ] #remove NA's

# get means and SD per plot to use as informed priors for the productivity mean
means.bm <- aggregate(biomass.g ~ plot, mean, data = latbiomass)[,2]
sd.bm    <- aggregate(biomass.g ~ plot, sd, data = latbiomass)[,2]       






# set the input data for the Stan model starting with:

# seeding only treatment: difference in species richness

standata.delta.rich.seed.only <- list(
          Nplot     = nrow(Xdelta15),             # number of plots the same everywhere
          
          Xdelta15     = Xdelta15,              # design matrix for first year
          NXdelta15    = ncol(Xdelta15),          # number of design matrix columns

          Xdelta19     = cbind(Xdelta15, delta.rich.seed.only15),  # design matrix for fifth year
          NXdelta19    = ncol(Xdelta15) + 1,        # number of design matrix columns
          
          delta15      = delta.rich.seed.only15, # response in first year
          delta19      = delta.rich.seed.only19, # response in fifth year
          
          Xbm = Xbm,                     # design matrix for latent productivity model
          NXbm = ncol(Xbm),              # number of design matrix columns
          
          # each five biomass measurements per grassland, with measurement error
          Nlatbm  = nrow(latbiomass),            # number of biomass observations (73 grasslands * 5 years minus 3 NAs)
          plotbm  = as.integer(factor(latbiomass$plot)), # integer coding of plots
          latbm   = latbiomass$biomass.g,        # biomass measurements on scale g/m²
          meansbm = means.bm,   # grassland biomass means as informed prior
          sdbm    = sd.bm       # standard deviations as informed prior for the mean
)




# copy the data for the three other models, 
# then replace the response variables

standata.delta.sown.seed.only            <- standata.delta.rich.seed.only
standata.delta.rich.seed.disturb         <- standata.delta.rich.seed.only
standata.delta.sown.seed.disturb         <- standata.delta.rich.seed.only


# Seeding only treatment: difference in sown species richness
standata.delta.sown.seed.only$delta15  <- delta.sown.seed.only15
standata.delta.sown.seed.only$delta19  <- delta.sown.seed.only19
standata.delta.sown.seed.only$Xdelta19 <- cbind(Xdelta15, delta.sown.seed.only15)

# Seeding and disturbance treatment: difference in species richness
standata.delta.rich.seed.disturb$delta15  <- delta.rich.seed.disturb15
standata.delta.rich.seed.disturb$delta19  <- delta.rich.seed.disturb19
standata.delta.rich.seed.disturb$Xdelta19 <- cbind(Xdelta15, delta.rich.seed.disturb15)

# Seeding and disturbance treatment: difference in sown species richness
standata.delta.sown.seed.disturb$delta15  <- delta.sown.seed.disturb15
standata.delta.sown.seed.disturb$delta19  <- delta.sown.seed.disturb19
standata.delta.sown.seed.disturb$Xdelta19 <- cbind(Xdelta15, delta.sown.seed.disturb15)




################################################################################
### 2.1 fit the multivariate model with Stan                               #####
################################################################################


# stan sampling options
rstan_options(auto_write = TRUE) # optimization


# list the parameters to save to save RAM space
save_pars <- c("c_Int_delta15", "b_delta15", "sigma_delta15",   # delta 15, first year model
               "c_Int_delta19", "b_delta19", "sigma_delta19",   # delta 19, fifth year model
               "Int_delta19",                                   # save for plotting later
               "Int_bm",        "b_bm",      "phi_bm",          # productivity model, phi equals variance
               "realbm",                     "sigma_bm",        # latent productivity model
               
               "mu_delta15",    "post_delta15", #mu and posterior predictions delta15, first year model
               "mu_delta19",    "post_delta19", #mu and posterior predictions delta19, fifth year model
               "mu_bm",         "post_bm",      #mu and posterior predictions productivity
               
               "sd_log_realbm",                 #save sd of log productivity for plotting later (backtransform scaled productivity)
               "log_lik")

chains <- 4
iter   <- 10000


# sample the Stan models
stanmod.delta.rich.seed.only <- stan(file = "Stan_files/01-2_delta-multivariate.stan", 
                                     data = standata.delta.rich.seed.only,
                                     chains = chains, cores = parallel::detectCores(), 
                                     iter = iter, thin = 1, 
                                     control = list(adapt_delta = 0.8, max_treedepth = 12),
                                     pars = save_pars, include = TRUE,
                                     seed = 123) 

stanmod.delta.sown.seed.only <- stan(file = "Stan_files/01-2_delta-multivariate.stan", 
                                     data = standata.delta.sown.seed.only,
                                     chains = chains, cores = parallel::detectCores(), 
                                     iter = iter, thin = 1, 
                                     control = list(adapt_delta = 0.8, max_treedepth = 12),
                                     pars = save_pars, include = TRUE,
                                     seed = 123) 

stanmod.delta.rich.seed.disturb <- stan(file = "Stan_files/01-2_delta-multivariate.stan", 
                                        data = standata.delta.rich.seed.disturb,
                                        chains = chains, cores = parallel::detectCores(), 
                                        iter = iter, thin = 1, 
                                        control = list(adapt_delta = 0.8, max_treedepth = 12),
                                        pars = save_pars, include = TRUE,
                                        seed = 123) 

stanmod.delta.sown.seed.disturb <- stan(file = "Stan_files/01-2_delta-multivariate.stan", 
                                        data = standata.delta.sown.seed.disturb,
                                        chains = chains, cores = parallel::detectCores(), 
                                        iter = iter, thin = 1, 
                                        control = list(adapt_delta = 0.8, max_treedepth = 12),
                                        pars = save_pars, include = TRUE,
                                        seed = 123) 




# save the model posteriors

save(stanmod.delta.rich.seed.only,    stanmod.delta.sown.seed.only,
     stanmod.delta.rich.seed.disturb, stanmod.delta.sown.seed.disturb,
     file="Rdata/01-2_stan-posteriors-delta-multivariate.RData")





# parameters to extract for table
pars.model <- c("c_Int_delta15", "b_delta15", "sigma_delta15",
                "c_Int_delta19", "b_delta19", "sigma_delta19", 
                "Int_bm", "b_bm", "phi_bm", "sigma_bm"
)

# get names of design matrices for responses
names.b.delta15 <- c(paste("delta15", c(colnames(Xdelta15), "logrealbm" ), sep ="_") )
names.b.delta19 <- c(paste("delta19", c(colnames(Xdelta15), "delta15", "logrealbm" ), sep ="_") )
names.b.bm      <- colnames(Xbm)

# and get the names for model parameters
names.pars.model <- c("c_Int_delta15", names.b.delta15, "sigma_delta15",
                      "c_Int_delta19", names.b.delta19, "sigma_delta19", 
                      "Int_bm", names.b.bm, "phi_bm", "sigma_bm"
)




# extract the posterior summaries, together with n_eff and Rhat

ci.delta.rich.seed.only    <- get_ci(stanmod.delta.rich.seed.only,
                                     pars = pars.model, names = names.pars.model)

ci.delta.sown.seed.only    <- get_ci(stanmod.delta.sown.seed.only, 
                                     pars = pars.model, names = names.pars.model)

ci.delta.rich.seed.disturb <- get_ci(stanmod.delta.rich.seed.disturb, 
                                     pars = pars.model, names = names.pars.model)

ci.delta.sown.seed.disturb <- get_ci(stanmod.delta.sown.seed.disturb, 
                                     pars = pars.model, names = names.pars.model)

cis.delta.rich <- merge(ci.delta.rich.seed.only, 
                        ci.delta.rich.seed.disturb,
                        by = "parameter")

cis.delta.sown <- merge(ci.delta.sown.seed.only, 
                        ci.delta.sown.seed.disturb, 
                        by = "parameter", all = TRUE)


cis.delta.sown <- cis.delta.sown[ order(match(cis.delta.sown$parameter, names.pars.model)) , ]
cis.delta.rich <- cis.delta.rich[ order(match(cis.delta.rich$parameter, names.pars.model)) , ]

# make nice interval labels
names(cis.delta.sown)[-1] <- rep(c("5%", "50%", "95%", "N_eff"), 2)
names(cis.delta.rich)[-1] <- rep(c("5%", "50%", "95%", "N_eff"), 2)


# print parameter, then replace by nice labels
cis.delta.rich$parameter

# the nice labels
parameter.names <- c("Intercept",
                     "Region HAI", "Region SCH", "Mowing intensity", 
                     "Grazing intensity (log.)", "Productivity (log.)", 
                     "$\\sigma$",
                     
                     "Intercept",
                     "Region HAI", "Region SCH", "Mowing intensity", 
                     "Grazing intensity (log.)", "$\\Delta$ richness 1st year", 
                     "Productivity (log.)", 
                     "$\\sigma$",
                     
                     "Intercept",
                     "Region HAI", "Region SCH", "Fertilization intensity (log.)",
                     "$\\phi$", 
                     "Measurement error $\\sigma$" )

cis.delta.rich$parameter <- parameter.names
cis.delta.sown$parameter <- parameter.names



# merge the credible intervals with the Bayes R2 later





################################################################################
### 2.2 calculate a Bayes R2, following Gelman et al. 2019               #######
################################################################################


# there are two implementations of the Bayesian R squared:
# the 2018 preprint version simply calculated variance_fit / (variance_observed_data). i.e. draws from residuals
# the published 2019 version used variance_fit / (variance_fit + variance_modelled), i.e. draws from sigma 
#
# brms and rstantools libraries use the version with draws from residuals, 
# which gives narrower R2 with small N, but almost the same results with plenty of data
# Ref: https://avehtari.github.io/bayes_R2/bayes_R2.html
#
# For the paper I used the first option (function bayesR2_res). For 
# comparison I also provide the other option below.


# used in the paper:

# delta richness seeding only treatment 
R2.rich.seed.only.delta15 <- bayesR2_res(fit = stanmod.delta.rich.seed.only, ypred = "mu_delta15", 
                                         y = standata.delta.rich.seed.only$delta15 )

R2.rich.seed.only.delta19 <- bayesR2_res(fit = stanmod.delta.rich.seed.only, ypred = "mu_delta19", 
                                         y = standata.delta.rich.seed.only$delta19 )

R2.rich.seed.only.bm <- bayesR2_res(fit = stanmod.delta.rich.seed.only, ypred = "mu_bm", 
                                    y = get_ci(stanmod.delta.rich.seed.only, pars = "realbm", 
                                               probs = 0.5)[ , 1] )

R2.rich.seed.only <- rbind(summarise_draws(R2.rich.seed.only.delta15, ~quantile(.x, probs = c(0.05, 0.50, 0.95))),
                           summarise_draws(R2.rich.seed.only.delta19, ~quantile(.x, probs = c(0.05, 0.50, 0.95))),
                           summarise_draws(R2.rich.seed.only.bm, ~quantile(.x, probs = c(0.05, 0.50, 0.95)))
)




# delta sown richness seeding only treatment 
R2.sown.seed.only.delta15 <- bayesR2_res(fit = stanmod.delta.sown.seed.only, ypred = "mu_delta15", 
                                         y = standata.delta.sown.seed.only$delta15 )

R2.sown.seed.only.delta19 <- bayesR2_res(fit = stanmod.delta.sown.seed.only, ypred = "mu_delta19", 
                                         y = standata.delta.sown.seed.only$delta19 )

R2.sown.seed.only.bm <- bayesR2_res(fit = stanmod.delta.sown.seed.only, ypred = "mu_bm", 
                                    y = get_ci(stanmod.delta.sown.seed.only, pars = "realbm", 
                                               probs = 0.5)[ , 1] )

R2.sown.seed.only <- rbind(summarise_draws(R2.sown.seed.only.delta15, ~quantile(.x, probs = c(0.05, 0.50, 0.95))),
                           summarise_draws(R2.sown.seed.only.delta19, ~quantile(.x, probs = c(0.05, 0.50, 0.95))),
                           summarise_draws(R2.sown.seed.only.bm, ~quantile(.x, probs = c(0.05, 0.50, 0.95)))
)






# delta richness seeding and disturbance 
R2.rich.seed.disturb.delta15 <- bayesR2_res(fit = stanmod.delta.rich.seed.disturb, ypred = "mu_delta15", 
                                            y = standata.delta.rich.seed.disturb$delta15 )

R2.rich.seed.disturb.delta19 <- bayesR2_res(fit = stanmod.delta.rich.seed.disturb, ypred = "mu_delta19", 
                                            y = standata.delta.rich.seed.disturb$delta19 )

R2.rich.seed.disturb.bm <- bayesR2_res(fit = stanmod.delta.rich.seed.disturb, ypred = "mu_bm", 
                                       y = get_ci(stanmod.delta.rich.seed.disturb, pars = "realbm", 
                                                  probs = 0.5)[ , 1] )

R2.rich.seed.disturb <- rbind(summarise_draws(R2.rich.seed.disturb.delta15, ~quantile(.x, probs = c(0.05, 0.50, 0.95))),
                              summarise_draws(R2.rich.seed.disturb.delta19, ~quantile(.x, probs = c(0.05, 0.50, 0.95))),
                              summarise_draws(R2.rich.seed.disturb.bm, ~quantile(.x, probs = c(0.05, 0.50, 0.95)))
)





# delta sown richness seeding and disturbance 
R2.sown.seed.disturb.delta15 <- bayesR2_res(fit = stanmod.delta.sown.seed.disturb, ypred = "mu_delta15", 
                                            y = standata.delta.sown.seed.disturb$delta15 )

R2.sown.seed.disturb.delta19 <- bayesR2_res(fit = stanmod.delta.sown.seed.disturb, ypred = "mu_delta19", 
                                            y = standata.delta.sown.seed.disturb$delta19 )

R2.sown.seed.disturb.bm <- bayesR2_res(fit = stanmod.delta.sown.seed.disturb, ypred = "mu_bm", 
                                       y = get_ci(stanmod.delta.sown.seed.disturb, pars = "realbm", 
                                                  probs = 0.5)[ , 1] )

R2.sown.seed.disturb <- rbind(summarise_draws(R2.sown.seed.disturb.delta15, ~quantile(.x, probs = c(0.05, 0.50, 0.95))),
                              summarise_draws(R2.sown.seed.disturb.delta19, ~quantile(.x, probs = c(0.05, 0.50, 0.95))),
                              summarise_draws(R2.sown.seed.disturb.bm, ~quantile(.x, probs = c(0.05, 0.50, 0.95)))
)




# and combine R2 for the two treatments and name first column 'parameter'
R2.rich <- cbind(R2.rich.seed.only, R2.rich.seed.disturb[, -1]) 
R2.sown <- cbind(R2.sown.seed.only, R2.sown.seed.disturb[, -1]) 

# put nice labels
R2.rich$variable <- c("$R^2$ $\\Delta$ richness 1st year", "$R^2$ $\\Delta$ richness 5th year", "$R^2$ productivity")
R2.sown$variable <- c("$R^2$ $\\Delta$ richness 1st year", "$R^2$ $\\Delta$ richness 5th year", "$R^2$ productivity")



# and add blank 'N_eff' column to merge with parameter credible intervals
R2.rich <- cbind(R2.rich[ , 1:4], N_eff = rep("", nrow(R2.rich)), 
                 R2.rich[ , 5:7], N_eff = rep("", nrow(R2.rich)) )

R2.sown <- cbind(R2.sown[ , 1:4], N_eff = rep("", nrow(R2.sown)), 
                 R2.sown[ , 5:7], N_eff = rep("", nrow(R2.sown)) )


# remove 'variable' name and make nice interval labels
names(R2.rich) <- c( "parameter", rep(c("5%", "50%", "95%", "N_eff"), 2) )
names(R2.sown) <- c( "parameter", rep(c("5%", "50%", "95%", "N_eff"), 2) )





# now rbind the Bayes R2 to the parameter estimates
cis.delta.rich <- rbind(cis.delta.rich, R2.rich)
cis.delta.sown <- rbind(cis.delta.sown, R2.sown)



# write to table    
write.table(cis.delta.sown, sep = "\t", dec = ".", 
            row.names = FALSE, quote = FALSE,
            "doc/tables/CIs-establishment-delta-sown-richness.txt")

write.table(cis.delta.rich, sep = "\t", dec = ".", 
            row.names = FALSE, quote = FALSE,
            "doc/tables/CIs-establishment-delta-richness.txt")





# just for for comparison, the second option to calculate Bayes R^2

R2_delta15 <- bayesR2(fit = stanmod.delta.rich.seed.only, ypred = "mu_delta15", 
                      error="sigma_delta15", family = "normal")
summarise_draws(R2_delta15, ~quantile(.x, probs = c(0.05, 0.50, 0.95)))

R2_delta19 <- bayesR2(fit = stanmod.delta.rich.seed.only, ypred = "mu_delta19", 
                      error="sigma_delta19", family = "normal")
summarise_draws(R2_delta19, ~quantile(.x, probs = c(0.05, 0.50, 0.95)))

R2_bm      <- bayesR2(fit = stanmod.delta.rich.seed.only, ypred = "mu_bm",     
                      error = "phi_bm",       family = "Gamma")
summarise_draws(R2_bm, ~quantile(.x, probs = c(0.05, 0.50, 0.95)))







################################################################################
### 2.1 Compare models with brms results                                 #######
################################################################################
###   
###   
###   # Get averaged biomass measurements to have one mean per grassland
###   meanbiomass  <- aggregate(biomass.g ~ plot, mean, data = latbiomass)[2]
###   meanbiomass <- data.frame(plot = levels(latbiomass$plot),
###                             meanbiomass = meanbiomass )
###   names(meanbiomass)[2] <- "meanbiomass"
###   
###   metadata <- merge(metadata, meanbiomass, by = "plot", all.x = TRUE)
###   
###   
###   # quickly scale land-use predictors 
###   metadata$mow.scl <- metadata$Mow / sd(metadata$Mow)
###   metadata$log.graz.scl <- log(metadata$Graz + 1) / sd(log(metadata$Graz + 1) )
###   metadata$log.fert.scl <- log(metadata$Fert + 1) / sd(log(metadata$Fert + 1) )
###   metadata$logmeanbiomass.scl <- log(metadata$meanbiomass) / sd(log(metadata$meanbiomass))
###   
###   
###   # set the data frames
###   brmdata.rich.seed.only <- cbind(subset(metadata, treatment == "seed" & year == 2019), 
###                                   delta15 = standata.delta.rich.seed.only$delta15, 
###                                   delta19 = standata.delta.rich.seed.only$delta19)
###   brmdata.sown.seed.only <- cbind(subset(metadata, treatment == "seed" & year == 2019), 
###                                   delta15 = standata.delta.sown.seed.only$delta15, 
###                                   delta19 = standata.delta.sown.seed.only$delta19)
###   brmdata.rich.seed.disturb <- cbind(subset(metadata, treatment == "seeddisturb" & year == 2019), 
###                                   delta15 = standata.delta.rich.seed.disturb$delta15, 
###                                   delta19 = standata.delta.rich.seed.disturb$delta19)
###   brmdata.sown.seed.disturb <- cbind(subset(metadata, treatment == "seeddisturb" & year == 2019), 
###                                   delta15 = standata.delta.sown.seed.disturb$delta15, 
###                                   delta19 = standata.delta.sown.seed.disturb$delta19)
###   
###   # define the model formulas
###   bf.delta15 <- bf(delta15 ~ region + mow.scl + log.graz.scl + logmeanbiomass.scl)
###   bf.delta19 <- bf(delta19 ~ region + mow.scl + log.graz.scl + logmeanbiomass.scl)
###   bf.biomass <- bf(meanbiomass ~ region + log.fert.scl, family = Gamma(link="log"))
###   
###   
###   # and run the models
###   
###   brms.rich.seed.only <- brm( bf.delta15 + bf.delta19 + bf.biomass, 
###                               data = brmdata.rich.seed.only, 
###                               chains = chains, iter = 2000)
###   
###   brms.sown.seed.only <- brm( bf.delta15 + bf.delta19 + bf.biomass, 
###                               data = brmdata.sown.seed.only, 
###                               chains = chains, iter = 2000)
###   
###   brms.rich.seed.disturb <- brm( bf.delta15 + bf.delta19 + bf.biomass, 
###                               data = brmdata.rich.seed.disturb, 
###                               chains = chains, iter = 2000)
###   
###   brms.sown.seed.disturb <- brm( bf.delta15 + bf.delta19 + bf.biomass, 
###                               data = brmdata.rich.seed.disturb, 
###                               chains = chains, iter = 2000)
###   
###   brms.rich.seed.only
###   brms.sown.seed.only
###   brms.rich.seed.disturb
###   brms.sown.seed.disturb
###   
###   
###   



################################################################################
### 3 predict number of established species along the productivity gradient ####
###   This is going to be lengthy, because there is no                      ####
###   'posterior_linpred' or predict-like function in stan                  ####
################################################################################


# source predict function
source("R/01-2_predict-delta.R")


predict.rich.seed.only.year5 <-    predict_established_year5(stanmod.delta.rich.seed.only, 
                                                             standata.delta.rich.seed.only$delta19, 
                                                             iter = iter, chains = chains)

predict.sown.seed.only.year5 <-    predict_established_year5(stanmod.delta.sown.seed.only, 
                                                            standata.delta.sown.seed.only$delta19, 
                                                            iter = iter, chains = chains)

predict.rich.seed.disturb.year5 <- predict_established_year5(stanmod.delta.rich.seed.disturb, 
                                                            standata.delta.rich.seed.disturb$delta19, 
                                                            iter = iter, chains = chains)

predict.sown.seed.disturb.year5 <- predict_established_year5(stanmod.delta.sown.seed.disturb, 
                                                            standata.delta.sown.seed.disturb$delta19, 
                                                            iter = iter, chains = chains)


save(predict.rich.seed.only.year5,
     predict.sown.seed.only.year5,
     predict.rich.seed.disturb.year5,
     predict.sown.seed.disturb.year5, 
     file="Rdata/01-2_stan-delta-richness-pred.RData")




################################################################################
### 2 Plot posterior predictive density overlay                          #######
################################################################################



# check posterior predictive distributions: did the model behave well, does it fit the data?
# extract posterior predictive samples


# seeding only treatment, delta richness
ppc.rich.seed.only.15 <- ppc_dens_overlay_stanfit(stanmod.delta.rich.seed.only, pars = "post_delta15", 
                                                  y = standata.delta.rich.seed.only$delta15)  +
        labs(title = expression(  Response ~ italic(y): Delta ~ richness ~ year ~ one, seeding ~ only ),
             subtitle = expression (Seeding ~ only ) )



ppc.rich.seed.only.19 <- ppc_dens_overlay_stanfit(stanmod.delta.rich.seed.only, pars = "post_delta19", 
                                                  y = standata.delta.rich.seed.only$delta19) +
        labs(title = expression(  Response ~ italic(y): Delta ~ richness ~ year ~ five) ,
             subtitle = expression (Seeding ~ only ) )



# seeding only treatment, delta sown richness
ppc.sown.seed.only.15 <- ppc_dens_overlay_stanfit(stanmod.delta.sown.seed.only, pars = "post_delta15", 
                                                  y = standata.delta.sown.seed.only$delta15) +
        labs(title = expression(  Response ~ italic(y): Delta ~ sown ~ richness ~ year ~ one ),
             subtitle = expression (Seeding ~ only ) )


ppc.sown.seed.only.19 <- ppc_dens_overlay_stanfit(stanmod.delta.sown.seed.only, pars = "post_delta19", 
                                                  y = standata.delta.sown.seed.only$delta19)+
        labs(title = expression(  Response ~ italic(y): Delta ~ sown ~ richness ~ year ~ five ),
             subtitle = expression (Seeding ~ only ) )



# seeding and disturbance treatment, delta richness
ppc.rich.seed.disturb.15 <- ppc_dens_overlay_stanfit(stanmod.delta.rich.seed.disturb, pars = "post_delta15", 
                                                     y = standata.delta.rich.seed.disturb$delta15)+
        labs(title = expression(  Response ~ italic(y): Delta ~ richness ~ year ~ one),
             subtitle = expression (Seeding ~ and ~ disturbance ) )


ppc.rich.seed.disturb.19 <- ppc_dens_overlay_stanfit(stanmod.delta.rich.seed.disturb, pars = "post_delta19", 
                                                     y = standata.delta.rich.seed.disturb$delta19)+
        labs(title = expression(  Response ~ italic(y): Delta ~ richness ~ year ~ five),
             subtitle = expression (Seeding ~ and ~ disturbance ) )


# seeding and disturbance treatment, delta sown richness
ppc.sown.seed.disturb.15 <- ppc_dens_overlay_stanfit(stanmod.delta.sown.seed.disturb, pars = "post_delta15", 
                                                     y = standata.delta.sown.seed.disturb$delta15) +
        labs(title = expression(  Response ~ italic(y): Delta ~ sown ~ richness ~ year ~ one),
             subtitle = expression (Seeding ~ and ~ disturbance ) )

ppc.sown.seed.disturb.19 <- ppc_dens_overlay_stanfit(stanmod.delta.sown.seed.disturb, pars = "post_delta19", 
                                                     y = standata.delta.sown.seed.disturb$delta19) +
        labs(title = expression(  Response ~ italic(y): Delta ~ sown ~ richness ~ year ~ five),
             subtitle = expression (Seeding ~ and ~ disturbance ) )


# and the productivity model - only once, as it was fitted the same way
# in all four models
# use the median of the moelled latent productivity as y data
ppc.productivity <- ppc_dens_overlay_stanfit(stanmod.delta.sown.seed.disturb, pars = "post_bm", 
                                             y = get_ci(stanmod.delta.sown.seed.disturb, pars = "realbm", 
                                                        probs = 0.5)[ , 1] ) +
        labs(title = expression(  Response ~ italic(y): Productivity ) ) 

# combine the plots
ggarrange(ppc.rich.seed.only.15, 
          ppc.rich.seed.only.19,
          ppc.rich.seed.disturb.15,
          ppc.rich.seed.disturb.19,
          ppc.sown.seed.only.15,
          ppc.sown.seed.only.19,
          ppc.sown.seed.disturb.15,
          ppc.sown.seed.disturb.19,
          ppc.productivity,
          labels = c("(a)", "(b)", "(c)", "(d)", 
                     "(d)", "(e)", "(f)", "(g)",
                     "(h)"), 
          nrow = 3, ncol = 3,
          hjust = c(0, 0),
          vjust = c(1.5, 1.5),
          font.label = list(size = 8),
          common.legend = TRUE, legend = "bottom")



ggsave("doc/figs/S1_pp-check-multivariate.pdf", plot = last_plot(), device = "pdf",
       width = 18, height = 18, units = "cm", dpi = 1200)

ggsave("doc/figs/S1_pp-check-multivariate.png", plot = last_plot(), device = "png",
       width = 18, height = 18, units = "cm", dpi = 1200)



