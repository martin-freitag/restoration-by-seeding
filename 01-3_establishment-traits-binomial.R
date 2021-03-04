

################################################################################
###   Do functional traits and productivity                            ######### 
###   predict establishment success?                                   #########
################################################################################


# NOTE: As of January 2020, the Stan program sometimes causes R to crash
#       ('C stack trace' error) when models are repeatedly fitted and saved 
#       to the same object (see also https://github.com/stan-dev/rstan/issues/844).
#       The script however runs fine from top to bottom when you source the 
#       whole script at once and wait until it is finished. 
#       
#       For testing purposes, I would recommend to set 'iter = 2000' to
#       reduce the runtime (and RAM consumption) of the models.



# some functions 
source("R/helpers.R")


# load packages
library(rstan)
library(brms)
library(loo)
# install.packages("posterior", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(posterior)
library(bayesplot)
library(ggpubr)
library(ggplot2)




# load data


# here, we source a script to prepare the 'newly.sown.long' data frame, 
# which contains species cover data of only newly sown species on the 
# two seeding treatments in 2019. Species that were already present on control 
# and disturbance treatments were removed.

source("R/00-1_sown-and-established-species-data-wrangling.R") 



# we are going to use the 'metadata' data, which contains info about  
# the experimental design and metadata like grazing and mowing intensities, 
# but also productivity (aboveground biomass) sampling data
metadata <- read.table("data/metadata.txt", header = TRUE,
                          sep = "\t", dec = ".", stringsAsFactors = TRUE)


# get LEDA traits
traits <- read.table("data/LEDA-traits.txt",
                     sep = "\t", header = TRUE, stringsAsFactors = TRUE)


# get the seeding density and germination rates (in percent), multiply the two to use as 
# 'live seeding density' predictor
seeddensity <- read.table("data/seeding-density-germination-rates.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = TRUE)

# we have germination rates for seeds from species separately for each exploratory,
# but we are missing germination rates for Hypericum perforatum in the ALB region
# take the mean from the other two regions

# calculate mean for Hypericum perforatum
mean.hyp.perf <- mean(seeddensity[ seeddensity$species == "Hypericum_perforatum" ,
                                   ]$germinationrate, 
                      na.rm = TRUE)

# replace missing Hypericum perforatum gerination rate
seeddensity[ seeddensity$species == "Hypericum_perforatum" 
             & seeddensity$region == "AEG", 
             ]$germinationrate <- mean.hyp.perf


seeddensity$live.seeddensity <- seeddensity$seeddensity.per.m2 * (seeddensity$germinationrate / 100)




# merge traits and seeding density into the newly sown species data frame
newly.sown.long <- merge(newly.sown.long, traits, by = "species", all.x = TRUE)
newly.sown.long <- merge(newly.sown.long, seeddensity, by = c("species", "region"), all.x = TRUE)



# order by plot and species, merging forced order by species (plot unordered)
newly.sown.long <- newly.sown.long[order(newly.sown.long$plotIdent, 
                                         newly.sown.long$species) , ]


# convert cover data observations to binary 0/1 observations
newly.sown.long$presence <- ifelse(newly.sown.long$cover > 0, 
                                   1, 0)



# get proportion of already present seeding species in regions
# Basically, it shouldn't matter which year to chose, but
# in 2018 there are some few missing observations in regions Alb and Hainich


1 - nrow(subset(newly.sown.long, region == "AEG" & treatment == "seed" & year == 2019)) / 
  nrow(subset(sown.long, region == "AEG" & treatment == "seed" & year == 2019))
# 37% of sown species observations deleted because species were already present


1 - nrow(subset(newly.sown.long, region == "HEG" & treatment == "seed" & year == 2019)) / 
  nrow(subset(sown.long, region == "HEG" & treatment == "seed" & year == 2019))
# 36% of sown species observations deleted because species were already present

1 - nrow(subset(newly.sown.long, region == "SEG" & treatment == "seed" & year == 2019)) / 
  nrow(subset(sown.long, region == "SEG" & treatment == "seed" & year == 2019))
# 21% of sown species observations deleted because species were already present


################################################################################
### 0.1 Prepare the productivity data and trait predictors              ########
###     Why the strange biomass data format?                            ########
###     We repeatedly measured productivity from 2015-2019,             ########
###     and we can use this data to model a latent (or 'true')          ########
###     productivity with measurement error                             ########
################################################################################



# take the biomass observations from 2015 to 2019 from the control treatment 

latbiomass  <- droplevels(subset(metadata, 
                                 treatment == "control", 
                                 select = c("plot", "biomass.g") ) )

latbiomass <- latbiomass[complete.cases(latbiomass) , ] #remove NA's




# Scale the trait and live seeding density predictors to unit standard deviation
# log height and seed mass, no need to log SLA (cf. pairs plot of predictors later)
newly.sown.long$logheight.scl         <- log(newly.sown.long$height) / sd(log(newly.sown.long$height))
newly.sown.long$logseedMass.scl       <- log(newly.sown.long$seedMass) / sd(log(newly.sown.long$seedMass))
newly.sown.long$SLA.scl               <- newly.sown.long$SLA / sd(newly.sown.long$SLA)
newly.sown.long$live.seeddensity.scl  <- newly.sown.long$live.seeddensity / sd(newly.sown.long$live.seeddensity)




# for comparisons with models accounting for trait-grazing and trait-mowing
# interactions, also add mowing and grazing intensities to the data
newly.sown.long <- merge(newly.sown.long, 
                         unique(metadata[c("plot", "Graz", "Mow")]), 
                         by  = "plot")


# scale also grazing and mowing intensities (used in later models)
newly.sown.long$log.graz.scl <- log(newly.sown.long$Graz + 1) / sd(log(newly.sown.long$Graz + 1))
newly.sown.long$mow.scl  <- newly.sown.long$Mow / sd(newly.sown.long$Mow)



# get averaged biomass to compare the Stan to brms model (verification and debugging)
meanbiomass <- aggregate(biomass.g ~ plot, mean, data = latbiomass)[2]
meanbiomass <- data.frame(plot = factor(levels(latbiomass$plot)),
                          logmeanbiomass = log(meanbiomass) )
names(meanbiomass)[2] <- "logmeanbiomass"


# add the logmeanbiomass variable to the species data frame, 
# and scale this variable (for use in brms models only)
newly.sown.long <- merge(newly.sown.long, meanbiomass, by  = "plot")
newly.sown.long$logmeanbiomass.scl <- newly.sown.long$logmeanbiomass / 
                                          sd(newly.sown.long$logmeanbiomass)






################################################################################
### 1.1 Establishment by species and plot on the seeding treatments     ########
###     Seeding only treatment                                          ########
################################################################################


# drop unused factor levels - necessary if you subset the data
newly.sown.long <- droplevels(newly.sown.long)


# first: subset establishment data to focal treatment and year 2019 
newly.sown.seed.only <- subset(newly.sown.long, treatment == "seed")




# define initial design matrix (interactions with productivity defined within stan model)
# note that predictors are scale to unit SD, but centered in the stan model
X <- with(newly.sown.seed.only,
         data.frame(
           cbind(
             logheight.scl = logheight.scl,
             logseedMass.scl = logseedMass.scl,
             SLA.scl = SLA.scl,
             live.seeddensity.scl = live.seeddensity.scl
             )
         )
)

# look at correlation/relation between trait predictors 
pairs(X, upper.panel = panel.cor, lower.panel = panel.smooth)


standata.traits.seed.only <- list(
  Nobs = nrow(newly.sown.seed.only),     # number of observations
  pres = newly.sown.seed.only$presence,  # binary response variable
  
  NX = ncol(X),  # number of predictors in initial design matrix
  Ntrait = 3,              # number of trait-productivity interactions (3 traits * productivity)
  X = X,         # initial design matrix
  
  # random intercept of plot (73 levels, ignore nesting within three regions)
  Nplot = length(levels(newly.sown.seed.only$plot)), # number of plot 
  Mplot = 1,              # number of effects per level, only intercept here
  plot = as.integer(newly.sown.seed.only$plot), # integer plot identifier
  Xplot = rep(1, nrow(newly.sown.seed.only)), # plot-level design matrix, vector of 1 for intercept only
  
  # random intercept of species (73 levels)
  Nspecies = length(levels(newly.sown.seed.only$species)), # number of species 
  Mspecies = 1,           # number of effects per level, only intercept here
  species = as.integer(newly.sown.seed.only$species), # integer species identifier
  Xspecies = rep(1, nrow(newly.sown.seed.only)), # group-level design matrix, vector of 1 for intercept
  
  # each five biomass observations per plot, with measurement error. Therefore, modelling of latent productivity
  Nlatbm = nrow(latbiomass),   # number of biomass observations (73 plots * 5 years minus 2 NAs)
  plotbm = as.integer(droplevels(latbiomass$plot)),
  latbm = latbiomass$biomass.g, # biomass observations with measurement error
  meansbm = aggregate(biomass.g ~ plot, mean, data = latbiomass)[,2],
  sdbm = aggregate(biomass.g ~ plot, sd, data = latbiomass)[,2]
  )




# second: subset establishment data to other treatment in year 2019 
newly.sown.seed.disturb <- subset(newly.sown.long, treatment == "seeddisturb" )


# the design matrix is identical to the matrix for the 
# seeding only treatment, only the response variable is different

standata.traits.seed.disturb <- standata.traits.seed.only
standata.traits.seed.disturb$pres <- newly.sown.seed.disturb$presence




################################################################################
### 1.2 Run the stan models                                             ########
################################################################################




# list the model parameters to save, because some parts (interaction with latent productivity) are ugly and memory-consuming
save.pars <- c("c_Int", "Intercept", "b", "r_plot", "r_species", "sd_plot", "sd_species", 
               "realbm", "sigma_bm", "sd_log_realbmX", "sd_X_interact", "log_lik", "post")
iter <- 10000
chains <- 4

# Stan sampling options
rstan_options(auto_write = TRUE)


# Sample the Stan models
stan.traits.seed.only <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
                              data = standata.traits.seed.only,
                              chains = chains, cores = chains, 
                              iter = iter, thin = 1, 
                              control = list(adapt_delta = 0.8, max_treedepth = 10),
                              pars = save.pars,
                              seed = 123) 

stan.traits.seed.disturb <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
                                 data = standata.traits.seed.disturb,
                                 chains = chains, cores = chains, 
                                 iter = iter, thin = 1, 
                                 control = list(adapt_delta = 0.8, max_treedepth = 10),
                                 pars = save.pars,
                                 seed = 123) 




# parameters to extract for table
pars_model <- c("c_Int", "b", "sd_plot", "sd_species", 
                "sigma_bm")

# get names of design matrices for responses
names_b      <- c(colnames(X), "log-productivity.scl",
                  paste("log-productivity_", colnames(X)[1:3], sep = "") )

# and get the names for model parameters
names_pars_model <- c("c_Int", names_b, "sd_plot", "sd_species", 
                      "sigma_bm")



# extract the posterior summaries, together with n_eff and Rhat

ci.traits.seed.only <- get_ci(stan.traits.seed.only, 
                              pars = pars_model, 
                              names = names_pars_model)

ci.traits.seed.disturb <- get_ci(stan.traits.seed.disturb, 
                                 pars = pars_model, 
                                 names = names_pars_model)



# Compute elpd and LOO IC to compare models

loo.seed.only     <- loo(stan.traits.seed.only)
loo.seed.disturb  <- loo(stan.traits.seed.disturb)




# save and remove models, they are quite large...

save(stan.traits.seed.only, stan.traits.seed.disturb,
     standata.traits.seed.only, standata.traits.seed.disturb,
     file = "Rdata/01-3_stan-posteriors-traits.Rdata")

rm(stan.traits.seed.only, 
   stan.traits.seed.disturb)




################################################################################
### 2 Check whether the results are robust to the distribution           #######
###   of the productivity (log.) predictor: there are three              #######
###   in the ALB region with very low productivity, which could          #######
###   bias the trait-productivity estimates                              #######
################################################################################



# which are the observations?

###   subset.obs <- which(newly.sown.seed.only$plot == "AEG07" | 
###                        newly.sown.seed.only$plot == "AEG09" | 
###                        newly.sown.seed.only$plot == "AEG26")
###   
###   newly.sown.seed.only.subset <- droplevels(newly.sown.seed.only[ - subset.obs , ])
###   
###   # also modify predictor matrix 
###   X.subset <- X[ - subset.obs , ]
###   
###   # ... and biomass data
###   latbiomass.subset <- droplevels(latbiomass[ - which(latbiomass$plot == "AEG07" | 
###                                                         latbiomass$plot == "AEG09" | 
###                                                         latbiomass$plot == "AEG26") , ])
###   
###     
###   # prepare input data
###   
###   standata.traits.seed.only.subset <- list(
###     Nobs = nrow(newly.sown.seed.only.subset),     # number of observations
###     pres = newly.sown.seed.only.subset$presence,  # binary response variable
###     
###     NX = ncol(X.subset),  # number of predictors in initial design matrix
###     Ntrait = 3,              # number of trait-productivity interactions (3 traits * productivity)
###     X = X.subset,         # initial design matrix
###     
###     # random intercept of plot (73 levels, ignore nesting within three regions)
###     Nplot = length(levels(newly.sown.seed.only.subset$plot)), # number of plot 
###     Mplot = 1,              # number of effects per level, only intercept here
###     plot = as.integer(newly.sown.seed.only.subset$plot), # integer plot identifier
###     Xplot = rep(1, nrow(newly.sown.seed.only.subset)), # plot-level design matrix, vector of 1 for intercept only
###     
###     # random intercept of species (73 levels)
###     Nspecies = length(levels(newly.sown.seed.only.subset$species)), # number of species 
###     Mspecies = 1,           # number of effects per level, only intercept here
###     species = as.integer(newly.sown.seed.only.subset$species), # integer species identifier
###     Xspecies = rep(1, nrow(newly.sown.seed.only.subset)), # group-level design matrix, vector of 1 for intercept
###     
###     # each five biomass observations per plot, with measurement error. Therefore, modelling of latent productivity
###     Nlatbm = nrow(latbiomass.subset),   # number of biomass observations (73 plots * 5 years minus 2 NAs)
###     plotbm = as.integer(droplevels(latbiomass.subset$plot)),
###     latbm = latbiomass.subset$biomass.g, # biomass observations with measurement error
###     meansbm = aggregate(biomass.g ~ plot, mean, data = latbiomass.subset)[,2],
###     sdbm = aggregate(biomass.g ~ plot, sd, data = latbiomass.subset)[,2]
###   )
###   
###   
###   
###   
###   # second: subset establishment data to other treatment in year 2019 
###   newly.sown.seed.disturb.subset <- newly.sown.seed.disturb[ - subset.obs , ]
###   
###   
###   # the design matrix is identical to the matrix for the 
###   # seeding only treatment, only the response variable is different
###   
###   standata.traits.seed.disturb.subset <- standata.traits.seed.only.subset
###   standata.traits.seed.disturb.subset$pres <- newly.sown.seed.disturb.subset$presence
###   
###   
###   
###   # sample the Stan models for subsetted data
###   
###   stan.traits.seed.only.subset <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
###                                        data = standata.traits.seed.only.subset,
###                                        chains = chains, cores = chains, 
###                                        iter = iter, thin = 1, 
###                                        control = list(adapt_delta = 0.8, max_treedepth = 10),
###                                        pars = save.pars,
###                                        seed = 123) 
###   
###   stan.traits.seed.disturb.subset <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
###                                           data = standata.traits.seed.disturb.subset,
###                                           chains = chains, cores = chains, 
###                                           iter = iter, thin = 1, 
###                                           control = list(adapt_delta = 0.8, max_treedepth = 10),
###                                           pars = save.pars,
###                                           seed = 123) 
###   
###   
###   
###    # Alternative possibility to check whether the grasslands with very low
###    # productivity are highly influential:
###    #
###    # Sample the Stan models, but with another .stan file without 
###    # grassland-specific informed priors for the latent productivity intercepts: 
###    # this model file uses global priors for mean productivity and hence 
###    # 'shrinks' extreme (here: low) productivity estimates towards the global mean
###    # 
###    # -> check if model results with informed priors are robust
###    
###    stan.traits.seed.only.shrink <- stan(file = "Stan_files/01-3_traits-binomial-productivity-shrink.stan", 
###                                         data = standata.traits.seed.only,
###                                         chains = chains, cores = chains, 
###                                         iter = iter, thin = 1, 
###                                         control = list(adapt_delta = 0.8, max_treedepth = 10),
###                                         pars = save.pars,
###                                         seed = 123) 
###    
###    stan.traits.seed.disturb.shrink <- stan(file = "Stan_files/01-3_traits-binomial-productivity-shrink.stan", 
###                                            data = standata.traits.seed.disturb,
###                                            chains = chains, cores = chains, 
###                                            iter = iter, thin = 1, 
###                                            control = list(adapt_delta = 0.8, max_treedepth = 10),
###                                            pars = save.pars,
###                                            seed = 123) 
###    
###    # Compute elpd and LOO IC to compare models
###    
###    loo.seed.only.subset     <- loo(stan.traits.seed.only.subset)
###    loo.seed.disturb.subset  <- loo(stan.traits.seed.disturb.subset)
###    
###    
###    # save and remove the models to save space
###    
###    save(stan.traits.seed.only.subset, stan.traits.seed.disturb.subset,
###         file = "Rdata/01-3_stan-posteriors-traits-productivity-subset.Rdata")
###    
###    rm(stan.traits.seed.only.subset, stan.traits.seed.disturb.subset)
###   
###   
###   
###   
###   # parameters to extract for table
###   pars_model <- c("c_Int", "b", "sd_plot", "sd_species", 
###                   "sigma_bm")
###   
###   # get names of design matrices for responses
###   names_b      <- c(colnames(X.subset), "log-productivity.scl",
###                     paste("log-productivity_", colnames(X.subset)[1:3], sep = "") )
###   
###   # and get the names for model parameters
###   names_pars_model <- c("c_Int", names_b, "sd_plot", "sd_species", 
###                         "sigma_bm")
###   
###   
###   
###   # extract the posterior summaries, together with n_eff and Rhat
###   
###   ci.traits.seed.only.subset <- get_ci(stan.traits.seed.only.subset, 
###                                       pars = pars_model, 
###                                       names = names_pars_model)
###   
###   ci.traits.seed.disturb.subset <- get_ci(stan.traits.seed.disturb.subset, 
###                                           pars = pars_model, 
###                                           names = names_pars_model)
###   
###   # Compute elpd and LOO IC to compare models
###   
###   loo.seed.only.subset     <- loo(stan.traits.seed.only.subset)
###   loo.seed.disturb.subset  <- loo(stan.traits.seed.disturb.subset)
###   
###   
###   # save and remove the models to save space
###   
###   save(stan.traits.seed.only.subset, stan.traits.seed.disturb.subset,
###        file = "Rdata/01-3_stan-posteriors-traits-productivity-subset.Rdata")
###   
###   rm(stan.traits.seed.only.subset, stan.traits.seed.disturb.subset)
###   


################################################################################
### 3 Compare the Stan model to brms and mean biomass                   ########
################################################################################


###   # seeding only   
###   
###   brms.traits.seed.only <- brm( presence ~  logheight.scl*logmeanbiomass.scl
###                          + SLA.scl*logmeanbiomass.scl 
###                          + logseedMass.scl*logmeanbiomass.scl 
###                          + live.seeddensity.scl
###                          + (1|plot) + (1|species) ,
###                          family = bernoulli, 
###                          data = newly.sown.seed.only, 
###                          cores = parallel::detectCores(), 
###                          chains = 4, iter = 10000,
###                          control = list(adapt_delta = 0.9, max_treedepth = 12),
###                          seed = 123, 
###                          prior = set_prior("normal(0, 5)", class = "b")
###   ) 
###   
###   brms.traits.seed.only
###   bayes_R2(brms.traits.seed.only)
###   
###   
###   # seeding and disturbance
###   
###   brms.traits.seed.disturb <- brm( presence ~  logheight.scl*logmeanbiomass.scl
###                                 + SLA.scl*logmeanbiomass.scl 
###                                 + logseedMass.scl*logmeanbiomass.scl
###                                 + live.seeddensity.scl
###                                 + (1|plot) + (1|species) ,
###                                 family = bernoulli, 
###                                 data = newly.sown.seed.disturb, 
###                                 cores = parallel::detectCores(), 
###                                 chains = 4, iter = 10000,
###                                 control = list(adapt_delta = 0.9, max_treedepth = 12),
###                                 seed = 123, 
###                                 prior = set_prior("normal(0, 5)", class = "b")
###   ) 
###   
###   brms.traits.seed.disturb
###   bayes_R2(brms.traits.seed.disturb)





################################################################################
### 4 Land use, i.e. grazing and mowing intensities could also have     ########
###   influenced establishment in addition to productivity              ########
################################################################################


# Compare the trait-productivity models with models including:
#   - trait-productivity and trait-grazing and trait-mowing 
#   - trait-productivity and trait-grazing 
#   - trait-productivity                   and trait-mowing 
#
# interactions, and use the LOO information criterion 
# to judge which model 'performs best': (http://mc-stan.org/loo/articles/loo2-with-rstan.html)




################################################################################
### 4.1 Prepare the predictor matrices X for the alternative models     ########
################################################################################


# seeding only treatment

standata.traits.seed.only.mow.graz <- standata.traits.seed.only
standata.traits.seed.only.graz     <- standata.traits.seed.only
standata.traits.seed.only.mow      <- standata.traits.seed.only


# add both mowing and grazing as well as their interactions with traits (cols 1-3 in X)
standata.traits.seed.only.mow.graz$X <- cbind(standata.traits.seed.only$X, 
                                              log.graz.scl = newly.sown.seed.only$log.graz.scl, 
                                              mow.scl = newly.sown.seed.only$mow.scl,
                                              log.graz = standata.traits.seed.only$X[, 1:3] * newly.sown.seed.only$log.graz.scl,
                                              mow = standata.traits.seed.only$X[, 1:3] * newly.sown.seed.only$mow.scl
                                              )

standata.traits.seed.only.graz$X <- cbind(standata.traits.seed.only$X, 
                                          log.graz.scl = newly.sown.seed.only$log.graz.scl,
                                          log.graz = standata.traits.seed.only$X[, 1:3] * newly.sown.seed.only$log.graz.scl
                                          )

standata.traits.seed.only.mow$X <- cbind(standata.traits.seed.only$X, 
                                         mow.scl  = newly.sown.seed.only$mow.scl,
                                         mow = standata.traits.seed.only$X[, 1:3] * newly.sown.seed.only$mow.scl
                                         )





# seeding and disturbance treatment, same as before

standata.traits.seed.disturb.mow.graz <- standata.traits.seed.disturb
standata.traits.seed.disturb.graz     <- standata.traits.seed.disturb
standata.traits.seed.disturb.mow      <- standata.traits.seed.disturb

standata.traits.seed.disturb.mow.graz$X <- standata.traits.seed.only.mow.graz$X
standata.traits.seed.disturb.graz$X     <- standata.traits.seed.only.graz$X
standata.traits.seed.disturb.mow$X      <- standata.traits.seed.only.mow$X



# also change 'NX', the number of predictors within the matrix
# NX is used in some parts of the .stan model file

# seeding only
standata.traits.seed.only.mow.graz$NX <- ncol(standata.traits.seed.only.mow.graz$X)
standata.traits.seed.only.graz$NX     <- ncol(standata.traits.seed.only.graz$X)
standata.traits.seed.only.mow$NX      <- ncol(standata.traits.seed.only.mow$X)

standata.traits.seed.disturb.mow.graz$NX <- ncol(standata.traits.seed.disturb.mow.graz$X)
standata.traits.seed.disturb.graz$NX     <- ncol(standata.traits.seed.disturb.graz$X)
standata.traits.seed.disturb.mow$NX      <- ncol(standata.traits.seed.disturb.mow$X)





################################################################################
### 4.2 Sample and summarize the alternative seedling only models       ########
################################################################################


stan.traits.seed.only.mow.graz <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
                                       data = standata.traits.seed.only.mow.graz,
                                       chains = chains, cores = chains, 
                                       iter = iter, thin = 1, 
                                       control = list(adapt_delta = 0.8, max_treedepth = 10),
                                       pars = save.pars,
                                       seed = 123) 

stan.traits.seed.only.graz <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
                                   data = standata.traits.seed.only.graz,
                                   chains = chains, cores = chains, 
                                   iter = iter, thin = 1, 
                                   control = list(adapt_delta = 0.8, max_treedepth = 10),
                                   pars = save.pars,
                                   seed = 123) 

stan.traits.seed.only.mow <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
                                  data = standata.traits.seed.only.mow,
                                  chains = chains, cores = chains, 
                                  iter = iter, thin = 1, 
                                  control = list(adapt_delta = 0.8, max_treedepth = 10),
                                  pars = save.pars,
                                  seed = 123) 

# extract the posterior summaries, together with n_eff and Rhat

# get names of design matrices for responses
names_b      <- c(colnames(standata.traits.seed.only.mow.graz$X), 
                  "log-productivity.scl",
                  paste("log-productivity_", colnames(X)[1:3], sep = "") )

# and get the names for model parameters
names_pars_model <- c("c_Int", names_b, "sd_plot", "sd_species", 
                      "sigma_bm")

ci.traits.seed.only.mow.graz <- get_ci(stan.traits.seed.only.mow.graz, 
                                       pars = pars_model, 
                                       names = names_pars_model)


# get names of design matrices for responses
names_b      <- c(colnames(standata.traits.seed.only.graz$X), 
                  "log-productivity.scl",
                  paste("log-productivity_", colnames(X)[1:3], sep = "") )

# and get the names for model parameters
names_pars_model <- c("c_Int", names_b, "sd_plot", "sd_species", 
                      "sigma_bm")

ci.traits.seed.only.graz <- get_ci(stan.traits.seed.only.graz, 
                                   pars = pars_model, 
                                   names = names_pars_model)


# get names of design matrices for responses
names_b      <- c(colnames(standata.traits.seed.only.mow$X), 
                  "log-productivity.scl",
                  paste("log-productivity_", colnames(X)[1:3], sep = "") )

# and get the names for model parameters
names_pars_model <- c("c_Int", names_b, "sd_plot", "sd_species", 
                      "sigma_bm")

ci.traits.seed.only.mow <- get_ci(stan.traits.seed.only.mow, 
                                  pars = pars_model, 
                                  names = names_pars_model)



# Compute elpd and LOO IC to compare models

loo.seed.only.mow.graz <- loo(stan.traits.seed.only.mow.graz)
loo.seed.only.graz     <- loo(stan.traits.seed.only.graz)
loo.seed.only.mow      <- loo(stan.traits.seed.only.mow)


# then save and delete the seeding models (too large files, RAM for loo calculations)
save(stan.traits.seed.only.mow.graz, stan.traits.seed.only.graz, stan.traits.seed.only.mow,
     file = "Rdata/01-3_stan-posteriors-traits-seed-only-mow-graz.Rdata")

rm(stan.traits.seed.only.mow.graz, 
   stan.traits.seed.only.graz, 
   stan.traits.seed.only.mow)







################################################################################
### 4.3 Sample and summarize the alternative seedling and               ########
###     disturbance models                                              ########
################################################################################



stan.traits.seed.disturb.mow.graz <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
                                     data = standata.traits.seed.disturb.mow.graz,
                                     chains = chains, cores = chains, 
                                     iter = iter, thin = 1, 
                                     control = list(adapt_delta = 0.8, max_treedepth = 10),
                                     pars = save.pars,
                                     seed = 123) 

stan.traits.seed.disturb.graz <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
                                      data = standata.traits.seed.disturb.graz,
                                      chains = chains, cores = chains, 
                                      iter = iter, thin = 1, 
                                      control = list(adapt_delta = 0.8, max_treedepth = 10),
                                      pars = save.pars,
                                      seed = 123) 

stan.traits.seed.disturb.mow <- stan(file = "Stan_files/01-3_traits-binomial.stan", 
                                     data = standata.traits.seed.disturb.mow,
                                     chains = chains, cores = chains, 
                                     iter = iter, thin = 1, 
                                     control = list(adapt_delta = 0.8, max_treedepth = 10),
                                     pars = save.pars,
                                     seed = 123) 



# extract the posterior summaries, together with n_eff and Rhat

# get names of design matrices for responses
names_b      <- c(colnames(standata.traits.seed.disturb.mow.graz$X), 
                  "log-productivity.scl",
                  paste("log-productivity_", colnames(X)[1:3], sep = "") )

# and get the names for model parameters
names_pars_model <- c("c_Int", names_b, "sd_plot", "sd_species", 
                      "sigma_bm")

ci.traits.seed.disturb.mow.graz <- get_ci(stan.traits.seed.disturb.mow.graz, 
                                          pars = pars_model, 
                                          names = names_pars_model)


# get names of design matrices for responses
names_b      <- c(colnames(standata.traits.seed.disturb.graz$X), 
                  "log-productivity.scl",
                  paste("log-productivity_", colnames(X)[1:3], sep = "") )

# and get the names for model parameters
names_pars_model <- c("c_Int", names_b, "sd_plot", "sd_species", 
                      "sigma_bm")

ci.traits.seed.disturb.graz <- get_ci(stan.traits.seed.disturb.graz, 
                                      pars = pars_model, 
                                      names = names_pars_model)


# get names of design matrices for responses
names_b      <- c(colnames(standata.traits.seed.disturb.mow$X), 
                  "log-productivity.scl",
                  paste("log-productivity_", colnames(X)[1:3], sep = "") )

# and get the names for model parameters
names_pars_model <- c("c_Int", names_b, "sd_plot", "sd_species", 
                      "sigma_bm")

ci.traits.seed.disturb.mow <- get_ci(stan.traits.seed.disturb.mow, 
                                     pars = pars_model, 
                                     names = names_pars_model)


# Compute elpd and LOO IC to compare models

loo.seed.disturb.graz     <- loo(stan.traits.seed.disturb.graz)
loo.seed.disturb.mow      <- loo(stan.traits.seed.disturb.mow)
loo.seed.disturb.mow.graz <- loo(stan.traits.seed.disturb.mow.graz)

# then save and delete the seeding and disturbance models (too large files, RAM for loo calculations)
save(stan.traits.seed.disturb.graz, stan.traits.seed.disturb.mow, stan.traits.seed.disturb.mow.graz, 
     file = "Rdata/01-3_stan-posteriors-traits-seed-disturb-mow-graz.Rdata")

rm(stan.traits.seed.disturb.graz, 
   stan.traits.seed.disturb.mow,
   stan.traits.seed.disturb.mow.graz)








################################################################################
### 4.4 Compare alternative models using leave-one-out                  ########
###     cross-validation                                                ########
################################################################################

# see also https://mc-stan.org/loo/reference/loo-glossary.html 
# and https://avehtari.github.io/modelselection/rats_kcv.html


# we used relative effective sample sizes, which allows for better estimates 
# of the PSIS effective sample sizes and Monte Carlo error


# compare models

model.names.seed.only <- c("seed.only", 
                           "seed.only.graz", 
                           "seed.only.mow",
                           "seed.only.mow.graz")

model.names.seed.disturb <- c("seed.disturb", 
                              "seed.disturb.graz", 
                              "seed.disturb.mow",
                              "seed.disturb.mow.graz")

comp.seed.only <- loo_compare(loo.seed.only, 
                              loo.seed.only.graz, 
                              loo.seed.only.mow,
                              loo.seed.only.mow.graz)
print(comp.seed.only)

comp.seed.disturb <- loo_compare(loo.seed.disturb, 
                                 loo.seed.disturb.graz, 
                                 loo.seed.disturb.mow,
                                 loo.seed.disturb.mow.graz)
print(comp.seed.disturb)


# extract information for table - ELPD_loo only
elpd.seed.only    <- comp.seed.only[order(rownames(comp.seed.only)),]
elpd.seed.disturb <- comp.seed.disturb[order(rownames(comp.seed.disturb)),]

elpd <- rbind(elpd.seed.only,
              elpd.seed.disturb)

elpd <- data.frame(elpd[ , c(3, 4, 1, 2)])
elpd$model <- c(model.names.seed.only, 
                model.names.seed.disturb)
elpd <- elpd[ , c(5, 1:4) ]

# make nice rownames
elpd$model <- rep(c("No grazing-/mowing-interactions",
                    "Grazing-interactions",
                    "Mowing-interactions",
                    "Grazing- and mowing-interactions"), 
                  2 )


# write to table
write.table(elpd, sep = "\t", dec = ".", 
            row.names = FALSE, quote = FALSE,
            "doc/tables/trait-comparison-elpd.txt")



rm(loo.seed.only, 
   loo.seed.only.graz,
   loo.seed.only.mow,
   loo.seed.only.mow.graz,
   #loo.seed.only.subset,
   loo.seed.disturb,
   loo.seed.disturb.graz,
   loo.seed.disturb.mow,
   loo.seed.disturb.mow.graz#,
   #loo.seed.disturb.subset
   )



################################################################################
### 5 Predict establishment by traits along the productivity gradient   ########
###   This is going to be lengthy, because there is no                  ########
###   predict-like function in stan                                     ########
################################################################################


# I have put the trait predictions in a separate R script to source 
# this script for both treatments. Script contains a single function
# to predict response-scale establishment probabilities for mean +- SD trait 
# values along the productivity gradient. Also calculates the Bayes R2 on the fly

source("R/01-3_predict-traits.R")


# load the simplest model with only trait-productivity interactions
load("Rdata/01-3_stan-posteriors-traits.Rdata")



pred.traits.seed.only <- predict_traits(model = stan.traits.seed.only,
                                        iter = iter, chains = chains,
                                        X = X,  
                                        newly.sown.data = newly.sown.seed.only)


pred.traits.seed.disturb <- predict_traits(model = stan.traits.seed.disturb,
                                          iter = iter, chains = chains,
                                          X = X,  
                                          newly.sown.data = newly.sown.seed.disturb)

# save predictions for trait values for plotting
save(pred.traits.seed.only, 
     pred.traits.seed.disturb, 
     file = "Rdata/01-3_stan-traits-pred.Rdata")


# get marginal and conditional Bayes R2 to report in the results section
R2.cond.seed.only <- round(quantile(pred.traits.seed.only$R2.cond, probs = c(0.05, 0.5, 0.95)),2)
R2.marg.seed.only <- round(quantile(pred.traits.seed.only$R2.marg, probs = c(0.05, 0.5, 0.95)),2)

R2.cond.seed.disturb <- round(quantile(pred.traits.seed.disturb$R2.cond, probs = c(0.05, 0.5, 0.95)),2)
R2.marg.seed.disturb <- round(quantile(pred.traits.seed.disturb$R2.marg, probs = c(0.05, 0.5, 0.95)),2)


# now merge R2's
R2.seed.only    <- rbind(R2.cond.seed.only, R2.marg.seed.only )
R2.seed.disturb <- rbind(R2.cond.seed.disturb, R2.marg.seed.disturb)


# merge and add empty N_eff column, as the credible intervals have this column
R2.traits <- cbind(parameter = c("R2.cond", "R2.marg") , 
                   R2.seed.only, N_eff = NA, 
                   R2.seed.disturb, N_eff = NA)




################################################################################
### 6 Save model parameter credible intervals as tables                 ########
################################################################################


# cbind the CIs for the two 'small' models  
cis.traits <- cbind(ci.traits.seed.only, ci.traits.seed.disturb[, -1])

# assign nice column names from R2.traits data frame, then rbind
colnames(cis.traits) <- colnames(R2.traits)
cis.traits <- rbind(cis.traits, R2.traits)


# assign nice row names
cis.traits$parameter <- c("Intercept", "Height (log.)", 
                          "Seed mass (log.)", "SLA", 
                          "Live seeding density",
                          "Productivity (log.)", 
                          "Height*Productivity", 
                          "Seed mass*Productivity",
                          "SLA*Productivity", 
                          "SD Grassland", "SD Species", 
                          "$\\sigma$ Productiviy",
                          "$R_{cond}^2$", 
                          "$R_{marg}^2$" )

# move live seeding density (row five) to row two
cis.traits <- cis.traits[ c(1, 5, 2:4, 6:14) , ]

# and save as table
write.table(cis.traits, sep = "\t", dec = ".", 
            row.names = FALSE, quote = FALSE,
            "doc/tables/trait-cis.txt")




###   # cbind the CIs for the two subsetted models  
###   cis.traits.subset <- cbind(ci.traits.seed.only.subset, ci.traits.seed.disturb.subset[, -1])
###   
###   # assign nice column names 
###   colnames(cis.traits.subset) <- colnames(R2.traits)
###   
###   # assign nice row names
###   cis.traits.subset$parameter <- c("Intercept", "Height (log.)", 
###                             "Seed mass (log.)", "SLA", 
###                             "Live seeding density", 
###                             "Productivity (log.)", 
###                             "Height*Productivity", 
###                             "Seed mass*Productivity",
###                             "SLA*Productivity", 
###                             "SD Grassland", "SD Species", 
###                             "$\\sigma$ Productiviy" )
###   
###   
###   # and save as tables
###   
###   write.table(cis.traits.subset, sep = "\t", dec = ".", 
###               row.names = FALSE, quote = FALSE,
###               "doc/tables/trait-cis-subset.txt")



# now the models with trait-mowing and trait-grazing interactions
# merge model credible intervals and replace NAs with empty cells
# column names are duplicated, but it's easy to tell
# by missing entries which model is which

cis.seed.only <- Reduce(function(x,y) merge(x = x, y = y, by = "parameter", all = TRUE), 
                        list(ci.traits.seed.only.graz, 
                             ci.traits.seed.only.mow, ci.traits.seed.only.mow.graz))

cis.seed.disturb <- Reduce(function(x,y) merge(x = x, y = y, by = "parameter", all = TRUE),
                           list(ci.traits.seed.disturb.graz, 
                           ci.traits.seed.disturb.mow, ci.traits.seed.disturb.mow.graz))



# and order rows right and nice... move grassland SD, species SD and biomass measurement error down
cis.seed.only    <- cis.seed.only[ c(1, 2, 11, 12, 20, 3:6, 9, 7, 8, 10, 15, 13, 14, 16:19) , ]
cis.seed.disturb <- cis.seed.disturb[ c(1, 2, 11, 12, 20, 3:6, 9, 7, 8, 10, 15, 13, 14, 16:19) , ]


# make nice interval labels
names(cis.seed.only)[-1] <- rep(c("5%", "50%", "95%", "N_eff"), 3)
names(cis.seed.disturb)[-1] <- rep(c("5%", "50%", "95%", "N_eff"), 3)


# and put nice row names
rownames(cis.seed.only) <- c("Intercept", 
                             "Live seeding density", 
                             "Height (log.)", 
                             "Seed mass (log.)", 
                             "SLA", 
                             "Productivity (log.)", 
                             "Height*Productivity", 
                             "Seed mass*Productivity", 
                             "SLA*Productivity", 
                             "Grazing (log.)", 
                             "Height*Grazing", 
                             "Seed mass*Grazing", 
                             "SLA*Grazing",
                             "Mowing",
                             "Height*Mowing", 
                             "Seed mass*Mowing", 
                             "SLA*Mowing", 
                             "SD Grassland",
                             "SD Species", 
                             "$\\sigma$ Productivity")

# different seeding treatment, but same parameters
rownames(cis.seed.disturb)<- rownames(cis.seed.only)


# and write to table
write.table(cis.seed.only, sep = "\t", dec = ".", 
            row.names = FALSE, quote = FALSE,
            "doc/tables/trait-mow-graz-seeding-only-ci.txt")

write.table(cis.seed.disturb, sep = "\t", dec = ".", 
            row.names = FALSE, quote = FALSE,
            "doc/tables/trait-mow-graz-seeding-disturbance-ci.txt")
