
################################################################################
# First of all, get the compiled data, that includes almost all data ###########
# measured in the Seed Addition and Disturbance Experiment SADE      ###########
################################################################################





# load some packages. Add version numbers as reported in the manuscript


library(brms)
library(bayesplot)
library(ggplot2)
library(ggpubr)
library(bayestestR) #function to extract effective sample sizes from brmsfit objects

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


# load the species cover data from our vegetation surveys, long format
species.long <- read.table("data/species-long.txt", header = TRUE,
                           sep = "\t", dec = ".", stringsAsFactors = TRUE)



# convert species data to wide format to calculate diversity indices
species.wide <- tidyr::pivot_wider(species.long, 
                                   names_from = species, 
                                   values_from = cover, 
                                   values_fill = 0)



# first make sure data frames in the same order
identical(metadata$plotIdent, species.wide$plotIdent)



# Three grasslands were not accessible in 2018 (3 grasslands * 4 treatments) 
# This gives 73 grasslands * 4 treatments * 5 years = 1,460 obs. - 12 missing obs.

nrow(metadata) #1,448 observations
nrow(species.wide) # same number of observations


# it is safe now to remove the 'plotIdent' column from the species data frame
species.wide <- species.wide[ , -1]

any(is.na(species.wide)) # no NA's in the data



###   # if you need species richness of sown species, take this one:
###   
###   # we source a script to prepare the 'sown.long' data frame, 
###   # which contains species cover data of only sown species. 
###   # Species that were already present on control and 
###   # disturbance treatments were removed. The full 'species.long' 
###   # data frame is sourced within the script as well.
###   
###   source("R/00-1_sown-and-established-species-data-wrangling.R") 
###   
###   sown.wide <- tidyr::pivot_wider(sown.long, 
###                                   names_from = species, 
###                                   values_from = cover, 
###                                   values_fill = 0) 
###   
###   # remove metadata here: plotIdent, plot, region, year and treatment
###   sown.wide <- sown.wide[ , -c(1:5)] 
###   
###   # replace species dataset
###   species.wide <- sown.wide


################################################################################
### 0.2 Calculate diversity metrics                                      #######
################################################################################


# we'll go with species richness and the Spie (effective number of species, 
# 1 - Simpson index (also called Gini-Simpson index), Jost 2006 and 
# Chase et al. 2018. Compare to formerly used Shannon entropy


# species richness
metadata$richness <- apply(species.wide, MARGIN = 1, function(x) sum(x > 0)) 



# Gini-Simpson index (probability of interspecific encounter, Jost 2006)
metadata$pie     <- apply(species.wide, MARGIN = 1, function(x) 
                               1 - sum( ( x / sum(x) )^2 )
                             # 1 - sum of    p_i^2
                               ) 

# Effective number of species S_pie (equivalent to 1/Simpson concentration)
metadata$spie  <- 1 / ( 1 - metadata$pie )



# plot the diversity indices
pairs(metadata[c("richness", "pie", "spie")], 
      upper.panel = panel.cor, lower.panel = panel.smooth)

# two S_pie values are very high: these are two very species-rich disturbance 
# treatment plots surveyed in 2015, 6 months after disturbance. The plots had 
# very low total cover, and almost all species were noted with 0.1% cover. 





################################################################################
### 1. To what degree do seeding, disturbance, their interaction,         ######
###    and their interactions with time affect plant diversity?           ######
################################################################################


# create the region:year variable to model the varying intercepts. Would
# also be possible to specify (1|region:year) in the model formula, but 
# the regionYear variable makes i easier to predict fitted responses later on

metadata$region.year <- interaction(metadata$region, metadata$year)

# and set year to factor
metadata$year        <- as.factor(metadata$year)



################################################################################
### 1.1 Multilevel model with species richness as the response            ######
################################################################################

chains <- 4
iter <- 10000


# set adapt_delta to a higher value, as in earlier models I had 
# some few divergent transitions
treat.eff.rich <- brm(richness ~ dist*seed*year + (1|plot) + (1|region.year),
                      data = metadata, 
                      family = poisson, cores = parallel::detectCores(), 
                      chains = chains, iter = iter,
                      control = list(adapt_delta = 0.99, max_treedepth = 12), 
                      seed = 123, 
                      prior = set_prior("normal(0,2)", class="b"))
treat.eff.rich


# get posterior summary (rounded) for 'fixed' effects
round(fixef(treat.eff.rich, probs = c(0.05, 0.5, 0.95)),2)
round(ranef(treat.eff.rich, probs = c(0.05, 0.5, 0.95)),2)

# and full posterior summary
round(posterior_summary(treat.eff.rich, probs = c(0.05, 0.5, 0.95)),2)


# check model performances with posterior predictive checks with 
# the bayesplot package

pp_check(treat.eff.rich, nsamples=50) 
pp_check(treat.eff.rich, nsamples = 50, type = "scatter_avg_grouped", 
         group = "region.year") 
pp_check(treat.eff.rich, type = "loo_intervals") 
pp_check(treat.eff.rich, type = "loo_pit") 
# posterior-predictive checks good, but slightly underdispersed (loo_pit)


# conditional Bayes R2 (Gelman et al. 2019)
R2.cond.rich <- brms::bayes_R2(treat.eff.rich, re_formula =  ~ (1|plot) + (1|region.year),
                               probs = c(0.05, 0.5, 0.95) ) 
# conditional R2 sounds high, but is close to the conditional R2 of the MuMIn 
# package (calculated using an equivalent model using the lme4 package).
# Most likely a result of the large number of varying intercepts, which
# explain a lot of the variation in the data

#marginal R2 (i.e. due to 'fixed' effects: treatments and year)
R2.marg.rich <- brms::bayes_R2(treat.eff.rich, re_formula = 1 ~ 1,
                               probs = c(0.05, 0.5, 0.95) ) 


R2.rich <- data.frame(rbind(R2.cond = R2.cond.rich[-c(1:2)],
                            R2.marg = R2.marg.rich[-c(1:2)] ) )





################################################################################
### 1.2 Multilevel model with the effective number of species (Spie)      ######
###     as the response                                                   ######
################################################################################


treat.eff.spie <- brm(spie ~ dist * seed * year + (1|plot) + (1|region.year) ,
                      data = metadata, 
                      family = Gamma(link = "log"), cores = parallel::detectCores(), 
                      chains = chains, iter = iter,
                      control = list(adapt_delta = 0.95, max_treedepth = 12), 
                      seed = 123, 
                      prior = set_prior("normal(0, 2)", class = "b"))
treat.eff.spie


# get posterior summary (rounded) for 'fixed' effects
round(fixef(treat.eff.spie, probs = c(0.05, 0.5, 0.95)),2)

# and full posterior summary
round(posterior_summary(treat.eff.spie, probs = c(0.05, 0.5, 0.95)),2)


# check model performances with posterior predictive checks with 
# the bayesplot package

pp_check(treat.eff.spie, nsamples = 100) 
pp_check(treat.eff.spie, nsamples = 100, type = "scatter_avg_grouped", 
         group = "region.year") 
pp_check(treat.eff.spie, type = "loo_intervals") 
pp_check(treat.eff.spie, type = "loo_pit") 


# conditional Bayes R2 (Gelman et al. 2019)
R2.cond.spie <- brms::bayes_R2(treat.eff.spie, re_formula = ~ (1|plot) + (1|region.year),
                               probs = c(0.05, 0.5, 0.95) ) 
# conditional R2 sounds high, but is close to the conditional R2 of the MuMIn 
# package (calculated using an equivalent model using the lme4 package).
# Most likely a result of the large number of varying intercepts, which
# explain a lot of the variation in the data

#marginal R2 (i.e. due to 'fixed' effects: treatments and year)
R2.marg.spie <- brms::bayes_R2(treat.eff.spie, re_formula = ~ 1,
                               probs = c(0.05, 0.5, 0.95) ) 


R2.spie <- data.frame(rbind(R2.cond = R2.cond.spie[-c(1:2)],
                            R2.marg = R2.marg.spie[-c(1:2)] ) )


################################################################################
### 2 save credible intervals                                             ######
################################################################################


# first the species richness model

# save only the first 22 parameters, omit posterior of individual varying intercepts
posterior.rich <- as.matrix(treat.eff.rich)[ , c(1:22)]

# and get 90% credible interval
posterior.rich <- apply(posterior.rich, MARGIN = 2, function(x) 
  quantile(x, probs = c(0.05, 0.5, 0.95)) )                        

posterior.rich <- t(posterior.rich)

posterior.rich.ess <- effective_sample(treat.eff.rich, effects = c("all"))[ c(1:22) ,  ]

posterior.rich <- cbind(posterior.rich, 
                        N_eff = c(posterior.rich.ess[ , 2 ])
                        )




# second  the S_PIE model

# only the first 22 parameters are interesting here
posterior.spie <- as.matrix(treat.eff.spie)[ , c(1:22)]

# and get 90% credible interval
posterior.spie <- apply(posterior.spie, MARGIN = 2, function(x) 
  quantile(x, probs = c(0.05, 0.5, 0.95)) )                     

posterior.spie <- t(posterior.spie)

posterior.spie.ess <- effective_sample(treat.eff.spie, effects = c("all"))[ c(1:22) , 2 ]

posterior.spie <- cbind(posterior.spie, 
                        N_eff = posterior.spie.ess)


cis.treatment <- data.frame(cbind(posterior.rich, posterior.spie ) )
names(cis.treatment) <- c(rep(c("5%", "50%", "95%", "N_eff"), 2) )


# add shape parameter for Gamma model, could not get the eff sample size
# of the estimate from the brms summary
cis.treatment <- rbind(cis.treatment, 
                       shape = c(NA, NA, NA, NA, 7.40, 7.88, 8.37, 14113))

# and add Bayes R2's
# add variable names to merge
R2.treatment <- cbind(R2.rich, NA, R2.spie, NA)
names(R2.treatment) <- c(rep(c("5%", "50%", "95%", "N_eff"), 2) )


# add 'parameter variable and rbind together
cis.treatment <- rbind(cis.treatment, R2.treatment)





################################################################################
### 3 save data and model posteriors      ######################################
################################################################################


# write credible intervals to table
write.table(cis.treatment, sep = "\t", dec = ".", 
            row.names = TRUE, quote = FALSE,
            "doc/tables/treatment-cis.txt")


save(metadata, treat.eff.rich, treat.eff.spie,
     file = "Rdata/01-1_brms_posteriors_treatment_effects.Rdata")




