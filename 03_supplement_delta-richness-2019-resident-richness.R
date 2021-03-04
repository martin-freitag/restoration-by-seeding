

# load some packages. Add version numbers as reported in the manuscript

library(data.table)
library(brms)
library(ggplot2)
library(ggpubr)


# functions for Bayes R2 and panel.cor correlation in pairs plots
source("R/helpers.R")





################################################################################
#   0.1 Load and prepare the data                                      #########
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
                                   values_from = cover)

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
###   0.2 Calculate diversity metrics                                    #######
################################################################################


# we'll go with overall species richness and richness of sown species


# species richness
metadata$richness <- apply(species.wide, MARGIN = 1, function(x) sum(x > 0)) 

metadata$sown.richness <- apply(sown.wide, MARGIN = 1, function(x) sum(x > 0)) 



# calculate differences ('Delta') in richness

setDT(metadata) #set to data.table format to use the diff-like function
metadata$plot.year <- interaction(metadata$plot, metadata$year) # identifier


# calculate differences between treatments within plot/year with 
# control as reference, and add 'resident richness' (control treatment)
metadata <- metadata[ , ":=" (delta.rich =  richness - richness[treatment == "control"]  ,
                                    delta.sown =  sown.richness - sown.richness[treatment == "control"],
                                    resident.rich = richness[treatment == "control"]
                                    ) ,
                            by = plot.year ] # differences within plot and year





################################################################################
###   1 Analyse the correlations between increases in species richness   #######
###     and resident species richness                                    #######
################################################################################



# We are going to fit four models:
#
# - seeding only treatment
#   + difference in      species richness between treatment and control
#   + difference in sown species richness between treatment and control
#
# - seeding and disturbance treatment
#   + difference in      species richness between treatment and control
#   + difference in sown species richness between treatment and control



# log-transform resident richness predictor
metadata$log.resident.rich <- log(metadata$resident.rich)


chains <- 4
iter <- 10000


# fit models with brms
fit.rich.seed.only <- brm(delta.rich ~ log.resident.rich + region, 
                       data = subset(metadata, treatment == "seed" & year == 2019), 
                       family = gaussian,
                       cores = chains, chains = chains, iter = iter,
                       control = list(adapt_delta = 0.95, max_treedepth = 12), seed = 123, 
                       prior = set_prior("normal(0, 10)", class="b")
                       )

fit.sown.seed.only <- brm(delta.sown ~ log.resident.rich + region, 
                          data = subset(metadata, treatment == "seed" & year == 2019), 
                          family = gaussian,
                          cores = chains, chains = chains, iter = iter,
                          control = list(adapt_delta = 0.95, max_treedepth = 12), seed = 123, 
                          prior = set_prior("normal(0, 10)", class="b")
                          )

fit.rich.seed.disturb <- brm(delta.rich ~ log.resident.rich + region, 
                          data = subset(metadata, treatment == "seeddisturb" & year == 2019), 
                          family = gaussian,
                          cores = chains, chains = chains, iter = iter,
                          control = list(adapt_delta = 0.95, max_treedepth = 12), seed = 123, 
                          prior = set_prior("normal(0, 10)", class="b")
                          )

fit.sown.seed.disturb <- brm(delta.sown ~ log.resident.rich + region, 
                          data = subset(metadata, treatment == "seeddisturb" & year == 2019), 
                          family = gaussian,
                          cores = chains, chains = chains, iter = iter,
                          control = list(adapt_delta = 0.95, max_treedepth = 12), seed = 123, 
                          prior = set_prior("normal(0, 10)", class="b")
                          )


fit.rich.seed.only
fit.sown.seed.only
fit.rich.seed.disturb
fit.sown.seed.disturb



################################################################################
### 3 save credible intervals                                             ######
################################################################################


# first the richness models

get_brms_ci <- function(brmsfit = brmsfit, probs = c(0.05, 0.5, 0.95) ) {

posterior <- as.matrix(brmsfit)
# save only the first 5 parameters, omit the log-posterior
posterior <- posterior[ , c(1:5)]

# and get 90% credible interval
posterior <- apply(posterior, MARGIN = 2, function(x) 
  quantile(x, probs = c(0.05, 0.5, 0.95) ) )                     
posterior <- t(posterior)

# effective_sample function from bayestestR package
posterior.ess.fix <- effective_sample(brmsfit, effects = c("all"))

posterior <- cbind(posterior, N_eff = c(posterior.ess.fix[ , 2 ], NA ) )

}

# get credible intervals of all models
ci.rich.seed.only    <- get_brms_ci(fit.rich.seed.only)
ci.rich.seed.disturb <- get_brms_ci(fit.rich.seed.disturb)
ci.sown.seed.only    <- get_brms_ci(fit.sown.seed.only)
ci.sown.seed.disturb <- get_brms_ci(fit.sown.seed.disturb)

R2.rich.seed.only    <- brms::bayes_R2(fit.rich.seed.only,    probs = c(0.05, 0.5, 0.95) )[ , c(3:5)]
R2.rich.seed.disturb <- brms::bayes_R2(fit.rich.seed.disturb, probs = c(0.05, 0.5, 0.95) )[ , c(3:5)]
R2.sown.seed.only    <- brms::bayes_R2(fit.sown.seed.only,    probs = c(0.05, 0.5, 0.95) )[ , c(3:5)]
R2.sown.seed.disturb <- brms::bayes_R2(fit.sown.seed.disturb, probs = c(0.05, 0.5, 0.95) )[ , c(3:5)]



# and combine with bayes R2 - first only richness, afterwards sown species richness
cis.rich <- data.frame(cbind(ci.rich.seed.only,
                             ci.rich.seed.disturb) )
names(cis.rich) <- c(rep(c("5%", "50%", "95%", "N_eff"), 2) )

R2.rich <- c(R2.rich.seed.only, NA, R2.rich.seed.disturb, NA)

# rbind together
cis.rich <- rbind(cis.rich, R2 = R2.rich)



# now the sown species richness
cis.sown <- data.frame(cbind(ci.sown.seed.only,
                             ci.sown.seed.disturb) )
names(cis.sown) <- c(rep(c("5%", "50%", "95%", "N_eff"), 2) )

R2.sown <- c(R2.sown.seed.only, NA, R2.sown.seed.disturb, NA)

# rbind together
cis.sown <- rbind(cis.sown, R2 = R2.sown)



# write credible intervals to table
write.table(cis.rich, sep = "\t", dec = ".", 
            row.names = TRUE, quote = FALSE,
            "doc/tables/delta-richness-resident-richness-cis.txt")

write.table(cis.sown, sep = "\t", dec = ".", 
            row.names = TRUE, quote = FALSE,
            "doc/tables/delta-sown-richness-resident-richness-cis.txt")



save(fit.rich.seed.only, fit.sown.seed.only,
     fit.rich.seed.disturb, fit.sown.seed.disturb, 
     file = "Rdata/03-1_brms-posteriors-resident-richness.Rdata")





################################################################################
### 4 plot observed and fitted delta richness                             ######
################################################################################




# plot observed and fitted delta richness, first for seeding only treatment


newrich <- seq(from = min(subset(metadata, treatment == "seed" & year == 2019)$log.resident.rich), 
               to = max(subset(metadata, treatment == "seed" & year == 2019)$log.resident.rich), 
               length.out = 200)
newdat <- data.frame(region = rep("ALB", 200),
                     log.resident.rich = newrich)


# first predict responses for the seeding only treatment
pred.rich.seed.only <- fitted(fit.rich.seed.only, newdata = newdat, 
                              probs = c(0.05, 0.50, 0.95))


# second, predict responses for the seeding and disturbance treatment
pred.rich.seed.disturb <- fitted(fit.rich.seed.disturb, newdata = newdat, 
                                 probs = c(0.05, 0.50, 0.95))


# add the mean here to plot on original scale


# cbind newdata to predicted repsonses and combine data frames
pred.rich.seed.only    <- cbind(newdat, pred.rich.seed.only)
pred.rich.seed.disturb <- cbind(newdat, pred.rich.seed.disturb)



p1 <- ggplot(data=pred.rich.seed.only, aes(x=log.resident.rich) ) +
  geom_ribbon(aes(ymin=Q5, ymax=Q95), linetype=0, alpha = 0.3) +
  geom_line(aes(y=Q50), size=0.5, linetype="solid", alpha = 0.8) +
  geom_point(data = subset(metadata, treatment == "seed" & year == 2019), 
             aes(y = delta.rich, x = log.resident.rich, colour = region)) +
  scale_x_continuous("Resident species richness [4 m²]" , breaks = log(c(10, 20, 40)), label = c(10, 20, 40) ) +
  scale_y_continuous(expression(paste("Predicted ", Delta, "richness [4 m²]")), 
                     limits = c(-10, 26)) +
  labs(title = "Seeding only") +  #specify both, because I specified both in first ggplot() call
  scale_colour_manual(values=c("royalblue4", "steelblue3", "brown3"), name = "Region" ) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(), 
                     legend.position="right", 
                     legend.title=element_text(size=9) ,
                     legend.text=element_text(size=8),
                     axis.text=element_text(size=8, colour="black"),
                     axis.title=element_text(size=9),
                     strip.background=element_blank() ) 
p1




p2 <- ggplot(data=pred.rich.seed.disturb, aes(x=log.resident.rich) ) +
  geom_ribbon(aes(ymin=Q5, ymax=Q95), linetype=0, alpha = 0.3) +
  geom_line(aes(y=Q50), size=0.5, linetype="solid", alpha = 0.8) +
  geom_point(data = subset(metadata, treatment == "seeddisturb" & year == 2019), 
             aes(y = delta.rich, x = log.resident.rich, colour = region)) +
  scale_x_continuous("Resident species richness [4 m²]" , breaks = log(c(10, 20, 40)), label = c(10, 20, 40) ) +
  scale_y_continuous("",
                     limits = c(-10, 26)) +
  labs(title = "Seeding and disturbance") +  
  scale_colour_manual(values=c("royalblue4", "steelblue3", "brown3"), name = "Region") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 10) ,
                     legend.text = element_text(size = 10),
                     axis.text = element_text(size = 10, colour = "black"),
                     axis.title = element_text(size = 10),
                     strip.background = element_blank()
                    ) 
p2

ggarrange(p1, p2, labels = c("(a)", "(b)"), 
          nrow=1, common.legend = TRUE, 
          legend = "bottom", 
          font.label = list(size = 15)
          )

ggsave("doc/figs/03-1_fitted_delta-richness-2019-resident-richness.pdf", plot = last_plot(), device = "pdf",
       width = 18, height = 10, units = "cm", dpi = 1200)

ggsave("doc/figs/03-1_fitted_delta-richness-2019-resident-richness.png", plot = last_plot(), device = "png",
       width = 18, height = 10, units = "cm", dpi = 1200)

