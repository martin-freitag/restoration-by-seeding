


library(ggplot2)
library(ggpubr)
library(GGally) # for multivariate pairs plots
library(rstan)
library(brms)




################################################################################
#   0.1 Load and prepare the data                                      #########
################################################################################




# we are going to use the 'metadata' data, which contains info about  
# the experimental design and metadata like grazing and mowing intensities, 
# but also productivity (aboveground biomass) sampling data
metadata <- read.table("data/metadata.txt", header = TRUE,
                       sep = "\t", dec = ".", stringsAsFactors = TRUE)


# reorder and rename treatments to get nice plot labels

metadata$treatment <- factor(metadata$treatment, 
                             levels = c("control", "seed", "seeddisturb", "disturb"), 
                             labels = c("Control", "Seeding",
                                        "Seeding and\ndisturbance", "Disturbance" ))



# here, we source a script to prepare the 'newly.sown.long' data frame, 
# which contains species cover data of only newly sown species on the 
# two seeding treatments in 2019. Species that were already present on control 
# and disturbance treatments were removed.

source("R/00-1_sown-and-established-species-data-wrangling.R") 


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
             & seeddensity$region == "AEG" ,
             ]$germinationrate <- mean.hyp.perf


seeddensity$live.seeddensity <- seeddensity$seeddensity.per.m2 * (seeddensity$germinationrate / 100)


################################################################################
#   1 Plot percent cover of litter, vegetation and bare soil           #########
################################################################################



# define some space between predicted lines
pd <- position_dodge(0.75)


# plot each response individually
# ... would be great to set shared plot style and just update the data... too late. 


# litter plot 

litter.plot <- ggplot(data = metadata, aes(x = as.factor(year), y = cover.litter.pc,
                                           fill = treatment)) +
  geom_boxplot(position = pd, show.legend = TRUE,
               outlier.size = 0.15, width = 0.7, size = 0.1) +
  scale_colour_viridis_d() + 
  facet_wrap(. ~ region ) + 
  scale_x_discrete("Study year", label = seq(1,5) ) +
  scale_y_continuous("Litter cover [%]") +
  labs(fill="Treatment") + 
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 10) ,
                     legend.text = element_text(size = 10),
                     axis.text = element_text(size = 10, colour = "black"),
                     axis.title = element_text(size = 10),
                     strip.text = element_text(size = 12), #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 15)
  )



herbs.plot <- ggplot(data = metadata, aes(x = as.factor(year), y = cover.herbs.pc,
                                          fill = treatment)) +
  geom_boxplot(position = pd, show.legend = TRUE,
               outlier.size = 0.15, width = 0.75, size = 0.1) +
  scale_colour_viridis_d() + 
  facet_wrap(. ~ region ) + 
  scale_x_discrete("Study year", label = seq(1,5) ) +
  scale_y_continuous("Vegetation cover [%]") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 10) ,
                     legend.text = element_text(size = 10),
                     axis.text = element_text(size = 10, colour = "black"),
                     axis.title = element_text(size = 10),
                     strip.text = element_text(size = 12), #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 15)
  )




soil.plot <- ggplot(data = metadata, aes(x = as.factor(year), y = cover.bare.soil.pc,
                                         fill = treatment)) +
  geom_boxplot(position = pd, show.legend = TRUE,
               outlier.size = 0.15, width = 0.75, size = 0.1) +
  scale_colour_viridis_d() + 
  facet_wrap(. ~ region ) + 
  scale_x_discrete("Study year", label = seq(1,5) ) +
  scale_y_continuous("Bare soil cover [%]") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 10) ,
                     legend.text = element_text(size = 10),
                     axis.text = element_text(size = 10, colour = "black"),
                     axis.title = element_text(size = 10),
                     strip.text = element_text(size = 12), #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 15)
  )




# combine the plots
ggarrange(herbs.plot, litter.plot, soil.plot,
          labels = c("(a) Vegetation", "(b) Litter", "(c) Bare soil" ), 
          nrow = 3, 
          hjust = c(-0.01, -0.05, -0.01),
          vjust = c(1.4, 1.4, 1.4),
          font.label = list(size = 12),
          common.legend = TRUE, legend = "bottom")



ggsave("doc/figs/03_litter-veg-soil-cover.pdf", plot = last_plot(), device = "pdf",
       width = 16, height = 16, units = "cm", dpi = 1200)

ggsave("doc/figs/03_litter-veg-soil-cover.png", plot = last_plot(), device = "png",
       width = 16, height = 16, units = "cm", dpi = 1200)



################################################################################
#   2 Plot fertilization, grazing, mowing                              #########
#     and productivity predictor distributions                         #########
################################################################################

# Refers to analysis part 2, where I analysed land-use effects on 
# the increase in sown species richness (delta richness relative to control)
# and used grazing, fertilization and mowing as predictors


# first get unique grassland land use information (i.e. get from 1448 obs to 73 grasslands)

grasslands <- unique(subset(metadata, treatment == "Control",
                            select = c("plot", "region", "Fert", "Graz", "Mow"))
)



# merge the latent productivity medians (from trait analyses) to the data
# productivity was modelled exactly the same in both analyses


# load the trait model posteriors
load("Rdata/01-3_stan-posteriors-traits.Rdata")

# get latent productivity data
latentbm <- as.data.frame(stan.traits.seed.only, pars = "realbm")

# calculate median and log-transform productivity for grasslands (N = 73)
# log-transform, as it was used as such as predictor in the analyses
log.median.bm <- apply(latentbm, MARGIN = 2, function(x) log( median(x) ) )


# productivity data from the trait establishment models
# are sorted by grassland ID ('plot', N = 73), therefore this data can be merged
# with the ordered grassland plot data
grasslands$log.median.bm <- log.median.bm


# transform land-use components and productivity as done in the analyses:
#   - Fertilization intensity in kg N ha^-1                       -> log(x + 1)
#   - Grazing intensity in livestock unit grazing days * ha^-1    -> log(x + 1)
#   - Mowing intensity in mean cuts per year                      ->     x
#   - Productivity modelled as aboveground biomss g * m^-2        -< log(x)

grasslands$log.Fert         <- log(grasslands$Fert + 1)
grasslands$log.Graz         <- log(grasslands$Graz + 1)

# plot distributions, pairs-plots and correlations of predictors

pl <- ggpairs(grasslands, mapping = aes(color = region, alpha = 0.5), 
              lower = list(continuous = wrap("points", alpha = 1, size=0.15) ), 
              columns = c("log.Fert", "log.Graz", "Mow", "log.median.bm"),
              columnLabels = c("Fertilization (log).", "Grazing (log).", 
                               "Mowing", "Productivity (log.)") ) +
  theme(strip.text = element_text(size = 8) ) #size of facet titles

pl



ggsave("doc/figs/03_land-use-productivity-distribution.pdf", 
       plot = last_plot(), device = "pdf",
       width = 16, height = 16, units = "cm", dpi = 1200)

ggsave("doc/figs/03_land-use-productivity-distribution.png", 
       plot = last_plot(), device = "png",
       width = 16, height = 16, units = "cm", dpi = 1200)



################################################################################
#   3 Plot distribution of productivity and trait distributions        #########
################################################################################



# merge traits into the newly sown species data frame
newly.sown.plot <- merge(newly.sown.long, traits, 
                         by = "species", all.x = TRUE)


# merge germination rates into the newly sown species data frame
newly.sown.plot <- merge(newly.sown.plot, seeddensity, 
                         by = c("region", "species"), all.x = TRUE)


# nicer labels
newly.sown.plot$region <- factor(newly.sown.plot$region, 
                                 levels = c("AEG", "HEG", "SEG"), 
                                 labels = c("ALB", "HAI", "SCH"))


# log-transform height and seed mass as done in the analysis
newly.sown.plot$log.height   <- log(newly.sown.plot$height)
newly.sown.plot$log.seedMass <- log(newly.sown.plot$seedMass)



# merge mean grasslands productivity (from section above) as well 
newly.sown.plot <- merge(newly.sown.plot, grasslands[ , c("plot", "log.median.bm")], 
                         by = "plot", all.x = TRUE)

# and select one of the two treatments, as the data is the same
newly.sown.plot <- subset(newly.sown.plot, treatment == "seeddisturb",
                          select = c("region", "log.height", 
                                     "log.seedMass", "SLA", 
                                     "live.seeddensity", "log.median.bm"))



pt <- ggpairs(newly.sown.plot, mapping = aes(color = region, alpha = 0.5), 
              lower = list(continuous = wrap("points",    size=0.1) ),
              columns = c("log.height", "log.seedMass", "SLA", 
                          "live.seeddensity", "log.median.bm"),
              columnLabels = c("Height (log).", "Seed mass (log).", 
                               "Specific leaf area", "Live seeding density", 
                               "Productivity (log.)")) +
  theme(strip.text = element_text(size = 8) ) #size of facet titles


pt


ggsave("doc/figs/03_traits-predictor-distribution.pdf", 
       plot = last_plot(), device = "pdf",
       width = 16, height = 16, units = "cm", dpi = 1200)

ggsave("doc/figs/03_traits-predictor-distribution.png", 
       plot = last_plot(), device = "png",
       width = 16, height = 16, units = "cm", dpi = 1200)



################################################################################
#   3 Plot bare soil cover along grazing and mowing intensities        #########
################################################################################


# make use of the 'grasslands' data frame used before
grasslands$cover.bare.soil.pc <- aggregate(cover.bare.soil.pc ~ plot, 
                                           data = metadata[metadata$treatment == "Control" , ],
                                           FUN = mean, na.action = na.omit)[,2]


# plot bare soil cover against grazing, mowing and productivity

soil.grazing <- ggplot(data = grasslands, aes(x = log.Graz, y = cover.bare.soil.pc,
                                              colour = region)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(show.legend = TRUE, size = 0.75) +
  scale_colour_manual(values = c("royalblue4", "steelblue3", "brown3") ,
                      name = "Region") + 
  scale_x_continuous("Grazing intensity (log.)" ) +
  scale_y_continuous("Bare soil cover [%]") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 10) ,
                     legend.text = element_text(size = 10),
                     axis.text = element_text(size = 10, colour = "black"),
                     axis.title = element_text(size = 10),
                     strip.text = element_text(size = 12), #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 15)
  )





soil.mowing <- ggplot(data = grasslands, aes(x = Mow, y = cover.bare.soil.pc,
                                             colour = region)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(show.legend = TRUE, size = 0.75) +
  scale_colour_manual(values = c("royalblue4", "steelblue3", "brown3") ,
                      name = "Region") + 
  scale_x_continuous("Mowing frequency" ) +
  scale_y_continuous("") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 10) ,
                     legend.text = element_text(size = 10),
                     axis.text = element_text(size = 10, colour = "black"),
                     axis.title = element_text(size = 10),
                     strip.text = element_text(size = 12), #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 15)
  )





# combine the plots
ggarrange(soil.grazing, soil.mowing,
          labels = c("(a)", "(b)"), 
          nrow = 1, ncol = 2, 
          hjust = c(-0.01, -0.05),
          vjust = c(1.4, 1.4),
          font.label = list(size = 15),
          common.legend = TRUE, legend = "bottom")


ggsave("doc/figs/03_soil-cover-grazing-mowing.pdf", 
       plot = last_plot(), device = "pdf",
       width = 12, height = 7, units = "cm", dpi = 1200)

ggsave("doc/figs/03_soil-cover-grazing-mowing.png", 
       plot = last_plot(), device = "png",
       width = 12, height = 7, units = "cm", dpi = 1200)

