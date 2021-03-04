



library(ggplot2)
library(ggpubr)
library(rstan)
library(dplyr)
library(purrr)


# some functions (credible intervals from Stan fit objects, ggplot inset function, ...)
source("R/helpers.R")


################################################################################
### 1 Conditional effects of functional traits along the Productivity   ########
###   gradient, including interactions with Productivity (log.)         ########
################################################################################


# we will predict establishment for low and high trait values along
# the productivity gradient, for three traits (i.e. three facets) in the two
# seeding treatments
#
# In each facet, we inset a small forest plot of the trait and
# trait*productivity parameter estimates, to see the effect sizes and 
# predictions at the same time
 



################################################################################
### 1.1 get posterior density intervals of parameter estimates          ########
###   for the two models                                                ########
################################################################################


# model posteriors
load("Rdata/01-3_stan-posteriors-traits.Rdata")


# I have saved the trait predictions in a separate RData file. Load here
# both treatments. File contains predicted response-scale establishment 
# probabilities for mean +- SD trait values along the productivity gradient. 

load("Rdata/01-3_stan-traits-pred.Rdata")


# extract posteriors of parameter estimates

# the order of traits was defined in the 'X' predictor matrix of 
# the trait analyses, plus productivity and trait*productivity parameters
# added in the .stan file
par_names <-  c("Height", "Seed mass", "SLA", "Live seeding density",
                "Productivity", "Height\n*Productivity",
                "Seed mass\n*Productivity", "SLA\n*Productivity")


probs.traits.seed.only <- get_ci(stan.traits.seed.only, pars = "b",
                                 names = par_names, 
                                 probs = c(0.05, 0.25, 0.50, 0.75, 0.95))

probs.traits.seed.disturb <- get_ci(stan.traits.seed.disturb, pars = "b",
                                 names = par_names, 
                                 probs = c(0.05, 0.25, 0.50, 0.75, 0.95))

# set specific order of levels for the inset plots
probs.traits.seed.only$parameter <- factor(probs.traits.seed.only$parameter, 
                                           levels = c("Productivity", 
                                                      "Height\n*Productivity",
                                                      "Height",
                                                      "Seed mass\n*Productivity",
                                                      "Seed mass",
                                                      "SLA\n*Productivity",
                                                      "SLA", 
                                                      "Live seeding density") )

probs.traits.seed.disturb$parameter <- factor(probs.traits.seed.disturb$parameter, 
                                           levels = c("Productivity", 
                                                      "Height\n*Productivity",
                                                      "Height",
                                                      "Seed mass\n*Productivity",
                                                      "Seed mass",
                                                      "SLA\n*Productivity",
                                                      "SLA", 
                                                      "Live seeding density") )


probs.traits.seed.only
probs.traits.seed.disturb

# get parameter estimates for each trait, with productivity estimate
probs.seed.only.height <- probs.traits.seed.only[ c(5, 6, 1) , ]
probs.seed.only.seed   <- probs.traits.seed.only[ c(5, 7, 2) , ]
probs.seed.only.SLA    <- probs.traits.seed.only[ c(5, 8, 3) , ]



probs.seed.disturb.height <- probs.traits.seed.disturb[ c(5, 6, 1) , ]
probs.seed.disturb.seed   <- probs.traits.seed.disturb[ c(5, 7, 2) , ]
probs.seed.disturb.SLA    <- probs.traits.seed.disturb[ c(5, 8, 3) , ]





################################################################################
### 1.2 Prepare the inset forest plots of parameters                    ########
################################################################################


# code inspired by https://clarewest.github.io/blog-posts/ggplotInset.html

get_inset <- function(data, ymin = NULL, ymax = NULL) {
p <- ggplot(data = data,
            aes(x = parameter, y = p0.5 ) ) +
  geom_hline(yintercept = 0, size = 0.2, linetype = "solid", colour = "black") +
  geom_linerange(aes(ymin =  p0.05, ymax = p0.95), size = 0.15,
                 show.legend = FALSE) +                 #plot 90% CI
  geom_linerange(aes(ymin = p0.25, ymax = p0.75), size = 0.4,
                 show.legend = FALSE) +                 #plot 50% CI
  geom_point(aes(y = p0.5), size = 0.5) + 
  scale_colour_manual(values = c("black", "gray50") ) + 
  scale_y_continuous("", breaks = c(-3, 0, 3),
                     limits = c(ymin, ymax) ) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_blank() ,
                     axis.title.x = element_blank() ,
                     axis.ticks.y = element_blank() , # this removes the ticks from the y axis, independent of later flip...
                     axis.ticks.x = element_line(colour = "black", size = 0.2) ,
                     axis.text   = element_text(size = 7, colour = "black"),
                     axis.text.y = element_text(hjust = 0),
                     panel.background = element_rect(fill = "transparent"),
                     plot.background = element_rect(fill = "transparent", color = NA), 
                     panel.border = element_rect(colour = "black", fill = NA, size = 0.2) ) + 
  coord_flip()
return(p)
}





# make individual inset plots
# yes, this is lengthy and I could have written a function... sorry :)


# seeding treatment 

# specify position in trait-productivity prediction plots

xmin <- log(60)
xmax <- log(600)
ymin <- 0.11
ymax <- 0.32

inset_seed.only.height <- get_inset(probs.seed.only.height, 
                                    ymin = min(probs.traits.seed.only$p0.05),  
                                    ymax = max(probs.traits.seed.only$p0.95))
inset_seed.only.seed   <- get_inset(probs.seed.only.seed, 
                                    ymin = min(probs.traits.seed.only$p0.05),  
                                    ymax = max(probs.traits.seed.only$p0.95))
inset_seed.only.SLA    <- get_inset(probs.seed.only.SLA, 
                                    ymin = min(probs.traits.seed.only$p0.05),  
                                    ymax = max(probs.traits.seed.only$p0.95))


# data = ...[[1]][1, ] is a row identifier: height data starts at 1 with length 2*200,
# seed mass data starts at 401 and so on
inset_seed.only.height <- annotation_custom2(ggplotGrob(inset_seed.only.height),
                            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                            data = pred.traits.seed.only[[1]][ 1 , ])

inset_seed.only.seed   <- annotation_custom2(ggplotGrob(inset_seed.only.seed), 
                            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                            data = pred.traits.seed.only[[1]][ 401 , ])

inset_seed.only.SLA    <- annotation_custom2(ggplotGrob(inset_seed.only.SLA), 
                            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                            data = pred.traits.seed.only[[1]][ 801 , ])




# seeding and disturbance treatment

# specify position in trait-productivity prediction plots

ymin <- 0.22
ymax <- 0.65


inset_seed.disturb.height <- get_inset(probs.seed.disturb.height, 
                                    ymin = min(probs.traits.seed.disturb$p0.05),  
                                    ymax = max(probs.traits.seed.disturb$p0.95))
inset_seed.disturb.seed   <- get_inset(probs.seed.disturb.seed, 
                                    ymin = min(probs.traits.seed.disturb$p0.05),  
                                    ymax = max(probs.traits.seed.disturb$p0.95))
inset_seed.disturb.SLA    <- get_inset(probs.seed.disturb.SLA, 
                                    ymin = min(probs.traits.seed.disturb$p0.05),  
                                    ymax = max(probs.traits.seed.disturb$p0.95))


# data = ...[[1]][1, ] is a row identifier: height data starts at 1 with length 2*200,
# seed mass data starts at 401 and so on
inset_seed.disturb.height <- annotation_custom2(ggplotGrob(inset_seed.disturb.height),
                                             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                             data = pred.traits.seed.disturb[[1]][ 1 , ])

inset_seed.disturb.seed   <- annotation_custom2(ggplotGrob(inset_seed.disturb.seed), 
                                             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                             data = pred.traits.seed.disturb[[1]][ 401 , ])

inset_seed.disturb.SLA    <- annotation_custom2(ggplotGrob(inset_seed.disturb.SLA), 
                                             xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                             data = pred.traits.seed.disturb[[1]][ 801 , ])




# prepare annotation data frame to label facets inside the box

labels.seed.only <- data.frame(trait = levels(pred.traits.seed.only[[1]]$trait),
                               label = c("(a) Height (log.)", "(b) Seed mass (log.)",
                                         "(c) Specific leaf area"),
                               xat = rep(log(25), 3), yat = rep(0.35, 3)
                               )

labels.seed.disturb <- data.frame(trait = levels(pred.traits.seed.only[[1]]$trait),
                                  label = c("(d) Height (log.)", "(e) Seed mass (log.)",
                                            "(f) Specific leaf area"),
                                  xat = rep(log(25), 3), yat = rep(0.7, 3)
                                  )


################################################################################
### 1.3 Plot the establishment predictions for trait-productivity       ########
###     interactions                                                    ########
################################################################################



# seeding only treatment

p.seed.only <- ggplot(data = pred.traits.seed.only[[1]], 
             aes(x = log(productivity)) ) +   
  facet_grid(. ~  trait) +
  geom_ribbon(aes(ymin = X5., ymax = X95., colour = value, fill = value), 
              linetype = 0, alpha = 0.3) +
  geom_line(aes(y = X50., colour = value), 
            size = 0.5, linetype = "solid", alpha = 0.8) +
  scale_x_continuous(expression("Productivity (g m" ^-2*")"), 
                     breaks = log(c(25, 100, 400)), label = c(25, 100, 400) ) +
  scale_y_continuous("Probability of establishment", limits = c(0, 0.37) ) +
  labs(title = "Seeding", colour = "Trait value", 
       fill = "Trait value") +  #specify both, because I specified both in first ggplot() call
  geom_text(data = labels.seed.only, aes(x = xat, y = yat, label = label),
            hjust = 0.1, fontface = "bold", size = 3.5) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 9) ,
                     legend.text = element_text(size = 9),
                     axis.text = element_text(size = 9, colour = "black"),
                     axis.title = element_text(size = 9),
                     strip.text = element_blank(),    #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 12)
                     ) 

# and add the insets
p.seed.only.inset <- p.seed.only + inset_seed.only.height + 
  inset_seed.only.seed + inset_seed.only.SLA




# seeding and disturbance treatment


p.seed.disturb <- ggplot(data = pred.traits.seed.disturb[[1]], 
             aes(x = log(productivity)) ) +   
  facet_grid(. ~  trait) +
  geom_ribbon(aes(ymin = X5., ymax = X95., colour = value, fill = value), 
              linetype = 0, alpha = 0.3) +
  geom_line(aes(y = X50., colour = value), 
            size = 0.5, linetype = "solid", alpha = 0.8) +
  scale_x_continuous(expression("Productivity (g m" ^-2*")"), 
                     breaks = log(c(25, 100, 400)), label = c(25, 100, 400) ) +
  scale_y_continuous("Probability of establishment", limits = c(0, 0.7)) +
  labs(title = "Seeding and disturbance", colour = "Trait value", 
       fill = "Trait value") +  #specify both, beause I specified both in first ggplot() call
  geom_text(data = labels.seed.disturb, aes(x = xat, y = yat, label = label),
            hjust = 0.1, fontface = "bold", size = 3.5) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 9) ,
                     legend.text = element_text(size = 9),
                     axis.text = element_text(size = 9, colour = "black"),
                     axis.title = element_text(size = 9),
                     strip.text = element_blank(),    #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 12) 
                     ) 

# and add the insets
p.seed.disturb.inset <- p.seed.disturb + inset_seed.disturb.height + 
  inset_seed.disturb.seed + inset_seed.disturb.SLA






# now combine both plots into one figure wit ggarrange from ggpubr
ggarrange(p.seed.only.inset, p.seed.disturb.inset, 
#          labels = c("(a)", "(b)"), 
          nrow = 2, 
          hjust = c(-0.5, -0.5),
          vjust = c(1.5, 1.5),
          font.label = list(size = 12),
          common.legend = TRUE, legend  = "bottom")



ggsave("doc/figs/01-4_traits-predict.pdf", plot = last_plot(), device = "pdf",
       width = 18, height = 15, units = "cm", dpi = 1200)

ggsave("doc/figs/01-4_traits-predict.png", plot = last_plot(), device = "png",
       width = 18, height = 15, units = "cm", dpi = 1200)





################################################################################
### 2 Plot posterior-predictive density overlay                         ########
################################################################################


# check posterior predictive distributions: did the model behave well, does it fit the data?
# extract posterior predictive samples

# use function from the bayesplot package to do so



# both treatments
ppc.traits.seed.only <- ppc_dens_overlay_stanfit(stan.traits.seed.only, pars = "post", 
                                                 y = standata.traits.seed.only$pres) +
  labs(title = expression(  Response ~ italic(y): Establishment ~ success ~ year ~ five),
       subtitle = expression (Seeding ~ only ) )

ppc.traits.seed.disturb <- ppc_dens_overlay_stanfit(stan.traits.seed.disturb, pars = "post", 
                                                    y = standata.traits.seed.disturb$pres) +
  labs(title = expression(  Response ~ italic(y): Establishment ~ success ~ year ~ five),
       subtitle = expression (Seeding ~ and ~ disturbance ) )


# combine plots
ggarrange(ppc.traits.seed.only, ppc.traits.seed.disturb, 
          labels = c("(a)", "(b)" ), 
          nrow = 1, 
          hjust = c(0, 0),
          vjust = c(1.5, 1.5),
          font.label = list(size = 10),
          common.legend = TRUE, legend = "bottom")



ggsave("doc/figs/S1_pp-check-traits-establishment.pdf", 
       plot = last_plot(), device = "pdf",
       width = 14, height = 8, units = "cm", dpi = 1200)

ggsave("doc/figs/S1_pp-check-traits-establishment.png", 
       plot = last_plot(), device = "png",
       width = 14, height = 8, units = "cm", dpi = 1200)





