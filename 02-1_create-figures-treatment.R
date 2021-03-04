





#################################################################################
### 1.1 Forest plots of predictors GLMM diversity models with   #################
###     three-way-interaction of seeding, disturbance and time  #################
#################################################################################



library(brms)
library(ggplot2)
library(ggpubr)


# load the models with posteriors and the data used to fit these models
load("Rdata/01-1_brms_posteriors_treatment_effects.Rdata") 



#first, get posterior density intervals for the species richness model
posterior.rich <- data.frame(fixef(treat.eff.rich, probs=c(0.05, 0.95, 0.5, 0.25, 0.75)))
posterior.rich$parameter <- factor(rownames(posterior.rich))
posterior.rich$parameter <- factor(posterior.rich$parameter, 
                                   levels=rev( c("Intercept", 
                                                 
                                                 "seed", 
                                                 "seed:year2016", 
                                                 "seed:year2017", 
                                                 "seed:year2018",
                                                 "seed:year2019", 
                                                 
                                                 "dist", 
                                                 "dist:year2016",
                                                 "dist:year2017",
                                                 "dist:year2018",
                                                 "dist:year2019",
                                                 
                                                 "dist:seed", 
                                                 "dist:seed:year2016",
                                                 "dist:seed:year2017",
                                                 "dist:seed:year2018",
                                                 "dist:seed:year2019",
                                                 
                                                 "year2016", "year2017", "year2018", "year2019"
                                   ) 
                                   ) ,
                                   
                                   labels=rev( c("Intercept", 
                                                 
                                                 "Seeding",   
                                                 "Seeding 2nd yr", 
                                                 "Seeding 3rd yr", 
                                                 "Seeding 4th yr", 
                                                 "Seeding 5th yr", 
                                                 
                                                 "Disturbance",
                                                 "Disturbance 2nd yr",
                                                 "Disturbance 3rd yr",
                                                 "Disturbance 4th yr",
                                                 "Disturbance 5th yr",
                                                 
                                                 "Seeding*\ndisturbance",
                                                 "Seeding*\ndisturbance 2nd yr",
                                                 "Seeding*\ndisturbance 3rd yr",
                                                 "Seeding*\ndisturbance 4th yr",
                                                 "Seeding*\ndisturbance 5th yr",
                                                 
                                                 "2nd year", "3rd year", "4th year", "5th year"
                                   ) 
                                   ) 
)



posterior.rich <- droplevels(posterior.rich[c(2,3,8:20) , ]) #remove intercept and years

posterior.rich$model <- as.factor(rep("Species richness", nrow(posterior.rich)))





# now get posteriors for the shannon model

posterior.spie <- data.frame(fixef(treat.eff.spie, probs=c(0.05, 0.95, 0.5, 0.25, 0.75)))
posterior.spie$parameter <- factor(rownames(posterior.spie))
posterior.spie$parameter <- factor(posterior.spie$parameter, 
                                   levels=rev( c("Intercept", 
                                                 
                                                 "seed", 
                                                 "seed:year2016", 
                                                 "seed:year2017", 
                                                 "seed:year2018",
                                                 "seed:year2019", 
                                                 
                                                 "dist", 
                                                 "dist:year2016",
                                                 "dist:year2017",
                                                 "dist:year2018",
                                                 "dist:year2019",
                                                 
                                                 "dist:seed", 
                                                 "dist:seed:year2016",
                                                 "dist:seed:year2017",
                                                 "dist:seed:year2018",
                                                 "dist:seed:year2019",
                                                 
                                                 "year2016", "year2017", "year2018", "year2019"
                                   ) 
                                   ) ,
                                   
                                   labels=rev( c("Intercept", 
                                                 
                                                 "Seeding",   
                                                 "Seeding 2nd yr", 
                                                 "Seeding 3rd yr", 
                                                 "Seeding 4th yr", 
                                                 "Seeding 5th yr", 
                                                 
                                                 "Disturbance",
                                                 "Disturbance 2nd yr",
                                                 "Disturbance 3rd yr",
                                                 "Disturbance 4th yr",
                                                 "Disturbance 5th yr",
                                                 
                                                 "Seeding*\ndisturbance",
                                                 "Seeding*\ndisturbance 2nd yr",
                                                 "Seeding*\ndisturbance 3rd yr",
                                                 "Seeding*\ndisturbance 4th yr",
                                                 "Seeding*\ndisturbance 5th yr",
                                                 
                                                 "2nd year", "3rd year", "4th year", "5th year"
                                   ) 
                                   ) 
)



posterior.spie <- droplevels(posterior.spie[c(2,3,8:20) , ]) #remove intercept and years

posterior.spie$model <- as.factor(rep("Spie", nrow(posterior.spie)))




#define labels for the treatments

labels_treat <- rev(expression(bold("Seeding"),   
                               "   * 2nd year", 
                               "   * 3rd year", 
                               "   * 4th year", 
                               "   * 5th year", 
                               
                               bold("Disturbance"),
                               "   * 2nd year", 
                               "   * 3rd year", 
                               "   * 4th year", 
                               "   * 5th year", 
                               
                               bold("\nSeeding\n* disturbance"),
                               "   * 2nd year", 
                               "   * 3rd year", 
                               "   * 4th year", 
                               "   * 5th year")
)



posterior.treat <- rbind(posterior.rich, posterior.spie)



# define background rectangles for the three treatments
rects <- data.frame(
  xmin = c(0, 10.5),
  xmax = c(5.5, 15.5),
  ymin = -Inf,
  ymax = Inf)


pd <- position_dodge(-0.4)


#plot the posteriors of both models together in one graph to make comparisons easy
params.plot <- ggplot(data=posterior.treat, aes(x=parameter, Q50, colour = model) ) + 
  # add background rectangles for the three treatments
  geom_rect(data=rects, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill="gray", alpha=0.3, inherit.aes = FALSE) +
  geom_hline(yintercept=0, size=0.25, linetype="solid", colour="black") +
  geom_linerange(aes(ymin=Q5, ymax=Q95), size=0.3, position=pd, 
                 show.legend = FALSE) + #plot 90% CI
  #geom_linerange(aes(ymin=Q25, ymax=Q75), size=0.9, position=pd, 
  #               show.legend = FALSE) + #plot 50% CI
  geom_point(aes(y=Q50), size=1.7, position=pd) + 
  scale_colour_manual(values=c("black", "gray50") ) + 
  scale_y_continuous("Effect sizes (log. scale)",
                     breaks = c(-0.2, 0, 0.2) ) +
  scale_x_discrete(limits = levels(posterior.rich$parameter), 
                   labels = labels_treat ) + 
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(), 
                     legend.title=element_blank(), 
                     legend.text=element_text(size=9) ,
                     legend.text.align = 0.5,
                     legend.position="bottom",
                     legend.direction="horizontal",
                     legend.background=element_rect(fill="white") , #alters background behind text and symbols
                     axis.title.y=element_blank() ,
                     axis.title.x=element_text(colour="black", size=9) ,
                     axis.ticks.y=element_blank() , # this removes the ticks from the y axis, independent of later flip...
                     axis.ticks.x=element_line(colour="black", size=0.5) ,
                     axis.text=element_text(size=9, colour="black"),
                     axis.text.y = element_text(hjust = 0),
                     panel.border=element_rect(colour="black", fill=NA) ) + 
  coord_flip()





################################################################################
### 1 Plot fitted responses of diversity to seeding, disturbance         #######
###   and their interaction over time                                    #######
################################################################################




################################################################################
### 1.1 create new data to predict                                      ########
################################################################################



# set empty data frame
newdata <- unique(metadata[,c("region", "dist", "seed", "year", 
                                  "region.year", "treatment")])
nrow(newdata)  # 4 treatments times 3 regions times 5 years -> 60 cases 



# reorder and rename treatments to get nice plot labels

newdata$treatment <- factor(newdata$treatment, 
                            levels = c("control", "disturb", "seed", "seeddisturb"), 
                            labels = c("Control", "(a) Disturbance", "(b) Seeding", 
                                     "(c) Seeding and\ndisturbance"))






################################################################################
### 1.2 Predict species richness on the response scale                  ########
################################################################################



# take region and year effects into account, but ignore grassland site effects

pred.rich <- posterior_epred(treat.eff.rich, newdata = newdata, 
                             scale = "response", 
                             re_formula = ~ (1|region.year), 
                             summary = FALSE
                             )

# cbind with newdata
pred.rich<- cbind(newdata, t(pred.rich))




# now substract the fitted values of the control treatment from the 
# fitted values of the three treatments, to have the control treatment as
# the baseline -> visualize the treatment effects relative to the control


# get only control predictions and prepare this data frame
ctrls.rich <- subset(pred.rich, treatment == "Control")

# repeat for the 3 other treatments
ctrls.rich <- ctrls.rich[rep(seq_len(nrow(ctrls.rich)), 3), ] 

#order by region, as treatments are ordered by region
ctrls.rich <- ctrls.rich[ order(ctrls.rich$region) , ] 

# and remove control treatments here
treatments.rich <- subset(pred.rich, treatment != "Control")

# substract control posteriors from treatment posteriors
delta.rich <- treatments.rich[ , -c(1:6) ] - ctrls.rich[ , -c(1:6) ]
delta.rich <- posterior_summary(t(delta.rich), probs = c(0.05, 0.95))

# combine with grouping variables
delta.rich <- cbind( treatments.rich[ , c(1:6) ] , delta.rich )

# rename region for plotting
levels(delta.rich$region) <- c("ALB", "HAI", "SCH")  

# year as numeric to draw lines between years in the plot
delta.rich$year <- as.numeric(as.character(delta.rich$year))





################################################################################
### 1.3 Predict S_pie, effective number of species                       ########
################################################################################


# take region and year effects into account, but ignore grassland site effects

pred.spie <- posterior_epred(treat.eff.spie, newdata = newdata, 
                             scale = "response", 
                             re_formula = ~ (1|region.year), 
                             summary = FALSE
)

# cbind with newdata
pred.spie<- cbind(newdata, t(pred.spie))




# now substract the fitted values of the control treatment from the 
# fitted values of the three treatments, to have the control treatment as
# the baseline -> visualize the treatment effects relative to the control


# get only control predictions and prepare this data frame
ctrls.spie <- subset(pred.spie, treatment == "Control")

# repeat for the 3 other treatments
ctrls.spie <- ctrls.spie[rep(seq_len(nrow(ctrls.spie)), 3), ] 

#order by region, as treatments are ordered by region
ctrls.spie <- ctrls.spie[ order(ctrls.spie$region) , ] 

# and remove control treatments here
treatments.spie <- subset(pred.spie, treatment != "Control")

# substract control posteriors from treatment posteriors
delta.spie <- treatments.spie[ , -c(1:6) ] - ctrls.spie[ , -c(1:6) ]
delta.spie <- posterior_summary(t(delta.spie), probs = c(0.05, 0.95))

# combine with grouping variables
delta.spie <- cbind( treatments.spie[ , c(1:6) ] , delta.spie )

# rename region for plotting
levels(delta.spie$region) <- c("ALB", "HAI", "SCH")  

# year as numeric to draw lines between years in the plot
delta.spie$year <- as.numeric(as.character(delta.spie$year))







# prepare annotation data frame to label facets inside the box

labels.rich <- data.frame(treatment = levels(droplevels(delta.rich$treatment)),
                               label = c("(b) Disturbance", "(c) Seeding", 
                                         "(d) Seeding and\n     disturbance"),
                               xat = rep(2015, 3), yat = rep(20, 3)
)

labels.spie <- data.frame(treatment = levels(droplevels(delta.spie$treatment)),
                                  label = c("(e) Disturbance", "(f) Seeding", 
                                            "(g) Seeding and\n      disturbance"),
                                  xat = rep(2015, 3), yat = rep(3.05, 3)
)






################################################################################
### 2. Create the plots                                                 ########
################################################################################



# define some space between predicted lines
pd <- position_dodge(0.4)


# plot each response individually
# ... would be great to set shared plot style and just update the data ... 


# Delta richness plot 

delta.rich.plot <- ggplot(data = delta.rich, aes(year, Estimate)) +
  geom_hline(yintercept = 0 , size = 0.5 , colour = "gray20") + 
  geom_hline(yintercept = 5 , size = 0.25, colour = "gray80") +
  geom_hline(yintercept = 10, size = 0.25, colour = "gray80") +
  geom_hline(yintercept = 15, size = 0.25, colour = "gray80") +
  geom_line(aes(x = year, y = Estimate, colour = region), 
            size = 0.375, linetype = "solid", position = pd, 
            show.legend = FALSE) +
  geom_linerange(aes(x = year, ymin = Q5, ymax = Q95, colour = region), 
                 size = 0.25, position = pd, 
                 show.legend = FALSE, alpha = 0.5 ) +
  geom_point(aes(x = year, y = Estimate, colour = region), 
             size = 1.5, position = pd) +
  scale_colour_manual(values = c("royalblue4", "steelblue3", "brown3") ,
                      name = "Region") + 
  facet_wrap(. ~ treatment ) + 
  scale_x_continuous("Study year", breaks = seq(2015, 2019, by = 1), 
                     label = seq(1,5) ) +
  scale_y_continuous(expression(atop(textstyle("Difference in species richness [4 m²]"),
                     textstyle("(Treatment - Control +/- 90% CI)") ) ),  
                     breaks = seq(0, 15, by = 5), limits =  c(-2, 22) ) +
  geom_text(data = labels.rich, aes(x = xat, y = yat, label = label),
            hjust = 0, fontface = "bold", size = 3.5) +
  # add title here instead of with ggarrange 'labels': can add expression with subscript
  labs(title = expression(bold( textstyle("Species richness") ) ) ) +     
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 9) ,
                     legend.text = element_text(size = 9),
                     axis.text = element_text(size = 9, colour = "black"),
                     axis.title = element_text(size = 9),
                     strip.text = element_blank(), #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 12)
                     )
                   


# Delta Spie plot

delta.spie.plot <- ggplot(data = delta.spie, aes(year, Estimate)) +
  geom_hline(yintercept = -1 , size = 0.25 , colour = "gray80") + 
  geom_hline(yintercept = 0 , size = 0.5 , colour = "gray20") + 
  geom_hline(yintercept = 1 , size = 0.25, colour = "gray80") +
  geom_hline(yintercept = 2, size = 0.25, colour = "gray80") +
  geom_line(aes(x = year, y = Estimate, colour = region), 
            size = 0.375, linetype = "solid", position = pd, 
            show.legend = FALSE) +
  geom_linerange(aes(x = year, ymin = Q5, ymax = Q95, colour = region), 
                 size = 0.25, position = pd, 
                 show.legend = FALSE, alpha = 0.5 ) +
  geom_point(aes(x = year, y = Estimate, colour = region), 
             size = 1.5, position = pd) +
  scale_colour_manual(values = c("royalblue4", "steelblue3", "brown3") ,
                      name = "Region") + 
  facet_wrap(. ~ treatment ) + 
  scale_x_continuous("Study year", breaks = seq(2015, 2019, by = 1), 
                     label = seq(1,5) ) +
  scale_y_continuous(expression(atop(textstyle("Difference in") ~ S[PIE],
                                paste("(Treatment - Control +/- 90% CI)")) ), 
                     breaks = seq(-1, 2, by = 1), limits = c(-1, 3.3)  ) +
  geom_text(data = labels.spie, aes(x = xat, y = yat, label = label),
            hjust = 0, fontface = "bold", size = 3.5) +
  # add title here instead of with ggarrange 'labels': can add expression with subscript
  labs(title = expression( bold(S[PIE]) ) ) +     
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     legend.position = "right", 
                     legend.title = element_text(size = 9) ,
                     legend.text = element_text(size = 9),
                     axis.text = element_text(size = 9, colour = "black"),
                     axis.title = element_text(size = 9),
                     strip.text = element_blank(), #size of facet titles
                     strip.background = element_blank(),
                     plot.title = element_text(face = "bold", size = 12)
                     ) 





# combine the plots
ggarrange(delta.rich.plot, delta.spie.plot, 
          nrow = 2, 
          common.legend = TRUE, legend = "bottom")




ggarrange(
  ggarrange(NULL, params.plot,
            nrow = 2,
            heights = c(0.05, 0.95) ),
  NULL, ggarrange(
    delta.rich.plot, delta.spie.plot, 
    nrow = 2, 
    common.legend = TRUE, legend = "bottom"),
  ncol = 3, labels = c("(a)        Parameters", "", ""),
  widths = c(0.32, 0.03, 0.65),
  hjust = -0.2,
  vjust = 1.5,
  common.legend = FALSE,
  font.label = list(size = 12)
)


ggsave("doc/figs/01-1_treatment-effects.pdf", plot = last_plot(), device = "pdf",
       width = 18, height = 15, units = "cm", dpi = 1200)

ggsave("doc/figs/01-1_treatment-effects.png", plot = last_plot(), device = "png",
       width = 18, height = 15, units = "cm", dpi = 1200)




################################################################################
### 3 Plot posterior-predictive density overlay                         ########
################################################################################


# use the 'pp_check' function from the bayesplot package to plot observed y vs.
# y_pred density

ppc.rich <- pp_check(treat.eff.rich, nsamples=50) +
  labs(title = expression(  Response ~ italic(y): Species ~ richness ) ) + 
  theme(plot.title = element_text(hjust=0.3, size = 9),
        plot.subtitle = element_text(hjust=0.2, size = 9) ) 

ppc.spie <- pp_check(treat.eff.spie, nsamples=50) +
  labs(title = expression(  Response ~ italic(y): bolditalic(S[PIE] ) ) ) + 
  theme(plot.title = element_text(hjust=0.3, size = 9),
        plot.subtitle = element_text(hjust=0.2, size = 9) ) 




ggarrange(ppc.rich, ppc.spie, 
          labels = c("(a)", "(b)" ), 
          nrow = 1, 
          hjust = c(0, 0),
          vjust = c(1.5, 1.5),
          font.label = list(size = 10),
          common.legend = TRUE, legend = "bottom")



ggsave("doc/figs/S1_pp-check-treatment-effects.pdf", plot = last_plot(), device = "pdf",
       width = 14, height = 8, units = "cm", dpi = 1200)

ggsave("doc/figs/S1_pp-check-treatment-effects.png", plot = last_plot(), device = "png",
       width = 14, height = 8, units = "cm", dpi = 1200)



