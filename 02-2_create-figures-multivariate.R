


library(ggplot2)
library(ggpubr)


################################################################################
### 0.1 Prepare scatter plots: delta richness vs. productivity           #######
################################################################################

# load the predicted data

load("Rdata/01-2_stan-delta-richness-pred.RData")


# add treatment variables and merge data frames


# focus on delta richness of sown species

predict.sown.seed.only.year5$pred.mean$treatment <- factor("(a) Seeding")
predict.sown.seed.only.year5$pred.regions.mean$treatment <- factor("(a) Seeding")
predict.sown.seed.only.year5$obs.delta$treatment <- factor("(a) Seeding")

predict.sown.seed.disturb.year5$pred.regions.mean$treatment <- factor("(b) Seeding and disturbance")
predict.sown.seed.disturb.year5$pred.mean$treatment <- factor("(b) Seeding and disturbance")
predict.sown.seed.disturb.year5$obs.delta$treatment <- factor("(b) Seeding and disturbance")




pred.mean.sown <- rbind(predict.sown.seed.only.year5$pred.mean,
                        predict.sown.seed.disturb.year5$pred.mean)

pred.regions.mean.sown <- rbind(predict.sown.seed.only.year5$pred.regions.mean,
                                predict.sown.seed.disturb.year5$pred.regions.mean)

obs.delta.sown <- rbind(predict.sown.seed.only.year5$obs.delta,
                        predict.sown.seed.disturb.year5$obs.delta)




# prepare annotation data frame to label facets inside the box

labels.sown <- data.frame(treatment = levels(pred.mean.sown$treatment),
                          label = levels(pred.mean.sown$treatment),
                          xat = rep(30, 2), yat = rep(29, 2)
)


################################################################################
### 1 Conditional effects of productivity (log.) on the                   ######
###   number of established species, relative to control                  ######
################################################################################

# I will plot the difference in (delta) absolute species richness
# as well as the difference in sown species richness (a subset of absolute richness)


# let's make it a function this year, this makes tweaking the graphics easier...
# plot title and axis labels can be changed later if needed

plot_delta <- function(data.mean = data.mean, 
                       data.regions = data.regions,
                       data.obs = data.obs){
p <- ggplot(data=data.mean, aes(x = productivity )) +
  geom_hline(yintercept=0 , size=0.3 , colour="gray20") + 
  geom_hline(yintercept=-10 , size=0.15, colour="gray80") +
  geom_hline(yintercept=10, size=0.15, colour="gray80") +
  geom_hline(yintercept=20, size=0.15, colour="gray80") +
  geom_ribbon(aes(ymin=X5., ymax=X95.), linetype=0, alpha = 0.3) +
  geom_line(aes(y = X50.), size=0.7, linetype="solid", alpha = 0.8) +
  geom_line(data = data.regions, aes(y = X50., colour = region), 
            size=0.3, linetype="solid", alpha = 0.8) +
  geom_point(data = data.obs, aes(y = delta, colour = region), 
             size = 0.5) +
  facet_wrap( . ~ treatment, nrow = 2) +
  scale_x_continuous(expression("Productivity (g m" ^-2*")"), 
                     breaks = c(100, 200, 300, 400), label = c(100, 200, 300, 400) ) +
  scale_y_continuous(limits = c(min(obs.delta.sown$delta), 30),
                     breaks = seq(-10, 20, by = 10)) +
  scale_colour_manual(values=c("royalblue4", "steelblue3", "brown3"), name = "Region") +
  ylab(expression(paste("Sown species ", Delta, "richness 5th year")) ) +
  theme_bw() + theme(panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank(), 
                     legend.position="bottom", 
                     legend.title=element_text(size=9) ,
                     legend.text=element_text(size=8),
                     axis.text=element_text(size=8, colour="black"),
                     axis.title=element_text(size=9),
                     strip.text=element_blank(),    #size of facet titles
                     strip.background=element_blank(),
                     plot.title = element_text(face = "bold", size = 12) )  
  }





# plots for sown species richness on the two seeding treatments
plot_delta(data.mean = pred.mean.sown, 
                        data.regions = pred.regions.mean.sown,
                        data.obs = obs.delta.sown)  +
  geom_text(data = labels.sown, aes(x = xat, y = yat, label = label),
            hjust = 0, fontface = "bold", size = 3.5) 



ggsave("doc/figs/01-3_delta-sown-richness-predict.pdf", plot = last_plot(), device = "pdf",
       width = 8, height = 16, units = "cm", dpi = 1200)

ggsave("doc/figs/01-3_delta-sown-richness-predict.png", plot = last_plot(), device = "png",
       width = 8, height = 16, units = "cm", dpi = 1200)



