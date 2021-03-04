

################################################################################
### Prepare data frames with only the sown species for the analysis of   #######
### 1. how productivity and land-use affect species establishment        #######
###    (analysis '01-2_delta-richness...R')                              #######
### 2. how plant functional traits affect establishment                  #######
###    (analysis '01-3_traits-binomial.R')                               #######
################################################################################


# load the species cover data from our vegetation surveys, long format
species.long <- read.table("data/species-long.txt", header = TRUE,
                           sep = "\t", dec = ".", stringsAsFactors = TRUE)

# the species.long data frame contains species cover estimates from our
# vegetation sampling campaigns. Absences in the dataset are 'true' zeros,
# i.e. species absences. There are no NAs in this data.


# convert species data to wide format and bck to long format 
# to produce the 'true' zeros / absences
species.wide <- tidyr::pivot_wider(species.long, 
                                   names_from = species, 
                                   values_from = cover,
                                   values_fill = 0)

# order columns by species, bu keep first column where it is
species.wide <- species.wide[ , c(1, 1 + order(colnames(species.wide[-1] ) ) ) ]


# and convert back to long format to go on
species.long <- tidyr::pivot_longer(species.wide, 
                                    cols = !plotIdent, # convert all but this column
                                    names_to = "species", 
                                    values_to = "cover")

#  make species column a vector again
species.long$species <- as.factor(species.long$species)

# get metadata from identifier 'plotIdent'
species.long$plot      <- as.factor(substr(species.long$plotIdent, 1,5))
species.long$region    <- as.factor(substr(species.long$plotIdent, 1,3))
species.long$year      <- as.factor(substr(species.long$plotIdent, 10,13))
species.long$treatment <- as.factor(substr(species.long$plotIdent, 7,8))


species.long$treatment <- plyr::revalue(species.long$treatment,
                                        c("S1" = "control", 
                                          "S2" = "seed", 
                                          "S3" = "seeddisturb", 
                                          "S4" = "disturb")
                                        )




# import list of sown species. This is a table with three columns,
# one column for each region

sown.species <- read.table("data/sown-species.txt", 
                           header = T, na.strings = "NA")

sownAEG <- sown.species[ , "AEG" ] # this region has no NAs in the sown species
sownHEG <- sown.species[ - which(is.na(sown.species[ , "HEG" ])) , "HEG" ]
sownSEG <- sown.species[ - which(is.na(sown.species[ , "SEG" ])) , "SEG" ]




### region by region, subset sown species --------------------------------------


longAEG <- subset(species.long, region=="AEG")
longAEG <- subset(longAEG, species %in% sownAEG)

longHEG <- subset(species.long, region=="HEG")
longHEG <- subset(longHEG, species %in% sownHEG)

longSEG <- subset(species.long, region=="SEG")
longSEG <- subset(longSEG, species %in% sownSEG)


# combine the data frames of sown species only for later use
sown.long <- rbind(longAEG, longHEG, longSEG)



################################################################################
### to be sure that sown species were not already present before, we remove ####
### observations where a species occurred on the control or disturbance     ####
### treatments during the five years of study                               ####
################################################################################

# we chose the year 2019 to analyse the establishment of sown species 
# depending on their traits and trait-environment interactions
# not so easy to get this for all years, because in 2018 we have missing 
# observations in the Alb and Hainich regions, this function would not work

remove_present <- function(long = long){    # species data frame in long format
    
    # get the species that were present
    present <- subset(long, treatment == "control" | treatment == "disturb")
    present <- droplevels(present)
    # aggregate results in data frame ordered by plot, then by species
    present <- aggregate(cover ~ species + plot, data = present,
                            FUN = sum, na.action = na.omit)
    #repeat for the two seeding treatments in the experiment
    present <- rbind(present, present)
    
    # and order by plot, as the 'long' data frame is ordered by plot
    present <- present[order(present$plot) , ]
    
    
    # get the observed species covers on seeding treatments
    long <- subset(long, (treatment == "seed"  |  treatment == "seeddisturb") &
                         year == 2019 )
    long <- droplevels(long)  # drop control and disturbance treatment levels

    if( identical(long$plot, present$plot) && identical(long$species, present$species) ) {
    # remove the recordings of species that were already present
    long <- long[ present$cover == 0 , ]
    } else {
        long <- NULL
    }
}


# now subset to species that were sown but not present
newly.sown.AEG  <- remove_present(long = longAEG)
newly.sown.HEG  <- remove_present(long = longHEG)
newly.sown.SEG  <- remove_present(long = longSEG)
newly.sown.long <- rbind(newly.sown.AEG, newly.sown.HEG, newly.sown.SEG)

# order by plot and species
newly.sown.long <- newly.sown.long[order(newly.sown.long$plotIdent, 
                                         newly.sown.long$species) , ]



###############################################################
### remove intermediate data                   ################
###############################################################


rm(sown.species, sownAEG, sownHEG, sownSEG,
   longAEG, longHEG, longSEG, 
   newly.sown.AEG, newly.sown.HEG, newly.sown.SEG,
   remove_present, species.long )

