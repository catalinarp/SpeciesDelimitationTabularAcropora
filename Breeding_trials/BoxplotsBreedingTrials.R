### Box-plots for breeding trials ###
# Fertilization success in cross-fertilization experiments

# 1. Set working directory
setwd("~/Projects/Speciesdelimitation_tabularAcropora/BreedingTrials/Analyses/")

# 2. Load packages
library(ggplot2)
library(RColorBrewer)
library(readr)

# 3. Load the data set
crosses <- read_delim("BreedingTrials.txt", 
                                    "\t", escape_double = FALSE, na = "na", 
                                    trim_ws = TRUE)
View(crosses)

# 4. Create factors

typecross <- factor(crosses$category, 
                  levels = c("control_self", "self",
                             "interspecific cross",
                             "intraspecific cross"),
                  labels = c("Eggs only (control)", 
                             "Self", "Between morphospecies",
                             "Within morphospecies"))
crosscat <- factor(crosses$cross, 
                   levels = c("control-bifurcata",
                              "control-cytherea",
                              "control-hyacinthus",
                              "bifurcata-cytherea",
                              "hyacinthus-bifurcata",
                              "hyacinthus-cytherea",
                              "bifurcata-bifurcata",
                              "cytherea-cytherea",
                              "hyacinthus-hyacinthus",
                              "self-bifurcata",
                              "self-cytherea",
                              "self-hyacinthus"),
                   labels = c("Eggs only (ACFbif)",
                              "Eggs only (ACFcyt)",
                              "Eggs only (ACFhya)",
                              "ACFbif x ACFcyt",
                              "ACFhya x ACFbif",
                              "ACFhya x ACFcyt",
                              "ACFbif x ACFbif",
                              "ACFcyt x ACFcyt",
                              "ACFhya x ACFhya",
                              "Self -  ACFbif",
                              "Self -  ACFcyt",
                              "Self -  ACFhya"))

# 4. Create box plots by type of cross (Embryos)

boxplot_width1 <- ggplot(crosses, aes(x = typecross, y= embryos)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() + xlab("Crosses category") +
  theme(axis.text.x = element_text(face = "italic")) + 
  ylab("Number of embryos (10h)") +
  geom_point()
plot(boxplot_width1)

# 5. Create box plots by type of cross (Embryos/Embryos+Unfertilized eggs)

boxplot_width2 <- ggplot(crosses, aes(x = typecross, y=fertratio)) + 
  stat_boxplot(geom ='errorbar') + 
   geom_boxplot() + xlab("Crosses") +
  theme_bw() + theme(panel.border = element_rect(size = 1)) +
  theme(plot.title =  element_blank()) + 
  ylab("Mean fertilization success (10 h)") +
  geom_point() +  theme(legend.position = "none") 

plot(boxplot_width2)

#6. Box Plots by cross category (Embryos/Embryos+Unfertilized eggs)

boxplot_width3 <- ggplot(crosses, aes(x = crosscat, y=fertratio)) + 
  stat_boxplot(geom ='errorbar') + 
   geom_boxplot() + xlab("Crosses") +
  theme_bw() + theme(panel.border = element_rect(size = 1)) +
  theme(plot.title =  element_blank()) +  
  ylab("Mean fertilization success (10 h)") +
  geom_point() +  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

plot(boxplot_width3)

