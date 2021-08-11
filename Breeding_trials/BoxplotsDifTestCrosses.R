#---------------------------------------------------------------------------------------
## Box-plots for cross-fertilization experiment
#---------------------------------------------------------------------------------------
#
# 1. Load required packages

library(ggplot2)
library(RColorBrewer)
library(readr)
library(bestNormalize)
library(PMCMRplus)
library(PMCMR)

# 2. Set working directory 

setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/Manuscript18Oki/BreedingData/")

# 3. Load dataset 
crosses <- read_delim("Breedingtrials18Oki.txt", 
                                    "\t", escape_double = FALSE, na = "na", 
                                    trim_ws = TRUE)
View(crosses)

# 3. Set factor: Type of experimental cross
typecross <- factor(crosses$Cross, 
                  levels = c("control_self", "self",
                             "interspecific cross",
                             "intraspecific cross"),
                  labels = c("Eggs only (control)", 
                             "Self", "Between morphospecies",
                             "Within morphospecies"))

# 4. Box Plots by cross types and fertilization ratio % (Embryos/Embryos+Unfertilized eggs)

boxplot_crosses <- ggplot(crosses, aes(x = typecross, y= FR_percentage)) + 
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot() + xlab("Crosses") +
  theme_bw() + theme(panel.border = element_rect(size = 1)) +
  theme(plot.title =  element_blank()) +  #coord_flip() +
  ylab("Mean fertilization success (%)") +
  geom_point() +  theme(legend.position = "none") 

plot(boxplot_crosses)

# 5. Test for significant differences in fertilization ratio between cross types:
# A) Check for normality and homocedasticity 

shapiro.test(crosses$Fertilization_ratio)
shapiro.test(crosses$FR_percentage)

leveneTest(crosses$Fertilization_ratio, typecross)
leveneTest(crosses$FR_percentage, typecross)

# B) Try to normalize data: Transform data

# Using Fertilization_ratio:
bestNormalize(crosses$Fertilization_ratio) #To see which transformation could work best
#Transform  according to the most frequent result obtained in the last step (run at least 5x)
x <- orderNorm(crosses$Fertilization_ratio)
x <- x$x.t
#Test for normality and homogeneity of variance:
shapiro.test(x) #Normality= NO
leveneTest(x, typecross) #Homogeneity of Variance= YES

# Using FR_percentage:
bestNormalize(crosses$FR_percentage) #To see which transformation could work best
#Transform:
x <- yeojohnson(crosses$FR_percentage)
x <- x$x.t
#Test for normality and homogeneity of variance:
shapiro.test(x) #Normality= NO
leveneTest(x, typecross) #Homogeneity of Variance= YES

# C) Non-parametric tests for significant difference 

kruskal.test(crosses$Fertilization_ratio ~ typecross,
             data = crosses) 
posthoc.kruskal.nemenyi.test(x=crosses$Fertilization_ratio,
                             g=typecross, 
                             dist="Chisquare")
posthoc.kruskal.dunn.test(x=crosses$Fertilization_ratio,
                          g=typecross,
                          p.adjust.method="bonferroni")
posthoc.kruskal.conover.test(x=crosses$Fertilization_ratio,
                             g=typecross,
                             p.adjust.method="bonferroni")
posthoc.vanWaerden.test(x=crosses$Fertilization_ratio,
                        g=typecross,
                        p.adjust.method="bonferroni")



