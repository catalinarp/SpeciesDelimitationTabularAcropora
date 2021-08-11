### Kruskal Wallis rank sum test (non-parametric test) ###
# Test for differences in fertilization success in breeding trials

# 1. Set working directory
setwd("~/Projects/Speciesdelimitation_tabularAcropora/BreedingTrials/Analyses/")

# 2. Load packages
library(PMCMRplus)
library(PMCMR)

# 3. Load the data set
crosses <- read_delim("BreedingTrials.txt", 
                      "\t", escape_double = FALSE, na = "na", 
                      trim_ws = TRUE)
View(crosses)

# 4. Create factors
typecross <- factor(crosses$category, 
                    levels = c("control_self", "interspecific cross",
                               "intraspecific cross", "self"),
                    labels = c("Control (only eggs)", 
                               "Interespecific", "Intraspecific",
                               "Self-crosses"))

# 5. Perform the test

# 5A. Using embryos
kruskal.test(crosses$embryos ~ typecross,
             data = crosses) 
posthoc.kruskal.nemenyi.test(x=crosses$embryos,
                             g=typecross, 
                             dist="Chisquare")
posthoc.kruskal.dunn.test(x=crosses$embryos,
                           g=typecross,
                          p.adjust.method="bonferroni")
posthoc.kruskal.conover.test(x=crosses$embryos,
                             g=typecross,
                             p.adjust.method="bonferroni")
posthoc.vanWaerden.test(x=crosses$embryos,
                        g=typecross,
                        p.adjust.method="bonferroni")

# 5B. Using fertilization success rate

kruskal.test(crosses$fertratio ~ typecross,
             data = crosses) 
posthoc.kruskal.nemenyi.test(x=crosses$embryos,
                             g=typecross, 
                             dist="Chisquare")
posthoc.kruskal.dunn.test(x=crosses$embryos,
                          g=typecross,
                          p.adjust.method="bonferroni")
posthoc.kruskal.conover.test(x=crosses$embryos,
                             g=typecross,
                             p.adjust.method="bonferroni")
posthoc.vanWaerden.test(x=crosses$embryos,
                        g=typecross,
                        p.adjust.method="bonferroni")



