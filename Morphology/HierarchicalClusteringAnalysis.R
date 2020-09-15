### Hierarchical Clustering Analysis ###
# Includes: Descriptive/qualitative and categorical data ###

# 1. Set working directory
setwd("~/Projects/Speciesdelimitation_tabularAcropora/MorphometryTaxonomy/Analyses/")

# 2. Load packages
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(NbClust)
library(pvclust)
library(nomclust)
library(purrr)

## 3. Load the datasets
# NOTE:Rows should contain observations (or data points) and columns should be variables.
# Check if your data has any missing values, if yes, remove or impute them: data <- na.omit(data)

# Qualitative variables
descriptivechars <- read_delim("QLDescriptiveData_recodedDEP.txt",
                               "\t", escape_double = FALSE,
                               na = "NA", trim_ws = TRUE)
descriptivechars <- descriptivechars %>% mutate_if(is.character,as.factor)

summary(descriptivechars) #Check the dataset characteristics such as type of variables etc.

# Quantitative variables that were transformed to categorical
categoricalchars <- read_delim("QTMorphometricDataCategorizedDEP.txt", 
                               "\t", escape_double = FALSE,
                               na = "NA", trim_ws = TRUE)
categoricalchars <- categoricalchars %>% mutate_if(is.character,as.factor)

summary(categoricalchars) #Check the dataset characteristics such as type of variables etc.

# 4. Merge the categorical and the qualitative datasets
merged <- merge(descriptivechars, categoricalchars,
                all= TRUE, all.x=TRUE, all.y=TRUE,
                no.dups = TRUE, sort = TRUE)

Allmorphochars <-data.frame(merged [,3:22],
                            row.names = (merged$Colony_ID),
                            stringsAsFactors=TRUE)
summary(Allmorphochars) #Check the dataset characteristics such as type of variables etc.

# 5. Compute dissimilarity matrix for categorical or mixed data, choose one:

# A. The general dissimilarity coefficient of Gower (1971): 
# Each variable (column) is first standardized by dividing each entry by the range of the corresponding variable, after subtracting the minimum value. 
# The rescaled variable has range from 0 to 1.
disimatrix <- daisy(Allmorphochars[,2:20], metric = c("gower"))

### OR:

# B. The simple matching coefficient (Sokal, 1958):
# Represents the simplest way of measuring similarity and it does not impose any weights. 
# For a given variable, it assigns the value 1 in case of match and value 0 otherwise.

disimatrix <- sm(Allmorphochars[2:20])

# 6. After computing the dissimilarity matrix by the chosen method, compare clustering methods using AGNES function.
# NOTE: With the agnes function you can get the agglomerative coefficient, which measures the amount of clustering structure found.

# Clustering methods to be compared (m)
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# Create a function for comparison:
ac <- function(x) {
  agnes(disimatrix, method = x)$ac
}

# Check the result of the agglomerative coefficient (values closer to 1 suggest stronger clustering)
map_dbl(m, ac)

# 7. Generate the clustering dendrogram with the clustering method chosen or compare those ones that gave best agglomerative coefficient values.

# Average method:
hc1 <- agnes(disimatrix, method = "average") # Calculates clustering
pltree(hc1, cex = 0.6, hang = -1, main = "Dendrogram", labels = descriptivechars$Morphosp) # plots dendrogram using morpho-species label
pltree(hc2, cex = 0.6, hang = -1, main = "Dendrogram", labels = descriptivechars$Colony_ID) # plots dendrogram using colony ID label

# Ward method:
hc2 <- agnes(disimatrix, method = "ward") # Calculates clustering
pltree(hc2, cex = 0.6, hang = -1, main = "Dendrogram", labels = descriptivechars$Morphosp) # plots dendrogram using morpho-species label
pltree(hc2, cex = 0.6, hang = -1, main = "Dendrogram", labels = descriptivechars$Colony_ID) # plots dendrogram using colony ID label

# Complete method:
hc3 <- agnes(disimatrix, method = "complete") # Calculates clustering
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram", labels = descriptivechars$Morphosp) # plots dendrogram using morpho-species label
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram", labels = descriptivechars$Colony_ID) # plots dendrogram using colony ID label
