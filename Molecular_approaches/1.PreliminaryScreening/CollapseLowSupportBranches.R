### Collapsing branches of low support in to polytomies in R ###
# Following: http://evoslav.blogspot.com/2015/01/how-to-collapse-unsupported-branches-of.html

# 1. Load package:
library(ape)

# NOTE: Run it for each marker separately.
# Then, clear objects from space and run the next:

## ==== AcroCR ==== 
#Set working directory:
setwd("~/Projects/Speciesdelimitation_tabularAcropora/MolecularAnalyses/Tree-based/AcroCR/") 

#Read one of the following tree files at the time:
tree <- read.tree("AcroCR_cps2.contree")
tree <- read.tree("AcroCR.recodedgaps_cps2.contree") # Long gaps re-coded to a single base gap

# Remove branches that have low support, in this case less than 85 bootstrap value
Badnodes <- which(as.numeric(tree$node.label) < 85) + length(tree$tip.label)

Badnodes_indexes <- c()
for(node in Badnodes){
  Badnodes_indexes <- c(Badnodes_indexes, which(tree$edge[,2] == node))
}

tree$edge.length[Badnodes_indexes] <- 0 # Make the branch length of the low support nodes = 0
tree_multi <- di2multi(tree) #Use di2multi function to convert the branches with length = 0, to multichotomies

#Save the resulting tree:
write.tree(tree_multi, file = "AcroCR_cps2.contree")
write.tree(tree_multi, file = "AcroCR_cps2.gapsrecoded.contree")

## ==== PMCA ==== 
#Set working directory:
setwd("~/Projects/Speciesdelimitation_tabularAcropora/MolecularAnalyses/Tree-based/PMCA/") 
#Read tree file:
tree4 <- read.tree("PMCA_alleles_cps2.contree")

# Remove branches that have low support, in this case less than 85 bootstrap value
Badnodes4 <- which(as.numeric(tree4$node.label) < 85) + length(tree4$tip.label)
Badnodes_indexes4 <- c()
for(node in Badnodes4){
  Badnodes_indexes4 <- c(Badnodes_indexes4, which(tree4$edge[,2] == node))
}

tree4$edge.length[Badnodes_indexes4] <- 0 # Make the branch length of the low support nodes = 0
tree_multi4 <- di2multi(tree4) #Use di2multi function to convert the branches with length = 0, to multichotomies

#Save the resulting tree:
write.tree(tree_multi4, file = "PMCA_alleles_cps2MOD.contree")

## ==== FZD ==== 
#Set working directory:
setwd("~/Projects/Speciesdelimitation_tabularAcropora/MolecularAnalyses/Tree-based/FZD/") 
#Read tree file:
tree5 <- read.tree("FZD_alleles_cps2.contree")

# Remove branches that have low support, in this case less than 85 bootstrap value
Badnodes5 <- which(as.numeric(tree5$node.label) < 85) + length(tree5$tip.label)
Badnodes_indexes5 <- c()
for(node in Badnodes5){
  Badnodes_indexes5 <- c(Badnodes_indexes5, which(tree5$edge[,2] == node))
}

tree5$edge.length[Badnodes_indexes5] <- 0 # Make the branch length of the low support nodes = 0
tree_multi5 <- di2multi(tree5) #Use di2multi function to convert the branches with length = 0, to multichotomies

#Save the resulting tree:
write.tree(tree_multi5, file = "FZD_alleles_cps2MOD.contree")

