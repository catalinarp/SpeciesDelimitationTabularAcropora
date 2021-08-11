
## Collapsing branches of low support in to polytomies in R ##
# Following: http://evoslav.blogspot.com/2015/01/how-to-collapse-unsupported-branches-of.html

library(ape) # Load the Ape library

#REMEMBER: Ran it for each marker separately:

## ==== AcroCR ==== 
#Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/18Oki_AcroCR/") 
#Read tree file:
tree <- read.tree("AcroCR_outgroup/AcroCR.mafft.fas_cps2.contree")
tree <- read.tree("AcroCR_outgroup_gapsrecoded/AcroCR.mafft.recodedgaps.fas_cps2.contree")

# Prune branches that have low support: Which nodes have low support?
Badnodes <- which(as.numeric(tree$node.label) < 85) + length(tree$tip.label)

# The number of those nodes are equal to the position in node.label plus the number of species (the first n nodes)
# The branch lengths INDEXES of each node are tree$edge[,1], and tree$edge[,2] are the nodes names
# The actual indexes of the nodes with low support in tree$edge
Badnodes_indexes <- c()
for(node in Badnodes){
  Badnodes_indexes <- c(Badnodes_indexes, which(tree$edge[,2] == node))
}

tree$edge.length[Badnodes_indexes] <- 0 # Make the branch length of the bad nodes = 0
tree_multi <- di2multi(tree) #Use di2multi function to convert the branches with lenght = 0, to multichotomies

write.tree(tree_multi, file = "AcroCR_cps2MOD.contree")
write.tree(tree_multi, file = "AcroCR_cps2MOD.gapsrecoded.contree")



## ==== TDH:E99029792 ==== 
#Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/Phased_outgroup/18Oki_E99029792/") 
#Read tree file:
tree1 <- read.tree("E99029792phased_cps2.contree")

# Prune branches that have low support: Which nodes have low support?
Badnodes <- which(as.numeric(tree1$node.label) < 85) + length(tree1$tip.label)

# The number of those nodes are equal to the position in node.label plus the number of species (the first n nodes)
# The branch lengths INDEXES of each node are tree$edge[,1], and tree$edge[,2] are the nodes names
# The actual indexes of the nodes with low support in tree$edge
Badnodes_indexes <- c()
for(node in Badnodes){
  Badnodes_indexes <- c(Badnodes_indexes, which(tree1$edge[,2] == node))
}

tree1$edge.length[Badnodes_indexes] <- 0 # Make the branch length of the bad nodes = 0
tree_multi <- di2multi(tree1) #Use di2multi function to convert the branches with lenght = 0, to multichotomies

write.tree(tree_multi, file = "E99029792phased_cps2MOD.contree")

## ==== DOPR: UCE111109 ==== 
#Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/Phased_outgroup/18Oki_UCE111109/") 
#Read tree file:
tree2 <- read.tree("UCE111109_k_alleles_cps2.contree")
Badnodes2 <- which(as.numeric(tree2$node.label) < 85) + length(tree2$tip.label)
Badnodes_indexes2 <- c()
for(node in Badnodes2){
  Badnodes_indexes2 <- c(Badnodes_indexes2, which(tree2$edge[,2] == node))
}

tree2$edge.length[Badnodes_indexes2] <- 0 # Make the branch length of the bad nodes = 0
tree_multi2 <- di2multi(tree2) #Use di2multi function to convert the branches with lenght = 0, to multichotomies

write.tree(tree_multi2, file = "UCE111109_k_alleles_cps2MOD.contree")

## ==== ASNA: E2711 ==== 
#Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/Phased_outgroup/18Oki_E2711/") 
#Read tree file:
tree3 <- read.tree("Exon2711_k_alleles_cps2.contree")
Badnodes3 <- which(as.numeric(tree3$node.label) < 85) + length(tree3$tip.label)
Badnodes_indexes3 <- c()
for(node in Badnodes3){
  Badnodes_indexes3 <- c(Badnodes_indexes3, which(tree3$edge[,2] == node))
}

tree3$edge.length[Badnodes_indexes3] <- 0 # Make the branch length of the bad nodes = 0
tree_multi3 <- di2multi(tree3) #Use di2multi function to convert the branches with lenght = 0, to multichotomies

write.tree(tree_multi3, file = "Exon2711_k_alleles_cps2MOD.contree")

## ==== PMCA ==== 
#Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/Phased_outgroup/18Oki_PMCA/") 
#Read tree file:
tree4 <- read.tree("PMCA_alleles_cps2.contree")
Badnodes4 <- which(as.numeric(tree4$node.label) < 85) + length(tree4$tip.label)
Badnodes_indexes4 <- c()
for(node in Badnodes4){
  Badnodes_indexes4 <- c(Badnodes_indexes4, which(tree4$edge[,2] == node))
}

tree4$edge.length[Badnodes_indexes4] <- 0 # Make the branch length of the bad nodes = 0
tree_multi4 <- di2multi(tree4) #Use di2multi function to convert the branches with lenght = 0, to multichotomies

write.tree(tree_multi4, file = "PMCA_alleles_cps2MOD.contree")

## ==== 5491L ==== 
#Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/Phased_outgroup/18Oki_5491L/") 
#Read tree file:
tree5 <- read.tree("5491L_alleles_cps2.contree")
Badnodes5 <- which(as.numeric(tree5$node.label) < 85) + length(tree5$tip.label)
Badnodes_indexes5 <- c()
for(node in Badnodes5){
  Badnodes_indexes5 <- c(Badnodes_indexes5, which(tree5$edge[,2] == node))
}

tree5$edge.length[Badnodes_indexes5] <- 0 # Make the branch length of the bad nodes = 0
tree_multi5 <- di2multi(tree5) #Use di2multi function to convert the branches with lenght = 0, to multichotomies

write.tree(tree_multi5, file = "5491L_alleles_cps2MOD.contree")

## ==== 3 TC loci:  ==== 
#Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/18Oki-Concatenated_outgroup/IUPAC_consensus/Loci_3nuTC/") 
#Read tree file:
tree <- read.tree("3nuTC_cps2.contree")
Badnodes <- which(as.numeric(tree$node.label) < 85) + length(tree$tip.label)
Badnodes_indexes <- c()
for(node in Badnodes){
  Badnodes_indexes <- c(Badnodes_indexes, which(tree$edge[,2] == node))
}

tree$edge.length[Badnodes_indexes] <- 0 # Make the branch length of the bad nodes = 0
tree_multi <- di2multi(tree) #Use di2multi function to convert the branches with lenght = 0, to multichotomies

write.tree(tree_multi, file = "3nuTC_cps2MOD.contree")

## ==== 5nu + 1 mt loci:  ==== 
#Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Tree-based/18Oki-Concatenated_outgroup/IUPAC_consensus/Loci_5nu1mt") 
#Read tree file:
tree <- read.tree("5nu1mt_cps2.contree")
Badnodes <- which(as.numeric(tree$node.label) < 85) + length(tree$tip.label)
Badnodes_indexes <- c()
for(node in Badnodes){
  Badnodes_indexes <- c(Badnodes_indexes, which(tree$edge[,2] == node))
}

tree$edge.length[Badnodes_indexes] <- 0 # Make the branch length of the bad nodes = 0
tree_multi <- di2multi(tree) #Use di2multi function to convert the branches with lenght = 0, to multichotomies

write.tree(tree_multi, file = "5nu1mt_cps2MOD.contree")


