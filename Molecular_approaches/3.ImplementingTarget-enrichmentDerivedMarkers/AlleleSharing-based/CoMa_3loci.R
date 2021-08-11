
## Set working directory

setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Allele_sharing-based/18Oki-3TCloci/")

## Load packages 
library(readr)
library(heatmap3)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(scales)

## Load data
nuclearloci3 <- read_delim("1.After_M13seq_3loci_known/CoMa.tsv",
                          "\t", escape_double = FALSE,trim_ws = TRUE)
full_matrix <- as.matrix(nuclearloci3[,-1])
Samples1 <- as.factor(x= nuclearloci3$X1)
rownames(full_matrix) <- Samples1
colnames(full_matrix) <- Samples1

# Color by morphogroup

morphogroup <- read.delim("individuals.txt",header = FALSE)
morphogroups <- factor(morphogroup$V2,
                       levels= c("Ahya", "Acyt", "Abif"),
                       labels= c("A. hyacinthus", "A.cytherea", "A. bifurcata"))
colormorpho <- c("#d5af29","#ce882a","#7e6729")[morphogroups]

# Draw heatmap

pdf(file='1.After_M13seq_3loci_known/CoMa_3loci-morphogroups_color2.pdf')
heatmap3(full_matrix, cexRow=0.9, cexCol=0.9, revC=TRUE, margins=c(6,9),
       col=colorRampPalette(c("#e8e7e7", "red"))(5000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=colormorpho, 
         ColSideLabs = "FFRs", 
         labCol=Samples1,labRow=Samples1,
         xlab="Samples")
dev.off()


## Load data with gaps as missing data
nuclearloci3gapsm <- read_delim("2.After_M13seq_known_gapsmiss/CoMa.tsv",
                           "\t", escape_double = FALSE,trim_ws = TRUE)
full_matrixgap <- as.matrix(nuclearloci3gapsm[,-1])
Samplesgap <- as.factor(x= nuclearloci3gapsm$X1)
rownames(full_matrixgap) <- Samplesgap
colnames(full_matrixgap) <- Samplesgap

# Draw heatmap

pdf(file='2.After_M13seq_known_gapsmiss/CoMa_3locigaps-morphogroups.pdf')
heatmap3(full_matrixgap, cexRow=0.9, cexCol=0.9, revC=TRUE, margins=c(6,9),
         col=colorRampPalette(c("#a5c3fe", "red"))(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=colormorpho, 
         ColSideLabs = "Morphological groups",
         labCol=Samplesgap,labRow=Samplesgap,
         xlab="Samples")
dev.off()




#Re-scale values

dataset_matrix <- melt(full_matrix)
summary(dataset_matrix)
qplot(dataset_matrix, geom = "histogram", bins = 50)
dataset_matrix_rescaled <- cbind(dataset_matrix,rescale(dataset_matrix$value, to = c(0,1), from = range(dataset_matrix$value)))
summary(dataset_matrix_rescaled)

full_matrix_rescaled <- rescale(full_matrix,to = c(0,1))

# New plot
pdf(file='CoMa_untrimmed-morphogroups_rescaled.pdf')
heatmap3(full_matrix_rescaled, cexRow=0.8, cexCol=0.8, revC=TRUE, margins=c(10,10),
                  col=colorRampPalette(c("beige","red"),
                                       interpolate=c("spline"), bias = 100)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=colormorpho, 
         ColSideLabs = "Morphological groups",
         labCol=Samples1,labRow=Samples1,
         xlab="Samples")
dev.off()

