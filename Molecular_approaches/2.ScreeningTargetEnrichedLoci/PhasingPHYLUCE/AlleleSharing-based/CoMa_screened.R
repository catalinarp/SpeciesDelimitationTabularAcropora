### Conspecificity matrix (CoMa) ###
# Includes: Output from HaplowebMaker (https://eeg-ebe.github.io/HaplowebMaker/)
# from Target-capture derived loci from the PHYLUCE phasing pipeline

# 1. Load packages:
library(readr)
library(heatmap3)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(scales)

# 2.Set working directory:
setwd("~/Projects/Speciesdelimitation_tabularAcropora/MolecularAnalyses/TargetCapture/Phyluce_Phasing/Allele_sharing-based/screen-fasta_80loci/")

# 3. Load dataset: Using Gaps as 5th character
gblocks <- read_delim("Gaps5thchar-all/CoMa.tsv",
                          "\t", escape_double = FALSE,trim_ws = TRUE)
matrix1 <- as.matrix(gblocks[,-1])
Samples1 <- factor(x= gblocks$X1,
                   levels= c("Acropora_CFhyacinthus1C282",
                             "Acropora_CFhyacinthus1C283","Acropora_CFhyacinthus1C284",
                             "Acropora_CFcytherea5C285", "Acropora_CFcytherea5C286",
                             "Acropora_CFcytherea5C287", "Acropora_CFbifurcataC288",
                             "Acropora_CFbifurcataC289", "Acropora_CFbifurcataC290"),
                   labels= c("18Oki21", "18Oki22","18Oki23","18Oki26", "18Oki27",
                             "18Oki29", "18Oki32","18Oki33","18Oki34"))
rownames(matrix1) <- Samples1
colnames(matrix1) <- Samples1

# 4. Set morpho-species as a factor to add color bar:

morphogroup <- read.delim("individuals.txt",header = FALSE)
morphogroups <- factor(morphogroup$V2,
                       levels= c("Ahya", "Abif", "Acyt"),
                       labels= c("A. cf. hyacinthus", "A. cf. bifurcata", "A.cf. cytherea"))
colormorpho <- c("#7e6729","#d5af29","#ce882a")[morphogroups]

# 5. Draw CoMa using heatmap package:

pdf(file='CoMa_Gaps5thchar_morphospp.pdf')
heatmap3(matrix1, cexRow=0.9, cexCol=0.9, revC=TRUE, margins=c(7,7),
         col=colorRampPalette(c("#ffffff", "red"))(5000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=colormorpho, 
         ColSideLabs = "FFRs",
         labCol=Samples1,labRow=Samples1,
         xlab="Samples")
dev.off()

# 6. Load dataset: Using Gaps as missing data
gblocks2 <- read_delim("GapsMissdata-all/CoMa.tsv",
                      "\t", escape_double = FALSE,trim_ws = TRUE)
matrix2 <- as.matrix(gblocks2[,-1])
Samples2 <- factor(x= gblocks2$X1,
                      levels= c("Acropora_CFhyacinthus1C282",
                      "Acropora_CFhyacinthus1C283","Acropora_CFhyacinthus1C284",
                      "Acropora_CFcytherea5C285", "Acropora_CFcytherea5C286",
                      "Acropora_CFcytherea5C287", "Acropora_CFbifurcataC288",
                      "Acropora_CFbifurcataC289", "Acropora_CFbifurcataC290"),
                      labels= c("18Oki21", "18Oki22","18Oki23","18Oki26", "18Oki27",
                                "18Oki29", "18Oki32","18Oki33","18Oki34"))
rownames(matrix2) <- Samples2
colnames(matrix2) <- Samples2

# 7. Draw CoMa using heatmap package:

pdf(file='CoMa_GapsMissdata_morphospp.pdf')
heatmap3(matrix2, cexRow=0.9, cexCol=0.9, revC=TRUE, margins=c(7,7),
         col=colorRampPalette(c("#ffffff", "red"))(5000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=colormorpho, 
         ColSideLabs = "FFRs",
         labCol=Samples2,labRow=Samples2,
         xlab="Samples")
dev.off()

