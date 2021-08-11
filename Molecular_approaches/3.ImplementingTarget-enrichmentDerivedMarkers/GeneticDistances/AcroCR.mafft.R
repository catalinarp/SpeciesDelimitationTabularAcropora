
#Distance-based- spp delimitation

## 0. Load packages
library(readr)
library(ade4)
library(adegenet)
library(ape)
library(reshape2)
library(otuSummary)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

## 1. Set working directory and load AcroCR genetic data

setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Distance-based/AcroCR/")

AcroCR_DNA <- fasta2DNAbin(file = "AcroCR.mafft.fasta")

AcroCR_morphogroups <- read_delim("AcroCR.mafft.morphogroups.txt",
                                        "\t", escape_double = FALSE,
                                        col_names = FALSE,
                                        trim_ws = TRUE)

# Subset individuals according to morphogroups
Ahya <-subset(AcroCR_morphogroups$X1,
              AcroCR_morphogroups$X2=="Ahya")
Acyt <-subset(AcroCR_morphogroups$X1,
              AcroCR_morphogroups$X2=="Acyt")
Abif <-subset(AcroCR_morphogroups$X1,
              AcroCR_morphogroups$X2=="Abif")

# Check if both the metadata and the sequences have the same dimensions and info:
all(AcroCR_morphogroups$X1==rownames(AcroCR_DNA))

## 2. Calculate the genetic distances according to the best model suported by AIC (Jmodeltest2):
AcroCRdist <- dist.dna(x = AcroCR_DNA,
                       model = "F81", 
                       pairwise.deletion = TRUE)

## 3. Convert lower triangle of distances to a dataframe

AcroCRpairwise <- matrixConvert(AcroCRdist,
                           colname = c("Ind1", "Ind2", "pwdist"))
  
summary(AcroCRpairwise$pwdist)

## 4. Add columns with pairwise qualifiers

AcroCRpwmorph1 <- mutate(AcroCRpairwise,
                 morpho1 = factor(case_when(Ind1 %in% Ahya ~ "Ahya",
                                            Ind1 %in% Acyt ~ "Acyt",
                                            Ind1 %in% Abif ~ "Abif",
                                             TRUE ~ NA_character_)))
AcroCRpwmorph2 <- mutate(AcroCRpairwise,
                 morpho2 = factor(case_when(Ind2 %in% Ahya ~ "Ahya",
                                            Ind2 %in% Acyt ~ "Acyt",
                                            Ind2 %in% Abif ~ "Abif",
                                              TRUE ~ NA_character_)))

AcroCRpwmorph <- merge(AcroCRpwmorph1,AcroCRpwmorph2)

AcroCRpwmorph$comp <- ifelse(AcroCRpwmorph$morpho1==AcroCRpwmorph$morpho2,
                              "intraspecific", "interspecific")

## 5. Perform histogram to detect potential barcode gaps

cols <- c("intraspecific" = "turquoise4",
          "interspecific" = "grey63")

a <- ggplot(AcroCRpwmorph, aes(x=pwdist, fill=factor(comp, levels=c("intraspecific","interspecific")))) +
  geom_histogram(bins = 20) + theme_classic() +
  scale_fill_manual(values = cols,
                    name= "Comparison",
                    labels = c("Intraspecific","Interspecific")) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size=12, vjust=-1),
        axis.text.x  = element_text(vjust=0.5, size=9),
        axis.title.y = element_text(size = 12, vjust=3),
        axis.text.y  = element_text(vjust=0.5, size=9)) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm")) +
  xlab("Pairwise genetic distance") +
  ylab ("Number of comparisons")
a

## 6. Perform histogram coloring by intraespecific groups according to morph

Markerpwdmorph1 <- mutate(AcroCRpairwise,
                          morpho1 = factor(case_when(Ind1 %in% Ahya ~ "Ahya",
                                                     Ind1 %in% Acyt ~ "Acyt",
                                                     Ind1 %in% Abif ~ "Abif",
                                                     TRUE ~ NA_character_)))
Markerpwdmorph2 <- mutate(AcroCRpairwise,
                          morpho2 = factor(case_when(Ind2 %in% Ahya ~ "Ahya",
                                                     Ind2 %in% Acyt ~ "Acyt",
                                                     Ind2 %in% Abif ~ "Abif",
                                                     TRUE ~ NA_character_)))

Markerpwmorphgroups <- merge(Markerpwdmorph1,Markerpwdmorph2)

Markerpwmorphgroups$comp <- ifelse(Markerpwmorphgroups$morpho1==Markerpwmorphgroups$morpho2,
                                   "intraspecific", "interspecific")

Markerpwdmorphgroup1 <- mutate(Markerpwmorphgroups,
                               comp = factor(case_when(morpho1 == "Ahya" & comp == "intraspecific" ~ "Ahya",
                                                       morpho1 == "Acyt" & comp == "intraspecific" ~ "Acyt",
                                                       morpho1 == "Abif" & comp == "intraspecific" ~ "Abif",
                                                       comp == "interspecific" ~ "interspecific",
                                                       TRUE ~ NA_character_)))
#Now create the histogram:

cols2 <- c("Ahya" = "#2c5ab4", #for colors according to genetic clustering
           "Acyt" = "#009933",
           "Abif" = "#FBA90A",
           "interspecific" = "grey63")

cols2 <- c("Ahya" = "#ce882a", #by morphogroup mod colors
           "Acyt" = "#7e6729",
           "Abif" = "#d5af29",
           "interspecific" = "grey63")

b <- ggplot(Markerpwdmorphgroup1,
            aes(x=pwdist, 
                fill=factor(comp, levels=c("Ahya","Acyt","Abif","interspecific")))) +
  geom_histogram(bins = 9) + theme_classic() +
  scale_fill_manual(values = cols2,
                    name= "Comparison",
                    labels = c( "A. cf. hyacinthus","A. cf. cytherea","A. cf. bifurcata",
                                "Interspecific")) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(size=12, vjust=-1),
        axis.text.x  = element_text(vjust=0.5, size=9),
        axis.title.y = element_text(size = 12, vjust=3),
        axis.text.y  = element_text(vjust=0.5, size=9)) +
  theme(plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm")) +
  theme(legend.position = "none") +
  xlab("Pairwise genetic distance") +
  ylab ("Number of comparisons")
b