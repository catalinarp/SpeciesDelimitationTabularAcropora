### Genetic distances estimation and histograms###
# Includes: Data from alignments: 
# -> Consensus sequences: mitochondrial AcroCR
# -> Phased alleles: FZD and PMCA performed in Mafft 
# NOTE: Remember to run the analyses one marker at a time

# 1. Set working directory
setwd("~/Projects/Speciesdelimitation_tabularAcropora/MolecularAnalyses/Distances")

# 2. Load packages:
library(readr)
library(ade4)
library(adegenet)
library(ape)
library(reshape2)
library(otuSummary)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# 3. Load genetic and morphological data from one marker at a time.
# Then, clear objects from space and run the next:

## AcroCR
DNA_aln <- fasta2DNAbin(file = "AcroCR.mafft.fasta")

Morphogroups <- read_delim("AcroCR.mafft.morphogroups.txt",
                                  "\t", escape_double = FALSE,
                                  col_names = FALSE,
                                  trim_ws = TRUE)

# NOTE: Check if both the metadata and the sequences have the same dimensions and info:
all(Morphogroups$X1==rownames(AcroCR_DNA))

## PMCA
DNA_aln <- fasta2DNAbin(file = "PMCA.alleles.fasta")

Morphogroups <- read_delim("morpho.PMCA_alleles.fasta.txt",
                           "\t", escape_double = FALSE,
                           col_names = FALSE,
                           trim_ws = TRUE)

# NOTE: Check if both the metadata and the sequences have the same dimensions and info:
all(Morphogroups$X1==rownames(DNA_aln))

## FZD
DNA_aln <- fasta2DNAbin(file = "FZD.mafft.alleles.fasta")


Morphogroups <- read_delim("morpho.FZD_alleles.fasta.txt",
                                        "\t", escape_double = FALSE,
                                        col_names = FALSE,
                                        trim_ws = TRUE)

# NOTE: Check if both the metadata and the sequences have the same dimensions and info:
all(Morphogroups$X1==rownames(DNA_aln))

# 4. Subset individuals according to morpho-species

Ahya <-subset(Morphogroups$X1,
              Morphogroups$X2=="Ahya")
Acyt <-subset(Morphogroups$X1,
              Morphogroups$X2=="Acyt")
Abif <-subset(Morphogroups$X1,
              Morphogroups$X2=="Abif")


# 5. Calculate the genetic distances according to the best model according to by BIC (ModelFinder Plus)
## NOTE: If the best supported is not available use the closest you can find available in dist.dna

DNA_aln_dist <- dist.dna(x = DNA_aln,
                       model = "TN93", #gamma=4, # Change model here
                       pairwise.deletion = TRUE)

# 6. Convert lower triangle of distances to a data frame:

DNA_aln_distpairwise <- matrixConvert(DNA_aln_dist,
                           colname = c("Ind1", "Ind2", "pwdist"))
  
summary(DNA_aln_distpairwise$pwdist)

# 7. Add columns with pairwise qualifiers:

Markerpwdmorph1 <- mutate(DNA_aln_distpairwise,
                 morpho1 = factor(case_when(Ind1 %in% Ahya ~ "Ahya",
                                            Ind1 %in% Acyt ~ "Acyt",
                                            Ind1 %in% Abif ~ "Abif",
                                             TRUE ~ NA_character_)))
Markerpwdmorph2 <- mutate(DNA_aln_distpairwise,
                 morpho2 = factor(case_when(Ind2 %in% Ahya ~ "Ahya",
                                            Ind2 %in% Acyt ~ "Acyt",
                                            Ind2 %in% Abif ~ "Abif",
                                              TRUE ~ NA_character_)))

Markerpwmorph <- merge(Markerpwdmorph1,Markerpwdmorph2)

Markerpwmorph$comp <- ifelse(Markerpwmorph$morpho1==Markerpwmorph$morpho2,
                              "intraspecific", "interspecific")

# 8. Perform histogram to detect potential barcoding gaps:

cols <- c("intraspecific" = "turquoise4",
          "interspecific" = "grey63")

a <- ggplot(Markerpwmorph, aes(x=pwdist, fill=factor(comp, levels=c("intraspecific","interspecific")))) +
  geom_histogram(bins = 9) + theme_classic() +
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

## 6. Perform histogram coloring by intraespecific groups according to morpho-species:

Markerpwdmorph1 <- mutate(DNA_aln_distpairwise,
                          morpho1 = factor(case_when(Ind1 %in% Ahya ~ "Ahya",
                                                     Ind1 %in% Acyt ~ "Acyt",
                                                     Ind1 %in% Abif ~ "Abif",
                                                     TRUE ~ NA_character_)))
Markerpwdmorph2 <- mutate(DNA_aln_distpairwise,
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

cols2 <- c("Ahya" = "#ce882a", #assess colors
           "Acyt" = "#7e6729",
           "Abif" = "#d5af29",
           "interspecific" = "grey63")

# Plot:

b <- ggplot(Markerpwdmorphgroup1,
            aes(x=pwdist, fill=factor(comp, levels=c("Ahya","Acyt","Abif","interspecific")))) +
  geom_histogram(bins = 5) + theme_classic() + #Bins can be readjusted if necessary  
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
