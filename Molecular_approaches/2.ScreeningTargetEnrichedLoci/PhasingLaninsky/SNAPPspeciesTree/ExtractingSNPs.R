### Extracting SNPs from DNA alignments ###
# Includes: 210 target-enriched loci that were recovered and aligned for all samples (n=9)

# 1. Load packages:
library(ape)
library(adegenet)
library(chopper)
library(phylotools)

# 2. Set working directory: 
setwd(dir ="~/Projects/Speciesdelimitation_tabularAcropora/MolecularAnalyses/TargetCapture/Laninsky_Phasing/Tree-based/working_fasta210/")

# 3. Importing fasta and pulling out SNPs simultaneously:

fl=list.files(pattern = ".fasta")

for (i in 1:length(fl)) assign(fl[i],
                               fasta2DNAbin(fl[i],snpOnly=T))
for(i in 1:length(fl)) {
  SNPs210 <- cbind (fasta2DNAbin(fl[i],snpOnly=T))
}

write.nexus.data(SNPs210, "../SNPsTC.nexus")

