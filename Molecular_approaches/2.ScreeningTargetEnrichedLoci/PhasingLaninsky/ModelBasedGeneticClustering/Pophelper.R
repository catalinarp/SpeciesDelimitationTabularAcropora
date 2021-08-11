### Model-based genetic clustering of available molecular Sanger Sequenced markers ###
# Following tutorial from: http://www.royfrancis.com/pophelper/articles/index.html#input-files
# Includes: Structure output for Target Capture-gblocks dataset from the Laninsky phasing pipeline

# 1. Load packages:
library(pophelper)
library(gridExtra)
library(strataG)
library(readr)
library(gplots)
library(RColorBrewer)

# 2.Set working directory:
setwd("~/Projects/Speciesdelimitation_tabularAcropora/MolecularAnalyses/TargetCapture/Laninsky_Phasing/Structure/STR_allk.Pophelper/")

# 3. Load dataset:
## STR1: Admixture and correlated frequencies w/o Morpho-species prior
filenames1 <- list.files(path = "STRauto1_allk/", pattern = "gblocks",full.names = TRUE)
slist1 <- readQ(files = filenames1, filetype = "structure")

# 4. Create a table and a summary
tr1 <- tabulateQ(qlist=slist1, writetable=TRUE,exportpath = dir.names[2])
sm1 <- summariseQ(tr1, writetable=TRUE, exportpath = dir.names[2])

# 5. Calculate best K according to Evanno (delta K) method
tr1 <- tabulateQ(qlist=slist1, writetable=FALSE)
sm1 <- summariseQ(tr1, writetable=FALSE)
p <- evannoMethodStructure(data=sm1,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p)

# 6. Use Pophelper + CLUMPP to combine, align and merge runs for each K:
clumpComb1 <- clumppExport(slist1) # To combine str files
# NOTE: Remember to move them to the appropriate directory in: STRauto1_allk

# NOTE: In the command line, run CLUMPP for each K you want to plot:
# Example: /Academic_software/CLUMPP_MacOSX.1.1.2/CLUMPP paramfile
# Now load data from CLUMPP output for the most likely best Ks:

clumpMerged1K2 <- readQClumpp("STRauto1_allk/pop_K2/pop_K2-combined-merged.txt")
clumpMerged1K3 <- readQClumpp("STRauto1_allk/pop_K3/pop_K3-combined-merged.txt")
clumpMerged1K4 <- readQClumpp("STRauto1_allk/pop_K4/pop_K4-combined-merged.txt")
clumpMerged1K5 <- readQClumpp("STRauto1_allk/pop_K5/pop_K5-combined-merged.txt")

# 7. Load metadata (individual labels):
metadata <- read.delim("LabelsInd.txt", sep=" ", header = FALSE,stringsAsFactors=F)

rownames(clumpMerged1K2[[1]]) <- metadata$V3
rownames(clumpMerged1K3[[1]]) <- metadata$V3
rownames(clumpMerged1K4[[1]]) <- metadata$V3
rownames(clumpMerged1K5[[1]]) <- metadata$V3

morphogroup <- metadata[,1,drop=F]

# NOTE: Build color sets according to the number of Ks.
# Examples:
colorset2 <- c("#ce882a","#7f682a")
colorset3 <- c("#7f682a","#ce882a","#d5af29")
colorset4 <- c("#7f682a","#ce882a","cyan4","#d5af29")
colorset5 <- c("cyan4","#0E682C","#ce882a","#7f682a","#d5af29")

# 8.Use the merged CLUMPP files to generate the final plots:
# Example with K=2
p1k2 <-  plotQMultiline(clumpMerged1K2,spl=9, height=8, width=14, useindlab=F,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, grplabsize = 14,
                        subtitlesize=14, clustercol = colorset2, showyaxis = T,
                        subtitlelab="k=2", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabangle=90, indlabvjust = 0.9,
                        indlabhjust = 0.5,indlabsize=12,returnplot=T,exportplot=F)
grid.arrange(p1k2$plot[[1]][[1]],nrow=1)

# Example with K=3
p1k3 <-  plotQMultiline(clumpMerged1K3,spl=9, height=8, width=14, useindlab=T,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, grplabsize = 14,
                        subtitlesize=14, clustercol = colorset3, showyaxis = T,
                        subtitlelab="k=3", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabangle=90, indlabvjust = 0.9,
                        indlabhjust = 0.5,indlabsize=12,returnplot=T,exportplot=F)
grid.arrange(p1k3$plot[[1]][[1]],nrow=1)

# Example with K=4
p1k4 <-  plotQMultiline(clumpMerged1K4,spl=9, height=8, width=14, useindlab=T,showindlab=T,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, showgrplab = T,
                        subtitlesize=8, clustercol = colorset4, showyaxis = T,
                        subtitlelab="k=4", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabsize=7,returnplot=T,exportplot=F)
grid.arrange(p1k4$plot[[1]][[1]],nrow=1)

# Example with K=5
p1k5 <-  plotQMultiline(clumpMerged1K5,spl=9, height=8, width=14, useindlab=T,showindlab=T,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, showgrplab = T,
                        subtitlesize=8, clustercol = colorset5, showyaxis = T,
                        subtitlelab="k=5", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabsize=7,returnplot=T,exportplot=F)
grid.arrange(p1k5$plot[[1]][[1]],nrow=1)


# Example of plot including k=3, k=4, and k=5
p1k3.2 <-  plotQMultiline(clumpMerged1K3,spl=9, height=8, width=14, useindlab=T,showindlab=F,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, grplabsize = 12,
                        subtitlesize=10, clustercol = colorset3, showyaxis = T,
                        subtitlelab="k=3", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabangle=90, indlabvjust = 0.9,indlabsize=10,
                        indlabhjust = 0.5,returnplot=T,exportplot=F)
p1k4.2 <-  plotQMultiline(clumpMerged1K4,spl=9, height=8, width=14, useindlab=T,showindlab=F,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, showgrplab = F,
                        subtitlesize=10, clustercol = colorset4, showyaxis = T,
                        subtitlelab="k=4", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabsize=10, returnplot=T,exportplot=F)
p1k5.2 <-  plotQMultiline(clumpMerged1K5,spl=9, height=8, width=14, useindlab=T,showindlab=T,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, showgrplab = F,
                        subtitlesize=10, clustercol = colorset5, showyaxis = T,
                        subtitlelab="k=5", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabsize=10,returnplot=T,exportplot=F)

grid.arrange(p1k3.2$plot[[1]][[1]],
             p1k4.2$plot[[1]][[1]],
             p1k5.2$plot[[1]][[1]],nrow=3)
