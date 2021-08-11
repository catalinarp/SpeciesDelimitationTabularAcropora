
## Pophelper with Sanger Sequenced markers - 4 loci dataset
## Following tutorial from: http://www.royfrancis.com/pophelper/articles/index.html#input-files

## 1. Load packages:
library(pophelper)
library(gridExtra)
library(strataG)
library(readr)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(pophelper)

## 2.Set working directory:
setwd("~/Dropbox/Research/1.ULB_PhD/Projects/Speciesdelimitation_tabularAcropora/MolecularLab/SangerSeqs/5.SpeciesDelimitation/Structure-based/Structure/")

## 3. Load dataset:
## STR1: Admixture and correlated frequencies w/o Morpho prior
filenames1 <- list.files(path = "STRauto1_5loci.out/STRauto1.allk/", pattern = "input5lociSS",full.names = TRUE)
slist1 <- readQ(files = filenames1, filetype = "structure")
dir.names <- list.dirs(full.names = TRUE)

#Create a table and a summary (this has to be done for each set of parameters if >1)
tr1 <- tabulateQ(qlist=slist1, writetable=TRUE,exportpath = dir.names[2])
sm1 <- summariseQ(tr1, writetable=TRUE, exportpath = dir.names[2])

## 4.Calculate best K according to Evanno (delta K) method
tr1 <- tabulateQ(qlist=slist1, writetable=F)
sm1 <- summariseQ(tr1, writetable=F)
p <- evannoMethodStructure(data=sm1,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p)

## 5. Create and align bar/distruct plots: for STR1
##    Note: Check tabulateQ file prev. created to id. the # runs when closer to "best K"
##    For STR1 k ~ 2-3, runs= 51-69 
#A) Using Pophelper ONLY (not sure if the alignment is correct) just to check stratification
slistaln1 <- alignK(slist1[c(51,53,55,58,60,63,65,67,70)])
p1 <- plotQ(slistaln1,imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11)
grid.arrange(p1$plot[[1]])

#B) Using Pophelper + CLUMPP to combine, align and merge runs for each K
clumpComb1 <- clumppExport(slist1) # To combine str files and generate clumpp input
# In the command line, run CLUMPP for each K you want to plot:
#Briefly: /Applications/Academic_software/CLUMPP_MacOSX.1.1.2/CLUMPP paramfile

clumpAln1K3 <- readQClumpp("STRauto1_5loci.out/clumppExport/pop_K3/pop_K3-combined-aligned.txt")
clumpAln1K2 <- readQClumpp("STRauto1_5loci.out/clumppExport/pop_K2/pop_K2-combined-aligned.txt")

metadata <- read.delim("LabelsInd.txt", sep="\t", header = FALSE,stringsAsFactors=F)
rownames(clumpAln1K3[[1]]) <- metadata$V1
rownames(clumpAln1K2[[1]]) <- metadata$V1
morphogroup <- metadata[,3,drop=F]

#Check the aligned plots first to see they all convey the same info. 
p1k2 <- plotQMultiline(clumpAln1K2,spl=36, height=8, width=14, useindlab=T,lpp=1,
                       sortind="Cluster2",barsize=1, grplab=morphogroup,
                       clustercol = brewer.pal(2,'YlGnBu'), ordergrp=TRUE)
p1k3 <- plotQMultiline(clumpAln1K3,spl=36, height=8, width=14, useindlab=T,lpp=1,
                       sortind="Cluster3",barsize=1, grplab=morphogroup,
                       clustercol = brewer.pal(3,'YlGnBu'), ordergrp=TRUE)

#Once that has been done, you can use the merged CLUMPP files to generate the final plots

clumpMer1K3 <- readQClumpp("STRauto1_5loci.out/clumppExport/pop_K3/pop_K3-combined-merged.txt")
clumpMer1K2 <- readQClumpp("STRauto1_5loci.out/clumppExport/pop_K2/pop_K2-combined-merged.txt")
rownames(clumpMer1K3[[1]]) <- metadata$V1
rownames(clumpMer1K2[[1]]) <- metadata$V1
morphogroup <- metadata[,3,drop=F]


#Plots: Build colorset factors according to Ks, make sure colors correspond
colorset2 <- c("goldenrod","royalblue4")
colorset2 <- c("#ce882a", "#7e6729")
colorset3 <- c("#FBA90A","#2c5ab4","#009933")
colorset3 <- c("#d5af29","#ce882a","#7e6729")

#colorset4 <- c("cyan4","goldenrod","forestgreen","royalblue4")

p1k2 <-  plotQMultiline(clumpMer1K2,spl=36, height=8, width=14, useindlab=T,showindlab=T,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, grplabsize = 10,
                        subtitlesize=8, clustercol = colorset2, showyaxis = T, ordergrp=TRUE,
                        subtitlelab="k=2", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabsize=7,returnplot=T,exportplot=F)
grid.arrange(p1k2$plot[[1]][[1]],nrow=1)


p1k3 <-  plotQMultiline(clumpMer1K3,spl=36, height=8, width=14, useindlab=T,subtitlelab="k=3",
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, showgrplab = T,
                        clustercol = colorset3, showyaxis = T, subtitlesize=8,ordergrp=TRUE,
                        barbordersize = 0.5, barbordercolour="white",subtitlehjust = 1,
                        indlabangle=90, indlabvjust = 0.9,indlabhjust = 0.5,indlabsize=7,
                        sortind ="label",returnplot=T,exportplot=F)
grid.arrange(p1k3$plot[[1]][[1]],nrow=1)

grid.arrange(p1k2$plot[[1]][[1]],p1k3$plot[[1]][[1]],nrow=2)

## Trying with all Ks for STR7
## STR7: NO admixture and uncorrelated frequencies w/o Morpho prior
filenames2 <- list.files(path = "STRauto7_5loci.out/all_k/", pattern = "input5lociSS",full.names = TRUE)
slist2 <- readQ(files = filenames2, filetype = "structure")
dir.names <- list.dirs(full.names = TRUE)

#Create a table and a summary (this has to be done for each set of parameters if >1)
tr2 <- tabulateQ(qlist=slist2, writetable=TRUE,exportpath = dir.names[2])
sm2 <- summariseQ(tr2, writetable=TRUE, exportpath = dir.names[2])

p2 <- evannoMethodStructure(data=sm2,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p2)

clumpComb2 <- clumppExport(slist2)

metadata <- read.delim("LabelsInd.txt", sep="\t", header = FALSE,stringsAsFactors=F)

clumpMer2K3 <- readQClumpp("STRauto7_5loci.out/all_k/pop_K3/pop_K3-combined-merged.txt")
clumpMer2K2 <- readQClumpp("STRauto7_5loci.out/all_k/pop_K2/pop_K2-combined-merged.txt")
clumpMer2K4 <- readQClumpp("STRauto7_5loci.out/all_k/pop_K4/pop_K4-combined-merged.txt")
rownames(clumpMer2K3[[1]]) <- metadata$V1
rownames(clumpMer2K2[[1]]) <- metadata$V1
rownames(clumpMer2K4[[1]]) <- metadata$V1
morphogroup <- metadata[,3,drop=F]

#Plots
dev.off()
colorset2 <- c("goldenrod","royalblue4")
colorset3 <- c("goldenrod","forestgreen","royalblue4")
colorset4 <- c("cyan4","forestgreen","royalblue4","goldenrod")

p7k2 <-  plotQMultiline(clumpMer2K2,spl=36, height=8, width=14, useindlab=T,showindlab=T,
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, grplabsize = 10,
                        subtitlesize=8, clustercol = colorset2, showyaxis = T, ordergrp=TRUE,
                        subtitlelab="k=2", barbordersize = 0.5, barbordercolour="white",
                        subtitlehjust = 1, indlabsize=7,returnplot=T,exportplot=F)

p7k3 <-  plotQMultiline(clumpMer2K3,spl=36, height=8, width=14, useindlab=T,subtitlelab="k=3",
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, showgrplab = T,
                        clustercol = colorset3, showyaxis = T, subtitlesize=8,ordergrp=TRUE,
                        barbordersize = 0.5, barbordercolour="white",subtitlehjust = 1,
                        indlabangle=90, indlabvjust = 0.9,indlabhjust = 0.5,indlabsize=7,
                        sortind ="label",returnplot=T,exportplot=F)

p7k4 <-  plotQMultiline(clumpMer2K4,spl=36, height=8, width=14, useindlab=T,subtitlelab="k=4",
                        lpp=1,barsize=1, grplab=morphogroup,showsubtitle=T, showgrplab = T,
                        clustercol = colorset4, showyaxis = T, subtitlesize=8,ordergrp=TRUE,
                        barbordersize = 0.5, barbordercolour="white",subtitlehjust = 1,
                        indlabangle=90, indlabvjust = 0.9,indlabhjust = 0.5,indlabsize=7,
                        sortind ="label",returnplot=T,exportplot=F)


grid.arrange(p7k2$plot[[1]][[1]],p7k3$plot[[1]][[1]],p7k4$plot[[1]][[1]], nrow=3)


