### Discriminant Analysis of Principal Components (DAPC) ###
# Data: Three Sanger Sequenced & target-capture derived markers

## 1. Load packages:
library(adegenet)
library(readr)

# 2.Set working directory:
setwd("~/Projects/Speciesdelimitation_tabularAcropora/MolecularAnalyses/SangerSeqs/5.SpeciesDelimitation/Structure-based/DAPC/")

## 3. Load data sets: STRUCTURE file to genind file

genotypedataframe <- read_delim("input3loci.mod.tab",
                           "\t", escape_double = FALSE, col_names = TRUE,
                           trim_ws = TRUE)
morphogroupsframe <- read_delim("pop.txt",
                                "\t", escape_double = FALSE, col_names = FALSE,
                                trim_ws = TRUE)

morphogroups <- factor(morphogroupsframe$X2,
                          levels = c(1,2,3),
                          labels=c("Ahya","Acyt","Abif"))
individuals <- factor(morphogroupsframe$X1)

objgenind <- df2genind(genotypedataframe, ploidy=2, sep="/", 
                 ind.names = individuals, pop = morphogroups)

# 4. Perform a basic PCA for exploration of genetic clustering:

sum(is.na(objgenind$tab)) # Check for missing data

X <- scaleGen(objgenind, NA.method="mean") # Replace missing data (NAs) and transform the allele frequencies

pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:5],main="PCA eigenvalues", col=heat.colors(5)) #Eigenvalues

s.class(pca1$li, pop(objgenind)) # simple PCA

col <- c("#2c5ab4", "#009933", "#FBA90A") # Adding some color
s.class(pca1$li, pop(objgenind),xax=1,yax=2, col=transp(col,.8),
        axesell=FALSE,
        cstar=1, cpoint=2, grid=FALSE)

# 5.Identify clusters: Options "kmeans" instead of "ward" and 'AIC' instead of 'BIC'

groups_3loci <- find.clusters(objgenind,method = "kmeans",
                           max.n.clust=10, stat = 'BIC')

# Check clusters
groups_3loci$size
groups_3loci$grp

# 6. Perform DAPC:
dapc1 <- dapc(objgenind, groups_3loci$grp, scale = TRUE)

scatter(dapc1) # basic scatterplot to check if number of PCs retained was enough 

# 7. Determine how many PCs needed to be retained
optim.a.score (dapc1)

# 8. Re-run DAPC after determining number of PCs to retain:
dapc2 <- dapc(objgenind, groups_3loci$grp, scale = TRUE)

# 9. Plot the DAPC:
col <- c("#d5af29","#ce882a","#7e6729")
scatter(dapc2, scree.da=FALSE, bg="white", pch=19,
        cell=2.5, cstar=1, col=col, solid=.7,
        cex=1.5,clab=0, leg=TRUE, posi.leg = NA)

# 10. Check if individual membership in each K cluster is specific:
assignplot(dapc2)

#------------------------------------------------------------------------------
# NOTE: Run Snapclust to explore Maximum-likelihood genetic clustering using EM algorithm

snapchoosek.aic <- snapclust.choose.k(10, objgenind,IC=AIC)
plot(1:10, snapchoosek.aic, xlab= "Number of clusters (K)",
     ylab = "AIC", type="b",pch=20, cex=3)

snapchoosek.kic <- snapclust.choose.k(10, objgenind,IC=KIC)
plot(1:10, snapchoosek.kic, xlab= "Number of clusters (K)",
     ylab = "KIC", type="b",pch=20, cex=3)

snapchoosek.bic <- snapclust.choose.k(10, objgenind,IC=BIC)
plot(1:10, snapchoosek.bic, xlab= "Number of clusters (K)",
     ylab = "BIC", type="b",pch=20, cex=3)

snapobj <- snapclust(objgenind, 3)
head(snapobj)
compoplot(snapobj)

