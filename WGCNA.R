
# Install Packages --------------------------------------------------------

list.of.packages <- c("matrixStats", "Hmisc", "splines", "foreach", "doParallel",
                      "fastcluster", "dynamicTreeCut", "survival", "BiocManager",
                      "cluster", "flashClust","reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

BiocManager::install(c("GO.db", "preprocessCore", "impute"))

install.packages("WGCNA")


# Library -----------------------------------------------------------------

library(WGCNA)
library(Seurat)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

dim(Retina)
# Get matrix from seurat object 
ret.expr = GetAssayData(object = Retina, slot = "data")
# change column names 
colnames(ret.expr) = Retina$Barcode
gene.names = row.names(Retina)
t.Retina = t(ret.expr)
