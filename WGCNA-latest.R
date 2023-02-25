
# Install Packages --------------------------------------------------------

list.of.packages <- c("matrixStats", "Hmisc", "splines", "foreach", "doParallel",
                      "fastcluster", "dynamicTreeCut", "survival", "BiocManager",
                      "cluster", "flashClust","reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

BiocManager::install(c("GO.db", "preprocessCore", "impute"))

install.packages("WGCNA")


# Library -----------------------------------------------------------------


suppressWarnings({lapply(c("Seurat",'WGCNA','ggplot2',
                           "dbplyr","cowplot","patchwork",
                           "AnnotationDbi", "GO.db", "preprocessCore", "impute", "igraph",
                           "tester","hdWGCNA","data.table","presto","dplyr"), library, character.only = T)})

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

# Subset with cell specific markers 
dim(t.Retina)

#genes = res.cone$feature
t.Retina.s = as.matrix(t.Retina)
celllist = Retina$Barcode[Retina$cell_types == "Cone bipolar cells" | Retina$cell_types =="Cone photoreceptor cells"]
t.Retina.d = t.Retina.s[celllist,]

save(res.wgcna, file = "cellspecificmarkerwgcna.RData")

# Checking for o

# # Re-cluster samples
# sampleTree2 = hclust(dist(datExpr), method = "average")
# # Convert traits to a color representation: white means low, red means high, grey means missing entry
# traitColors = numbers2colors(datTraits, signed = FALSE);
# # Plot the sample dendrogram and the colors underneath.
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = names(datTraits),
#                     main = "Sample dendrogram and trait heatmap")


# Choosing a soft-threshold to fit a scale-free topology to the network 
enableWGCNAThreads(nThreads = 8)

powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(t.Retina.d,dataIsExpr = TRUE,
                      powerVector = powers,
                      corFnc = cor,
                      corOptions = list(use = 'p'),
                      networkType = "signed hybrid")
  
  
# Visualise Result 

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Generating adjaceny and TOM similarity 

softPower = 2;

#calculate the adjacency matrix
adj= adjacency(t.Retina.d,type = "signed hybrid", power = softPower);

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(t.Retina.d,networkType = "signed hybrid", TOMType = "signed", power = softPower);

colnames(TOM) =rownames(TOM) = colnames(t.Retina.d)
dissTOM=1-TOM

# Module detection 
library(flashClust)
#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);

# Set the minimum module size
minModuleSize = 20;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)





plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                    main = "Retina Cone Cells (photoreceptors + bipolar cells), n= 18294")


#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(t.Retina.d[,restGenes], power = softPower)


colnames(diss1) =rownames(diss1) =colnames(t.Retina.d)[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, 
                    main = "Retina Cone Cells (photoreceptors + bipolar cells), n= 18294")


#set the diagonal of the dissimilarity to NA 
diag(diss1) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7,7)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))

SubGeneNames <- colnames(t.Retina.d)

#Extract modules

module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", 
              row.names=FALSE, col.names=FALSE,quote=FALSE)
}

#Quantify module similarity by eigengene correlation

MEList = moduleEigengenes(t.Retina.d, colors = dynamicColors)
save(MEList, file = "GeneModules.RData")
load("GeneModules.RData")

MEs = MEList$eigengenes
#plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
MEs = orderMEs(MEs)
meInfo<- data.frame(rownames(t.Retina.d), MEs)
colnames(meInfo)[1] = "SampleID"
# Intramodular connectivity 
#Calculation of (signed) eigengene-based connectivity, also known as module membership.
KMes <- signedKME(t.Retina.d, MEs, outputColumnName = "KME", corFnc = "bicor")

# complie into a module metadata table 
geneInfo = as.data.frame(cbind(colnames(t.Retina.d), dynamicColors, KMes))

# Numbers of modules 
nmodules = length(unique(dynamicColors))

# Merge gen symbole column

colnames(geneInfo)[1] = "GeneSymbol"
colnmaes(geneInfo)[2] = "Initially.Assigned.Module.Color"

# save data 

write.csv(geneInfo, file = "geneInfoSigned.csv")

#visualization
MEs -> PCvalues
plot_df = cbind(select(retina_cone@meta.data, c(phenotype, cell_types)), PCvalues)

plot_df<- reshape2::melt(plot_df, id.vars = c("phenotype", "cell_types"))

plot_df$cell_types <- factor(plot_df$cell_types, levels = c("Cone photoreceptor cells", "Cone bipolar cells"))
 colrs <- sub("ME", "", as.character(levels(plot_df$variable)))
 
 photo_cell = plot_df[plot_df$cell_types == "Cone photoreceptor cells",]
 
 bipolar_cell = plot_df[plot_df$cell_types == "Cone bipolar cells",]
 
 p <- ggplot(bipolar_cell, aes(x= variable, y = value, fill = phenotype)) +
   geom_boxplot(notch = FALSE) +
   RotatedAxis() + ylab("Module Eigengene") + xlab("") + 
   theme(
     axis.text.x = element_blank(),
     axis.ticks.x = element_blank()
     
   )
 w=2*nmodules; h=2*2;
 
 pdf('figures/ME_Plot_bipolar_condition.pdf',width=w,height=h,useDingbats=F)
 p + facet_wrap(cell_types~variable, scales='free', ncol=6 )
 dev.off()

# Write data to disk
write.table(t.Retina.s, file = "RetinaMatrix.txt", sep="\t", 
            row.names=TRUE, col.names=TRUE,quote=FALSE)










# Data Prepartion  --------------------------------------------------------

meanExpressionByArray = apply(t.Retina.d,1,mean, na.rm=T)
NumberMissingByArray = apply(is.na(data.frame(t.Retina.d)),1,sum)

sizeGrWindow(9, 5)
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        cex.names = 0.7)

# KeepArray= NumberMissingByArray<500
# table(KeepArray)
# datExpr=datExpr[KeepArray,]
# y=y[KeepArray]
# ArrayName[KeepArray]

y= Retina$phenotype

# Standard gene screening based on marginal correlation 
# Using Marginal Pearson to relate genes to phenotype 


