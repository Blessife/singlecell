
load(url("https://amdproject-1.eu-central-1.linodeobjects.com/Choriodplus.RData"))



# Install Packages --------------------------------------------------------

list.of.packages <- c("matrixStats", "Hmisc", "splines", "foreach", "doParallel",
                      "fastcluster", "dynamicTreeCut", "survival", "BiocManager",
                      "cluster", "flashClust","reshape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

BiocManager::install(c("GO.db", "preprocessCore", "impute"))

install.packages("WGCNA")
#install.packages("presto")

# Library -----------------------------------------------------------------


suppressWarnings({lapply(c("Seurat",'WGCNA','ggplot2',
                           "dbplyr","cowplot","patchwork",
                           "AnnotationDbi", "GO.db", "preprocessCore", "impute", "igraph",
                           "tester","data.table","dplyr"), library, character.only = T)})


suppressWarnings({lapply(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel",
                           "fastcluster", "dynamicTreeCut", "survival", "BiocManager",
                           "cluster", "flashClust","reshape"), library, character.only = T )})


#library(WGCNA)
#library(Seurat)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

dim(Choriod)
DefaultAssay(Choriod)<-"integrated"
# Get matrix from seurat object 
ret.expr = GetAssayData(object = Choriod, slot = "data")
# change column names 
colnames(ret.expr) = Choriod$Barcode
gene.names = row.names(Choriod)
t.Choriod = t(ret.expr)

# Subset with cell specific markers 
dim(t.Choriod)

#genes = res.cone$feature
t.Choriod.s = as.matrix(t.Choriod)
celllist = Choriod$Barcode[Choriod$cell_types == "Fibroblasts"]
t.Choriod.d = t.Choriod.s[celllist,]

#save(res.wgcna, file = "cellspecificmarkerwgcna.RData")




# Choosing a soft-threshold to fit a scale-free topology --------

 
#enableWGCNAThreads(nThreads = 8)

powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(t.Choriod.d,dataIsExpr = TRUE,
                      powerVector = powers,
                      corFnc = cor,
                      corOptions = list(use = 'p'),
                      networkType = "signed hybrid")


# Visualize Result 

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


# Generating adjaceny and TOM similarity  ---------------------------------



softPower = 2;

#calculate the adjacency matrix
adj= adjacency(t.Choriod.d,type = "signed hybrid", power = softPower);

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(t.Choriod.d,networkType = "signed hybrid", TOMType = "signed", power = softPower);

colnames(TOM) =rownames(TOM) = colnames(t.Choriod.d)
dissTOM=1-TOM

# Module detection 
library(flashClust)
#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);

# Set the minimum module size
minModuleSize = 30;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)





plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                    main = "Choriod Cone Cells (photoreceptors + bipolar cells), n= 18294")


#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(t.Choriod.d[,restGenes], power = softPower)


colnames(diss1) =rownames(diss1) =colnames(t.Choriod.d)[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )

dynamicMods1 = cutreeDynamic(dendro = hier1,  method="tree", minClusterSize = minModuleSize);
plotDendroAndColors(hier1, dynamicMods1, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, cex.axis = 3,
                    main = "Choriod Fibroblast Cell, n= 10,888", cex.main = 3)


#set the diagonal of the dissimilarity to NA 
diag(diss1) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7,7)
TOMplot(diss1, hier1, as.character(dynamicMods1))

SubGeneNames <- colnames(t.Choriod.d)

#Extract modules

module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", 
              row.names=FALSE, col.names=FALSE,quote=FALSE)
}

#Quantify module similarity by eigengene correlation

MEList = moduleEigengenes(t.Choriod.d, colors = dynamicColors)
save(MEList, file = "GeneModules.RData")
load("GeneModules.RData")

MEs = MEList$eigengenes
#plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
MEs = orderMEs(MEs)

meInfo<- data.frame(rownames(t.Choriod.d), MEs)
colnames(meInfo)[1] = "SampleID"
# Intramodular connectivity 
#Calculation of (signed) eigengene-based connectivity, also known as module membership.
KMes <- signedKME(t.Choriod.d, MEs, outputColumnName = "MM.", corFnc = "bicor")

# complie into a module metadata table 
geneInfo = as.data.frame(cbind(colnames(t.Choriod.d), dynamicColors, KMes))

# Numbers of modules 
nmodules = length(unique(dynamicColors))

# Merge gen symbole column

colnames(geneInfo)[1] = "GeneSymbol"
colnames(geneInfo)[2] = "Initially.Assigned.Module.Color"

# save data 

write.csv(geneInfo, file = "geneInfoSigned.csv")


choriod_sub = subset(x = Choriod, subset = cell_types == "Fibroblasts")


#visualization
MEs -> PCvalues
plot_df = cbind(select(choriod_sub@meta.data, c(phenotype, cell_types)), PCvalues)

plot_df<- reshape2::melt(plot_df, id.vars = c("phenotype", "cell_types"))

plot_df$phenotype <- factor(plot_df$phenotype, levels = c("Early AMD", "Normal"))
colrs <- sub("ME", "", as.character(levels(plot_df$variable)))


passive = c("MEgreenyellow", "MEsalmon","MEmagenta", "MEred","MEyellow","MEblack","MEtan","MEgrey")


# Declare a function 
`%!in%` <- Negate(`%in%`)
active_plotdf = subset(plot_df, plot_df[,"variable"] %!in% passive)



p <- ggplot(plot_df, aes(x= variable, y = value, fill = phenotype)) +
  geom_boxplot(notch = FALSE) +
  RotatedAxis() + ylab("Module Eigengene") + xlab("") + 
  theme(panel.background = element_blank(),
        legend.text = element_text(size=30),
        legend.key.size = unit(2, 'cm'),
        axis.text=element_text(size=25)
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank()
    )

p

plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
w=2*nmodules; h=2*2;

pdf('figures/ME_Plot_bipolar_condition.pdf',width=w,height=h,useDingbats=F)
p + facet_wrap(cell_types~variable, scales='free', ncol=6 )
dev.off()

# Write data to disk
write.table(t.Choriod.s, file = "ChoriodMatrix.txt", sep="\t", 
            row.names=TRUE, col.names=TRUE,quote=FALSE)




# Gene Significance  --------------------------------------------------------

y = as.numeric(as.factor(choriod_sub$phenotype))

GS1=as.numeric(cor(y,t.Choriod.d, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, dynamicColors, mean, na.rm=T)


sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,dynamicColors)

# Hub genes ---------------------------------------------------------------

chooseTopHubInEachModule(
  t.Choriod.d, 
  dynamicColors, 
  omitColors = "grey", 
  power = 2, 
  type = "signed")



# Get Hub genes -----------------------------------------------------------

signif(cor(y,MEs, use="p"),2)

cor.test(y, MEs$MEgreen)

p.values = corPvalueStudent(cor(y,MEs, use="p"), nSamples = length(y))



# Intramodular connectivity -----------------------------------------------

ADJ1=abs(cor(t.Choriod.d,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)
head(Alldegrees1)


# Finding genes with high gene significance and high intramodular  --------

FilterGenes= abs(GS1)>0.05 & abs(KMes$MM.pink)>0.6
table(FilterGenes)

dimnames(data.frame(t.Choriod.d))[[2]][FilterGenes]->hub_gene
hub_gene = hub_gene[!is.na(hub_gene)]
hub_gene

write.table(hub_gene, file = "hub_gene_wgcna_60.txt")


Seu_hub = c("ATF3","BTG2","DUSP1","EGR1","FOS","FOSB","IER2","JUN","JUNB","NR4A1")


# Compare hubgenes --------------------------------------------------------

intersect(hub_gene, Seu_hub)

library(VennDiagram)  
grid.newpage()  
draw.pairwise.venn(area1 = length(hub_gene),                        # Create pairwise venn diagram
                   area2 = length(Seu_hub),
                   cross.area = length(intersect(hub_gene,Seu_hub)),
                   fill = c("red", "blue"),
                   lty = "blank",
                   category = c("Hub genes in WGCNA", "Hub genes in PPI network"),
                   cat.pos = c(345, 165),
                   cex = 2,
                   cat.cex = 1.5)



