MEList = moduleEigengenes(t.Retina.d, colors = dynamicColors)
save(MEList, file = "GeneModules.RData")
load("GeneModules.RData")

MEs = MEList$eigengenes

datME <- MEs

dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")



# Diagnostic: displaying module heatmap and the eigengene -----------------

sizeGrWindow(8,7);
which.module="green"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")



# Finding the modules that relate to a clinical trait ---------------------
Trait = Retina$phenotype[Retina$cell_types == "Cone bipolar cells" | Retina$cell_types =="Cone photoreceptor cells"]

Trait = as.data.frame(Trait)      
Trait$y = 1
Trait$y[Trait$Trait == "Early AMD"] = 2
y = Trait$y

signif(cor(y,datME, use="p"),2)


# Module Significant ------------------------------------------------------

GS1=as.numeric(cor(y,t.Retina.d, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, dynamicColors, mean, na.rm=T)

sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,dynamicColors)



# Module Traits Correlation -----------------------------------------------

# convert sex to factor
SUB_RET_CELL <- subset(x = Retina, 
                       subset = cell_types == "Cone bipolar cells" | cell_types =="Cone photoreceptor cells")
# Convert category traits to factor
SUB_RET_CELL$phenotype <- as.factor(SUB_RET_CELL$phenotype)
SUB_RET_CELL$sex <- as.factor(SUB_RET_CELL$sex)

# list of traits to correlate
cur_traits <- c('nCount_RNA','nFeature_RNA','mitoRatio', "sex", "phenotype")

SUB_RET_CELL<- ModuleTraitCorrelation(
  SUB_RET_CELL,
  traits = cur_traits,
  group.by='cell_types'
)

