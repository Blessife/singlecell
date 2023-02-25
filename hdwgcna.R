load("Retinaplus.RData")
Retina -> retina 

DefaultAssay(retina) <- "RNA"
VariableFeatures(retina)<- rownames(retina[["RNA"]])


# Set up for wgcna --------------------------------------------------------


seurat_obj <- SetupForWGCNA(
  retina,
  gene_select = "variable", # the gene selection approach
  #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "retina" # the name of the hdWGCNA experiment
)


# Create Metacell ---------------------------------------------------------

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cell_types", "phenotype"), # specify the columns in seurat_obj@meta.data to group by
  k = 35, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cell_types' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

# Set up Expression martix ------------------------------------------------

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Cone photoreceptor cells","Cone bipolar cells"),
  # the name of the group of interest in the group.by column
  group.by='cell_types', # the metadata column containing the cell type info. 
  #This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)


# Select soft power threshold ---------------------------------------------

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)


# Construct co-expression Network  ----------------------------------------

seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=5,
  setDatExpr=FALSE,
  tom_name = 'Cone' # name of the topoligical overlap matrix written to disk
)


# Visualize ---------------------------------------------------------------

PlotDendrogram(seurat_obj, main='Cone Cells hdWGCNA Dendrogram')

TOM <- GetTOM(seurat_obj)



# Comute harmonise module eigengene ---------------------------------------

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))


save(seurat_obj, file = "scale_hdwgcna.RData")

load("scale_hdwgcna.RData")
# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="phenotype"
)


# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


# Compute connectivity  ---------------------------------------------------

library(qlcMatrix)
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_types', 
  group_name = c("Cone photoreceptor cells","Cone bipolar cells")
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Cone-M"
)

# get the module assignment table:
modules <- GetModules(seurat_obj)

# show the first 6 columns:
head(modules[,1:6])

# Get hub genes
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

save(hub_df, file = "hub_gene.RData")

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)


# Seurat plotting function 
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'cell_types')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') + 
  ggtitle("Module expression in cell types") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))

# plot output
p

save(seurat_obj, file = "ready_obj_visual.RData")

# Plot INH-M4 hME using Seurat VlnPlot function
p <- VlnPlot(
  seurat_obj,
  features = 'Cone-M7',
  group.by = 'cell_types',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p <- p + xlab('') + ylab('hME') + NoLegend()

# plot output
p



# Module trait correlation ------------------------------------------------

# convert sex to factor
seurat_obj$sex <- as.factor(seurat_obj$sex)

# convert age_death to numeric
seurat_obj$phenotype <- as.factor(seurat_obj$phenotype)

# list of traits to correlate
#cur_traits <- c('nCount_RNA','nFeature_RNA','mitoRatio', "sex", "phenotype")
cur_traits <- c("sex","phenotype")
seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by='cell_types'
)

# Result
# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(seurat_obj)

names(mt_cor)


names(mt_cor$cor)


PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 5,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)



# Differential Expressed  -------------------------------------------------


seurat_grp <- subset(seurat_obj, 
                     subset = cell_types == "Cone bipolar cells" | cell_types =="Cone photoreceptor cells")

seurat_grp@meta.data -> meta.data
meta.data2 = meta.data[!duplicated(meta.data$Barcode),]
rownames(meta.data2) <- meta.data2$Barcode
seurat_grp@meta.data <- meta.data2


# groupA <- seurat_grp@meta.data %>% 
#             subset(phenotype == "Early AMD")  
# groupB <- seurat_grp@meta.data %>% 
#               subset(phenotype == "Normal")
#   
# group1 = groupA$Barcode
# group2 = groupB$Barcode

group1 <- seurat_grp@meta.data %>% subset(phenotype == "Early AMD") %>% rownames
group2 <- seurat_grp@meta.data %>% subset(phenotype == "Normal") %>% rownames
 
head(group1)
head(group1)

DMEs <- FindDMEs(
  seurat_grp,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='Cone'
)

head(DMEs)

