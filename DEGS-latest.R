
# Installation  -----------------------------------------------------------
library(devtools)
install_github("immunogenomics/presto")

#collectGarbage()

# Load Library ------------------------------------------------------------

suppressWarnings({lapply(c("Seurat",'WGCNA','ggplot2',
                           "dbplyr","cowplot","patchwork",
                           "AnnotationDbi", "GO.db", "preprocessCore", "impute", "igraph",
                           "tester","hdWGCNA","data.table","presto","dplyr"), library, character.only = T)})

library(hdWGCNA)
library(presto)
library(data.table)
library(dplyr)
# Download seurat obj  -------------------------------------------------------
if (FALSE) {
  wget(url = c("https://amdproject-1.eu-central-1.linodeobjects.com/AMDIntregate.RData",
               "https://amdproject-1.eu-central-1.linodeobjects.com/Retina.RData"))
}

# Standard Single-cell Processing  -----------------------------------------

load("AMDIntregate.RData")
Retina <- seurat.integrated %>% 
          FindVariableFeatures(selection.method = "vst", nfeatures=2000) %>% 
          ScaleData() %>% 
          RunPCA() %>% 
          FindNeighbors(dims = 1:15) %>% 
          FindClusters(resolution = 0.02) %>% 
          RunUMAP(dims=1:10)

DimPlot(Retina, reduction = "umap")


# Cluster Annotation ------------------------------------------------------

library(HGNChelper)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")

# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# To prepare the gene set lets import the db

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Eye" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# Lets assign the cell  types to each cluster 

es.max = sctype_score(scRNAseqData = Retina@assays$integrated@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Merge by cluster

cL_results = do.call("rbind", lapply(unique(Retina@meta.data$seurat_clusters),
                                     function(cl){
                                       es.max.cl = sort(rowSums(es.max[,row.names(Retina@meta.data[Retina@meta.data$seurat_clusters==cl, ])]),
                                                        decreasing = !0)
                                       head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Retina@meta.data$seurat_clusters==cl)),
                                            10)
                                       
                                     }))

sctype_scores = cL_results %>% 
  group_by(cluster) %>% 
  top_n(n=1,wt =scores)

# set low-confident (low ScType score) clusters to "unknown"


sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells /
                     4] = "Unknown"
# View result
View(sctype_scores[, 1:3])


# Lets visualize the assigned cell 

Retina@meta.data$Cell_Identity = " "
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  Retina@meta.data$Cell_Identity[Retina@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(Retina, reduction = "umap", label = T, repel = T, group.by = "Cell_Identity") + ggtitle("Retina cell communities")
backretina <- Retina
Retina -> retina
# Create a vector list of the cell types

cell.types <- c("Rod bipolar cells","Muller cells","Cone bipolar cells","Horizontal cells",
                "Unknown","Cone photoreceptor cells","Microglial cells","Horizontal cells",
                "Pericytes","Astrocytes")

names(cell.types) <- levels(retina)

retina <- RenameIdents(retina, cell.types)
# Add the cell ident to the metadata
Idents(retina) -> cell_types

Retina <- retina
Retina$cell_types <- cell_types
save(Retina, file = "Retinaplus.RData")

# Custom function for visualisation

umap_theme <- function () {
  theme(axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), plot.title = element_text(hjust = 0.5))
}

DimPlot(Retina, reduction = "umap", label = T, repel = T, ) + umap_theme() + NoLegend() +
  ggtitle("Retina cell types")
  


# Perform DEGs with wilcoxauROC -------------------------------------------
#load("Retina.RData")
options(digits=2)
res <- wilcoxauc(Retina, 'cell_types' , seurat_assay = 'RNA')

# get the topmarker 

res %>% 
  filter(logFC > log(1.2) & pct_in > 50 & padj < 0.05) %>% 
  group_by(group) %>% 
  arrange(desc(logFC), .by_group=T) %>% 
  top_n(n=150, wt= logFC) %>%  
  filter(group == "Cone photoreceptor cells" | group == "Cone bipolar cells") -> res.cone

save(res.down, file = "cellspecificmarker.RData")

marker = top_markers(res, n=6, auc_min = 0.5, pct_in_min = 50)

all_markers<- marker %>%
                select(-rank) %>% 
                unclass() %>% 
                stack() %>%
                pull(values) %>%
                unique() %>%
                .[!is.na(.)]

p<- DotPlot(object = Retina, features = all_markers)
p

theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



save(Retina, file = "Retina.RData")

# Visualize specific genes 
p1= FeaturePlot(Retina, c("MEG3", "TTR", "LMO4"), ncol = 1)
p2 =VlnPlot(Retina, features = c("MEG3", "TTR", "LMO4"), pt.size = 0)

p1 | p2 


# Visualize top genes in each cluster -------------------------------------

colors <- c("red","blue")
res.dotplot$feature <- as.factor(res.dotplot$feature)
DotPlot(object = Retina, features = unique(res.dotplot$feature), cols = colors) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(size = 15)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  ggtitle("Cell Specific Marker genes")




# WGCNA -------------------------------------------------------------------


cellmarkers = res.down$feature

library(hdWGCNA)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 16)


# Set up seurat obj for WGCNA

seurat_obj <- SetupForWGCNA(
  Retina,
  features = cellmarkers,
  wgcna_name = "retina" #
)

# construct metacells  in each group
MetacellsByGroups(
  seurat_obj,
  group.by = c("cell_types","phenotype"),
  ident.group = "cell_types",
  k = 30,
  reduction = "umap",
  assay = "RNA",
  cells.use = NULL,
  slot = "data",
  mode = "average",
  min_cells = 100,
  max_shared = 15,
  target_metacells = 1000,
  max_iter = 5000,
  verbose = FALSE,
  wgcna_name = "retina"
)-> seurat_objs




# normalize metacell expression matrix:
seurat_objs <- NormalizeMetacells(seurat_objs)


# Visualize the metacell created  ------------------------------------------

# standard seurat work flow 

seurat_objs <- seurat_objs %>% 
                  NormalizeMetacells()
seurat_objs <- ScaleMetacells(seurat_objs, features=VariableFeatures(seurat_objs))
seurat_objs <- RunPCAMetacells(seurat_objs, features=VariableFeatures(seurat_objs))
seurat_objs <- RunHarmonyMetacells(seurat_objs, group.by.vars='phenotype')
seurat_objs <- RunUMAPMetacells(seurat_objs, reduction='harmony', dims=1:15)


p1 <-
  DimPlotMetacells(seurat_objs, group.by = 'cell_types') + umap_theme() +
  ggtitle("cell_types")
p2 <-
  DimPlotMetacells(seurat_objs, group.by = 'phenotype') + umap_theme() +
  ggtitle("phenotype")

p1 | p2



# Co-expression Network analysis ------------------------------------------

# set up expression matrix
seurat_objs <- SetDatExpr(
  seurat_objs,
  group_name = 'Cone bipolar cells', # the name of the group of interest in the group.by column
  group.by= 'cell_types', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

seurat_objs@assays$RNA@meta.features

seurat_objs <- TestSoftPowers(
  seurat_objs,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
