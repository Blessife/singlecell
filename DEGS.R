
# Installation  -----------------------------------------------------------
library(devtools)
install_github("immunogenomics/presto")

library(presto)
# Data Preparation  -------------------------------------------------------


# Standard Singlecell Processing  -----------------------------------------

Retina <- seurat.integrated %>% 
          FindVariableFeatures(selection.method = "vst", nfeatures=2000) %>% 
          ScaleData() %>% 
          RunPCA() %>% 
          FindNeighbors(dims = 1:15) %>% 
          FindClusters(resolution = 0.02) %>% 
          RunUMAP(dims=1:10)

DimPlot(Retina, reduction = "umap")

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

View(sctype_scores[, 1:3])


# Lets visualize the assigned cell 

Retina@meta.data$Cell_Identity = " "
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  Retina@meta.data$Cell_Identity[Retina@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(Retina, reduction = "tsne", label = T, repel = T, group.by = "Cell_Identity") + ggtitle("Retina cell communities")


# Create a vector list of the cell types

cell.types <- c("Rod bipolar cells","Muller cells","Cone bipolar cells","Horizontal cells",
                "Unknown","Cone photoreceptor cells","Microglial cells","Horizontal cells",
                "Pericytes","Astrocytes")

names(cell.types) <- levels(Retina)

Retina <- RenameIdents(Retina, cell.types)


# Custom function 

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
options(digits=2)
res <- wilcoxauc(Retina, 'seurat_clusters' , seurat_assay = 'RNA')

# get the topmarker 

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
