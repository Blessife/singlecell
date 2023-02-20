
# Load Library ------------------------------------------------------------


suppressWarnings({
  lapply(
    c(
      "dplyr",
      "Seurat",
      "HGNChelper",
      'GEOquery',
      'Matrix',
      'gridExtra',
      'tidyverse',
      'ggplot2'
    ),
    library,
    character.only = T
  )
})

# Fetch data  -------------------------------------------------------------

# set working directory 
workdir= "~/Retina" 
setwd(workdir)
accession_no = "GSE188280"
# download raw file to the directory
getGEOSuppFiles(accession_no)
# To open the tar file, access the folder with
paths=paste0(workdir,"/",accession_no)
setwd(paths)
# Open the tar file 
untar(paste0(accession_no,"_RAW.tar"))

# Create a vector list of all the samples
# samples=c("GSM5676873","GSM5676874","GSM5676875","GSM5676876","GSM5676877",
#           "GSM5676878","GSM5676879","GSM5676880","GSM5676881","GSM5676882",
#           "GSM5676883","GSM5676884")


retina_sample = c("GSM5676874","GSM5676876","GSM5676878","GSM5676880","GSM5676882","GSM5676884")
#choriod_sample = c("GSM5676873","GSM5676875","GSM5676877","GSM5676879","GSM5676881","GSM5676883")

samples = retina_sample
# Create sub folder  for each sample 


for (j in 1:length(samples)){
  folder<-dir.create(paste0(paths,'/',samples[j]))
}
# Move each sample file to the folder 

library(filesstrings)
dir.list= list.dirs(paths,recursive = FALSE)


for (i in 1:length(samples)){
  file_list<-list.files(paths, pattern  = paste0(samples[i],'_'))
  file.move(file_list, dir.list[i], overwrite=TRUE)
  
}

# Rename files in directory to fit into 10X format
for (i in dir.list){
  x=list.files(i, pattern = '.gz')
  for (a in x){
    b=gsub(".*_", "", a)
    setwd(i)
    file.rename(a,b)
  }
  
}

setwd(paths)

# Merge all together to prepare seurat object 

# get data location
dirs <- list.dirs(path = paths, recursive = F, full.names = F)

for(x in dirs){
  name <- x
  
  cts <- Read10X(data.dir = x)
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts, min.cells = 3, min.genes =200,
                                  project = 'AMD_retina'))
}



# Merge the dataset

merge_max = merge(GSM5676874, y=c(GSM5676876, GSM5676878, GSM5676880,GSM5676882,
                                  GSM5676884),
                    add.cell.ids=c("GSM5676874","GSM5676876","GSM5676878",
                     "GSM5676880","GSM5676882","GSM5676884"),
                   project = 'AMD_retina')


# # Save seurat object in file
#
save(merge_max,
      file=paste0(workdir,'/','merge_max.RData'))

# merge_max = merge(
#   GSM5676873,
#   y = c(
#     GSM5676874,
#     GSM5676875,
#     GSM5676876,
#     GSM5676877,
#     GSM5676878,
#     GSM5676879,
#     GSM5676880,
#     GSM5676881,
#     GSM5676882,
#     GSM5676883,
#     GSM5676884
#   ),
#   add.cell.ids = c(
#     "GSM5676873",
#     "GSM5676874",
#     "GSM5676875",
#     "GSM5676876",
#     "GSM5676877",
#     "GSM5676878",
#     "GSM5676879",
#     "GSM5676880",
#     "GSM5676881",
#     "GSM5676882",
#     "GSM5676883",
#     "GSM5676884"
#   ),
#   project = 'AMD'
# )



# save(merge_max, 
#      file=paste0(workdir,'/','mergeMax_obj.RData'))





# Perform Quality Check ---------------------------------------------------

# Prepare data for quality metrics

#load("mergeMax_obj.RData")

merge_max -> merge_max

View(merge_max@meta.data)

# create a sample column using the rowname of the metadata
merge_max$sample <- rownames(merge_max@meta.data)

# split sample column
merge_max@meta.data <-
  separate(
    merge_max@meta.data,
    col = 'sample',
    into = c('Sample_no', 'Barcode'),
    sep = '_'
  )
## Read phenotype data to R
phen=read.delim2('phendata.txt', sep = ',', header = T)

phen = phen %>%
  filter(Sample.Name %in% c("GSM5676874","GSM5676876","GSM5676878",
                            "GSM5676880","GSM5676882","GSM5676884"))

# integrate phendata with merge_max
for (y in 1:nrow(merge_max@meta.data)){
  z=merge_max@meta.data$Sample_no[y]
  merge_max@meta.data$sex[y] = phen$sex[phen$Sample.Name==z]
  merge_max@meta.data$source_name[y] = phen$source_name[phen$Sample.Name==z]
  merge_max@meta.data$phenotype[y] = phen$Phenotype[phen$Sample.Name==z]
}


view(merge_max@meta.data)
# compute quality assessment metrics  for the data

# calculate mitochondrial percentage
merge_max$mitoPercent <- PercentageFeatureSet(merge_max, pattern='^MT-')
merge_max$riboPercent <- PercentageFeatureSet(merge_max, pattern="^RP[LS]")

## Get a df from the slot meta data from the seurat object 

merge_max@meta.data -> meta.data


## Add some columns to the meta.data 

# Add number of genes per UMI for each cell to metadata
meta.data$log10GenesPerUMI <- log10(meta.data$nFeature_RNA) / log10(meta.data$nCount_RNA)

# Compute percent mito ratio
#meta.data$mitoRatio <- PercentageFeatureSet(object = merge_max, pattern = "^MT-")
meta.data$mitoRatio <- meta.data$mitoPercent / 100

merge_max@meta.data <- meta.data 

# save(merge_max,
#           file=paste0(workdir,'/','merge_max.RData'))


 
# Clone obj to perform some manipulation 
merge_max -> data.nomalat

#
#*****23762 genes and 30587 cells in 4 samples in its sparse matrix Only nonezeros values are stored
#Compute the percentage of the largest gene 
apply(data.nomalat@assays$RNA@counts,2,max) -> data.nomalat$largest_count

apply(data.nomalat@assays$RNA@counts,2,which.max)-> data.nomalat$largest_index

rownames(data.nomalat)[data.nomalat$largest_index] -> data.nomalat$largest_gene

100 * data.nomalat$largest_count / data.nomalat$nCount_RNA -> data.nomalat$percent.Largest.Gene

data.nomalat$largest_gene -> merge_max$largest_gene
data.nomalat$percent.Largest.Gene -> merge_max$percent.Largest.Gene
rm(data.nomalat)

merge_max@meta.data -> meta.data 

# Rename columns
meta.data <- meta.data %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
merge_max@meta.data <- meta.data


# Create .RData object to load at any time
save(merge_max, 
     file=paste0(paths,'/','merged_obj.RData'))

#################################################################
#load the merged_obj
library(ggplot2)

#load('merged_obj.RData')

merge_max@meta.data -> meta.data
## Exploring the quality of the data ..........................................


# Visualize the number of cell counts per sample
meta.data %>% 
  ggplot(aes(x=phenotype, fill=phenotype)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25)) +
  geom_text(aes(label= ..count..), stat = "count", vjust =2.0, size = 10)+
  ggtitle("Number of Retina cells")
# Visualize the number UMIs/transcripts per cell
meta.data %>% 
  ggplot(aes(color=phenotype, x=nUMI, fill= phenotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 23000)



# Visualize the distribution of genes detected per cell via histogram
meta.data %>% 
  ggplot(aes(color=phenotype, x=nGene, fill= phenotype)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


# filtering# Visualize the distribution of genes detected per cell via boxplot
meta.data %>% 
  ggplot(aes(x=phenotype, y=log10(nGene), fill=phenotype)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoPercent)) + 
  geom_point() + 
  scale_colour_gradient(low = "grey", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept =1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~phenotype)


meta.data %>% 
  ggplot(aes(color=phenotype, x=mitoPercent, fill=phenotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic()


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = phenotype, fill=phenotype)) +
  geom_density(alpha = 0.2) +
  theme_classic()

VlnPlot(merge_max, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))


## Cell level filtering 
merge_filtered <- subset(merge_max, subset = nUMI > 1000 &
                           nUMI < 15000 &
                           nGene > 500 &
                           mitoPercent > 0.5 &
                           log10GenesPerUMI > 0.85 &
                           riboPercent < 30)

# View quality after 
VlnPlot(merge_filtered, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))


merge_filtered

#merged_seurat



## gene level filtering 


# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = merge_filtered, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_AMDP <- CreateSeuratObject(filtered_counts, meta.data = merge_filtered@meta.data)



# Visualiz 
VlnPlot(filtered_AMDP, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))

#filtered_AMDP@meta.data->meta

#meta$source_name[grep('Chorod', meta$source_name)]='Choroid'

#merge_max@meta.data <- meta.data

filtered_AMDP@meta.data %>% 
  ggplot(aes(x=phenotype, fill=phenotype)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25)) +
  geom_text(aes(label= ..count..), stat = "count", vjust =2.0, size = 10)+
  ggtitle("Number of Retina cells")


# perform standard workflow steps to figure out if we see any batch effects --------

## Normalize the data 
merged_seurat_filtered <- NormalizeData(object = filtered_AMDP)
##Feature selection
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
## Scaling data
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
## Dimensional Reduction
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)

ElbowPlot(merged_seurat_filtered)
##Cell Sub population Identification using unsupervised clustering methods [!prior information based]
## Techniques for unsupervised clustering 
## (i) k-means; (ii) hierarchical clustering; (iii) density-based clustering; and (iv) graph-based clustering 
## 
## Other cluster methods are 
## single-cell consensus clustering (SC3) (Kiselev et al., 2017) 
## shared nearest neighbor (SNN) clustering algorithm (Satija et al., 2015) ....Seurat[findNeighbors]
##SC3 is an unsupervised approach that combines multiple clustering approaches, 
## which has a high accuracy and robustness in single-cell clustering.
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'phenotype') +
  ggtitle("Batch Effect")
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'source_name',
              cols = c('red','green'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)
#pca plot 
t1 <- DimPlot(merged_seurat_filtered, reduction = 'pca', group.by = 'phenotype')
t2 <- DimPlot(merged_seurat_filtered, reduction = 'pca', group.by = 'source_name',
              cols = c('red','green'))

## Get the list of most highly expressed genes
apply(merged_seurat_filtered@assays$RNA@data,1,mean) -> gene.expression

sort(gene.expression, decreasing = TRUE) -> gene.expression

head(gene.expression, n=10)

save(filtered_AMDP, file = 'AMD_Project.RData')



# Batch Effect Correction  ------------------------------------------------

#**************NORMALISATION & BATCH EFFECT CORRECTION*****************

# perform integration to correct for batch effects ------
#obj.list <- SplitObject(merged_seurat_filtered, split.by = 'phenotype')

obj.list <- SplitObject(filtered_AMDP, split.by = 'phenotype')

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}
# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)
# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'phenotype') + 
  ggtitle("Corrected Batch effect")
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'source_name',
              cols = c('red','green'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)


t3=DimPlot(seurat.integrated, reduction = 'pca', group.by = 'phenotype')
t4=DimPlot(seurat.integrated, reduction = 'pca', group.by = 'source_name',
           cols = c('red','green'))

## View all the pca plots
grid.arrange(t1, t3, t2, t4, ncol = 2, nrow = 2)

p1|p3
#save

save(seurat.integrated, file = 'AMDIntregate.RData')


# Dimension Reduction  ----------------------------------------------------
load("AMDIntregate.RData")


##################################################################
## Check the impact of cell cycle on the expression score 
##Predict the cell cycle of each cell 
## Get a new obj that allows for examine the cell cycle phase 

CellCycleScoring(seurat.integrated, 
                 s.features = cc.genes.updated.2019$s.genes, 
                 g2m.features = cc.genes.updated.2019$g2m.genes, 
                 set.ident = TRUE) -> cell.data

View(cell.data@meta.data)

as_tibble(cell.data[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()

## Another form of visual 
as_tibble(cell.data[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))


### Most of the cells in are in G1 phase which means they in cell growth phase. 

## Getting  the top 500 most varied genes 
FindVariableFeatures(
  cell.data, 
  selection.method = "vst", 
  nfeatures=500
) -> var.data


as_tibble(HVFInfo(var.data),rownames = "Gene") -> variance.data

variance.data %>% 
  mutate(hypervariable=Gene %in% VariableFeatures(var.data)
  ) -> variance.data

head(variance.data, n=10)

## Compare in graph
variance.data %>% 
  ggplot(aes(log(mean),log(variance),color=hypervariable)) + 
  geom_point() + 
  scale_color_manual(values=c("black","red"))


####**************************************DIMENSION REDUCTION ****************#######
## PCA to reduce the dimension 
RunPCA(cell.data,features=VariableFeatures(cell.data)) -> PCA.data

DimPlot(PCA.data, reduction = 'pca')

DimPlot(PCA.data,reduction="pca", group.by = "phenotype", 
        label = F, label.size = 3)

ElbowPlot(PCA.data)

DimHeatmap(PCA.data,dims=1:6, cells=500)


#### Performimg tSNE
8482 -> saved.seed
set.seed(saved.seed)



RunTSNE(
  PCA.data,
  dims=1:15,
  seed.use = saved.seed, 
  perplexity=100
) -> tdata

DimPlot(tdata,reduction = "tsne", pt.size = 1)+
  ggtitle("tSNE with Perplexity 100")


#save

save(tdata, file = "tdata.RData")


# Clustering --------------------------------------------------------------

##Defining Cluster using graph based methods
FindNeighbors(tdata,dims=1:10) -> kdata  ##--------kmeans = 20 
## Distance to the nearest 20 neigbours
kdata@graphs$integrated_snn[1:10,1:10]
##Segment the graph with FindCluster 
FindClusters(kdata,resolution = 0.01) -> fdata


###Cluster metadata 

DimPlot(fdata,reduction="pca",label = TRUE)+
  ggtitle("PC1 vs PC2 with Clusters")


DimPlot(fdata,reduction="pca", dims=c(4,9), label=TRUE)+
  ggtitle("PC4 vs PC9 with Clusters")


#### Performimg tSNE
8482 -> saved.seed
set.seed(saved.seed)


## check tSNE plot 
RunTSNE(
  fdata,
  dims=1:15,
  seed.use = saved.seed, 
  perplexity=100
) -> tfdata

DimPlot(tfdata,reduction="tsne",
        pt.size = 1, label = TRUE, label.size = 7)






## Compute the QC metric for thr clustered data 
VlnPlot(fdata,features="nCount_RNA")

VlnPlot(fdata,features="nFeature_RNA")

VlnPlot(fdata,features= 'mitoPercent', group.by = 'phenotype')

## Cell cycle per cluster 
fdata@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")



## Cell cycle per cluster - Normal
fdata@meta.data %>%
  filter(phenotype == "Normal") %>% 
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Normal Eye") -> cellCylNorm
## Cell cycle per cluster - Disease
fdata@meta.data %>%
  filter(phenotype == "Early AMD") %>% 
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Early AMD") -> cellcylAMD



cellCylNorm|cellcylAMD


fdata[[]] %>%
  group_by(seurat_clusters, phenotype) %>%
  count() %>%
  arrange(desc(n)) %>%
  group_by(seurat_clusters) %>%
  slice(1:2) %>%
  ungroup() %>%
  arrange(seurat_clusters, desc(n))


VlnPlot(fdata,features="MALAT1")

VlnPlot(fdata,features="percent.Largest.Gene", group.by = 'phenotype')

## Which gene is the largest 
fdata[[]] %>%
  filter(largest_gene != 'MALAT1') %>%
  group_by(seurat_clusters, largest_gene) %>%
  count() %>%
  arrange(desc(n)) %>%
  group_by(seurat_clusters) %>%
  slice(1:2) %>%
  ungroup() %>%
  arrange(seurat_clusters, desc(n)) -> data.gene



# tfdata@reductions$tsne@cell.embeddings %>%
#   as_tibble() %>%
#   add_column(seurat_clusters=tfdata$seurat_clusters, largest_gene=tfdata$largest_gene) %>%
#   filter(largest_gene %in% largest_genes_to_plot) %>%
#   ggplot(aes(x=tSNE_1, y=tSNE_2, colour=seurat_clusters)) +
#   geom_point() +
#   facet_wrap(vars(largest_gene))




# Cell Identification  ----------------------------------------------------


# Cell type assignment 
# Load the scType functions 
# load gene set preparation function

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

es.max = sctype_score(scRNAseqData = tfdata@assays$integrated@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Merge by cluster

cL_results = do.call("rbind", lapply(unique(tfdata@meta.data$seurat_clusters),
                                     function(cl){
                                       es.max.cl = sort(rowSums(es.max[,row.names(tfdata@meta.data[tfdata@meta.data$seurat_clusters==cl, ])]),
                                                        decreasing = !0)
                                       head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(tfdata@meta.data$seurat_clusters==cl)),
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

tfdata@meta.data$Cell_Identity = " "
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  tfdata@meta.data$Cell_Identity[tfdata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(tfdata, reduction = "tsne", label = T, repel = T, group.by = "Cell_Identity") + ggtitle("Retina cell communities")

DimPlot(tfdata, reduction = "umap", label = T, repel = T, group.by = c("Cell_Identity", 
                                                                       "phenotype"))

# custom function 
umap_theme <- function () {
  theme(axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), plot.title = element_text(hjust = 0.5))
        
}
 

DimPlot(tfdata, reduction = "umap", label = T, repel = T, 
        group.by = "Cell_Identity") + umap_theme() + NoLegend() +
  ggtitle("Retina cell types")
  
DimPlot(tfdata, reduction = "tsne", label = T, repel = T, 
        group.by = "Cell_Identity") + umap_theme()+
        ggtitle("Retina cell communities") + NoLegend()


# Findmarker --------------------------------------------------------------

## Identify differential expressed genes across normal and Earl AMD 

# Acomparative analysis  to look for the differences indused by EArly AMD 

#
DefaultAssay(tfdata) <- "RNA"

cone.markers <- FindConservedMarkers(tfdata, ident.1 = 7, grouping.var = "phenotype", verbose = FALSE)
head(cone.markers)

# Rename the ident for the obj based one cell type 

save(tfdata, file = "tfdata.RData")

Idents(tfdata) <- "seurat_clusters"


tfdatar <-
  RenameIdents(
    tfdata,
    `0` = "Rod bipolar cells",
    `1` = "Cone bipolar cells",
    `2` = "Horizontal cells",
    `3` = "Retinal ganglion cells",
    `4` = "Retinal ganglion cells",
    `5` = "Muller cells",
    `6` = "Rod bipolar cells",
    `7` = "Unknown",
    `8` = "Rod bipolar cells",
    `9` = "Microglial cells",
    `10` = "Glycinergic amacrine cells"
    
  )


DimPlot(tfdata, label = T)


DimPlot(tfdata, reduction = "tsne", label = T, repel = T, group.by = "Cell_Identity")

tfdatar$cell.phenotype <- paste(levels(tfdatar), tfdatar$phenotype, sep = "_")
tfdatar$cell.phenotype -> Idents(tfdatar)
Idents(tfdatar) <- "cell.phenotype"

# change the default assay to RNA

#DefaultAssay(tfdatar) = "RNA"
EarlyAMDrod <- FindMarkers(tfdatar, ident.1 = "Rod bipolar cells_Early AMD", ident.2 =  "Rod bipolar cells_Normal" , verbose = FALSE)
#test.use = "DESeq2")

EarlyAMDcone <- FindMarkers(tfdatar, ident.1 = "Cone bipolar cells_Early AMD", ident.2 =  "Cone bipolar cells_Normal" , verbose = FALSE)


head(EarlyAMDrod)

nrow(EarlyAMDrod)
head(EarlyAMDcone)





save(tfdatar, file = 'tfdatar.RData')


# co-expression network analysis packages with WGCNA -------------------------------------------------------------------

list.of.packages <- c("matrixStats", "Hmisc", "splines", "foreach", "doParallel",
                      "fastcluster", "dynamicTreeCut", "survival", "BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

BiocManager::install(c("GO.db", "preprocessCore", "impute"))

install.packages("WGCNA")
install.packages("devtools")

devtools::install_github("smorabit/hdWGCNA", ref="dev")


# Co-expression  ----------------------------------------------------------


# Load library ------------------------------------------------------------

suppressWarnings({lapply(c("Seurat",'WGCNA','ggplot2',
                           "dbplyr","cowplot","patchwork",
                           "AnnotationDbi", "GO.db", "preprocessCore", "impute", "igraph",
                           "tester"), library, character.only = T)})

library(hdWGCNA)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 16)

#load seurat obj

load("tfdata.RData")


# Set up seurat obj for WGCNA

seurat_obj <- SetupForWGCNA(
  tfdata,
  gene_select = "fraction",
  fraction = 0.05, # fraction of the cell that a gene need to be expressed to be included
  wgcna_name = "retina" #
)


