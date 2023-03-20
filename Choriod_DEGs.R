

# Load Seurat Obj ---------------------------------------------------------

load(url("https://amdproject-1.eu-central-1.linodeobjects.com/Choriodplus.RData"))

choriod_sub = subset(x = Choriod, subset = cell_types == "Fibroblasts")

# Rename the ident for the obj based one cell type 

choriod_sub$cell.phenotype <- paste(levels(choriod_sub), choriod_sub$phenotype, sep = "_")
choriod_sub$cell.phenotype -> Idents(choriod_sub)
Idents(choriod_sub) <- "cell.phenotype"

# change the default assay to RNA

DefaultAssay(choriod_sub) = "RNA"
choriod_fib <- FindMarkers(choriod_sub, ident.1 = "Fibroblasts_Early AMD", 
                           ident.2 =  "Fibroblasts_Normal" , verbose = FALSE)
#test.use = "DESeq2")

#Filter for DEGs
attach(choriod_fib)
Fib.DEGs <- row.names(choriod_fib)[abs(avg_log2FC)>2 & p_val_adj < 0.01]

write.table(Fib.DEGs, file = "Fib_DEGs.txt", sep = "/t")
# Compare with WGCNA 

#install.packages("ggVennDiagram")
library(ggVennDiagram)

# List of items
x <- list(A = hub_gene, B= Fib.DEGs)

# 2D Venn diagram
ggVennDiagram(x) 


intersect(hub_gene, Fib.DEGs)

library(VennDiagram)
#venn.diagram(DEG.list, filename = "VennDiagram.png")


display_venn <- function(x,...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, main.cex = 2.0, sub.cex = 2.0, main.fontface="bold", ...)
  grid.draw(venn_object)
}

DEG.list=list(hub_gene, Fib.DEGs)

display_venn(
  DEG.list,
  category.names = c("WGCNA_HubGenes", "Seurat_DEGs"),
  fill = c("#999999", "#E69F00"),
  main= "Gene Validation",
  cat.fontface = "bold",
  cex=2.0,
  cat.cex=2.0) 
