
# Load Seurat obj ---------------------------------------------------------


wget(url = c("https://amdproject-1.eu-central-1.linodeobjects.com/AMDIntregate.RData",
               "https://amdproject-1.eu-central-1.linodeobjects.com/Retina.RData"))


# Data preparation  -------------------------------------------------------
library('Matrix')
load("Retinaplus.RData")
# Save normalised counts
expr_mat <- as.matrix(Retina@assays$RNA@data)

writeMM(Retina@assays$RNA@data, file = "matrix.mtx")


write(x=rownames(Retina@assays$RNA@data), file = "features.tsv")
write(x = colnames(Retina@assays$RNA@data), file = "barcodes.tsv") 

# Metadata 


# Filter with 200 gene that are specific cell-type

#expr_mat = expr_mat[res.down$feature,]
# Select based on phenotype,
#fiter.vec = Retina$phenotype=="Early AMD"
#expr_mat = expr_mat[,fiter.vec]

# Subset of seurat object where phenotype is AMD 
SUB_OBJ_AMD <- subset(x = Retina, subset = phenotype == "Early AMD")

expr_mat <- as.matrix(SUB_OBJ_AMD@assays$RNA@data)
expr_mat -> b_expr
# save to file 
colnames(expr_mat) <- SUB_OBJ_AMD$Barcode
write.table(expr_mat, file="expr.txt", quote=F, sep="\t")
write.table(data.frame(Cell_bc = SUB_OBJ_AMD$Barcode, Cell_type = SUB_OBJ_AMD$cell_types),
            file="metadata.txt", row.names=F, quote=F, sep="\t")
