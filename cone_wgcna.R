
# Load libraries  ---------------------------------------------------------

devtools::install_github('smorabit/scWGCNA')

library(Seurat)
library(Matrix)
library(tidyverse)
library(rliger)
load("Retinaplus.RData")
library(scWGCNA)


# Subset cone cell --------------------------------------------------------


retina_cone <- subset(Retina, cell_types %in% c("Cone photoreceptor cells", ))