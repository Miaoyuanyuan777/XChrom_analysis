library(Seurat)
library(SeuratDisk)
rna <- readRDS('./COVID19_MHH50/scRNA/PBMC_scRNAseq.rds')
rna<-UpdateSeuratObject(rna)
rna_ms <- rna[,rna@meta.data$Severity == "mild"|rna@meta.data$Severity == "severe"]
rna_ms@assays$RNA@data <- rna_ms@assays$RNA@counts
SaveH5Seurat(rna_ms, filename = "./data/mild_severe_all.h5seurat")
Convert("./data/mild_severe_all.h5seurat", dest = "h5ad")
### 
# Adding scale.data from RNA as X
# Transfering meta.features to var
# Adding data from RNA as raw
# Transfering meta.features to raw/var

# When reading mild_severe_all.h5ad file in python, you need to execute ad=ad.raw. Only in this way can ad.x represent the original counts data; otherwise, it is considered scaled data