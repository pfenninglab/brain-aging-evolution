#!/usr/bin/env Rscript
R.version
library(Seurat)
library(SeuratDisk)

# just check which CON vs. SZ are in this seurat object....
base_dir <- "/singleCell/Ruzicka_snRNA-Seq"
data_dir <- paste0(base_dir, "/processed/")
seurat_object <- readRDS(paste0(data_dir, "21_10_11_seurat_raw_neurons_only.rds"))

#### FOLLOW NOTEBOOK AT: 21_10_04_Ruzicka_processing.ipynb
##### 1. remove NON NEURONS
toMatch <- c("Ex", "In")
neuron_names <- grep(paste(toMatch,collapse="|"),
                     unique(seurat_object$Celltype), value = TRUE)
neuron_names
seurat_object <- seurat_object[,seurat_object$Celltype %in% neuron_names]

print("Non neurons removed. Object dimensions:")
seurat_object

##### 2. remove PROBLEM PATIENTS
seurat_object <-seurat_object[,seurat_object$ID !="SZ57"]
pass <- which(table(seurat_object$ID) > 200)
print("length(pass)")
length(pass)
seurat_object <- seurat_object[,seurat_object$ID %in% names(pass)]
print("length(intersect(names(pass), seurat_object$ID))")
length(intersect(names(pass), seurat_object$ID))
seurat_object
gc()


#unique(seurat_object$Celltype)

SaveH5Seurat(seurat_object, 
    filename = paste0(base_dir,"/data/22_08_12_seurat_raw_FILTERED_neurons_only_controls_only.h5Seurat"))

Convert(paste0(base_dir,"/data/22_08_12_seurat_raw_FILTERED_neurons_only_controls_only.h5Seurat"), 
    dest = "h5ad")

### MAKE A BULK AND SAVE THAT
print("agg")
agg <- AggregateExpression(seurat_object,slot= "counts",
                          group.by = c("ID"),
                         return.seurat=TRUE) 
agg$ID <- rownames(agg[[]])
agg
SaveH5Seurat(agg, 
    filename = paste0(base_dir,
        "/data/22_08_12_seurat_raw_FILTERED_neurons_only_controls_only_pseudobulk.h5Seurat"))

Convert(paste0(base_dir,
        "/data/22_08_12_seurat_raw_FILTERED_neurons_only_controls_only_pseudobulk.h5Seurat"), 
    dest = "h5ad")

