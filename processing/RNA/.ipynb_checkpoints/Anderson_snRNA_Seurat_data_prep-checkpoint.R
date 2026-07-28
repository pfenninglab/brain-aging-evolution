library(Seurat)
library(SeuratDisk)
library(SoupX)
library(DropletUtils)
library(patchwork)
library(DoubletFinder)
library(dplyr)

### directories
raw_data_dir = '/Aging/Data/Anderson/Multiomics/Raw/'
data_dir = '/Aging/Data/Anderson/snRNA_Seq/'

# Load the obj dataset
obj.data <- Read10X_h5("/Aging/Data/Anderson/Multiomics/Raw/GSE214979_filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
#obj <- CreateSeuratObject(counts = obj.data$`Gene Expression`, project = "obj3k", min.cells = 3, min.features = 200)
obj <- CreateSeuratObject(counts = obj.data$`Gene Expression`)
obj

## Ambient RNA removal by SoupX
# list of samples 
dir = '/Aging/Data/Anderson/Multiomics/Raw/'
files = list.files(dir, pattern = "_feature_bc_matrix.h5")
print(files)

for(sample in files){
  ## read in the filtered counts and raw counts
  toc = Seurat::Read10X_h5(paste0(dir,"GSE214979_filtered_feature_bc_matrix.h5"))
  toc = toc$`Gene Expression`
  tod = Seurat::Read10X_h5(paste0(dir,"GSE214979_unfiltered_feature_bc_matrix.h5"))
  tod = tod$`Gene Expression`
  
  ## routine Seurat clustering for SoupX, will be redone later
  ## https://satijalab.org/seurat/articles/essential_commands.html
  print(paste('Processing counts for SoupX.'))
  obj <- CreateSeuratObject(counts = toc) %>%
    NormalizeData(verbose = F) %>%
    FindVariableFeatures(verbose = F) %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(algorithm = 2, resolution = 0.5, verbose = F)
  
  ## give SoupX the filtered counts, raw counts, and highly required basic clustering
  print(paste('Estimating ambient RNA w/ SoupX.'))
  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(obj$seurat_clusters, colnames(obj)))
  sc = autoEstCont(sc, doPlot = F)
  out = adjustCounts(sc, roundToInt=TRUE)
  
  print("Writing Output.")
    
  ## write 10X output
  dir_soupx = paste0(dir,"SoupX/")
  #write10xCounts(dir_soupx, out, barcodes = colnames(out),version = '3')
    
  obj = CreateSeuratObject(counts = out)
  #saveRDS(obj, paste0(dir_soupx,"GSE214979_feature_bc_matrix_seurat_soupx.rds"))
}
obj

## Doublet detection 
## QC - filtering cells 

# load dataset and create a list
obj.list <- c(obj = obj) 

for (p in 1:length(obj.list)) {
  print(names(obj.list[p]))
  obj = obj.list[[p]]
  
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  #seu_obj <- CreateSeuratObject(obj)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 5)
  obj <- NormalizeData(obj, verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F) %>%
    FindNeighbors(verbose = F) %>%
    FindClusters(algorithm = 2, resolution = 0.5, verbose = F) %>%
    RunUMAP(dims = 1:10, verbose = F)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:10, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep.stats_obj)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(obj@meta.data[["seurat_clusters"]])           ## ex: annotations <- obj@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  #obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  name <- names(obj.list[p])
  assign(name, obj)
}

#obj <- subset(obj, subset = DF.classifications_0.25_0.09_122 == "Singlet")
#obj@assays[["RNA"]]@data <- obj@assays[["RNA"]]@counts
#obj
saveRDS(obj, paste0(dir_soupx,"GSE214979_feature_bc_matrix_seurat_soupx_df.rds"))


####################################
## Add metadata 
### directories
raw_data_dir = '/Aging/Data/Anderson/Multiomics/Raw/'
dir = '/Aging/Data/Anderson/Multiomics/SoupX/'
obj <- readRDS(paste0(dir,"GSE214979_feature_bc_matrix_seurat_soupx_df.rds"))

# Add metadata 
metadata = read.csv(paste0(raw_data_dir,"GSE214979_cell_metadata.csv"))
metadata$barcode = metadata$X
obj@meta.data$barcode = rownames(obj@meta.data)
obj@meta.data = left_join(obj@meta.data, metadata, by = "barcode")
rownames(obj@meta.data) = obj@meta.data$barcode
#obj@meta.data
#colnames(obj@meta.data)
## Doublets removal 
obj <- subset(obj, DF.classifications_0.25_0.09_6948 == "Singlet")
#obj
# Subset neurons
obj <- subset(obj, subset = (predicted.id == "Excitatory" | predicted.id == "Inhibitory"))
obj@assays[["RNA"]]@data <- obj@assays[["RNA"]]@counts
VariableFeatures(obj) <- NULL
obj@reductions[["pca"]] <- NULL
obj@reductions[["umap"]] <- NULL
obj

#saveRDS(obj, paste0(data_dir,"GSE214979_filtered_feature_bc_matrix_processed_seurat.rds"))
#obj <- readRDS(paste0(data_dir,"GSE214979_filtered_feature_bc_matrix_processed_seurat.rds"))

####################################
## Processing 
# Subset Ctrl
obj_ctrl = subset(obj, subset = Status == "Ctrl")
obj_ctrl = NormalizeData(obj_ctrl, verbose = F) %>%
FindVariableFeatures(nfeatures = nrow(obj_ctrl), verbose = F) %>%
ScaleData(verbose = F) 
#%>%
#RunPCA(verbose = F) %>%
#FindNeighbors(verbose = F) %>%
#FindClusters(algorithm = 2, resolution = 0.5, verbose = F)
obj_ctrl

# Subset AD
obj_ad = subset(obj, subset = Status == "AD")
obj_ad = NormalizeData(obj_ad, verbose = F) %>%
FindVariableFeatures(nfeatures = nrow(obj_ad), verbose = F) %>%
ScaleData(verbose = F) 
#%>%
#RunPCA(verbose = F) %>%
#FindNeighbors(verbose = F) %>%
#FindClusters(algorithm = 2, resolution = 0.5, verbose = F)
obj_ad

## Converting to h5ad
data_dir = '/Aging/Data/Anderson/Multiomics/Processed/'
## Ctrl
# Converting factor variables to character [read in h5ad]
i <- sapply(obj_ctrl@meta.data, is.factor)
obj_ctrl@meta.data[i] <- lapply(obj_ctrl@meta.data[i], as.character)
setwd(data_dir)
SaveH5Seurat(obj_ctrl, filename = "GSE214979_filtered_feature_bc_matrix_processed_seurat_neurons_ctrl.h5Seurat")
Convert("GSE214979_filtered_feature_bc_matrix_processed_seurat_neurons_ctrl.h5Seurat", dest = "h5ad")
## AD
# Converting factor variables to character [read in h5ad]
i <- sapply(obj_ad@meta.data, is.factor)
obj_ad@meta.data[i] <- lapply(obj_ad@meta.data[i], as.character)
setwd(data_dir)
SaveH5Seurat(obj_ad, filename = "GSE214979_filtered_feature_bc_matrix_processed_seurat_neurons_AD.h5Seurat")
Convert("GSE214979_filtered_feature_bc_matrix_processed_seurat_neurons_AD.h5Seurat", dest = "h5ad")