library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
set.seed(1)

addArchRThreads(threads = 42) 
addArchRGenome("hg38")

setwd('/Aging/Data/Anderson/Multiomics/ATAC/ArchR/')
DIR <- '/Aging/Data/Anderson/Multiomics/ATAC/ArchR/Output/'

#saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRProject", load = FALSE)
proj <- loadArchRProject(path = "/Aging/Data/Anderson/Multiomics/ATAC/ArchR/ArchRProject_Neurons", force = FALSE, showLogo = FALSE)
getAvailableMatrices(proj)

# Subset 2 cells only -- we only care about the rowRanges
proj = proj[1:2,]

#--------------------------------------------------------
#--------------------------------------------------------

# Step 1a
print('Step 1a')
# Extract binarized PeakMatrix
PeakMatrix <- getMatrixFromProject(
ArchRProj = proj,
useMatrix = "PeakMatrix",
binarize=TRUE)
print(PeakMatrix)
rm(proj)

regions <- as.data.frame(rowRanges(PeakMatrix))            # Genomic coordinates
regions <- regions[,1:3]
# Save for bedtools intersect with annotations 
write.table(regions, '/Aging/Data/Anderson/Multiomics/ATAC/ArchR/Output/Anderson_PeakMatrix_Reads.bed', 
            row.names=FALSE,quote =FALSE, col.names=FALSE,sep = "\t") 

