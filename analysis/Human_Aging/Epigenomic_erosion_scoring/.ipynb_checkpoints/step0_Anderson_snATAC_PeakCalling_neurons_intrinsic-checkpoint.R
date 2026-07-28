#-----------------------------------------------------------------------------------------------------------------
# Load Required Libraries 
#-----------------------------------------------------------------------------------------------------------------
library(ArchR)          # Framework for single-cell ATAC-seq analysis
library(tidyverse)      # Collection of data manipulation packages
library(here)           # Project-relative file path management
library(limma)          # Differential expression analysis
library(rtracklayer)    # For reading and writing genomic data
library(parallelly)
library(future)
library(rhdf5)
#neg binom
library(MASS)
library(broom)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = availableCores())
addArchRGenome('hg38')

#future options
plan("multicore", workers = 8)
options(future.globals.maxSize = 300 * 1024^3)
options(future.rng.onMisuse = 'ignore')

#-----------------------------------------------------------------------------------------------------------------
# Load ArchR Project
#-----------------------------------------------------------------------------------------------------------------
#PROJDIR='/Aging/Data/Anderson/Multiomics/ATAC/ArchR/ArchRProject/'
#PROJDIR1='/Aging/Data/Anderson/Multiomics/ATAC/ArchR/ArchRProject_Neurons/'
PROJDIR2='/Aging/Data/Anderson/Multiomics/ATAC/ArchR/ArchRProject_Neurons_IntrinsicAging/'

## Load ArchR Project
setwd(PROJDIR2)
proj <- loadArchRProject(path = PROJDIR2, force = FALSE, showLogo = FALSE)
getAvailableMatrices(proj)

#-----------------------------------------------------------------------------------------------------------------
# Prep ArchR Project
#-----------------------------------------------------------------------------------------------------------------

## Adding coldata 
df = as.data.frame(proj$cellNames)
colnames(df) <- 'X'
# sort new column to the order of cellNames in ArchR project
metadata = read.csv('/home/gabdelha/Data/Anderson/Multiomics/Raw/GSE214979_cell_metadata.csv')
# AUCell obj 
data_dir2 <- '/home/gabdelha/Data/Anderson/snRNA_Seq/'
#query <- readRDS(paste0(data_dir2,"GSE214979_filtered_feature_bc_matrix_processed_seurat_neurons_only_with_Ruzicka_AUCell_M0_M9_M24_dd.rds"))
#metadata2 <- query@meta.data
metadata2 <- read.csv(paste0(data_dir2,"GSE214979_metadata_AUCell_all_scores.csv"), row.names = 1)
metadata <- left_join(metadata,metadata2[c(5,27,60:91)], by = 'X')

metadata$geneset0_status <- coalesce(metadata$geneset_pos0.y, metadata$geneset_neg0.y)
metadata$geneset24_status <- coalesce(metadata$geneset_pos24.y, metadata$geneset_neg24.y)

q1 = round(quantile(metadata$geneset_dd.x,0.25,na.rm=T),4)
q2 = round(quantile(metadata$geneset_dd.x,0.75,na.rm=T),4)
metadata = metadata %>% mutate(geneset_dd_status = case_when(geneset_dd.x <= q1 ~ "low",
                                                 geneset_dd.x >= q2 ~ "high"))

metadata$X = paste0('GSE214979_atac_fragments#',metadata$X)
metadata1 = left_join(df, metadata, by = 'X')

names = colnames(metadata1)[-1]
# loop through through columns
for (i in 1:length(names)) {
    print(names[i])
    proj <- addCellColData(ArchRProj = proj, data = paste0(metadata1[[names[i]]]),
                           cells = metadata1$X, name = names[i], force = TRUE)}

## Subset your data

idxSample <- BiocGenerics::which(proj$predicted.id.x %in% c("Excitatory","Inhibitory"))
table(proj$predicted.id.x[idxSample])
cellsSample <- proj$cellNames[idxSample]

proj <- proj[proj$id !='NA']
proj = proj[cellsSample,]

#proj = subsetArchRProject(
#  ArchRProj = proj,
#  cells = cellsSample,
#  outputDirectory = PROJDIR2,
#  force = TRUE
#)

proj
table(proj$predicted.id.x)
table(proj$predicted.id.y)

proj = saveArchRProject(ArchRProj = proj, outputDirectory= PROJDIR2)


#-----------------------------------------------------------------------------------------------------------------
# Peak Calling 
#-----------------------------------------------------------------------------------------------------------------

## Calling peaks
print('Calling Peaks')

# all nuclei -- predicted.id.x or subs
# neurons -- predicted.id.y
# neurons -- geneset24_status

proj <- addGroupCoverages(proj, groupBy = "geneset24_status", sampleLabels = 'id',
                          minCells = 40, maxCells = 2000, # per replicate
                          minReplicates = 2, maxReplicates = 12, #12 animal+region sample (celltype/sample pseudo-bulk replicates)
                          force=TRUE)

pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(proj, groupBy = "geneset24_status",
                                 pathToMacs2 = pathToMacs2,
                                 reproducibility = "(n+1)/2",
                                 genomeSize = 2.7e9,
                                 peaksPerCell = 500, maxPeaks = 300000,
                                 minCells = 25, cutOff = 0.1,
                                 promoterRegion = c(5000, 100),
                                 force=TRUE, plot = FALSE)


getPeakSet(proj)
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)

proj = saveArchRProject(ArchRProj = proj, outputDirectory= PROJDIR2)
