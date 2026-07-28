library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
set.seed(1)

addArchRThreads(threads = 42) 
addArchRGenome("hg38")

PROJ_DIR <- '/Aging/'
DIR <- paste0(PROJ_DIR,'Data/Anderson/Multiomics/ATAC/ArchR/Output/')
setwd(DIR)

## Loading ArchR Project + Adding Metadata
print('Loading ArchR Project + Adding Metadata')

#saveArchRProject(ArchRProj = proj, outputDirectory = "ArchRProject", load = FALSE)
proj <- loadArchRProject(path = "/Aging/Data/Anderson/Multiomics/ATAC/ArchR/ArchRProject_Neurons", force = FALSE, showLogo = FALSE)
getAvailableMatrices(proj)

## Adding coldata 
df = as.data.frame(proj$cellNames)
colnames(df) <- 'X'
# sort new column to the order of cellNames in ArchR project
metadata = read.csv(paste0(PROJ_DIR,'Data/Anderson/Multiomics/Raw/GSE214979_cell_metadata.csv'))
# AUCell obj 
data_dir2 <- paste0(PROJ_DIR,'Data/Anderson/snRNA_Seq/')
metadata2 <- read.csv(paste0(data_dir2,"GSE214979_metadata_AUCell_all_scores.csv"), row.names = 1)
metadata <- left_join(metadata,metadata2[c(5,27,60:91)], by = 'X')

metadata$geneset0_status <- coalesce(metadata$geneset_pos0.y, metadata$geneset_neg0.y)
metadata$geneset24_status <- coalesce(metadata$geneset_pos24.y, metadata$geneset_neg24.y)

metadata$X = paste0('GSE214979_atac_fragments#',metadata$X)
metadata1 = left_join(df, metadata, by = 'X')

names = colnames(metadata1)[-1]
# loop through through columns
for (i in 1:length(names)) {
    print(names[i])
    proj <- addCellColData(ArchRProj = proj, data = metadata1[[names[i]]],
                           cells = metadata1$X, name = names[i], force = TRUE)}
#getCellColData(proj)
#--------------------------------------------------------
#--------------------------------------------------------

# Subset Neurons only 

# Subset your data by the conditions you want to compare
idxSample1 <- BiocGenerics::which(proj$predicted.id.x %in% c("Excitatory"))
table(proj$predicted.id.x[idxSample1])
cellsSample1 <- proj$cellNames[idxSample1]

idxSample2 <- BiocGenerics::which(proj$Status %in% c("Ctrl"))
table(proj$Status[idxSample2])
cellsSample2 <- proj$cellNames[idxSample2]

idxSample <- intersect(idxSample1, idxSample2)
cellsSample <- intersect(cellsSample1, cellsSample2)

proj = proj[cellsSample,]
proj
table(proj$Status, proj$predicted.id.x)

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
#write.table(regions, paste0(DIR,'Anderson_PeakMatrix_Reads.bed'), 
#            row.names=FALSE,quote =FALSE, col.names=FALSE,sep = "\t") 


# Convert to PeakMatrix to data frame
# Add region identifier
Peak_df <- as.data.frame(as.matrix(assay(PeakMatrix)))   # columns = cells
print(dim(Peak_df))
rm(PeakMatrix)

combined_df <- cbind(regions, Peak_df)
combined_df$read <- paste0(combined_df$seqnames, ":", combined_df$start, "-", combined_df$end)
#write.table(combined_df, paste0('Anderson_PeakMatrix.tsv'), row.names=FALSE, sep = "\t")

# Stop here 
#----------------------------------------------------------

# Step 1b
# Intersect with annotations -- bedtools in bash 
#bedtools_genome_anno_intersect.sh

#----------------------------------------------------------

# Start again here 

# Step 2
print('Step 2')

# Merge annotations with reads 
# Annotate regions 
#regions <- read.table(paste0('Anderson_PeakMatrix_Reads.bed'))
colnames(regions) <- c('seqnames','start','end')
regions$read <- paste0(regions$seqnames,':',regions$start,'-',regions$end)
regions <- regions[,c('read'), drop=FALSE]

anno <- read.table(paste0('Anderson_PeakMatrix_Reads_rm_epigenomics_anno_intersect.txt'))
anno <- anno[,c(1:3,7)]
colnames(anno) <- c('seqnames','start','end','Description')
anno$read <- paste0(anno$seqnames,':',anno$start,'-',anno$end)
anno <- anno[,c('Description', 'read')]

regions <- left_join(regions, anno, by='read')
#write.table(regions, paste0('Anderson_PeakMatrix_Reads_Annotated.bed'), row.names=FALSE, sep = "\t")
#----------------------------------------------------------

# Step 3
print('Step 3')

#combined_df <- read.table(paste0(DIR,'Anderson_PeakMatrix.tsv'), header = T)
combined_df <- left_join(combined_df, regions, by = "read") %>% na.omit()
print(dim(combined_df))

#combined_df[,4:ncol(combined_df)-2] <- ((combined_df[,4:ncol(combined_df)-2])+1)/2

df <- combined_df[,4:ncol(combined_df)] %>%
  group_by(Description) %>%
  summarise(across(4:ncol(.)-2, sum)) 

df[2:ncol(df)] <- sweep(df[2:ncol(df)], 2, colSums(df[2:ncol(df)]), FUN = "/")
rownames(df) <- df$Description 
#df <- df[,-which(names(df) %in% 'Description')]
            
write.csv(df, paste0(DIR,'Anderson_Exc_Ctrl_PeakMatrix_reads_RoadMap_18States.csv'), row.names=F)

#----------------------------------------------------------

# Step 4
print('Step 4')

## Adding coldata 
df <- as.data.frame(t(df))
df <- df[-1,]
df$X <- rownames(df)
df = left_join(df,metadata1, by = 'X')
rownames(df) <- df$X

write.csv(tbl, paste0(DIR,'Anderson_Exc_Ctrl_PeakMatrix_reads_RoadMap_18States_plus_metadata.csv'))
