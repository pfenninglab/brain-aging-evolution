###############################################################################
# HPC Job Configuration
###############################################################################
# Request HPC resources: 1 node, 64 tasks, 8 hours runtime
# srun -p RM-shared -N 1 --ntasks-per-node=64 -t 8:00:00 --pty /bin/bash

###############################################################################
# Load Required Libraries
###############################################################################
library(ArchR)          # Framework for single-cell ATAC-seq analysis
library(tidyverse)      # Collection of data manipulation packages
library(here)           # Project-relative file path management
library(limma)          # Differential expression analysis
library(edgeR)          # Normalization of count data
library(swfdr)          # Multiple testing correction
library(rtracklayer)    # For reading and writing genomic data

###############################################################################
# Initial Data Processing and Cell Classification
###############################################################################
# Configure ArchR settings for parallel processing
#n_cores = Sys.getenv("SLURM_NTASKS_PER_NODE") %>% as.numeric()
n_cores = parallel::detectCores() %>% as.numeric()
addArchRThreads(threads = n_cores / 4 %>% round())

# Set up output directories
DATADIR='/Aging/Data/Anderson/Multiomics/ATAC/ArchR/'
PLOTDIR='/Aging/Aging_Clock/Plots/'
setwd(DATADIR)

# Define pre-seq variables
DIST_TO_TSS <- 10 * 10^3 # how far from TSS to annotate peak as "promoter"

dir.create(here(DATADIR, 'peaks'), showWarnings = F, recursive = T)
dir.create(here(DATADIR, 'rdas'), showWarnings = F, recursive = T)
dir.create(here(DATADIR, 'tables'), showWarnings = F, recursive = T)

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# String splitting utility
ss <- function(x, pattern, slot = 1, ...) {
  sapply(strsplit(x = x, split = pattern, ...), "[", slot)
}

# for narrowPeak r function utilities
#source(here('code/final_code/hal_scripts/narrowPeakFunctions.R'))
source('/singleCell/Stauffer_macaque_multiome_VentralPutamen/code/ldsc_multiomeATAC_Striatum/narrowPeakFunctions.R')

# this is a published R function that we will use for limma-voom DEA
source("https://github.com/YOU-k/voomByGroup/raw/main/voomByGroup.R")


###############################################################################
# Load project and add metadata
###############################################################################

proj <- loadArchRProject(path = "/Aging/Data/Anderson/Multiomics/ATAC/ArchR/ArchRProject_Neurons", force = FALSE, showLogo = FALSE)

## Adding coldata 
df = as.data.frame(proj$cellNames)
colnames(df) <- 'X'
# sort new column to the order of cellNames in ArchR project
metadata = read.csv('/Aging/Data/Anderson/Multiomics/Raw/GSE214979_cell_metadata.csv')
# AUCell obj 
data_dir2 <- '/Aging/Data/Anderson/snRNA_Seq/'
#query <- readRDS(paste0(data_dir2,"GSE214979_filtered_feature_bc_matrix_processed_seurat_neurons_only_with_Ruzicka_AUCell_M0_M9_M24_dd.rds"))
#metadata2 <- query@meta.data
metadata2 <- read.csv(paste0(data_dir2,"GSE214979_metadata_AUCell_all_scores.csv"), row.names = 1)
metadata <- left_join(metadata,metadata2[c(5,27,60:91)], by = 'X')

# Add geneset status columns 

# geneset_neg_status
metadata <- metadata %>% rowwise() %>% 
mutate(geneset_neg_status = {
      x <- c(geneset_neg24.x, geneset_neg0.x)
      nm <- c("geneset_neg24", "geneset_neg0")
      if (all(is.na(x))) NA_character_
      else if (sum(x == max(x, na.rm = TRUE), na.rm = TRUE) > 1) "tie"
      else nm[which.max(x)]
    }) %>% ungroup()

# geneset_pos_status
metadata <- metadata %>% rowwise() %>%
mutate(geneset_pos_status = {
      x <- c(geneset_pos24.x, geneset_pos0.x)
      nm <- c("geneset_pos24", "geneset_pos0")
      if (all(is.na(x))) NA_character_
      else if (sum(x == max(x, na.rm = TRUE), na.rm = TRUE) > 1) "tie"
      else nm[which.max(x)]
    }) %>% ungroup()

# geneset_intrinsic
metadata <- metadata %>% rowwise() %>% 
mutate(geneset_intrinsic = {
      x <- c(geneset_neg24.x, geneset_pos24.x)
      nm <- c("geneset_neg24", "geneset_pos24")
      if (all(is.na(x))) NA_character_
      else if (sum(x == max(x, na.rm = TRUE), na.rm = TRUE) > 1) "tie"
      else nm[which.max(x)]
    }) %>% ungroup()

# geneset_systemic
metadata <- metadata %>% rowwise() %>%
mutate(geneset_pos_status = {
      x <- c(geneset_neg0.x, geneset_pos0.x)
      nm <- c("geneset_neg0", "geneset_pos0")
      if (all(is.na(x))) NA_character_
      else if (sum(x == max(x, na.rm = TRUE), na.rm = TRUE) > 1) "tie"
      else nm[which.max(x)]
    }) %>% ungroup()


# geneset_dd_status
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

proj <- proj[proj$id !='NA']
proj <- proj[proj$predicted.id.y !='Ex-NRGN']


###############################################################################
# Subset your data by the conditions you want to compare -- FOR WITHIN COMPARISONS
###############################################################################

idxSample1 <- BiocGenerics::which(proj$predicted.id.x %in% c("Inhibitory"))
cellsSample1 <- proj$cellNames[idxSample1]
idxSample2 <- BiocGenerics::which(proj$Status %in% c("AD"))
cellsSample2 <- proj$cellNames[idxSample2]
idxSample <- intersect(idxSample1, idxSample2)
cellsSample <- intersect(cellsSample1, cellsSample2)
proj = proj[cellsSample,]

table(proj$Status, proj$geneset_neg0.y)


# one donor × intrinsic/systemic state
metadata=getCellColData(proj)
proj$Sample_Geneset <- paste(
  proj$id,
  proj$geneset_neg0.y,   # intrinsic/systemic
  sep = "#"
)


###############################################################################
# Process ATAC-seq Data and Perform Differential Analysis
###############################################################################

# 1. Get Peaks & Annotate
peak_gr = getPeakSet(proj)
peak_gr = nameNarrowPeakRanges(gr = peak_gr, genome = getGenome(proj))
names(peak_gr) = mcols(peak_gr)$name
peak_gr$name = names(peak_gr)
# Convert peak info to data frame
peak_df = data.frame(peak_gr) %>% mutate(name = names(peak_gr))

# 2. Build Pseudobulk (getGroupSE)
# Get pseudobulk accessibility matrices
se_peak = 
      # pseudo-bulk profiles by cell class
      getGroupSE(proj, useMatrix = 'PeakMatrix', 
                 groupBy = "Sample_Geneset", divideN = F)

# 3. Add Peak Metadata to se_peak
rowRanges(se_peak) = GRanges(rowData(se_peak))
rowRanges(se_peak) = nameNarrowPeakRanges(gr = rowRanges(se_peak), genome = getGenome(proj))
names(rowRanges(se_peak)) = mcols(rowRanges(se_peak))$name
rowData(se_peak) = peak_df[match(names(rowRanges(se_peak)), peak_df$name),]

# 4. Add Sample Metadata
se_peak$nCells = as.numeric(se_peak$nCells)
se_peak$Sample_Target = colnames(se_peak)
se_peak$Sample = ss(se_peak$Sample_Target, '#', 1)
se_peak$Target = ss(se_peak$Sample_Target, '#', 2) %>% make.names()
ind_match = match(se_peak$Sample, metadata$id)
se_peak$Sex = metadata$Sex[ind_match]
se_peak$Status = metadata$Status[ind_match]
se_peak$Batch = metadata$Sub.batch[ind_match]
se_peak$Age = as.numeric(metadata$Age[ind_match])
colData(se_peak)

# Save preprocessed se_peak
#se_peak_file = paste0('Anderson_PFC_se_peak_preprocessed.rds') %>%
#    here(DATADIR, 'rdas', .)
#saveRDS(se_peak, se_peak_file)
#se_peak = readRDS(se_peak_file)


###############################################################################
# Filter Peaks and Samples
###############################################################################
se_peak = se_peak[
    which(
      # Filter peaks (rows) on non-promoter and non-exonic peaks
      rowData(se_peak)$peakType %in% c('Distal', 'Intronic') & 
        # Filter peaks (rows) on too close to the TSS
        abs(rowData(se_peak)$distToTSS) > DIST_TO_TSS),
    # Filter samples (columns)
    #which(se_peak$nCells > 15)
  ]

# Remove NA
#keep <- !grepl("#NA$", colnames(se_peak))
#se_peak <- se_peak[, keep]        

# write out the peaks to predict on or use for HALPER
pred_peaks = rowRanges(se_peak)
pred_fn = here(DATADIR, 'peaks', 
             paste0('Anderson_PFC_all_peaks.narrowPeak.gz'))
if(!file.exists(pred_fn)){
pred_peaks = write_GRangesToNarrowPeak(
  gr = pred_peaks, file = pred_fn, genome = getGenome(proj))
}

bed_fn = here(DATADIR, 'peaks', 
             paste0('Anderson_PFC_all_peaks.bed'))
bed = as.data.frame(pred_peaks, row.names=NULL)[,c('seqnames','start','end')]
write.table(bed, bed_fn,quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

print("1. Normalize data using voom-limma")
#------------------------------------------------
y = DGEList(counts = assays(se_peak)[[1]]) %>% .[rowMeans(.$counts) > 1, , keep.lib.size = FALSE] %>% calcNormFactors()
dim(y)

print("2. Create design matrix with covariates")
#------------------------------------------------      
design <- model.matrix( ~ 0 + Target 
                       + nCells + ReadsInPeaks + FRIP,  # + Sex + Batch 
                       data = colData(se_peak))
colnames(design) <- make.names(colnames(design))  # replace invalid characters like ":" and "-"
dim(design)

print("3. Voom + duplicateCorrelation") 
#------------------------------------------------  
# (accounts for multiple pseudobulks per donor)
## voom precision weights and sample-wise quality weights normalization
print("Running voomWithQualityWeights...")
v <- voomWithQualityWeights(y, design)
print("Running duplicateCorrelation...")
cor <- duplicateCorrelation(v, design, block = se_peak$Sample) #what does this do - account for interrelatedness between celltype pseudobulks from the same sample
print(paste0("cor$consensus: ", cor$consensus)) # 0.272425023162038
## recalculate weights after adjusting for correlated samples from same subject
print("Rerunning voomWithQualityWeights including corr consensus...")
v <- voomWithQualityWeights(y, design, block = se_peak$Sample, 
          correlation = cor$consensus)
print("Rerunning duplicateCorrelation...")
cor <- duplicateCorrelation(v, design, block = se_peak$Sample)
  print(paste0("cor$consensus: ", cor$consensus)) # 0.2735817

print("4. Fit linear model")
#------------------------------------------------ 
vfit <- lmFit(v, design, block = se_peak$Sample, correlation = cor$consensus) %>% eBayes(., robust = TRUE)

print("5. Contrasts")
#------------------------------------------------ 
geneset_contrasts <- tribble(
  ~target,         ~off_target_class,   ~model,
  "geneset_neg0", "NA.",     "neg0_vs_all"
)

unique_geneset <- c(
  geneset_contrasts$target,
  unlist(strsplit(geneset_contrasts$off_target_class, "\\+"))
) %>% unique()

contrast.matrix = eval(parse(text = paste0(
"makeContrasts(",
paste0('Target', geneset_contrasts$target, "-", 'Target', 
       geneset_contrasts$off_target_class, collapse = ", "),
", levels = design)"
)))
colnames(contrast.matrix) = geneset_contrasts$model

print("6. Fit contrasts")
#------------------------------------------------ 
fit2 = contrasts.fit(vfit, contrast.matrix) %>% eBayes()

df_res = set_names(geneset_contrasts$model) %>% 
lapply(function(x) {
  topTable(fit2, coef=x, number=Inf, sort.by="P") %>%
  rownames_to_column('peak')
}) %>%
bind_rows(.id = 'model') %>%
# correct for multiple testing across all peaks by model
mutate(adj.P.Val = lm_qvalue(P.Value, X=AveExpr)$q)

# Save results
diff_results_file = 
paste0('Anderson_PFC_AD_Inh_neg0_diffPeakList_byModel.rds') %>% 
here(DATADIR, 'rdas', .)
saveRDS(split(df_res, f = df_res$model), diff_results_file)


# Save up & down-reg significant peaks
diff_gr_list = df_res %>% 
# Filter peaks with adjusted p-value < 0.05 and large effect size
filter(P.Value < 0.01, abs(logFC) > 1) %>% 
mutate(tmp = ifelse(logFC > 0, 
                    paste0(model, '.diff_up'),
                    paste0(model, '.diff_dn'))) %>%
split(f = .$tmp) %>% 
lapply(distinct) %>% 
lapply(function(x) peak_gr[x$peak])
lengths(diff_gr_list)

for(name in names(diff_gr_list)){
    diff_bed <- tibble(raw = diff_gr_list[[name]]$name) %>%
    separate(raw, into = c('prefix','chr','coord','width'), sep=':') %>%
    separate(coord, into = c('start','end'), sep='-') %>%
    mutate(
        start = as.integer(start),
        end = as.integer(end),
        width = as.integer(width)
    ) %>%
    select(chr, start, end)
    
    write.table(diff_bed, quote=F, row.names=F, col.names=F, 
            here(DATADIR, 'peaks', paste0('Anderson_PFC_AD_Inh_peaks_', name, '.bed')))
}
