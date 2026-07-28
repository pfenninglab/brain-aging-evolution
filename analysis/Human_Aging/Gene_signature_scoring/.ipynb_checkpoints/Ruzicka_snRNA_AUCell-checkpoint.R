library(Seurat)
library(SeuratDisk)
library(AUCell)
library(readxl)
library(patchwork)
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(GSEABase)
library(broom)
library(ggrepel)
library(gridExtra) # Or use the `patchwork` package for flexibility
library(ComplexHeatmap)
library(circlize)

proj_dir = '/Aging/'

### directories
data_dir = '/singleCell/Ruzicka_snRNA_Seq/data/'
data_dir2 = '/Aging/Data/Ruzicka_snRNA_Seq/'
out_dir2 = '/Aging/Aging_Clock/Output/'

PLOTDIR = '/Aging/Aging_Clock/Plots/'
setwd(PLOTDIR)

obj <- readRDS(paste0(data_dir,"22_08_12_seurat_raw_FILTERED_neurons_only_controls_only_lognormalized.rds"))
# exprMatrix
#obj = obj[VariableFeatures(obj),]
exprMatrix <- obj@assays[["RNA"]]@data
dim(exprMatrix)

# list of genesets 

# Lists of Genesets ---------------------------------------------------------------------------------------------------------------
# ------------- Aging Signatures -------------
geneset1 <- read.csv(paste0(proj_dir,"ENR_Age_Prediction/Output/24_10_02_model_M24_coef_gene_filtering_mean_age_corrected_pred_age.csv"), row.names = 1)
# M24
geneset_pos24 <- rownames(subset(geneset1, M_24 > 0))
geneset_neg24 <- rownames(subset(geneset1, M_24 < 0))
# MO -- subset top n number of genes -- n = number of genes from last model 
geneset_pos <- subset(geneset1, M_0 > 0)
geneset_pos0 <- rownames(geneset_pos[order(geneset_pos$M_0, decreasing = T),])[1:length(geneset_pos24)]
geneset_neg <- subset(geneset1, M_0 < 0)
geneset_neg0 <- rownames(geneset_neg[order(geneset_neg$M_0, decreasing = F),])[1:length(geneset_neg24)]

# ------------- DNA Damage Signature -------------
geneset_dd <- read.csv("/Aging/DNA_Damage/celltype_specific_gh2ax_plus_vs_NeuNplus_gh2ax_minus_l2fc1_signature.csv", row.names = 1)
geneset_dd <- unique(geneset_dd[geneset_dd$pvalue <0.05,]$hg_symbol)

# ------------- Senescence Signatures -------------
geneset_sen <- read_xlsx("/Aging/Senescence/43587_2021_142_MOESM3_ESM.xlsx")
can_sen <- na.omit(geneset_sen$`Canonical Senesce Pathway`[2:280])
sen_init <- na.omit(geneset_sen$`Senescence Initiators`[2:280])
sen_resp <- na.omit(geneset_sen$`Senescence Responses`[2:280])
cellage <- na.omit(geneset_sen$`CellAge`[2:280])


geneSets <- c(
  GeneSet(pos0, setName="geneset_pos0"),
  GeneSet(neg0, setName="geneset_neg0"),
  GeneSet(pos24, setName="geneset_pos24"),
  GeneSet(neg24, setName="geneset_neg24"),
  GeneSet(geneset_dd, setName="geneset_dd"),
  GeneSet(can_sen, setName="can_sen"),
  GeneSet(sen_init, setName="sen_init"),
  GeneSet(sen_resp, setName="sen_resp"),
  GeneSet(cellage, setName="cellage"))

geneSets
geneSets <- GeneSetCollection(geneSets)
geneSets

# ---------------------------------------------------------------------------------------------------------------------------------

## 1. Score gene signatures
# 1.1. Build gene-expression rankings for each cell

#plotGeneCount(exprMatrix)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
# 1.2. Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)

# ---------------------------------------------------------------------------------------------------------------------------------

# Add Scores to Metadata
attributes(cells_AUC)
auc_m <- as.data.frame(t(cells_AUC@assays@data$AUC))
auc_m$Id <- rownames(auc_m)
obj@meta.data = left_join(obj@meta.data, auc_m, by = "Id")

obj@meta.data$geneset_sc0 <- obj@meta.data$geneset_pos0 - obj@meta.data$geneset_neg0
obj@meta.data$geneset_sc24 <- obj@meta.data$geneset_pos24 - obj@meta.data$geneset_neg24

obj@meta.data = obj@meta.data %>% mutate(geneset_sc0_class = case_when(geneset_sc0 > 0 ~ "pos",
                                                 geneset_sc0 < 0 ~ "neg"))
obj@meta.data = obj@meta.data %>% mutate(geneset_sc24_class = case_when(geneset_sc24 > 0 ~ "pos",
                                                 geneset_sc24 < 0 ~ "neg"))
obj@meta.data = obj@meta.data %>% mutate(age_status = case_when(Age <= 52 ~ "young",
                                                 Age >= 69 ~ "old"))
#table(obj@meta.data$geneset_sc0_class, obj@meta.data$age_status)
#table(obj@meta.data$geneset_sc24_class, obj@meta.data$age_status)


# Aging Signature -------------------------------------------------------------------------
# Predicted Age Method --------------------------------------------------------------------

# M0
# geneset_pos
# Global k1
Global_k1 = cells_assignment$geneset_pos0$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("^geneset_pos0$", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$geneset_pos0 <- 'geneset_pos0'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
#------------------------------------------------------------------------------------------

# geneset_neg
# Global k1
Global_k1 = cells_assignment$geneset_neg0$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("^geneset_neg0$", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$geneset_neg0 <- 'geneset_neg0'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
#------------------------------------------------------------------------------------------

# M24
# geneset_pos
# Global k1
Global_k1 = cells_assignment$geneset_pos24$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("^geneset_pos24$", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$geneset_pos24 <- 'geneset_pos24'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
#------------------------------------------------------------------------------------------

# geneset_neg
# Global k1
Global_k1 = cells_assignment$geneset_neg24$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("^geneset_neg24$", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$geneset_neg24 <- 'geneset_neg24'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
#------------------------------------------------------------------------------------------

# DNA Damage Signature --------------------------------------------------------------------

# Global k1
Global_k1 = cells_assignment$geneset_dd$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("geneset_dd", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$geneset_dd <- 'geneset_dd'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
#------------------------------------------------------------------------------------------

# Senescenece Signature -------------------------------------------------------------------

# Canonical Senescence 
# Global k1
Global_k1 = cells_assignment$can_sen$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("^can_sen$", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$can_sen <- 'can_sen'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
#------------------------------------------------------------------------------------------

# Senescence Initiaiton 
# Global k1
Global_k1 = cells_assignment$sen_init$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("^sen_init$", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$sen_init <- 'sen_init'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
#------------------------------------------------------------------------------------------

# Senescence Response
# Global k1
Global_k1 = cells_assignment$sen_resp$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("^sen_resp$", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$sen_resp <- 'sen_resp'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
#------------------------------------------------------------------------------------------

# Cell Age
# Global k1
Global_k1 = cells_assignment$cellage$aucThr$thresholds[,'threshold'][['Global_k1']]
geneSetName <- rownames(cells_AUC)[grep("^cellage$", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=Global_k1)
abline(v=Global_k1)
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>Global_k1))
length(newSelectedCells)
assignmentTable <- as.data.frame(newSelectedCells)
assignmentTable$cellage <- 'cellage'
colnames(assignmentTable)[1] <- c('Id')
nrow(assignmentTable)
obj@meta.data = left_join(obj@meta.data, assignmentTable, by = "Id")
obj@meta.data
#------------------------------------------------------------------------------------------

#write.csv(obj@meta.data,paste0(data_dir2,"Ruzicka_metadata_AUCell_all_scores.csv"))