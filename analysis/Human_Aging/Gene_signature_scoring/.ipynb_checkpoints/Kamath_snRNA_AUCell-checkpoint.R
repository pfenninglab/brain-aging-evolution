library(Seurat)
library(patchwork)
library(dplyr)
library(AUCell)
library(GSEABase)
library(ggplot2)
library(reshape2)
library(readxl)

proj_dir = '/Aging/'
data_dir <- '/singleCell/Kamath2022_DAneuron/snRNA-seq/processed/'
data_dir2 <- '/Aging/Data/Kamath2022/'
data_dir3 <- '/Aging/Aging_Clock/Output/'


query <- readRDS(paste0(data_dir,'Kamath2022_DAneuron_snRNA_hs_processed.rds'))
# exprMatrix
exprMatrix <- query@assays$RNA$data
dim(exprMatrix)

# list of genesets 

# Lists of Genesets ---------------------------------------------------------------------------------------------------------------
# ------------- Aging Signatures -------------
geneset1 <- read.csv(paste0(proj_dir,"ENR_Age_Prediction/Output/24_10_02_model_M24_coef_gene_filtering_mean_age_corrected_pred_age.csv"), row.names = 1)
# M24
pos24 <- rownames(subset(geneset1, M_24 > 0))
neg24 <- rownames(subset(geneset1, M_24 < 0))
# MO -- subset top n number of genes -- n = number of genes from last model 
pos <- subset(geneset1, M_0 > 0)
pos0 <- rownames(pos[order(pos$M_0, decreasing = T),])[1:length(pos24)]
neg <- subset(geneset1, M_0 < 0)
neg0 <- rownames(neg[order(neg$M_0, decreasing = F),])[1:length(neg24)]

# ------------- DNA Damage Signature -------------
dna_damage <- read.csv(paste0(proj_dir,"DNA_Damage/celltype_specific_gh2ax_plus_vs_NeuNplus_gh2ax_minus_l2fc1_signature.csv"), row.names = 1)
dna_damage <- unique(dna_damage[dna_damage$pvalue <0.05,]$hg_symbol)

# ------------- Senescence Signatures -------------
sen <- read_xlsx(paste0(proj_dir,"Senescence/43587_2021_142_MOESM3_ESM.xlsx"))
can_sen <- na.omit(sen$`Canonical Senesce Pathway`[2:280])
sen_init <- na.omit(sen$`Senescence Initiators`[2:280])
sen_resp <- na.omit(sen$`Senescence Responses`[2:280])
cellage <- na.omit(sen$`CellAge`[2:280])


geneSets <- c(
  GeneSet(pos0, setName="pos0"),
  GeneSet(neg0, setName="neg0"),
  GeneSet(pos24, setName="pos24"),
  GeneSet(neg24, setName="neg24"),
  GeneSet(dna_damage, setName="dna_damage"),
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
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) # assign=TRUE


attributes(cells_AUC)
auc_m <- as.data.frame(t(cells_AUC@assays@data$AUC))
auc_m$Id <- rownames(auc_m)
query@meta.data = left_join(query@meta.data, auc_m, by = "Id")

query@meta.data$geneset_sc0 <- query@meta.data$pos0 - query@meta.data$neg0
query@meta.data$geneset_sc24 <- query@meta.data$pos24 - query@meta.data$neg24
query@meta.data = query@meta.data %>% mutate(geneset_sc0_class = case_when(geneset_sc0 > 0 ~ "pos",
                                                 geneset_sc0 < 0 ~ "neg"))
query@meta.data = query@meta.data %>% mutate(geneset_sc24_class = case_when(geneset_sc24 > 0 ~ "pos",
                                                 geneset_sc24 < 0 ~ "neg"))
rownames(query@meta.data) <- query@meta.data[,'Id']
query@meta.data


# Global k1
# canonical sensescence 
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

query@meta.data = left_join(query@meta.data, assignmentTable, by = "Id")

###################

# Global k1
# senesecence initiation 
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

query@meta.data = left_join(query@meta.data, assignmentTable, by = "Id")

###################

# Global k1
# senescence response 
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

query@meta.data = left_join(query@meta.data, assignmentTable, by = "Id")

###################

# Global k1
# cell age 
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

query@meta.data = left_join(query@meta.data, assignmentTable, by = "Id")
query@meta.data

###################
###################
###################

#saveRDS(query,paste0(data_dir2,"Kamath2022_DAneuron_snRNA_hs_processed_AUCell.rds"))