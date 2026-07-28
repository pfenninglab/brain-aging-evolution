# Tutorial https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")

library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggrepel)

# SET THE DESIRED ORGANISM HERE
# https://bioconductor.org/packages/release/BiocViews.html#___OrgDb
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
set.seed(123)

out_dir = '/Aging/Aging_Clock/Output/'
PLOTDIR <- '/Aging/Aging_Clock/Plots/'
setwd(PLOTDIR)

# Manuscript Figure 3 ---------------------------------------------------------
#------------------------------------------------------------------------------

# OLS -------------------------------------------------------------------------
# Sort based on Age coefficient 

# list of genesets 
geneset <- read.csv(paste0(out_dir,"25_02_09_OLS_Ruzicka_age_fullset.csv"), row.names = 1)


pdf('Dotplot_GSEA_ClusterProfiler_age_fullset.pdf', width=8,height=7)

geneset_sub <- subset(geneset, Age_pval < 0.05 & Age_coef != 0)
original_gene_list <- geneset_sub$Age_coef
names(original_gene_list) <- rownames(geneset_sub)
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = T)
length(gene_list)

# Check gene annotations for the chosen organism 
keytypes(org.Hs.eg.db)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             #nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH",
            seed = T)

# DotPlot
require(DOSE)
dotplot(gse, showCategory=10, split=".sign", x = 'enrichmentScore') + 
facet_grid(.~.sign) + 
scale_fill_gradientn(colours = c("#21FE80","#000435"), limits = c(min(gse$p.adjust),max(gse$p.adjust)))

emapplot(pairwise_termsim(gse), showCategory = 10) + 
scale_fill_gradientn(colours = c("#21FE80","#000435"), limits = c(min(gse$p.adjust),max(gse$p.adjust)))

dev.off()
######################


# OLS -------------------------------------------------------------------------
# Sort genes based on difference between predicted age 24 coef and Age coef

# list of genesets 
geneset1 <- read.csv(paste0(out_dir,"25_02_09_OLS_Ruzicka_age_fullset.csv"), row.names = 1)
#geneset1 <- subset(geneset1, Age_coef < 0)
geneset1 <- geneset1[,'Age_coef', drop=F]
geneset1$Gene <- rownames(geneset1)

geneset2 <- read.csv(paste0(out_dir,"25_02_09_OLS_Ruzicka_pred_age_24_fullset.csv"), row.names = 1)
geneset2 <- subset(geneset2, predicted_age_24_pval < 0.05 & predicted_age_24_coef != 0)
geneset2 <- geneset2[,'predicted_age_24_coef', drop = F]
geneset2$Gene <- rownames(geneset2)

geneset3 <- read.csv(paste0(out_dir,"25_02_09_OLS_Ruzicka_pred_age_0_fullset.csv"), row.names = 1)
#geneset3 <- subset(geneset3, predicted_age_0_pval < 0.05 & predicted_age_0_coef != 0)
geneset3 <- geneset3[,'predicted_age_0_coef', drop = F]
geneset3$Gene <- rownames(geneset3)

geneset <- left_join(geneset2, geneset1, by='Gene')
geneset$diff <- geneset$predicted_age_24_coef - geneset$Age_coef
geneset <- geneset[order(geneset$diff, decreasing = T),]
rownames(geneset) <- geneset$Gene

######################

pdf('Dotplot_GSEA_ClusterProfiler_sign_non-zero_predicted_age_24_age_diff_fullset.pdf', width=8,height=7)

geneset_sub <- geneset
original_gene_list <- geneset_sub$diff # Sort based on diff variable
names(original_gene_list) <- rownames(geneset_sub)
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = T)
length(gene_list)

# Check gene annotations for the chosen organism 
keytypes(org.Hs.eg.db)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             #nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH",
            seed = T)

# DotPlot
require(DOSE)
dotplot(gse, showCategory=10, split=".sign", x = 'enrichmentScore') + 
facet_grid(.~.sign) + 
scale_fill_gradientn(colours = c("#21FE80","#000435"), limits = c(min(gse$p.adjust),max(gse$p.adjust))) + 
scale_size_continuous(limits = c(5, 200))

emapplot(pairwise_termsim(gse), showCategory = 10) + 
scale_fill_gradientn(colours = c("#21FE80","#000435"), limits = c(min(gse$p.adjust),max(gse$p.adjust))) + 
scale_size_continuous(limits = c(5, 200))

dev.off()


#--------------------------------------------------------------------

# Plotting 

pdf('Scatterplot_predicted_age_24_age_coef_highlight_genes.pdf', width=5,height=5)

highlight_genes <- c('MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L', 'MT-ND5', 'MT-CO1', 'MT-CO2', 'MT-CYB', 'MT-ATP6',
                     "COX5B", "NDUFV1", 'CHCHD10','UQCRFS1')
ggplot(geneset, aes(x = Age_coef , y = diff)) +
  geom_point(fill = "#21FE80", color = "#000435", shape = 21, size = 3, stroke = 1) +
  geom_smooth(method = "lm", color = "#000435", fill = "#21FE80", linetype = "dashed", size = 1, se = TRUE) +  # Add regression line
  geom_hline(yintercept = 0, color = "#000435", linetype = "dashed") +  # Add x-axis at y=0
  geom_vline(xintercept = 0, color = "#000435", linetype = "dashed") +
  geom_text_repel(data = subset(geneset, rownames(geneset) %in% highlight_genes), # Use `data` here
                  aes(label = rownames(geneset)[rownames(geneset) %in% highlight_genes]),
                  size = 4, 
                  fontface = "bold",
                  color = "#000435", 
                  nudge_x = 0.2) +
  labs(x = "Age Coefficient",
       y = "Intrinsic Age Coefficient - Age Coefficient") + 
  #coord_cartesian(xlim = c(-0.05, 0.05), ylim = c(-0.12, 0.12)) + 
  theme_minimal()
dev.off()
######################