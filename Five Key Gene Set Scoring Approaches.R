library(AUCell)
library(GSVA)
library(clusterProfiler)
library(tidyverse)
library(xlsx)
library(Seurat)
library(ggpubr)
library(readxl)
library(reshape2)

# AUCell ------------------------------------------------------------------

obj <- readRDS("/home/data/zhuang/CT/3.1.PCa/PCa_merged_with_Ctrl_annotated.Rds") 

# 批量读取gmt制作list
all_data <- list()
gmt_files <- list.files("2.1.1PCa_by_c/打分/gmt/growth factor", pattern = "\\.gmt$", full.names = TRUE)
for(i in seq_along(gmt_files)) {
  gmt_df <- read.gmt(gmt_files[i])
  term_name <- as.character(unique(gmt_df$term))
  genes <- gmt_df$gene
  all_data[[term_name]] <- genes
}
max_length <- max(sapply(all_data, length))
all_data_padded <- lapply(all_data, function(x) {
  length(x) <- max_length
  return(x)
})
gene_data <- as_tibble(all_data_padded)
write.xlsx(gene_data, "gene_data.xlsx")

genesets <- lapply(gene_data, function(col) {
  col[!is.na(col) & col != ""] 
})

expMtx <- GetAssayData(obj, assay = "RNA", layer = "data")

cells_rankings <- AUCell_buildRankings(expMtx,
                                       splitByBlocks = TRUE, # 将细胞分组处理，减少内存使用
                                       plotStats = TRUE # 显示每个细胞中基因表达的分布统计
) 
cells_rankings

cells_AUC <- AUCell_calcAUC(genesets, 
                            cells_rankings, 
                            aucMaxRank = nrow(cells_rankings)*0.1 # 只考虑排名在前 aucMaxRank 位的基因
                            )
cells_AUC

## 小提琴图视化 -----

geneSet <- "GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
AUCell_auc
obj$AUCell <- AUCell_auc
head(obj@meta.data)

obj_meta <- obj@meta.data %>%
  dplyr::select(celltype, AUCell) %>%
  mutate(celltype = as.factor(celltype))

p <- ggplot(obj_meta, aes(x = celltype, y = AUCell, fill = celltype)) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(alpha = 0.1, color = "black", size = 0.05) +
  labs(title = paste0("GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION Score by Cell Type"),
       x = "Cell Type",
       y = "Module Score") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  stat_compare_means(
    comparisons = combn(levels(obj_meta$celltype), 2, simplify = FALSE),
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01,
    vjust = 0.5
  );p

## on/off -----

### 自动计算阈值 -----

cells_assignment <- AUCell_exploreThresholds(cells_AUC,     
                                             plotHist = TRUE,
                                             assignCells = TRUE # 是否根据AUC值自动分配细胞到不同的基因集
)
names(cells_assignment)
cells_assignment$GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION$aucThr$thresholds    
threshold    nCells        
tenPercentOfMax 0.07777778   1783    适合需要最大化检测率的场景 
Global_k1       0.36459154    942    适合需要高置信度的严格分析
minimumDens     0.17960426   1604    适合平衡准确性和完整性的常规分析

cells_assignment$GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION$aucThr$selected
length(cells_assignment$GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION$assignment)
head(cells_assignment$GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION$assignment)

### 手动设置阈值 -----

geneSetName <- rownames(cells_AUC)[grep("GOBP_FIBROBLAST_GROWTH_FACTOR_PRODUCTION", rownames(cells_AUC))]    
AUCell_plotHist(cells_AUC[geneSetName, ], aucThr = 0.25)
abline(v = 0.25)

### 可视化 -----

newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName, ] > 0.25))  
length(newSelectedCells)
head(newSelectedCells) 

par(mfrow = c(2,3))        
selectedThresholds <- getThresholdSelected(cells_assignment)    
nBreaks <- 5        
color_Neg <- grDevices::colorRampPalette(c("#1B3E6C", "#009DFB", "#75AADB"))(nBreaks)  
color_Pos <- grDevices::colorRampPalette(c("#FFC0CB", "#FF00FF", "#FF0000"))(nBreaks)            

cellsUmap <- obj@reductions$umap@cell.embeddings

for (geneSetName in names(selectedThresholds)) {     
  passThreshold <- getAUC(cells_AUC)[geneSetName, ] > selectedThresholds[geneSetName]  
  if (sum(passThreshold) > 0) {               
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)       
    cellColor <- c(setNames(color_Neg[cut(aucSplit[[1]], breaks = nBreaks)], names(aucSplit[[1]])),    
                   setNames(color_Pos[cut(aucSplit[[2]], breaks = nBreaks)], names(aucSplit[[2]])))    
    plot(cellsUmap, main=geneSetName,
         sub = "Pink/Red cells pass the threshold",
         col =cellColor[rownames(cellsUmap)], pch = 16)       
  }           
}    

selectedThresholds[1]
selectedThresholds[2]
par(mfrow = c(2,3))    
AUCell_plotTSNE(tSNE = cellsUmap, exprMat = GetAssayData(obj, layer = "data"),    
                cellsAUC = cells_AUC[1:2, ], 
                thresholds = selectedThresholds)

newThresholds <- selectedThresholds      
newThresholds[1] <- 0.1      
newThresholds[2] <- 0.25       
par(mfrow = c(2,3))       
AUCell_plotTSNE(tSNE = cellsUmap, exprMat = GetAssayData(obj, layer = "data"),       
                cellsAUC = cells_AUC[1:2, ], thresholds = newThresholds)

# GSVA ------------------------------------------------------------------

## 每个细胞打分 -----

obj <- readRDS("/home/data/zhuang/CT/3.1.PCa/PCa_merged_with_Ctrl_annotated.Rds") 
gene_data <- read_excel("2.1.1PCa_by_c/打分/gmt/heparin.xlsx", col_names = TRUE, sheet = "Sheet 1")
head(gene_data)
genesets <- lapply(gene_data, function(col) {
  col[!is.na(col) & col != ""] 
})
names(genesets)

expMtx <- GetAssayData(obj, assay = "RNA", slot = "data")

detectCores()
time <- system.time({
  GSVA <- gsva(expMtx, genesets, method = "gsva",
               kcdf = "Gaussian",
               verbose = TRUE, 
               parallel.sz = 20)
})
saveRDS(GSVA, "2.1.1PCa_by_c/打分/HEPARIN_GSVA.Rds")

GSVA <- readRDS("2.1.1PCa_by_c/打分/HEPARIN_GSVA.Rds") 
celltype <- obj$celltype %>% as.data.frame()
colnames(celltype) <- "celltype"
table(obj$celltype)
GSVA <- GSVA %>% as.data.frame() %>% t()
mtx <- cbind(GSVA, celltype)
mtx$celltype <- factor(mtx$celltype, levels = c("General", "Naive_Dom", "NCHT_General", "NCHT_nR_Dom", "NCHT_R_Dom"))
mtx <- mtx[order(mtx$celltype),]

data <- melt(mtx)

P = ggplot(data, aes(x = variable, y = value, fill = celltype)) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, colour="black", size = 12),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("HEPARIN pathways") +
  ylab("GSVA score") +
  xlab("");P

comparisons <- list(c("NCHT_nR_Dom", "NCHT_R_Dom"))
comparisons <- combn(levels(data$celltype), 2, simplify = FALSE)

ggplot(data, aes(x = celltype, y = value, fill = celltype)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    textsize = 2,
    vjust = 0.1  
  )

## pseudobulks -----

obj <- readRDS("/home/data/zhuang/CT/3.1.PCa/PCa_merged_with_Ctrl_annotated.Rds") 
gene_data <- read_excel("2.1.1PCa_by_c/打分/gmt/selected/Apoptosis  selected.xlsx", col_names = TRUE, sheet = "Sheet1")
head(gene_data)
genesets <- lapply(gene_data, function(col) {
  col[!is.na(col) & col != ""] 
})

Idents(obj) <- obj$celltype
expr <- AverageExpression(obj, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] 
expr <- as.matrix(expr)
head(expr)

gsvaPar <- gsvaParam(expr, genesets,maxDiff = TRUE)
gsvaPar 
gsva.res <- gsva(gsvaPar)
dim(gsva.res)

gsva_df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
rownames(gsva_df) <- gsva_df$Genesets
gsva_df <- gsva_df[, -1] 

gsva_df <- t(apply(gsva_df, 1, rescale, to=c(-1, 1)))
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))

new_column_order <- c("General", "Naive-Dom", "NCHT-General","NCHT-nR-Dom", "NCHT-R-Dom")
gsva_df <- gsva_df[, new_column_order]

pheatmap(gsva_df,
         show_colnames = T,
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = F,
         color = my.colors,
         border_color = "NA",
         fontsize = 8)

# ssGSEA ------------------------------------------------------------------

obj <- readRDS("/home/data/zhuang/CT/3.1.PCa/PCa_merged_with_Ctrl_annotated.Rds") 
gene.expr <-  GetAssayData(obj, assay = "RNA", slot = "data")
gene_data <- read_excel("2.1.1PCa_by_c/打分/gmt/heparin.xlsx", col_names = TRUE, sheet = "Sheet 1")
head(gene_data)
genesets <- lapply(gene_data, function(col) {
  col[!is.na(col) & col != ""] 
})
names(genesets)

ssGSEA <- GSVA::gsva(gene.expr, 
                     genesets,
                     method='ssgsea',
                     kcdf='Gaussian',
                     verbose = TRUE, 
                     parallel.sz = 20)
saveRDS(ssGSEA, "2.1.1PCa_by_c/打分/HEPARIN_ssGSEA.Rds")


celltype <- obj$celltype %>% as.data.frame()
colnames(celltype) <- "celltype"
table(obj$celltype)
ssGSEA <- ssGSEA %>% as.data.frame() %>% t()
mtx <- cbind(ssGSEA, celltype)
mtx$celltype <- factor(mtx$celltype, levels = c("General", "Naive_Dom", "NCHT_General", "NCHT_nR_Dom", "NCHT_R_Dom"))
mtx <- mtx[order(mtx$celltype),]

data <- melt(mtx)
data$celltype <- factor(data$celltype, 
                        levels = c("General", "Naive_Dom", "NCHT_General", "NCHT_nR_Dom", "NCHT_R_Dom"))

ggplot(data, aes(x = variable, y = value, fill = celltype)) +
  geom_boxplot(width = 0.6, outlier.colour = NA, position=position_dodge(0.8)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, colour="black", size = 8),
        axis.text.y=element_text(size=12, colour="black"), 
        axis.title.y=element_text(size = 14),
        panel.border = element_rect(fill=NA, color="black", linewidth=1, linetype="solid"),
        legend.text=element_text(colour="black", size=8),
        legend.title=element_text(colour="black", size=10),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("HEPARIN pathways") +
  ylab("ssGSEA score") +
  xlab("")

comparisons <- list(c("NCHT_nR_Dom", "NCHT_R_Dom"))

ggplot(data, aes(x = celltype, y = value, fill = celltype)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.5,
    textsize = 2,
    vjust = 0.1
  )

# AddModuleScore ----------------------------------------------------------

obj <- readRDS("3.3.Tcell/CD8/CD8_annotated.Rds")
gene_data <- read_excel("3.3.Tcell/gs.xlsx", col_names = TRUE, sheet = "CD8T")
gene_list <- lapply(gene_data, function(col) {
  col[!is.na(col) & col != ""] 
})

obj <- AddModuleScore(obj,
                      features = gene_list,
                      ctrl = 5,
                      name = "FunctionScore")
for(i in 1:length(gene_list)) {
  colnames(obj@meta.data)[colnames(obj@meta.data) == paste0("FunctionScore", i)] <- names(gene_list)[i]
}

## 热图 -----
Idents(obj) <- obj$celltype_minor
Differentiation <- c("Naive", "Activation:Effector function", "Exhaustion")
Function <- c("TCR Signaling", "Cytotoxicity", "Cytokine:Cytokine receptor",
              "Chemokine:Chemokine receptor", "Senescence", "Anergy",
              "NFKB Signaling", "Stress response", "MAPK Signaling", "Adhesion",
              "IFN Response")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(obj$celltype_minor)),
                              nrow = length(MarkerNameVector))
colnames(FunctionScoreMatrix) <- paste0(unique(obj$celltype_minor))
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)) {
  for(ri in 1:nrow(FunctionScoreMatrix)) {
    FunctionVec <- as_tibble(obj@meta.data) %>% pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[obj$celltype_minor == unique(obj$celltype_minor)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#6DCCFD", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#FD9AA0"))(length(my.breaks)/2))
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector
col_name_v <- paste0(sort(unique(obj$celltype_minor)))
FunctionScoreMatrix <- FunctionScoreMatrix[, col_name_v]

pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         annotation_row = signatureType_row,
         gaps_row = c(3, 14, 17),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 4.5,
         filename = file.path("3.3.Tcell/CD8/注释打分/heatmap.pdf"))

## 小提琴图 -----

X <- "REACTOME_PI3K_AKT_ACTIVATION"
module_score_df <- FetchData(
  obj,
  vars = c(paste0(X), "celltype")
)

module_score_df$celltype <- factor(module_score_df$celltype, 
                                   levels = c("General", "Naive_Dom", "NCHT_General", "NCHT_nR_Dom", "NCHT_R_Dom"))

p <- ggplot(module_score_df, aes(x = celltype, y = !!sym(X), fill = celltype)) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(alpha = 0.1, color = "black", size = 0.05) +
  labs(title = paste0(X," Score by Cell Type"),
       x = "Cell Type",
       y = "Module Score") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  stat_compare_means(
    comparisons = combn(levels(module_score_df$celltype), 2, simplify = FALSE),
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01,
    vjust = 0.5
  );p

# PercentageFeatureSet ----------------------------------------------------

obj <- readRDS("/home/data/zhuang/CT/3.1.PCa/PCa_merged_with_Ctrl_annotated.Rds") 
gene_data <- read_excel("2.1.1PCa_by_c/打分/gmt/selected/Apoptosis  selected.xlsx", col_names = TRUE, sheet = "Sheet1")

gene_columns <- names(gene_data)
celltype_levels <- c("General", "Naive_Dom", "NCHT_General", "NCHT_nR_Dom", "NCHT_R_Dom")
plot_list <- list()
for (col_name in gene_columns) {
  current_genes <- gene_data[[col_name]] %>%
    as.character() %>%
    .[!is.na(.)] %>%
    .[. != ""]
  matching_genes <- intersect(current_genes, rownames(obj))
  if (length(matching_genes) == 0) {
    message("No matching genes found for ", col_name)
    next
  }
  obj <- PercentageFeatureSet(obj, features = matching_genes, col.name = col_name)
  module_score_df <- FetchData(obj, vars = c(col_name, "celltype"))
  module_score_df$celltype <- factor(module_score_df$celltype, levels = celltype_levels)
  p <- ggplot(module_score_df, aes(x = celltype, y = !!sym(col_name), fill = celltype)) +
    geom_violin(scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA) +
    geom_jitter(alpha = 0.1, color = "black", size = 0.05) +
    labs(title = paste0(col_name, " Score by Cell Type"),
         x = "Cell Type",
         y = "Module Score") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none" 
    ) +
    stat_compare_means(
      comparisons = combn(levels(module_score_df$celltype), 2, simplify = FALSE),
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      vjust = 0.5
    )
  plot_list[[col_name]] <- p
}

print(names(plot_list))

pdf("all_gene_modules_plots.pdf", width = 10, height = 8)
for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
}
dev.off()




