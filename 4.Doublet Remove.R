library(DoubletFinder)

PCSC <- readRDS("1.data_processing/Annotation/NoInteg_PCSC_annotated.Rds") 

PCSC_split <- SplitObject(PCSC, split.by = "orig.ident")

for (i in 1:length(PCSC_split)) {
  # pK Identification
  sweep.res.list <- paramSweep(PCSC_split[[i]], PCs = 1:20)
  
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- PCSC_split[[i]]@meta.data$celltype
  homotypic.prop <- modelHomotypic(annotations)  
  DoubletRate = ncol(PCSC_split[[i]])*8*1e-6     #按每增加1000个细胞，双细胞比率增加千分之8来计算
  
  # 估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。
  nExp_poi <- round(DoubletRate*nrow(PCSC_split[[i]]@meta.data))  
  
  # 计算双细胞比例
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # doubletFinder
  PCSC_split[[i]] <- doubletFinder(PCSC_split[[i]], PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, sct = FALSE)
  PCSC_split[[i]] <- doubletFinder(PCSC_split[[i]], PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, sct = FALSE)
  
  # judge
  PCSC_split[[i]]@meta.data[,"Doublet"] <- PCSC_split[[i]]@meta.data[,length(colnames(PCSC_split[[i]]@meta.data))]
  
}

#merge seurat
PCSC$Doublet <- "Singlet"
for (i in 1:length(PCSC_split)) {
  PCSC@meta.data$Doublet[rownames(PCSC@meta.data) %in% 
                           rownames(PCSC_split[[i]]@meta.data)[PCSC_split[[i]]@meta.data$Doublet == "Doublet"]] <- "Doublet"
}


saveRDS(PCSC, "./1.data_processing/Annotation/NoInteg_PCSC_doublet_labeled.Rds")

# PCSC <- readRDS("./1.data_processing/Annotation/NoInteg_PCSC_doublet_labeled.Rds")

Idents(PCSC) <- "Doublet"

pdf("./1.data_processing/Annotation/NoInteg_Umap_Doublet.pdf", width = 5, height = 5)
DimPlot(PCSC,
        label = T,
        pt.size=0.5,
        shuffle = T,
        label.color = "black",
        label.box = T,
        label.size = 5,
        repel = T,
        raster=FALSE) +
  NoLegend() +
  theme(plot.title = element_blank())
dev.off()


doublet_stats <- PCSC@meta.data %>%
  group_by(celltype, Doublet) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(celltype) %>%
  mutate(total = sum(count),
         proportion = count / sum(count)) %>%
  ungroup()

plot_data <- dcast(doublet_stats, celltype ~ Doublet, value.var = "proportion", fill = 0)

plot_data_long <- melt(plot_data, id.vars = "celltype", variable.name = "Type", value.name = "Proportion")

pdf("./1.data_processing/Annotation/NoInteg_Doublet_proportion.pdf", width = 8, height = 5)
ggplot(plot_data_long, aes(x = celltype, y = Proportion, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Singlet vs Doublet Proportion by Cell Type",
       x = "Cell Type",
       y = "Proportion") +
  scale_fill_manual(values = c("Singlet" = "#66C2A5", "Doublet" = "#FC8D62"))
dev.off()

PCSC_doublet_Removed <- subset(PCSC, Doublet != "Doublet")
saveRDS(PCSC_doublet_Removed, "./1.data_processing/Annotation/NoInteg_PCSC_doublet_Removed.Rds")

PCSC_doublet <- subset(PCSC, Doublet != "Singlet")
saveRDS(PCSC_doublet, "./1.data_processing/Annotation/NoInteg_PCSC_doublet.Rds")

# PCSC_doublet <- readRDS("./1.data_processing/Annotation/NoInteg_PCSC_doublet.Rds") 
# table(PCSC_doublet$Group)
# dim(PCSC_doublet)
rm(list = ls());gc()

PCSC <- readRDS("./1.data_processing/Annotation/NoInteg_PCSC_doublet_Removed.Rds") 
# table(PCSC$Group)
# dim(PCSC)

df <- data.frame(Sample = names(table(PCSC$orig.ident)),
                 afterDoubletRemoved_Count = as.numeric(table(PCSC$orig.ident)))
write.xlsx(df, "./1.data_processing/Annotation/cell_counts_afterDoubletRemoved.xlsx")


color <- c("Epithelium"="#D1352B","T cell"="#7DBFA7","B cell"="#3C77AF",
           "Neutrophil"="#D0AFC4","Phagocyte"="#BBDD78","Mast cell"="#C6307C",
           "Endothelium"="#EE934E","Mesenchymal cell"="#c6b598","Neural cell"="#EBCC96")
pdf("./1.data_processing/Annotation/NoInteg_CellCluster_Umap_res1_celltype_doublet_removed.pdf", width = 6, height = 6)
DimPlot(PCSC, 
        label = T, 
        pt.size=0.5, 
        shuffle = T, 
        label.color = "black", 
        label.box = T, 
        label.size = 5,
        repel = T,
        group.by = "celltype",
        raster=FALSE,
        cols = color) + 
  NoLegend() +
  theme(
    plot.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  # 黑色边框
    panel.background = element_rect(fill = "white"),  # 确保背景为白色
    axis.line = element_blank(),         # 移除坐标轴线
    axis.ticks = element_blank(),        # 移除刻度线
    axis.text = element_blank()         # 移除刻度标签
  )
dev.off()

rm(list = ls());gc()

preQC_df <- read.xlsx("./1.data_processing/cell_counts_preQC.xlsx", sheet = "Sheet 1")
afterQC_df <- read.xlsx("./1.data_processing/QC/cell_counts_afterQC.xlsx", sheet = "Sheet 1")
afterDoubletRemoved_df <- read.xlsx("./1.data_processing/Annotation/cell_counts_afterDoubletRemoved.xlsx", sheet = "Sheet 1")

combined_df <- preQC_df %>%
  merge(afterQC_df, by = "Sample") %>%
  merge(afterDoubletRemoved_df, by = "Sample")

colnames(combined_df) <- c("Sample", "Pre_QC", "After_QC", "After_Doublet_Removed")

absolute_df <- combined_df %>%
  melt(id.vars = "Sample", 
       variable.name = "Stage", 
       value.name = "Cell_Count")

absolute_df$Stage <- factor(absolute_df$Stage, 
                            levels = c("Pre_QC", "After_QC", "After_Doublet_Removed"),
                            labels = c("Pre-QC", "After-QC", "After Doublet Removal"))

pdf("./1.data_processing/Annotation/QC_and_DoubletRemoved_Count_Analysis.pdf", width = 8, height = 5)
ggplot(absolute_df, aes(x = Sample, y = Cell_Count, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cell Counts at Different Processing Stages",
       x = "Sample",
       y = "Cell Count",
       fill = "Processing Stage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Pre-QC" = "#a4cde1", "After-QC" = "#f38989", "After Doublet Removal" = "#96cb8f"))
dev.off()

total_counts <- data.frame(
  Stage = c("Pre-QC", "After-QC", "After Doublet Removal"),
  Total_Count = c(sum(preQC_df$preQC_Count),
                  sum(afterQC_df$afterQC_Count),
                  sum(afterDoubletRemoved_df$afterDoubletRemoved_Count))
)

reduction_rates <- data.frame(
  Stage = c("QC Reduction Rate", "Doublet Reduction Rate"),
  Rate = c((total_counts$Total_Count[1] - total_counts$Total_Count[2]) / total_counts$Total_Count[1] * 100,
           (total_counts$Total_Count[2] - total_counts$Total_Count[3]) / total_counts$Total_Count[2] * 100)
)
