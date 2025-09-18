PCSC_preQC <- readRDS("./1.data_processing/PCSC_preQC.Rds")

folder <- "1.data_processing/QC"
if(!dir.exists(folder)){
  dir.create(folder)
}

# cell counts --------------------------------------------------------------------

before_sep <- df

pdf("./1.data_processing/QC/preQC_Count.pdf")
count_before <- PCSC_preQC@meta.data %>%
  ggplot(aes(x=orig.ident, fill=orig.ident)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells"); count_before
dev.off()

# nCount_RNA ----------------------------------------------------------------

pdf("./1.data_processing/QC/preQC_Count_Density.pdf")
count_den <- PCSC_preQC@meta.data %>%
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 1000); count_den
dev.off()

# nFeature_RNA----------------------------------------------------------

pdf("./1.data_processing/QC/preQC_Feature_Density.pdf")
feature_den <- PCSC_preQC@meta.data %>%
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ylab("Cell density") +
  scale_x_log10() +
  geom_vline(xintercept = 500); feature_den
dev.off()

# mito percentage, nCount vs. nFeature ----------------------------------------------------

PCSC_preQC[["percent.mt"]] <- PercentageFeatureSet(PCSC_preQC,pattern = "^MT-")

pdf("./1.data_processing/QC/preQC_Mt_Density.pdf")
mt_den <- PCSC_preQC@meta.data %>%
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ylab("Cell density")  +
  geom_vline(xintercept = 10); mt_den
dev.off()

pdf("./1.data_processing/QC/preQC_Feature_Count.pdf")
PCSC_preQC@meta.data %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~orig.ident)
dev.off()

# number of genes detected per UMI ----------------------------------------------------------------

PCSC_preQC$log10GenesPerUMI <- log10(PCSC_preQC$nFeature_RNA)/log10(PCSC_preQC$nCount_RNA)

pdf("./1.data_processing/QC/preQC_perUMI_Density.pdf")
perUMI_den <- PCSC_preQC@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 0.8); perUMI_den
dev.off()

library(gridExtra)
pdf("./1.data_processing/QC/preQC_Density.pdf", width = 12, height = 10)  
grid.arrange(count_den,feature_den,mt_den,perUMI_den)
dev.off() 

# VlnPlot ------------------------------------------------------------
P1 <- PCSC_preQC@meta.data %>%
  ggplot(aes(x=orig.ident, y=nFeature_RNA, fill=orig.ident)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature") +
  theme(legend.position = "none") 

P2 <- PCSC_preQC@meta.data %>%
  ggplot(aes(x=orig.ident, y=nCount_RNA, fill=orig.ident)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount") +
  theme(legend.position = "none") 

P3 <- PCSC_preQC@meta.data %>%
  ggplot(aes(x=orig.ident, y=percent.mt, fill=orig.ident)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("percent.mt") +
  theme(legend.position = "none") 

P4 <- PCSC_preQC@meta.data %>%
  ggplot(aes(x=orig.ident, y=log10GenesPerUMI, fill=orig.ident)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("log10GenesPerUMI") +
  theme(legend.position = "none")

pdf("./1.data_processing/QC/preQC_VlnPlot.pdf",width = 15,height = 8)  
grid.arrange(P1,P2,P3,P4,nrow=2,ncol=2)
dev.off() 

# distribution of nFeature, nCount, mito percentage and log10GenesPerUMI ---------------------------------------------------

# per 5 percentage
gene.freq <- do.call("cbind", tapply(PCSC_preQC@meta.data$nFeature_RNA,PCSC_preQC@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(PCSC_preQC@meta.data$nCount_RNA,PCSC_preQC@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(PCSC_preQC@meta.data$percent.mt,PCSC_preQC@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
GPU.freq <- do.call("cbind", tapply(PCSC_preQC@meta.data$log10GenesPerUMI,PCSC_preQC@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq,GPU.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"),
                            paste(colnames(GPU.freq),"Gene","per","UMI",sep = "_"))

write.xlsx(freq.combine, "./1.data_processing/QC/preQC_freq.xlsx")
freq.combine <- read.xlsx("./1.data_processing/QC/preQC_freq.xlsx", sheet = "Sheet 1")
rm(gene.freq,rna.freq,mt.freq,GPU.freq)
View(freq.combine)

# mito percentage vs. nCount & nGene vs. nCount

plot1 <- FeatureScatter(PCSC_preQC, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle = T, pt.size = 0.5) + 
  NoLegend()
plot2 <- FeatureScatter(PCSC_preQC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle = T, pt.size = 0.5) + 
  NoLegend()
pdf("./1.data_processing/QC/preQC_FeatureScatter.pdf", width = 6, height = 4.5)
grid.arrange(plot1,plot2,nrow=1,ncol=2)
dev.off()
rm(plot1,plot2)

# filter cells and genes -----------------------------------------------

before <- length(rownames(PCSC_preQC@meta.data))

PCSC_afterQC <- subset(PCSC_preQC, 
                       subset = 
                         nFeature_RNA > 500 & 
                         nCount_RNA > 1000 &
                         percent.mt < 10 &
                         log10GenesPerUMI > 0.80 &
                         nCount_RNA < 50000)

gene.freq <- do.call("cbind", tapply(PCSC_afterQC@meta.data$nFeature_RNA,PCSC_afterQC@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(PCSC_afterQC@meta.data$nCount_RNA,PCSC_afterQC@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
mt.freq <- do.call("cbind", tapply(PCSC_afterQC@meta.data$percent.mt,PCSC_afterQC@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
GPU.freq <- do.call("cbind", tapply(PCSC_afterQC@meta.data$log10GenesPerUMI,PCSC_afterQC@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,mt.freq,GPU.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(mt.freq),"MT",sep = "_"),
                            paste(colnames(GPU.freq),"Gene","per","UMI",sep = "_"))

write.xlsx(freq.combine, "./1.data_processing/QC/afterQC_freq.xlsx")
freq.combine <- read.xlsx("./1.data_processing/QC/afterQC_freq.xlsx", sheet = "Sheet 1")

after <- length(rownames(PCSC_afterQC@meta.data))

df <- data.frame(Sample = names(table(PCSC_afterQC$orig.ident)),
                 afterQC_Count = as.numeric(table(PCSC_afterQC$orig.ident)))
write.xlsx(df, "./1.data_processing/QC/cell_counts_afterQC.xlsx")
df <- read.xlsx("./1.data_processing/QC/cell_counts_afterQC.xlsx", sheet = "Sheet 1")

preQC_df <- read.xlsx("./1.data_processing/cell_counts_preQC.xlsx", sheet = "Sheet 1")
afterQC_df <- read.xlsx("./1.data_processing/QC/cell_counts_afterQC.xlsx", sheet = "Sheet 1")

combined_df <- merge(preQC_df, afterQC_df, by = "Sample")

combined_df <- combined_df %>%
  mutate(Change = afterQC_Count - preQC_Count,
         Change_Percent = round((Change / preQC_Count) * 100, 2))

write.xlsx(combined_df, "./1.data_processing/QC/QC_change.xlsx")
combined_df <- read.xlsx("./1.data_processing/QC/QC_change.xlsx", sheet = "Sheet 1")

combined_df_plot <- combined_df %>%
  mutate(Sample = reorder(Sample, Change))  # 保持排序一致



p1 <- ggplot(combined_df_plot, aes(x = Sample)) +
  geom_col(aes(y = Change, fill = ifelse(Change < 0, "Decrease", "Increase")), 
           alpha = 0.7) +
  scale_fill_manual(values = c("Decrease" = "coral", "Increase" = "lightblue"), 
                    name = "Change Direction") +
  geom_point(aes(y = rescale(Change_Percent, to = range(Change)), 
                 color = Change_Percent), 
             size = 2) +
  scale_color_gradient(low = "red", high = "blue", name = "Percent Change") +
  coord_flip() +
  labs(title = "QC Cell Count Analysis",
       x = "Sample",
       y = "Absolute Change in Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))

combined_long <- combined_df_plot %>%
  pivot_longer(cols = c(preQC_Count, afterQC_Count), 
               names_to = "QC_Status", 
               values_to = "Count") %>%
  mutate(QC_Status = ifelse(QC_Status == "preQC_Count", "Before QC", "After QC"))

p2 <- ggplot(combined_long, aes(x = Sample, y = Count, fill = QC_Status)) +
  geom_col(position = "dodge", width = 0.7) +
  labs(title = "Pre-QC vs After-QC Count Comparison",
       subtitle = "Comparison of read counts before and after quality control",
       x = "Sample",
       y = "Read Count",
       fill = "QC Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Before QC" = "#a4cde1", "After QC" = "#f38989"))

combined_long <- combined_df_plot %>%
  pivot_longer(cols = c(preQC_Count, afterQC_Count), 
               names_to = "QC_Status", 
               values_to = "Count") %>%
  mutate(QC_Status = ifelse(QC_Status == "preQC_Count", "Before QC", "After QC"))

pdf("./1.data_processing/QC/QC_Count_Analysis.pdf",width = 8,height = 10)  
grid.arrange(p1, p2, nrow = 2,ncol = 1)
dev.off() 

saveRDS(PCSC_afterQC, "./1.data_processing/QC/PCSC_afterQC.Rds")
PCSC_afterQC <- readRDS("./1.data_processing/QC/PCSC_afterQC.Rds")







