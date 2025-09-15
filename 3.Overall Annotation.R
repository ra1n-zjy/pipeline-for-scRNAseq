folder <- "1.data_processing/Annotation"
if(!dir.exists(folder)){
  dir.create(folder)
}

# Cell cycle scoring ---------------------------------------------------------

PCSC_NoInteg <- readRDS("./1.data_processing/QC/PCSC_afterQC.Rds")
table(PCSC_NoInteg$Group)
dim(PCSC_NoInteg)

PCSC_NoInteg <- NormalizeData(PCSC_NoInteg)
str(cc.genes)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
PCSC_NoInteg <- CellCycleScoring(PCSC_NoInteg,
                                 g2m.features = g2m.genes,
                                 s.features = s.genes)

# clustering and annotation without harmony ---------------------------------------------

## clustering ---------------------------------------------

PCSC_NoInteg <- FindVariableFeatures(PCSC_NoInteg,
                                     selection.method = "vst",
                                     nfeatures = 2000,
                                     verbose = FALSE) %>%
  ScaleData(vars.to.regress = c("percent.mt")) %>% 
  RunPCA(npcs = 50)

ElbowPlot(PCSC_NoInteg, ndims=50)
PCSC_NoInteg <- RunUMAP(PCSC_NoInteg, dims = 1:20)

pdf("./1.data_processing/Annotation/NoInteg_Umap_SampleGroup.pdf", width = 8, height = 7)
DimPlot(PCSC_NoInteg,
        reduction = "umap",
        group.by= "orig.ident")
dev.off()

PCSC_NoInteg <- FindNeighbors(PCSC_NoInteg, reduction = "pca", dims = 1:20  )
for (res in c(0.01, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1,1.2,1.5)) {
  PCSC_NoInteg <- FindClusters(PCSC_NoInteg, resolution = res)
}

apply(PCSC_NoInteg@meta.data[, grep("RNA_snn_res", colnames(PCSC_NoInteg@meta.data))],2,table)

library(clustree)

p2_tree=clustree(PCSC_NoInteg@meta.data, prefix = "RNA_snn_res.") 
ggsave(plot=p2_tree, filename="./1.data_processing/Annotation/NoInteg_Tree_diff_resolution.pdf",width = 10, height = 12)
p2_tree

res_columns <- grep("RNA_snn_res", colnames(PCSC_NoInteg@meta.data), value = TRUE)
plot_list <- lapply(res_columns, function(col) {
  DimPlot(PCSC_NoInteg, group.by = col, label = TRUE, raster = T) + 
    NoLegend() + 
    ggtitle(paste("Resolution:", gsub("RNA_snn_res\\.", "", col)))  # 添加分辨率标题
})
combined_plot <- wrap_plots(plot_list, ncol = 5) 
ggsave(plot = combined_plot, "./1.data_processing/Annotation/Nomix_DimPlot_cluster_res.pdf", width = 15, height = 7 )



## annotation --------------------------------------------------------------

Idents(PCSC_NoInteg) <- "RNA_snn_res.1"

pdf("./1.data_processing/Annotation/NoInteg_Cluster_Umap_res1.pdf",width = 7,height = 7)
DimPlot(PCSC_NoInteg, 
        reduction = "umap", 
        raster=FALSE,
        pt.size = 0.1, 
        label = T) + 
  NoLegend() +
  theme(plot.title = element_blank())
dev.off()

Allmarker=c("EPCAM","CDH1", 
            "VIM","PECAM1", "ENG", 
            "COL3A1", "FAP", "PDGFRA", 
            "THY1", "ACTA2", "MYL9",
            "PTPRC", 
            "CD3E", "CD3G", 
            "MS4A1", "CD79A", 
            "CD14", "FCGR3A", "CD68", "CD163", "LYZ",    
            "KIT", "MS4A2","TPSAB1","TPSB2","CPA3",
            "FCGR3B", "CSF3R","CXCR1",
            "S100B","PLP1","NRXN1",
            "TOP2A","MKI67")

library(scRNAtoolVis)
pdf("./1.data_processing/Annotation/NoInteg_AllmarkerBubble_res1.pdf",width = 10, height = 12)
jjDotPlot(PCSC_NoInteg,
          gene = Allmarker,
          id = "RNA_snn_res.1",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 6,
          x.text.vjust = 0.5)
dev.off()

celltype=data.frame(ClusterID=0:38, celltype='unknow')   
celltype[celltype$ClusterID %in% c(7,10,16,18,21,22,31),2]='Epithelium'
celltype[celltype$ClusterID %in% c(5,6,9,19,23,26,38),2]='Endothelium'
celltype[celltype$ClusterID %in% c(3,4,8,14,20),2]='Mesenchymal cell'
celltype[celltype$ClusterID %in% c(0,1,2,12,24,29),2]='T cell'  
celltype[celltype$ClusterID %in% c(17,30),2]='B cell' 
celltype[celltype$ClusterID %in% c(11,13,25,32,33,36),2]='Phagocyte' 
celltype[celltype$ClusterID %in% c(15,34,35,37),2]='Mast cell' 
celltype[celltype$ClusterID %in% c(27),2]='Neutrophil' 
celltype[celltype$ClusterID %in% c(28),2]='Neural cell'
head(celltype)

PCSC_NoInteg@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  PCSC_NoInteg@meta.data[which(PCSC_NoInteg@meta.data$RNA_snn_res.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(PCSC_NoInteg@meta.data$celltype)

Idents(PCSC_NoInteg) <- "celltype"

color <- c("Epithelium"="#D1352B","T cell"="#7DBFA7","B cell"="#3C77AF",
           "Neutrophil"="#D0AFC4","Phagocyte"="#BBDD78","Mast cell"="#C6307C",
           "Endothelium"="#EE934E","Mesenchymal cell"="#c6b598","Neural cell"="#EBCC96")

pdf("./1.data_processing/Annotation/NoInteg_CellCluster_Umap_res1_celltype.pdf", width = 5, height = 5)
DimPlot(PCSC_NoInteg, 
        label = T, 
        pt.size=0.5, 
        shuffle = T, 
        label.color = "black", 
        label.box = T,  
        label.size = 5,
        repel = T,
        raster=FALSE,
        cols = color) + 
  NoLegend() +
  theme(plot.title = element_blank())
dev.off()

PCSC_NoInteg$celltype <- factor(PCSC_NoInteg$celltype, levels = c('Epithelium','Endothelium','Mesenchymal cell','T cell','B cell','Phagocyte','Mast cell','Neutrophil','Neural cell'))

pdf("./1.data_processing/Annotation/NoInteg_AllmarkerBubble_celltype.pdf",width = 10, height = 8)
jjDotPlot(PCSC_NoInteg,
          gene = Allmarker,
          id = "celltype",
          dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "top",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 6,
          x.text.vjust = 0.5)
dev.off()

saveRDS(PCSC_NoInteg, "1.data_processing/Annotation/NoInteg_PCSC_annotated.Rds")


