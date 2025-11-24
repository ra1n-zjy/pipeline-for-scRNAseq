# common data base GSEA ---------------------------------------------------

## Group ---------------------------------------------------

folder <- "GSEA/Group"
if(!dir.exists(folder)){
  dir.create(folder, recursive = TRUE)
}

Idents(MyData) <- "Group"
table(MyData$Group)
markers <- FindMarkers(MyData,
                       assay = "RNA",
                       ident.1 = "AAAA",
                       ident.2 = "BBBB",
                       logfc.threshold = 0,
                       min.pct = 0, 
                       min.diff.pct = -Inf,
                       only.pos = FALSE)
markers <- rownames_to_column(markers, var = "SYMBOL")

write.xlsx(markers, 
           file = "GSEA/Group/markers.xlsx",
           row.names = FALSE)

colnames(markers)[c(1,3)] <- c("SYMBOL", "log2FC")
markers <- markers[,c(1,3)]

head(markers)
# SYMBOL     log2FC
# 1   Cox8a  0.7812288
# 2    Txn1  0.7264555
# 3   Prdx5  1.0425767
# 4 Gm42418 -0.4215925
# 5  Ifitm3  0.9824137
# 6   Fgfr2  0.6557044

markers_id <- bitr(markers$SYMBOL,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db") # Mm
markers_all <- merge(markers, markers_id, by = "SYMBOL", all = FALSE)
markers_all_sort <- markers_all[order(markers_all$log2FC, decreasing = TRUE), ]
genelist <- markers_all_sort$log2FC
names(genelist) <- markers_all_sort$ENTREZID

### GOBP ---------------------------------------------------

go <- gseGO(genelist,
            OrgDb = 'org.Hs.eg.db', # Mm
            ont = 'BP',
            pvalueCutoff = 1,
            minGSSize = 10,
            maxGSSize = 1000,
            eps = 1e-100,
            pAdjustMethod = 'BH',
            keyType = "ENTREZID")

go <- setReadable(go, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
data <- go@result

write.xlsx(data, 
           file = "GSEA/Group/GOBP.xlsx",
           row.names = FALSE)

### KEGG ---------------------------------------------------

kegg <- gseKEGG(genelist,
                organism = "hsa", # mmu
                nPerm = 1000,
                pvalueCutoff = 1,
                minGSSize = 10,
                maxGSSize = 1000,
                eps = 1e-100,
                pAdjustMethod = 'BH',
                keyType = "ncbi-geneid")

kegg <- setReadable(kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
data <- kegg@result

write.xlsx(data, 
           file = "GSEA/Group/KEGG.xlsx",
           row.names = FALSE)

### REACTOME ---------------------------------------------------

reactome <- gsePathway(genelist,
                       organism = "human", # mouse
                       minGSSize = 10,
                       maxGSSize = 1000,
                       pvalueCutoff = 1,
                       pAdjustMethod = "BH",
                       verbose = FALSE,
                       eps = 1e-100)
reactome <- setReadable(reactome, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
data <- reactome@result
write.xlsx(data, 
           file = "GSEA/Group/Reactome.xlsx",
           row.names = FALSE)

### DO ---------------------------------------------------

do <- gseDO(genelist,
            minGSSize = 10,
            maxGSSize = 1000,
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            verbose = FALSE,
            eps = 1e-100)

do <- setReadable(do, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
data <- do@result

write.xlsx(data, 
           file = "GSEA/Group/Do.xlsx",
           rowNames = FALSE)

### HALLMARK ---------------------------------------------------

markers_all_sort <- markers[order(markers$log2FC, decreasing = TRUE), ]
genelist <- markers_all_sort$log2FC
names(genelist) <- markers_all_sort$SYMBOL 
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>% # Mus musculus
  dplyr::select(gs_name, gene_symbol)

hallmark <- GSEA(
  geneList = genelist,
  TERM2GENE = hallmark_gene_sets,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 1000,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  eps = 1e-100,
  seed = 123
)
data <- hallmark@result
write.xlsx(data, 
           file = "GSEA/Group/HALLMARK.xlsx",
           row.names = FALSE)

### visualization, taking GOBP as an example ---------------------------------------------------

data <- read.xlsx("GSEA/Group/GOBP.xlsx")

#### Bar selected pathways ---------------------------------------------------

ID <- c(
  "GO:0045321",
  "GO:0002252",
  "GO:0006954",
  "GO:0002521",
  "GO:0042110",
  "GO:0050900",
  "GO:0042330",
  "GO:0001816",
  "GO:0046649",
  "GO:0050865"
)

colnames(data)
combined <- data[data$ID %in% ID,]
head(combined)

pdf("GSEA/Group/Bar_GOBP_selected.pdf",width = 6,height = 4)
ggplot(combined, aes(reorder(Description, NES), NES)) +
  geom_col(aes(fill = pvalue)) + 
  coord_flip() +
  scale_fill_gradientn(
    colours = c("#F2A49C", "#B9A3C6", "#8DAACB"),
    name = "pvalue",
    limits = c(0.01, 0.05),
    breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
    labels = c("0.01", "0.02", "0.03", "0.04", "0.05")
  )+
  labs(
    x = "",
    y = "Normalized Enrichment Score (NES)",
    title = "GO-BP GSEA"
  ) +
  scale_x_discrete(labels = function(dat) str_wrap(dat, width = 600)) + 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )
dev.off()

#### Bar top10 pathways ---------------------------------------------------

top10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05 & NES > 0) %>% head(n = 10)
bottom10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05 & NES < 0) %>% tail(n = 10)
combined <- bind_rows(top10, bottom10)

pdf("GSEA/Group/Bar_GOBP_nes10.pdf",width = 6,height = 8)
ggplot(combined, aes(reorder(Description, NES), NES)) +
  geom_col(aes(fill = pvalue)) + 
  coord_flip() +
  scale_fill_gradientn(
    colours = c("#F2A49C", "#B9A3C6", "#8DAACB"),  # 蓝-紫-红
    name = "pvalue",
    limits = c(0.01, 0.05),
    breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
    labels = c("0.01", "0.02", "0.03", "0.04", "0.05")
  )+
  labs(
    x = "",
    y = "Normalized Enrichment Score (NES)",
    title = "GO-BP GSEA"
  ) +
  scale_x_discrete(labels = function(dat) str_wrap(dat, width = 600)) + 
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )
dev.off()

#### GSEA enrichment plot ---------------------------------------------------

# selected
ID <- c(
  "GO:0045321",
  "GO:0002252",
  "GO:0006954",
  "GO:0002521",
  "GO:0042110",
  "GO:0050900",
  "GO:0042330",
  "GO:0001816",
  "GO:0046649",
  "GO:0050865"
)

# top
top10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05 & NES > 0) %>% head(n = 10)
ID <- top10$ID

pdf("GSEA/Group/enrichment plot.pdf",width = 4, height =6)

gseaNb(object = go,       
       geneSetID = ID,        
       newGsea = F,        
       addPval = T, 
       lineSize = 1,
       curveCol = c("#80B1D3", "#FF6B6B", "#B3DE69", "#FB8072", "#FCCDE5", "#D9D9D9", "#FFFFCC", "#BC80BD", "#A6CEE3", "#7FC97F"),
       htCol= c( "#80B1D3", "#FF6B6B"),  
       rankCol= c( "#1F618D", "white", "#C0392B"),
       pvalY = 0.7,  
       pvalX = 0.95, 
       markTopgene = F,       
       pHjust = 1, 
       subPlot = 3)

dev.off()

## celltype loop ---------------------------------------------------

celltypes <- unique(MyData$celltype_minor)
folder <- "GSEA/celltype_minor"
if(!dir.exists(folder)){
  dir.create(folder)
}

sub_folders <- c("GOBP", "KEGG", "Reactome", "DO", "HALLMARK")
for (sub_folder_name in sub_folders) {
  full_path <- file.path(folder, sub_folder_name)
  if(!dir.exists(full_path)){
    dir.create(full_path)
  } 
}

for (X in celltypes) {
  
  Idents(MyData) <- "celltype_minor"
  markers <- FindMarkers(MyData,
                         assay = "RNA",
                         ident.1 = X,
                         ident.2 = NULL,
                         logfc.threshold = 0,
                         min.pct = 0, 
                         min.diff.pct = -Inf,
                         only.pos = FALSE)
  markers <- rownames_to_column(markers, var = "SYMBOL")
  colnames(markers)[c(1,3)] <- c("SYMBOL", "log2FC")
  markers <- markers[,c(1,3)]
  markers_id <- bitr(markers$SYMBOL,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = "org.Hs.eg.db") # Mm
  markers_all <- merge(markers, markers_id, by = "SYMBOL", all = FALSE)
  markers_all_sort <- markers_all[order(markers_all$log2FC, decreasing = TRUE), ]
  genelist <- markers_all_sort$log2FC
  names(genelist) <- markers_all_sort$ENTREZID
  
  # GOBP
  go <- gseGO(genelist,
              OrgDb = 'org.Hs.eg.db', # Mm
              ont = 'BP',
              nPerm = 1000,
              pvalueCutoff = 1,
              minGSSize = 10,
              maxGSSize = 1000,
              eps = 1e-100,
              pAdjustMethod = 'BH',
              keyType = "ENTREZID")
  go <- setReadable(go, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
  data <- go@result
  write.xlsx(data, 
             file = paste0("GSEA/celltype_minor/GOBP/GOBP-", X, ".xlsx"),
             row.names = FALSE)
  cat("已完成", X, "的GOBP分析\n")  
  top10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% head(n = 15)
  bottom10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% tail(n = 15)
  combined <- bind_rows(top10, bottom10)
  p <- ggplot(combined, aes(reorder(Description, NES), NES)) +
    geom_col(aes(fill = pvalue)) + 
    coord_flip() +
    scale_fill_gradientn(
      colours = c("#F2A49C", "#B9A3C6", "#8DAACB"),  # 蓝-紫-红
      name = "pvalue",
      limits = c(0.01, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      labels = c("0.01", "0.02", "0.03", "0.04", "0.05")
    )+
    labs(
      x = "",
      y = "Normalized Enrichment Score (NES)",
      title = "GOBP GSEA"
    ) +
    scale_x_discrete(labels = function(dat) str_wrap(dat, width = 600)) + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  ggsave(plot=p, paste0("GSEA/celltype_minor/GOBP/Bar_GOBP_",X,".pdf"),width = 9,height = 8)
  
  # KEGG
  kegg <- gseKEGG(genelist,
                  organism = "hsa", # mmu
                  nPerm = 1000,
                  pvalueCutoff = 1,
                  minGSSize = 10,
                  maxGSSize = 1000,
                  eps = 1e-100,
                  pAdjustMethod = 'BH',
                  keyType = "ncbi-geneid")
  
  kegg <- setReadable(kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
  data <- kegg@result
  write.xlsx(data, 
             file = paste0("GSEA/celltype_minor/KEGG/KEGG-", X, ".xlsx"),
             row.names = FALSE)
  cat("已完成", X, "的KEGG分析\n")  
  top10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% head(n = 15)
  bottom10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% tail(n = 15)
  combined <- bind_rows(top10, bottom10)
  p <- ggplot(combined, aes(reorder(Description, NES), NES)) +
    geom_col(aes(fill = pvalue)) + 
    coord_flip() +
    scale_fill_gradientn(
      colours = c("#F2A49C", "#B9A3C6", "#8DAACB"),  # 蓝-紫-红
      name = "pvalue",
      limits = c(0.01, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      labels = c("0.01", "0.02", "0.03", "0.04", "0.05")
    )+
    labs(
      x = "",
      y = "Normalized Enrichment Score (NES)",
      title = "KEGG GSEA"
    ) +
    scale_x_discrete(labels = function(dat) str_wrap(dat, width = 600)) + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  ggsave(plot=p, paste0("GSEA/celltype_minor/KEGG/Bar_KEGG_",X,".pdf"),width = 9,height = 8)
  
  # Reactome
  reactome <- gsePathway(genelist,
                         organism = "human", # mouse
                         minGSSize = 10,
                         maxGSSize = 1000,
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         verbose = FALSE,
                         eps = 1e-100)
  reactome <- setReadable(reactome, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
  data <- reactome@result
  write.xlsx(data, 
             file = paste0("GSEA/celltype_minor/Reactome/Reactome-", X, ".xlsx"),
             row.names = FALSE)
  cat("已完成", X, "的Reactome分析\n")  
  top10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% head(n = 15)
  bottom10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% tail(n = 15)
  combined <- bind_rows(top10, bottom10)
  p <- ggplot(combined, aes(reorder(Description, NES), NES)) +
    geom_col(aes(fill = pvalue)) + 
    coord_flip() +
    scale_fill_gradientn(
      colours = c("#F2A49C", "#B9A3C6", "#8DAACB"),  # 蓝-紫-红
      name = "pvalue",
      limits = c(0.01, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      labels = c("0.01", "0.02", "0.03", "0.04", "0.05")
    )+
    labs(
      x = "",
      y = "Normalized Enrichment Score (NES)",
      title = "Reactome GSEA"
    ) +
    scale_x_discrete(labels = function(dat) str_wrap(dat, width = 600)) + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  ggsave(plot=p, paste0("GSEA/celltype_minor/Reactome/Bar_Reactome_",X,".pdf"),width = 9,height = 8)
  
  # DO
  do <- gseDO(genelist,
              minGSSize = 10,
              maxGSSize = 1000,
              pvalueCutoff = 1,
              pAdjustMethod = "BH",
              verbose = FALSE,
              eps = 1e-100)
  do <- setReadable(do, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
  data <- do@result
  write.xlsx(data, 
             file = paste0("GSEA/celltype_minor/DO/DO-", X, ".xlsx"),
             row.names = FALSE)
  cat("已完成", X, "的DO分析\n")  
  top10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% head(n = 15)
  bottom10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% tail(n = 15)
  combined <- bind_rows(top10, bottom10)
  p <- ggplot(combined, aes(reorder(Description, NES), NES)) +
    geom_col(aes(fill = pvalue)) + 
    coord_flip() +
    scale_fill_gradientn(
      colours = c("#F2A49C", "#B9A3C6", "#8DAACB"),  # 蓝-紫-红
      name = "pvalue",
      limits = c(0.01, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      labels = c("0.01", "0.02", "0.03", "0.04", "0.05")
    )+
    labs(
      x = "",
      y = "Normalized Enrichment Score (NES)",
      title = "DO GSEA"
    ) +
    scale_x_discrete(labels = function(dat) str_wrap(dat, width = 600)) + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  ggsave(plot=p, paste0("GSEA/celltype_minor/DO/Bar_DO_",X,".pdf"),width = 9,height = 8)
  
  # HALLMARK
  markers_all_sort <- markers[order(markers$log2FC, decreasing = TRUE), ]
  genelist <- markers_all_sort$log2FC
  names(genelist) <- markers_all_sort$SYMBOL 
  hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%  # Mus musculus
    dplyr::select(gs_name, gene_symbol)
  
  hallmark <- GSEA(
    geneList = genelist,
    TERM2GENE = hallmark_gene_sets,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 1000,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 1e-100,
    seed = 123
  )
  data <- hallmark@result
  write.xlsx(data, 
             file = paste0("GSEA/celltype_minor/HALLMARK/HALLMARK-", X, ".xlsx"),
             row.names = FALSE)
  cat("已完成", X, "的HALLMARK分析\n")  
  top10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% head(n = 15)
  bottom10 <- data %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pvalue < 0.05) %>% tail(n = 15)
  combined <- bind_rows(top10, bottom10)
  p <- ggplot(combined, aes(reorder(Description, NES), NES)) +
    geom_col(aes(fill = pvalue)) + 
    coord_flip() +
    scale_fill_gradientn(
      colours = c("#F2A49C", "#B9A3C6", "#8DAACB"),  # 蓝-紫-红
      name = "pvalue",
      limits = c(0.01, 0.05),
      breaks = c(0.01, 0.02, 0.03, 0.04, 0.05),
      labels = c("0.01", "0.02", "0.03", "0.04", "0.05")
    )+
    labs(
      x = "",
      y = "Normalized Enrichment Score (NES)",
      title = "HALLMARK GSEA"
    ) +
    scale_x_discrete(labels = function(dat) str_wrap(dat, width = 600)) + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  ggsave(plot=p, paste0("GSEA/celltype_minor/HALLMARK/Bar_HALLMARK_",X,".pdf"),width = 9,height = 8)
}

# single signature GSEA ---------------------------------------------------

Idents(MyData) <- "Group"
table(MyData$Group)

markers <- FindMarkers(
  object = MyData,
  ident.1 = "AAAA",
  ident.2 = "BBBB",
  only.pos = FALSE,  
  min.diff.pct = -Inf,
  min.pct = 0, 
  logfc.threshold = 0)

markers <- rownames_to_column(markers, var = "SYMBOL")
colnames(markers)[c(1,3)] <- c("SYMBOL", "log2FC")
markers <- markers[,c(1,3)]
markers_sort <- markers[order(markers$log2FC, decreasing = TRUE), ]
genelist <- markers_sort$log2FC
names(genelist) <- markers_sort$SYMBOL
head(genelist)

geneSet <- read.gmt("/home/data/zhuang/AFFAR_YY1_TARGETS_UP.v2025.1.Mm.gmt")

gsea <- GSEA(genelist, 
             TERM2GENE = geneSet, 
             exponent = 1, 
             minGSSize = 5, 
             maxGSSize = 1000,
             eps = 1e-10, 
             pvalueCutoff = 1, 
             pAdjustMethod = "BH",
             by = "fgsea"
)

pdf("AFFAR_YY1_TARGETS_UP.pdf",width = 5, height = 4)

gseaNb(object = gsea,       
       geneSetID = "AFFAR_YY1_TARGETS_UP",        
       newGsea = F,        
       addPval = T, 
       lineSize = 1,
       curveCol= c("#80B1D3", "#FF6B6B"),  
       htCol= c( "#80B1D3", "#FF6B6B"),  
       rankCol= c( "#1F618D", "white", "#C0392B"),
       pvalY = 0.7,  
       pvalX = 0.95, 
       markTopgene = F,       
       pHjust = 1, 
       subPlot = 3)

dev.off()















