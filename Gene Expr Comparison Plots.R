# violin + pie ------------------------------------------------------------

# function main

ks_VlnExp <- function(object,
                      group,
                      group_order,
                      features, # only supports a single gene
                      comparisons,
                      cols=NULL,
                      pieSize=NULL){
  
  #seurat obj
  meta <- object@meta.data
  if(group %in% colnames(meta)){
    
    Idents(object) <- group
    Idents(object) <- factor(Idents(object), levels = group_order)
    
  }else{
    
    print("Your group name does not exist, Please provide correct name in the colnames of metadata")
    
  }
  
  #gene expression percentage in cells
  
  # if(length(group_order) == 2){
  
  # gene_per <- FindMarkers(object, features = features,group.by = group,
  #                          ident.1 = group_order[1], ident.2 = group_order[2]) %>% .[,c("pct.1","pct.2")] %>% t()%>% as.data.frame()
  # 
  # colnames(gene_per) <- "gene"
  # 
  # gene_per$non_gene <- 1-gene_per$gene 
  # rownames(gene_per) <- group_order
  
  exp = FetchData(object, vars = features)
  colnames(exp) <- 'gene'
  exp$Group <- object@meta.data[,group]
  
  fearures_exp <- list()
  for (i in 1:length(group_order)) {
    
    exp1 = subset(exp, Group==group_order[i])
    colnames(exp1)[1] <- "gene"
    fearures_exp[[i]] <- exp1
    names(fearures_exp)[i] <- group_order[i]
    
  }
  
  pct_list <- list()
  
  for (i in 1:length(group_order)) {
    
    pct <- length(fearures_exp[[i]]$gene[fearures_exp[[i]]$gene>0]) / length(fearures_exp[[i]]$gene) %>% as.data.frame()
    pct_list[[i]] <- pct
    names(pct_list)[i] <- group_order[i]
    
  }
  
  exp_pct <- do.call(rbind, pct_list)
  colnames(exp_pct) <- "exp_gene"
  exp_pct$non_exp <- 1-exp_pct$exp_gene
  
  print(exp_pct)
  
  #plot pie
  
  if(is.null(cols)){
    
    colors_map = c("#FF5744","#208A42","#F98400", "#5BBCD6",'#7F3C8D' ,'#11A579', '#3969AC','#E73F74')
    cols <- colors_map[1:length(group_order)]
    
  }else{
    
    cols <- cols
  }
  
  plot_pie <- function(i) {
    
    df1 <- gather(exp_pct[i,], type, value, 1:2)
    df1$labels <- df1$value
    df1$labels <- scales::percent(df1$value,accuracy=0.1)
    df1$labels[2] <- ''
    
    ggplot(df1, aes(x= '', value, fill=type)) +
      geom_col(color='black') + #饼图边设置为黑色
      coord_polar(theta = 'y') +
      theme_void() + 
      theme_transparent() +
      theme(legend.position = "none")+
      geom_text(aes(label = labels), 
                position = position_stack(vjust = 0.5), size=2)+
      scale_fill_manual(values = c(cols[i],"grey80"))#修改饼图填充颜色
  }
  
  #data for pie
  
  exp_pct$pie <- lapply(1:nrow(exp_pct), plot_pie)
  
  exp_pct$x <- 1:length(group_order)
  exp_pct$y <- max(exp$gene)+2
  
  if(is.null(pieSize)){
    
    pieSize <- 1.5
  }else{
    pieSize <-pieSize
    
  }
  
  exp_pct$width <- pieSize
  exp_pct$height <-pieSize
  
  #plot Vln
  
  label.y.pos =  seq(max(exp$gene), max(exp$gene)+100, by = 0.3)
  
  p = VlnPlot(object, features = features)&
    theme_bw()&
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = 'black',face = "bold", size = 8),
          axis.text.y = element_text(color = 'black', face = "bold"),
          axis.title.y = element_text(color = 'black', face = "bold", size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black",size = 1.2, linetype="solid"),#修改边框大小
          panel.spacing = unit(0.12, "cm"),
          plot.title = element_text(hjust = 0.5, face = "bold.italic"),#字体加粗斜体
          legend.position = 'none')&
    stat_compare_means(method="t.test",hide.ns = F, #显著性检验
                       comparisons = comparisons,
                       label="p.signif",
                       bracket.size=0.8,
                       tip.length=0,
                       size=5,
                       vjust = 0.6,
                       label.y = label.y.pos)&
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))& #修改y轴，让其适应显著性检验标注
    scale_fill_manual(values = cols)&
    geom_subview(aes(x=x, y=y, 
                     subview=pie,width = width,
                     height = height), data=exp_pct)&
    ylim(0,max(exp$gene)+3)
  
  return(p)
  
}

my_comparisons <- list(c("NCHT", "Ctrl"))

pdf("XXX.pdf",width=3,height = 4)

ks_VlnExp(object = MyData, group = "Group",
          group_order = c("NCHT", "Ctrl"),
          features = "IL1RAP", comparisons = list())

dev.off()


# violin facet_wrap + statistical test ------------------------------------

markers <- c(
  "AAAA",
  "BBBB"
)
merged_df <- FetchData(
  MyData,
  vars = c(markers, "celltype_minor")
)

merged_df_long <- merged_df %>%
  pivot_longer(cols = -celltype_minor, names_to = "Gene", values_to = "Expression")

merged_df_long$Gene <- factor(
  merged_df_long$Gene,
  levels = markers 
)

color <- c(""="#D1352B", "" = "#f38989", "" = "#735e79", "" = "#555eAA", "Stress_Mast" = "#277fb8")

p <- ggplot(merged_df_long, aes(x = celltype_minor, y = Expression, fill = celltype_minor)) +
  geom_violin(scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.6, outlier.shape = NA) +
  # geom_jitter(alpha = 0.1, color = "black", size = 0.05) +
  scale_fill_manual(values = color) +
  labs(title = "Gene Expression by celltype",
       x = "Cell Type",
       y = "Expression") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank()
  ) +
  facet_wrap(~ Gene, scales = "free_y") +
  stat_compare_means(
    comparisons = combn(levels(merged_df_long$celltype_minor), 2, simplify = FALSE),
    method = "wilcox.test",  
    label = "p.signif",      
    tip.length = 0.01,
    vjust = 0.5
  );p

ggsave(p,file="XXX.pdf",width = 40,height = 20)


# dot ---------------------------------------------------------------------

jjDotPlot(MyData,
          gene = markers,
          id = "seurat_clusters",
          #dot.col = c("white", "#A52A2A"),
          ytree = F,
          xtree = F,
          legend.position = "bottom",
          plot.margin = c(0.5,0.5,0.5,0.5),
          dot.min = 0,
          dot.max = 6,
          x.text.vjust = 0.5)+ 
  ggsci::scale_fill_gsea()  # 添加GSEA配色



