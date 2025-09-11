# load data ------------------------------------------------------------------

folder <- "1.data_processing"
if(!dir.exists(folder)){
  dir.create(folder)
}

filelist <- c("X1","X2","X3_Y1")

for (file in filelist) {
  data <- Read10X(paste0("/home/data/matrix/", file))
  assign(file, data)
}
rm(data)

# add labels for samples and merge data ----------------------------------------------------------------

ls()
data <- X1
for (i in c(X2,X3_Y1)) {
  data <- merge(data, i, by = 0) %>% column_to_rownames("Row.names")
}

ncell <- ncol(X1)
for (i in c(X2,X3_Y1)) {
  ncell <- c(ncell, ncol(i))
}

ncell
sample <- c(rep("X1",ncell[1]),
            rep("X2",ncell[2]),
            rep("X3_Y1",ncell[3]))

colname <- paste0(colnames(data), "_", sample)
colnames(data) <- colname

# create Seurat objects based on merged data ----------------------------------------------------------

PCSC_preQC <- CreateSeuratObject(
  data,
  project = "PCSC_preQC", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")

PCSC_preQC$Group <- "X"
table(PCSC_preQC$orig.ident)
df <- data.frame(Sample = names(table(PCSC_preQC$orig.ident)), 
                 preQC_Count = as.numeric(table(PCSC_preQC$orig.ident)))

library(openxlsx)
write.xlsx(df, "./1.data_processing/cell_counts_preQC.xlsx")
getSheetNames("./1.data_processing/cell_counts_preQC.xlsx")
df <- read.xlsx("./1.data_processing/cell_counts_preQC.xlsx", sheet = "Sheet 1")

saveRDS(PCSC_preQC, file = "./1.data_processing/PCSC_preQC.Rds")

