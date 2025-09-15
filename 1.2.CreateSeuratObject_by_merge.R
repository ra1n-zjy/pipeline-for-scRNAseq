PCSC_X <- readRDS("/home/data/X/1.data_processing/PCSC_preQC.Rds")
table(PCSC_X$orig.ident)

PCSC_Y <- readRDS("/home/data/Y/1.data_processing/PCSC_preQC.Rds")
table(PCSC_Y$orig.ident)

PCSC_preQC <- merge(x = PCSC_X, y = PCSC_Y, merge.data = TRUE)
table(PCSC_preQC$orig.ident)

PCSC_preQC <- JoinLayers(PCSC_preQC)

folder <- "1.data_processing"
if(!dir.exists(folder)){
  dir.create(folder)
}

PCSC_preQC$Group <- recode(PCSC_preQC$orig.ident,
                           "" = "",
                           "" = "",
                           "" = "",
                           "" = "",
                           "" = ""
)

df <- data.frame(Sample = names(table(PCSC_preQC$orig.ident)), 
                 preQC_Count = as.numeric(table(PCSC_preQC$orig.ident)))

library(openxlsx)
write.xlsx(df, "./1.data_processing/cell_counts_preQC.xlsx")
getSheetNames("./1.data_processing/cell_counts_preQC.xlsx")
df <- read.xlsx("./1.data_processing/cell_counts_preQC.xlsx", sheet = "Sheet 1")


saveRDS(PCSC_preQC,"./1.data_processing/PCSC_preQC.Rds")

