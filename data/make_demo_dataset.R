
################################## Downloading TCGA COAD DATA #################################################

# Adapted from https://github.com/swayamjk10/TCGA-Data-Analysis-with-R/blob/main/TCGABiolinks_data_download_analysis_using_R.R
# References https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

# Set working directory to a very short path (must be close to root, or GDCdownload will fail)

# get list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-COAD')


# build a query to retrieve gene expression data ------------
query_RNA_TCGA <- GDCquery(project = 'TCGA-COAD',
                           data.category = 'Transcriptome Profiling',
                           access = 'open',
                           data.type = "Gene Expression Quantification")


unique(output_query_TCGA$analysis_workflow_type) 

query_miRNA_TCGA <- GDCquery(project = 'TCGA-COAD',
                             data.category = 'Transcriptome Profiling',
                             access = 'open',
                             data.type = "miRNA Expression Quantification")

# download data - GDCdownload
GDCdownload(query_RNA_TCGA)
GDCdownload(query_miRNA_TCGA)

# View data downloaded 
tcga_miRNA_data = GDCprepare(query_miRNA_TCGA)
tcga_RNA_data <- GDCprepare(query_RNA_TCGA)

### Save these objects ---

dim(assay(tcga_RNA_data))     # gene expression matrices. 524 samples and 60660 genes
head(assay(tcga_RNA_data)[,1:10]) # expression of first 6 genes and first 10 samples
head(rowData(tcga_RNA_data))

saveRDS(object = tcga_RNA_data,
        file = "tcga_coad_mRNA_data.RDS",
        compress = FALSE)

saveRDS(object = tcga_miRNA_data,
        file = "tcga_coad_miRNA_data.RDS",
        compress = FALSE)

### Load RDS objects ---
tcga_RNA_data <- readRDS("tcga_coad_mRNA_data.RDS")
tcga_miRNA_data <- readRDS("tcga_coad_miRNA_data.RDS")

### Wrangle metadata ---
colnames(colData(tcga_RNA_data))
table(tcga_RNA_data@colData$paper_vital_status)
table(tcga_RNA_data@colData$tissue_type)
table(tcga_RNA_data@colData$specimen_type)
tcga_RNA_data@colData$bcr_patient_barcode

# Extract metadata
metadata <- colData(tcga_RNA_data)

treatment_info <- data.frame(unlist(metadata$treatments))

metadata$treatments[[1]]$treatment_type
metadata$treatments[[300]]$treatment_type

specific_metadata <- data.frame(
  barcode = metadata$barcode,
  bcr_patient_barcode = metadata$bcr_patient_barcode,
  primary_diagnosis = metadata$primary_diagnosis,
  tissue_type = metadata$tissue_type,
  sex = metadata$gender
)

specific_metadata$base_id <- sapply(strsplit(specific_metadata$barcode, "-"), function(x) paste(x[1:3], collapse = "-"))


final_meta <- specific_metadata


### Wrangle count data --
mrna_df <- assay(tcga_RNA_data)
colnames(mrna_df) <- sapply(strsplit(colnames(mrna_df), "-"), function(x) paste(x[1:4], collapse = "-"))
colnames(mrna_df)
mrna_df <- mrna_df[,colnames(mrna_df) %in% final_meta$bcr_patient_barcode]
dim(mrna_df)


filtered_mirna_df <- tcga_miRNA_data %>% 
  column_to_rownames("miRNA_ID") %>%
  dplyr::select(starts_with("read_count"))

colnames(filtered_mirna_df) <- gsub("read_count_","",colnames(filtered_mirna_df))
colnames(filtered_mirna_df) <- sapply(strsplit(colnames(filtered_mirna_df), "-"), function(x) paste(x[1:4], collapse = "-"))
filtered_mirna_df <- filtered_mirna_df[,colnames(filtered_mirna_df) %in% final_meta$bcr_patient_barcode]

### one last filter to match mrna and mirna samples 
mrna_df <- mrna_df[,intersect(colnames(mrna_df),colnames(filtered_mirna_df))]
dim(mrna_df)

#length(unique(colnames(mrna_df)))

mrna_unique <- mrna_df
dim(mrna_unique)

filtered_mirna_df <- filtered_mirna_df[,colnames(filtered_mirna_df) %in% colnames(mrna_unique)]
filtered_mirna_df <- filtered_mirna_df[,colnames(mrna_unique)]

# check for equality 
colnames(filtered_mirna_df) == colnames(mrna_unique)


# ensure metadata is also in the same order 
final_meta <- final_meta %>% filter(bcr_patient_barcode %in% colnames(filtered_mirna_df)) %>% 
  dplyr::select(-c("barcode")) %>% 
  unique()
final_meta$bcr_patient_barcode <- factor(final_meta$bcr_patient_barcode, levels =colnames(filtered_mirna_df))



# Reorder the data frame based on the ordered SampleID
metadata_ordered <- final_meta[order(final_meta$bcr_patient_barcode), ]
metadata_ordered$SampleID <- gsub("-",".",metadata_ordered$bcr_patient_barcode)
metadata_ordered <- metadata_ordered %>% dplyr::select(SampleID, everything())

colnames(filtered_mirna_df) <- gsub("-", ".", colnames(filtered_mirna_df))
colnames(mrna_unique) <- gsub("-",".", colnames(mrna_unique))

# Add sequencing depth and study id as a variable 
metadata_ordered <- metadata_ordered %>% mutate(Sequencing_Depth = colSums(filtered_mirna_df))
metadata_ordered$TCGA <- "TCGA-COAD"

# Save wrangled count data 
write.csv(mrna_unique, "Tumor_vs_Normal_TCGA_COAD_MRNA_counts.csv")
write.csv(filtered_mirna_df, "Tumor_vs_Normal_TCGA_COAD_miRNA_counts.csv")
write.csv(metadata_ordered, "Tumor_vs_Normal_TCGA_COAD_metadata.csv",row.names=FALSE)


################################## Downloading TCGA READ DATA #################################################

# Adapted from https://github.com/swayamjk10/TCGA-Data-Analysis-with-R/blob/main/TCGABiolinks_data_download_analysis_using_R.R
# References https://www.costalab.org/wp-content/uploads/2020/11/R_class_D3.html
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

#Set working directory to a very short path (must be close to root, or GDCdownload will fail)
getwd()

# get list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-READ')


# build a query to retrieve gene expression data ------------
query_RNA_TCGA <- GDCquery(project = 'TCGA-READ',
                           data.category = 'Transcriptome Profiling',
                           access = 'open',
                           data.type = "Gene Expression Quantification")


query_miRNA_TCGA <- GDCquery(project = 'TCGA-READ',
                             data.category = 'Transcriptome Profiling',
                             access = 'open',
                             data.type = "miRNA Expression Quantification")

# download data - GDCdownload
GDCdownload(query_RNA_TCGA)
GDCdownload(query_miRNA_TCGA)

# View data downloaded 
tcga_miRNA_data = GDCprepare(query_miRNA_TCGA)
tcga_RNA_data <- GDCprepare(query_RNA_TCGA)

### Save these objects ---

dim(assay(tcga_RNA_data))     # gene expression matrices. 524 samples and 60660 genes
head(assay(tcga_RNA_data)[,1:10]) # expression of first 6 genes and first 10 samples
head(rowData(tcga_RNA_data))

saveRDS(object = tcga_RNA_data,
        file = "tcga_read_mRNA_data.RDS",
        compress = FALSE)

saveRDS(object = tcga_miRNA_data,
        file = "tcga_read_miRNA_data.RDS",
        compress = FALSE)

### Load RDS objects if starting from here  ---
tcga_RNA_data <- readRDS("tcga_read_mRNA_data.RDS")
tcga_miRNA_data <- readRDS("tcga_read_miRNA_data.RDS")

### Wrangle metadata ---
colnames(colData(tcga_RNA_data))
table(tcga_RNA_data@colData$paper_vital_status)
table(tcga_RNA_data@colData$tissue_type)
table(tcga_RNA_data@colData$specimen_type)
tcga_RNA_data@colData$bcr_patient_barcode

# Extract metadata
metadata <- colData(tcga_RNA_data)
treatment_info <- data.frame(unlist(metadata$gender))


specific_metadata <- data.frame(
  barcode = metadata$barcode,
  bcr_patient_barcode = metadata$bcr_patient_barcode,
  primary_diagnosis = metadata$primary_diagnosis,
  tissue_type = metadata$tissue_type,
  sex = metadata$gender
)

specific_metadata$base_id <- sapply(strsplit(specific_metadata$barcode, "-"), function(x) paste(x[1:3], collapse = "-"))

final_meta <- specific_metadata



### Wrangle count data --
mrna_df <- assay(tcga_RNA_data)
colnames(mrna_df) <- sapply(strsplit(colnames(mrna_df), "-"), function(x) paste(x[1:4], collapse = "-"))
colnames(mrna_df)

mrna_df <- mrna_df[,colnames(mrna_df) %in% final_meta$bcr_patient_barcode] %>%
  as.data.frame()
dim(mrna_df)

filtered_mirna_df <- tcga_miRNA_data %>% 
  column_to_rownames("miRNA_ID") %>%
  dplyr::select(starts_with("read_count"))

colnames(filtered_mirna_df) <- gsub("read_count_","",colnames(filtered_mirna_df))
colnames(filtered_mirna_df) <- sapply(strsplit(colnames(filtered_mirna_df), "-"), function(x) paste(x[1:4], collapse = "-"))
filtered_mirna_df <- filtered_mirna_df[,colnames(filtered_mirna_df) %in% final_meta$bcr_patient_barcode]

# one last filter to match mrna and mirna samples 
names <- colnames(filtered_mirna_df)
dim(mrna_df)

mrna_unique <- mrna_df[,colnames(mrna_df) %in% names]

filtered_mirna_df <- filtered_mirna_df[,colnames(filtered_mirna_df) %in% colnames(mrna_unique)]
filtered_mirna_df <- filtered_mirna_df[,colnames(mrna_unique)]

# check for equality 
all(colnames(filtered_mirna_df) == colnames(mrna_unique))


# ensure metadata is also in the same order 
final_meta <- final_meta %>% filter(bcr_patient_barcode %in% colnames(filtered_mirna_df))
final_meta$bcr_patient_barcode <- factor(final_meta$bcr_patient_barcode, levels =colnames(filtered_mirna_df))


# Reorder the data frame based on the ordered SampleID
metadata_ordered <- final_meta[order(final_meta$bcr_patient_barcode), ]
metadata_ordered$SampleID <- gsub("-",".",metadata_ordered$bcr_patient_barcode)
metadata_ordered <- metadata_ordered %>% dplyr::select(SampleID, everything())
colnames(filtered_mirna_df) <- gsub("-", ".", colnames(filtered_mirna_df))
colnames(mrna_unique) <- gsub("-",".", colnames(mrna_unique))

# Add sequencing depth as a variable 
metadata_ordered$TCGA <- "TCGA-READ"
metadata_ordered <- metadata_ordered %>% mutate(Sequencing_Depth = colSums(filtered_mirna_df))

# Save wrangled count data 
write.csv(mrna_unique, "Tumor_vs_Normal_TCGA_READ_MRNA_counts.csv")
write.csv(filtered_mirna_df, "Tumor_vs_Normal_TCGA_READ_miRNA_counts.csv")
write.csv(metadata_ordered, "Tumor_vs_Normal_TCGA_READ_metadata.csv",row.names=FALSE)

################################## Combining TCGA COAD and TCGA READ to create COREAD #################################################

library(caret)
library(scutr)
library(tidyverse)
library(miRBaseConverter)

mirna_coad <- read.csv("Tumor_vs_Normal_TCGA_COAD_miRNA_counts.csv",row.names = 1)
mirna_read <- read.csv("Tumor_vs_Normal_TCGA_READ_miRNA_counts.csv", row.names =1)

mirna_coread <- cbind(mirna_coad,mirna_read)
write.csv(mirna_coread, "COREAD_Tumor_vs_Normal_miRNA.csv")

coad_meta <- read.csv("Tumor_vs_Normal_TCGA_COAD_metadata.csv")
read_meta <- read.csv("Tumor_vs_Normal_TCGA_READ_metadata.csv")

meta_coread <- dplyr::bind_rows(coad_meta,read_meta) 
write.csv(meta_coread,"COREAD_Tumor_vs_Normal_Metadata.csv",row.names=FALSE)

## Extract dimension 1 coordinates - downloaded from miRQuest app
mds <- read.csv("mds_xy_2025-05-08.csv")
names(mds)
mds <- mds %>% select(c("SampleID", "Dim1", "Dim2", "tissue_type")) %>% 
  column_to_rownames("SampleID")
set.seed(123)
undersamp <- undersample_tomek(mds, "Tumor", "tissue_type", 11, tomek = "diff", force_m = TRUE)


## Grab tomek samples --
normal <- meta_coread %>% filter(tissue_type=="Normal") %>% pull(SampleID)

samples <- c(normal,row.names(undersamp))

## Pull downsampled samples from mirna_coread_counts --
downsampled_coread_mirna <- mirna_coread[,names(mirna_coread) %in% samples]

## Convert miRBase ID --
version=checkMiRNAVersion(row.names(downsampled_coread_mirna), verbose = TRUE)
converted_miRNA_names <- miRNAVersionConvert(row.names(downsampled_coread_mirna),targetVersion = "v22",exact = TRUE)
all(converted_miRNA_names$OriginalName==converted_miRNA_names$TargetName)

mature_miRNA_names <- miRNA_PrecursorToMature(row.names(downsampled_coread_mirna))
downsampled_coread_mirna <- downsampled_coread_mirna %>% rownames_to_column("OriginalName") 

downsampled_coread_mirna <- left_join(downsampled_coread_mirna, mature_miRNA_names, by="OriginalName")

exclude_columns <- c("Mature1","Mature2", "OriginalName")
downsampled_coread_mirna_sum <- downsampled_coread_mirna %>% group_by(Mature1) %>%
  summarise(across(setdiff(names(downsampled_coread_mirna), exclude_columns), sum, na.rm = TRUE)) %>% 
  column_to_rownames("Mature1")
write.csv(downsampled_coread_mirna_sum, "COREAD_Downsampled_Tumor_vs_Normal_Mature_miRNA.csv")

## Match mRNA Samples to miRNA and write out --
mrna_coad <- read.csv("Tumor_vs_Normal_TCGA_COAD_MRNA_counts.csv", row.names = 1, check.names = FALSE)
mrna_read <- read.csv("Tumor_vs_Normal_TCGA_READ_MRNA_counts.csv", row.names = 1, check.names = FALSE)

mirna_coread <- read.csv("COREAD_Downsampled_Tumor_vs_Normal_Mature_miRNA.csv", row.names = 1, check.names = FALSE)

# Keep only mRNA samples that are present in miRNA COREAD
mrna_read <- mrna_read[, intersect(colnames(mrna_read), colnames(mirna_coread)), drop = FALSE]
mrna_coad <- mrna_coad[, intersect(colnames(mrna_coad), colnames(mirna_coread)), drop = FALSE]

# Pool genes safely
common_genes <- intersect(rownames(mrna_coad), rownames(mrna_read))
mrna_coread <- cbind(mrna_coad[common_genes, , drop = FALSE],
                     mrna_read[common_genes, , drop = FALSE])

write.csv(mrna_coread, "COREAD_Downsampled_Tumor_vs_Normal_mRNA.csv")