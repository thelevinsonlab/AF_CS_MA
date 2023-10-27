library(limma)
library(edgeR)
library(stringr)
library(tidyverse)
library(biomaRt)

Herrera_rivero_primeros_12 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Herrera-Rivero/Herrera_rivero_primeros_12.tabular")
Herrera_rivero_segundos_12 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Herrera-Rivero/Herrera_rivero_segundos_12.tabular")

df <- merge(Herrera_rivero_primeros_12, Herrera_rivero_segundos_12, by = "Geneid")
colnames(df)

rownames(df) <- df[,1]
df[,1] <- NULL

metadata <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/data_processing/scripts/Herrera-Rivero/metadata.txt")

#Targets file creation with only LA samples
targets <- metadata %>% filter(Site == "RA") %>% dplyr::select(c(sample, condition))

save(targets, file = "./Herrera_Rivero_targets.ro")

# Get the column names from df that are in the 'sample' column of 'targets'
cols_to_keep <- intersect(names(df), targets$sample)

# Subset 'df' to only keep the matching columns
df <- df[, cols_to_keep]

cols_order <- c(targets$sample)

# Reorder the columns
df <- df %>% dplyr::select(all_of(cols_order))

#Checking columns order
colnames(df) == targets$sample

#create DGE class object
group <- targets$condition
dge <- DGEList(counts=df, group=group)

#detect and remove low expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Apply normalization method TMM (trimmed mean of M)
dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)

Herrera_Rivero_count= v$E

Herrera_Rivero_count <- as.data.frame(Herrera_Rivero_count)

Herrera_Rivero_count <- Herrera_Rivero_count %>% rownames_to_column("Gene")

# Define the BioMart object
mart <- useMart(biomart = "ensembl", 
                dataset = "hsapiens_gene_ensembl", 
                host = 'useast.ensembl.org')

# Query BioMart for the gene symbols corresponding to your Ensembl IDs
# Here, 'Zhu_counts$Gene' refers to the column of your data frame that contains the Ensembl IDs
# You may need to replace 'Zhu_counts$Gene' with the actual column name in your data frame
IDs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
             filters = "ensembl_gene_id",
             values = Herrera_Rivero_count$Gene,
             mart = mart)

Herrera_Rivero_count$GeneSymbol <- IDs$hgnc_symbol[match(Herrera_Rivero_count$Gene, IDs$ensembl_gene_id)]

# Remove rows with NA or empty in 'GeneSymbol'
Herrera_Rivero_count <- Herrera_Rivero_count[!is.na(Herrera_Rivero_count$GeneSymbol), ]
Herrera_Rivero_count <- Herrera_Rivero_count[Herrera_Rivero_count$GeneSymbol != "", ]

# Create a table of duplicated gene symbols
dupes <- table(Herrera_Rivero_count$GeneSymbol)
dupes <- dupes[dupes > 1]

# Loop through each duplicated gene symbol
for (gene in names(dupes)) {
  # Get the indices of the rows with this gene symbol
  indices <- which(Herrera_Rivero_count$GeneSymbol == gene)
  
  # Append a unique identifier to each duplicate
  Herrera_Rivero_count$GeneSymbol[indices] <- paste(Herrera_Rivero_count$GeneSymbol[indices], "_", seq_along(indices), sep = "")
}

# Now set the rownames and remove the columns as before
rownames(Herrera_Rivero_count) <- Herrera_Rivero_count$GeneSymbol
Herrera_Rivero_count <- subset(Herrera_Rivero_count, select = -c(Gene, GeneSymbol))

save(Herrera_Rivero_count, file = "./Herrera_Rivero_count.ro")
