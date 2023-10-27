library(limma)
library(edgeR)
library(stringr)
library(tidyverse)
library(biomaRt)


HS13 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Santos_2020/HS13.tabular")
HS1415 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Santos_2020/HS14-15.tabular")
HS1617 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Santos_2020/HS16-17.tabular")
HS1819 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Santos_2020/HS18-19.txt")
HS2021 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Santos_2020/HS20-21.tabular")
HS2223 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Santos_2020/HS22-23.tabular")


df <- merge(HS13, HS1415, by = "Geneid")
df <- merge(df, HS1617, by = "Geneid")
df <- merge(df, HS1819, by = "Geneid")
df <- merge(df, HS2021, by = "Geneid")
df <- merge(df, HS2223, by = "Geneid")
colnames(df)

rownames(df) <- df[,1]
df[,1] <- NULL

targets <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/data_processing/scripts/Santos/metadata.txt")

save(targets, file = "./Santos_targets.ro")

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

Santos_count= v$E

Santos_count <- as.data.frame(Santos_count)

Santos_count <- Santos_count %>% rownames_to_column("Gene")

# Define the BioMart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Query BioMart for the gene symbols corresponding to your Ensembl IDs
# Here, 'Santos_counts$Gene' refers to the column of your data frame that contains the Ensembl IDs
# You may need to replace 'Santos_counts$Gene' with the actual column name in your data frame
IDs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
             filters = "ensembl_gene_id",
             values = Santos_count$Gene,
             mart = mart)

Santos_count$GeneSymbol <- IDs$hgnc_symbol[match(Santos_count$Gene, IDs$ensembl_gene_id)]

# Remove rows with NA or empty in 'GeneSymbol'
Santos_count <- Santos_count[!is.na(Santos_count$GeneSymbol), ]
Santos_count <- Santos_count[Santos_count$GeneSymbol != "", ]

# Create a table of duplicated gene symbols
dupes <- table(Santos_count$GeneSymbol)
dupes <- dupes[dupes > 1]

# Loop through each duplicated gene symbol
for (gene in names(dupes)) {
  # Get the indices of the rows with this gene symbol
  indices <- which(Santos_count$GeneSymbol == gene)
  
  # Append a unique identifier to each duplicate
  Santos_count$GeneSymbol[indices] <- paste(Santos_count$GeneSymbol[indices], "_", seq_along(indices), sep = "")
}

# Now set the rownames and remove the columns as before
rownames(Santos_count) <- Santos_count$GeneSymbol
Santos_count <- subset(Santos_count, select = -c(Gene, GeneSymbol))

save(Santos_count, file = "./Santos_count.ro")
