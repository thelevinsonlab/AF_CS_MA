library(limma)
library(edgeR)
library(stringr)
library(tidyverse)
library(biomaRt)


df <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/data_processing/scripts/Chiang/Chiang.txt")

colnames(df)

rownames(df) <- df[,1]
df[,1] <- NULL

targets <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/data_processing/scripts/Chiang/metadata.txt")

save(targets, file = "./Chiang_targets.ro")

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

Chiang_count= v$E

Chiang_count <- as.data.frame(Chiang_count)

Chiang_count <- Chiang_count %>% rownames_to_column("Gene")

# Define the BioMart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Query BioMart for the gene symbols corresponding to your Ensembl IDs
# Here, 'Chiang_counts$Gene' refers to the column of your data frame that contains the Ensembl IDs
# You may need to replace 'Chiang_counts$Gene' with the actual column name in your data frame
IDs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
             filters = "ensembl_gene_id",
             values = Chiang_count$Gene,
             mart = mart)

Chiang_count$GeneSymbol <- IDs$hgnc_symbol[match(Chiang_count$Gene, IDs$ensembl_gene_id)]

# Remove rows with NA or empty in 'GeneSymbol'
Chiang_count <- Chiang_count[!is.na(Chiang_count$GeneSymbol), ]
Chiang_count <- Chiang_count[Chiang_count$GeneSymbol != "", ]

# Create a table of duplicated gene symbols
dupes <- table(Chiang_count$GeneSymbol)
dupes <- dupes[dupes > 1]

# Loop through each duplicated gene symbol
for (gene in names(dupes)) {
  # Get the indices of the rows with this gene symbol
  indices <- which(Chiang_count$GeneSymbol == gene)
  
  # Append a unique identifier to each duplicate
  Chiang_count$GeneSymbol[indices] <- paste(Chiang_count$GeneSymbol[indices], "_", seq_along(indices), sep = "")
}

# Now set the rownames and remove the columns as before
rownames(Chiang_count) <- Chiang_count$GeneSymbol
Chiang_count <- subset(Chiang_count, select = -c(Gene, GeneSymbol))

save(Chiang_count, file = "./Chiang_count.ro")