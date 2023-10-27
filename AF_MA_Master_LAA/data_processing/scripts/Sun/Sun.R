library(limma)
library(edgeR)
library(stringr)
library(tidyverse)
library(biomaRt)

df <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Sun/sun_counts.tabular")


colnames(df)

rownames(df) <- df[,1]
df[,1] <- NULL

metadata <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/data_processing/scripts/Sun/metadata.txt")
targets <- dplyr::select(metadata, sample, condition)

save(targets, file = "./Sun_targets.ro")

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

Sun_count= v$E

Sun_count <- as.data.frame(Sun_count)

Sun_count <- Sun_count %>% rownames_to_column("Gene")

# Define the BioMart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Query BioMart for the gene symbols corresponding to your Ensembl IDs
# Here, 'Sun_counts$Gene' refers to the column of your data frame that contains the Ensembl IDs
# You may need to replace 'Sun_counts$Gene' with the actual column name in your data frame
IDs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
             filters = "ensembl_gene_id",
             values = Sun_count$Gene,
             mart = mart)

Sun_count$GeneSymbol <- IDs$hgnc_symbol[match(Sun_count$Gene, IDs$ensembl_gene_id)]

# Remove rows with NA or empty in 'GeneSymbol'
Sun_count <- Sun_count[!is.na(Sun_count$GeneSymbol), ]
Sun_count <- Sun_count[Sun_count$GeneSymbol != "", ]

# Create a table of duplicated gene symbols
dupes <- table(Sun_count$GeneSymbol)
dupes <- dupes[dupes > 1]

# Loop through each duplicated gene symbol
for (gene in names(dupes)) {
  # Get the indices of the rows with this gene symbol
  indices <- which(Sun_count$GeneSymbol == gene)
  
  # Append a unique identifier to each duplicate
  Sun_count$GeneSymbol[indices] <- paste(Sun_count$GeneSymbol[indices], "_", seq_along(indices), sep = "")
}

# Now set the rownames and remove the columns as before
rownames(Sun_count) <- Sun_count$GeneSymbol
Sun_count <- subset(Sun_count, select = -c(Gene, GeneSymbol))

save(Sun_count, file = "./Sun_count.ro")
