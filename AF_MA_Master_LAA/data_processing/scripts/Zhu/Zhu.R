library(limma)
library(edgeR)
library(stringr)
library(tidyverse)
library(biomaRt)


Zhu_primeros_10 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Zhu_X_2020/Zhu_primeros_10.tabular")
Zhu_segundos_10 <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/New_datasets_alignment/Zhu_X_2020/Zhu_segundos_10.tabular")

df <- merge(Zhu_primeros_10, Zhu_segundos_10, by = "Geneid")
colnames(df)

rownames(df) <- df[,1]
df[,1] <- NULL

metadata <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/data_processing/scripts/Zhu/metadata.txt")

#Targets file refining (it already includes only LA samples)
targets <- metadata %>% dplyr::select(c(sample, condition))

save(targets, file = "./Zhu_targets.ro")

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

Zhu_count= v$E

Zhu_count <- as.data.frame(Zhu_count)

Zhu_count <- Zhu_count %>% rownames_to_column("Gene")

# Define the BioMart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Query BioMart for the gene symbols corresponding to your Ensembl IDs
# Here, 'Zhu_counts$Gene' refers to the column of your data frame that contains the Ensembl IDs
# You may need to replace 'Zhu_counts$Gene' with the actual column name in your data frame
IDs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
             filters = "ensembl_gene_id",
             values = Zhu_count$Gene,
             mart = mart)

Zhu_count$GeneSymbol <- IDs$hgnc_symbol[match(Zhu_count$Gene, IDs$ensembl_gene_id)]

# Remove rows with NA or empty in 'GeneSymbol'
Zhu_count <- Zhu_count[!is.na(Zhu_count$GeneSymbol), ]
Zhu_count <- Zhu_count[Zhu_count$GeneSymbol != "", ]

# Create a table of duplicated gene symbols
dupes <- table(Zhu_count$GeneSymbol)
dupes <- dupes[dupes > 1]

# Loop through each duplicated gene symbol
for (gene in names(dupes)) {
  # Get the indices of the rows with this gene symbol
  indices <- which(Zhu_count$GeneSymbol == gene)
  
  # Append a unique identifier to each duplicate
  Zhu_count$GeneSymbol[indices] <- paste(Zhu_count$GeneSymbol[indices], "_", seq_along(indices), sep = "")
}

# Now set the rownames and remove the columns as before
rownames(Zhu_count) <- Zhu_count$GeneSymbol
Zhu_count <- subset(Zhu_count, select = -c(Gene, GeneSymbol))

save(Zhu_count, file = "./Zhu_count.ro")
