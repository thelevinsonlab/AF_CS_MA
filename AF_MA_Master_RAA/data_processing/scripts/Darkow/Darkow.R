library(limma)
library(edgeR)
library(stringr)
library(tidyverse)
library(biomaRt)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library("DESeq2")


df <- read.delim("data_processing/scripts/Darkow/rm_darkow.tabular")

rownames(df) <- df[,1]
df[,1] <- NULL

#Targets creation

targets <- as.data.frame(matrix(NA,length(names(df)),2))
names(targets) <- c("sample","condition")
targets$sample <- names(df)
targets$condition = str_sub(targets$sample, start = 2, end = 3)
targets <- targets %>% mutate(condition = case_when(
  condition > 42 & condition < 75 ~ "AF",
  condition < 42 & condition > 20 ~ "SR"))

save(targets, file = "./Darkow_targets.ro")

# Get the column names from df that are in the 'sample' column of 'targets'
cols_to_keep <- intersect(names(df), Darkow_targets$sample)

# Subset 'df' to only keep the matching columns
df <- df[, cols_to_keep]

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

Darkow_count= v$E

Darkow_count <- as.data.frame(Darkow_count)

Darkow_count <- Darkow_count %>% rownames_to_column("Gene")

# Define the BioMart object
mart <- useMart(biomart = "ensembl", 
                dataset = "hsapiens_gene_ensembl", 
                host = 'useast.ensembl.org')

# Query BioMart for the gene symbols corresponding to your Ensembl IDs
# Here, 'Zhu_counts$Gene' refers to the column of your data frame that contains the Ensembl IDs
# You may need to replace 'Zhu_counts$Gene' with the actual column name in your data frame
IDs <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
             filters = "ensembl_gene_id",
             values = Darkow_count$Gene,
             mart = mart)

Darkow_count$GeneSymbol <- IDs$hgnc_symbol[match(Darkow_count$Gene, IDs$ensembl_gene_id)]

# Remove rows with NA or empty in 'GeneSymbol'
Darkow_count <- Darkow_count[!is.na(Darkow_count$GeneSymbol), ]
Darkow_count <- Darkow_count[Darkow_count$GeneSymbol != "", ]

# Create a table of duplicated gene symbols
dupes <- table(Darkow_count$GeneSymbol)
dupes <- dupes[dupes > 1]

# Loop through each duplicated gene symbol
for (gene in names(dupes)) {
  # Get the indices of the rows with this gene symbol
  indices <- which(Darkow_count$GeneSymbol == gene)
  
  # Append a unique identifier to each duplicate
  Darkow_count$GeneSymbol[indices] <- paste(Darkow_count$GeneSymbol[indices], "_", seq_along(indices), sep = "")
}

# Now set the rownames and remove the columns as before
rownames(Darkow_count) <- Darkow_count$GeneSymbol
Darkow_count <- subset(Darkow_count, select = -c(Gene, GeneSymbol))

Darkow_count <- as.matrix(Darkow_count)

save(Darkow_count, file = "./Darkow_count.ro")



































# Assuming 'af_types' is your variable
af_types <- c("Permanent", "Not reported", "Permanent", "Persistent", 
              "Paroxysmal", "Permanent", "Persistent", "Paroxysmal;Permanent", 
              "Paroxysmal;Persistent;Permanent", "Paroxysmal;Persistent", 
              "Not reported", "Persistent;Permanent", "Not reported", 
              "Paroxysmal;Persistent", "Persistent", "Permanent", "Not reported", 
              "Permanent", "Not reported", "Permanent", "Not reported", 
              "Persistent;Permanent", "Paroxysmal;Persistent", "Permanent", 
              "Persistent", "Paroxysmal;Persistent", "Persistent", "Not reported", 
              "Persistent", "Paroxysmal;Persistent;Permanent", "Persistent", 
              "Persistent", "Permanent", "Persistent")

# Split the 'af_types' by ';'
split_types <- strsplit(af_types, split = ";")

# Unlist the result to get a single vector
unlisted_types <- unlist(split_types)

# Use the table function to count the occurrences of each category
type_counts <- table(unlisted_types)

# Print the result
print(type_counts)

# Convert the table to a dataframe for easier manipulation
df <- as.data.frame(type_counts)
names(df) <- c("Type", "Count")

# Calculate proportions
df$Proportion <- df$Count / sum(df$Count)

# First, rearrange the data in descending order
df <- df[order(-df$Count), ]

# Next, set the levels of the factor variable Type to match the order in df
df$Type <- factor(df$Type, levels = df$Type)

# Create the bar plot
barplot <- ggplot(df, aes(x=Type, y=Count, fill=Type)) +
  geom_bar(stat="identity", width=0.7, colour="black") +
  geom_text(aes(label=Count), vjust=-0.3, size=3.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Type", y = "Count", title = "Barplot of AF Types") 

# Save the plot
ggsave("af_types_barplot.png", plot = barplot, width = 10, height = 8, dpi = 600)


# Create a data frame
df <- data.frame(
  Method = c("Microarrays", "RNA-Seq"),
  Percentage = c(44, 56)
)

# Define darker colors
colors <- c("darkblue", "darkgreen")

# Create a pie plot
ggplot(df, aes(x="", y=Percentage, fill=Method)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position="right") +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), color = "white") +
  labs(title = "Method Distribution")





tab_information <- read.delim2("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/tab_information.txt")

head(tab_information)

# Make sure you have tidyverse installed
# install.packages("tidyverse")

# Convert the data to long format
tab_information_long <- tab_information %>%
  pivot_longer(cols = -First.author.name, names_to = "Variable", values_to = "Type")

# Generate the graph
ggplot(tab_information_long, aes(x = First.author.name, y = Variable, fill = Type)) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Group" = "#DDE91A", "Individual" = "#5DC207", "NR" = "#CF3F17")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE))

tab_information <- read.delim2("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/tab_information.txt")

# Convert the data to long format
tab_information_long <- tab_information %>%
  pivot_longer(cols = -First.author.name, names_to = "Variable", values_to = "Type")

# Define the desired order
levels_order <- c("Age", "Sex", "Ethnicity", "Comorbidities", "Treatment", "LVEF")

# Convert 'Variable' to a factor with the specified level order
tab_information_long$Variable <- factor(tab_information_long$Variable, levels = levels_order)

# Generate the graph
final_plot <- ggplot(tab_information_long, aes(x = reorder(First.author.name, First.author.name), y = Variable, fill = Type)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Type == "Group", "Gr", ifelse(Type == "Individual", "Ind", "NR"))), 
            size = 4, color = "black") +
  scale_fill_manual(values = c("Group" = "#DDE91A", "Individual" = "#5DC207", "NR" = "#CF3F17")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE))

# Save the plot
ggsave("info_plot.png", plot = final_plot, width = 6, height = 8, dpi = 600)

sample_information <- read.delim("~/PhD Heidelberg/Projects/AF_MA_Master_LAA/sample_information.txt")
head(sample_information)

sample_information_long <- sample_information %>%
  gather(key = "Sample Type", value = "Count", 
         c("Total_AF_samples", "Total_SR_samples")) 

# Make a variable to order the bars
sample_information <- sample_information %>% 
  mutate(Order = row_number(Total_samples))

# Transform the data to long format
sample_information_long <- sample_information %>%
  gather(key = "Sample Type", value = "Count", 
         c("Total_AF_samples", "Total_SR_samples")) 

sample_information_long <- sample_information_long %>% mutate(`Sample Type`=case_when(`Sample Type` == "Total_AF_samples" ~ "AF",
                                                                                      `Sample Type` == "Total_SR_samples" ~ "SR"))


# Create the plot
sample_plot <- ggplot(sample_information_long, 
       aes(x = reorder(First.author.name, Order), y = Count, fill = `Sample Type`)) +
  geom_bar(stat = "identity", position = 'stack', color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Study", y = "Sample Count", fill = "Sample Type") +
  coord_flip() + scale_y_continuous(breaks = seq(0, 300, by = 20))

# Save the plot
ggsave("sample_plot.png", plot = sample_plot, width = 10, height = 8, dpi = 600)