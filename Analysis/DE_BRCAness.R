# Load the necessary libraries
library("DESeq2")

# Set the seed for reproducibility
set.seed(42)

# Read in count data from a tab-separated file
counts = read.csv("Counts_for_DEseq2.tab", 
                  sep = '\t', header = T, row.names = 1, check.names = FALSE)

# Read in metadata which contains the experimental design and sample information
meta = read.csv('Mata.csv', 
                header = TRUE, row.names = 1, check.names = FALSE)

# Prepare the metadata by matching the row names with the sample order
meta$Sample = rownames(meta)
meta <- meta[order(meta$Sample),]

# Ensure that the count data columns match the order of samples in the metadata
counts <- counts[,sort(names(counts))]

# Create a DESeq2 dataset object with the count data and the design formula
dds = DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~BRCAness)

# Relevel the factor to use "noBRCAness" as the reference level
dds$BRCAness <- relevel(dds$BRCAness, "noBRCAness")

# Run the DESeq pipeline for differential expression analysis
dds = DESeq(dds)

# Extract results and convert them to a DataFrame for easier handling
res1 = results(dds)
res_df1 = as.data.frame(res1)


# Write the results of the differential expression analysis to a CSV file
write.csv(res_df1,'DE_BRCAness.csv')
