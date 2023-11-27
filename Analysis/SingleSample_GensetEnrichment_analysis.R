# Load the GSVA library for gene set variation analysis
library(GSVA)

# Read in a CSV file containing gene lists for a method called avex
genelists_avex <- read.csv('ImmuneSignaturesShort.csv', check.names = FALSE)
# Get column names from the avex gene lists
cols_avex <- colnames(genelists_avex)

# Read in a CSV file containing gene lists for ssGSEA analysis
genelists_ssgsea <- read.csv('ImmuneSignatures_ssGSEA.csv', check.names = FALSE)
# Convert the gene lists for ssGSEA into a list format
genelists_ssgsea = as.list(genelists_ssgsea)

# Read in expression data from a CSV file, setting the first column as row names
data <- read.csv('log2TPM1.csv', row.names = 1, check.names = FALSE)

# Transpose the expression data to have genes as columns and samples as rows
expression <- t(data)

# Perform ssGSEA using the expression data and gene lists
SSGSEA <- gsva(expression,
               genelists_ssgsea,
               method = 'ssgsea',
               ssgsea.norm=FALSE) # no normalization in ssGSEA

# Get column names from the expression data to use as column names for the avex matrix
names = colnames(expression)
# Create an empty data frame to store avex results with the same column names as 'expression'
avex = data.frame(matrix(ncol = 321, nrow = 0))
colnames(avex) = names

# Loop over each gene list in avex
for (i in cols_avex){
  print(i) # Print the current gene list being processed
  # Extract the gene list from the genelists_avex dataframe
  genlist = genelists_avex[,i]
  # Remove empty entries from the gene list
  genlist = genlist[genlist != ""]
  # Subset the expression data to include only genes in the current gene list
  expression_genes = expression[genlist,]
  # Calculate the average expression for each gene list across all samples
  av_expression = colMeans(expression_genes)
  # Bind the average expression as a new row in the avex dataframe
  avex = rbind(avex, as.data.frame(t(av_expression)))
}
# Set the row names of the avex dataframe to the gene list names
rownames(avex) = cols_avex

# Convert the SSgsea results into a dataframe
ssGSEA = as.data.frame(SSGSEA)

# Combine the ssGSEA and avex results into a single dataframe
signatures = rbind(ssGSEA, avex)

# Write the combined signature results to a CSV file
write.csv(signatures, 'GSEA_avex_Signature.csv')
