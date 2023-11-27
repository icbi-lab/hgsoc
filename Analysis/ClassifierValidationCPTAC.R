# Load required libraries
library(ggplot2)
library(dplyr)

# Load combined data from a CSV file
data = read.csv('CombinedData_matrix.csv')

# Load expression data from a text file
expression_cptac = read.csv('CPTAC_OV_expression.txt', sep = '\t', check.names = FALSE)
# Remove duplicate entries based on the 'GeneID' column
expression_cptac1 <- expression_cptac[!duplicated(expression_cptac$GeneID),]
# Set row names to be the gene IDs
rownames(expression_cptac1) = expression_cptac1$GeneID
# Remove the GeneID column and transpose the data
expression_cptac1 = expression_cptac1[,-1]
expression_cptac1 = as.data.frame(t(expression_cptac1))

# Define gene sets for HRD positive and negative signatures
HRD_plus = list(c('A2M', 'MRPL46', ... 'ZNF7'))
HRD_neg = list(c('ALDH4A1', 'ALDOC', ... 'ZNF84'))

# Gene sets for BRCAness positive and negative
genesetBRCAness_neg1 = list(c('TSPAN1', 'TCEA3', ... 'LTA4H'))
genesetBRCAness_pos1 = list(c('MT1M', 'PLGRKT', ... 'TMEM38B'))

# Perform ssGSEA for HRD positive and negative signatures
HRD_plus_ssgsea = gsva(t(expression_cptac1), HRD_plus, method=c("ssgsea"))
HRD_neg_ssgsea = gsva(t(expression_cptac1), HRD_neg, method=c("ssgsea"))

# Convert the results to data frames and rename the columns
HRD_neg_ssgsea_res = as.data.frame(t(HRD_neg_ssgsea))
colnames(HRD_neg_ssgsea_res)[1] = 'HRD_neg_ssgsea'
HRD_plus_ssgsea_res = as.data.frame(t(HRD_plus_ssgsea))
colnames(HRD_plus_ssgsea_res)[1] = 'HRD_plus_ssgsea'
# Combine both positive and negative HRD ssGSEA results into one data frame
HRD_signature_df = cbind(HRD_neg_ssgsea_res, HRD_plus_ssgsea_res)
# Create a combined HRD signature score by subtracting negative from positive
HRD_signature_df['HRD_signature'] = HRD_signature_df$HRD_plus_ssgsea - HRD_signature_df$HRD_neg_ssgsea

# Perform ssGSEA for BRCAness positive and negative gene sets
gs_BRCAness_neg = gsva(t(expression_cptac1), genesetBRCAness_neg1, method=c("ssgsea"))
gs_BRCAness_pos = gsva(t(expression_cptac1), genesetBRCAness_pos1, method=c("ssgsea"))
# Convert the results to data frames and rename the columns
Sig_brcaness_neg = as.data.frame(t(gs_BRCAness_neg))
colnames(Sig_brcaness_neg)[1] = 'Sig_brcaness_neg'
Sig_brcaness_pos = as.data.frame(t(gs_BRCAness_pos))
colnames(Sig_brcaness_pos)[1] = 'Sig_brcaness_pos'
# Combine both positive and negative BRCAness ssGSEA results into one data frame
BRCAnessScore = cbind(Sig_brcaness_neg, Sig_brcaness_pos)
# Create a combined BRCAness score by subtracting negative from positive
BRCAnessScore['BRCAnessScore'] = BRCAnessScore$Sig_brcaness_pos - BRCAnessScore$Sig_brcaness_neg

# Merge the combined data with BRCAness scores
data = merge(data, BRCAnessScore, all.x=TRUE, by.x = 8, by.y = 0)
# Merge the combined data with HRD signatures
data = merge(data, HRD_signature_df, all.x=TRUE, by.x = 1, by.y = 0)

# Plotting and correlation analysis between SigMA signature and BRCAness probability
ggplot(data, aes(x=Signature_3_c, y=BRCAness_prob)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  xlab("SigMA") + 
  ylab("BRCAness prob") + 
  theme(axis.title=element_text(size=30))
# Calculate Spearman correlation
cor(data$Signature_3_c, data$BRCAness_prob, method = 'spearman')

# Plotting and correlation analysis between SigMA signature and BRCAness score
ggplot(data, aes(x=Signature_3_c, y=BRCAnessScore)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  xlab("SigMA") + 
  ylab("BRCAness signature") +
  theme(axis.text = element_text(size = 30)) +
  theme(axis.title=element_text(size=30))
# Perform Spearman correlation test
cor.test(data$Signature_3_c, data$BRCAnessScore, method = 'spearman')

# Plotting and correlation analysis between SigMA signature and HRD signature
ggplot(data, aes(x=Signature_3_c, y=HRD_signature)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  xlab("SigMA") + 
  ylab("BRCAness prob") +
  theme(axis.title=element_text(size=30))
# Calculate Spearman correlation
cor(data$Signature_3_c, data$HRD_signature, method = 'spearman')
