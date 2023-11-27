# Load the necessary libraries
library(psych)         # For statistical operations like correlations
library(ComplexHeatmap) # For creating complex heatmaps for data visualization
library(gdata)         # For data manipulation, especially reading and writing files

# Read in the main data set, gene set enrichment analysis (GSEA) signatures, and immunophenogram scores
data <- read.csv('/path/to/MASTERFILE.csv', check.names = FALSE)
data_signatures = read.csv('/path/to/GSEA_avex_Signature.csv', row.names = 1, check.names = FALSE)
data_signatures = as.data.frame(t(data_signatures)) # Transpose the signatures data for compatibility
data_ips = read.csv('/path/to/IPS.txt', row.names = 1, check.names = FALSE, sep = '\t')

# Combine all datasets into a single data frame for analysis
data = cbind(data, data_signatures)
data = merge(data, data_ips, by = "row.names", all.x = TRUE)

# Define the columns that will be used for correlation analysis
columns = c('NeoAG_Load', 'TMB', ..., 'CD8 EXHAUSTION')

# Subset the data to include only the columns of interest
data_correlations = data[, columns]

data_correlations = data[,columns]

# Define the continuous variables that will be analyzed
continous = c('NeoAG_Load','TMB','HRD_Score','SIG3RATIO','BRCA1_METHYLATION_BETA_VALUE','BRCA2_METHYLATION_BETA_VALUE','CD8 Jiang','IFNG Ayers',
              'CTL Jiang','CYT Rooney','INFLAM_SprangerssGSEA','EXPAND_IMM_AyersssGSEA','T_EXCL_jerby_anonssGSEA','CD8 EXHAUSTION')

# Rename the variables for readability in the final output
final_names = c('NeoAG load','TMB','HRDscore','SIG3 ratio','BRCA1 meth','BRCA2 meth','CD8 sign','IFNG sign',
                'CTL sign','CYT sign','INFLAM sign','EXPAND IMM sign','T EXCLUSION sign','CD8 EXHAUSTION sign','BRCAness','HRR mut')

# Initialize data frames to store correlation coefficients and p-values
corrDF_pointbiserial = data.frame(matrix(ncol = 14, nrow = 0))
corrSig_pointbiserial = data.frame(matrix(ncol = 14, nrow = 0))

# Initialize vectors to hold the correlation results
p = c()
coef = c()


# Perform point biserial correlation for 'BRCANESS' against continuous variables
for (i in continuous) {
  x <- ifelse(data_correlations[,'BRCANESS'] == 'BRCAness', 1, 0)
  y <- data_correlations[,i]
  res <- cor.test(x, y)
  corrDF_pointbiserial <- rbind(corrDF_pointbiserial, res$estimate)
  corrSig_pointbiserial <- rbind(corrSig_pointbiserial, res$p.value)
}


# Continue the same process for 'HRR_MUT'
# The code performs similar operations for 'HRR_MUT' variable against continuous variables
p = c()
coef = c()

for (i in continous) {
  x= data_correlations[,'HRR_MUT']
  y = data_correlations[,i]
  res = cor.test(x, y)
  res
  p = c(p,res$p.value)
  coef = c(coef,res$estimate)
}

corrDF_pointbiserial = rbind(corrDF_pointbiserial,coef)
corrSig_pointbiserial = rbind(corrSig_pointbiserial,p)

rownames(corrDF_pointbiserial) = c('BRCANESS','HRR_MUT')
rownames(corrSig_pointbiserial) = c('BRCANESS','HRR_MUT')
colnames(corrDF_pointbiserial) = continous
colnames(corrSig_pointbiserial) = continous

# Prepare a table and perform a chi-squared test for 'HRR_MUT' and 'BRCANESS'
t <- table(hrr_mut, brcaness)
chisq_res <- chisq.test(hrr_mut, brcaness, correct=FALSE)
phi_coef <- phi(t)  # Calculate phi coefficient for the chi-squared test


# Calculate correlations between all continuous variables
# The code constructs a correlation matrix for continuous variables

data_frame_continous = data_correlations[,continous]

cont_corr = data.frame(nrow = 14)
cont_corr_p = data.frame(nrow = 14)

for (i in 1:length(data_frame_continous)){
  col_con = c()
  col_p = c()
  for (x in 1:length(data_frame_continous)){
    correlation = cor.test(data_frame_continous[,i],data_frame_continous[,x])
    coef = cor(data_frame_continous[,i],data_frame_continous[,x],use='complete.obs')
    p = correlation$p.value
    col_con = c(col_con,coef)
    col_p = c(col_p,p)
  }
  cont_corr = cbind(cont_corr,col_con)
  cont_corr_p = cbind(cont_corr_p,col_p)
}

cont_corr = cont_corr[-1]
cont_corr_p = cont_corr_p[-1]

rownames(cont_corr_p) = colnames(data_frame_continous)
colnames(cont_corr_p) = colnames(data_frame_continous)

rownames(cont_corr) = colnames(data_frame_continous)
colnames(cont_corr) = colnames(data_frame_continous)

# Add BRCANESS and HRR_MUT data to the correlation matrix
# The code creates final correlation and p-value matrices including BRCANESS and HRR_MUT

es_df = rbind(cont_corr,corrDF_pointbiserial)

brcaness_col = as.numeric(as.vector(corrDF_pointbiserial['BRCANESS',]))
hrr_mut_col = as.numeric(as.vector(corrDF_pointbiserial['HRR_MUT',]))

BRCANESS = c(brcaness_col,1,phi_coef)
HRR_MUT = c(hrr_mut_col,phi_coef,1)

res_df = cbind(res_df,BRCANESS)
res_df = cbind(res_df,HRR_MUT)



sig_df = rbind(cont_corr_p, corrSig_pointbiserial)

brcaness_col_sig = as.numeric(as.vector(corrSig_pointbiserial['BRCANESS',]))
hrr_mut_col_sig = as.numeric(as.vector(corrSig_pointbiserial['HRR_MUT',]))

BRCANESS = c(brcaness_col_sig,0,chisq_res$p.value)
HRR_MUT = c(hrr_mut_col_sig,chisq_res$p.value,0)


sig_df = cbind(sig_df,BRCANESS)
sig_df = cbind(sig_df,HRR_MUT)

p_vals = c() 


for (i in 1:length(sig_df)){
  s = as.numeric(as.vector(sig_df[,i]))
  p = first(s,n=i)
  p_vals = c(p_vals,p)
}

adj_pvals = p.adjust(p_vals,method = 'BH')

len_m = length(sig_df)+1

b= matrix(0, len_m, len_m)
b[upper.tri(b, diag=FALSE)] <- adj_pvals 
b = b[,-1]
b = b[-len_m,]
#b<-round(b,4)

f <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

b = f(b)


b

adj_p = as.data.frame(b)
colnames(adj_p) = final_names
rownames(adj_p) = final_names
colnames(res_df) = final_names
rownames(res_df) = final_names

# Create a heatmap of the correlation matrix with annotations for significance levels
ht = Heatmap(as.matrix(res_df),width = unit(10, "cm"), height = unit(10, "cm"),column_names_rot = 45,heatmap_legend_param=list(
  title = "coef", at = c(-1, 0, 1), 
  labels = c(-1, 0, 1)),cell_fun = function(j, i, x, y, w, h, fill) {
    # The code inside this function adds circles of varying size to the heatmap cells based on significance levels
      if(as.matrix(sig_df)[i, j] < 0.001) {
      r = min(unit.c(w, h))*0.20
      grid.circle(x, y, r, gp = gpar(col = "white"))
    }
    else if(as.matrix(sig_df)[i, j] < 0.01) {
      r = min(unit.c(w, h))*0.15
      grid.circle(x, y, r, gp = gpar(col = "white"))
    }
    else if(as.matrix(sig_df)[i, j] < 0.1) {
      r = min(unit.c(w, h))*0.10
      grid.circle(x, y, r, gp = gpar(col = "white"))
    }
    else if(as.matrix(sig_df)[i, j] > 0.1) {
      grid.text("", x, y)
    }
  })

# Add a legend to the heatmap
lgd <- Legend(at = 1:3, title = "FDR", labels = c("<0.001", "<0.01",'<0.1'), legend_gp = gpar(col = 1:3))

# Output the heatmap to an SVG file
svg('IMM_BRCAness_Correlation.svg')
draw(ht, heatmap_legend_list = lgd)
dev.off()  # Close the SVG device

