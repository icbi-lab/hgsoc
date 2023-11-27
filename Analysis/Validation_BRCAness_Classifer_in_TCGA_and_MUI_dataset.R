
Here is the R script with detailed comments explaining each part of the code:

R
Copy code
# Load the required library for Gene Set Variation Analysis
library(GSVA)

# Read the TCGA data matrix from a CSV file, setting the first column as row names
data_matrix_TCGA = read.csv('TCGA_DataMatrix.csv', row.names = 1)
# Read the TCGA expression data from a CSV file, setting the first column as row names
exp_data_TCGA = read.csv('MUI_log2TPMs.csv', row.names = 1)

# Define lists of genes associated with HRD (homologous recombination deficiency)
# HRD positive gene list
HRD_plus = list(c('A2M','MRPL46','ABCD1','MRPS11','ADCK2','MT1E','ANKRD49','MTM1','ANXA1','MYL5','ASB7','NAB1','ASNSD1','NADSYN1','ATAD2','NARS2','ATP6V1E1','NDRG1','B3GALT4','NPEPL1',
                  'BBOX1','NQO2','BCAP31','NSDHL','BCL2L13','NUBPL','BID','OGFRL1','BOP1','P2RY6','BRCC3','PCGF3','C11orf71','PDE6B','C8orf33','PHF20L1','CCDC90B','PHKA1','CD83','PIR','CEP57','PNPLA4','CHN1','POLD4',
                  'CHORDC1','POLE3','COQ10B','PRCP','CPSF1','PRPS2','CRABP2','PRSS23','CREB3','PTK2','CRTC3','PUF60','CXCL10','PYCR1','CXCL11','RAB40B','CYC1','RAB9A','DCUN1D4','RABEPK','DCXR','RAP2C','DERL1',
                  'RECQL4','DGAT1','RMI1','DUSP22','RNF139','EHD1','S100A9','EIF1AX','SAC3D1','EML1','SAT1','ENDOD1','SCD','ERG','SCML1','ETFB','SDF2L1','EXOSC4','SECTM1','FAM3A','SFXN1','FANCI','SGPP1','FBXL6',
                  'SHARPIN','FJX1','SIRT5','FMR1','SLBP','FRYL','SLC10A3','FZD4','SLC25A1','GABRE','SLC39A4','GABRP','SLC3A2','GADD45B','SLC9A6','GBP1','SNRPA1','GBP2','SQLE','GCH1','ST6GALNAC5','GEMIN8',
                  'STEAP1','GPAA1','STK17B','GRINA','SYNGR2','GRPEL1','TANK','GSTK1','TFAP2C','HCCS','TINF2','HLA-F','TK1','HMGB3','TM2D3','HOPX','TMEM126B','HPRT1','TMEM135','HSF1','TMEM165','HSPA1A',
                  'TMEM187','IDH2','TMEM38B','IDH3A','TNFRSF21','IER3','TPM1','IL13RA1','TRMT12','IMPA2','UBB','INSIG1','UBL4A','JRKL','VBP1','KLK8','VPS28','KRT6A','YIF1A','LDB2','YIPF6','LRRC14','ZC3H3',
                  'LRRK1','ZDHHC24','MAPK9','ZNF16','MED17','ZNF250','MED7','ZNF267','MEF2A','ZNF34','METRN','ZNF623','MRPL13','ZNF7','MRPL40'))

# HRD negative gene list
HRD_neg = list(c('ALDH4A1','ALDOC','APP','ATXN2','BBS9','BMI1','BRCA1',
                 'C11orf49','CCDC130','CCDC92','CDK7','CDKN1C','CHRNB1','CHST12','CLASP1','CLIP3','CLSTN3','CPOX','DUOX1','ECSIT','EIF2AK1','EIF4G3','ERAL1','EXOSC10','FBXL18','FBXO2','FLOT2','GIPC1','GLTP','ICA1','LAMB2',
                 'LGALS8','LTA4H','LTBP4','MMS19','MRPS27','MYL6B','NAGPA','NAP1L1','NR2F6','PHACTR4','PHF21A','PKN1','PLTP','PMS2CL','POLR3B','POP4','PRDM4','PRKCSH','PTPN14','RAD17','RCAN3','RCBTB1','RCOR3',
                 'RPS6KC1','RRP15','SARS2','SMYD2','SMYD3','STX10','STX1A','SUPT6H','TEAD1','TP53BP2','TRAF5','TRIM8','TSC2','TTC31','USP48','UTP15','WIPI2','ZFYVE16','ZMYND8','ZNF10','ZNF12','ZNF84',''))

# Define gene sets for BRCAness negative and positive signatures
# BRCAness negative gene list
genesetBRCAness_neg1 = list(c('TSPAN1','TCEA3','FGGY','ZSWIM4','DCAF15','PTDSS1','DAPK3','NUDT15','CCAR2','RAD17','LTA4H'))
# BRCAness positive gene list
genesetBRCAness_pos1 = list(c('MT1M','PLGRKT','SNRPA1','DNAJC25','MRS2','USP28','FZD4','CCDC90B','PRCP','DERL3','NAPRT','GPAA1','CRABP2'))

# Perform GSVA for the BRCAness gene sets on the TCGA expression data using the ssGSEA method
gs_BRCAness_neg = gsva(t(exp_data_TCGA), genesetBRCAness_neg1, method = c("ssgsea"))
gs_BRCAness_pos = gsva(t(exp_data_TCGA), genesetBRCAness_pos1, method = c("ssgsea"))

# Convert the GSVA results to data frames and set the column names
Sig_brcaness_neg = as.data.frame(t(gs_BRCAness_neg))
colnames(Sig_brcaness_neg)[1] = 'Sig_brcaness_neg'
Sig_brcaness_pos = as.data.frame(t(gs_BRCAness_pos))
colnames(Sig_brcaness_pos)[1] = 'Sig_brcaness_pos'

# Combine the GSVA results into a single data frame and calculate the BRCAness Score
BRCAnessScore = cbind(Sig_brcaness_neg, Sig_brcaness_pos)
BRCAnessScore['BRCAnessScore'] = BRCAnessScore$Sig_brcaness_pos - BRCAnessScore$Sig_brcaness_neg

# Perform GSVA for the HRD gene sets on the TCGA expression data using the ssGSEA method
HRD_plus_ssgsea = gsva(t(exp_data_TCGA), HRD_plus, method = c("ssgsea"))
HRD_neg_ssgsea = gsva(t(exp_data_TCGA), HRD_neg, method = c("ssgsea"))

# Convert the GSVA results to data frames and set the column names
HRD_neg_ssgsea_res = as.data.frame(t(HRD_neg_ssgsea))
colnames(HRD_neg_ssgsea_res)[1] = 'HRD_neg_ssgsea'
HRD_plus_ssgsea_res = as.data.frame(t(HRD_plus_ssgsea))
colnames(HRD_plus_ssgsea_res)[1] = 'HRD_plus_ssgsea'

# Combine the HRD GSVA results into a single data frame and calculate the HRD signature score
HRD_signature_df = cbind(HRD_neg_ssgsea_res, HRD_plus_ssgsea_res)
HRD_signature_df['HRD_signature'] = HRD_signature_df$HRD_plus_ssgsea - HRD_signature_df$HRD_neg_ssgsea

# Merge the HRD signature results with the TCGA data matrix
data_matrix_TCGA_1 = merge(data_matrix_TCGA, HRD_signature_df, by = "row.names", all = TRUE)
# Merge the BRCAness Score with the TCGA data matrix
data_matrix_TCGA_1 = merge(data_matrix_TCGA_1, BRCAnessScore, by = "row.names", all = TRUE)


# Plot the relationship between HRD signature and BRCAness Score
ggplot(data_matrix_TCGA_1, aes(x = HRD_signature, y = BRCAnessScore)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  xlab("HRD signature") + 
  ylab("BRCAness signature") + 
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30))

# Plot the relationship between HRD signature and BRCAness probability
ggplot(data_matrix_TCGA_1, aes(x = HRD_signature, y = BRCAness_prob)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  xlab("HRD signature") + 
  ylab("BRCAness probability") + 
  theme(axis.text = element_text(size = 30), axis.title = element_text(size = 30))

# Perform Spearman correlation tests between HRD signature and BRCAness Score/Probability
cor.test(data_matrix_TCGA_1$HRD_signature, data_matrix_TCGA_1$BRCAnessScore, method = 'spearman')
cor.test(data_matrix_TCGA_1$HRD_signature, data_matrix_TCGA_1$BRCAness_prob, method = 'spearman')


# Read the MUI data matrix and expression data from CSV files
data_matrix_Mui = read.csv('/MUI_DataMatrix.csv',row.names = 1)
exp_data_Mui = read.csv('/MUI_log2_TPMs.csv',row.names = 1)

# (Repeat the steps above for the MUI dataset)
HRD_plus_ssgsea = gsva(as.matrix(exp_data_Mui),HRD_plus,method=c("ssgsea"))
HRD_neg_ssgsea = gsva(as.matrix(exp_data_Mui),HRD_neg,method=c("ssgsea"))

HRD_neg_ssgsea_res = as.data.frame(t(HRD_neg_ssgsea))
colnames(HRD_neg_ssgsea_res)[1] = 'HRD_neg_ssgsea'
HRD_plus_ssgsea_res = as.data.frame(t(HRD_plus_ssgsea))
colnames(HRD_plus_ssgsea_res)[1] = 'HRD_plus_ssgsea'
HRD_signature_df = cbind(HRD_neg_ssgsea_res,HRD_plus_ssgsea_res)

HRD_signature_df['HRD_signature'] = HRD_signature_df$HRD_plus_ssgsea-HRD_signature_df$HRD_neg_ssgsea
data_matrix_Mui_1 = merge(data_matrix_Mui,HRD_signature_df,all.x=TRUE,by.x = 0, by.y = 0)
rownames(data_matrix_Mui_1) = data_matrix_Mui_1$Row.names

gs_BRCAness_neg = gsva(as.matrix(exp_data_Mui),genesetBRCAness_neg1,method=c("ssgsea"))
gs_BRCAness_pos = gsva(as.matrix(exp_data_Mui),genesetBRCAness_pos1,method=c("ssgsea"))
Sig_brcaness_neg = as.data.frame(t(gs_BRCAness_neg))
colnames(Sig_brcaness_neg)[1] = 'Sig_brcaness_neg'
Sig_brcaness_pos = as.data.frame(t(gs_BRCAness_pos))
colnames(Sig_brcaness_pos)[1] = 'Sig_brcaness_pos'
BRCAnessScore = cbind(Sig_brcaness_neg,Sig_brcaness_pos)
BRCAnessScore['BRCAnessScore'] = BRCAnessScore$Sig_brcaness_pos-BRCAnessScore$Sig_brcaness_neg



data_matrix_Mui_1 = merge(data_matrix_Mui_1,BRCAnessScore,all.x=TRUE,by.x = 0, by.y = 0)
data_matrix_Mui_1 = data_matrix_Mui_1[, -1]



data_matrix_Mui_1$BRCAness_prob
data_matrix_Mui_1$HRD_signature

ggplot(data_matrix_Mui_1, aes(x = HRD_signature, y = BRCAness_prob)) +
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x)+
  xlab("HRD signature") + 
  ylab("BRCAness probability")+
  theme(axis.text = element_text(size = 30))+
  theme(axis.title=element_text(size=30))

ggplot(data_matrix_Mui_1, aes(x = HRD_signature, y = BRCAnessScore)) +
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x)+
  xlab("HRD signature") + 
  ylab("BRCAness signature")+ 
  theme(axis.text = element_text(size = 30))+
  theme(axis.title=element_text(size=30))

cor.test(data_matrix_Mui_1$HRD_signature,data_matrix_Mui_1$BRCAness_prob,method = 'spearman')
cor.test(data_matrix_Mui_1$HRD_signature,data_matrix_Mui_1$BRCAnessScore,method = 'spearman')


