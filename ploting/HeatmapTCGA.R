# Load necessary libraries for creating complex heatmaps and circular visualizations
library(ComplexHeatmap)
library(circlize)

# Read various data sets and annotations from CSV files into data frames, with the first column as row names
data_deconv <- read.csv('DataHeatmapImm/deconvHeatmap.csv',row.names = 1)
annotations <-  read.csv('DataHeatmapImm/annotHeatmap.csv',row.names = 1)
data_antigenpres_proc <-  read.csv('DataHeatmapImm/antigenpres_proc.csv',row.names = 1)
data_tgfb_caf <-  read.csv('DataHeatmapImm/tgfb.csv',row.names = 1)
data_checkpoints <- read.csv('DataHeatmapImm/checkpoints.csv',row.names = 1)
data_inf <- read.csv('DataHeatmapImm/inf.csv',row.names = 1)
data_cytoxEFF <- read.csv('DataHeatmapImm/cytoxeff.csv',row.names = 1)
data_T_exhaustion <- read.csv('DataHeatmapImm/T_exhaustion.csv',row.names = 1)

# Extract and convert annotations to character vectors
BRCAness <- as.character(as.vector(annotations['BRCANESS',]))
TumorImmFeno <- as.character(as.vector(annotations['TUMOUR_IMMUNE_PHENOTYPE',]))
MolSubtype <- as.character(as.vector(annotations['MOLECULAR_SUBTYPES',]))
HRR_MUT <- as.character(as.vector(annotations['HRR_MUT',]))
BRCA_MUT <- as.character(as.vector(annotations['BRCA_MUT',]))
BRCA_Methyl <- as.character(as.vector(annotations['BRCA_METHYL',]))

# Configure heatmap annotations with colors and labels
ha <- HeatmapAnnotation(BRCAness = BRCAness,
                        Tumor_immune_Phenotype= TumorImmFeno,
                        Molecular_Subtype = MolSubtype,
                        BRCA_mutation = BRCA_MUT,
                        col = list(BRCAness = c('BRCAness' = 'gray10','noBRCAness'='gray60'),
                                   Tumor_immune_Phenotype = c('Desert'='blue','Excluded'='darkorange','Infiltrated'='red3','unclassified'='gray50'),
                                   Molecular_Subtype = c('DIF'='red','IMR'='darkgreen','MES'='blue','PRO'='darkmagenta'),
                                   BRCA_mutation = c('no Mutation'='gray60','Mutation'='gray10')
                                   ))

# Define empty annotations for additional heatmap rows
ha1 = rowAnnotation(data_cytoxEFF = anno_empty(border = TRUE),width = unit(0.5, "cm"))
ha2 = rowAnnotation(data_imresp = anno_empty(border = TRUE),width = unit(0.5, "cm"))
ha3 = rowAnnotation(data_antigenpres_proc = anno_empty(border = TRUE),width = unit(0.5, "cm"))
ha4 = rowAnnotation(data_checkpoints = anno_empty(border = TRUE),width = unit(0.5, "cm"))
ha5 = rowAnnotation(data_inf = anno_empty(border = TRUE),width = unit(0.5, "cm"))
ha6 = rowAnnotation(data_T_exhaustion = anno_empty(border = TRUE),width = unit(0.5, "cm"))
ha7 = rowAnnotation(data_tgfb_caf = anno_empty(border = TRUE),width = unit(0.5, "cm"))

# Define the color gradient for quantitative immune infiltration
col_qinf = colorRamp2(c(0, 0.15), c("white", "tan3"))

# Create individual heatmaps for each data set without clustering and with defined styles
h1 = Heatmap(data_cytoxEFF,cluster_rows = FALSE,cluster_columns = FALSE, rect_gp = gpar(col = "grey",  lwd = 0.2),row_names_gp = grid::gpar(fontsize = 7),column_names_gp = grid::gpar(fontsize = 10),left_annotation = ha1)
h2 = Heatmap(data_imresp,cluster_rows = FALSE,cluster_columns = FALSE, rect_gp = gpar(col = "grey",  lwd = 0.2),row_names_gp = grid::gpar(fontsize = 7),column_names_gp = grid::gpar(fontsize = 10),left_annotation = ha2)
h3 = Heatmap(data_antigenpres_proc,cluster_rows = FALSE,cluster_columns = FALSE, rect_gp = gpar(col = "grey",  lwd = 0.2),row_names_gp = grid::gpar(fontsize = 7),column_names_gp = grid::gpar(fontsize = 10),left_annotation = ha3)
h4 = Heatmap(data_checkpoints,cluster_rows = FALSE,cluster_columns = FALSE, rect_gp = gpar(col = "grey",  lwd = 0.2),row_names_gp = grid::gpar(fontsize = 7),column_names_gp = grid::gpar(fontsize = 10),left_annotation = ha4)
h5 = Heatmap(data_inf,cluster_rows = FALSE,cluster_columns = FALSE, rect_gp = gpar(col = "grey",  lwd = 0.2),row_names_gp = grid::gpar(fontsize = 7),column_names_gp = grid::gpar(fontsize = 10),left_annotation = ha5)
h6 = Heatmap(data_T_exhaustion,cluster_rows = FALSE,cluster_columns = FALSE, rect_gp = gpar(col = "grey",  lwd = 0.2),row_names_gp = grid::gpar(fontsize = 7),column_names_gp = grid::gpar(fontsize = 10),left_annotation = ha6)
h7 = Heatmap(data_tgfb_caf,cluster_rows = FALSE,cluster_columns = FALSE, rect_gp = gpar(col = "grey",  lwd = 0.2),row_names_gp = grid::gpar(fontsize = 7),column_names_gp = grid::gpar(fontsize = 10),left_annotation = ha7)
h8=Heatmap(data_deconv,cluster_rows = FALSE,cluster_columns = FALSE, rect_gp = gpar(col = "grey",  lwd = 0.2),row_names_gp = grid::gpar(fontsize = 7),column_names_gp = grid::gpar(fontsize = 10),left_annotation = ha12,col = col_qinf)


# Combine all heatmaps and annotations vertically
ht_list = ha %v% h1 %v% h2 %v% h3 %v% h4 %v% h5 %v% h6 %v% h7 %v% h8

# Save the combined heatmap to an SVG file
svg(filename = 'Heatmap/HeatmapTCGA.svg', width = 21, height = 29.7)
draw(ht_list)
dev.off()

# Save the combined heatmap to a PNG file
png(filename = 'Heatmap/HeatmapTCGA.png', width = 2100, height = 2970)
draw(ht_list)
dev.off()
