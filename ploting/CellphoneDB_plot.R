# Load the required libraries for heatmap and data manipulation
library(ComplexHeatmap)
library(dplyr)
library(circlize)

# Read in receptor data, order by interaction ID, and select relevant columns
data_cellPDB_receptor = read.csv('/path/to/ReceptorData.txt', sep = '\t', check.names = FALSE)
data_cellPDB_receptor = data_cellPDB_receptor[order(data_cellPDB_receptor$id_cp_interaction),]
receptor = data_cellPDB_receptor[,c(2,7:12,14:25)]
# Filter out rows where all selected columns have numeric values greater than 3
receptor = receptor %>% filter_all(any_vars(is.numeric(.) & . > 3))
receptor = receptor[order(receptor$id_cp_interaction),]
# Extract gene names from the receptor data
receptor_genenames = receptor$gene_name
# Finalize receptor data by selecting a range of columns
receptor_final = receptor[,3:19]
# Name the columns of the final receptor data
names(receptor_final) = c("B cells","CD169 macrophages", ... , "Regulatory T cells")

# Read in p-value data and order by interaction ID
data_cellPDB_p_value = read.csv('/path/to/PValueData.txt', sep = '\t')
data_cellPDB_p_value = data_cellPDB_p_value[order(data_cellPDB_p_value$id_cp_interaction),]
# Select columns and filter rows based on receptor interaction IDs
sig_df = data_cellPDB_p_value[,c(2,13:30)]
sig_df <- sig_df[sig_df$id_cp_interaction %in% receptor$id_cp_interaction, ]
sig_df = sig_df[order(sig_df$id_cp_interaction),]
# Finalize the p-value data by selecting and sorting columns
sig_df_final = sig_df[,c(2:6,8:19)]
sig_df_final = sig_df_final[,sort(colnames(sig_df_final))]
# Ensure that all values in the data frame are numeric
sig_df_final <- mutate_all(sig_df_final, function(x) as.numeric(as.character(x)))

# Read in ligand data, order by interaction ID, and filter for matching receptor interactions
data_cellPDB_ligand = read.csv('/path/to/LigandData.txt', sep = '\t')
ligand = data_cellPDB_ligand[,c(2,7,13)]
ligand = ligand[order(ligand$id_cp_interaction),]
ligand = ligand[ligand$id_cp_interaction %in% receptor$id_cp_interaction, ]
# Extract gene names from the ligand data
ligand_genenames = ligand$gene_name
# Select the final column for the ligand data
ligand_final = ligand[3]

# Define color functions for the heatmap
col_fun = colorRamp2(c(0, 8), c("white", "red"))
col_fun1 = colorRamp2(c(0, 8), c("white", "blue"))

# Set the row and column names for the significance data frame
rownames(sig_df_final) = rownames(receptor_final)
colnames(sig_df_final) = colnames(receptor_final)


# Create heatmaps for ligand and receptor expression
# Ligand heatmap
ht1 = Heatmap(t(ligand_final), width = unit(33, "cm"), height = unit(1, "cm"),column_names_gp = grid::gpar(fontsize = 20),
              row_names_gp = grid::gpar(fontsize = 20),cluster_columns =  FALSE , col = col_fun1, column_labels=ligand_genenames, column_names_side = "top", 
              heatmap_legend_param=list(title = 'logCPM', at = c(0,4,8)))

# Receptor heatmap with significance overlay
ht2 = Heatmap(t(receptor_final),width = unit(33, "cm"), height = unit(17, "cm"),column_names_gp = grid::gpar(fontsize = 20),
              row_names_gp = grid::gpar(fontsize = 20),cluster_columns =  FALSE , column_names_rot = 45, col = col_fun ,column_labels=receptor_genenames,
             heatmap_legend_param=list(title = 'logCPM', at = c(0,4,8)),
             cell_fun = function(j, i, x, y, w, h, fill) {
               # Overlay circles based on the p-value significance
             if(t(sig_df_final)[i, j] < 0.001) {
              r = min(unit.c(w, h))*0.20
              grid.circle(x, y, r, gp = gpar(col = "white"))
             }
             else if(t(sig_df_final)[i, j] < 0.01) {
              r = min(unit.c(w, h))*0.15
              grid.circle(x, y, r, gp = gpar(col = "white"))
             }
             else if(t(sig_df_final)[i, j] <= 0.1) {
              r = min(unit.c(w, h))*0.10
              grid.circle(x, y, r, gp = gpar(col = "white"))
             }
             else if(t(sig_df_final)[i, j] > 0.1) {
              grid.text("", x, y)
            }
          })
          
# Create a legend for the significance levels
lgd = Legend(at = 1:3, title = "FDR", labels = c("<0.001", "<0.01", '<0.1'), ...)

# Combine the two heatmaps vertically
htlist = ht1 %v% ht2

# Save the combined heatmap to a PDF file
pdf('/path/to/output/HeatmapCellPhoneDB.svg', width = 35, height = 20)
draw(htlist, heatmap_legend_list = lgd)
dev.off()
