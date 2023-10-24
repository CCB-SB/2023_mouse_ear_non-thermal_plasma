suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))


set.seed(snakemake@params$parameters_props$set_seed)

print("heatmap_corr")


#---------------------------------- Functions ----------------------------------
plot_heatmap_annot = function(df, row_names, plots_props, prop, annotation_row, colors, filename) {
  
  th = 0.3
  df = t(df)
  df[abs(df) <= th] = 0
  df[df < -th] = -0.5
  df[df > th] = 0.5
  
  color_bar_name = "Correlation (Spearman)" 
  cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  
  col_val = 0.75
  col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  col_legend = c(cmap(-col_val), cmap(col_val))
  names(col_legend) = c("negative", "positive")
  
  annotation_colors = list()
  ucolors <- unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_row[[prop]] %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  colnames(df) = row_names
  column_ha = rowAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"), 
                            annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
                            show_legend=FALSE)
  
  p = ComplexHeatmap::Heatmap(df, col = col_fun,
                              height = unit(3, "cm"), width = unit(6.25, "cm"),
                              show_column_names = FALSE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family), 
                              cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              column_dend_reorder = TRUE, 
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              left_annotation = column_ha
  )

  legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
                                  direction="horizontal",
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(4, "cm"))
  
  png(sprintf("%s.png",filename), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(p, padding = unit(c(0.25, 0.25, -1.5, 0.25), "cm"))  # bottom, left, top and right margins
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s.svg",filename), width = plots_props$image_width / 2.54, height = plots_props$image_height / 2.54)
  draw(p, padding = unit(c(0.25, 0.25, -1.5, 0.25), "cm"))  # bottom, left, top and right margins
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
}

plot_heatmap_annot_filtered = function(df, row_names, plots_props, prop, annotation_row, miRNA_list, colors, filename) {
  
  th = 0.3
  df = t(df)
  df[abs(df) <= th] = 0
  df[df < -th] = -0.5
  df[df > th] = 0.5
  
  row_names = row_names[colSums(df != 0) > 0]
  tmp = !(row_names %in% miRNA_list$miRNA)
  row_names[tmp] = ""
  df_label = matrix("*", nrow(df), ncol(df))
  df_label[,tmp] = ""
  df = df[, colSums(df != 0) > 0]
  
  color_bar_name = "Correlation (Spearman)" 
  cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  
  col_val = 0.75
  col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  col_legend = c(cmap(-col_val), cmap(col_val))
  names(col_legend) = c("negative", "positive")
  
  annotation_colors = list()
  ucolors <- unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_row[[prop]] %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  colnames(df) = row_names
  column_ha = rowAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"), 
                            annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
                            show_legend=FALSE)
  cell_fun = function(j, i, x, y, width, height, fill) {grid.text(df_label[i, j], x, y, gp = gpar(fontsize = 6))}
  
  p = ComplexHeatmap::Heatmap(df, col = col_fun,
                              cell_fun = cell_fun,
                              height = unit(1.5, "cm"), width = unit(6.25, "cm"),
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family), 
                              column_names_side = "bottom", column_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              column_dend_reorder = TRUE, 
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              left_annotation = column_ha
  )

  legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
                                  direction="horizontal",
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(4, "cm"))
  
  png(sprintf("%s.png",filename), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(p, padding = unit(c(-0.25, 0.25, -1.5, 0.25), "cm"))  # bottom, left, top and right margins
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s.svg",filename), width = plots_props$image_width / 2.54, height = plots_props$image_height / 2.54)
  draw(p, padding = unit(c(-0.5, 0.25, -1.5, 0.25), "cm"))  
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
}


#---------------------------------- Load data ---------------------------------- 
expr = fread(snakemake@input$expr, sep='\t')
annot = fread(snakemake@input$annot, sep='\t')
miRNA_list = fread("data/miRNA_list.csv")

plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

output_folder = sprintf("%s/heatmap_complex/corr_plots", results_folder)
dir.create(output_folder, recursive=TRUE)
supp_table_output_folder = sprintf("%s/matrices/supp_corr", results_folder)
dir.create(supp_table_output_folder, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
# corr between miRNA and Timepoints for every Status
correlation_Status_Timepoints = c()
for(group in unique(annot$Status)){
  groups_ID = annot[annot$Status == group,]$ID
  tmp = colnames(expr) %in% groups_ID
  expr_sort = data.frame(miRNA=expr$miRNA, expr[,..tmp], check.names = FALSE)
  timepoints = annot[annot$ID %in% groups_ID,]$Time
  timepoints = as.numeric(gsub("m", "", timepoints))
  
  correlation_Status_Timepoints[[group]] = as.data.frame(cor(t(expr_sort[,2:dim(expr_sort)[2]]), timepoints, method ="spearman"))
}

correlation_Status_Timepoints = as.data.frame(correlation_Status_Timepoints)
colnames(correlation_Status_Timepoints) = unique(annot$Status)
annot_row = data.frame(Status=c("Control", "Treatment"))
rownames(annot_row) = c("Control", "Treatment")
# all miRNAs
plot_heatmap_annot(correlation_Status_Timepoints, expr_sort$miRNA , plots_props, "Status", annot_row, colors, sprintf("%s/Time_Status", output_folder))
# only miRNAs for which one group fulfill the threshold
plot_heatmap_annot_filtered(correlation_Status_Timepoints, expr_sort$miRNA , plots_props, "Status", annot_row, miRNA_list, colors, sprintf("%s/Time_Status_filtered", output_folder))

# save
rownames(correlation_Status_Timepoints) = expr$miRNA
fwrite(correlation_Status_Timepoints, sprintf("%s/supp_table_Timepoint_Status.csv", supp_table_output_folder), sep = "\t", row.names = TRUE)
write.xlsx(correlation_Status_Timepoints, sprintf("%s/supp_table_Timepoint_Status.xlsx", supp_table_output_folder), colNames = TRUE, rowNames = TRUE, append = FALSE)

# corr between miRNA and Timepoints for every ear group
correlation_Group_Timepoints = c()
for(group in unique(annot$Group)){
  groups_ID = annot[annot$Group == group,]$ID
  tmp = colnames(expr) %in% groups_ID
  expr_sort = data.frame(miRNA=expr$miRNA, expr[,..tmp], check.names = FALSE)
  timepoints = annot[annot$ID %in% groups_ID,]$Time
  timepoints = as.numeric(gsub("m", "", timepoints))
  correlation_Group_Timepoints[[group]] = as.data.frame(cor(t(expr_sort[,2:dim(expr_sort)[2]]), timepoints, method ="spearman"))
}

correlation_Group_Timepoints = as.data.frame(correlation_Group_Timepoints)
colnames(correlation_Group_Timepoints) = unique(annot$Group)
annot_row = data.frame(Group=c("Control_left ear", "Treatment_left ear", "Treatment_right ear", "Control_right ear"))
rownames(annot_row) = c("Control_left ear", "Treatment_left ear", "Treatment_right ear", "Control_right ear")
# all miRNAs
plot_heatmap_annot(correlation_Group_Timepoints, expr$miRNA , plots_props, "Group", annot_row, colors, sprintf("%s/Timepoint_Group", output_folder))
# only miRNAs for which one group fulfill the threshold
plot_heatmap_annot_filtered(correlation_Group_Timepoints, expr$miRNA , plots_props, "Group", annot_row, miRNA_list, colors, sprintf("%s/Timepoint_Group_filtered", output_folder))

# save
rownames(correlation_Group_Timepoints) = expr$miRNA
fwrite(correlation_Group_Timepoints, sprintf("%s/supp_table_Timepoint_Group.csv", supp_table_output_folder), sep = "\t", row.names = TRUE)
write.xlsx(correlation_Group_Timepoints, sprintf("%s/supp_table_Timepoint_Group.xlsx", supp_table_output_folder), colNames = TRUE, rowNames = TRUE, append = FALSE)

# combine all corr to one dataframe
combined_corr = cbind(correlation_Status_Timepoints,correlation_Group_Timepoints)

# save 
rownames(combined_corr) = expr$miRNA
fwrite(combined_corr, sprintf("%s/supp_table_all.csv", supp_table_output_folder), sep = "\t", row.names = TRUE)
write.xlsx(combined_corr, sprintf("%s/supp_table_all.xlsx", supp_table_output_folder), colNames = TRUE, rowNames = TRUE, append = FALSE)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
