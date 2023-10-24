suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))

data.table::setDTthreads(1)

set.seed(snakemake@params$parameters_porps$set_seed)

print("heatmap pathways")

label_split = function(string, collapse = "\n", length=25) {
  if (nchar(string) > length) {
    return(paste(strwrap(string, width = length), collapse = collapse))
  } else {
    return(string)
  }
}

plot_heatmap = function(df, row_names, plot_props, width, height, filename){
  df = as.matrix(df)
  rownames_without_species = row_names
  rownames(df) = rownames_without_species
  
  legend_ticks = seq(from = 0, to = max(t(df)), length.out = 5)
  legend_tick_labels = round(legend_ticks, digits = 2)

  maxi = max(df)
  col_fun = colorRamp2(c(0,maxi), hcl_palette = "Reds 3", reverse = TRUE)
  col_fun_legend = col_fun
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (df[i, j] >= -log10(0.001)){ 
      grid.text("***", x, y, gp = gpar(fontsize = 6, col = "white"))
    } else if (df[i, j] >= -log10(0.01)){
      grid.text("**", x, y, gp = gpar(fontsize = 6))
    } else if (df[i, j] >= -log10(0.05)){
      grid.text("*", x, y, gp = gpar(fontsize = 6))
    }
  }

  p = ComplexHeatmap::Heatmap(df,
                              col = col_fun,
                              rect_gp = gpar(col = "white", lwd = .3),
                              cell_fun = cell_fun,
                              height = unit(height - 5.5, "cm"), width = unit(width - 2, "cm"),
                              row_names_side = "left", row_names_gp = gpar(fontsize = plot_props$font_size*2/3, fontfamily=plot_props$font_family), 
                              column_names_side = "bottom", column_names_gp = gpar(fontsize = plot_props$font_size*2/3, fontfamily=plot_props$font_family),
                              column_dend_reorder = FALSE, cluster_columns = FALSE,
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
  )
  legend = ComplexHeatmap::Legend(title = "Adj. p-value from Hypergeometric test (-log10)", col_fun = col_fun_legend, at = legend_tick_labels, direction="horizontal", #col_fun = circlize::colorRamp2(legend_ticks, viridis(5))
                                  title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
  
  png(sprintf("%s.png",filename), width = width, height = height, units = plot_props$image_units, res = plot_props$dpi)
  draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, -1, 0.25), "cm"))
  draw(legend, x = unit(2.9, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s.svg",filename), width = width / 2.54, height = height / 2.54)
  draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, -2, 0.25), "cm"))
  draw(legend, x = unit(2.9, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
}

tbl = fread("data/external_data/mirpathdb/KEGG_prediction_union.csv", sep=',', header=T)

min_pathway = snakemake@params$filter_params$min_pathway
min_mirna = snakemake@params$filter_params$min_mirna
sig_th = snakemake@params$filter_params$sig_th

mirna_list = interesting_miRNA = fread(snakemake@input$miRNA_list)
params = snakemake@params$params
plot_props = snakemake@params$plot_props

results_folder = snakemake@params$results_folder
output_folder = sprintf("%s/heatmap_complex/pathway_pvalues", results_folder)
dir.create(output_folder, recursive=TRUE)


# filtering only interesting RNAs
tbl_filtered = tbl[tbl$RNA %in% mirna_list$miRNA, ]
# make binary
tbl_filtered_bin = as.matrix(tbl_filtered[, -1]) < sig_th
# sum over columns and remove all columns < min_mirna
tmp = colSums(tbl_filtered_bin) >= min_mirna
# get the column names
sig_pathways = names(tmp)[tmp == TRUE]

tbl_plot = tbl_filtered[, ..sig_pathways]
# log10
tbl_plot = -log10(tbl_plot)

pw_full = colnames(tbl_plot)

colnames(tbl_plot) = pw_full

# plotting 
plot_heatmap(tbl_plot, tbl_filtered$RNA, plot_props, width = 11, height = plot_props$image_height*1.3, filename = sprintf("%s/kegg_prediction_union", output_folder))


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
